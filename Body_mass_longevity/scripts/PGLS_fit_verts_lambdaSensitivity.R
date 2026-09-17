rm(list = ls())

library(dplyr)
library(ggplot2)
library(ape)
library(caper)

# ============================================================
# Input files and settings
# ============================================================

anage <- read.delim(
  file = "data/anage_data.txt",
  stringsAsFactors = FALSE
)

tree_files <- c(
  Amphibia = "output/ams.species.nwk",
  Aves = "output/aves.species.nwk",
  Reptilia = "output/rept.species.nwk",
  Mammalia = "output/mamms.species.nwk"
)

species_files <- c(
  Amphibia = "output/ams.species.txt",
  Aves = "output/aves.species.txt",
  Reptilia = "output/rept.species.txt",
  Mammalia = "output/mamms.species.txt"
)

classes <- names(tree_files)

lambda_near0_value <- 0.001

all_pgls_results <- list()
all_model_comparisons <- list()

# ============================================================
# Helper functions
# ============================================================

fit_pgls_safely <- function(formula, comp, lambda_value, class_name) {
  
  tryCatch(
    caper::pgls(
      formula = formula,
      data = comp,
      lambda = lambda_value
    ),
    error = function(e) {
      message(
        class_name,
        ", lambda = ",
        lambda_value,
        " PGLS fit failed: ",
        e$message
      )
      NULL
    }
  )
}

extract_model_row <- function(
    fit,
    class_name,
    model_name,
    lambda_value,
    n_species
) {
  
  fit_summary <- summary(fit)
  
  data.frame(
    class = class_name,
    model = model_name,
    lambda = lambda_value,
    n_species = n_species,
    intercept = unname(coef(fit)["(Intercept)"]),
    slope = unname(coef(fit)["log_adult_weight"]),
    slope_SE = fit_summary$coefficients[
      "log_adult_weight",
      "Std. Error"
    ],
    t_value = fit_summary$coefficients[
      "log_adult_weight",
      "t value"
    ],
    p_value = fit_summary$coefficients[
      "log_adult_weight",
      "Pr(>|t|)"
    ],
    AIC = AIC(fit),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Main analysis loop
# ============================================================

for (a in classes) {
  
  cat("\n====================================================\n")
  cat("Analyzing:", a, "\n")
  cat("====================================================\n")
  
  # ----------------------------------------------------------
  # Prepare ANAGE data
  # ----------------------------------------------------------
  
  tmp <- anage %>%
    dplyr::filter(
      Class == a,
      !is.na(Maximum.longevity..yrs.),
      !is.na(Adult.weight..g.),
      Maximum.longevity..yrs. > 0,
      Adult.weight..g. > 0
    ) %>%
    dplyr::mutate(
      Species = paste(Genus, Species, sep = "_"),
      log_adult_weight = log10(Adult.weight..g.),
      log_max_long = log10(Maximum.longevity..yrs.)
    ) %>%
    dplyr::select(
      Genus,
      Species,
      log_adult_weight,
      log_max_long
    )
  
  # Write original genus and species epithet for external tools.
  tmp.species <- tmp %>%
    dplyr::mutate(
      Species_epithet = sub("^[^_]+_", "", Species)
    ) %>%
    dplyr::select(Genus, Species_epithet)
  
  write.table(
    tmp.species,
    file = species_files[a],
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Read tree and match species
  # ----------------------------------------------------------
  
  if (!file.exists(tree_files[a])) {
    warning(
      a,
      ": tree file not found: ",
      tree_files[a],
      ". Skipping."
    )
    next
  }
  
  tree <- ape::read.tree(tree_files[a])
  
  if (is.null(tree$tip.label) || length(tree$tip.label) < 4) {
    warning(a, ": tree has fewer than four tips. Skipping.")
    next
  }
  
  cat("Tree file:", tree_files[a], "\n")
  cat("ANAGE observations:", nrow(tmp), "\n")
  cat("Tree tips:", length(tree$tip.label), "\n")
  cat(
    "ANAGE species missing from tree:",
    sum(!tmp$Species %in% tree$tip.label),
    "\n"
  )
  cat(
    "Tree tips missing from ANAGE:",
    sum(!tree$tip.label %in% tmp$Species),
    "\n"
  )
  
  shared_species <- intersect(tmp$Species, tree$tip.label)
  
  if (length(shared_species) < 4) {
    warning(a, ": fewer than four matched species. Skipping.")
    next
  }
  
  # Retain one record for each matched species.
  tmp_matched <- tmp %>%
    dplyr::filter(Species %in% shared_species) %>%
    dplyr::distinct(Species, .keep_all = TRUE) %>%
    dplyr::select(Species, log_adult_weight, log_max_long)
  
  tree_matched <- ape::keep.tip(
    phy = tree,
    tip = tmp_matched$Species
  )
  
  # Order data in exactly the same order as tree tips.
  tmp_matched <- tmp_matched[
    match(tree_matched$tip.label, tmp_matched$Species),
    ,
    drop = FALSE
  ]
  
  rownames(tmp_matched) <- tmp_matched$Species
  
  cat("Matched taxa:", length(tree_matched$tip.label), "\n")
  
  if (length(tree_matched$tip.label) < 4) {
    warning(a, ": fewer than four tree tips remain after pruning. Skipping.")
    next
  }
  
  if (anyDuplicated(tmp_matched$Species) > 0) {
    warning(a, ": duplicated species in matched data. Skipping.")
    next
  }
  
  if (anyDuplicated(tree_matched$tip.label) > 0) {
    warning(a, ": duplicated tree-tip labels. Skipping.")
    next
  }
  
  # ----------------------------------------------------------
  # Repair topology and branch lengths
  # ----------------------------------------------------------
  
  used_grafen <- FALSE
  
  if (is.null(tree_matched$edge.length)) {
    warning(
      a,
      ": tree has no branch lengths; assigning Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
  }
  
  if (any(!is.finite(tree_matched$edge.length))) {
    warning(
      a,
      ": tree has non-finite branch lengths; assigning Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
  }
  
  # Turn polytomies into deterministic bifurcations.
  if (!ape::is.binary.tree(tree_matched)) {
    cat("Resolving polytomies in", a, "tree.\n")
    
    tree_matched <- ape::multi2di(
      phy = tree_matched,
      random = FALSE
    )
  }
  
  # Replace nonpositive branches with a small positive value.
  positive_lengths <- tree_matched$edge.length[
    tree_matched$edge.length > 0 &
      is.finite(tree_matched$edge.length)
  ]
  
  if (length(positive_lengths) == 0) {
    warning(
      a,
      ": no valid positive branch lengths; assigning Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
    
    positive_lengths <- tree_matched$edge.length[
      tree_matched$edge.length > 0 &
        is.finite(tree_matched$edge.length)
    ]
  }
  
  min_positive <- min(positive_lengths)
  
  tree_matched$edge.length[
    !is.finite(tree_matched$edge.length) |
      tree_matched$edge.length <= 0
  ] <- min_positive * 1e-6
  
  cat("Tree is binary:", ape::is.binary.tree(tree_matched), "\n")
  cat(
    "Nonpositive edges after repair:",
    sum(tree_matched$edge.length <= 0),
    "\n"
  )
  cat("Used Grafen branch lengths:", used_grafen, "\n")
  
  # ----------------------------------------------------------
  # Create comparative object and validate covariance matrix
  # ----------------------------------------------------------
  
  comp <- tryCatch(
    caper::comparative.data(
      phy = tree_matched,
      data = tmp_matched,
      names.col = "Species",
      vcv = TRUE,
      warn.dropped = TRUE
    ),
    error = function(e) {
      message(a, ": comparative.data failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(comp)) {
    next
  }
  
  V <- comp$vcv
  V_rank <- qr(V)$rank
  V_size <- nrow(V)
  
  cat("VCV matrix:", V_size, "x", V_size, "\n")
  cat("VCV rank:", V_rank, "\n")
  cat("VCV rank deficiency:", V_size - V_rank, "\n")
  
  # Retry with Grafen branch lengths if original lengths give singular VCV.
  if (V_rank < V_size) {
    
    warning(
      a,
      ": singular covariance matrix; retrying with Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
    
    comp <- tryCatch(
      caper::comparative.data(
        phy = tree_matched,
        data = tmp_matched,
        names.col = "Species",
        vcv = TRUE,
        warn.dropped = TRUE
      ),
      error = function(e) {
        message(
          a,
          ": comparative.data failed after Grafen lengths: ",
          e$message
        )
        NULL
      }
    )
    
    if (is.null(comp)) {
      next
    }
    
    V <- comp$vcv
    V_rank <- qr(V)$rank
    V_size <- nrow(V)
    
    cat("Grafen VCV rank:", V_rank, "\n")
    cat("Grafen VCV rank deficiency:", V_size - V_rank, "\n")
    
    if (V_rank < V_size) {
      warning(
        a,
        ": covariance matrix remains singular after Grafen transformation. ",
        "Skipping."
      )
      next
    }
  }
  
  # ----------------------------------------------------------
  # Fit OLS and PGLS models
  # ----------------------------------------------------------
  
  fit_ols <- lm(
    log_max_long ~ log_adult_weight,
    data = tmp_matched
  )
  
  fit_lambda_near0 <- fit_pgls_safely(
    formula = log_max_long ~ log_adult_weight,
    comp = comp,
    lambda_value = lambda_near0_value,
    class_name = a
  )
  
  fit_lambda_half <- fit_pgls_safely(
    formula = log_max_long ~ log_adult_weight,
    comp = comp,
    lambda_value = 0.5,
    class_name = a
  )
  
  fit_lambda1 <- fit_pgls_safely(
    formula = log_max_long ~ log_adult_weight,
    comp = comp,
    lambda_value = 1,
    class_name = a
  )
  
  # Optional: estimate lambda from the data by maximum likelihood.
  fit_lambda_ml <- tryCatch(
    caper::pgls(
      formula = log_max_long ~ log_adult_weight,
      data = comp,
      lambda = "ML",
      bounds = list(lambda = c(0.001, 1))
    ),
    error = function(e) {
      message(a, ", lambda = ML PGLS fit failed: ", e$message)
      NULL
    }
  )
  
  # ----------------------------------------------------------
  # Store models and identify successful PGLS fits
  # ----------------------------------------------------------
  
  model_list <- list(
    lambda_near0 = fit_lambda_near0,
    lambda_half = fit_lambda_half,
    lambda_1 = fit_lambda1,
    lambda_ML = fit_lambda_ml
  )
  
  fixed_lambda_values <- c(
    lambda_near0 = lambda_near0_value,
    lambda_half = 0.5,
    lambda_1 = 1
  )
  
  successful_models <- !vapply(
    model_list,
    is.null,
    logical(1)
  )
  
  cat(
    "Successful PGLS models:",
    paste(names(model_list)[successful_models], collapse = ", "),
    "\n"
  )
  
  # ----------------------------------------------------------
  # Build OLS + successful-PGLS comparison table
  # ----------------------------------------------------------
  
  ols_summary <- summary(fit_ols)
  
  comparison <- data.frame(
    class = a,
    branch_lengths = ifelse(
      used_grafen,
      "Grafen_topology_derived",
      "source_tree"
    ),
    model = "OLS",
    lambda = NA_real_,
    n_species = nrow(tmp_matched),
    intercept = unname(coef(fit_ols)["(Intercept)"]),
    slope = unname(coef(fit_ols)["log_adult_weight"]),
    slope_SE = ols_summary$coefficients[
      "log_adult_weight",
      "Std. Error"
    ],
    t_value = ols_summary$coefficients[
      "log_adult_weight",
      "t value"
    ],
    p_value = ols_summary$coefficients[
      "log_adult_weight",
      "Pr(>|t|)"
    ],
    AIC = AIC(fit_ols),
    stringsAsFactors = FALSE
  )
  
  if (any(successful_models)) {
    
    for (model_name in names(model_list)[successful_models]) {
      
      fit <- model_list[[model_name]]
      fit_summary <- summary(fit)
      
      lambda_for_table <- if (model_name == "lambda_ML") {
        unname(fit$param["lambda"])
      } else {
        unname(fixed_lambda_values[model_name])
      }
      
      comparison <- rbind(
        comparison,
        data.frame(
          class = a,
          branch_lengths = ifelse(
            used_grafen,
            "Grafen_topology_derived",
            "source_tree"
          ),
          model = model_name,
          lambda = lambda_for_table,
          n_species = length(comp$phy$tip.label),
          intercept = unname(coef(fit)["(Intercept)"]),
          slope = unname(coef(fit)["log_adult_weight"]),
          slope_SE = fit_summary$coefficients[
            "log_adult_weight",
            "Std. Error"
          ],
          t_value = fit_summary$coefficients[
            "log_adult_weight",
            "t value"
          ],
          p_value = fit_summary$coefficients[
            "log_adult_weight",
            "Pr(>|t|)"
          ],
          AIC = AIC(fit),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  print(comparison)
  
  all_model_comparisons[[a]] <- comparison
  
  # Keep a PGLS-only table for downstream analyses.
  if (any(successful_models)) {
    all_pgls_results[[a]] <- comparison %>%
      dplyr::filter(model != "OLS")
  }
  
  # ----------------------------------------------------------
  # Plot raw data plus successful model lines
  # ----------------------------------------------------------
  
  x_values <- seq(
    from = min(tmp_matched$log_adult_weight),
    to = max(tmp_matched$log_adult_weight),
    length.out = 200
  )
  
  prediction_data <- data.frame(
    log_adult_weight = x_values
  )
  
  # Start with OLS predictions to establish y-axis range.
  y_lines <- predict(
    fit_ols,
    newdata = prediction_data
  )
  
  # Add successful PGLS prediction ranges to y-axis calculation.
  for (model_name in names(model_list)[successful_models]) {
    y_lines <- c(
      y_lines,
      predict(
        model_list[[model_name]],
        newdata = prediction_data
      )
    )
  }
  
  plot(
    tmp_matched$log_adult_weight,
    tmp_matched$log_max_long,
    main = paste(a, "- longevity vs adult mass"),
    xlab = "log10 adult weight (g)",
    ylab = "log10 maximum longevity (years)",
    ylim = range(c(tmp_matched$log_max_long, y_lines), finite = TRUE),
    pch = 16,
    col = grDevices::adjustcolor("gray25", alpha.f = 0.55)
  )
  
  # OLS: black dashed line.
  lines(
    x_values,
    predict(fit_ols, newdata = prediction_data),
    col = "black",
    lwd = 2,
    lty = 2
  )
  
  # PGLS line colors.
  model_colors <- c(
    lambda_near0 = "steelblue",
    lambda_half = "darkorange3",
    lambda_1 = "firebrick",
    lambda_ML = "darkgreen"
  )
  
  model_line_widths <- c(
    lambda_near0 = 2,
    lambda_half = 2,
    lambda_1 = 3,
    lambda_ML = 2
  )
  
  # Add only fitted PGLS models.
  for (model_name in names(model_list)[successful_models]) {
    lines(
      x_values,
      predict(
        model_list[[model_name]],
        newdata = prediction_data
      ),
      col = model_colors[model_name],
      lwd = model_line_widths[model_name],
      lty = 1
    )
  }
  
  # Build legend only from models that were successfully fitted.
  legend_labels <- c("OLS")
  
  if ("lambda_near0" %in% names(model_list)[successful_models]) {
    legend_labels <- c(
      legend_labels,
      paste0("PGLS: lambda = ", lambda_near0_value)
    )
  }
  
  if ("lambda_half" %in% names(model_list)[successful_models]) {
    legend_labels <- c(
      legend_labels,
      "PGLS: lambda = 0.5"
    )
  }
  
  if ("lambda_1" %in% names(model_list)[successful_models]) {
    legend_labels <- c(
      legend_labels,
      "PGLS: lambda = 1"
    )
  }
  
  if ("lambda_ML" %in% names(model_list)[successful_models]) {
    lambda_ml_value <- unname(fit_lambda_ml$param["lambda"])
    
    legend_labels <- c(
      legend_labels,
      paste0(
        "PGLS: lambda ML = ",
        round(lambda_ml_value, 3)
      )
    )
  }
  
  legend_colors <- c("black")
  
  if ("lambda_near0" %in% names(model_list)[successful_models]) {
    legend_colors <- c(legend_colors, "steelblue")
  }
  
  if ("lambda_half" %in% names(model_list)[successful_models]) {
    legend_colors <- c(legend_colors, "darkorange3")
  }
  
  if ("lambda_1" %in% names(model_list)[successful_models]) {
    legend_colors <- c(legend_colors, "firebrick")
  }
  
  if ("lambda_ML" %in% names(model_list)[successful_models]) {
    legend_colors <- c(legend_colors, "darkgreen")
  }
  
  legend_lty <- c(2, rep(1, length(legend_labels) - 1))
  legend_lwd <- c(2, rep(2, length(legend_labels) - 1))
  
  if ("lambda_1" %in% names(model_list)[successful_models]) {
    lambda1_position <- which(legend_labels == "PGLS: lambda = 1")
    
    legend_lwd[lambda1_position] <- 3
  }
  
  legend(
    "topleft",
    legend = legend_labels,
    col = legend_colors,
    lty = legend_lty,
    lwd = legend_lwd,
    bty = "n",
    cex = 0.8
  )
}

# ============================================================
# Combine and save final results
# ============================================================

if (length(all_model_comparisons) > 0) {
  
  all_comparison_df <- dplyr::bind_rows(all_model_comparisons)
  
  write.csv(
    all_comparison_df,
    file = "output/OLS_and_PGLS_longevity_adult_weight_comparison.csv",
    row.names = FALSE
  )
  
  print(all_comparison_df)
  
} else {
  
  warning("No model-comparison results were generated.")
}

if (length(all_pgls_results) > 0) {
  
  all_pgls_df <- dplyr::bind_rows(all_pgls_results)
  
  write.csv(
    all_pgls_df,
    file = "output/PGLS_longevity_adult_weight_results.csv",
    row.names = FALSE
  )
  
  print(all_pgls_df)
  
} else {
  
  warning("No successful PGLS models were generated.")
}