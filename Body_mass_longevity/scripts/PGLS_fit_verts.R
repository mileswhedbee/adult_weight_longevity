rm(list = ls())

library(dplyr)
library(ggplot2)
library(ape)
library(caper)

# ============================================================
# Input data and file paths
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

# Results containers
all_results <- list()
all_missing_species <- list()

# ============================================================
# Helper function: safely fit ML-lambda PGLS
# ============================================================

fit_ml_pgls <- function(formula, comp, class_name) {
  
  tryCatch(
    caper::pgls(
      formula = formula,
      data = comp,
      lambda = "ML",
      bounds = list(lambda = c(0.001, 1))
    ),
    error = function(e) {
      message(
        class_name,
        ": ML-lambda PGLS fit failed: ",
        e$message
      )
      NULL
    }
  )
}

# ============================================================
# Main loop: one analysis per vertebrate class
# ============================================================

for (a in classes) {
  
  cat("\n====================================================\n")
  cat("Analyzing:", a, "\n")
  cat("====================================================\n")
  
  # ----------------------------------------------------------
  # Prepare lifespan/body-mass data
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
  
  # Write a two-column genus/species file for phylogeny tools.
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
  # Read matching tree
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
  
  # ----------------------------------------------------------
  # Match tree tips and ANAGE species
  # ----------------------------------------------------------
  
  shared_species <- intersect(tmp$Species, tree$tip.label)
  
  missing_from_tree <- setdiff(tmp$Species, tree$tip.label)
  missing_from_data <- setdiff(tree$tip.label, tmp$Species)
  
  all_missing_species[[a]] <- list(
    data_species_missing_from_tree = missing_from_tree,
    tree_tips_missing_from_data = missing_from_data
  )
  
  cat("Tree file:", tree_files[a], "\n")
  cat("ANAGE records before matching:", nrow(tmp), "\n")
  cat("Tree tips before matching:", length(tree$tip.label), "\n")
  cat("Matched species:", length(shared_species), "\n")
  cat("ANAGE species absent from tree:", length(missing_from_tree), "\n")
  cat("Tree tips absent from ANAGE:", length(missing_from_data), "\n")
  
  if (length(shared_species) < 4) {
    warning(a, ": fewer than four matched species. Skipping.")
    next
  }
  
  # Retain a single ANAGE row per species.
  tmp_matched <- tmp %>%
    dplyr::filter(Species %in% shared_species) %>%
    dplyr::distinct(Species, .keep_all = TRUE) %>%
    dplyr::select(
      Species,
      log_adult_weight,
      log_max_long
    )
  
  # Prune phylogeny to exactly the species retained in the data.
  tree_matched <- ape::keep.tip(
    phy = tree,
    tip = tmp_matched$Species
  )
  
  # Put data in tree-tip order.
  tmp_matched <- tmp_matched[
    match(tree_matched$tip.label, tmp_matched$Species),
    ,
    drop = FALSE
  ]
  
  rownames(tmp_matched) <- tmp_matched$Species
  
  if (anyDuplicated(tmp_matched$Species) > 0) {
    warning(a, ": duplicate species in matched data. Skipping.")
    next
  }
  
  if (anyDuplicated(tree_matched$tip.label) > 0) {
    warning(a, ": duplicate labels in tree. Skipping.")
    next
  }
  
  if (
    any(!is.finite(tmp_matched$log_adult_weight)) ||
    any(!is.finite(tmp_matched$log_max_long))
  ) {
    warning(a, ": non-finite trait values after transformation. Skipping.")
    next
  }
  
  cat("Matched taxa used:", nrow(tmp_matched), "\n")
  
  # ----------------------------------------------------------
  # Validate tree topology and branch lengths
  # ----------------------------------------------------------
  
  used_grafen <- FALSE
  
  # Supply topology-based branch lengths if the tree has none.
  if (is.null(tree_matched$edge.length)) {
    
    warning(
      a,
      ": tree has no branch lengths; using Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
  }
  
  # Replace non-finite tree lengths entirely, rather than selectively.
  if (any(!is.finite(tree_matched$edge.length))) {
    
    warning(
      a,
      ": tree has non-finite branch lengths; using Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
  }
  
  # Resolve polytomies so that the tree is bifurcating.
  if (!ape::is.binary.tree(tree_matched)) {
    
    cat("Resolving polytomies in", a, "tree.\n")
    
    tree_matched <- ape::multi2di(
      phy = tree_matched,
      random = FALSE
    )
  }
  
  # Assign tiny positive lengths to zero or negative edges.
  positive_lengths <- tree_matched$edge.length[
    is.finite(tree_matched$edge.length) &
      tree_matched$edge.length > 0
  ]
  
  if (length(positive_lengths) == 0) {
    
    warning(
      a,
      ": no positive branch lengths; using Grafen branch lengths."
    )
    
    tree_matched <- ape::compute.brlen(
      phy = tree_matched,
      method = "Grafen",
      power = 1
    )
    
    used_grafen <- TRUE
    
    positive_lengths <- tree_matched$edge.length[
      tree_matched$edge.length > 0
    ]
  }
  
  min_positive <- min(positive_lengths)
  
  tree_matched$edge.length[
    !is.finite(tree_matched$edge.length) |
      tree_matched$edge.length <= 0
  ] <- min_positive * 1e-6
  
  cat("Tree is binary:", ape::is.binary.tree(tree_matched), "\n")
  cat(
    "Nonpositive edge lengths after repair:",
    sum(tree_matched$edge.length <= 0),
    "\n"
  )
  cat("Used Grafen branch lengths:", used_grafen, "\n")
  
  # ----------------------------------------------------------
  # Build comparative data and check phylogenetic covariance
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
  
  cat("VCV dimensions:", V_size, "x", V_size, "\n")
  cat("VCV rank:", V_rank, "\n")
  cat("VCV rank deficiency:", V_size - V_rank, "\n")
  
  # If source branch lengths create a singular VCV, retry using Grafen lengths.
  if (V_rank < V_size) {
    
    warning(
      a,
      ": source-tree covariance is singular; retrying with Grafen lengths."
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
        ": VCV is still singular after Grafen transformation. Skipping."
      )
      next
    }
  }
  
  # ----------------------------------------------------------
  # Fit models: OLS plus ML-lambda PGLS
  # ----------------------------------------------------------
  
  fit_ols <- lm(
    log_max_long ~ log_adult_weight,
    data = tmp_matched
  )
  
  fit_lambda_ml <- fit_ml_pgls(
    formula = log_max_long ~ log_adult_weight,
    comp = comp,
    class_name = a
  )
  
  # ----------------------------------------------------------
  # Extract model summaries
  # ----------------------------------------------------------
  
  ols_summary <- summary(fit_ols)
  
  results <- data.frame(
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
  
  if (!is.null(fit_lambda_ml)) {
    
    ml_summary <- summary(fit_lambda_ml)
    
    results <- rbind(
      results,
      data.frame(
        class = a,
        branch_lengths = ifelse(
          used_grafen,
          "Grafen_topology_derived",
          "source_tree"
        ),
        model = "PGLS_lambda_ML",
        lambda = unname(fit_lambda_ml$param["lambda"]),
        n_species = length(comp$phy$tip.label),
        intercept = unname(coef(fit_lambda_ml)["(Intercept)"]),
        slope = unname(coef(fit_lambda_ml)["log_adult_weight"]),
        slope_SE = ml_summary$coefficients[
          "log_adult_weight",
          "Std. Error"
        ],
        t_value = ml_summary$coefficients[
          "log_adult_weight",
          "t value"
        ],
        p_value = ml_summary$coefficients[
          "log_adult_weight",
          "Pr(>|t|)"
        ],
        AIC = AIC(fit_lambda_ml),
        stringsAsFactors = FALSE
      )
    )
  }
  
  print(results)
  
  all_results[[a]] <- results
  
  # ----------------------------------------------------------
  # Plot raw matched data, OLS, and ML-lambda PGLS only
  # ----------------------------------------------------------
  
  x_values <- seq(
    from = min(tmp_matched$log_adult_weight),
    to = max(tmp_matched$log_adult_weight),
    length.out = 200
  )
  
  prediction_data <- data.frame(
    log_adult_weight = x_values
  )
  
  ols_predictions <- predict(
    fit_ols,
    newdata = prediction_data
  )
  
  if (!is.null(fit_lambda_ml)) {
    
    ml_predictions <- predict(
      fit_lambda_ml,
      newdata = prediction_data
    )
    
    lambda_ml_value <- unname(
      fit_lambda_ml$param["lambda"]
    )
    
    y_limits <- range(
      c(
        tmp_matched$log_max_long,
        ols_predictions,
        ml_predictions
      ),
      finite = TRUE
    )
    
  } else {
    
    ml_predictions <- NULL
    lambda_ml_value <- NA_real_
    
    y_limits <- range(
      c(
        tmp_matched$log_max_long,
        ols_predictions
      ),
      finite = TRUE
    )
  }
  
  plot(
    tmp_matched$log_adult_weight,
    tmp_matched$log_max_long,
    main = paste(a, "- longevity vs adult mass"),
    xlab = "log10 adult weight (g)",
    ylab = "log10 maximum longevity (years)",
    xlim = range(x_values),
    ylim = y_limits,
    pch = 16,
    col = grDevices::adjustcolor(
      "gray25",
      alpha.f = 0.55
    )
  )
  
  # OLS line: black dashed.
  lines(
    x_values,
    ols_predictions,
    col = "black",
    lwd = 2,
    lty = 2
  )
  
  # ML-lambda PGLS line: dark green solid.
  if (!is.null(fit_lambda_ml)) {
    
    lines(
      x_values,
      ml_predictions,
      col = "darkgreen",
      lwd = 3,
      lty = 1
    )
    
    legend(
      "topleft",
      legend = c(
        "OLS",
        paste0(
          "PGLS: lambda ML = ",
          round(lambda_ml_value, 3)
        )
      ),
      col = c("black", "darkgreen"),
      lty = c(2, 1),
      lwd = c(2, 3),
      bty = "n",
      cex = 0.85
    )
    
  } else {
    
    legend(
      "topleft",
      legend = c(
        "OLS",
        "PGLS: lambda ML failed"
      ),
      col = c("black", "gray50"),
      lty = c(2, 1),
      lwd = c(2, 2),
      bty = "n",
      cex = 0.85
    )
  }
}

# ============================================================
# Combine and save results
# ============================================================

if (length(all_results) > 0) {
  
  final_results <- dplyr::bind_rows(all_results)
  
  print(final_results)
  
  write.csv(
    final_results,
    file = "output/OLS_and_ML_PGLS_longevity_adult_weight_results.csv",
    row.names = FALSE
  )
  
} else {
  
  warning("No OLS/PGLS results were generated.")
}