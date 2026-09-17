rm(list = ls())
library(dplyr)
#install.packages("caper")
library(ggplot2)
library(ape)
library(caper)
# read in from tab delim total anage database file
read.delim(file = "data/anage_data.txt") -> anage


for (a in c("Amphibia", "Aves", "Reptilia","Mammalia")) {
  print(a)
  anage %>% filter(Class == a) -> tmp
  
  tmp[which(!is.na(tmp$Maximum.longevity..yrs.)), ] -> tmp
  tmp[which(!is.na(tmp$Adult.weight..g.)), ] -> tmp
  
  tmp[,c("Genus","Species")] -> tmp.species
  
  write.table(tmp.species, file = "output/mamms.species.txt",
              quote = F, col.names = F, row.names = F)

  for (clade in c("output/ams.species.nwk","output/aves.species.nwk",
                  "output/rept.species.nwk","output/mamms.species.nwk")) {
    
    tree <- read.tree(clade)
  }
 
  
  log10(tmp$Adult.weight..g.) -> tmp$log_adult_weight
  log10(tmp$Maximum.longevity..yrs.) -> tmp$log_max_long
  
  tmp[,c("Genus","Species","log_adult_weight","log_max_long")] -> tmp
  
  tree$tip.label
  
  paste(tmp$Genus, tmp$Species, sep = "_") -> tmp$Species
  
  tmp[,c("Species","log_adult_weight","log_max_long")] -> tmp
  
  comp <- comparative.data(
    phy = tree,
    data = tmp,
    names.col = "Species",
    vcv = TRUE,
    warn.dropped = TRUE
  )
  
  fit_lambda_near0 <- pgls(
    log_max_long ~ log_adult_weight,
    data = comp,
    lambda = 1e-6
  )
  
  fit_lambda_half <- pgls(
    log_max_long ~ log_adult_weight,
    data = comp,
    lambda = 0.5
  )
  
  fit_lambda1 <- pgls(
    log_max_long ~ log_adult_weight,
    data = comp,
    lambda = 1
  )
  
  summary(fit_lambda_near0)
  summary(fit_lambda_half)
  summary(fit_lambda1)
  
 
  results <- data.frame(
    lambda = c(1e-6, 0.5, 1),
    slope = c(
      coef(fit_lambda_near0)["log_adult_weight"],
      coef(fit_lambda_half)["log_adult_weight"],
      coef(fit_lambda1)["log_adult_weight"]
    ),
    slope_SE = c(
      summary(fit_lambda_near0)$coefficients[
        "log_adult_weight", "Std. Error"
      ],
      summary(fit_lambda_half)$coefficients[
        "log_adult_weight", "Std. Error"
      ],
      summary(fit_lambda1)$coefficients[
        "log_adult_weight", "Std. Error"
      ]
    ),
    p_value = c(
      summary(fit_lambda_near0)$coefficients[
        "log_adult_weight", "Pr(>|t|)"
      ],
      summary(fit_lambda_half)$coefficients[
        "log_adult_weight", "Pr(>|t|)"
      ],
      summary(fit_lambda1)$coefficients[
        "log_adult_weight", "Pr(>|t|)"
      ]
    ),
    AIC = c(
      AIC(fit_lambda_near0),
      AIC(fit_lambda_half),
      AIC(fit_lambda1)
    )
  )
  
  results

  plot(tmp$log_adult_weight, tmp$log_max_long, main = a)
  abline(fit_lambda1, col = "firebrick", lwd = 2)
  
}

