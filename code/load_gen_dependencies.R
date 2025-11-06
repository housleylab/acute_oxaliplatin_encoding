########################### load general dependencies ########################### 
packagesS<-c("devtools",
             "dplyr",
             "parallel",
             "ggplot2",
             "readxl",
             "crayon",
             "PerformanceAnalytics",
             "Hmisc",
             "tidyr",
             "GGally",
             "tibble",
             "ggplot2",
             "leaps",
             "caret",
             "rstan",
             "rstanarm",
             "coda",
             "broom",
             "tidybayes",
             "emmeans",
             "loo",
             "bayesplot",
             "magrittr",
             "ggpubr",
             "modelr",
             "broom.mixed",
             "EnvStats",
             "ggrepel",
             "plotrix",
             "R.matlab",
             "data.table",
             "FactoMineR",
             "factoextra",
             "corrplot",
             "forcats",
             "tidyquant",
             "readr",
             "PK",
             "smplot2",
             "cowplot",
             "MASS",
             "stringr",
             "multcomp",
             "raster",
             "barcode",
             "forcats"
             
             
             
)

invisible(suppressWarnings(suppressMessages(package.check <- lapply(
  packagesS,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
))))



## do not exclude these options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rm(packagesS)
