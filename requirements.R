# Choose repository
repos <- "https://cloud.r-project.org/"

# List of required packages
package_deps <- c("argparser",
                  "CARBayes",
                  "devtools",
                  "dplyr",
                  "ggplot2",
                  "ggmap",
                  "ggrepel",
                  "ggridges",
                  "gridExtra",
                  "loo",
                  "metRology",
                  "pROC",
                  "reshape2",
                  "rjags",
                  "RProtoBuf",
                  "sf",
                  "sn",
                  "spdep")

# Install packages from CRAN (in parallel)
num_cores <- parallel::detectCores() - 1
install.packages(package_deps, repos = repos, Ncpus = num_cores)
