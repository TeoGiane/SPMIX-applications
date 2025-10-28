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
                  "metRology",
                  "reshape2",
                  "rjags",
                  "RProtoBuf",
                  "sf",
                  "sn",
                  "spdep")

# Install packages from CRAN
install.packages(package_deps, repos = repos)

# Install SPMIX from GitHub
devtools::install_github("TeoGiane/SPMIX")
