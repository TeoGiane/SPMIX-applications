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
                  "reshape2",
                  "rjags",
                  "RProtoBuf",
                  "sf",
                  "spdep")

# Install packages from CRAN
install.packages(package_deps, repos = repos)

# Register stadia map API key for ggmap
ggmap::register_stadiamaps(key = "1f3e613a-2adc-4c22-8bc9-901f1c33e05b",
                           write = TRUE)

# Install SPMIX from GitHub
devtools::install_github("TeoGiane/SPMIX")
