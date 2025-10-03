# Choose repository
repos <- "https://cloud.r-project.org/"

# List of required packages
package_deps <- c("devtools", "sf", "spdep", "CARBayes", "argparser", "rjags", "reshape2", "RProtoBuf")

# Install packages from CRAN
install.packages(package_deps, repos = repos)
# install.packages("CARBayes", repos = repos)
# install.packages("argparser", repos = repos, dependencies = FALSE)
# install.packages("rjags", repos = repos, dependencies = FALSE)
# install.packages("RProtoBuf", repos = repos, dependencies = FALSE)

# Install SPMIX from GitHub
devtools::install_github("TeoGiane/SPMIX")