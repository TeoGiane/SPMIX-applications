# # ---- Comparative Study on Synthetic Data - SPMIX Sampler ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_SPMIX", hide.opts = TRUE,
                         description = "Runs SPMIX sampler on synthetic data")
opt_parser <- add_argument(opt_parser, arg = "--input-dir", type = "character", default = "input",
                           help = "Relative path to the input data directory")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "output",
                           help = "Relative path to the output file")
opt_parser <- add_argument(opt_parser, arg = "--a-beta", type = "numeric", default = 2, 
                           help = "Beta prior parameter a")
opt_parser <- add_argument(opt_parser, arg = "--b-beta", type = "numeric", default = 36, 
                           help = "Beta prior parameter b")
extra_args <- parse_args(opt_parser)

# Set working directory relative to the script location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    # Running in RStudio
    setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
} else {
    # Running from command line
    initial.options <- commandArgs(trailingOnly = FALSE)
    script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
    setwd(dirname(dirname(script.name)))
}
cat("Setting working directory to: ", getwd(), "\n")

# Check if input directory exists
input_dir <- file.path(getwd(), extra_args$input_dir)
if (!dir.exists(input_dir)) {
  stop("Input directory: ", input_dir, " does not exist.")
}

# Generate output directory if it does not exist
output_dir <- file.path(getwd(), extra_args$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory: ", output_dir, "\n")
}

# Required Libraries
suppressMessages(library("SPMIX"))
suppressMessages(library("sf"))
suppressMessages(library("rjags"))

# Load data
load(file.path(input_dir, "data.dat"))

# Load shapefile
grid_sf <- st_read(file.path(input_dir, "shp", "grid.shp"), quiet = TRUE)

# Compute summary statistics for MCAR model
Y <- c(sapply(data, quantile, 0.05),
       sapply(data, quantile, 0.25),
       sapply(data, quantile, 0.5),
       sapply(data, quantile, 0.75),
       sapply(data, quantile, 0.95))
names(Y) <- NULL

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid_sf, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Prepare data for JAGS
naiveMCAR_data <- list("Y" = Y,
                       "I" = length(data),
                       "rho" = 0.95,
                       "W" = W,
                       "Ncov" = 5,
                       "mu0" = mean(Y),
                       "abeta" = extra_args$a_beta,
                       "bbeta" = extra_args$b_beta)

# Compile model
naiveMCAR_model <- jags.model(file.path(getwd(), "src", "bug", "naiveMCAR_marginalized.bug"),
                              data = naiveMCAR_data, n.chains=1, n.adapt = 0)

# Run MCMC
naiveMCAR_params <- c("G", "p", "one_over_sigmasq", "one_over_tausq")
naiveMCAR_fit <- jags.samples(naiveMCAR_model, variable.names = naiveMCAR_params, n.iter = 1000, thin = 1)

# Store naiveMCAR_fit
if (exists("naiveMCAR_fit")) {
  filename <- file.path(output_dir, "naiveMCAR-fit.dat")
  save(naiveMCAR_fit, file = filename)
  cat(sprintf("Saved naiveMCAR chain at: %s\n", filename))
}
