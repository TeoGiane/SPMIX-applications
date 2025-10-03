# # ---- Comparative Study on Synthetic Data - CARBayes Sampler ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_CARBayes", hide.opts = TRUE,
                         description = "Runs CARBayes sampler on synthetic data")
opt_parser <- add_argument(opt_parser, arg = "--input-dir", type = "character", default = "input",
                           help = "Relative path to the input data directory")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "output",
                           help = "Relative path to the output file")
# opt_parser <- add_argument(opt_parser, arg = "--n-cuts", type = "integer", default = 10,
#                            help = "Number of cuts for SKATER algorithm")
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
suppressMessages(library("CARBayes"))

# Load data
load(file.path(input_dir, "data.dat"))

# Load shapefile
grid_sf <- st_read(file.path(input_dir, "shp", "grid.shp"), quiet = TRUE)

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid_sf, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Compute summary statistic for CARBayes model
CARBayes_data <- data.frame("median" = sapply(data, quantile, 0.5))

# Compute list of dissimilarities for CARBayes model
Z <- list("Z1" = as.matrix(dist(sapply(data, quantile, 0.05))),
          "Z2" = as.matrix(dist(sapply(data, quantile, 0.25))),
          "Z3" = as.matrix(dist(sapply(data, quantile, 0.75))),
          "Z4" = as.matrix(dist(sapply(data, quantile, 0.95))))

# Run CARBayes model
CARBayes_fit <- S.CARdissimilarity(median ~ 1, data = CARBayes_data, family = "gaussian", W = W, Z = Z, W.binary = TRUE,
                                   n.sample = 10000, burnin = 5000, thin = 1, verbose = TRUE)

# Store CARBayes_fit
if (exists("CARBayes_fit")) {
  filename <- file.path(output_dir, "CARBayes-fit.dat")
  save(CARBayes_fit, file = filename)
  cat(sprintf("Saved CARBayes chain at: %s\n", filename))
}
