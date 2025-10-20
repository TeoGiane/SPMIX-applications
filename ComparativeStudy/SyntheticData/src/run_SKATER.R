# # ---- Comparative Study on Synthetic Data - SKATER Algorithm ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_SKATER", hide.opts = TRUE,
                         description = "Runs SKATER algorithm on synthetic data")
opt_parser <- add_argument(opt_parser, arg = "--input-dir", type = "character", default = "input",
                           help = "Relative path to the input data directory")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "output",
                           help = "Relative path to the output file")
opt_parser <- add_argument(opt_parser, arg = "--n-cuts", type = "integer", default = 10,
                           help = "Number of cuts for SKATER algorithm")
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
suppressMessages(library("spdep"))

# Load data
load(file.path(input_dir, "data.dat"))

# Load shapefile
grid_sf <- st_read(file.path(input_dir, "shp", "grid.shp"), quiet = TRUE)

# Compute edges matrix for SKATER algorithm
adj_list <- poly2nb(grid_sf, queen = FALSE)
edges <- mstree(nb2listw(adj_list, style = "B"), 1)[, 1:2]

# Compute summary statistics for SKATER algorithm
Y <- c(sapply(data, quantile, 0.05),
       sapply(data, quantile, 0.25),
       sapply(data, quantile, 0.5),
       sapply(data, quantile, 0.75),
       sapply(data, quantile, 0.95))
names(Y) <- NULL

# Run SKATER algorithm
SKATER_fit <- skater(edges = edges, data = Y, ncuts = extra_args$n_cuts, crit = 2)

# Store SKATER_fit
if (exists("SKATER_fit")) {
    filename <- file.path(output_dir, "SKATER-fit.dat")
    save(SKATER_fit, file = filename)
    cat(sprintf("Saved SKATER output at: %s\n", filename))
}
