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

# Load data
load(file.path(input_dir, "data.dat"))

# Load shapefile
grid_sf <- st_read(file.path(input_dir, "shp", "grid.shp"), quiet = TRUE)

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid_sf, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Setting MCMC parameters
burnin <- 5000
niter <- 5000
thin <- 1

# Set sampler parameters
params_template <-
  "
  num_components: 10

  p0_params {
    mu0: 0
    a: 2
    b: 2
    lam_: 0.1
  }

  rho {
    fixed: 0.95
  }

  sigma {
    inv_gamma_prior {
      alpha: 2
      beta: 2
    }
  }

  graph_params {
    beta_prior {
      a: %g
      b: %g
    }
  }
  "

# Specify beta parameters
params <- sprintf(params_template, extra_args$a_beta, extra_args$b_beta)

# Run SPMIX Sampler
SPMIX_fit <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params, type = "rjmcmc")

# Store SPMIX_fit
if (exists("SPMIX_fit")) {
  filename <- file.path(output_dir, "SPMIX-fit.dat")
  save(SPMIX_fit, file = filename)
  cat(sprintf("Saved SPMIX chain at: %s\n", filename))
}
