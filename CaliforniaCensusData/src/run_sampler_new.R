# # ---- CALIFORNIA CENSUS DATA - RUN SAMPLER ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_sampler", hide.opts = TRUE,
                         description = "Run SPMIX::Sampler.BoundaryDetection on the California Census dataset")
opt_parser <- add_argument(opt_parser, arg = "input-file", type = "character",
                           help = "Relative path to input data file")
opt_parser <- add_argument(opt_parser, arg = "--num-components-prior", default = "shifted_poisson_prior { rate: 1.0 }",
                           help = "Prior on the number of components (in ASCII protocol buffer format)")
opt_parser <- add_argument(opt_parser, arg = "--rho-prior", default = "fixed: 0.95",
                           help = "Prior on the rho parameter (in ASCII protocol buffer format)")
opt_parser <- add_argument(opt_parser, arg = "--sigma-prior", default = "inv_gamma_prior { alpha: 6, beta: 4 }",
                           help = "Prior on the sigma parameter (in ASCII protocol buffer format)")
opt_parser <- add_argument(opt_parser, arg = "--graph-prior", default = "beta_prior { a: 2, b: 93 }",
                           help = "Prior on the graph parameters (in ASCII protocol buffer format)")
opt_parser <- add_argument(opt_parser, arg = "--output-file", type = "character", default = "./output/chain.dat",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)


# Preliminary checks ------------------------------------------------------

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

# Check if data file exists
data_file <- file.path(getwd(), extra_args$input_file)
if(!file.exists(data_file)){
  stop(sprintf("%s does not exist", data_file))
}
data_file <- normalizePath(data_file)
cat(sprintf("Data file: %s\n", data_file)) # Log

# Check if Adj. Matrix file exists
adj_file <- file.path(getwd(), dirname(extra_args$input_file), "adj_matrix.dat")
if(!file.exists(adj_file)){
  stop(sprintf("%s does not exist", adj_file))
}
adj_file <- normalizePath(adj_file)
cat(sprintf("Adj. Matrix file: %s\n", adj_file)) # Log

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log


# Main code ---------------------------------------------------------------

# Import required packages
suppressMessages(library("SPMIX"))

# Load input files
load(data_file)
load(adj_file)

# Setting MCMC parameters
burnin <- 30000
niter <- 35000
thin <- 1

# Set sampler parameters template
params_template =
  "
  num_components { %s }

  p0_params {
    mu0: 10
    a: 4
    b: 4
    lam_: 0.1
  }

  rho { %s }

  sigma { %s }

  graph_params { %s }
  "

# Set sampler parameter
params <- sprintf(params_template, extra_args$num_components_prior, extra_args$rho_prior, extra_args$sigma_prior, extra_args$graph_prior)

# Run Spatial sampler
SPMIX_fit <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)
if (exists("SPMIX_fit")) {
  save(SPMIX_fit, file = out_file)
}
