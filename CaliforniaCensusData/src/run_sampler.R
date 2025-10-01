# # ---- CALIFORNIA CENSUS DATA - RUN SAMPLER ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_sampler", hide.opts = TRUE,
                         description = "Run SPMIX::Sampler.BoundaryDetection on the California Census dataset")
opt_parser <- add_argument(opt_parser, arg = "input-file", type = "character",
                           help = "Relative path to input data file")
opt_parser <- add_argument(opt_parser, arg = "--rho", type = "double", default = 0.95,
                           help = "Value of 'rho' parameter, fixed in (0,1)")
opt_parser <- add_argument(opt_parser, arg = "--num-components", type = "character",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--output-file", type = "character", default = "./output/chain.dat",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)


# Preliminary checks ------------------------------------------------------

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

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

# Set algo type according to --num-components
if (extra_args$num_components == "RJ") {
  algo_type <- "rjmcmc"
  H <- 10L
} else {
  algo_type <- "no_rjmcmc"
  H <- as.integer(extra_args$num_components)
}
cat(sprintf("Algorithm Type: %s\n", algo_type)) # Log
cat(sprintf("N° of Components: %g\n", H)) # Log
cat(sprintf("rho: %g\n", extra_args$rho)) # Log

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
burnin <- 20000
niter <- 40000
thin <- 1

# Set sampler parameters template
params_template =
  "
  num_components: %g

  p0_params {
    mu0: 10
    a: 4
    b: 4
    lam_: 0.1
  }

  rho {
    fixed: %g
  }

  sigma {
    inv_gamma_prior {
      alpha: 4
      beta: 4
    }
  }

  graph_params {
    beta_prior {
      a: 2
      b: 93
    }
  }
  "

# Set sampler parameter
params <- sprintf(params_template, H, extra_args$rho)

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params, type = algo_type)
if (exists("out")) {
  save(out, file = out_file)
}
