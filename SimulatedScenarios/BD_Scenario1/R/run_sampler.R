# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_sampler", hide.opts = TRUE,
                         description = "Run SPMIX::Sampler.BoundaryDetection on the data stored in 'file'")
opt_parser <- add_argument(opt_parser, arg = "input-file", type = "character",
                           help = "Relative path to input data file")
opt_parser <- add_argument(opt_parser, arg = "--rho", type = "double", default = 0.95,
                           help = "Value of 'rho' parameter, fixed in (0,1)")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--output-file", type = "character", default = "./output/chain.dat",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check if input-file exits
data_file <- file.path(getwd(), extra_args$input_file)
if(!file.exists(data_file)){
  stop(sprintf("%s does not exist", data_file))
}
data_file <- normalizePath(data_file)
cat(sprintf("Data file: %s\n", data_file)) # Log

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


# Required libraries
suppressMessages(library("SPMIX"))
suppressMessages(library("sf"))

# Load data
load(data_file)

# Generate shapefile from scratch
numGroups <- 36
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 1

# Set sampler parameters template
params_template =
  "
  num_components: %g

  p0_params {
    mu0: 0
    a: 2
    b: 2
    lam_: 0.1
  }

  rho {
    fixed: %g
  }

  sigma {
    inv_gamma_prior {
      alpha: 6
      beta: 6
    }
  }

  graph_params {
    beta_prior {
      a: 2
      b: 36
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
