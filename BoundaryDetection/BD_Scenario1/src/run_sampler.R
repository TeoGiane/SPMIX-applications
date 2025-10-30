# # ---- BD SCENARIO 1 - RUN SAMPLER ---- # #

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
opt_parser <- add_argument(opt_parser, arg = "--poisson-rate", type = "double", default = NULL,
                           help = "Rate parameter for the shifted Poisson prior on the number of components (used only if --num-components is 'RJ')")
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

# Check if input-file exits
data_file <- file.path(getwd(), extra_args$input_file)
if(!file.exists(data_file)){
  stop(sprintf("%s does not exist", data_file))
}
data_file <- normalizePath(data_file)
cat(sprintf("Data file: %s\n", data_file)) # Log

# Set prior on H according to --num-components
if (extra_args$num_components == "RJ") {
  if (is.null(extra_args$poisson_rate)) {
    stop("Please provide a value for --poisson-rate when using RJMCMC")
  }
  H <- sprintf("{ shifted_poisson_prior { rate: %g } }", extra_args$poisson_rate)
} else {
  H <- sprintf("{ fixed: %g }", as.integer(extra_args$num_components))
}
cat(sprintf("N° of Components: %s\n", H)) # Log
cat(sprintf("rho: %g\n", extra_args$rho)) # Log

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log


# Main code ---------------------------------------------------------------

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
  num_components %s

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
      alpha: 4
      beta: 4
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
SPMIX_fit <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)
if (exists("SPMIX_fit")) {
  save(SPMIX_fit, file = out_file)
}
