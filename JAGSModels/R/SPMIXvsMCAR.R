# Extra option parser -----------------------------------------------------

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "SPMIXvsMCAR", hide.opts = TRUE,
                         description = "Compares SPMIX::Sampler.BoundaryDetection and MCAR model on summary statistics on a simulated scenario")
opt_parser <- add_argument(opt_parser, arg = "--a-beta", type = "double", default = 2,
                           help = "Value of 'a-beta' parameter, fixed in (0,1)")
opt_parser <- add_argument(opt_parser, arg = "--b-beta", type = "double", default = 36,
                           help = "Value of 'b-beta' parameter, fixed in (0,1)")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "output",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)


# Preliminary Checks ------------------------------------------------------

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Create directory for output if does not exist
out_dir <- file.path(getwd(), extra_args$output_dir)
if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
cat(sprintf("Output Directory: %s\n", normalizePath(out_dir))) # Log


# Main Code ---------------------------------------------------------------

# Required libraries
suppressMessages(library("rjags"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))
suppressMessages(library("SPMIX"))

# Data generation ---------------------------------------------------------

# Generate shape file from scratch
numGroups <- 36
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Generate distribution picker label
dist_picker <- c(rep(cbind(rep(1,3),rep(2,3)),3), rep(cbind(rep(2,3),rep(1,3)),3))

# Generate sf object
sf_grid <- st_sf(data.frame("Group" = as.factor(dist_picker)), geometry = grid)

# Aux quantities
mean_g1 <- c(-2,2); std_dev_g1 <- c(1,1); clus_allocs_g1 <- c(rep(1,50), rep(2,50))
mean_g2 <- 0; std_dev_g2 <- sqrt(5)

# Generate data
set.seed(230196); data <- list()
for (i in 1:numGroups) {
  if (dist_picker[i] == 1) {
    data[[i]] <- rnorm(200, mean_g1[clus_allocs_g1], std_dev_g1[clus_allocs_g1])
  } else {
    data[[i]] <- rnorm(200, mean_g2, std_dev_g2)
    attributes(data[[i]]) <- NULL
  }
}

# Compute summary statistics for MCAR model
Y = c(sapply(data, quantile, 0.05),
      sapply(data, quantile, 0.25),
      sapply(data, quantile, 0.5),
      sapply(data, quantile, 0.75),
      sapply(data, quantile, 0.95)); names(Y) = NULL


# SPMIX sampler -----------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 1

# Set sampler parameters
params_template =
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

# Run Sampler and store result

out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params, type = "rjmcmc")
if (exists("out")) {
  filename <- sprintf("%s/SPMIXvsMCAR-SPMIX-chain-%s.dat", out_dir, format(Sys.time(), format = "%Y%m%d-%H%M"))
  save(out, file = filename)
  cat(sprintf("Saved SPMIX chain at: %s\n", filename))
}


# # MCAR sampler ------------------------------------------------------------

# Prepare data for JAGS
MCAR_data <- list("Y"=Y,
                  "I"=numGroups,
                  "rho"=0.95,
                  "W"=W,
                  "Ncov"=5,
                  "mu0"=mean(Y),
                  "abeta"=extra_args$a_beta,
                  "bbeta"=extra_args$b_beta)

# Compile model
MCAR_model <- jags.model("bug/multivariate_CAR_marginalized.bug", data=MCAR_data, n.chains=1, n.adapt = 0)

# Run MCMC
MCAR_params <- c("G", "p", "one_over_sigmasq", "one_over_tausq")
MCAR_fit <- jags.samples(MCAR_model, variable.names = MCAR_params, n.iter = 1000, thin = 1)

# Store MCAR_fit
filename <- sprintf("%s/SPMIXvsMCAR-MCAR-chain-%s.dat", out_dir, format(Sys.time(), format = "%Y%m%d-%H%M"))
save(MCAR_fit, file = filename)

