# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compare_densities", hide.opts = TRUE,
                         description = "For each simulation in a given scenario, compare the estimated and true densities via L1 Distance")
opt_parser <- add_argument(opt_parser, arg = "--num-datasets", type = "integer", short = "-d",
                           help = "Number of datasets to consider in the given scenario")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--rho", short = "-r",
                           help = "Value of 'rho' parameter, fixed in (0,1)")
opt_parser <- add_argument(opt_parser, arg = "output-file",
                           help = "Relative path to the output .csv file")
extra_args <- parse_args(opt_parser)

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(basedir)
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check input parameters
if(is.na(extra_args$num_datasets)) {stop("Input parameter '--num-datasets' not specified")}
if(is.na(extra_args$num_components)) {stop("Input parameter '--num-components' not specified")}
if(is.na(extra_args$rho)){stop("Input parameter '--rho' not specified")}

# Set input parameters
num_datasets <- as.integer(extra_args$num_datasets)
rho <- extra_args$rho
H <- extra_args$num_components

# Check if deduced folders exist
data_folder <- file.path(getwd(), "input")
if(!dir.exists(data_folder)){stop(sprintf("'%s' does not exists", data_folder))}
cat(sprintf("Data Folder: %s\n", data_folder)) # Log
chain_folder <- file.path(getwd(), "output", sprintf("H_%s",H), sprintf("rho_%s",rho))
if(!dir.exists(chain_folder)){stop(sprintf("'%s' does not exists", chain_folder))}
cat(sprintf("Chain Folder: %s\n", chain_folder)) # Log


if(length(list.files(chain_folder)) < num_datasets){
  stop("Number of available chain files are less that required")
}

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log

# --- End of checks --- #

# Load required packages
suppressMessages(library("SPMIX"))
suppressMessages(library("sf"))

# Function that computes the WAIC distance
WAIC_distance <- function(true, est){
  if(is.null(nrow(true))){
    # Replicate true to estimate waic using loo package
    true <- t(replicate(nrow(est), true))
  }
  # Diff between estimated and true WAIC
  return(loo::waic(log(abs(est-true)))$estimates["waic", "Estimate"])
}

# Define numGroups
numGroups <- 36

# Load true densities and common x_grid
x_grid <- as.numeric(ReadMatrixFromCSV("truth/common_grid.csv"))
true_dens <- ReadMatrixFromCSV("truth/true_densities.csv")

# Compute mean L1 distances for each simulated dataset
df <- data.frame("MeanWAIC"=rep(NA,num_datasets))
for (n in 1:num_datasets) {
  # Load data from file
  data_file <- file.path(data_folder, sprintf("data_%03d.dat", n))
  load(data_file)
  # Load chain from file
  chain_file <- file.path(chain_folder, sprintf("chain_%03d.dat", n))
  load(chain_file)
  # Deserialize chain
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  # Compute point estimate of posterior densities in each area
  est_dens <- ComputeDensities(chains, x_grid, verbose = T)
  # Compute mean L1 distance
  df[n,] <- suppressWarnings(mean(sapply(1:numGroups, function(a){WAIC_distance(true_dens[a,], est_dens[[a]])})))
  cat(sprintf("\rProcessed file: %d / %d", n, num_datasets))
}
cat("\n")

# Write csv to file
write.table(df, file=out_file, sep=",", row.names=F)
