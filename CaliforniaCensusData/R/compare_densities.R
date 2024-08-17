# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compare_densities", hide.opts = TRUE,
                         description = "Compare the estimated densities via L1 Distance for the California Census Dataset")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "output-file",
                           help = "Relative path to the output .csv file")
extra_args <- parse_args(opt_parser)

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check input parameters
if(is.na(extra_args$num_components)) {stop("Input parameter '--num-components' not specified")}

# Set input parameters
H <- extra_args$num_components

# Check if deduced folders exist
# Data
data_folder <- file.path(getwd(), "data")
if(!dir.exists(data_folder)){stop(sprintf("'%s' does not exists", data_folder))}
cat(sprintf("Data Folder: %s\n", data_folder)) # Log
# Chains Parent folder
chains_folder <- file.path(getwd(), "output", sprintf("H_%s", H))
if(!dir.exists(chains_folder)){stop(sprintf("'%s' does not exists", chains_folder))}
cat(sprintf("Chains Parent Folder: %s\n", chains_folder)) # Log
# Get all subfolders
chains_folder <- list.dirs(chains_folder, recursive = FALSE)
# Reference chain folder
ref_chain_folder <- file.path(getwd(), "output", "H_RJ", "rho_0.95")
if(!dir.exists(ref_chain_folder)){stop(sprintf("'%s' does not exists", ref_chain_folder))}
cat(sprintf("Reference Chain Folder: %s\n", ref_chain_folder)) # Log

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

# Function that computes the L1 distance
L1_distance <- function(y1, y2, x){
  I <- diff(range(x))
  return(I * mean(abs(y1 - y2)))
}

# Load true densities and common x_grid
load("data/data_001.dat")
x_grid <- seq(range(data)[1], range(data)[2], length.out = 500)

# Define numGroups
numGroups <- length(data)

# Load reference chain from file
ref_chain_file <- file.path(ref_chain_folder, "chain_001.dat")
load(ref_chain_file)
# Deserialize chain
ref_chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
# Compute point estimate of posterior densities in each area
ref_est_dens <- ComputeDensities(ref_chains, x_grid, verbose = F)


# Compute mean L1 distances
df <- data.frame("MeanL1" = numeric(0))
for (curr_chain_folder in chains_folder) {
  # Load current chain from file
  load(file.path(curr_chain_folder, "chain.dat"))
  # Deserialize chain
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  # Compute point estimate of posterior densities in each area
  est_dens <- ComputeDensities(chains, x_grid, verbose = F)
  # Compute mean L1 distance
  curr_meanL1 <- data.frame("MeanL1" = mean(sapply(1:numGroups, function(a){L1_distance(colMeans(est_dens[[a]]), colMeans(ref_est_dens[[a]]), x_grid)})) )
  df <- rbind(df, curr_meanL1)
  cat(sprintf("Processed folder: %s\n", curr_chain_folder))
}

# Write csv to file
write.table(df, file=out_file, sep=",", row.names=F)
