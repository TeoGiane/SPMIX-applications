###########################################################################
# Command line input options via argparser --------------------------------

# Parse argument from terminal
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compute_estimated_densities", hide.opts = TRUE,
                         description = "Compute posterior estimated mean densities from a posterior MCMC chain for the California Census Dataset")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c", default = "RJ",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--rho", short = "-r", default = "0.95",
                           help = "Value for rho parameter")
opt_parser <- add_argument(opt_parser, arg = "--output-folder", short = "-o", default = "summary/mean_est_dens",
                           help = "Relative path to the output folder")
opt_parser <- add_argument(opt_parser, arg = "num-sample",
                           help = "The number of subsample to consider")
extra_args <- parse_args(opt_parser)

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Set input parameters
H <- extra_args$num_components
rho <- extra_args$rho
num_sample <- as.integer(extra_args$num_sample)

# Check if data file exists
data_file <- file.path(getwd(), "data", sprintf("data_%03d.dat",num_sample))
if(!file.exists(data_file)){stop(sprintf("'%s' does not exists", data_file))}
data_file <- normalizePath(data_file)
cat(sprintf("Data File: %s\n", data_file)) # Log

# Check if chain file exists
chain_file <- file.path(getwd(), "output", sprintf("H_%s",H), sprintf("rho_%s",rho), sprintf("chain_%03d.dat",num_sample))
if(!file.exists(chain_file)){stop(sprintf("'%s' does not exists", chain_file))}
chain_file <- normalizePath(chain_file)
cat(sprintf("Chain File: %s\n", chain_file)) # Log

# Create output directory if does not exist
out_folder <- file.path(getwd(), extra_args$output_folder)
if(!dir.exists(out_folder)){dir.create(out_folder, recursive = TRUE)}
out_folder <- normalizePath(out_folder)
cat(sprintf("Output Folder: %s\n", out_folder)) # Log

cat("\n") # Blank line

###########################################################################

###########################################################################
# Main code ---------------------------------------------------------------

# Required libraries
suppressMessages(library("SPMIX"))
suppressMessages(library("sf"))

# Load Data File
load(data_file)
cat(sprintf("Loaded data: %s\n", data_file)) # Log

# Load and deserialse chain File
load(chain_file)
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
cat(sprintf("Loaded chain: %s\n", chain_file)) # Log

# Compute estimated density
x <- seq(range(data)[1], range(data)[2], length.out = 500)
estimated_densities <- ComputeDensities(chains, x, verbose = TRUE)
mean_est_dens <- t(sapply(estimated_densities, colMeans))

# Save mean_est_dens to file
output_file <- file.path(out_folder, sprintf("mean_est_dens_%03d.csv",num_sample))
write.table(mean_est_dens, file = output_file, sep=",", col.names = F, row.names = F)
cat(sprintf("Successfully write 'mean_est_dens' to: %s\n", output_file)) # Log

###########################################################################

