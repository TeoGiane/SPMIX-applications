# # ---- BOUNDARY DETECTION SCENARIO 1: COMPUTE WAIC ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compute_WAIC", hide.opts = TRUE,
                         description = "For each simulation in a given scenario, compute the WAIC for the estimated densities")
opt_parser <- add_argument(opt_parser, arg = "--num-datasets", type = "integer", short = "-d",
                           help = "Number of datasets to consider in the given scenario")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--rho", short = "-r",
                           help = "Value of 'rho' parameter, fixed in (0,1)")
opt_parser <- add_argument(opt_parser, arg = "output-file",
                           help = "Relative path to the output .csv file")
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

# Check input parameters
if(is.na(extra_args$num_datasets)) {stop("Input parameter '--num-datasets' not specified")}
if(is.na(extra_args$num_components)) {stop("Input parameter '--num-components' not specified")}
if(is.na(extra_args$rho)){stop("Input parameter '--rho' not specified")}

# Set input parameters
num_datasets <- as.integer(extra_args$num_datasets)
rho <- extra_args$rho
H <- extra_args$num_components

# Deduce input folder
data_folder <- file.path(getwd(), "input")
if(!dir.exists(data_folder)){
  stop(sprintf("'%s' does not exists", data_folder))
}
cat(sprintf("Data Folder: %s\n", data_folder)) # Log

# Deduce chains folder
chains_folder <- file.path(getwd(), "output", sprintf("H%s",H), sprintf("rho%s",rho))
if(!dir.exists(chains_folder)){
  stop(sprintf("'%s' does not exists", chains_folder))
}
cat(sprintf("Chain Folder: %s\n", chains_folder)) # Log

# Check if too many datasets were required
if(length(list.files(chains_folder)) < num_datasets){
  stop("Number of available chain files are less that required")
}

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log


# Main code ---------------------------------------------------------------

# Load required packages
suppressMessages(library("SPMIX"))
suppressMessages(library("sf"))
suppressMessages(library("loo"))
suppressMessages(library("parallel"))

# Function to process a single dataset
process_dataset <- function(id) {
  # Load data from file
  data_file <- file.path(data_folder, sprintf("data_%03d.dat", id))
  load(data_file)
  # Load chain from file
  chain_file <- file.path(chains_folder, sprintf("chain_%03d.dat", id))
  load(chain_file)
  # Deserialize chain
  chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("UnivariateState",x))
  # Compute point estimate of posterior densities in each area
  post_lpdf <- ComputePosteriorLPDF(data, chains, verbose = F)
  # Compute WAIC
  suppressWarnings(waic(post_lpdf)$estimates["waic", "Estimate"])
}

# Process datasets in parallel
num_cores <- min(num_datasets, detectCores() - 1)
cat("Processing datasets in parallel using", num_cores, "cores... ") # Log
cl <- makeCluster(num_cores)
clusterExport(cl, c("data_folder", "chains_folder", "DeserializeSPMIXProto", "ComputePosteriorLPDF", "waic"))
results <- parSapply(cl, 1:num_datasets, process_dataset)
df <- data.frame("WAIC" = as.numeric(results))
stopCluster(cl)
cat("Done!\n") # Log

# Write csv to file
write.table(df, file=out_file, sep=",", row.names=F)
