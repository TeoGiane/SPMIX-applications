###########################################################################
# Command line input options via argparser --------------------------------

# Parse argument from terminal
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compute_plinks_matrix", hide.opts = TRUE,
                         description = "Compute posterior probabiliy of inclusion matrix from a posterior MCMC chain for the California Census Dataset")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c", default = "RJ",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--rho", short = "-r", default = "0.95",
                           help = "Value for rho parameter")
opt_parser <- add_argument(opt_parser, arg = "--output-folder", short = "-o", default = "summary/plinks",
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

# Check if shape file exists
sf_file <- file.path(getwd(),"data","counties-pumas","counties-pumas.shp")
if(!file.exists(sf_file)){stop(sprintf("'%s' does not exists", sf_file))}
sf_file <- normalizePath(sf_file)
cat(sprintf("Shape File: %s\n", sf_file)) # Log

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
suppressMessages(library("spdep"))

# Load Data File
sf_counties <- read_sf(sf_file)
cat(sprintf("Loaded shapefile: %s\n", sf_file)) # Log

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Load and deserialize chain file
load(chain_file)
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
cat(sprintf("Loaded chain: %s\n", chain_file)) # Log

# Compute plinks matrix
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Save plinks to file
output_file <- file.path(out_folder, sprintf("plinks_%03d.csv",num_sample))
write.table(plinks, file = output_file, sep=",", col.names = F, row.names = F)
cat(sprintf("Successfully write 'plinks' to: %s\n", output_file)) # Log

###########################################################################
