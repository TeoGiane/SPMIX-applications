###########################################################################
# Command line input options via argparser --------------------------------

# Parse argument from terminal
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compute_global_comparison", hide.opts = TRUE,
                         description = "Compute global comparison of L1 distances from a posterior MCMC chain for the California Census Dataset")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c", default = "RJ",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--rho", short = "-r", default = "0.95",
                           help = "Value for rho parameter")
opt_parser <- add_argument(opt_parser, arg = "--output-folder", short = "-o", default = "summary/L1_Glob",
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
sf_file <- file.path(getwd(), "data", "counties-pumas", "counties-pumas.shp")
if(!file.exists(sf_file)){stop(sprintf("File: '%s' does not exists", sf_file))}
sf_file <- normalizePath(sf_file)
cat(sprintf("File '%s' exists\n", sf_file)) # Log

# Check if data file exists
data_file <- file.path(getwd(), "data", sprintf("data_%03d.dat",num_sample))
if(!file.exists(data_file)){stop(sprintf("File: '%s' does not exists", data_file))}
data_file <- normalizePath(data_file)
cat(sprintf("File '%s' exists\n", data_file)) # Log

# Check if plinks file exists
plinks_file <- file.path(getwd(), "summary", "plinks", sprintf("plinks_%03d.csv",num_sample))
if(!file.exists(plinks_file)){stop(sprintf("File: '%s' does not exists", plinks_file))}
plinks_file <- normalizePath(plinks_file)
cat(sprintf("File '%s' exists\n", plinks_file)) # Log

# Check if estimated densitites file exists
mean_est_dens_file <- file.path(getwd(), "summary", "mean_est_dens", sprintf("mean_est_dens_%03d.csv",num_sample))
if(!file.exists(mean_est_dens_file)){stop(sprintf("File: '%s' does not exists", mean_est_dens_file))}
mean_est_dens_file <- normalizePath(mean_est_dens_file)
cat(sprintf("File '%s' exists\n", mean_est_dens_file)) # Log
 
# Create output directory if does not exist
out_folder <- file.path(getwd(), extra_args$output_folder)
if(!dir.exists(out_folder)){dir.create(out_folder, recursive = TRUE)}
out_folder <- normalizePath(out_folder)
cat(sprintf("Output Folder: '%s'\n", out_folder)) # Log

cat("\n") # Blank line

###########################################################################

###########################################################################
# Main code ---------------------------------------------------------------

# Required libraries
suppressMessages(library("SPMIX"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))

# FUNCTION - Compute L1 distance
L1_distance <- function(y1, y2, x) {
  I <- diff(range(x))
  return(I * mean(abs(y1 - y2)))
}

# Load Data File
load(data_file)
cat(sprintf("Loaded: %s\n", data_file)) # Log

# Load shape file
sf_counties <- read_sf(sf_file)
cat(sprintf("Loaded: %s\n", sf_file)) # Log

# Load plinks file
plinks <- as.matrix(read.csv(plinks_file, header = F))
cat(sprintf("Loaded: %s\n", plinks_file)) # Log

# Load estimated densities file
mean_est_dens <- as.matrix(read.csv(mean_est_dens_file, header = F))
cat(sprintf("Loaded: %s\n", mean_est_dens_file)) # Log

# Compute common grid
x <- seq(range(data)[1], range(data)[2], length.out = 500)

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Compute admissible edges
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute Neighbouring pairs
Gn <- matrix(NA, nrow(plinks), ncol(plinks))
Gn[Eadj] <- ifelse(plinks[Eadj] >= 0.5, 1, NA)
Gn[lower.tri(Gn)] <- NA
neigh_pairs <- which(Gn == 1, arr.ind = T)

# Compute Boundary pairs
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)
Gb[lower.tri(Gb)] <- NA
bound_pairs <- which(Gb == 1, arr.ind = T)

# Compute L1 distance between all neighbouring pairs (GLOBAL)
L1_neigh <- data.frame("Type" = rep("Neigh", nrow(neigh_pairs)), "Dist" = rep(NA, nrow(neigh_pairs)))
for (i in 1:nrow(neigh_pairs)) {
  L1_neigh$Dist[i] <- L1_distance(mean_est_dens[neigh_pairs[i,1],], mean_est_dens[neigh_pairs[i,2],], x)
}

# Compute L1 distance between all boundary pairs (GLOBAL)
L1_bound <- data.frame("Type" = rep("Bound", nrow(bound_pairs)), "Dist" = rep(NA, nrow(bound_pairs)))
for (i in 1:nrow(bound_pairs)) {
  L1_bound$Dist[i] <- L1_distance(mean_est_dens[bound_pairs[i,1],], mean_est_dens[bound_pairs[i,2],], x)
}

# Bind by rows
L1_Glob <- rbind(L1_neigh, L1_bound); L1_Glob$Dataset <- rep(num_sample, nrow(L1_Glob))

# Save plinks to file
output_file <- file.path(out_folder, sprintf("L1_Glob-%03d.csv",num_sample))
write.table(L1_Glob, file = output_file, sep=",", col.names = T, row.names = F)
cat(sprintf("Successfully write L1_Glob to '%s'\n", output_file)) # Log

###########################################################################
