# # ---- BD SCENARIO 1 - COMPUTE CONFUSION MATRICES ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compare_boundary_graph", hide.opts = TRUE,
                         description = "For each simulation in a given scenario, compare the estimated and true graph via Standardized Hamming Distance")
opt_parser <- add_argument(opt_parser, arg = "--num-datasets", type = "integer",
                           help = "Number of datasets to consider in the given scenario")
opt_parser <- add_argument(opt_parser, arg = "--num-components",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--poisson-rate", default = NULL,
                           help = "Rate parameter for the shifted Poisson prior on the number of components (used only if --num-components is 'RJ')")
opt_parser <- add_argument(opt_parser, arg = "--rho",
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
if(extra_args$num_components == "RJ" && is.null(extra_args$poisson_rate)){
  stop("Please provide a value for --poisson-rate when using RJMCMC")
}
if(is.na(extra_args$rho)){stop("Input parameter '--rho' not specified")}

# Set input parameters
num_datasets <- as.integer(extra_args$num_datasets)
rho <- extra_args$rho
H <- extra_args$num_components
poisson_rate <- extra_args$poisson_rate

# Deduce input folder
data_folder <- file.path(getwd(), "input")
if(!dir.exists(data_folder)){
  stop(sprintf("'%s' does not exists", data_folder))
}
cat(sprintf("Data Folder: %s\n", data_folder)) # Log

# Deduce chains folder
if(H == "RJ") { H <- sprintf("RJ/poisson%s", extra_args$poisson_rate) }
chains_folder <- file.path(getwd(), "output", sprintf("H%s",H), sprintf("rho%s",rho))
if(!dir.exists(chains_folder)){s
  top(sprintf("'%s' does not exists", chains_folder))
}
cat(sprintf("Chains Folder: %s\n", chains_folder)) # Log

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
suppressMessages(library("parallel"))

# Function that computes confusion matrix
confusion_df <- function(G_est, G_true) {
  # Check data type
  if(!is.matrix(G_est)) {stop("'G_est' must be a matrix")}
  if(!is.matrix(G_true)) {stop("'G_true' must be a matrix")}
  # Check dimension
  if((nrow(G_est) != nrow(G_true)) || (ncol(G_est) != ncol(G_true))) {
    stop("Dimension mismatch between 'G_est' and 'G_true'")
  }
  # Take only upper diagonal part (since simmetric)
  G_est <- factor(G_est[upper.tri(G_est, diag = F)], levels = c(0,1))
  G_true <- factor(G_true[upper.tri(G_true, diag = F)], levels = c(0,1))
  # Compute confusion matrix
  out_table <- table(G_true, G_est)
  # Store confusion matrix in vector form
  out <- c("TP"=out_table[2,2], "TN"=out_table[1,1],
           "FP"=out_table[1,2], "FN"=out_table[2,1])
}

# Function to process a single dataset
process_dataset <- function(id) {
  # Load data from file
  data_file <- file.path(data_folder, sprintf("data%03d.dat", id))
  load(data_file)
  # Load chain from file
  chain_file <- file.path(chains_folder, sprintf("chain%03d.dat", id))
  load(chain_file)
  # Deserialize chain
  chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("spmix.UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  # Compute point estimate of posterior boundary graph
  plinks <- Reduce('+', G_chain)/length(G_chain)
  Gb_est <- matrix(NA, nrow(plinks), ncol(plinks))
  Gb_est[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, 0)
  # Compute standardized hamming distance
  return(confusion_df(Gb_est, Gb_true))
}

# Generate shapefile
numGroups <- 36
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Compute adjacency matrix and related quantities
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")
Eadj <- which(W != 0, arr.ind = T)
Ena <- which(W == 0, arr.ind = T)

# Load true graph
Gb_true <- ReadMatrixFromCSV("truth/true_graph.csv")
Gb_true[Ena] <- NA

# Process datasets in parallel
num_cores <- min(num_datasets, detectCores() - 1)
cat("Processing datasets in parallel using", num_cores, "cores... ") # Log
cl <- makeCluster(num_cores)
clusterExport(cl, c("data_folder", "chains_folder", "Gb_true", "Eadj", "confusion_df", "DeserializeSPMIXProto"))
results <- parSapply(cl, 1:num_datasets, process_dataset)
stopCluster(cl)
cat("Done!\n") # Log

# Convert results to data frame (NON CREDO SIA CORRETTO QUI)
df <- as.data.frame(t(results))

# Write csv to file
write.table(df, file=out_file, sep=",", row.names=F)
cat(sprintf("Confusion matrices written to: %s\n", normalizePath(out_file))) # Log
