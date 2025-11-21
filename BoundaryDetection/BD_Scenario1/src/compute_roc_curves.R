# # ---- BOUNDARY DETECTION SCENARIO 1: COMPUTE ROC CURVES ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "compute_roc_curves", hide.opts = TRUE,
                         description = "For each simulation in a given scenario, compute the ROC curve for all possible boundary detection thresholds")
opt_parser <- add_argument(opt_parser, arg = "--num-datasets", type = "integer", short = "-d",
                           help = "Number of datasets to consider in the given scenario")
opt_parser <- add_argument(opt_parser, arg = "--num-components", short = "-c",
                           help = "Value for the number of components, or 'RJ' if the reverisble jump sampler is considered")
opt_parser <- add_argument(opt_parser, arg = "--poisson-rate", default = NULL,
                           help = "Rate parameter for the shifted Poisson prior on the number of components (used only if --num-components is 'RJ')")
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
if(extra_args$num_components == "RJ" && is.null(extra_args$poisson_rate)){
  stop("Please provide a value for --poisson-rate when using RJMCMC")
}
if(is.na(extra_args$num_components)) {stop("Input parameter '--num-components' not specified")}
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
suppressMessages(library("pROC"))
suppressMessages(library("parallel"))

# Generate shapefile
numGroups <- 36
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Compute adjacency matrix and related quantities
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")
Eadj <- which(W != 0, arr.ind = T)
Ena <- which(W == 0, arr.ind = T)

# Load true boundary graph
Gb_true <- ReadMatrixFromCSV("truth/true_graph.csv")
Gb_true[Ena] <- NA

compute_roc_object <- function(true_graph, plinks) {
  # Compute posterior probability of a boundary in vector form
  predicted_boundaries <- na.omit(1.0 - plinks[upper.tri(plinks)])
  # Compute actual boundaries in vector form
  actual_boundaries <- na.omit(true_graph[upper.tri(true_graph)])
  return(roc(actual_boundaries, predicted_boundaries))
}

# Function to process a single dataset
process_dataset <- function(id) {
  # Load data from file
  data_file <- file.path(data_folder, sprintf("data%03d.dat", id))
  load(data_file)
  # Load chain from file
  chain_file <- file.path(chains_folder, sprintf("chain%03d.dat", id))
  load(chain_file)
  # Deserialize chain and compute plinks matrix
  chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("spmix.UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)
  plinks[Ena] <- NA
  # Compute ROC values and store in proper dataset structure
  roc_obj <- compute_roc_object(Gb_true, plinks)
  roc_df <- coords(roc_obj, x = seq(0, 1, by = 0.01), ret = c("threshold", "tn", "tp", "fn", "fp"))
  roc_df$auc <- rep(roc_obj$auc, nrow(roc_df))
  roc_df$dataset_id <- rep(id, nrow(roc_df))
  return(roc_df)
}

# Process datasets in parallel
num_cores <- min(num_datasets, detectCores() - 1)
cat("Processing datasets in parallel using", num_cores, "cores... ") # Log
cl <- makeCluster(num_cores)
clusterExport(cl, c("data_folder", "chains_folder", "DeserializeSPMIXProto", "compute_roc_object", "roc", "coords", "Gb_true", "Ena"))
results <- parLapply(cl, 1:num_datasets, process_dataset)
stopCluster(cl)
cat("Done!\n") # Log

# Combine all results into a single data frame
df <- do.call(rbind, results)

# Write csv to file
write.table(df, file = out_file, sep = ",", row.names = FALSE)
cat(sprintf("ROC curves written to file: %s\n", out_file)) # Log
