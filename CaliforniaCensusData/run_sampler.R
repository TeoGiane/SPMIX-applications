# Import required packages
library("optparse")
suppressMessages(library("SPMIX"))

# Argument parser
option_list <- list (
  make_option(c("--sim_name"), type = "character", default = NULL,
              help = "Name of the simulation to run. If the file 'input/params_<sim_name>.asciipb' exists,
                      the sampler will the prior hyperparameters specified in that file,
                      otherwise default values will be used.",
              metavar = "character")
)
opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

# log
cat(sprintf("Current Directory: %s\n", getwd()))

# Create directory to store files
sim_folder <- sprintf("%s/output/%s", getwd(), args$sim_name)
dir.create(sim_folder, showWarnings = F)

# Select parameter file to use
params_file <- ifelse(file.exists(sprintf("%s/input/params_%s.asciipb", getwd(), args$sim_name)),
                      sprintf("%s/input/params_%s.asciipb", getwd(), args$sim_name),
                      sprintf("%s/input/params_default.asciipb", getwd()))

# log
cat(sprintf("Using prior hyperparameters in file: %s\n", params_file))

# Copy parameter file into output folder
success <- file.copy(params_file, sim_folder, overwrite = TRUE)

#log
if(success) {
  cat(sprintf("Prior hyperparameters file has been copied in: %s\n", sim_folder))
}

# Load data and adjacency matrix
load(sprintf("%s/data/clean_data.dat", getwd()))
load(sprintf("%s/data/adj_matrix.dat", getwd()))

# log
cat("Data have been loaded\n")

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 1

# log
cat("Sampler is about to start:\n\n")

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params_file)

# Save output
if (exists("out")) {
  filename <- sprintf("%s/chain_%s.dat", sim_folder, args$sim_name)
  save(out, file = filename)
}