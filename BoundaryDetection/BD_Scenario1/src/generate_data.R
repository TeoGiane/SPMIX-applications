# # ---- BOUNDARY DETECTION SCENARIO 1: GENERATE DATA ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_data", hide.opts = TRUE,
                         description = "Generate N syntethic datasets for the simulation study")
opt_parser <- add_argument(opt_parser, arg = "--num-datasets", type = "integer", default = 50,
                           help = "Number of datasets to generate")
opt_parser <- add_argument(opt_parser, arg = "--dest-dir", type = "character", default = "input",
                           help = "Directory to save generated datasets")
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

# Create input directory if not present
datasets_dir <- file.path(getwd(), extra_args$dest_dir)
if(!dir.exists(datasets_dir)){
  dir.create(datasets_dir, recursive = TRUE)
}
datasets_dir <- normalizePath(datasets_dir)
cat(sprintf("Destination Directory: %s\n", datasets_dir))


# Main code ---------------------------------------------------------------

# Required libraries
suppressMessages(library("SPMIX"))
# suppressMessages(library("sf"))

# Required quantities
numGroups <- 36
dist_picker <- c(rep(cbind(rep(1,3),rep(2,3)),3), rep(cbind(rep(2,3),rep(1,3)),3))

# Set seed
set.seed(230196)

# Generate data
for (sim in 1:extra_args$num_datasets) {
  # Clean data list
  data <- list()
  # Generate data
  for (i in 1:numGroups) {
    if (dist_picker[i] == 1) {
      data[[i]] <- metRology::rt.scaled(n=100, df=6, mean=4, sd=1.5)
    } else {
      data[[i]] <- sn::rsn(n=100, xi=4, omega=1.3, alpha=-3)
      attributes(data[[i]]) <- NULL
    }
  }
  # Save data to file
  if (exists("data")) {
    filename <- file.path(datasets_dir, sprintf("data%03d.dat", sim))
    save(data, file = filename)
  }
}
