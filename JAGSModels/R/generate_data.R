# # ---- CALIFORNIA CENSUS DATASET - GENERATE DATA ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_data", hide.opts = TRUE,
                         description = "Generate N subsampled datasets from the California Census Dataset
                                        and stores them in the data/ directory, along with files required in simulations")
opt_parser <- add_argument(opt_parser, arg = "--num-datasets", type = "integer", default = 10,
                           help = "Number of datasets to generate")
opt_parser <- add_argument(opt_parser, arg = "--subsample-size", type = "integer", default = 100,
                           help = "Number of sub-sampled data in each area")
extra_args <- parse_args(opt_parser)


# Preliminary checks ------------------------------------------------------

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Create data directory if not present
data_dir <- file.path(getwd(), "data")
if(!dir.exists(data_dir)){
  dir.create(data_dir, recursive = TRUE)
}
data_dir <- normalizePath(data_dir)
cat(sprintf("Data Directory: %s\n", data_dir)) # Log


# Main code ---------------------------------------------------------------

# Import packages
suppressMessages(library("sf"))
suppressMessages(library("dplyr"))
suppressMessages(library("spdep"))

# Import usa shapefiles
sf_usa <- read_sf(file.path(getwd(),"raw","us-pumas","us-pumas.shp"))

# Extract LA, Orange & Ventura counties (sort PUMAs in ascending order)
sf_counties <- sf_usa %>%
  filter(State == 'California', grepl("Los Angeles County|Ventura County|Orange County", Name)) %>%
  select(PUMA, Name, geometry) %>%
  arrange(PUMA)

# Compute adjacency matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Import raw data
raw_data <- read.csv(file.path(getwd(),"raw","psam_p06.csv"))

# Clean data
clean_data <- raw_data %>%
  mutate(PUMA = sprintf("%05d", PUMA)) %>%
  filter(PINCP > 0, PUMA %in% sf_counties$PUMA) %>%
  mutate(LPINCP = log(PINCP)) %>%
  select(PUMA, LPINCP) %>%
  group_by(PUMA) %>%
  group_split() %>%
  as.list()

# Check if subsample size is valid
if(extra_args$subsample_size > min(sapply(clean_data, nrow))){
  msg <- sprintf("'--subsample-size' parameter not valid. Set a value <= %d",
                 min(sapply(clean_data, nrow)))
  stop(msg)
}

# Define subsample_data function
subsample_data <- function(df){
  out <- df$LPINCP[sample(1:nrow(df), extra_args$subsample_size)]
}

# Set seed
set.seed(230196)

# Generate data
for (sim in 1:extra_args$num_datasets) {
  # Random subsample (n = 100)
  data <- lapply(clean_data, subsample_data)
  names(data) <- sf_counties$PUMA
  # Save data to file
  if (exists("data")) {
    filename <- file.path(data_dir,sprintf("data_%03d.dat",sim))
    save(data, file = filename)
  }
}

# Save files in the data/ directory
shp_dir <- file.path(data_dir, "counties-pumas")
if(!dir.exists(shp_dir)){
  dir.create(shp_dir, recursive = TRUE)
}
st_write(sf_counties, file.path(shp_dir,"counties-pumas.shp"))
save(W, file = file.path(data_dir,"adj_matrix.dat"))
