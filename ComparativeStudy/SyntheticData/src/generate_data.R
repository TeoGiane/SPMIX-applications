# # ---- Comparative Study on Synthetic Data - Data generation ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_data", hide.opts = TRUE,
                         description = "Generates synthetic data for comparative study")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "input",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)

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
cat("Setting working directory to:", getwd(), "\n")

# Generate output directory if it does not exist
output_dir <- file.path(getwd(), extra_args$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Required libraries
library("sf")

# Generate shape file from scratch
numGroups <- 36
box <- st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Generate distribution picker label
dist_picker <- c(rep(cbind(rep(1, 3),rep(2, 3)), 3), rep(cbind(rep(2, 3),rep(1, 3)), 3))

# Generate sf object
sf_grid <- st_sf(data.frame("Group" = as.factor(dist_picker)), geometry = grid, crs = 4326)

# Aux quantities
mean_g1 <- c(-2, 2); std_dev_g1 <- c(1,1); clus_allocs_g1 <- c(rep(1,50), rep(2,50))
mean_g2 <- 0; std_dev_g2 <- sqrt(5)

# Generate data
set.seed(230196); data <- list()
for (i in 1:numGroups) {
  if (dist_picker[i] == 1) {
    data[[i]] <- rnorm(200, mean_g1[clus_allocs_g1], std_dev_g1[clus_allocs_g1])
  } else {
    data[[i]] <- rnorm(200, mean_g2, std_dev_g2)
    attributes(data[[i]]) <- NULL
  }
}

# Create shapefile directory
shapefile_dir <- file.path(extra_args$output_dir, "shp")
if (!dir.exists(shapefile_dir)) {
  dir.create(shapefile_dir, recursive = TRUE)
  cat("Created shapefile directory:", shapefile_dir, "\n")
}

# Save shapefile
write_sf(sf_grid, file.path(shapefile_dir, "grid.shp"))
cat("Shapefile saved to:", file.path(shapefile_dir, "grid.shp"), "\n")

# Save data
save(data, file = file.path(extra_args$output_dir, "data.dat"))
cat("Data saved to:", file.path(extra_args$output_dir, "data.dat"), "\n")
