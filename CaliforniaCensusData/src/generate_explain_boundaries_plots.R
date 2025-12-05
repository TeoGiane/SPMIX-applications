# # ---- CALIFORNIA CENSUS DATA - GENERATE EXPLAIN BOUNDARIES PLOTS ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_explain_boundaries_plots", hide.opts = TRUE,
                         description = "Generate explain boundaries plots for California Census Dataset.")
opt_parser <- add_argument(opt_parser, arg = "--input-dir", type = "character", default = NULL,
                           help = "Relative path to the input directory.")
opt_parser <- add_argument(opt_parser, arg = "--sim-file", type = "character", default = NULL,
                           help = "Relative path to the simulation file to use.")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "plots",
                           help = "Relative path to save the generated plots")
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
cat("Setting working directory to: ", getwd(), "\n")

# Check if input directory exists
input_dir <- file.path(getwd(), extra_args$input_dir)
if (!dir.exists(input_dir)) {
  stop(sprintf("Input directory '%s' does not exists.", input_dir))
}

# Check if input/explain_boundaries exists
explain_boundaries_dir <- file.path(input_dir, "explain_boundaries")
if (!dir.exists(explain_boundaries_dir)) {
  stop(sprintf("Input explain boundaries directory '%s' does not exists.", explain_boundaries_dir))
}

# Check if simulation requested exists
sim_file <- file.path(getwd(), extra_args$sim_file)
if(!file.exists(sim_file)){
  stop(sprintf("Simulation '%s' does not exists.", sim_file))
}

# Generate output directory if it does not exist
output_dir <- file.path(getwd(), extra_args$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory: ", output_dir, "\n")
}


# Main code ---------------------------------------------------------------

# Load required libraries
suppressMessages(library("dplyr"))
suppressMessages(library("ggmap"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))
suppressMessages(library("SPMIX"))

# Register Stadia Map API key
register_stadiamaps(key = "1f3e613a-2adc-4c22-8bc9-901f1c33e05b")
if(!has_stadiamaps_key()){
  stop("Stadia Maps API key is not set")
}

# Functions
sf_ggmap <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, and convert back to a bbox
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  
  # Return
  return(map)
}

# Import counties shapefile
sf_path <- file.path(input_dir, "counties-pumas", "counties-pumas.shp")
counties_sf <- st_read(sf_path)
adj_list <- poly2nb(counties_sf, queen = FALSE)
W <- nb2mat(adj_list, style = "B")
cat(sprintf("Shapefiles imported from: %s\n", sf_path)) # log

# Import crimes_sf
crimes_sf_path <- file.path(explain_boundaries_dir, "crimes_sf.dat")
load(crimes_sf_path)
cat(sprintf("Crimes shapefile imported from: %s\n", crimes_sf_path)) # log

# Import health_insurance_sf
health_insurance_sf_path <- file.path(explain_boundaries_dir, "health_insurance_sf.dat")
load(health_insurance_sf_path)
cat(sprintf("Health insurance shapefiles imported from: %s\n", health_insurance_sf_path)) # log

# Perform primary operations
intersections <- st_intersects(counties_sf, health_insurance_sf)
compute_mean_uninsured_pct <- function(idx) {
  intersected_geoms <- st_intersection(counties_sf$geometry[idx], health_insurance_sf$geometry)
  weights <- (st_area(intersected_geoms) / st_area(health_insurance_sf$geometry[intersections[[idx]]]))
  mean_uninsured_pct <- weighted.mean(health_insurance_sf$uninsured_pct[intersections[[idx]]], w = as.numeric(weights))
  return(mean_uninsured_pct)
}

# Add explain boundaries information to lacounty_sf
lacounty_idxs <- grep("Los Angeles County", counties_sf$Name)
lacounty_sf <- counties_sf[lacounty_idxs, ] %>%
  mutate(NumCrimes = sapply(st_contains(., crimes_sf), length)) %>%
  rowwise() %>%
  mutate(MeanUninsuredPct = compute_mean_uninsured_pct(cur_group_id()))

# Import final simulation and compute the boundary geometry
load(sim_file)
cat(sprintf("MCMC chain imported from: %s\n", sim_file)) # log
SPMIX_fit <- SPMIX_fit[35001:40000]
chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("spmix.UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute adjacency matrix
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute PPI matrix (remove not admissible edges)
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)

# Restrict to LA County
W_la <- W[lacounty_idxs, lacounty_idxs]
Gb_la <- Gb[lacounty_idxs, lacounty_idxs]

# Compute boundary geometry
if(!all(is.na(Gb_la))){
  bound_list <- apply(Gb_la, 1, function(x){which(x == 1)})
  bound_sf <- BoundaryGeometry(bound_list, lacounty_sf)
} else {
  cat("No Boundaries have been found\n")
}

# Get maps from stadia
la_bbox <- unname(st_bbox(st_transform(lacounty_sf, 4326)))
la_bbox[2] <- la_bbox[4] - (la_bbox[3] - la_bbox[1])
lacounty_map <- sf_ggmap(get_map(la_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))
lacounty_sf_3857 <- st_transform(lacounty_sf, 3857)
if(!all(is.na(Gb))){
  bound_sf_3857 <- st_transform(bound_sf, 3857)
}

# Plot - Number of crimes heatmap with detected boundaries in red
plt_numcrimes <- ggmap(lacounty_map) +
  geom_sf(data = lacounty_sf_3857, aes(fill = NumCrimes), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("N° of Crimes in 2020", direction = "horizontal", barwidth=unit(3,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
if(!all(is.na(Gb_la))){
   plt_numcrimes <- plt_numcrimes +
    geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = F)
}
pdf(file.path(output_dir, "plt_numcrimes.pdf"), height = 4, width = 4); print(plt_numcrimes); dev.off()

# Plot - Mean Precentage of uninsured population heatmap with detected boundaries in red
plt_uninsuredpct <- ggmap(lacounty_map) +
  geom_sf(data = lacounty_sf_3857, aes(fill = MeanUninsuredPct), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("Population w/o Health Insurance (%)", direction = "horizontal", barwidth=unit(3,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
if(!all(is.na(Gb_la))){
  plt_uninsuredpct <- plt_uninsuredpct +
    geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = F)
}
pdf(file.path(output_dir, "plt_uninsuredpct.pdf"), height = 4, width = 4); print(plt_uninsuredpct); dev.off()


# # ---- END OF SCRIPT ---- # #
