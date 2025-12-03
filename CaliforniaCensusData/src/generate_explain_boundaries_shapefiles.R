# # ---- CALIFORNIA CENSUS DATA - GENERATE EXPLAIN BOUNDARIES SHAPEFILES ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_explain_boundaries_shapefiles", hide.opts = TRUE,
                         description = "Generate shapefiles used for explain boundaries section in the manuscript.")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", default = "input/explain_boundaries",
                           help = "Output directory for generated shapefiles.")
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

# Create output directory if it doesn't exist
if (!dir.exists(extra_args$output_dir)) {
  dir.create(extra_args$output_dir, recursive = TRUE)
  cat("Created output directory:", extra_args$output_dir, "\n")
}


# Main code ---------------------------------------------------------------

# Import required packages
suppressMessages(library("httr"))
suppressMessages(library("jsonlite"))
suppressMessages(library("geojson"))
suppressMessages(library("geojsonsf"))
suppressMessages(library("sf"))
suppressMessages(library("dplyr"))

# Function to fetch from arcGIS
fetch_arcgis_features <- function(base_url, query_list, requested_page_size = 2000) {
  
  # Initialize
  all_features <- list()
  offset <- 0
  actual_page_size <- requested_page_size
  
  repeat {
    cat("Requesting", requested_page_size, "records from offset", offset, "...\n")
    
    # Build the query with paging fields
    extra_query_args <- list("resultRecordCount" = requested_page_size,
                             "resultOffset" = offset,
                             "f" = "geojson")
    full_query <- c(query_list, extra_query_args)
    
    # Perform request
    res <- GET(base_url, query = full_query)
    stop_for_status(res)
    curr_features <- geojson_sf(as.geojson(content(res, "text")))
    returned <- nrow(curr_features)
    cat(" → Server returned", returned, "features.\n")
    
    # Stop if there are no other records to dump
    if (returned == 0) {
      cat("No more records.\n")
      break
    }
    
    # Concatenate new features
    all_features <- append(all_features, list(curr_features))
    
    # Detect smaller server-imposed page size
    if (returned < requested_page_size) {
      actual_page_size <- returned
      cat(" → Adjusting future page size to", actual_page_size, "\n")
      requested_page_size <- actual_page_size
    }
    
    # Update offset
    offset <- offset + returned
  }
  
  # Bind all batches and return
  all_features <- do.call(rbind, all_features)
  cat("Total features downloaded:", nrow(all_features), "\n")
  return(all_features)
}

# Get Recorded Crimes Shapefile -------------------------------------------

# Define query
base_url <- "https://services.arcgis.com/RmCCgQtiZLDCtblq/arcgis/rest/services/PART_I_AND_II_CRIMES-HISTORICAL/FeatureServer/0/query"
query = list("where" = "INCIDENT_DATE >= DATE '2020-01-01 00:00:00' AND INCIDENT_DATE < DATE '2021-01-01 00:00:00'",
             "outFields" = "*",
             "returnGeometry" = "true")

# Get content as shapefile
crimes_sf_raw <- fetch_arcgis_features(base_url, query)

# Remove empty geometries
crimes_sf <- crimes_sf_raw %>%
  filter(!st_is_empty(geometry)) %>%
  st_transform(st_crs(st_read("input/counties-pumas/counties-pumas.shp")))

# Save output as .dat file
dir.create("input/explain_boundaries", recursive = T, showWarnings = F)
save(crimes_sf, file = "input/explain_boundaries/crimes_sf.dat")


# Get people without health insurance shapefile ---------------------------

# Define query
base_url <- "https://services.arcgis.com/RmCCgQtiZLDCtblq/arcgis/rest/services/Health_Insurance_tract/FeatureServer/0/query"
query = list("where" = "1=1",
             "outFields" = "*",
             "returnGeometry" = "true")

# Get content as shapefile
health_insurance_sf_raw <- fetch_arcgis_features(base_url, query)

# Remove empty geometries and transform
health_insurance_sf <- health_insurance_sf_raw %>%
  filter(!st_is_empty(geometry)) %>%
  st_transform(st_crs(st_read("input/counties-pumas/counties-pumas.shp")))

# Save output as .dat file
dir.create("input/explain_boundaries", recursive = T, showWarnings = F)
save(health_insurance_sf, file = "input/explain_boundaries/health_insurance_sf.dat")


# # ---- END OF SCRIPT ---- #

# counties_sf <- st_read("input/counties-pumas/counties-pumas.shp")
# Generate plots comes after
# num_crimes_sf <- st_read("input/counties-pumas/counties-pumas.shp") %>%
#   mutate(NumCrimes = sapply(st_contains(., crimes_sf), length))


