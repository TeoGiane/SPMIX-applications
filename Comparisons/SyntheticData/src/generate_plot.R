# # ---- Comparative Study on Synthetic Data - Generate Plot ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_plot", hide.opts = TRUE,
                         description = "Generates plots for synthetic data study")
opt_parser <- add_argument(opt_parser, arg = "--input-dir", type = "character", default = "input",
                           help = "Relative path to the input data directory")
opt_parser <- add_argument(opt_parser, arg = "--chains-dir", type = "character", default = "output",
                           help = "Relative path to the input data directory where the output from competitor models are stored")
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
  stop("Input directory: ", input_dir, " does not exist.")
}

# Check if chains directory exists
chains_dir <- file.path(getwd(), extra_args$chains_dir)
if (!dir.exists(chains_dir)) {
  stop("Chain directory: ", chains_dir, " does not exist.")
}

# Generate output directory if it does not exist
output_dir <- file.path(getwd(), extra_args$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory: ", output_dir, "\n")
}

# Required Libraries
suppressMessages(library("SPMIX"))
suppressMessages(library("CARBayes"))
suppressMessages(library("rjags"))
suppressMessages(library("spdep"))
suppressMessages(library("sf"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gridExtra"))

# Auxiliary functions
boundary_geometry <- function(boundary_list, sf_geometry) {
  if (!inherits(sf_geometry, "sf")) {
    stop("'sf_geometry' must be an sf object")
  }
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }

  # Create empty list
  geom_bdd <- list()

  for(i in 1:nrow(sf_geometry)) {
    # Get current area and its boundaries
    if (length(boundary_list[[i]]) > 0) {
      sel_geom <- sf_geometry[c(i, boundary_list[[i]]), ]

      # Compute geometry of boundary
      bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
      bounds <- st_geometry(bounds[bounds$id != bounds$id.1, ])

      # Add to list
      geom_bdd[[i]] <- st_sf(geometry = bounds)
    }
  }
  geom_bdd <- do.call(rbind, geom_bdd)

  # Drop points if present
  points <- which(attr(geom_bdd$geometry, "classes") == "POINT")
  if(length(points) > 0){
    geom_bdd <- geom_bdd[-points, ]
  }

  # Condense everything into a unique sf object and return
  return(geom_bdd)
}

# Load data
load(file.path(input_dir, "data.dat"))

# Load shapefile
grid_sf <- st_read(file.path(input_dir, "shp", "grid.shp"), quiet = TRUE)
st_crs(grid_sf) <- NA

# Compute adjacency list and matrix
adj_list <- poly2nb(grid_sf, queen = FALSE)
W <- nb2mat(adj_list, style = "B")
Eadj <- which(W == 1, arr.ind = TRUE)

# Models to compare
models <- c("CARBayes", "naiveMCAR", "SKATER", "SPMIX")

for(model in models){
  cat("Generating plot for model:", model, "\n")
  switch(model,
         "CARBayes" = {#
            load(file.path(chains_dir, "CARBayes-fit.dat"))
            plinks <- 1 - CARBayes_fit$localised.structure$W.posterior
            Gb <- ifelse(CARBayes_fit$localised.structure$W.posterior == 1, 1, NA)
          },
         "naiveMCAR" = {# 
            load(file.path(chains_dir, "naiveMCAR-fit.dat"))
            plinks <- summary(naiveMCAR_fit$G, FUN = mean)$stat
            plinks[which(W == 0, arr.ind = TRUE)] <- NA
            Gb <- matrix(NA, nrow(plinks), ncol(plinks))
            Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)
          },
         "SKATER" = {#
            load(file.path(chains_dir, "SKATER-fit.dat"))
            Gb <- matrix(NA, length(SKATER_fit$groups), length(SKATER_fit$groups))
            Gb[Eadj] <- ifelse(SKATER_fit$groups[Eadj[,1]] != SKATER_fit$groups[Eadj[,2]], 1, NA)
          },
         "SPMIX" = {#
            load(file.path(chains_dir, "SPMIX-fit.dat"))
            chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("UnivariateState",x))
            G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
            plinks <- Reduce('+', G_chain)/length(G_chain)
            plinks[which(W == 0, arr.ind = TRUE)] <- NA
            Gb <- matrix(NA, nrow(plinks), ncol(plinks))
            Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)
          },
         stop("Unknown model: ", model)
  )


# CARBayes
# load("output/CARBayes-chain-20250929-2320.dat")
# plinks <- 1 - CARBayes_fit$localised.structure$W.posterior
# Gb <- ifelse(CARBayes_fit$localised.structure$W.posterior == 1, 1, NA)

# # naiveMCAR
# load("output/naiveMCAR-chain-20250929-2334.dat")
# plinks <- summary(naiveMCAR_fit$G, FUN = mean)$stat
# plinks[which(W == 0, arr.ind = TRUE)] <- NA
# Gb <- matrix(NA, nrow(plinks), ncol(plinks))
# Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)

# # SKATER
# load("output/SKATER-output-20251002-1636.dat") # NON RICORDO UN CAZZO, CHE DUE COGLIONI
# Gb <- matrix(NA, length(SKATER_fit$groups), length(SKATER_fit$groups))
# Gb[Eadj] <- ifelse(SKATER_fit$groups[Eadj[,1]] != SKATER_fit$groups[Eadj[,2]], 1, NA)

# # SPMIX
# load("output/SPMIX-chain-20250929-2353.dat")
# chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
# G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
# plinks <- Reduce('+', G_chain)/length(G_chain)
# plinks[which(W == 0, arr.ind = TRUE)] <- NA
# Gb <- matrix(NA, nrow(plinks), ncol(plinks))
# Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)


  # Parte in comune
  # Compute boundary geometry
  if(!all(is.na(Gb))){
    bound_list <- apply(Gb, 1, function(x){which(x == 1)})
    bound_sf <- boundary_geometry(bound_list, grid_sf)
  } else {
    cat("No Boundaries have been found\n")
  }

  # PLOT - boundaries on the map in red + group
  plt_boundaries <- ggplot() +
    geom_sf(data = grid_sf, aes(fill=Group), col='gray25', alpha = 0.6, linewidth=0.4) +
    scale_fill_manual(values = c("steelblue","darkorange"), labels = c("Student's t", "Skew Normal"),
                      guide = guide_legend(title = "Data Generating Distribution", title.hjust = 0.5, label.position = "bottom",
                                          direction = "horizontal", title.position = "bottom", keywidth = unit(0.8,"in"))) +
    theme_void() + theme(legend.position = "bottom")
  if(!all(is.na(Gb))) {
    plt_boundaries <- plt_boundaries +
      geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2)
  }
  pdf(file.path(output_dir, paste0("plt_boundaries-", model, ".pdf")), height = 4, width = 4); print(plt_boundaries); dev.off()
  
  # PLOT - plinks matrix and with boundary edges in red
  if (model != "SKATER") {
    plinks_df <- reshape2::melt(plinks, c("x","y"), value.name="PPI") %>% na.omit()
    Gb_df <- reshape2::melt(Gb, c("x","y"), value.name="Gb") %>% na.omit()
    plt_plinks <- ggplot() +
      geom_tile(data = plinks_df, aes(x=x, y=y, fill=PPI), width=1, height=1) +
      geom_rect(aes(xmin=0.5, xmax=nrow(plinks)+0.5, ymin=0.5, ymax=nrow(plinks)+0.5), fill=NA, color="gray25", linewidth=0.5) +
      scale_fill_gradient2(low='steelblue4', mid = "white", high = 'darkorange', midpoint = 0.5, na.value = 'white',
                          guide = guide_colorbar("Post. Prob. of Inclusion", position = "bottom", direction = "horizontal", barwidth=unit(1.75,"in"),
                                                  title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
      theme_void() + theme(legend.position = "bottom") + coord_equal()
    if(!all(is.na(Gb))){
      plt_plinks <- plt_plinks +
        geom_tile(data = Gb_df, aes(x=x,y=y), fill=NA, col='darkred', linewidth=0.5)
    }
    pdf(file.path(output_dir, paste0("plt_plinks-", model, ".pdf")), height = 4, width = 4); print(plt_plinks); dev.off()
  }

}
