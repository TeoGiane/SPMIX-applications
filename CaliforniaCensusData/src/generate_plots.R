# # ---- CALIFORNIA CENSUS DATASET - GENERATE PLOTS ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_plot", hide.opts = TRUE,
                         description = "Generate plots for the analysis of the California Census Dataset")
opt_parser <- add_argument(opt_parser, arg = "--data-file", type = "character", default = NULL,
                           help = "Relative path to the data file to use.")
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
input_dir <- file.path(getwd(), "input")
if (!dir.exists(input_dir)) {
  stop(sprintf("Input directory '%s' does not exists.", input_dir))
}

# Check if data file exists
data_file <- file.path(getwd(), extra_args$data_file)
if(!file.exists(data_file)){
  stop(sprintf("Data file '%s' does not exists.", data_file))
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

# Import required packages
suppressMessages(library("ggplot2"))
suppressMessages(library("ggridges"))
suppressMessages(library("ggrepel"))
suppressMessages(library("ggmap"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))
suppressMessages(library("SPMIX"))
suppressMessages(library("dplyr"))

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

boundary_geometry <- function(boundary_list, sf_geometry) {
  if (!inherits(sf_geometry, "sf")) {
    stop("'sf_geometry' must be an sf object")
  }
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }

  # Create empty list
  geom_bdd <- list()
  count <- 1
  # Populate
  for(i in 1:nrow(sf_geometry)) {
    # Get current area and its boundaries
    if (length(boundary_list[[i]]) > 0) {
      for (j in boundary_list[[i]]) {
        # Compute geometry of boundary
        sel_geom <- sf_geometry[c(i, j),"id"]
        bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
        bounds <- st_union(st_geometry(bounds[bounds$id != bounds$id.1, ]))
        # Add to list
        geom_bdd[[count]] <- st_sf(geometry = bounds)
        count <- count + 1
      }
    }
  }
  # Bind all objects
  geom_bdd <- do.call(rbind, geom_bdd)

  # Drop points if present
  points <- which(attr(geom_bdd$geometry, "classes") == "POINT")
  if(length(points) > 0){
    geom_bdd <- geom_bdd[-points, ]
  }

  # Condense everything into a unique sf object and return
  return(geom_bdd)
}

L1_distance <- function(y1, y2, x) {
  I <- diff(range(x))
  return(I * mean(abs(y1 - y2)))
}

# Import shapefile
sf_file <- file.path(input_dir, "counties-pumas", "counties-pumas.shp")
sf_counties <- read_sf(sf_file)
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")
cat(sprintf("Shapefiles imported from: %s\n", sf_file)) # log

# Import data
load(data_file)
cat(sprintf("Data imported from: %s\n", data_file)) # log

# Load output
load(sim_file)
cat(sprintf("MCMC chain imported from: %s\n", sim_file)) # log

# Deserialization
SPMIX_fit <- SPMIX_fit[35001:40000]
chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("spmix.UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
Nedge_chain <- sapply(G_chain, function(x){sum(x[upper.tri(x)])})
sigma_chain <- sapply(chains, function(x){x$Sigma$data[1]})
p_chain <- sapply(chains, function(x){x$p})

# Compute posterior mean and variance for each PUMA
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$mean))
vars_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$stdev^2))
weights_chain <- lapply(chains, function(x) t(sapply(x$groupParams, function(y) y$weights)))
post_means <- matrix(nrow = nrow(sf_counties), ncol=length(chains))
post_vars <- matrix(nrow = nrow(sf_counties), ncol=length(chains))
for (j in 1:length(chains)) {
  post_means[,j] <- weights_chain[[j]] %*% means_chain[[j]]
  second_moment <- as.vector(weights_chain[[j]] %*% (vars_chain[[j]] + means_chain[[j]]^2))
  post_vars[,j] <- second_moment - post_means[,j]^2
}

# Add to shapefile dataframe
sf_counties$emp_mean <- sapply(data, mean)
sf_counties$post_mean <- apply(post_means, 1, mean)
sf_counties$emp_var <- sapply(data, var)
sf_counties$post_var <- apply(post_vars, 1, mean)

# Compute estimated density
x <- seq(range(data)[1], range(data)[2], length.out = 500)
estimated_densities <- lapply(ComputePredictiveLPDFs(SPMIX_fit, x), exp)

# Compute admissible edges
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute PPI matrix (remove not admissible edges)
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Compute neighbouring graph
Gn <- matrix(NA, nrow(plinks), ncol(plinks))
Gn[Eadj] <- ifelse(plinks[Eadj] >= 0.5, 1, NA)
neigh_list <- apply(Gn, 1, function(x){which(x==1)})

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1, NA)

# Compute boundary geometry
if(!all(is.na(Gb))){
  bound_list <- apply(Gb, 1, function(x){which(x == 1)})
  bound_sf <- boundary_geometry(bound_list, sf_counties)
} else {
  cat("No Boundaries have been found\n")
}

# L1 distance between areas - Local Comparison
L1 <- data.frame("bounds" = rep(NA, nrow(sf_counties)),
                 "no_bounds" = rep(NA, nrow(sf_counties)))
for (i in 1:nrow(sf_counties)) {
  p <- colMeans(estimated_densities[[i]])
  q_b <- lapply(estimated_densities[bound_list[[i]]], function(x) { colMeans(x) })
  q_nb <- lapply(estimated_densities[neigh_list[[i]]], function(x) { colMeans(x) })
  if(length(q_b) > 0) {
    L1[i, "bounds"] <- mean(sapply(q_b, function(y){ L1_distance(p, y, x) }))
  }
  if(length(q_nb) > 0) {
    L1[i, "no_bounds"] <- mean(sapply(q_nb, function(y){ L1_distance(p, y, x) }))
  }
}

# PLOT - boxplot comparison (local)
L1loc <- rbind(data.frame("Dist" = L1$no_bounds, "Type" = "Neigh"),
               data.frame("Dist" = L1$bounds, "Type" = "Bound"))
L1loc$Type <- as.factor(L1loc$Type)
L1loc <- na.omit(L1loc)
plt_L1loc <- ggplot() +
  geom_boxplot(data = L1loc, aes(x=Type, y=Dist)) +
  geom_boxplot(data = L1loc, aes(x=Type, y=Dist, fill=Type, color=Type), staplewidth = 0.3, alpha = 0.3, show.legend = F) +
  scale_x_discrete(labels = c(bquote(d[hat(BE)[~loc]]), bquote(d[hat(NE)[~loc]]))) + labs(x=NULL,y=NULL) +
  scale_fill_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  scale_color_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  theme(text = element_text(size = 14))
pdf(file.path(output_dir, "plt_L1loc.pdf"), height = 3, width = 4); print(plt_L1loc); dev.off()

# L1 distance between areas - Global Comparison
Gn_up <- Gn; Gn_up[lower.tri(Gn_up)] <- NA
neigh_pairs <- which(Gn_up == 1, arr.ind = T)
Gb_up <- Gb; Gb_up[lower.tri(Gb_up)] <- NA
bound_pairs <- which(Gb_up == 1, arr.ind = T)
L1_neigh <- data.frame("Type" = rep("Neigh", nrow(neigh_pairs)),
                       "Dist" = rep(NA, nrow(neigh_pairs)))
for (i in 1:nrow(neigh_pairs)) {
  L1_neigh$Dist[i] <- L1_distance(colMeans(estimated_densities[[neigh_pairs[i,1]]]),
                                  colMeans(estimated_densities[[neigh_pairs[i,2]]]), x)
}
L1_bound <- data.frame("Type" = rep("Bound", nrow(bound_pairs)),
                       "Dist" = rep(NA, nrow(bound_pairs)))
for (i in 1:nrow(bound_pairs)) {
  L1_bound$Dist[i] <- L1_distance(colMeans(estimated_densities[[bound_pairs[i,1]]]),
                                  colMeans(estimated_densities[[bound_pairs[i,2]]]), x)
}

# PLOT - boxplot comparison (global)
L1glob <- rbind(L1_neigh, L1_bound); rm(L1_neigh, L1_bound)
L1glob$Type <- as.factor(L1glob$Type)
plt_L1glob <- ggplot() +
  geom_boxplot(data = L1glob, aes(x=Type, y=Dist)) +
  geom_boxplot(data = L1glob, aes(x=Type, y=Dist, fill=Type, color=Type), staplewidth = 0.3, alpha = 0.3, show.legend = F) +
  scale_x_discrete(labels = c(bquote(d[hat(BE)]), bquote(d[hat(NE)]))) + labs(x=NULL,y=NULL) +
  scale_fill_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  scale_color_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  theme(text = element_text(size = 14))
pdf(file.path(output_dir, "plt_L1glob.pdf"), height = 3, width = 4); print(plt_L1glob); dev.off()

# PLOT - Traceplot of |G|
df <- data.frame("Iteration" = 1:length(Nedge_chain), "Value" = Nedge_chain)
plt_Nedge_chain <- ggplot(data = df) + 
  geom_line(aes(x=Iteration, y=Value)) + ylab("|G|")
pdf(file.path(output_dir, "plt_Nedge_chain.pdf"), height = 4, width = 4); print(plt_Nedge_chain); dev.off()

# PLOT - Traceplot of sigma^2
df <- data.frame("Iteration" = 1:length(sigma_chain), "Value" = sigma_chain)
plt_sigma_chain <- ggplot(data = df) + 
  geom_line(aes(x=Iteration, y=Value)) + ylab(bquote(sigma^2))
pdf(file.path(output_dir, "plt_sigma_chain.pdf"), height = 4, width = 4); print(plt_sigma_chain); dev.off()

# PLOT - traceplot of p
df <- data.frame("Iteration" = 1:length(p_chain), "Value" = p_chain)
plt_p_chain <- ggplot(data = df) + 
  geom_line(aes(x=Iteration, y=Value)) + ylab("p")
pdf(file.path(output_dir, "plt_p_chain.pdf"), height = 4, width = 4); print(plt_p_chain); dev.off()

# PLOT - Posterior distribution of H
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plt_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  xlab("N° of Components")
pdf(file.path(output_dir, "plt_postH.pdf"), height = 4, width = 4); print(plt_postH); dev.off()

# PLOT - Traceplot of H
df <- data.frame("Iteration"=1:length(chains), "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plt_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
pdf(file.path(output_dir, "plt_traceH.pdf"), height = 4, width = 4); print(plt_traceH); dev.off()

# PLOT - plinks matrix and with boundary edges in red
plinks_df <- reshape2::melt(plinks, c("x", "y"), value.name="PPI") %>% na.omit()
Gb_df <- reshape2::melt(Gb, c("x", "y"), value.name="Gb") %>% na.omit()
plt_plinks <- ggplot() +
  geom_tile(data = plinks_df, aes(x=x, y=y, fill=PPI), width=1, height=1) +
  geom_rect(aes(xmin=0.5, xmax=nrow(plinks)+0.5, ymin=0.5, ymax=nrow(plinks)+0.5), fill=NA, color="gray25", linewidth=0.5) +
  scale_fill_gradient2(low='steelblue4', mid = "white", high = 'darkorange', midpoint = 0.5, na.value = 'white',
                       guide = guide_colorbar("Post. Prob. of Inclusion", position = "bottom", direction = "horizontal", barwidth=unit(2.5,"in"),
                                              title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom") + coord_equal()
if(!all(is.na(Gb))){
  plt_plinks <- plt_plinks + geom_tile(data = Gb_df, aes(x=x,y=y), fill=NA, col='darkred', linewidth=0.5)
}
pdf(file.path(output_dir, "plt_plinks.pdf"), height = 4, width = 4); print(plt_plinks); dev.off()

# PLOT - Posterior Summary graph
# df <- reshape2::melt(G_sum, c("x","y"), value.name = "Val"); df[which(df$Val == 0), "Val"] <- NA; df$Val <- as.factor(df$Val)
# # Generate
# plt_Gsum <- ggplot() +
#   geom_tile(data = df, aes(x=x, y=y, fill=Val), width=1, height=1) +
#   geom_rect(aes(xmin=0.5, xmax=nrow(G_est)+0.5, ymin=0.5, ymax=nrow(G_est)+0.5), fill=NA, color='gray25', linewidth=0.5) +
#   scale_fill_manual(labels = c("Boundary edge", "Graph edge", NULL),
#                     values = c("darkred", "darkorange"),
#                     na.translate = F,
#                     guide = guide_legend(title="", direction = "horizontal", title.position = "bottom",
#                                          label.position = "bottom", keywidth = unit(3,"cm"), )) +
#   theme_void() + theme(legend.position = "bottom") + coord_equal()
# # Show
# plt_Gsum

# # PLOT - Posterior boundary graph
# df <- reshape2::melt(G_b, c("x","y"), value.name = "Val")
# # Generate
# plt_Gbound <- ggplot() +
#   geom_tile(data = df, aes(x=x, y=y, fill=Val), width=1, height=1) +
#   geom_rect(aes(xmin=0.5, xmax=nrow(G_est)+0.5, ymin=0.5, ymax=nrow(G_est)+0.5), fill=NA, color='gray25', linewidth=0.5) +
#   scale_fill_gradient(low='white', high='darkred',
#                       guide = guide_legend(title="", direction = "horizontal", title.position = "bottom",
#                                            label.position = "bottom", override.aes = list(alpha = 0))) +
#   theme_void() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent")) + coord_equal()
# # Show
# plt_Gbound

# Get maps from stadia
counties_bbox <- unname(st_bbox(st_transform(sf_counties, 4326)))
counties_map <- sf_ggmap(get_map(counties_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))
sf_counties_3857 <- st_transform(sf_counties, 3857)
if(!all(is.na(Gb))){
  bound_sf_3857 <- st_transform(bound_sf, 3857)
}

# PLOT - Posterior mean heatmap + detected boundaries
plt_boundaries_mean <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("Post. Mean", direction = "horizontal", barwidth=unit(3,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
if(!all(is.na(Gb))){
  plt_boundaries_mean <- plt_boundaries_mean +
    geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = FALSE)
}
pdf(file.path(output_dir, "plt_boundaries_mean.pdf"), height = 4, width = 4); print(plt_boundaries_mean); dev.off()

# PLOT - Posterior variance heatmap + detected boundaries
plt_boundaries_var <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_var), col='gray25', alpha = 0.6, inherit.aes = FALSE) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("Post. Variance", direction = "horizontal", barwidth=unit(3,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
if(!all(is.na(Gb))){
  plt_boundaries_var <- plt_boundaries_var +
    geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = FALSE)
}
pdf(file.path(output_dir, "plt_boundaries_var.pdf"), height = 4, width = 4); print(plt_boundaries_var); dev.off()

# PLOT - Empirical mean in each PUMA on the map
plt_emp_mean <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=emp_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colourbar(title = "Empirical Mean", direction = "horizontal", barwidth = unit(3, "in"),
                                              title.position = "bottom", title.hjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
pdf(file.path(output_dir, "plt_emp_mean.pdf"), height = 4, width = 4); print(plt_emp_mean); dev.off()

# PLOT - Empirical variance in each PUMA on the map
plt_emp_var <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=emp_var), col='gray25', alpha = 0.6, inherit.aes = FALSE) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colourbar(title = "Empirical Variance", direction = "horizontal", barwidth = unit(3, "in"),
                                              title.position = "bottom", title.hjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
pdf(file.path(output_dir, "plt_emp_var.pdf"), height = 4, width = 4); print(plt_emp_var); dev.off()

# PLOT - Empirical density histogram in bordering areas
areas <- c(58,60,48)
areas_names <- c("Gardena,\nLawndale\n& West Athens", "Redondo Beach,\nManhattan Beach\n& Hermosa Beach", "Marina del Rey,\nWestchester\n& Cluver City")
df_hist <- data.frame("logPINCP"=numeric(0), "PUMA"=numeric(0))
df_dens <- data.frame("x"=numeric(0), "y"=numeric(0), "PUMA"=numeric(0))
for (i in 1:length(areas)) {
  to_add <- data.frame("logPINCP" = data[[areas[i]]],
                       "PUMA" = rep(areas_names[i], length(data[[areas[i]]])))
  df_hist <- rbind(df_hist, to_add)
  to_add <- data.frame("x" = x,
                       "y" = colMeans(estimated_densities[[areas[i]]]),
                       "PUMA" = rep(areas_names[i], length(x)))
  df_dens <- rbind(df_dens, to_add)
}
df_hist$PUMA <- factor(df_hist$PUMA, levels = areas_names)
df_dens$PUMA <- factor(df_dens$PUMA, levels = areas_names)
plt_DensCompare <- ggplot() +
  geom_density_ridges(data = df_hist, aes(x=logPINCP, y=PUMA, height=after_stat(density), fill=PUMA, color=PUMA), stat="binline", bins=50, alpha=0.4, show.legend = F) +
  geom_ridgeline(data=df_dens, aes(x=x, y=PUMA, height=y, color=PUMA), scale=4.5, fill=NA, linewidth = 1.2, show.legend = F) +
  scale_y_discrete(expand = c(0.05,0,1,0)) +
  scale_color_manual(NULL, values = c('darkorange', "steelblue", 'forestgreen')) +
  scale_fill_manual(NULL, values = c('darkorange', "steelblue", 'forestgreen')) +
  xlab("log(PINCP)") + ylab("Density") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf(file.path(output_dir, "plt_DensCompare.pdf"), height = 4, width = 4); print(plt_DensCompare); dev.off()

# PLOT - Zoom on the map
zoom_bbox <- unname(st_bbox(st_transform(sf_counties[areas,], 4326)))
zoom_map <- sf_ggmap(get_map(zoom_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))
names_zoom <- data.frame("PUMA" = areas_names, st_coordinates(st_centroid(st_geometry(sf_counties_3857[areas, ]))))
names_zoom$PUMA <- factor(names_zoom$PUMA, levels = areas_names)#, levels = c("U.S.C. & Exposition Park", "Hancock Park & Mid-Wilshire", "West Hollywood & Beverly Hills"))
plt_zoom <- ggmap(zoom_map) +
  geom_sf(data = sf_counties_3857, fill=NA, col='gray25', linewidth=0.7, inherit.aes = F) +
  geom_sf(data = sf_counties_3857[areas, ], aes(fill=PUMA), col='gray25', alpha = 0.5, linewidth = 0.7, show.legend = F, inherit.aes = F) +
  geom_sf(data = bound_sf_3857, fill=NA, col='darkred', linewidth=1.4, inherit.aes = F) +
  scale_fill_manual(values = c("03758"='darkorange',"03760"='steelblue',"03748"='forestgreen')) +
  geom_label_repel(data = names_zoom, aes(x=X, y=Y, label=PUMA, color=PUMA), inherit.aes = F, show.legend = F) +
  scale_color_manual(values = c('darkorange', 'steelblue', 'forestgreen')) + theme_void()
pdf(file.path(output_dir, "plt_zoom.pdf"), height = 5, width = 5); print(plt_zoom); dev.off()

# Final log
cat(sprintf("All plots have been generated and saved in %s directory\n", output_dir))

# # ---- END OF SCRIPT ---- # #
