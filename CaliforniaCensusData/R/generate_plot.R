# Import required packages
suppressMessages(library("ggplot2"))
suppressMessages(library("ggridges"))
suppressMessages(library("ggrepel"))
suppressMessages(library("ggmap"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))
suppressMessages(library("SPMIX"))
# suppressMessages(library("dplyr"))

# setwd("~/Documents/GitHub/SPMIX-applications/CaliforniaCensusData")


# Argument parser to select from terminal the simulation to parse
# library("optparse")
# option_list <- list (
#   make_option(c("--sim_name"), type = "character", default = NULL,
#               help = "Name of the simulation to use. If the directory 'output/<sim_name>' does not exist,
#                       it will rise an error.",
#               metavar="character")
# )
# opt_parser <- OptionParser(option_list=option_list)
# args <- parse_args(opt_parser)
# sim_name <- args$sim_name
# sim_name <- "sim28"
# sim_name <- "rho0"

sim_folder <- file.path(getwd(), "output/H_RJ/rho_0.95")

# log
cat(sprintf("Current Directory: %s\n", getwd()))

# Check is sim_folder exists
# sim_folder <- sprintf("%s/output/%s", getwd(), sim_name)
# if(!dir.exists(sim_folder)){
#   stop(sprintf("'%s' folder does not exists.", sim_folder))
# }

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
  if(!inherits(sf_geometry, "sf")) { stop("'sf_geometry' must be an sf object") }
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }
  
  # Create empty list
  geom_bdd <- list()
  
  for(i in 1:nrow(sf_geometry)) {
    # Get current area and its boundaries
    sel_geom <- sf_geometry[c(i, boundary_list[[i]]), c("id", "geometry")]
    
    # Compute geometry of boundary
    bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
    bounds <- st_geometry(bounds[bounds$id != bounds$id.1, ])
    
    # Add to list
    geom_bdd[[i]] <- st_sf(geometry = bounds)
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

# bottom_colorbar <- function(title, length){
#   out <- guide_colorbar(title, direction = "horizontal", barwidth = unit(length, "in"), barheight = unit(0.25, "in"),
#                         title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)
#   return(out)
# }
# 
# bottom_legend <- function(title, nrow){
#   out <- guide_legend(title, direction = "horizontal", title.position = "bottom",
#                       title.hjust = 0.5, nrow = nrow, override.aes = list(size=2))
#   return(out)
# }

L1_distance <- function(y1, y2, x) {
  I <- diff(range(x))
  return(I * mean(abs(y1 - y2)))
}

# Import shapefile
sf_counties <- read_sf("data/counties-pumas/counties-pumas.shp")

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Import data
load("data/data_001.dat")




# Import usa shapefile
# sf_usa <- read_sf("shp/us-pumas/us-pumas.shp")
# sf_usa$id <- row.names(sf_usa)
# sf_usa$PUMA <- as.character(as.numeric(sf_usa$PUMA))

# Select only LA + Ventura + Orange counties, California
# sf_cali <- sf_usa[sf_usa$State == "California", ]
# sel_county <- c("Los Angeles County", "Ventura County", "Orange County")
# sf_counties <- sf_cali[grep(paste(sel_county, collapse = "|"), sf_cali$Name), ]
# pumas <- sf_counties$PUMA

# Import data and adjacency matrix
# load("data/clean_data.dat")
# load("data/adj_matrix.dat")
# adj_list <- apply(W, 1, function(x){which(x==1)})
cat("Shapefiles have been imported\n") # log

# Load output
out_file = file.path(sim_folder, "chain_001.dat")
load(out_file)
cat(sprintf("Loaded chain: %s\n", out_file)) # log

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Some chains to inspect
Nedges_chain <- sapply(G_chain, sum); mcmcse::ess(Nedges_chain)
sigma_chain <- sapply(chains, function(x){x$Sigma$data[1]}); mcmcse::ess(sigma_chain)
p_chain <- sapply(chains, function(x){x$p}); mcmcse::ess(p_chain)

# Visualization of traceplots - |G|
df <- data.frame("Iteration"=1:length(Nedges_chain), "Value"=Nedges_chain)
# Generate
plt_Nedges <- ggplot(data = df) + 
  geom_line(aes(x=Iteration, y=Value)) + ylab("|G|")
# Show / Save
x11(height = 3, width = 3); plt_Nedges
# pdf("plots/plt_Nedges.pdf", height = 3, width = 3); plt_Nedges; dev.off()

# Visualization of traceplots - sigma
df <- data.frame("Iteration"=1:length(sigma_chain), "Value"=sigma_chain)
# Generate
plt_sigma <- ggplot(data = df) + 
  geom_line(aes(x=Iteration, y=Value)) + ylab(bquote(sigma^2))
# Show / Save
x11(height = 3, width = 3); plt_sigma
# pdf("plots/plt_sigma.pdf", height = 3, width = 3); plt_sigma; dev.off()

# Visualization of traceplots - p
df <- data.frame("Iteration"=1:length(p_chain), "Value"=p_chain)
# Generate
plt_p <- ggplot(data = df) + 
  geom_line(aes(x=Iteration, y=Value)) + ylab("p")
# Show / Save
x11(height = 3, width = 3); plt_p
# pdf("plots/plt_p.pdf", height = 3, width = 3); plt_p; dev.off()

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
estimated_densities <- ComputeDensities(chains, x, verbose = TRUE)

# Compute admissible edges
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute PPI matrix (remove not admissible edges)
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Threshold -> posterior median of p
gamma_graph <- 0.5

# Compute neighbouring graph
Gn <- matrix(NA, nrow(plinks), ncol(plinks))
Gn[Eadj] <- ifelse(plinks[Eadj] >= gamma_graph, 1, NA)

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < gamma_graph, 1, NA)

# Compute neighbouring and boundary adjacency lists
neigh_list <- apply(Gn, 1, function(x){which(x==1)})
bound_list <- apply(Gb, 1, function(x){which(x==1)})

# Compute boundary geometry
bound_sf <- boundary_geometry(bound_list, sf_counties)


# L1 distance between areas - V1 ------------------------------------------

# Prepare buffer
L1 <- data.frame("bounds" = rep(NA, nrow(sf_counties)),
                 "no_bounds" = rep(NA, nrow(sf_counties)))

# Compute L1 distances between densities
for (i in 1:nrow(sf_counties)) {
  p <- colMeans(estimated_densities[[i]])
  q <- lapply(estimated_densities[adj_list[[i]]], function(x) { colMeans(x) })
  q_b <- lapply(estimated_densities[bound_list[[i]]], function(x) { colMeans(x) })
  q_nb <- lapply(estimated_densities[neigh_list[[i]]], function(x) { colMeans(x) })
  mean_dist <- mean(sapply(q, function(y) { L1_distance(p, y, x) }))#, na.rm = T)
  #print(mean_dist)
  # q_b <- lapply(estimated_densities[bound_list[[i]]], function(x) {x['est',]})
  if(length(q_b) > 0) {
    L1[i, "bounds"] <- mean(sapply(q_b, function(y){ L1_distance(p, y, x) }))
      # abs(mean(sapply(q_b, function(x) { sum(abs(x - p)) }), na.rm = T) - mean_dist)
  }
  # q_nb <- lapply(estimated_densities[nobound_list[[i]]], function(x) {x['est',]})
  if(length(q_nb) > 0) {
    L1[i, "no_bounds"] <- mean(sapply(q_nb, function(y){ L1_distance(p, y, x) }))
      # abs(mean(sapply(q_nb, function(x) { sum(abs(x - p)) }), na.rm = T) - mean_dist)
  }
}

# Visualization - Paired boxplots
L1_v1 <- rbind(data.frame("Dist" = L1$no_bounds, "Type" = "Neigh"),
               data.frame("Dist" = L1$bounds, "Type" = "Bound")); rm(L1)
L1_v1$Type <- as.factor(L1_v1$Type)
L1_v1 <- na.omit(L1_v1)
# Generate plot
plt_L1v1 <- ggplot() +
  geom_boxplot(data = L1_v1, aes(x=Type, y=Dist)) +
  geom_boxplot(data = L1_v1, aes(x=Type, y=Dist, fill=Type, color=Type), staplewidth = 0.3, alpha = 0.3, show.legend = F) +
  scale_x_discrete(labels = c(bquote(d[FM[i]]), bquote(d[TM[i]]))) + labs(x=NULL,y=NULL) +
  scale_fill_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  scale_color_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  theme(text = element_text(size = 14))
# Show
x11(height = 3, width = 4); plt_L1v1
# Save
# pdf("plots/plt_L1v1.pdf", height = 3, width = 4); plt_L1v1; dev.off()
  

# L1 distance between areas - V2 ------------------------------------------

# Neighbouring pairs
Gn_up <- Gn; Gn_up[lower.tri(Gn_up)] <- NA
neigh_pairs <- which(Gn_up == 1, arr.ind = T)

# Boundary pairs
Gb_up <- Gb; Gb_up[lower.tri(Gb_up)] <- NA
bound_pairs <- which(Gb_up == 1, arr.ind = T)

# Compute L1 distance between all neighbouring pairs
L1_neigh <- data.frame("Type" = rep("Neigh", nrow(neigh_pairs)),
                       "Dist" = rep(NA, nrow(neigh_pairs)))
for (i in 1:nrow(neigh_pairs)) {
  L1_neigh$Dist[i] <- L1_distance(colMeans(estimated_densities[[neigh_pairs[i,1]]]),
                                  colMeans(estimated_densities[[neigh_pairs[i,2]]]), x)
}

# Compute L1 distance between all boundary pairs
L1_bound <- data.frame("Type" = rep("Bound", nrow(bound_pairs)),
                       "Dist" = rep(NA, nrow(bound_pairs)))
for (i in 1:nrow(bound_pairs)) {
  L1_bound$Dist[i] <- L1_distance(colMeans(estimated_densities[[bound_pairs[i,1]]]),
                                  colMeans(estimated_densities[[bound_pairs[i,2]]]), x)
}


# Visualization - Paired boxplots
L1_v2 <- rbind(L1_neigh, L1_bound); rm(L1_neigh, L1_bound)
L1_v2$Type <- as.factor(L1_v2$Type)
# Generate plot
plt_L1v2 <- ggplot() +
  geom_boxplot(data = L1_v2, aes(x=Type, y=Dist)) +
  geom_boxplot(data = L1_v2, aes(x=Type, y=Dist, fill=Type, color=Type), staplewidth = 0.3, alpha = 0.3, show.legend = F) +
  scale_x_discrete(labels = c(bquote(d[FM]), bquote(d[TM]))) + labs(x=NULL,y=NULL) +
  scale_fill_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  scale_color_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  theme(text = element_text(size = 14))
# Show
x11(height = 3, width = 4); plt_L1v2
# Save
# pdf("plots/plt_L1v2.pdf", height = 3, width = 4); plt_L1v2; dev.off()

cat("Posterior analysis has been completed\n")

# PLOT - Posterior distribution of H
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
# Generate
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components")
# Show
x11(height = 4, width = 4); plot_postH


# PLOT - PUMAs with names (if I need to select specific areas)
text <- data.frame("PUMA" = sf_counties$PUMA,  st_coordinates(st_centroid(st_geometry(sf_counties))))
# Generate
plt_pumanames <- ggplot() + 
  geom_sf(data = sf_counties) + 
  geom_text(data = text, aes(x = X, y = Y, label=PUMA), check_overlap = TRUE)
# Show
x11(height = 4, width = 4); plt_pumanames


# PLOT - plinks matrix and with boundary edges in red
plinks_df <- reshape2::melt(plinks, c("x", "y"), value.name="PPI") %>% na.omit()
Gb_df <- reshape2::melt(Gb, c("x","y"), value.name="Gb") %>% na.omit()
# Generate
plt_plinks <- ggplot() +
  geom_tile(data = plinks_df, aes(x=x, y=y, fill=PPI), width=1, height=1) +
  geom_rect(aes(xmin=0.5, xmax=nrow(plinks)+0.5, ymin=0.5, ymax=nrow(plinks)+0.5), fill=NA, color="gray25", linewidth=0.5) +
  scale_fill_gradient2(low='steelblue4', mid = "white", high = 'darkorange', midpoint = 0.5, na.value = 'white',
                       guide = guide_colorbar("Post. Prob. of Inclusion", position = "bottom", direction = "horizontal", barwidth=unit(2.5,"in"),
                                              title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  geom_tile(data = Gb_df, aes(x=x,y=y), fill=NA, col='darkred', linewidth=0.5) +
  theme_void() + theme(legend.position = "bottom") + coord_equal()
# Show
x11(height = 4, width = 4); plt_plinks
# Save
# pdf("plots/plt_plinks.pdf", height = 4, width = 4); plt_plinks; dev.off()

# # PLOT - Posterior Summary graph
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
bound_sf_3857 <- st_transform(bound_sf, 3857)

# PLOT - Posterior mean heatmap + detected boundaries
# Generate
plt_boundaries_mean <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("Post. Mean", direction = "horizontal", barwidth=unit(2.5,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = FALSE) +
  theme_void() + theme(legend.position = "bottom")
# Show
x11(height = 4, width = 4); plt_boundaries_mean
# Save
# pdf("plots/plt_boundaries_mean.pdf", height = 4, width = 4); plt_boundaries_mean; dev.off()


# PLOT - Posterior variance heatmap + detected boundaries
# Generate
plt_boundaries_var <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_var), col='gray25', alpha = 0.6, inherit.aes = FALSE) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("Post. Variance", direction = "horizontal", barwidth=unit(2.5,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = FALSE) +
  theme_void() + theme(legend.position = "bottom")
# Show
x11(height = 4, width = 4); plt_boundaries_var
# Save
# pdf("plots/plt_boundaries_var.pdf", height = 4, width = 4); plt_boundaries_var; dev.off()


# PLOT - Empirical mean in each PUMA on the map
plt_emp_mean <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=emp_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colourbar(title = "Empirical Mean", direction = "horizontal", barwidth = unit(2.5, "in"),
                                              title.position = "bottom", title.hjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
# Show
x11(height = 4, width = 4); plt_emp_mean
# Save
pdf("plots/plt_emp_mean.pdf", height = 4, width = 4); plt_emp_mean; dev.off()

# PLOT - Empirical variance in each PUMA on the map
plt_emp_var <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=emp_var), col='gray25', alpha = 0.6, inherit.aes = FALSE) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colourbar(title = "Empirical Variance", direction = "horizontal", barwidth = unit(3, "in"),
                                              title.position = "bottom", title.hjust = 0.5)) +
  theme_void() + theme(legend.position = "bottom")
# Show
x11(height = 4, width = 4); plt_emp_var
# Save
pdf("plots/plt_emp_var.pdf", height = 4, width = 4); plt_emp_var; dev.off()


# PLOT - Empirical density histogram in bordering areas
areas <- c(30,46,31); areas_names <- c("Hancock Park & Mid-Wilshire", "U.S.C. & Exposition Park", "West Hollywood & Beverly Hills")
# Auxiliary dataframes
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
df_hist$PUMA <- factor(df_hist$PUMA, levels = c("U.S.C. & Exposition Park", "Hancock Park & Mid-Wilshire", "West Hollywood & Beverly Hills"))
df_dens$PUMA <- factor(df_dens$PUMA, levels = c("U.S.C. & Exposition Park", "Hancock Park & Mid-Wilshire", "West Hollywood & Beverly Hills"))
# Generate plot
plt_DensCompare <- ggplot() +
  geom_density_ridges(data = df_hist, aes(x=logPINCP, y=PUMA, height=after_stat(ndensity), fill=PUMA, color=PUMA, scale=1.5), stat="binline", bins=10, alpha=0.4, show.legend = F) +
  geom_ridgeline(data=df_dens, aes(x=x, y=PUMA, height=y, color=PUMA), scale=4.5, fill=NA, linewidth = 1.2, show.legend = F) +#, size=1.2, show.legend = F) +
  scale_y_discrete(expand = c(0.05,0,0.9,0)) +
  scale_color_manual(NULL, values = c('darkorange', "steelblue", 'forestgreen')) +
  scale_fill_manual(NULL, values = c('darkorange', "steelblue", 'forestgreen')) +
  xlab("log(PINCP)") + ylab("Density") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# Show
x11(height = 4, width = 4); plt_DensCompare
# Save
# pdf("plots/plt_DensCompare.pdf", height = 4, width = 4); plt_DensCompare; dev.off()


# plt_areas <- list(); nbins <- 9
# for (i in 1:length(areas)) {
#   df <- data.frame("logPINCP" = data[[areas[i]]])
#   plt_areas[[i]] <- ggplot(data = df, aes(x=logPINCP, y = after_stat(density))) +
#     geom_histogram(col=NA, fill='white', bins = nbins) +
#     geom_histogram(col='steelblue', fill='steelblue', alpha = 0.4, bins = nbins) +
#     xlab("log(PINCP)") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
# }
# # Show plot
# plt_areas[[1]]; plt_areas[[2]]
# Save plot
# pdf("output/plt_areas_empirical.pdf", height = 3, width = 6)
# gridExtra::grid.arrange(grobs = plt_areas, ncol=2)
# dev.off()

# PLOT - Empirical density histogram + Estimated density
# plt_areasdens <- list()
# for (i in 1:length(areas)) {
#   df <- data.frame("x" = seq(data_range[1,areas[i]], data_range[2,areas[i]], length.out = Npoints),
#                    "y" = estimated_densities[[areas[i]]]['est', ],
#                    "ymin" = estimated_densities[[areas[i]]]['low', ],
#                    "ymax" = estimated_densities[[areas[i]]]['up', ])
#   plt_areasdens[[i]] <- plt_areas[[i]] +
#     geom_ribbon(data = df, aes(x=x, ymin=ymin, ymax=ymax), fill='orange', alpha = 0.2, inherit.aes = F) +
#     geom_line(data = df, aes(x=x, y=y), col='darkorange', linewidth=1.2)
#     
# }
# # Show plot
# plt_areasdens[[1]]; plt_areasdens[[2]]


# log
cat(sprintf("Plots have been generated in file: %s", pdf_out))

# area <- 93
# plot(seq(data_range[area,1], data_range[area,2], length.out = Npoints), estimated_densities[[area]]['est',], type='l')
# hist(data[[area]], probability = T, add = T)
# lines(seq(data_range[area,1], data_range[area,2], length.out = Npoints), estimated_densities[[area]]['est',], lwd=2)
# 

# Map
zoom_bbox <- unname(st_bbox(st_transform(sf_counties[areas,], 4326)))
zoom_map <- sf_ggmap(get_map(zoom_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))
# Labels
names_zoom <- data.frame("PUMA" = areas_names, st_coordinates(st_centroid(st_geometry(sf_counties_3857[areas, ]))))
names_zoom$PUMA <- factor(names_zoom$PUMA, levels = c("U.S.C. & Exposition Park", "Hancock Park & Mid-Wilshire", "West Hollywood & Beverly Hills"))
# Generate
plt_zoom <- ggmap(zoom_map) +
  geom_sf(data = sf_counties_3857, fill=NA, col='gray25', linewidth=0.7, inherit.aes = F) +
  geom_sf(data = sf_counties_3857[areas, ], aes(fill=PUMA), col='gray25', alpha = 0.5, linewidth = 1.4, show.legend = F, inherit.aes = F) +
  geom_sf(data = bound_sf_3857, fill=NA, col='darkred', linewidth=1.4, inherit.aes = F) +
  scale_fill_manual(values = c('steelblue', 'forestgreen', 'darkorange')) +
  geom_label_repel(data = names_zoom, aes(x=X, y=Y, label=PUMA, color=PUMA), inherit.aes = F, show.legend = F) +
  scale_color_manual(values = c('darkorange', 'steelblue', 'forestgreen')) + theme_void()
# Show
x11(height = 4, width = 4); plt_zoom
# Save
# pdf("plots/plt_zoom.pdf", height = 4, width = 4); plt_zoom; dev.off()

# 1 - Info interessanti sui crimini e i boundaries che ho trovato ----

sf_la <- sf_counties[grep("Los Angeles County", sf_counties$Name), ]
sf_la_3857 <- st_transform(sf_la, 3857)
la_bbox <- unname(st_bbox(st_transform(sf_la, 4326)))
# Crop islands
la_bbox[2] <- la_bbox[4] - (la_bbox[3] - la_bbox[1])
la_map <- sf_ggmap(get_map(la_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))

# Select W for LA county
selected_rowcols <- which(sf_counties$Name %in% sf_la$Name)
W_la <- W[selected_rowcols, ]
W_la <- W_la[, selected_rowcols]

# Compute not admissible edges
Eadj_la <- which(W_la == 1, arr.ind = TRUE)
Eb_la <- which(W_la == 0, arr.ind = TRUE)

# Select plinks for LA county
plinks_la <- plinks[selected_rowcols, ]
plinks_la <- plinks_la[, selected_rowcols]

# Compute neighbouring graph for LA
Gn_la <- matrix(0, nrow(plinks_la), ncol(plinks_la))
Gn_la[Eadj_la] <- ifelse(plinks_la[Eadj_la] >= 0.5, 1, 0)

# Compute boundary graph for LA
Gb_la <- matrix(0, nrow(plinks_la), ncol(plinks_la))
Gb_la[Eadj_la] <- ifelse(plinks_la[Eadj_la] < 0.5, 1, 0)

# Compute boundary adjacency list and geometry for LA
bound_list_la <- spdep::mat2listw(Gb_la, style = "B")$neighbours
bound_sf_la <- boundary_geometry(bound_list_la, sf_la)
bound_sf_la_3857 <- st_transform(bound_sf_la, 3857)

# Drop zeros
plinks_la[Eb_la] <- NA
Gn_la[which(Gn_la == 0, arr.ind=T)] <- NA
Gb_la[which(Gb_la == 0, arr.ind=T)] <- NA


df <- read.csv("data/explain_boundaries/2020_LACounty_Crimes_dataset.csv")
sf <- st_as_sf(df, coords = c("LONGITUDE", "LATITUDE")) %>%
  st_set_crs(4326) %>% st_transform(3857) %>%
  mutate(CATEGORY = as.factor(CATEGORY))

# plt_boundaries_mean <- ggplot() +
#   geom_sf(data = sf_counties_3857, aes(fill=post_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
#   scale_fill_gradient(low = 'steelblue', high = 'darkorange') +
#   geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = F) +
#   theme_void() + theme(legend.position = "bottom") + 
#   guides(fill = guide_colourbar(title = "Post. Mean", direction = "horizontal", barwidth = unit(2.5, "in"),
#                                 title.position = "bottom", title.hjust = 0.5))

categories <- c("DISORDERLY CONDUCT", "LIQUOR LAWS", "RECEIVING STOLEN PROPERTY", "VAGRANCY")
sf_crimes <- sf[sf$CATEGORY %in% categories, ]

plt_crimescat <- ggmap(la_map) +
  geom_sf(data = sf_la_3857, col='gray25', linewidth=0.25, fill='white', alpha=0.6, inherit.aes=F) +
  geom_sf(data = sf_crimes, aes(color=CATEGORY), size=0.75, alpha = 0.5, inherit.aes=F) +
  geom_sf(data = bound_sf_la_3857, col='darkred', linewidth=0.5, inherit.aes=F) +
  scale_color_manual(NULL,
                     values = c("steelblue", 'darkorange', 'lightgreen', 'salmon'),
                     labels = c("Disorderly Conduct", "Liquor Laws", "Recieving Stolen Property", "Vagrancy"),
                     guide = bottom_legend(NULL,2)) +
  theme_void() + theme(legend.position = "bottom")
  # scale_color_manual(NULL, values = c("white",'darkred')) + guides(color = guide_legend(override.aes = list(position="none")))
  # guides(color = bottom_legend(NULL, 2)) + theme_void() + theme(legend.position = "bottom")
pdf("plots/plt_CrimesCategories.pdf", height = 4, width = 3.5); plt_crimescat; dev.off()
# x11(height = 4, width = 3.5); plt_crimescat


for (cat in levels(sf$CATEGORY)) {
  print(plt_boundaries_mean + 
    geom_sf(data = sf[sf$CATEGORY == cat, ], inherit.aes = F) +
    geom_sf(data = bound_sf, col='darkred', linewidth = 0.5, inherit.aes = F) +
    ggtitle(cat))
}

plt_gangcrimes <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, col='gray25', linewidth=0.25, fill='white', alpha=0.6, inherit.aes=F) +
  geom_sf(data = sf[sf$GANG_RELATED == "YES", ], aes(col=GANG_RELATED), inherit.aes=F) +
  geom_sf(data = bound_sf, col='darkred', linewidth = 0.5, inherit.aes=F) +
  guides(color = bottom_legend(NULL, 2)) + theme_void() + theme(legend.position = "bottom")
pdf("plots/plt_GangCrimes.pdf", height = 5, width = 4); plt_gangcrimes; dev.off()
# x11(height = 5, width = 4); plt_gangcrimes

# 2 - Percentuale di persone senza assicurazione medica

sf <- read_sf("data/explain_boundaries/Health_Insurance_(census_tract)/Health_Insurance_(census_tract).shp") %>%
  st_transform(3857) %>%
  group_by(csa) %>%
  summarise(uninsured = mean(uninsure_1, na.rm = T))

bottom_colorbar <- function(title, length){
  out <- guide_colorbar(title, direction = "horizontal", barwidth = unit(length, "in"), barheight = unit(0.15, "in"),
                        title.position = "bottom", title.vjust = 0.5, label.vjust = 0.5)
  return(out)
}

plt_insurance <- ggmap(la_map) +
  geom_sf(data = sf_la_3857, col=NA, fill='white', alpha=0.6, inherit.aes=F) +
  geom_sf(data = sf, aes(fill=uninsured), col=NA, inherit.aes=F) +
  geom_sf(data = sf_la_3857, col="gray25", fill=NA, linewidth=0.25, inherit.aes=F) +
  geom_sf(data = bound_sf_la_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
  scale_fill_gradient(low="steelblue", high="darkorange",
                      guide = bottom_colorbar("Population w/o Health Insurance (%)", 2.5)) +
  theme_void() + theme(legend.position = "bottom")
pdf("plots/plt_NoInsurance.pdf", height = 4, width = 3.5); plt_insurance; dev.off()
# x11(height = 4, width = 3.5); plt_insurance
#png("people_without_insurance.png", height = 5, width = 5, units = "in", res = 200); plt_insurance; dev.off()

x11(height = 4, width = 7); gridExtra::grid.arrange(plt_crimescat, plt_insurance, ncol=2)

# 3 - Impatto delle incarcerazioni: numero di arresti della LAPD e della LASD nel 2020

sf <- read_sf("data/Incarceration_Impact_(neighborhood)/Incarceration_Impact__neighborhood_.shp") %>%
  st_zm() %>% st_transform(3857)

plt_arrests <- ggmap(counties_map) +
  geom_sf(data = sf, aes(fill=Combo_Arre), col=NA, inherit.aes=F) +
  geom_sf(data = sf_counties_3857, col = "gray25", fill=NA, linewidth=0.25, inherit.aes=F) +
  geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
  scale_fill_gradient(low="steelblue", high="darkorange",
                      guide = bottom_colorbar("N° of Total Arrests", 3)) +
  theme_void() + theme(legend.position = "bottom")
png("total_arrests.png", height = 5, width = 5, units = "in", res = 200); plt_arrests; dev.off()

# 4 - Potrebbe essere interessante: below Count and Percent below federal poverty level (fpl) and below 200 percent fpl.

sf <- read_sf("data/Below_Poverty_(census_tract)/Below_Poverty__census_tract_.shp") %>%
  st_transform(3857) %>%
  group_by(csa) %>%
  summarise(below_poverty = mean(below_fpl_, na.rm = T),
            below_200poverty = mean(below_20_1, na.rm = T))
  
plt_fpl <- ggmap(counties_map) +
  geom_sf(data = sf, aes(fill=below_poverty), col=NA, inherit.aes=F) +
  geom_sf(data = sf_counties_3857, col = "gray25", linewidth=0.25, fill=NA, inherit.aes=F) +
  geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
  scale_fill_gradient(low="steelblue", high="darkorange",
                      guide = bottom_colorbar("People Below Federal Poverty Level (%)", 3)) +
  theme_void() + theme(legend.position = "bottom")
png("people_below_fpl.png", height = 5, width = 5, units = "in", res = 200); plt_fpl; dev.off()

plt_200fpl <- ggmap(counties_map) +
  geom_sf(data = sf, aes(fill=below_200poverty), col=NA, inherit.aes=F) +
  geom_sf(data = sf_counties_3857, col = "gray25", linewidth=0.25, fill=NA, inherit.aes=F) +
  geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
  scale_fill_gradient(low="steelblue", high="darkorange",
                      guide = bottom_colorbar("People 200% Below Federal Poverty Level (%)", 3)) +
  theme_void() + theme(legend.position = "bottom")
png("people_below_200fpl.png", height = 5, width = 5, units = "in", res = 200); plt_200fpl; dev.off()
