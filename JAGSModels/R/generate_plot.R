library("ggplot2")
library("rjags")
library("sf")
library("spdep")

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

boundary_geometry <- function(boundary_graph, sf_geometry) {
  # Check
  if(!inherits(sf_geometry, "sf")) { stop("'sf_geometry' must be an sf object") }
  
  # Add id column
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }
  
  # Generate boundary list using upper-triangular view of boundary_graph
  boundary_graph[lower.tri(boundary_graph)] <- NA
  boundary_list <- apply(boundary_graph, 1, function(x){which(x==1)})
  
  # Create empty list
  geom_bdd <- list(); k <- 1
  for(i in 1:nrow(sf_geometry)) {
    for (j in (boundary_list[[i]])) {
      # Compute boundary geometry between i and j
      sel_geom <- sf_geometry[c(i,j), c("id", "geometry")]
      bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
      bounds <- st_geometry(bounds[bounds$id != bounds$id.1, ])
      # Append to list
      geom_bdd[[k]] <- st_sf(geometry = bounds); k <- k+1
    }
    
    # Get current area and its boundaries
    # sel_geom <- sf_geometry[c(i, boundary_list[[i]]), c("id", "geometry")]
    
    # Compute geometry of boundary
    # bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
    # bounds <- st_geometry(bounds[bounds$id != bounds$id.1, ])
    
    # Add to list
    # geom_bdd[[i]] <- st_sf(geometry = bounds)
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


# CODE --------------------------------------------------------------------

# Import shapefile
sf_counties <- read_sf("data/counties-pumas/counties-pumas.shp")

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Import data
load("data/data_001.dat")

# Set beta parameters
abeta <- 2
bbeta <- length(data)

# Import chain
filename <- file.path(getwd(), "output", sprintf("MCAR_fit-a%gb%g.dat", abeta, bbeta))
load(filename)

# 100 Iterations of burnin
burnin <- 100
p_chain <- lapply(as.mcmc.list(MCAR_fit$p), function(chain){ as.vector(chain[(burnin+1):dim(chain)[1],1:dim(chain)[2]]) })
p_chain <- do.call(c, p_chain)
sigmasq_chain <- lapply(as.mcmc.list(MCAR_fit$one_over_sigmasq), function(chain){ as.vector(chain[(burnin+1):dim(chain)[1],1:dim(chain)[2]]) })
sigmasq_chain <- 1 / do.call(c, sigmasq_chain)
tausq_chain <- lapply(as.mcmc.list(MCAR_fit$one_over_tausq), function(chain){ as.vector(chain[(burnin+1):dim(chain)[1],1:dim(chain)[2]]) })
tausq_chain <- 1 / do.call(c, tausq_chain)
psi_chain <- lapply(as.mcmc.list(MCAR_fit$psi), function(chain){ as.matrix(chain[(burnin+1):dim(chain)[1],1:dim(chain)[2]]) })
psi_chain <- do.call(rbind, psi_chain)

# Get posterior median in each area
post_quantiles <- matrix(nrow = length(data), ncol = 5)
for (q in 1:5) {
  post_quantiles[,q] <- colMeans(psi_chain[,(length(data)*(q-1)+1):(length(data)*q)])
}

# Compute admissible edges
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute PPI matrix (remove not admissible edges)
plinks <- summary(MCAR_fit$G, FUN = mean)$stat # + abs(jitter(matrix(0,length(data),length(data))))
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Threshold -> posterior median of p
gamma_graph <- 0.5

# Compute neighbouring graph
Gn <- matrix(NA, nrow(plinks), ncol(plinks))
Gn[Eadj] <- ifelse(plinks[Eadj] >= gamma_graph, 1, NA)

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < gamma_graph, 1, NA)

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

# PLOT - Traceplot of p
p_df <- data.frame("Iteration"=1:length(p_chain), "p"=p_chain)
# Generate
plt_p <- ggplot() +
  geom_line(data = p_df, aes(x=Iteration, y=p), linewidth=1.2)

# PLOT - traceplot of sigma^2
sigma2_df <- data.frame("Iteration"=1:length(sigmasq_chain), "sigma2"=sigmasq_chain)
# Generate
plt_sigma2 <- ggplot() +
  geom_line(data = sigma2_df, aes(x=Iteration, y=sigma2), linewidth=1.2) +
  ylab(bquote(sigma^2))

# PLOT - traceplot of tau^2
tau2_df <- data.frame("Iteration"=1:length(tausq_chain), "tau2"=tausq_chain)
# Generate
plt_tau2 <- ggplot() +
  geom_line(data = tau2_df, aes(x=Iteration, y=tau2), linewidth=1.2) +
  ylab(bquote(tau^2))

# Show all plots together
x11(width = 12, height = 8); gridExtra::grid.arrange(plt_plinks, plt_boundaries_median, plt_L1glob, plt_p, plt_sigma2, plt_tau2, ncol = 3, nrow = 2)
# Save all plots together
filename <- file.path(getwd(),"plots",sprintf("plt_BDSumStats-a%gb%g.pdf", abeta, bbeta))
pdf(filename, width = 12, height = 4); gridExtra::grid.arrange(plt_plinks, plt_p, plt_sigma2, ncol = 3); dev.off()

# quick plot of psis
boxplot(psi_chain); abline(v = 93*(1:4), col='red', lwd=2, lty=1)

# Show
# x11(height = 4, width = 4); plt_plinks
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_plinks-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_plinks; dev.off()

# Show
# x11(height = 4, width = 4); plt_p
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_p-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_p; dev.off()

# Show
# x11(height = 4, width = 4); plt_sigma2
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_sigma2-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_sigma2; dev.off()

# Show
x11(height = 4, width = 4); plt_tau2
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_tau2-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_tau2; dev.off()


# Compute boundary geometry
sf_boundary <- boundary_geometry(Gb, sf_counties)

# Get maps from stadia
counties_bbox <- unname(st_bbox(st_transform(sf_counties, 4326)))
counties_map <- sf_ggmap(get_map(counties_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))

# PLOT - Posterior mean heatmap + detected boundaries
sf_counties$post_median <- post_quantiles[,3]
sf_counties_3857 <- st_transform(sf_counties, 3857)
sf_boundary_3857 <- st_transform(sf_boundary, 3857)
# Generate
plt_boundaries_median <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_median), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange',
                      guide = guide_colorbar("Post. Median", direction = "horizontal", barwidth=unit(3,"in"),
                                             title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  geom_sf(data = sf_boundary_3857, col='darkred', linewidth = 0.3, inherit.aes = FALSE) +
  theme_void() + theme(legend.position = "bottom")
# Show / Save
x11(height = 4, width = 4); plt_boundaries_median


# L1 distance between areas - Global Comparison ---------------------------

# Neighbouring pairs
Gn_up <- Gn; Gn_up[lower.tri(Gn_up)] <- NA
neigh_pairs <- which(Gn_up == 1, arr.ind = T)

# Boundary pairs
Gb_up <- Gb; Gb_up[lower.tri(Gb_up)] <- NA
bound_pairs <- which(Gb_up == 1, arr.ind = T)

# Compute L1 distance between all neighbouring pairs (GLOBAL)
L1_neigh <- data.frame("Type" = rep("Neigh", nrow(neigh_pairs)),
                       "Dist" = rep(NA, nrow(neigh_pairs)))
for (i in 1:nrow(neigh_pairs)) {
  L1_neigh$Dist[i] <- sum(abs(post_quantiles[neigh_pairs[i,1],] - post_quantiles[neigh_pairs[i,2],]))
}

# Compute L1 distance between all boundary pairs (GLOBAL)
L1_bound <- data.frame("Type" = rep("Bound", nrow(bound_pairs)),
                       "Dist" = rep(NA, nrow(bound_pairs)))
for (i in 1:nrow(bound_pairs)) {
  L1_bound$Dist[i] <- sum(abs(post_quantiles[bound_pairs[i,1],] - post_quantiles[bound_pairs[i,2],]))
}

# Visualization - boxplot comparison
L1glob <- rbind(L1_neigh, L1_bound); rm(L1_neigh, L1_bound)
L1glob$Type <- as.factor(L1glob$Type)
# Generate plot
plt_L1glob <- ggplot() +
  geom_boxplot(data = L1glob, aes(x=Type, y=Dist)) +
  geom_boxplot(data = L1glob, aes(x=Type, y=Dist, fill=Type, color=Type), staplewidth = 0.3, alpha = 0.3, show.legend = F) +
  scale_x_discrete(labels = c(bquote(d[FM]), bquote(d[TM]))) + labs(x=NULL,y=NULL) +
  scale_fill_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  scale_color_manual(values = c("Neigh"="gray25", "Bound"="darkred")) +
  theme(text = element_text(size = 14))
# Show / Save
x11(height = 3, width = 4); plt_L1glob
