# Required libraries
library("SPMIX")
library("sf")
library("ggplot2")

boundary_geometry <- function(boundary_list, sf_geometry) {
  if(!inherits(sf_geometry, "sf")) { stop("'sf_geometry' must be an sf object") }
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }
  
  # Create empty list
  geom_bdd <- list()
  
  for(i in 1:nrow(sf_geometry)) {
    # Get current area and its boundaries
    sel_geom <- sf_geometry[c(i, boundary_list[[i]]), ]
    
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

# Generate shape file from scratch
numGroups <- 36
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Generate distribution picker label
dist_picker <- c(rep(cbind(rep(1,3),rep(2,3)),3), rep(cbind(rep(2,3),rep(1,3)),3))

# Generate sf object
sf_grid <- st_sf(data.frame("Group" = as.factor(dist_picker)), geometry = grid)

# Generate data
set.seed(230196); data <- list(); Ndata = 3000
for (i in 1:numGroups) {
  if (dist_picker[i] == 1) {
    tmp_mean <- rnorm(4, 0.12)
    data[[i]] <- metRology::rt.scaled(n=Ndata, df=6, mean=tmp_mean, sd=1.5)
  } else {
    data[[i]] <- sn::rsn(n=Ndata, xi=4, omega=1.3, alpha=-3)
    attributes(data[[i]]) <- NULL
  }
}

# # Setting MCMC parameters
# burnin = 5000
# niter = 5000
# thin = 1
# 
# # Set sampler parameters
# params =
#   "
#   num_components: 10
# 
#   p0_params {
#     mu0: 0
#     a: 2
#     b: 2
#     lam_: 0.1
#   }
# 
#   rho {
#     fixed: 0.95
#   }
# 
#   sigma {
#     inv_gamma_prior {
#       alpha: 2
#       beta: 2
#     }
#   }
# 
#   graph_params {
#     beta_prior {
#       a: 2
#       b: 36
#     }
#   }
#   "
# 
# # Sparse inducing prior --> a = 1, b = (2*I - 2) / 3 -1 (see Paci and Consonni (2020))
# 
# # Run Spatial sampler
# out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params, type = "rjmcmc")
# if (exists("out")) {
#   filename <- sprintf("BD_Scenario2_chain_%s.dat", format(Sys.time(), format = "%Y%m%d-%H%M"))
#   save(out, file = filename)
# }

# Load chain
load("BD_Scenario2_chain_.dat")

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute posterior mean and variance for each area
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$mean))
vars_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$stdev^2))
weights_chain <- lapply(chains, function(x) t(sapply(x$groupParams, function(y) y$weights)))
post_means <- matrix(nrow = numGroups, ncol=length(chains))
post_vars <- matrix(nrow = numGroups, ncol=length(chains))
for (j in 1:length(chains)) {
  post_means[,j] <- weights_chain[[j]] %*% means_chain[[j]]
  second_moment <- as.vector(weights_chain[[j]] %*% (vars_chain[[j]] + means_chain[[j]]^2))
  post_vars[,j] <- second_moment - post_means[,j]^2
}

# Add to shapefile dataframe
sf_grid$post_mean <- apply(post_means, 1, mean)
sf_grid$post_var <- apply(post_vars, 1, mean)

# Compute estimated density
estimated_densities <- ComputeDensities(chains, seq(range(data)[1],range(data)[2],length.out=500), verbose = T)

# Compute admissible edges
admissible_edges <- which(W != 0, arr.ind = T)

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = T)] <- NA

# Compute median graph
G_est <- matrix(NA, nrow(plinks), ncol(plinks))
G_est[admissible_edges] <- ifelse(plinks[admissible_edges] >= 0.5, 1, NA)

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[admissible_edges] <- as.factor(ifelse(plinks[admissible_edges] < 0.5, 1, NA))

# Compute boundary geometry
if(!all(is.na(Gb))){
  bound_list <- apply(Gb, 1, function(x){which(x == 1)})
  bound_sf <- boundary_geometry(bound_list, sf_grid)
} else {
  cat("No Boundaries have been found\n")
}

# Posterior of H - barplot
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plt_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components")
# Show plot
pdf("plt_postH.pdf", height = 4, width = 4); print(plt_postH); dev.off()

# Posterior of H - Traceplot
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plt_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# Show plot
pdf("plt_traceH.pdf", height = 4, width = 4); print(plt_traceH); dev.off()

# Plot plinks matrix
df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
plt_plinks <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val)) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=numGroups+0.5, ymax=numGroups+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient2(low='steelblue', mid = "white", high = 'darkorange', midpoint = 0.5, na.value = 'white',
                       guide = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                               barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5)) +
  coord_equal() + theme_void() + theme(legend.position = "bottom")

# Plot estimated graph
df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
plt_Gest <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient(low = 'darkorange', high = 'darkorange', na.value = "white",
                      guide = guide_legend(override.aes = list(alpha = 0), title="",
                                           title.position="bottom", label.position="bottom")) +
  theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))
# Show plots
pdf("plt_plinks.pdf", height = 4, width = 4); print(plt_plinks); dev.off()
pdf("plt_Gest.pdf", height = 4, width = 4); print(plt_Gest); dev.off()
# gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2)

# Plot boundaries on the grid + group
plt_boundaries <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=Group), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_manual(values = c("steelblue","darkorange"), labels = c("Student's t", "Skew Normal"),
                    guide = guide_legend(title = NULL, title.hjust = 0.5, label.position = "bottom",
                                         direction = "horizontal", title.position = "bottom", keywidth = unit(1,"in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
pdf("plt_boundaries.pdf", height = 4, width = 4); print(plt_boundaries); dev.off()
# Save
# pdf("plt_BDgroup.pdf", height = 4, width = 4); plt_boundaries; dev.off()

# Plot boundaries on the grid + post_mean
plt_postmean <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=post_mean), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(title = "Post. Mean", title.hjust = 0.5, title.position = "bottom",
                                             direction = "horizontal", barwidth = unit(2.5, "in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
pdf("plt_postmean.pdf", height = 4, width = 4); print(plt_postmean); dev.off()
# Save
# pdf("plt_BDpostMean.pdf", height = 4, width = 4); plt_postmean; dev.off()

# Plot boundaries on the grid + post_var
plt_postvar <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=post_var), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(title = "Post. Variance", title.hjust = 0.5, title.position = "bottom",
                                             direction = "horizontal", barwidth = unit(2.5, "in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
pdf("plt_postvar.pdf", height = 4, width = 4); print(plt_postvar); dev.off()

# # Save
# # pdf("plt_BDpostVar.pdf", height = 4, width = 4); plt_postvar; dev.off()
# 
# # Plot - Empirical density histogram in bordering areas
# areas <- c(3,4)
# plt_areas <- list()
# for (i in 1:length(areas)) {
#   df <- data.frame("x" = data[[areas[i]]])
#   plt_areas[[i]] <- ggplot(data = df, aes(x=x, y = after_stat(density))) +
#     geom_histogram(col=NA, fill='white', bins = 10) +
#     geom_histogram(col='steelblue', fill='steelblue', alpha = 0.4, bins = 10) +
#     xlab("Data") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
# }
# # Show plot
# plt_areas[[1]]; plt_areas[[2]]
# # Save plot
# # pdf("output/plt_areas_empirical_3.pdf", height = 3, width = 3); plt_areas[[1]]; dev.off()
# # pdf("output/plt_areas_empirical_4.pdf", height = 3, width = 3); plt_areas[[2]]; dev.off()
# 
# # Plot - Empirical density histogram + Estimated density
# plt_areasdens <- list()
# for (i in 1:length(areas)) {
#   df <- data.frame("x" = seq(data_ranges[1,areas[i]], data_ranges[2,areas[i]], length.out = Npoints),
#                    "y" = estimated_densities[[areas[i]]]['est', ],
#                    "ymin" = estimated_densities[[areas[i]]]['low', ],
#                    "ymax" = estimated_densities[[areas[i]]]['up', ])
#   plt_areasdens[[i]] <- plt_areas[[i]] +
#     geom_ribbon(data = df, aes(x=x, ymin=ymin, ymax=ymax), fill='orange', alpha = 0.2, inherit.aes = F) +
#     geom_line(data = df, aes(x=x, y=y), col='darkorange', linewidth=1.2) +
#     ylim(c(0,0.6))
#   
# }
# # Show plot
# plt_areasdens[[1]]; plt_areasdens[[2]]
# # Save plot
# # pdf("output/plt_areas_est_3.pdf", height = 3, width = 3); plt_areasdens[[1]]; dev.off()
# # pdf("output/plt_areas_est_4.pdf", height = 3, width = 3); plt_areasdens[[2]]; dev.off()
# 
# 
# # ESS, WAIC, robe per review ----------------------------------------------
# 
# # PLOT - Traceplot of p
# p_chain <- sapply(chains, function(x){x$p})
# df <- data.frame("Iter"=1:length(p_chain), "Value"=p_chain)
# # Generate plot
# plt_p_chain <- ggplot() +
#   geom_line(data = df, aes(x = Iter, y = Value)) +
#   xlab('Iteration') + ylab('p')
# # Show / Save
# # x11(height = 4, width = 4); plt_p_chain
# pdf("plots/plt_p_chain.pdf", height = 4, width = 4); plt_p_chain; dev.off()
# 
# # ESS of p
# mcmcse::ess(p_chain)
# 
# # PLOT - Traceplot of |G|
# Nedge_chain <- sapply(G_chain, function(x){sum(x[upper.tri(x)])})
# df <- data.frame("Iter"=1:length(Nedge_chain), "Value"=Nedge_chain)
# # Generate plot
# plt_Nedge_chain <- ggplot() +
#   geom_line(data = df, aes(x = Iter, y = Value)) +
#   xlab('Iteration') + ylab('|G|')
# # Show / Save
# # x11(height = 4, width = 4); plt_Nedge_chain
# pdf("plots/plt_Nedge_chain.pdf", height = 4, width = 4); plt_Nedge_chain; dev.off()
# 
# # ESS of |G|
# mcmcse::ess(Nedge_chain)
# 
# # PLOT - Traceplot of sigma^2
# sigma_chain <- sapply(chains, function(x){x$Sigma$data[1]})
# df <- data.frame("Iter"=1:length(sigma_chain), "Value"=sigma_chain)
# # Generate plot
# plt_sigma_chain <- ggplot() +
#   geom_line(data = df, aes(x = Iter, y = Value)) +
#   xlab('Iteration') + ylab(bquote(sigma^2))
# # Show / Save
# # x11(height = 4, width = 4); plt_sigma_chain
# pdf("plots/plt_sigma_chain.pdf", height = 4, width = 4); plt_sigma_chain; dev.off()
# 
# # ESS of sigma^2
# mcmcse::ess(sigma_chain)
# 
# # TODO: improve doc for WAIC
# tmp <- ComputePosteriorLPDF(data, chains, verbose = T)
# 
# # Plot - |G| with different initializations (DA SISTEMARE)
# load("output/Nedge_chain_startempty.dat")
# df <- data.frame("Iter" = 1:length(Nedge_chain), "Value"=Nedge_chain, "Group"=rep("Empty", length(Nedge_chain)))
# load("output/Nedge_chain_startfull.dat")
# df <- rbind(df, data.frame("Iter" = 1:length(Nedge_chain), "Value"=Nedge_chain, "Group"=rep("Full", length(Nedge_chain))))
# load("output/Nedge_chain_startrandom.dat")
# df <- rbind(df, data.frame("Iter" = 1:length(Nedge_chain), "Value"=Nedge_chain, "Group"=rep("Random", length(Nedge_chain))))
# #Generate plot
# plt_Nedge_chain_diffInit <- ggplot() +
#   geom_line(data=df, aes(x=Iter, y=Value, color=Group)) +
#   scale_color_manual(values = c("Full"="darkorange", "Empty"='steelblue', "Random"='forestgreen')) +
#   xlab('Iteration') + ylab('|G|') +
#   guides(color = guide_legend("Initial Value", direction = 'horizontal', position = 'bottom'))
# # Show / Save
# # x11(height = 4, width = 4); plt_Nedge_chain_diffInit
# pdf("plots/plt_Nedge_chain_diffInit.pdf", height = 4, width = 4); plt_Nedge_chain_diffInit; dev.off()
