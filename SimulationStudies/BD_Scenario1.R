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
set.seed(230196); data <- list()
for (i in 1:numGroups) {
  if (dist_picker[i] == 1) {
    data[[i]] <- metRology::rt.scaled(n=100, df=6, mean=4, sd=1.5)
  } else {
    data[[i]] <- sn::rsn(n=100, xi=4, omega=1.3, alpha=-3)
    attributes(data[[i]]) <- NULL
  }
}

# Plot - empirical density histograms in areas
# areas <- c(1,4); plt_areas <- list(); k <- 1
# for (i in areas) {
#   # Generate plot
#   tmp <- ggplot() +
#     geom_histogram(data = data.frame("Data" = data[[i]]), aes(x=Data, y = after_stat(density)),
#                    binwidth = function(x) { diff(pretty(x,nclass.Sturges(x)))[1] }, col='steelblue', fill='lightblue') +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     ylab("Density") + ggtitle(sprintf("Area %d", i))
#   # Save plot
#   pdf(sprintf("plt_histArea_%d.pdf", i), height = 4, width = 4); print(tmp); dev.off()
# }


# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 1

# Set sampler parameters
params =
  "
  num_components: 8

  p0_params {
    mu0: 0
    a: 2
    b: 2
    lam_: 0.1
  }

  rho {
    fixed: 0.99
  }

  sigma {
    inv_gamma_prior {
      alpha: 2
      beta: 2
    }
  }

  graph_params {
    beta_prior {
      a: 1
      b: 23
    }
  }
  "

# Sparse inducing prior --> a = 1, b = (2*I - 2) / 3 -1 (see Paci and Consonni (2020))

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

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
  post_vars[,j] <- weights_chain[[j]]^2 %*% vars_chain[[j]]
}

# Add to shapefile dataframe
sf_grid$post_mean <- apply(post_means, 1, mean)
sf_grid$post_var <- apply(post_vars, 1, mean)

# Compute estimated density
data_ranges <- sapply(data, range); Npoints <- 500
estimated_densities <- ComputeDensities(chains, Npoints, data_ranges, alpha = 0.05)

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

# Compute boundary matrix, boundary adj list and geometry
bound_matrix <- W - G_est
bound_list <- spdep::mat2listw(bound_matrix)$neighbours
bound_sf <- boundary_geometry(bound_list, sf_grid)

# Posterior of H - barplot
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plt_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components")
# Show plot
plt_postH

# Posterior of H - Traceplot
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plt_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# Show plot
plt_traceH

# Plot plinks matrix
plinks[which(plinks == 0, arr.ind = T)] <- NA
df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
plt_plinks <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val)) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=numGroups+0.5, ymax=numGroups+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey',
                       guide = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                               barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5)) +
  coord_equal() + theme_void() + theme(legend.position = "bottom")
# Plot estimated graph
df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
plt_Gest <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient(low = 'white', high = 'darkorange',
                      guide = guide_legend(override.aes = list(alpha = 0), title="",
                                           title.position="bottom", label.position="bottom")) +
  theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))
# Show plots
gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2)

# Plot boundaries on the grid + group
plt_boundaries <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=Group), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_manual(values = c("steelblue","darkorange"), labels = c("Student's t", "Skew Normal"),
                    guide = guide_legend(title = "", title.hjust = 0.5, label.position = "bottom",
                                         direction = "horizontal", title.position = "bottom", keywidth = unit(1,"in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom", legend.title = element_text(colour = "transparent"))
# Show plot
plt_boundaries
# Save
pdf("plt_BDgroup.pdf", height = 4, width = 4); plt_boundaries; dev.off()

# Plot boundaries on the grid + post_mean
plt_postmean <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=post_mean), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(title = "Post. Mean", title.hjust = 0.5, title.position = "bottom",
                                             direction = "horizontal", barwidth = unit(2.5, "in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
plt_postmean
# Save
pdf("plt_BDpostMean.pdf", height = 4, width = 4); plt_postmean; dev.off()

# Plot boundaries on the grid + post_var
plt_postvar <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=post_var), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(title = "Post. Variance", title.hjust = 0.5, title.position = "bottom",
                                             direction = "horizontal", barwidth = unit(2.5, "in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
plt_postvar
# Save
pdf("plt_BDpostVar.pdf", height = 4, width = 4); plt_postvar; dev.off()

