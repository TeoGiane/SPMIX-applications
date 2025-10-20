## Second test for the Reversible Jump Sampler ##
# The scenario is described in Sec. 6.2 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")
library("sf")

###########################################################################
# Data Generation ---------------------------------------------------------

# Dimensions
numGroups <- 9
numComponents <- 3

# Generate shape file from scratch
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))

# Generate true weights in each area
centers <- st_coordinates(st_centroid(grid))
transformed_weights <- matrix(0, nrow = numGroups, ncol = numComponents-1)
transformed_weights[,1] <- + 3*(centers[,1] - 0.5) + 3*(centers[,2] - 0.5)
transformed_weights[,2] <- - 3*(centers[,1] - 0.5) - 3*(centers[,2] - 0.5)
weights <- t(apply(transformed_weights, 1, InvAlr))

# Generate sf object
sf_grid <- st_sf(data.frame(weights), geometry = grid)

# Compute adjacency matrix
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Generate data
set.seed(230196)
means <- c(-5,0,5); sds <- c(1,1,1); Ns <- rep(100, numGroups)
data <- list()
for (i in 1:numGroups) {
  cluster_alloc <- sample(1:numComponents, prob = weights[i,], size = Ns[i], replace = T)
  data[[i]] <- rnorm(Ns[i], mean = means[cluster_alloc], sd = sds[cluster_alloc])
}

###########################################################################

###########################################################################
# Sampler Execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Default parameters
params =
  "
  num_components: 5

  p0_params {
    mu0: 0
    a: 2
    b: 2
    lam_: 0.1
  }

  rho {
    beta_prior {
      a: 2
      b: 2
    }
  }

  sigma {
    inv_gamma_prior {
      alpha: 2
      beta: 2
    }
  }
  "

# Run Spatial sampler
out <- Sampler.DensityEstimation(burnin, niter, thin, data, W, params, type = "rjmcmc")

###########################################################################

###########################################################################
# Posterior analysis ------------------------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$mean))
stdev_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$stdev))

# Computing estimated densities
estimated_densities <- ComputeDensities(chains, seq(range(data)[1], range(data)[2], length.out=500), verbose = T)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:numGroups) {
  x_grid <- seq(data_ranges[1,i], data_ranges[2,i], length.out = Npoints)
  xgrid_expand <- t(rbind(replicate(numComponents, x_grid, simplify = "matrix")))
  true_dens <- t(as.matrix(weights[i,])) %*% dnorm(xgrid_expand, means, sds)
  true_densities[[i]] <- true_dens
}

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Weights distribution on the spatial grid - 1st component
df_text <- data.frame("Name" = row.names(sf_grid),
                      st_coordinates(st_centroid(st_geometry(sf_grid))))
plt_w1 <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=X1), col=NA) +
  geom_text(data = df_text, aes(x=X,y=Y,label=Name), col='gray25', size = 10) +
  geom_rect(aes(xmin=0, ymin=0, xmax=1, ymax=1), fill=NA, col='gray25', linewidth=1) +
  scale_fill_gradient(low="steelblue", high="orange") +
  theme_void() + theme(legend.position = "none")
# Show plot
plt_w1
# Save
pdf("plots/DE_Scenario2-plt_w1.pdf", height = 4, width = 4); print(plt_w1); dev.off()

# Weights distribution on the spatial grid - 2nd component
plt_w2 <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=X2), col=NA) +
  # geom_text(data = df_text[c(1,7),], aes(x=X,y=Y,label=Name), col='gray25', size = 10) +
  geom_rect(aes(xmin=0, ymin=0, xmax=1, ymax=1), fill=NA, col='gray25', linewidth=1) +
  scale_fill_gradient(low="steelblue", high="orange") +
  theme_void() + theme(legend.position = "none")
# Show plot
plt_w2
# Save
pdf("plots/DE_Scenario2-plt_w2.pdf", height = 4, width = 4); print(plt_w2); dev.off()

# Adjacency matrix
df <- reshape2::melt(W, c("x", "y"), value.name = "value")
df_text <- data.frame("Name" = 1:numGroups, "X" = df$x[1:numGroups], "Y" = df$y[1:numGroups])
plt_W <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=value)) +
  geom_text(data = df_text, aes(x=X,y=Y,label=Name), col='gray25', size=5) +
  geom_text(data = df_text, aes(x=Y,y=X,label=Name), col='gray25', size=5) +
  scale_fill_gradient(low="white", high="orange") +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(W)+0.5, ymax=ncol(W)+0.5),
            fill=NA, col='gray25', linewidth=1) +
  theme_void() + coord_equal() + theme(legend.position = "none")
# Show plot
plt_W
# Save
pdf("plots/DE_Scenario2-plt_G.pdf", height = 4, width = 4); print(plt_W); dev.off()

# Posterior of H - Barplot
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plt_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") + xlab("N° of Components") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# Show plot
plt_postH
# Save
# pdf("DE_Scenario2/plt_postH.pdf", height = 4, width = 4); print(plt_postH); dev.off()

# Posterior of H - Traceplot
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plt_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# Show plot
plt_traceH
# Save
# pdf("DE_Scenario2/plt_traceH.pdf", height = 4, width = 4); print(plt_traceH); dev.off()

# Comparison plots between estimated and true densities in i-th area + bands
plts_area <- list()
for (i in 1:numGroups) {
  # Auxiliary dataframe
  df <- data.frame('grid'=seq(data_ranges[1,i], data_ranges[2,i], length.out=Npoints),
                   t(estimated_densities[[i]]), 'true'=t(true_densities[[i]]))
  # Generate plot
  tmp <- ggplot(data = df, aes(x=grid)) +
    geom_line(aes(y=est, color="Estimated"), linewidth = 1) +
    geom_line(aes(y=true, color="True"), linewidth = 1) +
    scale_color_manual(breaks=c("Estimated","True"), values=c("darkorange", "steelblue"),
                       guide = guide_none()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Grid") + ylab("Density") + ggtitle(paste0("Area ", i))
  # Add credibility band if present
  if (dim(estimated_densities[[i]])[1] > 1)
    tmp <- tmp + geom_ribbon(aes(ymax=up, ymin=low), fill="orange", alpha=0.3)
  # Save plot and clean useless variables
  plts_area[[i]] <- tmp; rm(list=c('df','tmp'))
}

# Printing plots
gridExtra::grid.arrange(grobs=plts_area, nrow=3, ncol=3)
# Save
# pdf("DE_Scenario2/plt_densCompare.pdf", height = 8, width = 8); gridExtra::grid.arrange(grobs=plts_area, nrow=3, ncol=3); dev.off()
# pdf("DE_Scenario2/plt_area147.pdf", height = 4, width = 12); gridExtra::grid.arrange(plts_area[[1]], plts_area[[4]], plts_area[[7]], ncol=3); dev.off()

###########################################################################
