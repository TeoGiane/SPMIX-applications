## Second test for the Reversible Jump Sampler ##
# The scenario is described in Sec. 6.2 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")
library("sp") # also rgdal, rgeos, maptools have been installed

# Helper functions
BuildLattice <- function(coords) {
  numGroups <- dim(coords)[1]

  # Build spatial grid from coordinates
  spPixels <- sp::SpatialPixels(sp::SpatialPoints(coords))
  spatialGrid <- sp::as.SpatialPolygons.GridTopology(spPixels@grid)

  # Compute proximity matrix
  nblist <- sp::gridIndex2nb(spPixels, maxdist = 1)
  W <- matrix(0,numGroups,numGroups)
  for (i in 1:numGroups) {
    row_ind <- i; col_inds <- unlist(nblist[[i]]); names(col_inds) <- NULL
    W[row_ind,col_inds] <- 1; W[i,i] <- 0
  }

  # Return output list
  out <- list(); out[[1]] <- spatialGrid; out[[2]] <- W
  names(out) <- c("spatialGrid", "proxMatrix")
  return(out)
}

###########################################################################
# Data Generation ---------------------------------------------------------

# Build spatial grid and proximity matrix from coordinates
numGroups <- 9; numComponents <- 3
coords <- expand.grid(x = seq(0,1,length.out = sqrt(numGroups)),
                      y = seq(0,1,length.out = sqrt(numGroups)))
xbar <- ybar <- 0.5

spatial_grid <- BuildLattice(coords)$spatialGrid
W <- BuildLattice(coords)$proxMatrix
labels <- row.names(coordinates(spatial_grid))
rm(list='coords')

# Generating weights matrices
trasformed_weights <- matrix(0,nrow=numGroups,ncol=numComponents-1)
weights <- matrix(0,nrow=numGroups,ncol=numComponents)
for (i in 1:numGroups) {
  center <- coordinates(spatial_grid)[i,]
  trasformed_weights[i,1] <- 3*(center[1]-xbar)+3*(center[2]-ybar)
  trasformed_weights[i,2] <- -3*(center[1]-xbar)-3*(center[2]-ybar)
  weights[i,] <- InvAlr(trasformed_weights[i,])
}
rm(list=c('i','center'))
row.names(trasformed_weights) <- labels
row.names(weights) <- labels

# Generate data
set.seed(230196)
means <- c(-5,0,5); sds <- c(1,1,1); Ns <- rep(25, numGroups)
data <- list()
for (i in 1:numGroups) {
  cluster_alloc <- sample(1:numComponents, prob = weights[i,], size = Ns[i], replace = T)
  data[[i]] <- rnorm(Ns[i], mean = means[cluster_alloc], sd = sds[cluster_alloc])
}
names(data) <- labels
rm(list=c('cluster_alloc','i'))

###########################################################################

###########################################################################
# Sampler Execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 2000
niter = 2000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- Sampler.DensityEstimation(burnin,niter,thin,data,W,params_filename,type = "rjmcmc")

###########################################################################

###########################################################################
# Sampler Execution (full run) --------------------------------------------

# Setting MCMC parameters
burnin = 0
niter = 10000
thin = 1

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- Sampler.DensityEstimation(burnin,niter,thin,data,W,params_filename,type = "rjmcmc")

###########################################################################

###########################################################################
# Posterior analysis ------------------------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$mean))
stdev_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$stdev))

# Computing estimated densities
data_ranges <- sapply(data, range); Npoints <- 500
estimated_densities <- ComputeDensities(chains, Npoints, data_ranges, alpha=0.05, names=labels)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:numGroups) {
  x_grid <- seq(data_ranges[1,i], data_ranges[2,i], length.out = Npoints)
  xgrid_expand <- t(rbind(replicate(numComponents, x_grid, simplify = "matrix")))
  true_dens <- t(as.matrix(weights[i,])) %*% dnorm(xgrid_expand, means, sds)
  true_densities[[i]] <- true_dens
}
names(true_densities) <- labels

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Weights distribution on the spatial grid - Plot
df <- fortify(spatial_grid)
for (i in 1:length(labels)) {
  df[which(df$id==labels[i]),'w_i1'] <- weights[i,1]
  df[which(df$id==labels[i]),'w_i2'] <- weights[i,2]
}
plot_weightsi1 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i1)) +
  geom_polygon() + scale_fill_gradient(low="steelblue", high="orange") +
  xlab("") + ylab("") + theme(legend.position = "none")
plot_weightsi2 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i2)) +
  geom_polygon() + scale_fill_gradient(low="steelblue", high="orange") +
  xlab("") + ylab("") + theme(legend.position = "none")
rm(list='df')

# Adjacency matrix - plot
df <- reshape2::melt(W)
plot_W <- ggplot(data=df, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + scale_fill_gradient(low="white", high="orange") +
  xlab("") + ylab("") + theme(legend.position = "none")
rm(list='df')

# Posterior of H - Barplot
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") + ylab("N° of Components") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #+
  #ggtitle("Posterior of H")
rm(list='df')

# Posterior of H - Traceplot
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plot_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #+
  #ggtitle("Traceplot of H")
rm(list='df')

# Comparison plots between estimated and true densities in i-th area + bands
plots_area <- list()
for (i in 1:numGroups) {
  # Auxiliary dataframe
  df <- data.frame('grid'=seq(data_ranges[1,i], data_ranges[2,i], length.out=Npoints),
                   t(estimated_densities[[i]]),
                   'true'=t(true_densities[[i]]))
  # Generate plot
  tmp <- ggplot(data = df, aes(x=grid)) +
    geom_line(aes(y=est, color="Estimated"), lwd=1) +
    geom_line(aes(y=true, color="True"), lwd=1) +
    scale_color_manual("", breaks=c("Estimated","True"), values=c("Estimated"="darkorange", "True"="steelblue")) +
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none") +
    xlab("Grid") + ylab("Density") + ggtitle(paste0("Area ", i))
  # Add credibility band if present
  if (dim(estimated_densities[[i]])[1] > 1)
    tmp <- tmp + geom_ribbon(aes(ymax=up, ymin=low), fill="orange", alpha=0.3)
  # Save plot and clean useless variables
  plots_area[[i]] <- tmp; rm(list=c('df','tmp'))
}
names(plots_area) <- labels

# Printing plots
x11(height = 4, width = 4); plot_weightsi1
x11(height = 4, width = 4); plot_weightsi2
x11(height = 4, width = 4); plot_W
x11(height = 4, width = 4); plot_postH
x11(height = 4, width = 4); plot_traceH
x11(height = 8.27, width = 8.27); gridExtra::grid.arrange(grobs=plots_area, nrow=3, ncol=3)

###########################################################################
