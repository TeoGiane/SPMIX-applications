## First test for the Reversible Jump Sampler ##

# Required libraries
library("SPMIX")
library("ggplot2")

###########################################################################
# Data Simulation ---------------------------------------------------------

# Generate data (1 location, mixture of 3 normals)
set.seed(230196)
ngroups <- 1; ncomponents <- 4; N <- 1000
means <- c(-4,1,5,8); std_devs <- c(1,1,1,1); weights <- matrix(c(2/8,3/8,1/8,2/8),ngroups,ncomponents,T)
cluster_alloc <- sample(1:ncomponents, prob = weights, size = N, replace = T)
data <- list(); data[[1]] <- rnorm(N, mean = means[cluster_alloc], sd = std_devs[cluster_alloc])

# Generate Matrix W
W <- matrix(0, nrow = 1, ncol = 1, byrow = T)

###########################################################################

###########################################################################
# Sampler Run -------------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- Sampler.DensityEstimation(burnin, niter, thin, data, W, params_filename, type="rjmcmc")

###########################################################################

###########################################################################
# Analysis and Visualization ----------------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)

# Barplot of the estimated posterior for H
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components") #+ ggtitle("Posterior of H")
rm(list='df')

x11(height = 4, width = 4); plot_postH

# Plotting the average density over iterations and compare with true curve
# Compute estimated density
data_ranges <- sapply(data, range); Npoints <- 500
estimated_densities <- ComputeDensities(chains, Npoints, data_ranges, alpha = 0.05)

# Compute true densities
true_densities <- list()
for (i in 1:ngroups) {
  x_grid <- seq(data_ranges[1,i], data_ranges[2,i], length.out = Npoints)
  xgrid_expand <- t(rbind(replicate(ncomponents, x_grid, simplify = "matrix")))
  true_dens <- t(as.matrix(weights[i,])) %*% dnorm(xgrid_expand, means, std_devs)
  true_densities[[i]] <- true_dens
}

# Density comparison - Plot
# Auxiliary dataframe
df <- data.frame('grid'=seq(data_ranges[1,1], data_ranges[2,1], length.out=Npoints),
                 t(estimated_densities[[1]]),
                 'true'=t(true_densities[[1]]))
# Generate plot
plot_densCompare <- ggplot(data = df, aes(x=grid)) +
  geom_line(aes(y=est, color="Estimated"), linewidth = 1) +
  geom_line(aes(y=true, color="True"), linewidth = 1) +
  scale_color_manual("", breaks=c("Estimated","True"), values=c("Estimated"="darkorange", "True"="steelblue")) +
  theme(plot.title = element_text(face="bold", hjust = 0.5)) +
  theme(legend.title=element_blank(), legend.position="bottom") +
  xlab("Grid") + ylab("Density") + ggtitle(paste0("Area ", i))
# Add credibility band if present
if (dim(estimated_densities[[i]])[1] > 1) {
  plot_densCompare <- plot_densCompare + geom_ribbon(aes(ymax=up, ymin=low), fill="orange", alpha=0.3)
}
# Clean useless variables
rm(list='df')

x11(width = 6, height = 4); plot_densCompare

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
out <- Sampler.DensityEstimation(burnin, niter, thin, data, W, params_filename, type = "rjmcmc")

###########################################################################

###########################################################################
# Analyses and Visualization (full run) -----------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)

# Traceplot for the whole chain (no burnin, no thinning)
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plot_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.2) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
rm(list='df')

x11(height = 4, width = 4); plot_traceH

###########################################################################

###########################################################################
# Sensitivity w.r.t. number of observations -------------------------------

# Clean console and set number of observations
rm(list=ls()); set.seed(230196)
Ns <- c(100,250,500,1000)
W <- matrix(0, nrow = 1, ncol = 1, byrow = T)
plots <- list()

# Sampler run loop
for (i in 1:length(Ns)) {
  # Generate data (1 location, mixture of 3 normals)
  ngroups <- 1; ncomponents <- 3; N <- Ns[i]
  means <- c(-4,1,5); std_devs <- c(1,1,1); weights <- c(1/6,3/6,2/6)
  cluster_alloc <- sample(1:ncomponents, prob = weights, size = N, replace = T)
  data <- list(); data[[1]] <- rnorm(N, mean = means[cluster_alloc], sd = std_devs[cluster_alloc])

  # Grab input filenames
  params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

  # Run Spatial sampler
  out <- Sampler.DensityEstimation(0, 10000, 1, data, W, params_filename, type = "rjmcmc")

  # Deserialization
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  H_chain <- sapply(chains, function(x) x$num_components)

  # Traceplot for the whole chain (no burnin, no thinning)
  df <- data.frame("Iteration"=1:10000, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
  tmp <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
    ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(lwd=0.1) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    ggtitle(paste0("Traceplot of H - ", Ns[i], " obs."))
  rm(list='df')

  # Save plot
  plots[[i]] <- tmp; rm(list='tmp')
}

# Visualization
x11(height = 6, width = 6); gridExtra::grid.arrange(grobs=plots, nrow=2, ncol=2)

###########################################################################
