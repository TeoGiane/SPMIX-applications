## Test for boundary detection for the Reversible Jump Sampler ##
# The scenario is described in Sec. 6.1 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")

###########################################################################
# Data Generation - Scenario I --------------------------------------------

# Generating data
set.seed(230196)
I <- 6; Ns <- rep(100,I); means <- c(rep(-5,2),rep(0,2),rep(5,2))
data <- list(); labels <- sprintf("g%d", 1:I)
for (i in 1:I) {
  data[[i]] <- rnorm(Ns[i], means[i], 1)
}
rm(list = 'i')
names(data) <- labels

# Setting initial W
W <- matrix(1,I,I)

###########################################################################

###########################################################################
# Data Generation - Scenario II -------------------------------------------

# Generating data
set.seed(230196)
I <- 6; Ns <- rep(100,I)
data <- list(); labels <- sprintf("g%d", 1:I)
for (i in 1:I) {
  if (i %in% c(1,2)){
    data[[i]] <- metRology::rt.scaled(n=Ns[i], df=6, mean=-4, sd=1)
  }
  if (i %in% c(3,4)){
    data[[i]] <- sn::rsn(n=Ns[i], xi=4, omega=4, alpha=1); attributes(data[[i]]) <- NULL
  }
  if (i %in% c(5,6)) {
    data[[i]] <- rchisq(n=Ns[i], df=3)
  }
}
rm(list = 'i')
names(data) <- labels

# Setting initial W
W <- matrix(1,I,I)

###########################################################################

###########################################################################
# Sampler execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params_filename)

###########################################################################

###########################################################################
# Posterior Analysis - Scenario I -----------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

# Computing estimated densities
data_ranges <- sapply(data, range); Npoints <- 500
estimated_densities <- ComputeDensities(chains, Npoints, data_ranges, names=labels)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:I) {
  x <- seq(data_ranges[1,i], data_ranges[2,i], length.out=Npoints)
  true_densities[[i]] <- t(dnorm(x, means[i], 1))
}
rm(list = c('x','i'))
names(true_densities) <- labels

###########################################################################

###########################################################################
# Posterior Analysis - Scenario II ----------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

# Computing estimated densities
data_ranges <- sapply(data, range); Npoints <- 500
estimated_densities <- ComputeDensities(chains, Npoints, data_ranges, names=labels)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:I) {
  x <- seq(data_ranges[1,i], data_ranges[2,i], length.out=Npoints)
  if (i %in% c(1,2)){
    true_densities[[i]] <- t(metRology::dt.scaled(x, df=6, mean=-4, sd=1))
  }
  if (i %in% c(3,4)){
    true_densities[[i]] <- t(sn::dsn(x, xi=4, omega=4, alpha=1))
  }
  if (i %in% c(5,6)) {
    true_densities[[i]] <- t(dchisq(x, df=3))
  }
}
rm(list = c('x','i'))
names(true_densities) <- labels

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Comparison plots between estimated and true densities in i-th area
plots_area <- list()
for (i in 1:I) {
  # Auxiliary dataframe
  df <- data.frame('grid'=seq(data_ranges[1,i], data_ranges[2,i], length.out=Npoints),
                   t(estimated_densities[[i]]),
                   'true'=t(true_densities[[i]]))
  # Generate plot
  tmp <- ggplot(data = df, aes(x=grid)) +
    geom_line(aes(y=est, color="Estimated"), linewidth=1) +
    geom_line(aes(y=true, color="True"), linewidth=1) +
    scale_color_manual("", breaks=c("Estimated","True"), values=c("Estimated"="orange", "True"="steelblue")) +
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none") +
    xlab("Grid") + ylab("Density") + ggtitle(paste0("Area ", i))
  # Add credibility band if present
  if (dim(estimated_densities[[i]])[1] > 1)
    tmp <- tmp + geom_ribbon(aes(ymax=up, ymin=low), fill="orange", alpha=0.3)
  # Save plot and clean useless variables
  plots_area[[i]] <- tmp; rm(list=c('df','tmp'))
}

# Plot plinks matrix
plinks[which(plinks == 0, arr.ind = T)] <- NA
df <- reshape2::melt(plinks); names(df) <- c("x", "y", "val")

plt_plinks <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val)) +
  scale_fill_gradient(low='steelblue', high = 'orange', na.value = 'lightgrey') +
  coord_equal() + theme_void() + theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

# Visualization
x11(height = 6, width = 8.27); gridExtra::grid.arrange(grobs=plots_area, nrow=2, ncol=3)
x11(height = 4, width = 4); plt_plinks

# Save pdf
dev.copy2pdf(device=x11, file="BD_fixed_p05_plinks.pdf"); dev.off()
dev.copy2pdf(device=x11, file="BD_fixed_p05_Densities.pdf"); dev.off()

###########################################################################
