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

# Set sampler parameters
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
    b: 5
  }
}
"

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

###########################################################################

###########################################################################
# Plinks under different prior specs --------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Set parameters template (p fixed)
params_template =
"
num_components: 5

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
  fixed: %g
}
"

# Set values for p
p_fixed <- c(0.2, 0.3, 0.5)

# Sampler runs
for (i in 1:length(p_fixed)) {

  # Set parameters
  params <- sprintf(params_template, p_fixed[i])

  # Run Spatial sampler
  out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

  # Compute plinks
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)

  # Generate plinks plot
  plinks[which(plinks == 0, arr.ind = T)] <- NA
  df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
  plt_plinks <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val)) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey') +
    coord_equal() + theme_void() + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                  barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

  # Save plot
  filename <- sprintf("plt_plinks_fixed_%g.pdf", p_fixed[i])
  pdf(filename, height = 4, width = 4); print(plt_plinks); dev.off()
}

# Set parameters template (p beta distributed)
params_template =
"
num_components: 5

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
    a: %g
    b: %g
  }
}
"

# Set values for p
beta_params <- rbind(c(2,2), c(1,5))

# Sampler runs
for (i in 1:nrow(beta_params)) {

  # Set parameters
  params <- sprintf(params_template, beta_params[i,1], beta_params[i,2])

  # Run Spatial sampler
  out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

  # Compute plinks
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)

  # Generate plinks plot
  plinks[which(plinks == 0, arr.ind = T)] <- NA
  df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
  plt_plinks <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val)) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey') +
    coord_equal() + theme_void() + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                  barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

  # Save plot
  filename <- sprintf("plt_plinks_beta_a%g_b%g.pdf", beta_params[i,1], beta_params[i,2])
  pdf(filename, height = 4, width = 4); print(plt_plinks); dev.off()
}



###########################################################################

###########################################################################
# Prior specs on rho and sigma --------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Set parameters template (rho)
params_template =
"
num_components: 5

p0_params {
  mu0: 0
  a: 2
  b: 2
  lam_: 0.1
}

rho {
  fixed: %g
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
    b: 5
  }
}
"

rho_vals = c(0.9, 0.95, 0.99)

# Sampler runs
for (i in 1:length(rho_vals)) {

  # Set parameters
  params <- sprintf(params_template, rho_vals[i])

  # Run Spatial sampler
  out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

  # Compute plinks
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)

  # Generate plinks plot
  plinks[which(plinks == 0, arr.ind = T)] <- NA
  df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
  plt_plinks <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val)) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey') +
    coord_equal() + theme_void() + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                  barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

  # Save plot
  filename <- sprintf("plt_plinks_rho_%g.pdf", rho_vals[i])
  pdf(filename, height = 4, width = 4); print(plt_plinks); dev.off()
}

# Set parameters template (sigma)
params_template =
"
num_components: 5

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
    alpha: %g
    beta: %g
  }
}

graph_params {
  beta_prior {
    a: 1
    b: 5
  }
}
"

# Values for m
ms <- c(1, 2, 5)

# Sampler runs
for (i in 1:length(ms)) {

  # Set parameters
  params <- sprintf(params_template, 2*(ms[i]^2+1), ms[i]*(2*ms[i]^2+1))

  # Run Spatial sampler
  out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

  # Compute plinks
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)

  # Generate plinks plot
  plinks[which(plinks == 0, arr.ind = T)] <- NA
  df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
  plt_plinks <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val)) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey') +
    coord_equal() + theme_void() + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                  barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

  # Save plot
  filename <- sprintf("plt_plinks_m_%g.pdf", ms[i])
  pdf(filename, height = 4, width = 4); print(plt_plinks); dev.off()
}

# Values for v
vs <- c(0.5, 1, 2)

# Sampler runs
for (i in 1:length(vs)) {

  # Set parameters
  params <- sprintf(params_template, 1/vs[i]+2, 1/vs[i]+1)

  # Run Spatial sampler
  out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

  # Compute plinks
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)

  # Generate plinks plot
  plinks[which(plinks == 0, arr.ind = T)] <- NA
  df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
  plt_plinks <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val)) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey') +
    coord_equal() + theme_void() + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                  barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

  # Save plot
  filename <- sprintf("plt_plinks_v_%g.pdf", vs[i])
  pdf(filename, height = 4, width = 4); print(plt_plinks); dev.off()
}

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
df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")

plt_plinks <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val)) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient2(low='steelblue', mid = "lightgrey", high = 'darkorange', midpoint = 0.5, na.value = 'lightgrey') +
  coord_equal() + theme_void() + theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5))

# Plot estimated graph
df <- reshape2::melt(G_est, c("x","y"), value.name = "val")

plt_Gest <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient(low = 'white', high = 'darkorange',
                      guide = guide_legend(override.aes = list(alpha = 0), title="", title.position="bottom", label.position="bottom")) +
  theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))

# Visualization
x11(height = 6, width = 8.27); gridExtra::grid.arrange(grobs=plots_area, nrow=2, ncol=3)
x11(height = 4, width = 8); gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2)

# Save pdf
# dev.copy2pdf(device=x11, file="BD_fixed_p05_plinks.pdf"); dev.off()
# dev.copy2pdf(device=x11, file="BD_fixed_p05_Densities.pdf"); dev.off()

###########################################################################
