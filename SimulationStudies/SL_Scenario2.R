# Required libraries
library("SPMIX")
library("ggplot2")

# Generating data
set.seed(230196)
I <- 6; Ns <- rep(100, I)
data <- list()
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

# Setting initial W
W <- matrix(1,I,I)

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 1

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
      a: 2
      b: 2
    }
  }
  "

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

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
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=I+0.5, ymax=I+0.5), col='gray25', fill=NA, linewidth=1) +
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


# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 1

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
p_fixed <- c(0.1, 0.2, 0.3)

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
  G_est <- ifelse(plinks > 0.5, 1, 0)

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

  # Generate Gest plot
  df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
  plt_Gest <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient(low = 'white', high = 'darkorange',
                        guide = guide_legend(override.aes = list(alpha = 0), title="",
                                             title.position="bottom", label.position="bottom")) +
    theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))

  # Save plot
  filename <- sprintf("plt_plinksGest_fixed_%g.pdf", p_fixed[i])
  pdf(filename, height = 4, width = 6); gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2); dev.off()
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
beta_params <- rbind(c(1,5), c(2,5), c(2,2))

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
  G_est <- ifelse(plinks > 0.5, 1, 0)

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
  # Generate Gest plot
  df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
  plt_Gest <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient(low = 'white', high = 'darkorange',
                        guide = guide_legend(override.aes = list(alpha = 0), title="",
                                             title.position="bottom", label.position="bottom")) +
    theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))

  # Save plot
  filename <- sprintf("plt_plinksGest_beta_a%g_b%g.pdf", beta_params[i,1], beta_params[i,2])
  pdf(filename, height = 4, width = 6); gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2); dev.off()
}

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

  # Compute plinks and Gest
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)
  G_est <- ifelse(plinks > 0.5, 1, 0)

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

  # Generate Gest plot
  df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
  plt_Gest <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient(low = 'white', high = 'darkorange',
                        guide = guide_legend(override.aes = list(alpha = 0), title="",
                                             title.position="bottom", label.position="bottom")) +
    theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))

  # Save plot
  filename <- sprintf("plt_plinksGest_rho_%g.pdf", rho_vals[i])
  pdf(filename, height = 4, width = 6); gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2); dev.off()
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
  G_est <- ifelse(plinks > 0.5, 1, 0)

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

  # Generate Gest plot
  df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
  plt_Gest <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient(low = 'white', high = 'darkorange',
                        guide = guide_legend(override.aes = list(alpha = 0), title="",
                                             title.position="bottom", label.position="bottom")) +
    theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))

  # Save plot
  filename <- sprintf("plt_plinksGest_m_%g.pdf", ms[i])
  pdf(filename, height = 4, width = 6); gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2); dev.off()
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
  G_est <- ifelse(plinks > 0.5, 1, 0)

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

  # Generate Gest plot
  df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
  plt_Gest <- ggplot() +
    geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
    geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
    scale_fill_gradient(low = 'white', high = 'darkorange',
                        guide = guide_legend(override.aes = list(alpha = 0), title="",
                                             title.position="bottom", label.position="bottom")) +
    theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))

  # Save plot
  filename <- sprintf("plt_plinksGest_v_%g.pdf", vs[i])
  pdf(filename, height = 4, width = 6); gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2); dev.off()
}
