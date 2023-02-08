## First test for the Reversible Jump Sampler ##

# Required libraries
library("SPMIX")
library("ggplot2")

###########################################################################
# Data Simulation ---------------------------------------------------------

# Generate data (1 location, mixture of 3 normals)
set.seed(230196)
ngroups <- 1; ncomponents <- 3; N <- 1000
means <- c(-4,1,5); std_devs <- c(1,1,1)
weights <- matrix(c(2/6,3/6,1/6), ngroups, ncomponents, TRUE)
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
thin = 1

# Sampler parameters
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
      a: 1
      b: 5
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
# Analysis and Visualization ----------------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState", x))
H_chain <- sapply(chains, function(x) x$num_components)

# Barplot of the estimated posterior for H
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components") #+ ggtitle("Posterior of H")
# Show plot
x11(height = 4, width = 4); plot_postH

# Traceplot for the whole chain (no burnin, no thinning)
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plot_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.2) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# Show plot
x11(height = 4, width = 4); plot_traceH

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
# Show plot
x11(width = 6, height = 4); plot_densCompare

###########################################################################

###########################################################################
# Sensitivity w.r.t. initial value of H -----------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Sampler parameters
params_template =
  "
num_components: %d

p0_params {
  mu0: 0
  a: 2
  b: 2
  lam_: 0.1
}

rho {
  beta_prior {
    a: 1
    b: 5
  }
}

sigma {
	inv_gamma_prior {
		alpha: 2
		beta: 2
	}
}
"

# Specify range of parameters and run simulations
H_init = c(2,4,8,10)
for (i in 1:length(H_init)) {
  # Run Spatial sampler
  params <- sprintf(params_template, H_init[i])
  out <- Sampler.DensityEstimation(burnin, niter, thin, data, W, params, type="rjmcmc")

  # Deserialize chain
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  H_chain <- sapply(chains, function(x) x$num_components)

  # Plot
  df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
  tmp <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
    geom_bar(stat="identity", color="steelblue", fill="lightblue") +
    ylim(c(0,1)) + xlab("N° of Components")

  # Save to pdf
  filename <- sprintf("plt_postH_Hinit%d.pdf", H_init[i])
  pdf(filename, height = 2, width = 2); print(tmp); dev.off()
}

###########################################################################

###########################################################################
# Sampler Execution (full run) --------------------------------------------

# Setting MCMC parameters
burnin = 0
niter = 10000
thin = 1

# Sampler parameters
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
    a: 1
    b: 5
  }
}

sigma {
	inv_gamma_prior {
		alpha: 2
		beta: 2
	}
}

mtilde_sigmasq: 5
"

# Run Spatial sampler
out <- Sampler.DensityEstimation(burnin, niter, thin, data, W, params, type = "rjmcmc")

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
# Show plot
x11(height = 4, width = 4); plot_traceH

###########################################################################

###########################################################################
# Sensitivity w.r.t. number of observations -------------------------------

# Sampler parameters
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
    a: 1
    b: 5
  }
}

sigma {
	inv_gamma_prior {
		alpha: 2
		beta: 2
	}
}
"

# Set number of observations
Ns <- c(100,250,500,1000)

# Run simulations
for (i in 1:length(Ns)) {

  # Generate data (1 location, mixture of 3 normals)
  set.seed(230196)
  ngroups <- 1; ncomponents <- 3; N <- Ns[i]
  means <- c(-4,1,5); std_devs <- c(1,1,1); weights <- c(1/6,3/6,2/6)
  cluster_alloc <- sample(1:ncomponents, prob = weights, size = N, replace = T)
  data <- list(); data[[1]] <- rnorm(N, mean = means[cluster_alloc], sd = std_devs[cluster_alloc])

  # Run Spatial sampler
  out <- Sampler.DensityEstimation(5000, 5000, 2, data, W, params, type = "rjmcmc")

  # Deserialization
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  H_chain <- sapply(chains, function(x) x$num_components)

  # Traceplot for the whole chain (no burnin, no thinning)
  df <- data.frame("Iteration"=1:length(H_chain), "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
  tmp <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
    ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1)
  rm(list='df')

  # Save pdf
  filename <- sprintf("plt_traceH_N%d.pdf", N)
  pdf(filename, height = 2, width = 2); print(tmp); dev.off()
}

# Visualization
x11(height = 6, width = 6); gridExtra::grid.arrange(grobs=plots, nrow=2, ncol=2)

###########################################################################
