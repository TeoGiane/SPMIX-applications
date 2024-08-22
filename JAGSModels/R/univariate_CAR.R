library("rjags")
library("sf")

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

# Prepare data for JAGS
rho = 0.95
Y = sapply(data, median)
CAR_data <- list("Y"=Y, "n"=length(Y), "rho"=0.95, "Gadj"=W)
CAR_model <- jags.model("univariate_CAR.bug", data=CAR_data, n.chains=2)

# Run MCMC
CAR_params <- c("G", "p", "tau")
CAR_fit <- jags.samples(CAR_model, variable.names = CAR_params, n.iter = 1000, thin = 2)

# Get posterior median graph
G_post <- summary(CAR_fit$G, FUN = mean)$stat
Gb <- W - G_post
fields::image.plot(Gb)

# Traceplots for p and sigma^2
plot(CAR_fit$p, type='l', lwd=1.5)
plot(1/CAR_fit$tau, type='l', lwd=1.5)

