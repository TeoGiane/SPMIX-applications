# Preliminary Checks ------------------------------------------------------

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(dirname(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Create directory for output if does not exist
out_dir <- file.path(getwd(), "output")
if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(out_dir))) # Log


# Main Code ---------------------------------------------------------------

# Required libraries
suppressMessages(library("rjags"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))

# Import shapefile
sf_counties <- read_sf("data/counties-pumas/counties-pumas.shp")

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Import data
load("data/full_dataset.dat")

# Get number of areas and prior hypers
numGroups <- length(data)
abeta <- 2
bbeta <- numGroups

# # Generate shape file from scratch
# numGroups <- 36
# box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
# grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))
# 
# # Compute adjacency matrix
# adj_list <- spdep::poly2nb(grid, queen = FALSE)
# W <- spdep::nb2mat(adj_list, style = "B")
# 
# # Generate distribution picker label
# dist_picker <- c(rep(cbind(rep(1,3),rep(2,3)),3), rep(cbind(rep(2,3),rep(1,3)),3))
# 
# # Generate sf object
# sf_grid <- st_sf(data.frame("Group" = as.factor(dist_picker)), geometry = grid)
# 
# # Generate data
# set.seed(230196); data <- list()
# for (i in 1:numGroups) {
#   if (dist_picker[i] == 1) {
#     data[[i]] <- metRology::rt.scaled(n=100, df=6, mean=4, sd=1.5)
#   } else {
#     data[[i]] <- sn::rsn(n=100, xi=4, omega=1.3, alpha=-3)
#     attributes(data[[i]]) <- NULL
#   }
# }

# Prepare data for JAGS
Y = c(sapply(data, quantile, 0.05),
      sapply(data, quantile, 0.25),
      sapply(data, quantile, 0.5),
      sapply(data, quantile, 0.75),
      sapply(data, quantile, 0.95)); names(Y) = NULL
MCAR_data <- list("Y"=Y,
                  "I"=numGroups,
                  "rho"=0.95,
                  "W"=W,
                  "Ncov"=5,
                  "mu0"=mean(Y),
                  "abeta"=abeta,
                  "bbeta"=bbeta)
MCAR_model <- jags.model("bug/multivariate_CAR_marginalized.bug", data=MCAR_data, n.chains=1, n.adapt = 0)

# Run MCMC
MCAR_params <- c("G", "p", "one_over_sigmasq", "one_over_tausq")
MCAR_fit <- jags.samples(MCAR_model, variable.names = MCAR_params, n.iter = 100, thin = 1)

# Store MCAR_fit
filename <- sprintf("output/MCAR_marginalized_fit-a%gb%g.dat", abeta, bbeta)
save(MCAR_fit, file = filename)

# # Get posterior median graph
# G_post <- summary(MCAR_fit$G, FUN = mean)$stat
# Gb <- W - G_post
# fields::image.plot(Gb)
# 
# # Traceplots for p and sigma^2
# plot(MCAR_fit$p, type='l', lwd=1.5)
# plot(1/MCAR_fit$tau, type='l', lwd=1.5)
# 
