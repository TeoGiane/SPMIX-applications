# Import required packages
library("ggplot2")
library("ggmap")
library("sf")
library("dplyr")
library("SPMIX")

# Functions
sf_ggmap <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, and convert back to a bbox
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  
  # Return
  return(map)
}
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

# Import usa shapefile
sf_usa <- read_sf("shp/us-pumas/us-pumas.shp")
sf_usa$id <- row.names(sf_usa)
sf_usa$PUMA <- as.character(as.numeric(sf_usa$PUMA))

# Select only LA + Ventura + Orange counties, California
sf_cali <- sf_usa[sf_usa$State == "California", ]
sel_county <- c("Los Angeles County", "Ventura County", "Orange County")
sf_counties <- sf_cali[grep(paste(sel_county, collapse = "|"), sf_cali$Name), ]

# Import raw data
rdata <- read.csv("data/psam_p06.csv")
rdata <- rdata[which(rdata$PINCP > 0), ]
rdata$PUMA <- as.character(rdata$PUMA)

# Subsample data for SPMIX
pumas <- sf_counties$PUMA
data <- list()
for (i in 1:length(pumas)) {
  set.seed(230196)
  tmp <- rdata$PINCP[which(rdata$PUMA == pumas[i])]
  data[[i]] <- log(tmp[sample(1:length(tmp), 100)])
}

# Compute adjacency matrix
adj_list <- spdep::poly2nb(sf_counties, row.names = pumas, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")


###########################################################################
# Sampler run -------------------------------------------------------------

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
    fixed: 0.2
  }
  "

# Run Spatial sampler
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params)

# Save output
if (exists("out")) {
  filename <- sprintf("output/chain_%s_burn%d_iter%d.dat", format(Sys.time(), "%Y%m%d_%H%M"), burnin, niter)
  save(out, file = filename)
}

###########################################################################

###########################################################################
# Posterior Analysis ------------------------------------------------------

# Load output
load("output/chain_20230208_2225_burn5000_iter5000.dat")

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute posterior mean and variance for each PUMA
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$mean))
vars_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$stdev^2))
weights_chain <- lapply(chains, function(x) t(sapply(x$groupParams, function(y) y$weights)))
post_means <- matrix(nrow = length(pumas), ncol=length(chains))
post_vars <- matrix(nrow = length(pumas), ncol=length(chains))
for (j in 1:length(chains)) {
  post_means[,j] <- weights_chain[[j]] %*% means_chain[[j]]
  post_vars[,j] <- weights_chain[[j]]^2 %*% vars_chain[[j]]
}

# Add to shapefile dataframe
sf_counties$post_mean <- apply(post_means, 1, mean)
sf_counties$post_var <- apply(post_vars, 1, mean)

# Compute estimated density
data_ranges <- sapply(data, range); Npoints <- 500
estimated_densities <- ComputeDensities(chains, Npoints, data_ranges, alpha = 0.05)

# Compute plinks and median graph estimate
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

# Compute boundary matrix, boundary adj list and geometry
bound_matrix <- W - G_est
bound_list <- spdep::mat2listw(bound_matrix)$neighbours
bound_sf <- boundary_geometry(bound_list, sf_counties)

###########################################################################

###########################################################################
# PLOT - Posterior distribution of H --------------------------------------

# Make auxiliary dataframe
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
# Generate plot
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components")
# Show plot
plot_postH

###########################################################################

###########################################################################
# PLOT - Empirical vs. Estimated density ----------------------------------

# Make auxiliary dataframes
area <- 60
df_hist <- data.frame("Data" = data[[area]])
df_dens <- data.frame("x" = seq(data_ranges[1,area], data_ranges[2,area], length.out = Npoints),
                      t(estimated_densities[[area]]))
# Generate plot
plt_denscompare <- ggplot() +
  geom_histogram(data = df_hist, aes(x = Data, after_stat(density)), color="steelblue", fill="lightblue") +
  geom_ribbon(data = df_dens, aes(x = x, ymin = low, ymax = up), fill='darkorange', alpha = 0.35) +
  geom_line(data = df_dens, aes(x = x, y = est), col='darkorange', linewidth = 1) +
  ylab("Density")
# Clean auxiliary dataframes
rm(list = c("df_hist", "df_dens"))

# Show
x11(height = 4, width = 4); plt_denscompare

###########################################################################

###########################################################################
# PLOT - plinks matrix ----------------------------------------------------

# Make auxiliary dataframe
df <- reshape2::melt(plinks, c("x", "y"), value.name = "PPI")
df[which(df$PPI == 0), "PPI"] <- NA
# Generate plot
plt_plinks <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=PPI), width=1, height=1) +
  geom_rect(aes(xmin=0.5, xmax=nrow(plinks)+0.5, ymin=0.5, ymax=nrow(plinks)+0.5), fill=NA, color="grey25", linewidth=0.5) +
  scale_fill_gradient(low='steelblue', high = "darkorange", na.value = "white") +
  theme_void() + theme(legend.position = "bottom") + coord_equal() +
  guides(fill = guide_colorbar(title = "Post. Prob. of Link Inclusion",
                               title.position = "bottom", title.hjust = 0.5,
                               barwidth = unit(2.5,"in")))
# Clean auxiliary dataframe
rm('df')

# Show
plt_plinks

###########################################################################

###########################################################################
# PLOT - Detected boundaries on the map -----------------------------------

# Get map with proper attributes
counties_bbox <- unname(st_bbox(st_transform(sf_counties, 4326)))
counties_map <- sf_ggmap(get_map(counties_bbox, source = "stamen", crop = FALSE))
# CRS conversions (for plotting)
sf_counties_3857 <- st_transform(sf_counties, 3857)
bound_sf_3857 <- st_transform(bound_sf, 3857)
# Generate plot
plt_boundaries_mean <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange') +
  geom_sf(data = bound_sf, col='darkred', lwd=0.8, inherit.aes = FALSE) +
  theme_void() + theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title = "Post. Mean", direction = "horizontal", barwidth = unit(3, "in"),
                                title.position = "bottom", title.hjust = 0.5))
plt_boundaries_var <- ggmap(counties_map) +
  geom_sf(data = sf_counties_3857, aes(fill=post_var), col='gray25', alpha = 0.6, inherit.aes = FALSE) +
  scale_fill_gradient(low = 'steelblue', high = 'darkorange') +
  geom_sf(data = bound_sf, col='darkred', lwd=0.8, inherit.aes = FALSE) +
  theme_void() + theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(title = "Post. Variance", direction = "horizontal", barwidth = unit(3, "in"),
                                title.position = "bottom", title.hjust = 0.5))

# View
x11(height = 7, width = 14); gridExtra::grid.arrange(grobs = list(plt_boundaries_mean, plt_boundaries_var), ncol=2)

###########################################################################

# Quick comparison
x11(); par(mfrow=c(1,2))
area <- 4
hist(data[[area]], probability = T, ylim = c(0,0.5))
lines(seq(data_ranges[1,area], data_ranges[2,area], length.out = Npoints), estimated_densities[[area]]["est",], col='blue', lwd=2)
area <- 11
hist(data[[area]], probability = T, ylim = c(0,0.5))
lines(seq(data_ranges[1,area], data_ranges[2,area], length.out = Npoints), estimated_densities[[area]]["est",], col='blue', lwd=2)



# Script di supporto per aggiornare primi_test.R nel plotting
# # Import usa shape file
# sf_usa <- read_sf("shapefiles/us-pumas_boundaries/ipums_puma_2010.shp"); sf_usa$id <- row.names(sf_usa)
# 
# # Select only LA + Ventura + Orange counties, CA
# sf_cali <- sf_usa[sf_usa$State=="California",]
# sel_county <- c("Los Angeles County", "Ventura County", "Orange County")
# sf_counties <- sf_cali[grep(paste(sel_county, collapse = "|"), sf_cali$Name), ]
# 
# # Find internal borders (prima prova per fare il plot automatico dei boundaries)
# borders <- st_intersection(sf_counties, sf_counties)
# inner_borders <- st_geometry(borders[borders$id != borders$id.1, ])
# 
# # Get map with proper attributes
# counties_bbox <- unname(st_bbox(st_transform(sf_counties, 4326)))
# counties_map <- sf_ggmap(get_map(counties_bbox, source = "stamen", crop = FALSE))
# 
# # CRS conversions (for plotting)
# sf_counties_3857 <- st_transform(sf_counties, 3857)
# inner_borders_3857 <- st_transform(inner_borders, 3857)
# 
# # Plot
# ggmap(counties_map) +
#   geom_sf(data = sf_counties_3857, fill="white", alpha=0.5, col='gray25', inherit.aes = FALSE) +
#   geom_sf(data = inner_borders_3857, col='red', lwd=1.2, inherit.aes = FALSE) +
#   theme_void()
