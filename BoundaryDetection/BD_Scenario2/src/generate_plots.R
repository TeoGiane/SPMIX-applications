# # ---- BOUNDAY DETECTION HIGH DIMENSION: GENERATE PLOTS ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_plot", hide.opts = TRUE,
                         description = "Generate plots for the analysis of the California Census Dataset")
opt_parser <- add_argument(opt_parser, arg = "--data-file", type = "character", default = NULL,
                           help = "Relative path to the data file to use.")
opt_parser <- add_argument(opt_parser, arg = "--sim-file", type = "character", default = NULL,
                           help = "Relative path to the simulation file to use.")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "plots",
                           help = "Relative path to save the generated plots")
extra_args <- parse_args(opt_parser)

# Set working directory relative to the script location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # Running in RStudio
  setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
} else {
  # Running from command line
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  setwd(dirname(dirname(script.name)))
}
cat("Setting working directory to: ", getwd(), "\n")

# Check if data file exists
data_file <- file.path(getwd(), extra_args$data_file)
if(!file.exists(data_file)){
  stop(sprintf("Data file '%s' does not exists.", data_file))
}

# Check if simulation requested exists
sim_file <- file.path(getwd(), extra_args$sim_file)
if(!file.exists(sim_file)){
  stop(sprintf("Simulation '%s' does not exists.", sim_file))
}

# Generate output directory if it does not exist
output_dir <- file.path(getwd(), extra_args$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory: ", output_dir, "\n")
}


# Main code ---------------------------------------------------------------

# Required libraries
suppressMessages(library("ggplot2"))
# suppressMessages(library("ggridges"))
# suppressMessages(library("ggrepel"))
# suppressMessages(library("ggmap"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))
suppressMessages(library("SPMIX"))
# suppressMessages(library("dplyr"))

boundary_geometry <- function(boundary_list, sf_geometry) {
  if (!inherits(sf_geometry, "sf")) {
    stop("'sf_geometry' must be an sf object")
  }
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }

  # Create empty list
  geom_bdd <- list()
  count <- 1
  # Populate
  for(i in 1:nrow(sf_geometry)) {
    # Get current area and its boundaries
    if (length(boundary_list[[i]]) > 0) {
      for (j in boundary_list[[i]]) {
        # Compute geometry of boundary
        sel_geom <- sf_geometry[c(i, j),"id"]
        bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
        bounds <- st_union(st_geometry(bounds[bounds$id != bounds$id.1, ]))
        # Add to list
        geom_bdd[[count]] <- st_sf(geometry = bounds)
        count <- count + 1
      }
    }
  }
  # Bind all objects
  geom_bdd <- do.call(rbind, geom_bdd)

  # Drop points if present
  points <- which(attr(geom_bdd$geometry, "classes") == "POINT")
  if(length(points) > 0){
    geom_bdd <- geom_bdd[-points, ]
  }

  # Condense everything into a unique sf object and return
  return(geom_bdd)
}

L1_distance <- function(y1, y2, x) {
  I <- diff(range(x))
  return(I * mean(abs(y1 - y2)))
}

# Generate shapefile and adjacency matrix from scratch
numGroups <- 36
box <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
grid <- st_make_grid(box, cellsize = 1/sqrt(numGroups))
adj_list <- spdep::poly2nb(grid, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")
dist_picker <- c(rep(cbind(rep(1,3),rep(2,3)),3), rep(cbind(rep(2,3),rep(1,3)),3))
sf_grid <- st_sf(data.frame("Group" = as.factor(dist_picker)), geometry = grid)
cat("Generated shapefile and adjacency matrix from scratch\n")

# Import data
load(data_file)
cat(sprintf("Data imported from: %s\n", data_file)) # log

# Import output
load(sim_file)
cat(sprintf("MCMC chain imported from: %s\n", sim_file)) # log

# Deserialization
chains <- sapply(SPMIX_fit, function(x) DeserializeSPMIXProto("spmix.UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute posterior mean and variance for each area
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$mean))
vars_chain <- lapply(chains, function(x) sapply(x$atoms, function(y) y$stdev^2))
weights_chain <- lapply(chains, function(x) t(sapply(x$groupParams, function(y) y$weights)))
post_means <- matrix(nrow = numGroups, ncol=length(chains))
post_vars <- matrix(nrow = numGroups, ncol=length(chains))
for (j in 1:length(chains)) {
  post_means[,j] <- weights_chain[[j]] %*% means_chain[[j]]
  second_moment <- as.vector(weights_chain[[j]] %*% (vars_chain[[j]] + means_chain[[j]]^2))
  post_vars[,j] <- second_moment - post_means[,j]^2
}

# Add to shapefile dataframe
sf_grid$post_mean <- apply(post_means, 1, mean)
sf_grid$post_var <- apply(post_vars, 1, mean)

# Compute estimated density
x <- seq(range(data)[1], range(data)[2], length.out = 500)
estimated_densities <- lapply(ComputePredictiveLPDFs(SPMIX_fit, x), exp)

# Compute admissible and non-admissible edges
admissible_edges <- which(W != 0, arr.ind = T)
non_admissible_edges <- which(W == 0, arr.ind = T)

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = T)] <- NA

# Compute median graph
G_est <- matrix(NA, nrow(plinks), ncol(plinks))
G_est[admissible_edges] <- ifelse(plinks[admissible_edges] >= 0.5, 1, NA)

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[admissible_edges] <- as.factor(ifelse(plinks[admissible_edges] < 0.5, 1, NA))

# Compute boundary geometry
if(!all(is.na(Gb))){
  bound_list <- apply(Gb, 1, function(x){which(x == 1)})
  bound_sf <- boundary_geometry(bound_list, sf_grid)
} else {
  cat("No Boundaries have been found\n")
}

# Posterior of H - barplot
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plt_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("N° of Components")
# Show plot
pdf("plt_postH.pdf", height = 4, width = 4); print(plt_postH); dev.off()

# Posterior of H - Traceplot
df <- data.frame("Iteration"=1:length(chains), "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plt_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("N° of Components") + geom_segment(linewidth = 0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# Show plot
pdf("plt_traceH.pdf", height = 4, width = 4); print(plt_traceH); dev.off()

# Plot plinks matrix
plinks_df <- reshape2::melt(plinks, c("x", "y"), value.name = "val")
Gb_df <- reshape2::melt(Gb, c("x", "y"), value.name="Gb") %>% na.omit()
plt_plinks <- ggplot() +
  geom_tile(data = plinks_df, aes(x=x, y=y, fill=val)) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=numGroups+0.5, ymax=numGroups+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient2(low='steelblue', mid = "white", high = 'darkorange', midpoint = 0.5, na.value = 'white',
                       guide = guide_colourbar(title = "Post. Prob. of Inclusion", direction = "horizontal",
                                               barwidth = unit(2.5, "in"), title.position = "bottom", title.hjust = 0.5)) +
  coord_equal() + theme_void() + theme(legend.position = "bottom")
if(!all(is.na(Gb))){
  plt_plinks <- plt_plinks + geom_tile(data = Gb_df, aes(x=x,y=y), fill=NA, col='darkred', linewidth=0.5)
}

# Plot estimated graph
df <- reshape2::melt(G_est, c("x","y"), value.name = "val")
plt_Gest <- ggplot() +
  geom_tile(data = df, aes(x=x, y=y, fill=val), width=1, height=1) +
  geom_rect(aes(xmin=0.5, ymin=0.5, xmax=nrow(G_est)+0.5, ymax=ncol(G_est)+0.5), col='gray25', fill=NA, linewidth=1) +
  scale_fill_gradient(low = 'darkorange', high = 'darkorange', na.value = "white",
                      guide = guide_legend(override.aes = list(alpha = 0), title="",
                                           title.position="bottom", label.position="bottom")) +
  theme_void() + coord_equal() + theme(legend.position = "bottom", legend.text = element_text(colour = "transparent"))
# Show plots
pdf("plt_plinks.pdf", height = 4, width = 4); print(plt_plinks); dev.off()
pdf("plt_Gest.pdf", height = 4, width = 4); print(plt_Gest); dev.off()
# gridExtra::grid.arrange(plt_plinks, plt_Gest, ncol=2)

# Plot boundaries on the grid + group
plt_boundaries <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=Group), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_manual(values = c("steelblue","darkorange"), labels = c("Student's t", "Skew Normal"),
                    guide = guide_legend(title = NULL, title.hjust = 0.5, label.position = "bottom",
                                         direction = "horizontal", title.position = "bottom", keywidth = unit(1,"in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
pdf("plt_boundaries.pdf", height = 4, width = 4); print(plt_boundaries); dev.off()
# Save
# pdf("plt_BDgroup.pdf", height = 4, width = 4); plt_boundaries; dev.off()

# Plot boundaries on the grid + post_mean
plt_postmean <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=post_mean), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(title = "Post. Mean", title.hjust = 0.5, title.position = "bottom",
                                             direction = "horizontal", barwidth = unit(2.5, "in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Show plot
pdf("plt_postmean.pdf", height = 4, width = 4); print(plt_postmean); dev.off()
# Save
# pdf("plt_BDpostMean.pdf", height = 4, width = 4); plt_postmean; dev.off()

# Plot boundaries on the grid + post_var
plt_postvar <- ggplot() +
  geom_sf(data = sf_grid, aes(fill=post_var), col='gray25', alpha = 0.6, linewidth=0.4) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(title = "Post. Variance", title.hjust = 0.5, title.position = "bottom",
                                             direction = "horizontal", barwidth = unit(2.5, "in"))) +
  geom_sf(data = bound_sf, fill=NA, col='darkred', linewidth=1.2) +
  theme_void() + theme(legend.position = "bottom")
# Save plot
pdf("plt_postvar.pdf", height = 4, width = 4); print(plt_postvar); dev.off()

# Plot - Empirical density histogram in bordering areas
areas <- c(3,4)
plt_areas <- list()
for (i in 1:length(areas)) {
  df <- data.frame("x" = data[[areas[i]]])
  plt_areas[[i]] <- ggplot(data = df, aes(x=x, y = after_stat(density))) +
    geom_histogram(col=NA, fill='white', bins = 50) +
    geom_histogram(col='steelblue', fill='steelblue', alpha = 0.4, bins = 50) +
    xlim(range(data[areas])) + ylim(0, 0.55) +
    xlab("Data") + ylab("Density") + theme(plot.title = element_text(hjust = 0.5))
}
# Save plots
pdf(sprintf("plt_areas%d.pdf", areas[1]), height = 4, width = 4); print(plt_areas[[1]]); dev.off()
pdf(sprintf("plt_areas%d.pdf", areas[2]), height = 4, width = 4); print(plt_areas[[2]]); dev.off()

# # Plot - Empirical density histogram + Estimated density
plt_areasdens <- list()
for (i in 1:length(areas)) {
  df <- data.frame("x" = x, "y" = colMeans(estimated_densities[[areas[i]]]))
  plt_areasdens[[i]] <- plt_areas[[i]] +
    geom_line(data = df, aes(x=x, y=y), col='darkorange', linewidth=1.2) +
    ylim(0, 0.55) + theme(plot.title = element_text(hjust = 0.5))
}
# Save plots
pdf(sprintf("plt_areasdens_%d.pdf", areas[1]), height = 4, width = 4); print(plt_areasdens[[1]]); dev.off()
pdf(sprintf("plt_areasdens_%d.pdf", areas[2]), height = 4, width = 4); print(plt_areasdens[[2]]); dev.off()

# Plot - Traceplot of p
p_chain <- sapply(chains, function(x){x$p})
df <- data.frame("Iter"=1:length(p_chain), "Value"=p_chain)
plt_p_chain <- ggplot() +
  geom_line(data = df, aes(x = Iter, y = Value)) +
  xlab('Iteration') + ylab('p')
# Save plot
pdf("plt_p_chain.pdf", height = 4, width = 4); print(plt_p_chain); dev.off()

# Plot - Traceplot of |G|
Nedge_chain <- sapply(G_chain, function(x){sum(x[upper.tri(x)])})
df <- data.frame("Iter"=1:length(Nedge_chain), "Value"=Nedge_chain)
plt_Nedge_chain <- ggplot() +
  geom_line(data = df, aes(x = Iter, y = Value)) +
  xlab('Iteration') + ylab('|G|')
# Save plot
pdf("plt_Nedge_chain.pdf", height = 4, width = 4); print(plt_Nedge_chain); dev.off()

# Plot - Traceplot of sigma^2
sigma_chain <- sapply(chains, function(x){x$Sigma$data[1]})
df <- data.frame("Iter"=1:length(sigma_chain), "Value"=sigma_chain)
plt_sigma_chain <- ggplot() +
  geom_line(data = df, aes(x = Iter, y = Value)) +
  xlab('Iteration') + ylab(bquote(sigma^2))
# Save plot
pdf("plt_sigma_chain.pdf", height = 4, width = 4); print(plt_sigma_chain); dev.off()

cat(sprintf("All plots have been generated and saved in: %s\n", output_dir)) # log


# # ---- END OF SCRIPT ---- # #

# ESS, WAIC, robe per review ----------------------------------------------

# # ESS of p
# mcmcse::ess(p_chain)
# # ESS of |G|
# mcmcse::ess(Nedge_chain)
# ESS of sigma^2
# mcmcse::ess(sigma_chain)

# # TODO: improve doc for WAIC
# tmp <- ComputePosteriorLPDF(data, chains, verbose = T)
#
# # Plot - |G| with different initializations (DA SISTEMARE)
# load("output/Nedge_chain_startempty.dat")
# df <- data.frame("Iter" = 1:length(Nedge_chain), "Value"=Nedge_chain, "Group"=rep("Empty", length(Nedge_chain)))
# load("output/Nedge_chain_startfull.dat")
# df <- rbind(df, data.frame("Iter" = 1:length(Nedge_chain), "Value"=Nedge_chain, "Group"=rep("Full", length(Nedge_chain))))
# load("output/Nedge_chain_startrandom.dat")
# df <- rbind(df, data.frame("Iter" = 1:length(Nedge_chain), "Value"=Nedge_chain, "Group"=rep("Random", length(Nedge_chain))))
# #Generate plot
# plt_Nedge_chain_diffInit <- ggplot() +
#   geom_line(data=df, aes(x=Iter, y=Value, color=Group)) +
#   scale_color_manual(values = c("Full"="darkorange", "Empty"='steelblue', "Random"='forestgreen')) +
#   xlab('Iteration') + ylab('|G|') +
#   guides(color = guide_legend("Initial Value", direction = 'horizontal', position = 'bottom'))
# # Show / Save
# # x11(height = 4, width = 4); plt_Nedge_chain_diffInit
# pdf("plots/plt_Nedge_chain_diffInit.pdf", height = 4, width = 4); plt_Nedge_chain_diffInit; dev.off()

