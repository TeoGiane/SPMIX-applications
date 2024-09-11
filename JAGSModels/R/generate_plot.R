library("ggplot2")
library("rjags")
library("sf")
library("spdep")

# Import shapefile
sf_counties <- read_sf("data/counties-pumas/counties-pumas.shp")

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Import data
load("data/data_001.dat")

# Set beta parameters
abeta <- 2
bbeta <- 2

# Import chain
filename <- file.path(getwd(), "output", sprintf("MCAR_fit-a%gb%g.dat", abeta, bbeta))
load(filename)

# 100 Iterations of burnin
burnin <- 100
p_chain <- lapply(as.mcmc.list(MCAR_fit$p), function(chain){ as.vector(chain[(burnin+1):dim(chain)[1],1:dim(chain)[2]]) })
p_chain <- do.call(c, p_chain)
tau_chain <- lapply(as.mcmc.list(MCAR_fit$tau), function(chain){ as.vector(chain[(burnin+1):dim(chain)[1],1:dim(chain)[2]]) })
tau_chain <- do.call(c, tau_chain)

# Compute admissible edges
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute PPI matrix (remove not admissible edges)
plinks <- summary(MCAR_fit$G, FUN = mean)$stat + abs(jitter(matrix(0,length(data),length(data))))
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Threshold -> posterior median of p
gamma_graph <- 0.5

# Compute neighbouring graph
Gn <- matrix(NA, nrow(plinks), ncol(plinks))
Gn[Eadj] <- ifelse(plinks[Eadj] >= gamma_graph, 1, NA)

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < gamma_graph, 1, NA)

# PLOT - plinks matrix and with boundary edges in red
plinks_df <- reshape2::melt(plinks, c("x", "y"), value.name="PPI") %>% na.omit()
Gb_df <- reshape2::melt(Gb, c("x","y"), value.name="Gb") %>% na.omit()
# Generate
plt_plinks <- ggplot() +
  geom_tile(data = plinks_df, aes(x=x, y=y, fill=PPI), width=1, height=1) +
  geom_rect(aes(xmin=0.5, xmax=nrow(plinks)+0.5, ymin=0.5, ymax=nrow(plinks)+0.5), fill=NA, color="gray25", linewidth=0.5) +
  scale_fill_gradient2(low='steelblue4', mid = "white", high = 'darkorange', midpoint = 0.5, na.value = 'white',
                       guide = guide_colorbar("Post. Prob. of Inclusion", position = "bottom", direction = "horizontal", barwidth=unit(2.5,"in"),
                                              title.position = "bottom", title.hjust = 0.5, label.vjust = 0.5)) +
  geom_tile(data = Gb_df, aes(x=x,y=y), fill=NA, col='darkred', linewidth=0.5) +
  theme_void() + theme(legend.position = "bottom") + coord_equal()

# PLOT - Traceplot of p
p_df <- data.frame("Iteration"=1:length(p_chain), "p"=p_chain)
# Generate
plt_p <- ggplot() +
  geom_line(data = p_df, aes(x=Iteration, y=p), linewidth=1.2)

# PLOT - traceplot of sigma^2
sigma2_df <- data.frame("Iteration"=1:length(tau_chain), "sigma2"=1/tau_chain)
# Generate
plt_sigma2 <- ggplot() +
  geom_line(data = sigma2_df, aes(x=Iteration, y=sigma2), linewidth=1.2) +
  ylab(bquote(sigma^2))

# Show all plots together
x11(width = 12, height = 4); gridExtra::grid.arrange(plt_plinks, plt_p, plt_sigma2, ncol = 3)
# Save all plots together
filename <- file.path(getwd(),"plots",sprintf("plt_BDSumStats-a%gb%g.pdf", abeta, bbeta))
pdf(filename, width = 12, height = 4); gridExtra::grid.arrange(plt_plinks, plt_p, plt_sigma2, ncol = 3); dev.off()

# Show
# x11(height = 4, width = 4); plt_plinks
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_plinks-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_plinks; dev.off()

# Show
# x11(height = 4, width = 4); plt_p
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_p-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_p; dev.off()

# Show
# x11(height = 4, width = 4); plt_sigma2
# Save
# filename <- file.path(getwd(),"plots",sprintf("plt_sigma2-a%gb%g.pdf", abeta, bbeta))
# pdf(filename, height = 4, width = 4); plt_sigma2; dev.off()
