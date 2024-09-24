# Import required packages
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggmap"))
suppressMessages(library("sf"))
suppressMessages(library("spdep"))
suppressMessages(library("SPMIX"))

# Select sim_folder
sim_folder <- file.path(getwd(), "output/H_RJ/rho_0.95")

# log
cat(sprintf("Current Directory: %s\n", getwd()))

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

boundary_geometry <- function(boundary_graph, sf_geometry) {
  # Check
  if(!inherits(sf_geometry, "sf")) { stop("'sf_geometry' must be an sf object") }
  
  # Add id column
  if(!("id" %in% names(sf_geometry))){
    sf_geometry$id <- 1:nrow(sf_geometry)
  }
  
  # Generate boundary list using upper-triangular view of boundary_graph
  boundary_graph[lower.tri(boundary_graph)] <- NA
  boundary_list <- apply(boundary_graph, 1, function(x){which(x==1)})
  
  # Create empty list
  geom_bdd <- list(); k <- 1
  for(i in 1:nrow(sf_geometry)) {
    for (j in (boundary_list[[i]])) {
      # Compute boundary geometry between i and j
      sel_geom <- sf_geometry[c(i,j), c("id", "geometry")]
      bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
      bounds <- st_geometry(bounds[bounds$id != bounds$id.1, ])
      # Append to list
      geom_bdd[[k]] <- st_sf(geometry = bounds); k <- k+1
    }
    
    # Get current area and its boundaries
    # sel_geom <- sf_geometry[c(i, boundary_list[[i]]), c("id", "geometry")]
    
    # Compute geometry of boundary
    # bounds <- suppressWarnings(st_intersection(sel_geom, sel_geom))
    # bounds <- st_geometry(bounds[bounds$id != bounds$id.1, ])
    
    # Add to list
    # geom_bdd[[i]] <- st_sf(geometry = bounds)
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

# Import shapefile
sf_counties <- read_sf("data/counties-pumas/counties-pumas.shp")

# Import data
load("data/data_001.dat")

# Load output
load(file.path(sim_folder, "chain_001.dat"))
# cat(sprintf("Loaded chain: %s\n", file.path(sim_folder, "chain_001.dat"))) # log

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Select only LA county
la_county <- grepl("Los Angeles", sf_counties$Name)
sf_la <- sf_counties %>% filter(la_county) %>% st_transform(5070)
data <- data[which(la_county)]
G_chain <- lapply(G_chain, function(G){G[which(la_county),which(la_county)]})

# Compute adjacency list for LA county
adj_list <- poly2nb(sf_la, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# Compute admissible edges
Eadj <- which(W == 1, arr.ind = TRUE)

# Compute PPI matrix (remove not admissible edges)
plinks <- Reduce('+', G_chain)/length(G_chain)
plinks[which(W == 0, arr.ind = TRUE)] <- NA

# Threshold -> posterior median of p
gamma_graph <- 0.5

# Compute neighbouring graph
Gn <- matrix(NA, nrow(plinks), ncol(plinks))
Gn[Eadj] <- ifelse(plinks[Eadj] >= gamma_graph, 1, NA)

# Compute boundary graph
Gb <- matrix(NA, nrow(plinks), ncol(plinks))
Gb[Eadj] <- ifelse(plinks[Eadj] < gamma_graph, 1, NA)

# Compute boundary geometry
bound_sf <- boundary_geometry(Gb, sf_la)

# Compute neighbouring and boundary adjacency lists
neigh_list <- apply(Gn, 1, function(x){which(x==1)})
bound_list <- apply(Gb, 1, function(x){which(x==1)})

# Get LA map and Crop Islands
la_bbox <- unname(st_bbox(st_transform(sf_la, 4326)))
la_bbox[2] <- la_bbox[4] - (la_bbox[3] - la_bbox[1])
la_map <- sf_ggmap(get_map(la_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))

###########################################################################
# Recorded Crimes in LA during 2020 ---------------------------------------

# Get raw data for recorded crimes in LA county
df <- read.csv("data/explain_boundaries/2020_LACounty_Crimes_dataset.csv")

# Select crimes categories
categories <- c("DISORDERLY CONDUCT", "LIQUOR LAWS", "RECEIVING STOLEN PROPERTY", "VAGRANCY")

# Create sf dataset
sf_crimes <- st_as_sf(df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
  filter(CATEGORY %in% categories) %>%
  mutate(CATEGORY = as.factor(CATEGORY)) %>%
  select(CATEGORY, geometry)

# PLOT - Recorded non-violent crimes during 2020
sf_la_3857 <- st_transform(sf_la, 3857)
sf_crimes_3857 <- st_transform(sf_crimes, 3857)
bound_sf_3857 <- st_transform(bound_sf, 3857)
# Generate
plt_crimescat <- ggmap(la_map) +
  geom_sf(data = sf_la_3857, col='gray25', linewidth=0.25, fill='white', alpha=0.6, inherit.aes=F) +
  geom_sf(data = sf_crimes_3857, aes(color=CATEGORY), size=1, inherit.aes=F) +
  geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, inherit.aes=F) +
  scale_color_manual(values = c("steelblue", 'darkorange', 'lightgreen', 'salmon'),
                    labels = c("Disorderly Conduct", "Liquor Laws", "Recieving Stolen Property", "Vagrancy")) +
  guides(color = guide_legend(NULL, position = "bottom", direction = "horizontal", title.position = "bottom",
                             title.hjust = 0.5, nrow = 2, override.aes = list(size=2.5))) +
  theme_void()
# View / Save
x11(height = 4, width = 4); plt_crimescat
# pdf("plots/plt_CrimesCategories.pdf", height = 4, width = 4); plt_crimescat; dev.off()

###########################################################################

###########################################################################
# % of people without medical insurance -----------------------------------

sf_insurance <- read_sf("data/explain_boundaries/Health_Insurance_(census_tract)/Health_Insurance_(census_tract).shp") %>%
  group_by(csa) %>%
  summarise(UNINSURED = mean(uninsure_1, na.rm = T))

# PLOT- % of people without medical insurance
sf_la_3857 <- st_transform(sf_la, 3857)
sf_insurance_3857 <- st_transform(sf_insurance, 3857)
bound_sf_3857 <- st_transform(bound_sf, 3857)
# Generate
plt_insurance <- ggmap(la_map) +
  geom_sf(data = sf_la_3857, col=NA, fill='white', alpha=0.6, inherit.aes=F) +
  geom_sf(data = sf_insurance_3857, aes(fill=UNINSURED, color=UNINSURED), inherit.aes=F) +
  geom_sf(data = sf_la_3857, col="gray25", fill=NA, linewidth=0.25, inherit.aes=F) +
  geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
  scale_fill_gradient(low="steelblue", high="darkorange") +
  scale_color_gradient(low="steelblue", high="darkorange") +
  guides(fill=guide_colorbar("Population w/o Health Insurance (%)", position = "bottom", direction = "horizontal",
                             title.position = "bottom", title.hjust = 0.5, barwidth = unit(3,"in"), barheight = unit(0.12, "in")),
         color = "none") + theme_void()
pdf("plots/plt_NoInsurance.pdf", height = 4, width = 4); plt_insurance; dev.off()

###########################################################################

###########################################################################
# OLD STUFF TO REMOVE -----------------------------------------------------

# # 1 - Info interessanti sui crimini e i boundaries che ho trovato
# 
# sf_la <- sf_counties[grep("Los Angeles County", sf_counties$Name), ]
# sf_la_3857 <- st_transform(sf_la, 3857)
# la_bbox <- unname(st_bbox(st_transform(sf_la, 4326)))
# # Crop islands
# la_bbox[2] <- la_bbox[4] - (la_bbox[3] - la_bbox[1])
# la_map <- sf_ggmap(get_map(la_bbox, maptype = "stamen_terrain", source = "stadia", crop = F))
# 
# # Select W for LA county
# selected_rowcols <- which(sf_counties$Name %in% sf_la$Name)
# W_la <- W[selected_rowcols, ]
# W_la <- W_la[, selected_rowcols]
# 
# # Compute not admissible edges
# Eadj_la <- which(W_la == 1, arr.ind = TRUE)
# Eb_la <- which(W_la == 0, arr.ind = TRUE)
# 
# # Select plinks for LA county
# plinks_la <- plinks[selected_rowcols, ]
# plinks_la <- plinks_la[, selected_rowcols]
# 
# # Compute neighbouring graph for LA
# Gn_la <- matrix(0, nrow(plinks_la), ncol(plinks_la))
# Gn_la[Eadj_la] <- ifelse(plinks_la[Eadj_la] >= 0.5, 1, 0)
# 
# # Compute boundary graph for LA
# Gb_la <- matrix(0, nrow(plinks_la), ncol(plinks_la))
# Gb_la[Eadj_la] <- ifelse(plinks_la[Eadj_la] < 0.5, 1, 0)
# 
# # Compute boundary adjacency list and geometry for LA
# bound_list_la <- spdep::mat2listw(Gb_la, style = "B")$neighbours
# bound_sf_la <- boundary_geometry(bound_list_la, sf_la)
# bound_sf_la_3857 <- st_transform(bound_sf_la, 3857)
# 
# # Drop zeros
# plinks_la[Eb_la] <- NA
# Gn_la[which(Gn_la == 0, arr.ind=T)] <- NA
# Gb_la[which(Gb_la == 0, arr.ind=T)] <- NA
# 
# 
# df <- read.csv("data/explain_boundaries/2020_LACounty_Crimes_dataset.csv")
# sf <- st_as_sf(df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
#   st_transform(3857) %>% mutate(CATEGORY = as.factor(CATEGORY))
# 
# # plt_boundaries_mean <- ggplot() +
# #   geom_sf(data = sf_counties_3857, aes(fill=post_mean), col='gray25', alpha = 0.6, inherit.aes = F) +
# #   scale_fill_gradient(low = 'steelblue', high = 'darkorange') +
# #   geom_sf(data = bound_sf, col='darkred', linewidth = 0.3, inherit.aes = F) +
# #   theme_void() + theme(legend.position = "bottom") + 
# #   guides(fill = guide_colourbar(title = "Post. Mean", direction = "horizontal", barwidth = unit(2.5, "in"),
# #                                 title.position = "bottom", title.hjust = 0.5))
# 
# categories <- c("DISORDERLY CONDUCT", "LIQUOR LAWS", "RECEIVING STOLEN PROPERTY", "VAGRANCY")
# sf_crimes <- sf[sf$CATEGORY %in% categories, ]
# 
# plt_crimescat <- ggmap(la_map) +
#   geom_sf(data = sf_la_3857, col='gray25', linewidth=0.25, fill='white', alpha=0.6, inherit.aes=F) +
#   geom_sf(data = sf_crimes, aes(color=CATEGORY), size=0.75, alpha = 0.5, inherit.aes=F) +
#   geom_sf(data = bound_sf_la_3857, col='darkred', linewidth=0.5, inherit.aes=F) +
#   scale_color_manual(NULL,
#                      values = c("steelblue", 'darkorange', 'lightgreen', 'salmon'),
#                      labels = c("Disorderly Conduct", "Liquor Laws", "Recieving Stolen Property", "Vagrancy"),
#                      guide = bottom_legend(NULL,2)) +
#   theme_void() + theme(legend.position = "bottom")
#   # scale_color_manual(NULL, values = c("white",'darkred')) + guides(color = guide_legend(override.aes = list(position="none")))
#   # guides(color = bottom_legend(NULL, 2)) + theme_void() + theme(legend.position = "bottom")
# pdf("plots/plt_CrimesCategories.pdf", height = 4, width = 3.5); plt_crimescat; dev.off()
# # x11(height = 4, width = 3.5); plt_crimescat
# 
# 
# for (cat in levels(sf$CATEGORY)) {
#   print(plt_boundaries_mean + 
#     geom_sf(data = sf[sf$CATEGORY == cat, ], inherit.aes = F) +
#     geom_sf(data = bound_sf, col='darkred', linewidth = 0.5, inherit.aes = F) +
#     ggtitle(cat))
# }
# 
# plt_gangcrimes <- ggmap(counties_map) +
#   geom_sf(data = sf_counties_3857, col='gray25', linewidth=0.25, fill='white', alpha=0.6, inherit.aes=F) +
#   geom_sf(data = sf[sf$GANG_RELATED == "YES", ], aes(col=GANG_RELATED), inherit.aes=F) +
#   geom_sf(data = bound_sf, col='darkred', linewidth = 0.5, inherit.aes=F) +
#   guides(color = bottom_legend(NULL, 2)) + theme_void() + theme(legend.position = "bottom")
# pdf("plots/plt_GangCrimes.pdf", height = 5, width = 4); plt_gangcrimes; dev.off()
# # x11(height = 5, width = 4); plt_gangcrimes
# 
# # 2 - Percentuale di persone senza assicurazione medica
# 
# sf <- read_sf("data/explain_boundaries/Health_Insurance_(census_tract)/Health_Insurance_(census_tract).shp") %>%
#   st_transform(3857) %>%
#   group_by(csa) %>%
#   summarise(uninsured = mean(uninsure_1, na.rm = T))
# 
# bottom_colorbar <- function(title, length){
#   out <- guide_colorbar(title, direction = "horizontal", barwidth = unit(length, "in"), barheight = unit(0.15, "in"),
#                         title.position = "bottom", title.vjust = 0.5, label.vjust = 0.5)
#   return(out)
# }
# 
# plt_insurance <- ggmap(la_map) +
#   geom_sf(data = sf_la_3857, col=NA, fill='white', alpha=0.6, inherit.aes=F) +
#   geom_sf(data = sf, aes(fill=uninsured), col=NA, inherit.aes=F) +
#   geom_sf(data = sf_la_3857, col="gray25", fill=NA, linewidth=0.25, inherit.aes=F) +
#   geom_sf(data = bound_sf_la_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
#   scale_fill_gradient(low="steelblue", high="darkorange",
#                       guide = bottom_colorbar("Population w/o Health Insurance (%)", 2.5)) +
#   theme_void() + theme(legend.position = "bottom")
# pdf("plots/plt_NoInsurance.pdf", height = 4, width = 3.5); plt_insurance; dev.off()
# # x11(height = 4, width = 3.5); plt_insurance
# #png("people_without_insurance.png", height = 5, width = 5, units = "in", res = 200); plt_insurance; dev.off()
# 
# x11(height = 4, width = 7); gridExtra::grid.arrange(plt_crimescat, plt_insurance, ncol=2)
# 
# # 3 - Impatto delle incarcerazioni: numero di arresti della LAPD e della LASD nel 2020
# 
# sf <- read_sf("data/Incarceration_Impact_(neighborhood)/Incarceration_Impact__neighborhood_.shp") %>%
#   st_zm() %>% st_transform(3857)
# 
# plt_arrests <- ggmap(counties_map) +
#   geom_sf(data = sf, aes(fill=Combo_Arre), col=NA, inherit.aes=F) +
#   geom_sf(data = sf_counties_3857, col = "gray25", fill=NA, linewidth=0.25, inherit.aes=F) +
#   geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
#   scale_fill_gradient(low="steelblue", high="darkorange",
#                       guide = bottom_colorbar("N° of Total Arrests", 3)) +
#   theme_void() + theme(legend.position = "bottom")
# png("total_arrests.png", height = 5, width = 5, units = "in", res = 200); plt_arrests; dev.off()
# 
# # 4 - Potrebbe essere interessante: below Count and Percent below federal poverty level (fpl) and below 200 percent fpl.
# 
# sf <- read_sf("data/Below_Poverty_(census_tract)/Below_Poverty__census_tract_.shp") %>%
#   st_transform(3857) %>%
#   group_by(csa) %>%
#   summarise(below_poverty = mean(below_fpl_, na.rm = T),
#             below_200poverty = mean(below_20_1, na.rm = T))
#   
# plt_fpl <- ggmap(counties_map) +
#   geom_sf(data = sf, aes(fill=below_poverty), col=NA, inherit.aes=F) +
#   geom_sf(data = sf_counties_3857, col = "gray25", linewidth=0.25, fill=NA, inherit.aes=F) +
#   geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
#   scale_fill_gradient(low="steelblue", high="darkorange",
#                       guide = bottom_colorbar("People Below Federal Poverty Level (%)", 3)) +
#   theme_void() + theme(legend.position = "bottom")
# png("people_below_fpl.png", height = 5, width = 5, units = "in", res = 200); plt_fpl; dev.off()
# 
# plt_200fpl <- ggmap(counties_map) +
#   geom_sf(data = sf, aes(fill=below_200poverty), col=NA, inherit.aes=F) +
#   geom_sf(data = sf_counties_3857, col = "gray25", linewidth=0.25, fill=NA, inherit.aes=F) +
#   geom_sf(data = bound_sf_3857, col='darkred', linewidth=0.5, fill=NA, inherit.aes=F) +
#   scale_fill_gradient(low="steelblue", high="darkorange",
#                       guide = bottom_colorbar("People 200% Below Federal Poverty Level (%)", 3)) +
#   theme_void() + theme(legend.position = "bottom")
# png("people_below_200fpl.png", height = 5, width = 5, units = "in", res = 200); plt_200fpl; dev.off()

###########################################################################
