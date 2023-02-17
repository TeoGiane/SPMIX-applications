# Import required packages
library("sf")

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

# Save data in file
save(data, file = "data/clean_data.dat")

# Compute adjacency matrix
adj_list <- spdep::poly2nb(sf_counties, row.names = pumas, queen = FALSE)
W <- spdep::nb2mat(adj_list, style = "B")

# Save W in file
save(W, file = "data/adj_matrix.dat")