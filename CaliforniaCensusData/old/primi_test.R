# Import required packages
library("ggplot2")
library("ggmap")
library("sf")
library("sp")
library("dplyr")

# Import dataset
data <- read.csv("csv_pca/psam_p06.csv")
data <- data[which(data$PINCP > 0),]
data$PUMA <- as.character(data$PUMA)

sum_stats <- data %>%
  group_by(PUMA) %>%
  summarise(n = length(PINCP),
            n_na = sum(is.na(PINCP)),
            mean = mean(PINCP, na.rm = TRUE),
            sd = sd(PINCP, na.rm = TRUE))

sum_stats <- data %>%
  group_by(PUMA) %>%
  summarise(n = length(WAGP),
            n_na = sum(is.na(WAGP)),
            mean = mean(WAGP, na.rm = TRUE),
            sd = sd(WAGP, na.rm = TRUE))


# Define data transformation
tr <- function(x) { return(x) } #return(log(x)) }

# Import SpatialPolygonsDataFram and proper setting
spdf_usa <- as_Spatial(read_sf("us-pumas_boundaries/ipums_puma_2010.shp"))
spdf_usa <- spTransform(spdf_usa, CRS("+proj=longlat +datum=WGS84"))
spdf_usa$id <- row.names(spdf_usa)
spdf_usa$PUMA <- as.character(as.numeric(spdf_usa$PUMA))
spdf_usa@data <- left_join(spdf_usa@data, sum_stats, by = "PUMA")

# Import counties boundaries
cali_counties <- as_Spatial(read_sf("ca-county_boundaries/CA_Counties_TIGER2016.shp"))
cali_counties <- spTransform(cali_counties, CRS("+proj=longlat +datum=WGS84"))
cali_counties$id <- row.names(cali_counties)

# Show California with PUMAs and Counties + Mean view on the map
spdf_cali <- spdf_usa[which(spdf_usa$State == "California"), ]
# Datasets for ggplot
df_cali <- broom::tidy(spdf_cali) %>% left_join(spdf_cali@data, by = "id")
df_counties <- broom::tidy(cali_counties) %>% left_join(cali_counties@data, by = "id")
cnames <- aggregate(cbind(long, lat) ~ NAME, data=df_counties, FUN=function(x)mean(range(x)))
# Plot
plt_cali <- ggmap(get_map(spdf_cali@bbox, crop = FALSE)) +
  geom_polygon(data = df_cali, aes(x=long, y=lat, group=group, fill = tr(mean)), col="gray25", alpha = 0.75) +
  geom_polygon(data = df_counties, aes(x=long, y=lat, group=group), fill= NA, col="darkred") +
  geom_label(data = cnames, mapping = aes(x = long, y = lat, label = NAME)) +
  scale_fill_gradient(low = "steelblue", high = "darkorange") +
  theme_void()
# Show
x11(height = 5, width = 5); plt_cali

# LA + Ventura + Orange County - Mean view on the map
sel_county <- c("Los Angeles County", "Ventura County", "Orange County")
spdf <- spdf_cali[grep(paste(sel_county, collapse = "|"), spdf_cali$Name), ]
# Dataset for ggplot
df <- broom::tidy(spdf) %>% left_join(spdf@data, by = "id")
# Plot
plt_mean <- ggmap(get_map(spdf@bbox, crop = FALSE)) +
  geom_polygon(data = df, aes(x=long, y=lat, group=group, fill = tr(mean)), col="gray25", alpha = 0.75) +
  scale_fill_gradient(low = "steelblue", high = "darkorange") +
  theme_void()
plt_n <- ggmap(get_map(spdf@bbox, crop = FALSE)) +
  geom_polygon(data = df, aes(x=long, y=lat, group=group, fill = n), col="gray25", alpha = 0.75) +
  scale_fill_gradient(low = "steelblue", high = "darkorange") +
  theme_void()
# Show
x11(height = 5, width = 10); gridExtra::grid.arrange(grobs = list(plt_mean, plt_n), ncol = 2)

# Select subset of data in considered counties
data_sub <- data[which(data$PUMA %in% spdf$PUMA), ]

plt_facet <- ggplot(data = data_sub, aes(x = log(WAGP), col = PUMA, fill = PUMA)) +
  geom_density(alpha = 0.5, show.legend = FALSE) + facet_wrap(~ PUMA) +
  geom_vline(xintercept = 2.5, col='darkred', lty=1) +
  geom_vline(xintercept = 5, col='darkred', lty=1) + 
  geom_vline(xintercept = 7.5, col='darkred', lty=1) +
  geom_vline(xintercept = 10, col='darkred', lty=1)
x11(height = 10, width = 10); plt_facet

# Suitable transformation
data_sub$PINCP <- tr(data_sub$PINCP)

# Compute residual kde density in each PUMA
Npoints <- 1000
mean_dens <- density(data_sub$PINCP, n = Npoints)
res_kde_dens <- sapply(unique(data_sub$PUMA),
                       function(x) density(data_sub[which(data_sub$PUMA == x), "PINCP"], n = Npoints)$y - mean_dens$y)
rownames(res_kde_dens) <- mean_dens$x

# Dataframe for ggplot
df <- as.data.frame.table(res_kde_dens)
names(df) <- c("tr(PINCP)", "PUMA", "Res. Density")
df$`tr(PINCP)` <- as.numeric(levels(df$`tr(PINCP)`))
df$PUMA <- as.character(df$PUMA)
df <- left_join(df, spdf[, c("PUMA", "Name")], by = "PUMA", copy = TRUE)
# Plot
plt_resdens <- ggplot(df, aes(x = `tr(PINCP)`, y = `Res. Density`, group = PUMA)) +
  geom_line( col = 'steelblue', alpha = 0.8, show.legend = FALSE) + 
  geom_abline(intercept = 0, slope = 0, col='darkorange', lwd=1.2, lty=1)
# View
x11(height = 5, width = 5); plt_resdens


# Diamo ora un'occhiata alla baia di san francisco e dintroni
sel_county <- c("San Francisco County", "San Mateo County", "Santa Clara County",
                "Alameda County", "Contra Costa County", "Solano County", "Napa County",
                "Sonoma County", "Marin County", "Sacramento County")
spdf <- spdf_cali[grep(paste(sel_county, collapse = "|"), spdf_cali$Name), ]
# Dataset for ggplot
df <- broom::tidy(spdf) %>% left_join(spdf@data, by = "id")
# Plot
plt_mean <- ggmap(get_map(spdf@bbox, crop = FALSE)) +
  geom_polygon(data = df, aes(x=long, y=lat, group=group, fill = tr(mean)), col="gray25", alpha = 0.75) +
  scale_fill_gradient(low = "steelblue", high = "darkorange") +
  theme_void()
plt_n <- ggmap(get_map(spdf@bbox, crop = FALSE)) +
  geom_polygon(data = df, aes(x=long, y=lat, group=group, fill = n), col="gray25", alpha = 0.75) +
  scale_fill_gradient(low = "steelblue", high = "darkorange") +
  theme_void()
# Show
x11(height = 5, width = 10); gridExtra::grid.arrange(grobs = list(plt_mean, plt_n), ncol = 2)

# Select subset of data in considered counties
data_sub <- data[which(data$PUMA %in% spdf$PUMA), ]

# Suitable transformation
data_sub$PINCP <- tr(data_sub$PINCP)

# Compute residual kde density in each PUMA
Npoints <- 1000
mean_dens <- density(data_sub$PINCP, n = Npoints)
res_kde_dens <- sapply(unique(data_sub$PUMA),
                       function(x) density(data_sub[which(data_sub$PUMA == x), "PINCP"], n = Npoints)$y - mean_dens$y)
rownames(res_kde_dens) <- mean_dens$x

# Dataframe for ggplot
df <- as.data.frame.table(res_kde_dens)
names(df) <- c("tr(PINCP)", "PUMA", "Res. Density")
df$`tr(PINCP)` <- as.numeric(levels(df$`tr(PINCP)`))
df$PUMA <- as.character(df$PUMA)
df <- left_join(df, spdf[, c("PUMA", "Name")], by = "PUMA", copy = TRUE)
# Plot
plt_resdens <- ggplot(df, aes(x = `tr(PINCP)`, y = `Res. Density`, group = PUMA)) +
  geom_line( col = 'steelblue', alpha = 0.8, show.legend = FALSE) + 
  geom_abline(intercept = 0, slope = 0, col='darkorange', lwd=1.2, lty=1)
# View
x11(height = 5, width = 5); plt_resdens


