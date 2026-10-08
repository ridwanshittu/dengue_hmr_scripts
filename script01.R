# Dengue Modeling Pipeline - Module 01: Data Cleaning and Preparation ----
# Author: Ridwan A. Shittu
  
#### load package ----
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(blockCV)
  library(SDMtune)
  library(dplyr)
  library(sdm)
  library(fs)
  library(kuenm)
  library(ENMeval)
  library(rJava)
  library(readxl)
  library(dplyr)
  library(CoordinateCleaner)
  library(rnaturalearthdata)
  library(gtools)
  library(geodata)
  library(spThin)
  library(corrplot)
  library(usdm)
})

#set paths and seed ----

set.seed(24)
p.results <- "H:/Dengue_bioclimatic/results"
#dir.create(p.results, showWarnings = TRUE, recursive = FALSE)
p.dbio <- "H:/Dengue_bioclimatic/data"

# LOAD AND CLEAN DENGUE OCCURRENCE DATA
# Load datasets
Ms19_pts <- read.csv(file.path(p.dbio, "occurrence_records/DENGV/Mesina_2019/points_standard_checkedEC.csv"))
Ms19_ply <- read.csv(file.path(p.dbio, "occurrence_records/DENGV/Mesina_2019/poly_standard_checkedEC.csv"))
bio_df <- read_excel(file.path(p.dbio, "occurrence_records/DENGV/biogeo/DENV_2018-11-22.xlsx"), sheet = 2)
bio_df20 <- read.csv(file.path(p.dbio, "Join_thinned_data/dengue_geocode/dengue_geocode_ctrycode.csv"))
bio_df20 <- subset(bio_df20, Status != "imported")

opendengue_org_cleaned_path<- file.path(p.dbio,"opendengue/OpenDengue_geocoded_gt1990_cleaned.csv")
org_df <- read.csv(opendengue_org_cleaned_path)


# Standardize and clean
bio_df$longitude <- as.numeric(bio_df$longitude)
bio_df$latitude <- as.numeric(bio_df$latitude)

clean_df <- function(df, lon_col, lat_col) {
  df <- df[!is.na(df[[lon_col]]) & !is.na(df[[lat_col]]), ]
  df <- df[!duplicated(df[, c(lon_col, lat_col)]), ]
  return(df)
}

ms19pt <- clean_df(Ms19_pts, "Longitude", "Latitude")
ms19ply <- clean_df(Ms19_ply, "Longitude", "Latitude")
bio_df <- clean_df(bio_df, "longitude", "latitude")
bio_df20 <- clean_df(bio_df20, "lon", "lat")
org_df <- clean_df(org_df, "lon", "lat")

# Merge all datasets into one unified data frame
dengv_ds <- data.frame(
  dataset = c(
    rep("Messina Points", nrow(ms19pt)),
    rep("Messina Poly", nrow(ms19ply)),
    rep("biogeo", nrow(bio_df)),
    rep("biogeo2019 above", nrow(bio_df20)),
    rep("dengue.org", nrow(org_df))
  ),
  year = c(ms19pt$Year, ms19ply$Year, bio_df$first_year, bio_df20$Year, org_df$Year),
  longitude = c(ms19pt$Longitude, ms19ply$Longitude, bio_df$longitude, bio_df20$lon, org_df$lon),
  latitude = c(ms19pt$Latitude, ms19ply$Latitude, bio_df$latitude, bio_df20$lat, org_df$lat)
)

dim(dengv_ds) # 15206     4

# Filter for records from 1970 onwards
dengv_ds <- dengv_ds[!duplicated(dengv_ds[, c("longitude", "latitude")]), ]
dim(dengv_ds) # 15080     4
dengv_ds1970 <- subset(dengv_ds, year >= 1970)

dengv_ds1990 <- subset(dengv_ds, year>=1990)
dim(dengv_ds1990) ####   14688     4
table(dengv_ds1990$year)
dengv_by_year <- data.frame(table(dengv_ds1990$year))
names(dengv_by_year) <- c("name", "unique_count")

### remove bad Australia records ----
dg_rem <- read.csv(file.path(p.dbio, "Qgis/04_records_to_remove.csv"))
dg_rem_coord <- dg_rem[, c("longitude", "latitude")]
dengv <- dengv_ds1990[!(paste(dengv_ds1990$longitude, dengv_ds1990$latitude) %in%
                          paste(dg_rem_coord$longitude, dg_rem_coord$latitude)), ]

dim(dengv)  ###14684     4 # remove bad records in Australia

### spatially join with country shapefile  ----
world_sf <- sf::st_read(file.path(p.dbio, "NaturalEarth/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp"))
world <- terra::vect(file.path(p.dbio, "NaturalEarth/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp"))
dengv_points <- terra::vect(dengv, geom = c("longitude", "latitude"), crs = crs(world))
dengv_points <- terra::project(dengv_points, crs(world))
dengv_data <- terra::extract(world, dengv_points)
head(dengv_data)
dim(dengv_data) ###  14684   169

#
dengv_iso3 <- cbind(dengv, dengv_data[, -1])
dengv_iso3 <- dengv_iso3[!is.na(dengv_iso3$ISO_A3), ]

head(dengv_iso3)
dim(dengv_iso3) ###  14378   172

# coordinate cleaning ----

dengv_iso3$id <- 1:nrow(dengv_iso3)
dengv_iso3$species <- "dengue"
flags <- clean_coordinates(
  x = dengv_iso3,
  lon = "longitude", lat = "latitude",
  countries = "ISO_A3", country_refcol = "ISO_A3_EH",
  country_ref = world_sf, species = "species",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros", "countries")
)

#Exclude problematic records
dengv_cl <- dengv_iso3[flags$.summary, ]
dengv_cl$Country <- dengv_cl$FORMAL_EN
dim(dengv_cl) ## 13510   175

#The flagged records
dengv_fl <- dengv_iso3[!flags$.summary,]
dim(dengv_fl) #868 174

summary(dengv_cl$year)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1990    2001    2007    2007    2012    2024  

# Output cleaned data
dengue_cleaned_path <- file.path(p.results, "Dengue_cleaned_1990_2024_clean2025.csv")
dengue_flagged_path <- file.path(p.results, "Dengue_cleaned_1990_2024_flagged2025.csv")
write.csv(dengv_cl, dengue_cleaned_path, row.names = FALSE)
write.csv(dengv_fl, dengue_flagged_path, row.names = FALSE)

## bioclimatic Data Preparation -----

bio_curr_path <- file.path(p.dbio,"climate_data/BIO_1990_2020_climatology_renamed2.tif")
bio_curr <- terra::rast(bio_curr_path)
names(bio_curr)

# Remove cell duplicates in occurrence ----

dengue_cleaned_path <- file.path(p.results, "Dengue_cleaned_1990_2024_clean2025.csv")
dengv_cl <- read.csv(dengue_cleaned_path)

dengue_cells <- terra::extract(bio_curr[[1]], dengv_cl[, c("longitude", "latitude")], cellnumbers = TRUE, ID = FALSE)
cell_duplicates <- duplicated(dengue_cells[, 1])
dengv_cl_unique <- dengv_cl[!cell_duplicates, ]
dim(dengv_cl_unique) #6431  175

# Save cleaned unique dataset
dengue_unique_path <- file.path(p.results, "Dengue_cleaned_unique_cells2025.csv")
write.csv(dengv_cl_unique, dengue_unique_path, row.names = FALSE)

## Spatial Thinning ----
gc()
# Run spatial thinning 
# Try a sequence of thinning distances to identify optimal reduction
spThin.pars <- seq(from = 5, to = 20, by = 1)
thinned_results <- list()

for (i in seq_along(spThin.pars)) {
  set.seed(42)
  thinned_output <- thin(
    loc.data = dengv_cl_unique,
    lat.col = "latitude",
    long.col = "longitude",
    spec.col = "species",
    thin.par = spThin.pars[i],
    reps = 1,
    locs.thinned.list.return = TRUE,
    write.files = FALSE
  )
  thinned_results[[i]] <- thinned_output[[1]]
  cat("Thinning distance:", spThin.pars[i], "- Records retained:", nrow(thinned_output[[1]]), "\n")
}

### finalize thinning -select optimal ----

# Choose the result from thin.par = 20 based on visual/manual assessment
final_thinned <- thinned_results[[which(spThin.pars == 20)]]
final_thinned$isin <- 1

# Merge back with original dataset to annotate retained records
dengv_thinned_merge <- merge(
  x = dengv_cl_unique,
  y = final_thinned,
  by.x = c("longitude", "latitude"),
  by.y = c("Longitude", "Latitude"),
  all.x = TRUE
)

dengv_thinned <- subset(dengv_thinned_merge, !is.na(dengv_thinned_merge$isin))

dim(dengv_thinned)
#  4007  176

### Save thinned dataset ----
thinned_merge_csv <- file.path(p.results, "Dengue_cleaned_thinned_merge_2025.csv")
thinned_csv <-  file.path(p.results,"Dengv_cleaned_thinned_2025.csv")
thinned_rds <-  file.path(p.results,"Dengv_cleaned_thinned_2025.RData")

write.csv(dengv_thinned_merge, thinned_merge_csv, row.names = FALSE)
write.csv(dengv_thinned, thinned_csv, row.names = FALSE)
save(dengv_thinned, file = thinned_rds)

#### for comparison
dim(dengv_cl_unique) # 6431  175
dim(dengv_thinned_merge) # 6431  175
dim(dengv_thinned) # 4007  176

.jinit(parameters = "-Xmx56g")

### Background sampling and masking ----

pres_points <- terra::vect(dengv_thinned, geom = c("longitude", "latitude"),
                           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

buffer_area <- terra::buffer(pres_points, width = 1500000)  # 1000km
buffer_union <- terra::aggregate(buffer_area)

# Crop the bioironment rasters to buffer extent
bio_buff <- terra::crop(bio_curr, buffer_union, mask = TRUE)
bio_buff_pth <- file.path(p.results, "dengv_bio_buffer_2025.tif")
writeRaster(bio_buff, bio_buff_pth, overwrite=TRUE)

# CREATE BACKGROUND MASK (EXCLUDE PRESENCE CELLS)

grayscale <- bio_curr[[1]]
values(grayscale)[!is.na(values(grayscale))] <- 1
bio_mask <- terra::mask(grayscale, buffer_union)
#plot(bio_mask)
bio_mask_pth <- file.path(p.results, "dengv_bio_mask_buffered.tif")
writeRaster(bio_mask, bio_mask_pth, overwrite=TRUE)

# Remove presence cells
pres_cells <- terra::extract(bio_mask, pres_points, cells = TRUE)$cell
values(bio_mask)[pres_cells] <- NA
plot(bio_mask)

# SAMPLE BACKGROUND POINTS

set.seed(24)
bg_points <- terra::spatSample(bio_mask, size = 10000, method = "random", na.rm = TRUE,
                               as.points = TRUE, exhaustive = TRUE)
dim(bg_points)
# Save sampled background mask and points
mask_path <- file.path(p.results,"bio_background_mask_2025.tif")

bg_points_df_csv <- file.path(p.results,"background_points_2025.csv")
bg_points_df <- data.frame(terra::geom(bg_points)[, c("x", "y")])
colnames(bg_points_df) <- c("longitude", "latitude")

terra::writeRaster(bio_mask, filename = mask_path, overwrite = TRUE)
write.csv(bg_points_df, file = bg_points_df_csv, row.names = FALSE)


# CHECK FOR MULTICOLLINEARITY (AFTER EXPERT FILTERING)

# Extract bioclimatic values for occurrence points
dengv_coords <- data.frame(dengv_thinned[, c("longitude", "latitude")])
occs.bio <- terra::extract(bio_curr, dengv_coords, ID = FALSE)
dengv_coords <- dengv_coords[complete.cases(occs.bio), ]
occs.bio <- occs.bio[complete.cases(occs.bio), ]
dim(occs.bio)

dengv_coords_complete <- file.path(p.results, "dengv_coords_complete_2025.csv")
write.csv(dengv_coords, file = dengv_coords_complete, row.names = FALSE)



