# Dengue Modeling Pipeline - Module 02: Model tunning and predictions ----
# Author: Ridwan A. Shittu

# load library

library(ENMeval)
library(tidyverse)
library(sdm)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(stringr)

set.seed(24)
p.results <- "H:/Dengue_manuscript_2026/results"
# filter the lon lat column
dengv_coords <- data.frame(dengv_thinned[, c("longitude", "latitude")])

#First we randomly split the dengue occurrence data into two, training set (cross validation approach) and held-out test set.

train_idx <- sample(seq_len(nrow(dengv_coords)), size=round(0.8*nrow(dengv_coords)))
dg_train <- dengv_coords[train_idx, ]
dg_test <- dengv_coords[-train_idx,]
dim(dg_test)
dim(dg_train)

#create a 10 folds cross validation on the train datasets.
occs.bio <- terra::extract(bio_curr, dg_train, xy=TRUE,ID = FALSE)
dim(occs.bio)
occs.bio <- occs.bio[complete.cases(occs.bio), ]
dim(occs.bio)
colnames(occs.bio)
colnames(occs.bio)[20:21] <- c("longitude", "latitude")

# load the buffered rasterfiles
bio_buff_pth <- file.path(p.results, "dengv_bio_buffer_2025.tif")
bio_buff <- rast(bio_buff_pth)
#plot(bio_buff[[1]])

dg_train_points <- terra::vect(dg_train, geom = c("longitude", "latitude"),
                           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Remove presence cells
bio_mask_pth <- file.path(p.results, "dengv_bio_mask_buffered.tif")
bio_mask <- rast(bio_mask_pth)

bg_points <- terra::spatSample(bio_mask, size = 10000, method = "random", na.rm = TRUE,
                               as.points = TRUE, exhaustive = TRUE)

dim(bg_points)
dg_train_bg <- data.frame(terra::geom(bg_points)[, c("x", "y")])

bg.bio <- terra::extract(bio_curr, dg_train_bg , xy=TRUE, ID = FALSE)
bg.bio <- bg.bio[complete.cases(bg.bio), ]
dim(bg.bio)
colnames(bg.bio)
colnames(bg.bio)[20:21] <- c("longitude", "latitude")

# Combine and clean
occbg_bio <- rbind(occs.bio, bg.bio)
dim(occbg_bio) 
occbg_bio <- na.omit(occbg_bio)
dim(occbg_bio)

# Expert filtering before VIF 
vars_to_remove <- c("bio_2", "bio_3", "bio_4", "bio_8", "bio_9",
                    "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "longitude", "latitude")
occbg_bio_filtered <- occbg_bio[, !(names(occbg_bio) %in% vars_to_remove)]
dim(occbg_bio_filtered) 

# Correlation and VIF check
cor_mat <- cor(occbg_bio_filtered, method = "spearman")
corrplot(cor_mat, method = "number", col = "black", cl.pos = "n", type = "lower", tl.cex = 0.6)

#for all bioironmental variable
vif_bio <- vifcor(occbg_bio_filtered, th = 0.7, method='spearman', keep = c("bio_1", "bio_12"))

# Final selected variables (based on expert + VIF)
sel_bio <- c("bio_1", "bio_5", "bio_12", "bio_19")

p_coords <- occs.bio[, c("longitude", "latitude")]
bg_coords <- bg.bio[, c("longitude", "latitude")]

pts <- rbind(
  data.frame(longitude =p_coords$longitude, latitude=p_coords$latitude, pa=1),
  data.frame(longitude = bg_coords$longitude, latitude=bg_coords$latitude, pa=0)
)

### BlockCV Spatial partitioning ----
pts_sf <- sf::st_as_sf(pts, coords=c("longitude", "latitude"), crs = 4326, remove = FALSE)

bio_buff_sel <- bio_buff[[sel_bio]]

# Ensure CRS match
crs(pts_sf)
crs(bio_buff_sel)

stopifnot(!is.na(terra::crs(bio_buff_sel)))
pts_sf <- sf::st_transform(pts_sf, terra::crs(bio_buff_sel))

set.seed(24)
#create 10 spatial blocks
folds10 <- blockCV::cv_spatial(
  x=pts_sf,
  r = bio_buff_sel,
  column = "pa",
  k=10,
  rows_cols = c(12,12),
  selection = 'random',
  hexagon = FALSE,
  extend = 0,
  iteration = 100
  
) 

occs.grp <- folds10$folds_ids[1:nrow(p_coords)]
bg.grp <- folds10$folds_ids[(nrow(p_coords)+1):length(folds10$folds_ids)]

user.grp <- list(occs.grp = occs.grp, bg.grp = bg.grp)

cat("\nENMeval training folds (1..10) presence counts:\n")
print(table(user.grp$occs.grp, useNA = "ifany"))
print(table(user.grp$bg.grp, useNA = "ifany"))

# Maxent
pROC <- function(vars) {
  pROC <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(pROC_auc_ratio = pROC$pROC_summary[1], 
                    pROC_pval = pROC$pROC_summary[2], row.names = NULL)
  return(out)
}

cat("\n--- Tuning Maxent (ENMeval + spatial user partitions) ---\n")

tune.args <- list(
  fc = c( "LQ", "LQP", "LQH"),
  rm = 1:10
)

gc()

### Run ENMeval ----
.jinit(parameters = "-Xmx56g")

sel_bio_lonlat <- c("longitude", "latitude", sel_bio)

e_mx_bio10<- ENMevaluate(
  occs=occs.bio[sel_bio_lonlat],  
  bg=bg.bio[sel_bio_lonlat],
  algorithm = "maxent.jar",
  partitions = "user",
  user.grp = user.grp,
  tune.args = tune.args,
  user.eval = pROC,
  numCores = 8) 

saveRDS(e_mx_bio10, file.path(p.results, "ENMeval_bioclimatic_spatialfolds10_split10kbg_2026.rds"))
res_mx_bio <- eval.results(e_mx_bio10)

# Choose best model 
# res_bio <- res_mx_bio %>%
#   filter(delta.AICc == min(delta.AICc)) %>%
#   filter(pROC_pval.avg == min(pROC_pval.avg)) %>%
#   filter(or.10p.avg == min(or.10p.avg))
# res_bio

res_bio <- res_mx_bio %>%
  filter(pROC_pval.avg == min(pROC_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg)) %>%
  filter(delta.AICc == min(delta.AICc))
res_bio

### Final model using Sdm package
train_pres_df <- data.frame(pa = 1, occs.bio[sel_bio_lonlat])
train_bg_df   <- data.frame(pa = 0, bg.bio[sel_bio_lonlat])

sdm_train <- rbind(train_pres_df, train_bg_df)
sdm_train <- sdm_train[complete.cases(sdm_train), , drop = FALSE]

#prepare test data
test.bio <- terra::extract(bio_curr[[sel_bio]], dg_test, xy=TRUE,ID = FALSE)
dim(test.bio)
test.bio <- test.bio[complete.cases(test.bio), ]
colnames(test.bio)
colnames(test.bio)[5:6] <- c("longitude", "latitude")
test_pres_df <-  test.bio[, c("longitude", "latitude", sel_bio)]
head(test_pres_df)
colnames(test_pres_df)

set.seed(20)
test_bg_points <- terra::spatSample(bio_buff_sel, size = 5000, method = "random", na.rm = TRUE,
                                    xy = TRUE, as.df=TRUE, exhaustive = TRUE)

colnames(test_bg_points)
colnames(test_bg_points)[1:2] <- c("longitude", "latitude")

test_pres_df <- data.frame(pa = 1, test.bio[sel_bio_lonlat])
test_bg_df <- data.frame(pa = 0, test_bg_points[sel_bio_lonlat])

sdm_test <- rbind(test_pres_df, test_bg_df)
sdm_test <- sdm_test[complete.cases(sdm_test), , drop = FALSE]

# sdm object with test_data independent
sdm_data <- sdmData(pa ~ .+coords(longitude+latitude), train = sdm_train, test = sdm_test)

sdm_data
gc()

m_bio_sdm <- sdm(
  pa ~ .,
  data = sdm_data,
  methods = c("maxent"),
  #replication=c('boot'),n=100,
  modelSettings = list(
    maxent =list(betamultiplier=1, hinge=TRUE, threshold=FALSE, 
                 linear=TRUE, quadratic=TRUE,
                 product=FALSE, 
                 responsecurves=TRUE,
                 maximumiterations=5000)
  ),
  parallelSetting = list(ncore = 8)
)

m_bio_sdm@setting
m_bio_sdm@models


#### plot variable importance and rcurves
sdm::rcurve(m_bio_sdm)
sdm::roc(x=m_bio_sdm)
plot(getVarImp(m_bio_sdm))
getModelInfo(m_bio_sdm)

### plot out the model graphs
# Open the PNG graphics device and define parameters

png(filename = "H:/Dengue_manuscript_2026/figures/FigS1d.png", width = 12.8, height = 8.09, units = "in", res = 600)

#  Execute plotting code
sdm::roc(x=m_bio_sdm)

# Close the device to finalize and save the file
dev.off()

png(filename = "H:/Dengue_manuscript_2026/figures/FigS1e.png", width = 12.8, height = 8.09, units = "in", res = 600)

plot(getVarImp(m_bio_sdm), 'auc', col="#feb351")

dev.off()

png(filename = "H:/Dengue_manuscript_2026/figures/FigS1f.png", width = 12.8, height = 8.09, units = "in", res = 600)

sdm::rcurve(m_bio_sdm, gg=FALSE)

dev.off()

#prediction on buffered area
pred_buff <- predict(
  object   = m_bio_sdm,
  newdata  = bio_buff_sel,
  parallelSetting = list(ncore = 8)
)

plot(pred_buff)

#prediction on global area 
pred_glo <- predict(
  object   = m_bio_sdm,
  newdata  = bio_curr[[sel_bio]],
  parallelSetting = list(ncore = 8))

plot(pred_glo)

out_glo <- file.path(p.results,"baseline")
dir.create(out_glo, recursive = TRUE)

writeRaster(
  pred_glo,
  names="baseline2020",
  filename = file.path(out_glo, "prediction_global_baseline.tif" ),
  overwrite = TRUE
)


#### future predictions ----
fut_clim_dir <- "H:/future_climate"
out_fut <- file.path(p.results, "fut_predictions")
dir.create(out_fut, recursive = TRUE)
years <- c("2030", "2050", "2070")
ssps  <- c("ssp126", "ssp245", "ssp370", "ssp585")
mask_r <- rast("H:/DengueENM/results/maxent/map_mask.tif")
plot(mask_r)

#sel_bio

bionames <- paste0("bio_", 1:19)
ncore = 8 #for parallel computing

for (yr in years) {
  for (ssp in ssps) {
    
    message("Processing ", yr, " / ", ssp)
    in_dir <- file.path(fut_clim_dir, yr, ssp)
    
    if (!dir.exists(in_dir)) next
    
    gcm_files <- list.files(in_dir, pattern = "\\.tif$", full.names = TRUE)
    if (length(gcm_files) == 0) next
    
    out_ssp_yr <- file.path(out_fut, yr, ssp)
    dir.create(out_ssp_yr, showWarnings = FALSE, recursive = TRUE)
    final_file <- file.path(out_ssp_yr, paste0(yr, "_", ssp, "_GCM_sd.tif"))
    
    if (file.exists(final_file)) {
      message("Skipping ", yr, " / ", ssp, " (already completed)")
      next
    }
    
    gcm_means <- list()
    
    
    # loop over GCMs
    
    for (gcm_file in gcm_files) {
      try({
        
        gcm_name <- sub(
          "wc2\\.1_2\\.5m_bioc_|_ssp[0-9]+_.*\\.tif",
          "",
          basename(gcm_file)
        )
        
        message("  → GCM: ", gcm_name)
        
        bio_future <- rast(gcm_file)
        
        # ensure standard WorldClim names
        if (nlyr(bio_future) != 19) {
          warning("Skipping ", gcm_name, ": expected 19 bioclim layers")
          next
        }
        names(bio_future) <- bionames
        
        # Now safely subset
        bio_future <- bio_future[[sel_bio]]
        
        # apply mask
        bio_future <- terra::crop(bio_future, mask_r, mask=TRUE)
        
        out_gcm <- file.path(out_ssp_yr, gcm_name)
        dir.create(out_gcm, showWarnings = FALSE, recursive = TRUE)
        
        # raw individual model prediction
        
        pred_fut <- predict(
          object   = m_bio_sdm,
          newdata  = bio_future,
          overwrite = TRUE,
          parallelSetting = list(ncore = ncore)
        )
        
        writeRaster(
          pred_fut,
          filename = file.path(out_gcm, paste0(gcm_name, "_prediction.tif")),
          overwrite = TRUE
        )
        
        gcm_means[[gcm_name]] <- pred_fut
        
        
      })
    }
    
    # SSP–YEAR aggregation (across GCMs)
    
    if (length(gcm_means) > 0) {
      
      gcm_stack <- rast(gcm_means)
      
      mean(
        gcm_stack,
        na.rm = TRUE,
        filename = file.path(out_ssp_yr, paste0(yr, "_", ssp, "_GCM_mean.tif")),
        overwrite = TRUE
      )
      
      app(
        gcm_stack,
        median,
        na.rm = TRUE,
        filename = file.path(out_ssp_yr, paste0(yr, "_", ssp, "_GCM_median.tif")),
        overwrite = TRUE
      )
      
      stdev(
        gcm_stack,
        na.rm = TRUE,
        filename = file.path(out_ssp_yr, paste0(yr, "_", ssp, "_GCM_sd.tif")),
        overwrite = TRUE
      )
    }
  }
}


#create binary files
binary_dir <- file.path(p.results, "binary_maps")
change_dir <- file.path(p.results, "change_maps")
dir.create(binary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(change_dir, showWarnings = FALSE, recursive = TRUE)

# Occurrence data (all occurrence records to compute thresholds)
dengv_coords_complete <- file.path(p.results, "dengv_coords_complete_2025.csv")
p_coords <- read.csv(dengv_coords_complete)
any(is.na(p_coords))

occ <- p_coords

#Threshold function (MTP, 5%, 10%) ----
sdm_threshold <- function(suit_raster, occ_xy, type = "mtp") {
  
  occ_vals <- terra::extract(suit_raster, occ_xy)[,2]
  occ_vals <- na.omit(occ_vals)
  
  if (length(occ_vals) == 0)
    stop("No occurrence values extracted — check CRS and extent.")
  
  if (type == "mtp") {
    return(min(occ_vals))
  }
  
  if (type == "p05") {
    return(quantile(occ_vals, probs = 0.05))
  }
  
  if (type == "p10") {
    return(quantile(occ_vals, probs = 0.10))
  }
  
  stop("Unknown threshold type")
}

# Compute thresholds  ----
#Use current ensemble mean.
# Load current ensemble
ens_glo_bio <- rast(file.path(out_glo, "prediction_global_baseline.tif"))

# Compute thresholds
th_mtp <- sdm_threshold(ens_glo_bio, occ, "mtp")
th_p05 <- sdm_threshold(ens_glo_bio, occ, "p05")
th_p10 <- sdm_threshold(ens_glo_bio, occ, "p10")

thresholds <- c(MTP = th_mtp, P05 = th_p05, P10 = th_p10)
print(thresholds)
names(thresholds) <- c("MTP", "P05", "P10")
saveRDS(thresholds, file.path(p.results, "thresholds_current_bio.rds"))

#Compute thresholds from current predictions
#Apply thresholds to future  prediction
#We apply only to future scenario mean & median (never to SD).

apply_threshold <- function(r, th) {
  out <- r
  out[out < th] <- 0
  out[out >= th] <- 1
  return(out)
}

# apply threshold to current prediction.

glo_bio_mtp <- apply_threshold(ens_glo_bio, th_mtp)
plot(glo_bio_mtp)
glo_bio_p05 <- apply_threshold(ens_glo_bio, th_p05)
plot(glo_bio_p05)
glo_bio_p10 <- apply_threshold(ens_glo_bio, th_p10)
plot(glo_bio_p10)

# Save outputs
writeRaster(
  glo_bio_mtp,
  filename = file.path(binary_dir,"global_current_bio_MTP_binary.tif"),
  overwrite = TRUE
)

writeRaster(
  glo_bio_p05,
  filename = file.path(binary_dir,"global_current_bio_P05_binary.tif"),
  overwrite = TRUE
)

writeRaster(
  glo_bio_p10,
  filename = file.path(binary_dir,"global_current_bio_P10_binary.tif"),
  overwrite = TRUE
)

#only p05 and p10
thresholds <- c(P05 = th_p05, P10 = th_p10)
names(thresholds) <- c("P05", "P10")

#Loop through future ensemble outputs
years <- c("2030", "2050", "2070")
ssps  <- c("ssp126", "ssp245", "ssp370", "ssp585")

for (yr in years) {
  for (ssp in ssps) {
    
    message("Thresholding ", yr, " / ", ssp)
    
    ens_mean_file   <- file.path(out_fut, yr, ssp, paste0(yr, "_", ssp, "_GCM_mean.tif"))
    ens_median_file <- file.path(out_fut, yr, ssp, paste0(yr, "_", ssp, "_GCM_median.tif"))
    
    if (!file.exists(ens_mean_file)) next
    
    ens_mean   <- rast(ens_mean_file)
    ens_median <- rast(ens_median_file)
    
    # Apply thresholds
    for (th_name in names(thresholds)) {
      
      th_val <- thresholds[th_name]
      
      mean_bin <- apply_threshold(ens_mean, th_val)
      median_bin <- apply_threshold(ens_median, th_val)
      
      writeRaster(
        mean_bin,
        file.path(binary_dir, paste0(yr, "_", ssp, "_GCM_mean_bio_", th_name, "_binary.tif")),
        overwrite = TRUE
      )
      
      writeRaster(
        median_bin,
        file.path(binary_dir,  paste0(yr, "_", ssp, "_GCM_median_bio_", th_name, "_binary.tif")),
        overwrite = TRUE
      )
    }
  }
}

