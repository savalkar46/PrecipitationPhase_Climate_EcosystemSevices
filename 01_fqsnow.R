# =============================================================================
# 01_fqsnow.R   --   fQ,snow: CALCULATION (Part A) and FIGURES (Part B) in one script
## Author:      Supriya Savalkar
## Created:     2025-10-09
## Last updated: 2026-10-05
## Code cleanup and documentation assisted by Claude (Anthropic)
# PAPER SECTIONS : 3.1, 3.2; Figs 3, 4, 5, A2, A3, A4; CRB-level fQ,snow values
#
# PART A  CALCULATION  (slow; set RUN_CALCULATION <- FALSE to skip it)
# WHAT IT DOES   : Computes fQ,snow -- the fraction of streamflow that came from
#                  snowmelt -- for every VIC-CropSyst grid cell, separately for the two
#                  rain-snow partitioning methods:
#                      DRS = TRH-RSP  (temperature + relative humidity)
#                      SRS = T-RSP    (temperature only)
#                  for the historical period and (optionally) for every future
#                  GCM x RCP x timeframe.
#
# HOW fQ,snow IS CALCULATED (per cell, from the DAILY VIC-CropSyst flux file)
#   1. Daily water input is SnowMelt + Rain. p_melt and p_rain are the melt and
#      rain shares of that input (0 on days with no input).
#   2. Root-zone water is sm = sum of the 17 soil-moisture layers. Infiltration is
#      the water needed to close the soil balance:
#          infil = max(0, d(sm) + Baseflow + Evap)
#   3. f_so = fraction of soil water that originated as snowmelt. It is updated
#      every day with a mixing equation (inflow of melt-derived water, minus the
#      melt-fraction share of what leaves by Evap + Baseflow):
#          f_so[t] = ( f_so[t-1]*sm[t-1] + infil[t]*p_melt[t]
#                      - (Evap[t] + Baseflow[t]) * f_so[t-1] ) / sm[t]
#      (the original code also tracked the rain complement g_so; it never enters
#       fQ,snow, so it is not computed here -- results are identical.)
#   4. Input that does not leave as streamflow is surf_loss = melt + rain -
#      (Runoff + Baseflow), shared out pro rata:
#          tot_rsn = SnowMelt - surf_loss * p_melt      (snow-derived quick flow)
#          tot_rp  = Rain     - surf_loss * p_rain      (rain-derived quick flow)
#   5. Per water year (1 Oct - 30 Sep):
#          fQ,snow = sum(tot_rsn + Baseflow*f_so) / sum(tot_rp + tot_rsn + Baseflow)
#      rounded to 2 decimals, then averaged over the water years of the period.
#
#   Starting values: f_so = 0 on the first day of the period (no spin-up).
#
# meanRunoff     : mean annual (Runoff + Baseflow), mm per water year. It is the
#                  weight used later to aggregate cells to sub-basins (Part B).
#                  >>> CONFIRM this matches the runoff variable you used for your
#                  >>> existing gridlevel_hist_*.rds files (see README, "Verify").
#
# INPUTS         : raw VIC-CropSyst daily flux files, one per cell, no header,
#                  tab-delimited, 32 columns, named flux_<lat>_<lon>
#                  (historical: paths$vic_flux_dir("DRS"/"SRS");
#                   future: see future_flux_dir() below)
# OUTPUTS        : historical -> <paths$fq_hist_dir>/gridlevel_hist_{DRS,SRS}.rds
#                  future     -> <paths$fq_future_dir>/<GCM>_<RCP>_<DRS|SRS>.rds
#                  columns: lat_lon, lat, lon, meanfqsnow, medfqsnow, meanRunoff
#                           (+ Timeframe, Partition, GCM, RCP for the future files)
#
#
# PART B  FIGURES  (set RUN_FIGURES <- FALSE to skip it)
#   Reads the grid-level fQ,snow written by Part A, aggregates to HUC6 sub-basins (runoff-
#   weighted), clips everything to the Columbia River Basin and draws the figures:
#   Fig 3   grid-level historical fQ,snow, T-RSP vs TRH-RSP
#   Fig 4   HUC6 historical fQ,snow: T-RSP | TRH-RSP | CHANGE (TRH-RSP minus T-RSP)
#   Fig 5   dumbbell plot, Historical -> 2040s -> 2080s, RCP8.5   (A3 = RCP4.5)
#   Fig A2  grids contributing 80% of the runoff in each sub-watershed
#   Fig A4  projected change in fQ,snow by sub-basin;   crb_level_fqsnow.csv = CRB-level value
#
# Everything is limited to the CRB boundary, so the coastal watersheds are gone:
#   * grid cells are kept only if their CENTRE lies inside the CRB polygon
#   * HUC6 sub-basins are dropped if less than MIN_INSIDE of their area is inside
#     the CRB (or if listed in DROP_BASINS) and the rest are clipped to the boundary
#   * basin values are recomputed from the kept cells (switch: RECOMPUTE_FROM_CLIPPED)
#
# Also written to OUT_DIR (check these before editing the paper):
#   subbasin_inside_crb.csv        area fraction inside the CRB per sub-basin, kept/dropped
#   basin_fqsnow_all_scenarios.csv basin values used in the figures
#   basin_fqsnow_clip_effect.csv   basin values with vs without the clip (historical)
#
# Part B INPUTS  : output of Part A, HUC6 + CRB shapefiles, HUC6 lookup CSVs
# Part B OUTPUTS : outputs/fq_snow_CRB/*.png and *.csv
# =============================================================================
source("R/setup.R")
use_packages(c("dplyr", "stringr", "parallel", "tidyr", "purrr", "ggplot2", "sf", "terra",
               "patchwork", "ggtext", "scales"))
sf_use_s2(FALSE)


## ---- 0. SETTINGS -----------------------------------------------------------------
RUN_CALCULATION <- TRUE        # PART A (slow). FALSE = reuse the .rds files from an earlier run
RUN_FIGURES     <- TRUE        # PART B (figures and tables)

PARTITIONS     <- c("DRS", "SRS")
RUN_HISTORICAL <- TRUE
RUN_FUTURE     <- FALSE        # TRUE after you have filled in FUTURE_RUNS / future_flux_dir() below
MAX_CELLS      <- Inf          # set e.g. 50 for a quick test run
N_CORES        <- max(1, parallel::detectCores() - 1)    # forced to 1 on Windows

## Historical period: water years 1986-2015 (1 Oct 1985 - 30 Sep 2015)
HIST_START <- as.Date("1985-10-01")
HIST_END   <- as.Date("2015-09-30")

## Folder with the historical flux files of one partition.
hist_flux_dir <- function(part) paths$vic_flux_dir(part)

## Future runs. One row per GCM x RCP. Fill in your GCM names.
FUTURE_RUNS <- expand.grid(GCM = c("<GCM_1>", "<GCM_2>"),      # <-- replace
                           RCP = c("RCP45", "RCP85"),
                           stringsAsFactors = FALSE)
## Future timeframes (water years). >>> CONFIRM these match the windows used in the paper
## (2040s = 2026-2055, 2080s = 2065-2094).
FUTURE_PERIODS <- data.frame(Timeframe = c("2040s", "2080s"),
                             start = as.Date(c("2025-10-01", "2064-10-01")),
                             end   = as.Date(c("2055-09-30", "2094-09-30")))
## Folder with the future flux files of one partition x GCM x RCP.  <-- EDIT to your layout
future_flux_dir <- function(part, gcm, rcp)
  file.path(cfg$data_root, "05_Analysis/Q1/Future", part, paste(gcm, rcp, sep = "_"))

VIC_FLUX_COLS <- c("Year", "Month", "Day", "Precip", "Evap", "Runoff", "Baseflow", "Rain",
                   "Snow", "Snow_Depth", "SWE", "Delta_SWE", "SnowMelt", "SnowEvap", "AirTemp",
                   paste0("SM", 1:17))                          # SM1-SM17 = columns 16-32
if (.Platform$OS.type == "windows") N_CORES <- 1

## ---- PART A. fQ,snow CALCULATION -- 1. Helpers -------------------------------------------------------------------
## Water year starting 1 October (Oct-Dec belong to the next water year)
water_year <- function(dates, start_month = 10) {
  as.integer(format(dates, "%Y")) + as.integer(as.integer(format(dates, "%m")) >= start_month)
}

read_flux <- function(file) {
  d <- tryCatch(
    if (requireNamespace("data.table", quietly = TRUE))
      as.data.frame(data.table::fread(file, header = FALSE, sep = "\t"))
    else read.delim(file, header = FALSE),
    error = function(e) NULL)
  if (is.null(d) || ncol(d) != length(VIC_FLUX_COLS)) return(NULL)
  names(d) <- VIC_FLUX_COLS
  d
}

## fQ,snow of ONE grid cell for ONE period. Returns one row (or NULL if the file is unusable).
fqsnow_cell <- function(file, start_date, end_date) {
  tmp <- read_flux(file)
  if (is.null(tmp)) return(NULL)

  ll <- str_match(basename(file), "flux_([0-9.]+)_(-[0-9.]+)")
  latitude <- as.numeric(ll[2]); longitude <- as.numeric(ll[3])

  v <- tmp |>
    mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"), format = "%Y-%m-%d"),
           WaterYear = water_year(Date)) |>
    filter(Date >= start_date, Date <= end_date) |>
    arrange(Date)
  n <- nrow(v)
  if (n < 2) return(NULL)

  sm      <- rowSums(v[, paste0("SM", 1:17)])          # root-zone soil water
  inflow  <- v$SnowMelt + v$Rain
  p_melt  <- ifelse(inflow == 0, 0, v$SnowMelt / inflow)
  p_rain  <- ifelse(inflow == 0, 0, v$Rain / inflow)

  infil <- c(0, diff(sm)) + v$Baseflow + v$Evap
  infil[infil < 0] <- 0
  infil[1] <- 0

  f_so <- numeric(n)                                    # snowmelt fraction of soil water
  for (d in 2:n) {
    f_so[d] <- if (inflow[d] == 0) f_so[d - 1] else
      (f_so[d - 1] * sm[d - 1] + infil[d] * p_melt[d] - (v$Evap[d] + v$Baseflow[d]) * f_so[d - 1]) / sm[d]
  }

  surf_loss <- inflow - (v$Runoff + v$Baseflow)
  tot_rsn   <- v$SnowMelt - surf_loss * p_melt
  tot_rp    <- v$Rain     - surf_loss * p_rain

  v$numerator   <- tot_rsn + v$Baseflow * f_so
  v$denominator <- tot_rp + tot_rsn + v$Baseflow
  v$runoff_tot  <- v$Runoff + v$Baseflow

  yearly <- v |>
    group_by(WaterYear) |>
    summarise(fqsnow = round(sum(numerator) / sum(denominator), 2),
              runoff = sum(runoff_tot), .groups = "drop")

  tibble::tibble(lat_lon = paste(latitude, longitude, sep = "_"), lat = latitude, lon = longitude,
                 meanfqsnow = round(mean(yearly$fqsnow), 2),
                 medfqsnow  = median(yearly$fqsnow),
                 meanRunoff = round(mean(yearly$runoff), 2))
}

## All cells in a folder
fqsnow_dir <- function(dir, start_date, end_date, label = basename(dir)) {
  if (!dir.exists(dir)) stop("Flux folder not found: ", dir, call. = FALSE)
  files <- list.files(dir, pattern = "^flux_", full.names = TRUE)
  if (!length(files)) stop("No flux_* files in ", dir, call. = FALSE)
  if (is.finite(MAX_CELLS)) files <- files[seq_len(min(length(files), MAX_CELLS))]
  message(sprintf("  %s: %d cells", label, length(files)))
  t0 <- Sys.time()
  res <- bind_rows(parallel::mclapply(files, fqsnow_cell, start_date = start_date,
                                      end_date = end_date, mc.cores = N_CORES))
  message("    done in ", format(Sys.time() - t0))
  res
}

## ---- PART A -- 2. Historical, both partitions ---------------------------------------------------
if (RUN_CALCULATION && RUN_HISTORICAL) {
  dir.create(paths$fq_hist_dir, recursive = TRUE, showWarnings = FALSE)
  for (part in PARTITIONS) {
    message("[historical ", part, "]")
    res <- fqsnow_dir(hist_flux_dir(part), HIST_START, HIST_END, part)
    saveRDS(res, file.path(paths$fq_hist_dir, sprintf("gridlevel_hist_%s.rds", part)))
  }
}

## ---- PART A -- 3. Future, both partitions x GCM x RCP x timeframe ---------------------------------
if (RUN_CALCULATION && RUN_FUTURE) {
  dir.create(paths$fq_future_dir, recursive = TRUE, showWarnings = FALSE)
  for (part in PARTITIONS) for (r in seq_len(nrow(FUTURE_RUNS))) {
    gcm <- FUTURE_RUNS$GCM[r]; rcp <- FUTURE_RUNS$RCP[r]
    message("[future ", part, " ", gcm, " ", rcp, "]")
    res <- bind_rows(lapply(seq_len(nrow(FUTURE_PERIODS)), function(k)
      fqsnow_dir(future_flux_dir(part, gcm, rcp), FUTURE_PERIODS$start[k], FUTURE_PERIODS$end[k],
                 paste(part, gcm, rcp, FUTURE_PERIODS$Timeframe[k])) |>
        mutate(Timeframe = FUTURE_PERIODS$Timeframe[k])))
    res$Partition <- part; res$GCM <- gcm; res$RCP <- rcp
    saveRDS(res, file.path(paths$fq_future_dir, sprintf("%s_%s_%s.rds", gcm, rcp, part)))
  }
}

## =============================================================================
## PART B. FIGURES AND TABLES (CRB only)
## =============================================================================
if (RUN_FIGURES) {

## ---- PART B. FIGURES -- settings (paths come from config/config.R via R/setup.R) -----
CRB_SHP      <- cfg$crb_shp                                    # CRB boundary shapefile
HUC6_SHP     <- cfg$huc6_shp                                   # HUC6 sub-basin polygons
HUC6_LOOKUP  <- paths$huc6_lookup                              # lat_lon, lat, lon, Name
HUC6_CLASS   <- paths$huc6_class                               # Name, Classification
HIST_DIR     <- paths$fq_hist_dir                              # gridlevel_hist_{DRS,SRS}.rds   (from Part A)
FUTURE_DIR   <- paths$fq_future_dir                            # one .rds per GCM x RCP x method (from Part A)
OUT_DIR      <- out_dir("fq_snow_CRB")
MIN_INSIDE   <- 0.5                     # drop a sub-basin if < 50% of its area is inside the CRB
DROP_BASINS  <- c("WashingtonCoastal", "SouthernOregonCoastal", "NorthernOregonCoastal")  # always dropped
KEEP_BASINS  <- character(0)            # force-keep (overrides the area rule)
RECOMPUTE_FROM_CLIPPED <- TRUE          # TRUE: basin values from cells inside the CRB only
USE_CACHED_FILTER <- TRUE               # reuse the saved grid filter (delete the CSV to rebuild)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
GRID_FILTER_CSV <- file.path(OUT_DIR, "crb_grid_filter.csv")   # lat_lon, lat, lon, inside
save_png <- function(p, f, w, h) { ggsave(file.path(OUT_DIR, f), p, width = w, height = h, units = "cm", dpi = 300); message("Saved ", f) }

## ---- 1. Historical grid data (both methods) and the lookup of cells -> sub-basins ----
read_hist <- function(code, label)
  readRDS(file.path(HIST_DIR, sprintf("gridlevel_hist_%s.rds", code))) |> ungroup() |> mutate(Partition = label)
hist_grid <- bind_rows(read_hist("DRS", "TRH-RSP"), read_hist("SRS", "T-RSP")) |>
  mutate(Partition = factor(Partition, levels = c("T-RSP", "TRH-RSP")))

huc6 <- read.csv(HUC6_LOOKUP)                   # lat_lon, lat, lon, Name
huc6$Name <- gsub("\\s+", "", huc6$Name)

## ---- 2. CRB boundary and the GRID FILTER (built once, then reused) --------------------
## Rasterise the CRB polygon on the model grid and keep the cells whose CENTRE falls inside.
## This is fast (no point-in-polygon tests) and is saved to GRID_FILTER_CSV, so later runs
## just read it. Delete that file to rebuild (e.g. after changing the boundary).
crb <- st_sf(geometry = st_union(st_make_valid(st_read(CRB_SHP, quiet = TRUE))))

domain_cells <- bind_rows(distinct(hist_grid, lat_lon, lat, lon), distinct(huc6, lat_lon, lat, lon)) |>
  distinct() |> ungroup()

build_grid_filter <- function(cells) {
  xres <- median(diff(sort(unique(cells$lon)))); yres <- median(diff(sort(unique(cells$lat))))
  r <- terra::rast(xmin = min(cells$lon) - xres / 2, xmax = max(cells$lon) + xres / 2,
                   ymin = min(cells$lat) - yres / 2, ymax = max(cells$lat) + yres / 2,
                   resolution = c(xres, yres), crs = "EPSG:4326")
  m <- terra::rasterize(terra::vect(crb), r, field = 1, touches = FALSE)   # touches = FALSE: centre rule
  hit <- terra::extract(m, as.matrix(cells[, c("lon", "lat")]))[[1]]
  cells |> mutate(inside = !is.na(hit))
}

grid_filter <- if (USE_CACHED_FILTER && file.exists(GRID_FILTER_CSV)) read.csv(GRID_FILTER_CSV) else NULL
if (is.null(grid_filter) || !all(domain_cells$lat_lon %in% grid_filter$lat_lon)) {
  grid_filter <- build_grid_filter(domain_cells)
  write.csv(grid_filter, GRID_FILTER_CSV, row.names = FALSE)
}
cat("\nGrid cells in the domain:", nrow(grid_filter), " inside the CRB:", sum(grid_filter$inside), "\n")

cells_in_crb <- function(df)                    # df needs a lat_lon column
  semi_join(df, filter(grid_filter, inside), by = "lat_lon")

## ---- 3. Sub-basins that stay, clipped to the boundary -------------------------------------
subbasin_all <- st_read(HUC6_SHP, quiet = TRUE) |>
  st_make_valid() |> mutate(Name = gsub("\\s+", "", Name))

pieces <- suppressWarnings(st_intersection(subbasin_all["Name"], crb)) |>      # ONE vectorised clip
  st_make_valid()
area_km2 <- function(x) as.numeric(st_area(st_transform(x, 5070))) / 1e6      # equal-area CRS
basin_inside <- tibble(Name = subbasin_all$Name, area_km2 = area_km2(subbasin_all)) |>
  left_join(tibble(Name = pieces$Name, inside_km2 = area_km2(pieces)) |>
              group_by(Name) |> summarise(inside_km2 = sum(inside_km2)), by = "Name") |>
  mutate(inside_km2 = coalesce(inside_km2, 0),
         frac_inside = inside_km2 / area_km2,
         keep = (frac_inside >= MIN_INSIDE | Name %in% KEEP_BASINS) & !(Name %in% DROP_BASINS))
write.csv(basin_inside, file.path(OUT_DIR, "subbasin_inside_crb.csv"), row.names = FALSE)
cat("\nSub-basins vs the CRB boundary (REVIEW THIS):\n")
print(basin_inside |> mutate(across(where(is.numeric), \(x) round(x, 2))) |> arrange(frac_inside), n = Inf)
KEEP <- basin_inside$Name[basin_inside$keep]

classification <- read.csv(HUC6_CLASS)
subbasin <- pieces |> filter(Name %in% KEEP) |>
  group_by(Name) |> summarise(geometry = suppressWarnings(st_union(geometry)), .groups = "drop") |>
  left_join(distinct(classification, Name, Classification), by = "Name")

huc6_crb <- cells_in_crb(huc6) |> filter(Name %in% KEEP)
cat("Grid cells with a sub-basin: ", nrow(huc6), " -> inside CRB and kept: ", nrow(huc6_crb), "\n")

## ---- 4. Basin values (runoff-weighted) ----------------------------------------------
basin_value <- function(df, cells)              # df: grid rows with meanfqsnow, meanRunoff
  df |> inner_join(cells |> dplyr::select(lat_lon, lat, lon, Name), by = c("lat_lon", "lat", "lon")) |>
    mutate(Snow_runoff = meanfqsnow * meanRunoff)

hist_basin <- function(cells)
  basin_value(hist_grid, cells) |>
    group_by(Name, Partition) |>
    summarise(basin_fqsnow = round(sum(Snow_runoff) / sum(meanRunoff), 2), .groups = "drop") |>
    mutate(Timeframe = "Historical", RCP = "Historical", Partition = as.character(Partition))

future_files <- list.files(FUTURE_DIR, full.names = TRUE)
future_grid <- map_dfr(future_files, \(f) {
  d <- ungroup(readRDS(f))      # the saved files can be grouped data frames
  if (!"RCP" %in% names(d)) d$RCP <- str_extract(basename(f), "RCP[0-9.]+")
  d
}) |> mutate(RCP = gsub("\\.", "", RCP),
             Partition = recode(Partition, DRS = "TRH-RSP", SRS = "T-RSP"))

future_basin <- function(cells)
  future_grid |>
    group_by(lat_lon, lat, lon, Timeframe, Partition, RCP) |>
    summarise(fq = round(mean(meanfqsnow), 2), ro = round(mean(meanRunoff), 2), .groups = "drop") |>
    inner_join(cells |> dplyr::select(lat_lon, lat, lon, Name), by = c("lat_lon", "lat", "lon")) |>
    group_by(Timeframe, Partition, RCP, Name) |>
    summarise(basin_fqsnow = round(sum(fq * ro) / sum(ro), 2), .groups = "drop")

cells_use <- if (RECOMPUTE_FROM_CLIPPED) huc6_crb else huc6 |> filter(Name %in% KEEP)
basin_hist <- hist_basin(cells_use)
basin_fut  <- future_basin(cells_use) |> filter(Timeframe != "2060s")

write.csv(bind_rows(basin_hist, basin_fut), file.path(OUT_DIR, "basin_fqsnow_all_scenarios.csv"), row.names = FALSE)

## how much does the clip change the historical basin values?
effect <- full_join(hist_basin(huc6 |> filter(Name %in% KEEP)) |> rename(unclipped = basin_fqsnow),
                    hist_basin(huc6_crb) |> rename(clipped = basin_fqsnow),
                    by = c("Name", "Partition", "Timeframe", "RCP")) |>
  mutate(change = clipped - unclipped) |> arrange(desc(abs(change)))
write.csv(effect, file.path(OUT_DIR, "basin_fqsnow_clip_effect.csv"), row.names = FALSE)
cat("\nLargest changes in historical basin fQ,snow caused by the clip:\n"); print(head(effect, 8))

## ---- 5. Shared plot style -------------------------------------------------------------
fq_name <- expression(f[Q*","*snow])
theme_fq <- theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")
colourbar <- guide_colorbar(barwidth = unit(8, "cm"), barheight = unit(0.8, "cm"),
                            label.theme = element_text(size = 12, face = "bold"))
crb_outline <- geom_sf(data = crb, fill = NA, colour = "grey25", linewidth = 0.3, inherit.aes = FALSE)
lonlat_axes <- list(scale_x_continuous(labels = number_format(accuracy = 1)),
                    scale_y_continuous(labels = number_format(accuracy = 1)))

## =============================================================================
## FIGURE 3  grid-level historical fQ,snow, clipped to the CRB
## =============================================================================
g3 <- cells_in_crb(hist_grid)
res <- median(diff(sort(unique(g3$lon))))                         # grid spacing from the data
fig3 <- ggplot(g3, aes(lon, lat, fill = meanfqsnow)) +
  geom_tile(width = res, height = res) +
  scale_fill_viridis_c(option = "viridis", direction = -1, limits = c(0, 1), name = fq_name) +
  coord_sf(crs = 4326) + lonlat_axes + guides(fill = colourbar) +
  labs(x = NULL, y = NULL) + theme_fq + facet_grid(~Partition) +
  theme(legend.title = element_text(size = 20, face = "bold"))
save_png(fig3, "Fig3_fqsnow_grid_historical_CRB.png", 26, 17)

## =============================================================================
## FIGURE 4  HUC6 historical: T-RSP | TRH-RSP | CHANGE
## =============================================================================
labels_at <- function(sf_obj) {
  pts <- suppressWarnings(st_point_on_surface(sf_obj))
  bind_cols(st_drop_geometry(sf_obj), as_tibble(st_coordinates(pts)) |> rename(lon = X, lat = Y))
}

hist_sf   <- subbasin |> inner_join(basin_hist, by = "Name")
hist_lab  <- labels_at(hist_sf)
p_two <- ggplot() +
  geom_sf(data = hist_sf, aes(fill = basin_fqsnow), colour = "grey35", linewidth = 0.2) +
  geom_text(data = hist_lab, aes(lon, lat, label = round(basin_fqsnow, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis_c(option = "viridis", direction = -1, limits = c(0, 1),
                       name = expression("Basin level " * f[q*","*snow] * " value")) +
  lonlat_axes + guides(fill = colourbar) + labs(x = NULL, y = NULL) +
  theme_fq + theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  facet_wrap(~factor(Partition, levels = c("T-RSP", "TRH-RSP")))

change_sf <- subbasin |>
  inner_join(basin_hist |> dplyr::select(Name, Partition, basin_fqsnow) |>
               pivot_wider(names_from = Partition, values_from = basin_fqsnow) |>
               mutate(Difference = `TRH-RSP` - `T-RSP`), by = "Name")
change_lab <- labels_at(change_sf)
p_change <- ggplot() +
  geom_sf(data = change_sf, aes(fill = Difference), colour = "grey35", linewidth = 0.2) +
  geom_text(data = change_lab, aes(lon, lat, label = round(Difference, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis_c(option = "viridis", direction = -1,
                       name = expression("Difference in " * f[q*","*snow] * " between T-RSP and TRH-RSP")) +
  lonlat_axes + guides(fill = guide_colorbar(barwidth = unit(6, "cm"), barheight = unit(0.8, "cm"),
                                             title.position = "bottom", title.hjust = 0.5)) +
  labs(title = "CHANGE", x = NULL, y = NULL) + theme_fq +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        plot.title = element_text(face = "bold", size = 18),
        legend.title = element_text(size = 11, face = "bold"))
fig4 <- (p_two | p_change) + plot_layout(widths = c(2, 1))
save_png(fig4, "Fig4_fqsnow_historical_difference_CRB.png", 40, 17)

## =============================================================================
## FIGURE 5 / A3  dumbbell: Historical -> 2040s -> 2080s per sub-basin
## =============================================================================
class_order  <- c("LowerCRB", "CentralCRB", "RockiesCRB", "CascadeCRB", "HeadwaterCRB")  # bottom -> top
class_colors <- c(HeadwaterCRB = "#E6AB02", CascadeCRB = "#66A61E", RockiesCRB = "#E7298A",
                  CentralCRB = "#7570B3", LowerCRB = "#D95F02")
gradient_df <- tibble(xmin = seq(0, 0.9, 0.1), xmax = seq(0.1, 1, 0.1),
                      fill = colorRampPalette(c("lightblue", "grey"))(10))

cls <- st_drop_geometry(subbasin) |> distinct(Name, Classification) |>
  mutate(Classification = factor(Classification, levels = class_order)) |>
  arrange(Classification, Name)
y_labels <- setNames(sprintf("<span style='color:%s'>%s</span>", class_colors[as.character(cls$Classification)], cls$Name), cls$Name)

dumbbell <- function(emission_code, emission_label) {
  d <- bind_rows(basin_hist, basin_fut) |>
    filter(RCP %in% c("Historical", emission_code), Name %in% cls$Name) |>
    mutate(scenario = factor(Timeframe, levels = c("Historical", "2040s", "2080s")),
           Name = factor(Name, levels = cls$Name),
           Partition = factor(Partition, levels = c("T-RSP", "TRH-RSP")))
  ggplot(d) +
    geom_rect(data = gradient_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill), alpha = 0.4) +
    geom_vline(xintercept = 0.5, colour = "darkred", linetype = "dashed", linewidth = 1) +
    geom_line(aes(basin_fqsnow, Name, group = Name), colour = "gray70", linewidth = 2) +
    geom_point(aes(basin_fqsnow, Name, colour = scenario), size = 8, shape = 18) +
    scale_colour_manual(values = c(Historical = "#0072B2", `2040s` = "#56B4E9", `2080s` = "#009E73")) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_y_discrete(labels = y_labels) +
    annotate("text", x = 0, y = Inf, label = "Rain-dominant", vjust = 2, hjust = 0, colour = "grey40", fontface = "italic", size = 4) +
    annotate("text", x = 1, y = Inf, label = "Snow-dominant", vjust = 2, hjust = 1, colour = "steelblue", fontface = "italic", size = 4) +
    labs(x = NULL, y = NULL,
         title = paste0("Change in fq,snow value from Historical to 2040s to 2080s for ", emission_label),
         subtitle = "Background gradient shows range of fq,snow from low (rain-dominant) to high (snow-dominant)") +
    theme_fq + theme(axis.text.y = element_markdown(size = 16), axis.text.x = element_text(face = "plain"),
                     legend.text = element_text(size = 16, face = "bold"), legend.title = element_blank(),
                     plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12)) +
    facet_grid(~Partition)
}

inset <- ggplot(subbasin |> mutate(Classification = factor(Classification, levels = rev(class_order)))) +
  geom_sf(aes(fill = Classification), colour = "red", linewidth = 0.3) +
  scale_fill_manual(values = class_colors, breaks = rev(class_order), name = NULL) +
  theme_void() + theme(legend.position = "left", legend.text = element_text(size = 11))

save_png(inset + dumbbell("RCP85", "RCP8.5") + plot_layout(widths = c(1, 2.6)),
         "Fig5_fqsnow_dumbbell_RCP85_CRB.png", 36, 22)
save_png(inset + dumbbell("RCP45", "RCP4.5") + plot_layout(widths = c(1, 2.6)),
         "FigA3_fqsnow_dumbbell_RCP45_CRB.png", 36, 22)

## =============================================================================
## FIGURE A4  projected decrease in fQ,snow (historical minus future), per sub-basin
##            rows: T-RSP / TRH-RSP;  columns: RCP4.5 2040s | 2080s, RCP8.5 2040s | 2080s
## =============================================================================
hist_ref <- basin_hist |> dplyr::select(Name, Partition, historical = basin_fqsnow)
change_a4 <- basin_fut |>
  inner_join(hist_ref, by = c("Name", "Partition")) |>
  mutate(change = historical - basin_fqsnow,
         RCP = recode(RCP, RCP45 = "RCP4.5", RCP85 = "RCP8.5"),
         Partition = factor(Partition, levels = c("T-RSP", "TRH-RSP")))
write.csv(change_a4, file.path(OUT_DIR, "basin_fqsnow_change_future.csv"), row.names = FALSE)

a4_sf  <- subbasin |> inner_join(change_a4, by = "Name")
a4_lab <- labels_at(a4_sf)
figA4 <- ggplot() +
  geom_sf(data = a4_sf, aes(fill = change), colour = "grey35", linewidth = 0.2) +
  geom_text(data = a4_lab, aes(lon, lat, label = round(change, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis_c(option = "viridis", direction = -1,
                       name = expression("Difference in " * f[q*","*snow] * " between historical and future")) +
  lonlat_axes + guides(fill = colourbar) + labs(x = NULL, y = NULL) +
  theme_fq + theme(axis.text.x = element_text(angle = 90, hjust = 0.5), strip.text = element_text(size = 14)) +
  facet_grid(Partition ~ RCP + Timeframe)
save_png(figA4, "FigA4_fqsnow_change_future_CRB.png", 45, 25)

## =============================================================================
## FIGURE A2  grids contributing 80% of the runoff in each sub-watershed (CRB only)
##   Only cells inside the CRB and in the kept sub-basins are ranked, so the 80% is
##   80% of the CRB-clipped runoff. Label = % of a sub-basin's grids needed to reach 80%.
## =============================================================================
ranked <- hist_grid |>
  inner_join(dplyr::select(huc6_crb, lat_lon, lat, lon, Name), by = c("lat_lon", "lat", "lon")) |>
  group_by(Partition, Name) |>
  arrange(desc(meanRunoff), .by_group = TRUE) |>
  mutate(percent = cumsum(meanRunoff) / sum(meanRunoff) * 100) |>
  ungroup()

A2_PARTITION <- "TRH-RSP"                       # set to NULL to facet both methods
cells80 <- ranked |> filter(percent <= 80)

share80 <- ranked |> group_by(Name, Partition) |>
  summarise(total_grids = n(), grids_80 = sum(percent <= 80),
            percent_grids = round(grids_80 / total_grids * 100, 2), .groups = "drop")
write.csv(share80, file.path(OUT_DIR, "FigA2_grids_80pct_runoff.csv"), row.names = FALSE)

a2_lab <- labels_at(inner_join(subbasin, share80, by = "Name"))   # label points + percent_grids
if (!is.null(A2_PARTITION)) {
  cells80 <- filter(cells80, Partition == A2_PARTITION)
  a2_lab  <- filter(a2_lab,  Partition == A2_PARTITION)
}
res_a2 <- median(diff(sort(unique(hist_grid$lon))))
figA2 <- ggplot() +
  geom_tile(data = cells80, aes(lon, lat), width = res_a2, height = res_a2, fill = "steelblue") +
  geom_sf(data = subbasin, fill = NA, colour = "gray40", linewidth = 0.3, inherit.aes = FALSE) +
  crb_outline +
  geom_text(data = a2_lab, aes(lon, lat, label = round(percent_grids, 2)),
            fontface = "bold", size = 4) +
  coord_sf(crs = 4326) + lonlat_axes + labs(x = NULL, y = NULL) + theme_fq
if (is.null(A2_PARTITION)) figA2 <- figA2 + facet_wrap(~Partition)
save_png(figA2, "FigA2_grids_80pct_runoff_CRB.png", 20, 20)

## =============================================================================
## CRB-LEVEL fQ,snow (runoff-weighted), both methods
##   = sum(fQsnow x runoff) / sum(runoff) over all cells inside the CRB and in the kept
##   sub-basins (the same cells as the basin values above). Historical and future.
## =============================================================================
crb_hist <- basin_value(hist_grid, cells_use) |>
  group_by(Partition) |>
  summarise(crb_fqsnow = sum(Snow_runoff) / sum(meanRunoff), n_cells = n(), .groups = "drop") |>
  mutate(Timeframe = "Historical", RCP = "Historical", Partition = as.character(Partition))

crb_fut <- future_grid |> filter(Timeframe != "2060s") |>
  group_by(lat_lon, lat, lon, Timeframe, Partition, RCP) |>          # GCM mean per cell
  summarise(fq = mean(meanfqsnow), ro = mean(meanRunoff), .groups = "drop") |>
  inner_join(dplyr::select(cells_use, lat_lon, lat, lon, Name), by = c("lat_lon", "lat", "lon")) |>
  group_by(Timeframe, Partition, RCP) |>
  summarise(crb_fqsnow = sum(fq * ro) / sum(ro), n_cells = n(), .groups = "drop")

crb_level <- bind_rows(crb_hist, crb_fut) |> mutate(crb_fqsnow = round(crb_fqsnow, 3))
write.csv(crb_level, file.path(OUT_DIR, "crb_level_fqsnow.csv"), row.names = FALSE)
cat("\nCRB-level runoff-weighted fQ,snow:\n"); print(crb_level, n = Inf)

message("Done. Figures and CSVs in: ", OUT_DIR)

}  # end RUN_FIGURES
