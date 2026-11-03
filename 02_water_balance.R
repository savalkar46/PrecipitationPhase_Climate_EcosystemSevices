# =============================================================================
# 02_water_balance.R   --   mass-balance comparison, T-RSP vs TRH-RSP, in one script
## Author:      Supriya Savalkar
## Created:     2025-10-09
## Last updated: 2026-10-05
## Code cleanup and documentation assisted by Claude (Anthropic)
# PAPER SECTIONS : 2.9 and Appendix A (Table A1, Figure A1)
#   DRS = TRH-RSP (temperature + relative humidity)    SRS = T-RSP (temperature only)
#   Differences are T-RSP minus TRH-RSP.
#
# PART A  COLLATE  (slow: one file read per grid cell, twice; RUN_COLLATE <- FALSE to skip)
#   Reads the raw VIC-CropSyst daily flux file of every grid cell, for each method, and
#   collapses it to ANNUAL (water-year) totals plus snowpack metrics:
#     SumPrecip, SumEvap, SumRunoff, SumBase (mm / water year); Evap_perc, Runoff_perc,
#     Base_perc (% of precipitation); PeakSWE, PeakSWE_Day, MeltOutDay (day of water year).
#   Then writes 30-year-mean summaries per cell.
#
# PART B  TABLE A1 AND FIGURE A1, CRB ONLY  (RUN_TABLE_FIGURE)
#   Keeps the cells whose centre lies inside the CRB polygon (fast raster filter) and
#   recomputes the domain-mean fluxes (Table A1, also all-cells for comparison) and the
#   four-panel Figure A1 (a ET difference, b streamflow difference, c divergence by
#   precipitation decile, d mean fluxes).
#
# PART C  ANNUAL-RESOLUTION CHECK  (RUN_INTERANNUAL_CHECK; all cells, NOT clipped to the CRB)
#   Bias check (mean of annual ratios vs sum of sums, the ~0.13 pp in the text), year-by-
#   year series and trends, closure residual by year, zero-baseflow diagnostics, 4-panel figure.
#
# INPUTS  : paths$vic_flux_dir("DRS"/"SRS") (flux_<lat>_<lon>, tab-delimited, 32 columns); CRB shapefile
# OUTPUTS : paths$wb_annual/wb_summary/wb_summary_swe(part)  Waterbalance_<DRS|SRS>*.rds   (Part A)
#           outputs/water_balance/: TableA1_flux_table_CRB.csv, TableA1_all_vs_CRB.csv,
#           FigA1_rsp_waterbalance_CRB.png, precip_gradient_CRB.csv, diagnostics_CRB.txt (Part B);
#           bias_check.csv, flux_table_direct.csv, yearly_series.csv, closure_by_year.csv,
#           diagnostics_interannual.txt, rsp_interannual.png (Part C)
#
# Changes from the original collation loops: both partitions and both outputs come from one
# pass; outputs are NOT written into the VIC flux folder; MAX_CELLS defaults to ALL cells
# (the original loop was `for (i in 1:100)`).
# =============================================================================
source("R/setup.R")
use_packages(c("dplyr", "tidyr", "stringr", "ggplot2", "patchwork", "sf", "terra"))
sf_use_s2(FALSE)

## ---- 0. SETTINGS -----------------------------------------------------------------
RUN_COLLATE            <- TRUE     # Part A (slow). FALSE = reuse the Waterbalance_*.rds files
RUN_TABLE_FIGURE       <- TRUE     # Part B
RUN_INTERANNUAL_CHECK  <- TRUE     # Part C
CRB_SHP <- cfg$crb_shp
WB_DRS  <- paths$wb_summary("DRS")      # Waterbalance_DRS.rds  (written by Part A)
WB_SRS  <- paths$wb_summary("SRS")      # Waterbalance_SRS.rds  (written by Part A)
OUT_DIR <- out_dir("water_balance")

## ---- PART A. COLLATE -- settings and helpers ---------------------------------------
PARTITIONS    <- c("DRS", "SRS")
WY_START      <- as.Date("1985-10-01")   # water years 1986-2015
WY_END        <- as.Date("2015-09-30")
MAX_CELLS     <- Inf                     # set e.g. 100 for a quick test run
MELT_SWE_MM   <- 1                       # pack counts as melted out at SWE <= 1 mm

VIC_FLUX_COLS <- c("Year", "Month", "Day", "Precip", "Evap", "Runoff", "Baseflow", "Rain",
                   "Snow", "Snow_Depth", "SWE", "Delta_SWE", "SnowMelt", "SnowEvap", "AirTemp",
                   paste0("SM", 1:17))

## Water year starting 1 October (months Oct-Dec belong to the next water year)
water_year <- function(dates, start_month = 10) {
  as.integer(format(dates, "%Y")) + as.integer(as.integer(format(dates, "%m")) >= start_month)
}

## Peak SWE, its day, and the first day after the peak with SWE <= MELT_SWE_MM.
## NOTE: a cell-year with no snow at all returns PeakSWE = 0, PeakSWE_Day = 1 and
## MeltOutDay = 2. A cell-year that never melts out returns NA.
snow_metrics <- function(swe, dowy) {
  i     <- which.max(swe)
  after <- which(dowy > dowy[i] & swe <= MELT_SWE_MM)
  list(PeakSWE = swe[i], PeakSWE_Day = dowy[i],
       MeltOutDay = if (length(after)) dowy[after[1]] else NA_integer_)
}

process_cell <- function(file) {
  tmp <- tryCatch(read.delim(file, header = FALSE), error = function(e) NULL)
  if (is.null(tmp) || ncol(tmp) != length(VIC_FLUX_COLS)) return(NULL)
  names(tmp) <- VIC_FLUX_COLS

  ll <- str_match(basename(file), "flux_([0-9.]+)_(-[0-9.]+)")
  latitude  <- as.numeric(ll[2])
  longitude <- as.numeric(ll[3])

  tmp |>
    mutate(Date      = as.Date(paste(Year, Month, Day, sep = "-"), format = "%Y-%m-%d"),
           WaterYear = water_year(Date)) |>
    filter(Date >= WY_START, Date <= WY_END) |>
    mutate(DOWY = as.integer(Date - as.Date(paste0(WaterYear - 1, "-10-01"))) + 1L) |>
    arrange(Date) |>
    group_by(WaterYear) |>
    group_modify(function(d, key) {
      sm <- snow_metrics(d$SWE, d$DOWY)
      tibble::tibble(
        SumPrecip   = sum(d$Precip),
        SumEvap     = sum(d$Evap),
        SumRunoff   = sum(d$Runoff),
        SumBase     = sum(d$Baseflow),
        Runoff_Base = sum(d$Runoff) + sum(d$Baseflow),
        WaterBal_Q  = sum(d$Precip) - sum(d$Evap),
        Evap_perc   = round(sum(d$Evap)      / sum(d$Precip) * 100, 2),
        Runoff_perc = round(sum(d$Runoff)    / sum(d$Precip) * 100, 2),
        Base_perc   = round(sum(d$Baseflow)  / sum(d$Precip) * 100, 2),
        PeakSWE     = sm$PeakSWE,
        PeakSWE_Day = sm$PeakSWE_Day,
        MeltOutDay  = sm$MeltOutDay)
    }) |>
    ungroup() |>
    mutate(lat = latitude, lon = longitude)
}

if (RUN_COLLATE) {

for (part in PARTITIONS) {
  files <- list.files(paths$vic_flux_dir(part), full.names = TRUE)
  if (is.finite(MAX_CELLS)) files <- files[seq_len(min(length(files), MAX_CELLS))]
  message(sprintf("[%s] %d grid cells", part, length(files)))

  t0 <- Sys.time()
  annual <- bind_rows(lapply(seq_along(files), function(i) {
    if (i %% 250 == 0) message(sprintf("  [%s] %d / %d", part, i, length(files)))
    process_cell(files[i])
  }))
  message("  elapsed: ", format(Sys.time() - t0))

  saveRDS(annual, paths$wb_annual(part))

  ## 30-year-mean summary per cell; column suffix = partition code (Pmean_DRS, ...)
  sfx <- paste0("_", part)
  flux_summary <- annual |>
    group_by(lat, lon) |>
    summarise(Pmean = mean(SumPrecip), Eperc = mean(Evap_perc),
              Rperc = mean(Runoff_perc), Bperc = mean(Base_perc), .groups = "drop")
  saveRDS(rename_with(flux_summary, ~ paste0(.x, sfx), c(Pmean, Eperc, Rperc, Bperc)),
          paths$wb_summary(part))

  snow_summary <- annual |>
    group_by(lat, lon) |>
    summarise(PeakSWE = mean(PeakSWE), MeltOutDay = mean(MeltOutDay, na.rm = TRUE),
              .groups = "drop")
  saveRDS(left_join(flux_summary, snow_summary, by = c("lat", "lon")) |>
            rename_with(~ paste0(.x, sfx), c(Pmean, Eperc, Rperc, Bperc, PeakSWE, MeltOutDay)),
          paths$wb_summary_swe(part))
}

}  # end RUN_COLLATE

## =============================================================================
## PART B. TABLE A1 AND FIGURE A1, CRB ONLY
## =============================================================================
if (RUN_TABLE_FIGURE) {

## ---- 1. Load, join, and flag the cells inside the CRB -------------------------------
drs <- readRDS(WB_DRS) |> ungroup()
srs <- readRDS(WB_SRS) |> ungroup()
m_all <- inner_join(drs, srs, by = c("lat", "lon")) |> na.omit()

crb <- st_sf(geometry = st_union(st_make_valid(st_read(CRB_SHP, quiet = TRUE))))
xres <- median(diff(sort(unique(m_all$lon)))); yres <- median(diff(sort(unique(m_all$lat))))
r <- terra::rast(xmin = min(m_all$lon) - xres / 2, xmax = max(m_all$lon) + xres / 2,
                 ymin = min(m_all$lat) - yres / 2, ymax = max(m_all$lat) + yres / 2,
                 resolution = c(xres, yres), crs = "EPSG:4326")
mask <- terra::rasterize(terra::vect(crb), r, field = 1, touches = FALSE)      # centre rule
m_all$inside <- !is.na(terra::extract(mask, as.matrix(m_all[, c("lon", "lat")]))[[1]])
cat("cells in the file:", nrow(m_all), "  inside the CRB:", sum(m_all$inside), "\n")

## ---- 2. Differences (T-RSP minus TRH-RSP) and fluxes in mm/yr ------------------------
add_fluxes <- function(d) d |>
  mutate(Ediff = Eperc_SRS - Eperc_DRS, Rdiff = Rperc_SRS - Rperc_DRS, Bdiff = Bperc_SRS - Bperc_DRS,
         Qdiff = Rdiff + Bdiff,
         E_D_mm = Pmean_DRS * Eperc_DRS / 100, R_D_mm = Pmean_DRS * Rperc_DRS / 100, B_D_mm = Pmean_DRS * Bperc_DRS / 100,
         E_S_mm = Pmean_SRS * Eperc_SRS / 100, R_S_mm = Pmean_SRS * Rperc_SRS / 100, B_S_mm = Pmean_SRS * Bperc_SRS / 100,
         Q_D_mm = R_D_mm + B_D_mm, Q_S_mm = R_S_mm + B_S_mm,
         close_D = Eperc_DRS + Rperc_DRS + Bperc_DRS, close_S = Eperc_SRS + Rperc_SRS + Bperc_SRS)
m_all <- add_fluxes(m_all)
m <- filter(m_all, inside)                      # CRB only

## ---- 3. Table A1 ---------------------------------------------------------------------
flux_table <- function(d) tibble(
  Flux   = c("Precipitation", "Evapotranspiration", "Surface runoff", "Baseflow", "Total streamflow"),
  TRH_RSP = c(mean(d$Pmean_DRS), mean(d$E_D_mm), mean(d$R_D_mm), mean(d$B_D_mm), mean(d$Q_D_mm)),
  T_RSP   = c(mean(d$Pmean_SRS), mean(d$E_S_mm), mean(d$R_S_mm), mean(d$B_S_mm), mean(d$Q_S_mm))) |>
  mutate(Diff = T_RSP - TRH_RSP, Pct = 100 * Diff / TRH_RSP) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))

tab_crb <- flux_table(m); tab_all <- flux_table(m_all)
write.csv(tab_crb, file.path(OUT_DIR, "TableA1_flux_table_CRB.csv"), row.names = FALSE)
write.csv(bind_rows(mutate(tab_all, Domain = "All cells"), mutate(tab_crb, Domain = "CRB only")),
          file.path(OUT_DIR, "TableA1_all_vs_CRB.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "diagnostics_CRB.txt"), split = TRUE)
cat("Table A1, CRB only (mm/yr):\n"); print(tab_crb)
cat("\nTable A1, all cells (for comparison):\n"); print(tab_all)
cat("\ncells CRB:", nrow(m), "  all:", nrow(m_all), "\n")
cat("max |P_DRS - P_SRS| =", max(abs(m$Pmean_DRS - m$Pmean_SRS)), "mm\n")
cat("frac cells T-RSP has more ET:", round(mean(m$Ediff > 0), 4), "\n")
cat("frac cells T-RSP has less Q :", round(mean(m$Qdiff < 0), 4), "\n")
cat("baseflow share of total Q loss:",
    round(100 * (mean(m$B_D_mm) - mean(m$B_S_mm)) / (mean(m$Q_D_mm) - mean(m$Q_S_mm)), 1), "%\n")
cat("mean closure TRH-RSP / T-RSP (%):", round(mean(m$close_D), 2), "/", round(mean(m$close_S), 2), "\n")

## divergence vs wetness (deciles of mean annual precipitation)
pgrad <- m |> mutate(Pbin = ntile(Pmean_DRS, 10)) |> group_by(Pbin) |>
  summarise(P = median(Pmean_DRS), ET = mean(Ediff), Runoff = mean(Rdiff), Baseflow = mean(Bdiff), .groups = "drop")
print(pgrad |> mutate(across(where(is.numeric), \(x) round(x, 2))))
write.csv(pgrad, file.path(OUT_DIR, "precip_gradient_CRB.csv"), row.names = FALSE)
sink()

## ---- 4. Figure A1 --------------------------------------------------------------------
base_theme <- theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 11),
        axis.text = element_text(face = "bold", colour = "black"),
        plot.title = element_text(face = "bold", size = 11), legend.position = "right")
outline <- geom_sf(data = crb, fill = NA, colour = "grey25", linewidth = 0.25, inherit.aes = FALSE)

p1 <- ggplot(m, aes(lon, lat, fill = Ediff)) +
  geom_tile(width = xres, height = yres) + outline +
  scale_fill_viridis_c(option = "C", direction = -1, limits = c(0, 20), oob = scales::squish) +
  coord_sf(crs = 4326) + base_theme +
  labs(x = NULL, y = NULL, title = "ET difference (T-RSP \u2013 TRH-RSP), % of P", fill = NULL)

p2 <- ggplot(m, aes(lon, lat, fill = Qdiff)) +
  geom_tile(width = xres, height = yres) + outline +
  scale_fill_viridis_c(option = "D", limits = c(-18, 0), oob = scales::squish) +
  coord_sf(crs = 4326) + base_theme +
  labs(x = NULL, y = NULL, title = "Total streamflow difference (T-RSP \u2013 TRH-RSP), % of P", fill = NULL)

p3 <- pgrad |> pivot_longer(c(ET, Runoff, Baseflow)) |>
  ggplot(aes(P, value, colour = name)) +
  geom_hline(yintercept = 0, linewidth = 0.4) + geom_line() + geom_point() +
  scale_x_log10() +
  scale_colour_manual(values = c(ET = "#CC6677", Runoff = "#332288", Baseflow = "#44AA99")) +
  base_theme + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1)) +
  labs(x = "Mean annual precipitation (mm)", y = "T-RSP \u2212 TRH-RSP (% of P)",
       title = "Divergence grows with wetness", colour = NULL)

bars <- tab_crb |> filter(Flux %in% c("Evapotranspiration", "Surface runoff", "Baseflow")) |>
  mutate(Flux = recode(Flux, Evapotranspiration = "ET", `Surface runoff` = "Runoff")) |>
  mutate(Flux = factor(Flux, levels = c("ET", "Runoff", "Baseflow")))
p4 <- bars |> pivot_longer(c(TRH_RSP, T_RSP), names_to = "method", values_to = "value") |>
  mutate(method = factor(recode(method, TRH_RSP = "TRH-RSP", T_RSP = "T-RSP"), levels = c("TRH-RSP", "T-RSP"))) |>
  ggplot(aes(Flux, value, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(data = bars, inherit.aes = FALSE,
            aes(Flux, pmax(TRH_RSP, T_RSP) + 12, label = sprintf("%+.0f mm", Diff)),
            fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c(`TRH-RSP` = "#332288", `T-RSP` = "#888888")) +
  base_theme + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1)) +
  labs(x = NULL, y = "mm/yr", fill = NULL,
       title = sprintf("Domain mean fluxes (P = %.0f mm/yr)", tab_crb$TRH_RSP[tab_crb$Flux == "Precipitation"]))

figA1 <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(file.path(OUT_DIR, "FigA1_rsp_waterbalance_CRB.png"), figA1, width = 13, height = 10, dpi = 300)
message("Table A1 and Figure A1 written to: ", OUT_DIR)

}  # end RUN_TABLE_FIGURE

## =============================================================================
## PART C. ANNUAL-RESOLUTION CHECK (all cells, not clipped to the CRB)
## =============================================================================
if (RUN_INTERANNUAL_CHECK) {

out <- OUT_DIR
sink(file.path(out, "diagnostics_interannual.txt"), split = TRUE)

## ---- 1. Load and join on WaterYear + lat + lon --------------------------------
## Only the water-balance columns are used, so cells that never melt out
## (MeltOutDay = NA) are not dropped by na.omit() below.
wb_cols <- c("WaterYear", "lat", "lon", "SumPrecip", "SumEvap", "SumRunoff", "SumBase",
             "Runoff_Base", "WaterBal_Q", "Evap_perc", "Runoff_perc", "Base_perc")
drs <- readRDS(paths$wb_annual("DRS")) |> dplyr::select(all_of(wb_cols))
srs <- readRDS(paths$wb_annual("SRS")) |> dplyr::select(all_of(wb_cols))

m <- left_join(drs, srs, by = c("WaterYear", "lat", "lon"), suffix = c("_D", "_S")) |> na.omit()

cat("cells:", m |> distinct(lat, lon) |> nrow(), " years:", n_distinct(m$WaterYear),
    " (", min(m$WaterYear), "-", max(m$WaterYear), ")\n")
cat("lat range:", range(m$lat), "  lon range:", range(m$lon), "\n")
cat("P identical?", isTRUE(all.equal(m$SumPrecip_D, m$SumPrecip_S)), "\n")

## ---- 2. Bias check: mean-of-annual-ratios vs sum-of-sums ----------------------
cell_old <- m |>
  group_by(lat, lon) |>
  summarise(P = mean(SumPrecip_D),
            E_old_D = mean(Evap_perc_D),   E_old_S = mean(Evap_perc_S),
            R_old_D = mean(Runoff_perc_D), R_old_S = mean(Runoff_perc_S),
            B_old_D = mean(Base_perc_D),   B_old_S = mean(Base_perc_S),
            sP   = sum(SumPrecip_D),
            sE_D = sum(SumEvap_D),   sE_S = sum(SumEvap_S),
            sR_D = sum(SumRunoff_D), sR_S = sum(SumRunoff_S),
            sB_D = sum(SumBase_D),   sB_S = sum(SumBase_S),
            .groups = "drop") |>
  mutate(E_new_D = 100 * sE_D / sP, E_new_S = 100 * sE_S / sP,
         R_new_D = 100 * sR_D / sP, R_new_S = 100 * sR_S / sP,
         B_new_D = 100 * sB_D / sP, B_new_S = 100 * sB_S / sP,
         Qdiff_old = (R_old_S + B_old_S) - (R_old_D + B_old_D),
         Qdiff_new = (R_new_S + B_new_S) - (R_new_D + B_new_D),
         Ediff_old = E_old_S - E_old_D,
         Ediff_new = E_new_S - E_new_D)

bias_check <- tibble(
  quantity = c("ET difference (pp of P)", "Total Q difference (pp of P)"),
  mean_of_annual_ratios = c(mean(cell_old$Ediff_old), mean(cell_old$Qdiff_old)),
  sum_of_sums           = c(mean(cell_old$Ediff_new), mean(cell_old$Qdiff_new))) |>
  mutate(bias = mean_of_annual_ratios - sum_of_sums, across(where(is.numeric), \(x) round(x, 3)))
cat("\n--- bias check (mean-of-ratios vs sum-of-sums) ---\n")
print(bias_check)
write.csv(bias_check, file.path(out, "bias_check.csv"), row.names = FALSE)

## ---- 3. Domain-mean absolute fluxes (mm/yr), computed directly -----------------
m <- m |>
  mutate(Qmm_D = SumRunoff_D + SumBase_D,
         Qmm_S = SumRunoff_S + SumBase_S,
         Qdiff_mm = Qmm_S - Qmm_D,
         Ediff = Evap_perc_S - Evap_perc_D,
         close_D = Evap_perc_D + Runoff_perc_D + Base_perc_D,
         close_S = Evap_perc_S + Runoff_perc_S + Base_perc_S)

flux_tab <- m |>
  summarise(P = mean(SumPrecip_D),
            E_D = mean(SumEvap_D),   E_S = mean(SumEvap_S),
            R_D = mean(SumRunoff_D), R_S = mean(SumRunoff_S),
            B_D = mean(SumBase_D),   B_S = mean(SumBase_S),
            Q_D = mean(Qmm_D),       Q_S = mean(Qmm_S)) |>
  mutate(Q_diff_pct = 100 * (Q_S - Q_D) / Q_D) |>
  mutate(across(everything(), \(x) round(x, 1)))
print(flux_tab)
write.csv(flux_tab, file.path(out, "flux_table_direct.csv"), row.names = FALSE)

## ---- 4. Cell-year consistency -------------------------------------------------
cat("\n--- cell-year consistency ---\n")
cat("frac cell-years SRS has more ET:", round(mean(m$Ediff > 0), 3), "\n")
cat("frac cell-years SRS has less Q :", round(mean(m$Qdiff_mm < 0), 3), "\n")

## ---- 5. Year-by-year series ---------------------------------------------------
yr <- m |>
  group_by(WaterYear) |>
  summarise(P = mean(SumPrecip_D), Q_D = mean(Qmm_D), Q_S = mean(Qmm_S), .groups = "drop") |>
  mutate(dQ = Q_S - Q_D, dQpct = 100 * dQ / Q_D, Panom = 100 * (P - mean(P)) / mean(P))
write.csv(yr, file.path(out, "yearly_series.csv"), row.names = FALSE)

cat("\n--- every year, DRS wetter? ---\n")
cat("years where SRS >= DRS:", sum(yr$dQ >= 0), "of", nrow(yr), "\n")
cat("range of relative loss:", round(range(yr$dQpct), 1), "%\n")
cat("corr(P, dQ mm):", round(cor(yr$P, yr$dQ), 3), "\n")
cat("corr(P, dQ %) :", round(cor(yr$P, yr$dQpct), 3), "\n")

wet <- slice_max(yr, P, n = 8)
dry <- slice_min(yr, P, n = 8)
cat("wettest 8 yrs: mean dQ% =", round(mean(wet$dQpct), 1), "\n")
cat("driest  8 yrs: mean dQ% =", round(mean(dry$dQpct), 1), "\n")

## simple linear trend (swap in trend::mk.test for a Mann-Kendall test if preferred)
trend_test <- function(x, y) summary(lm(y ~ x))$coefficients[2, c(1, 4)]
cat("\n--- trends (slope, p-value) ---\n")
print(rbind(Q_DRS = trend_test(yr$WaterYear, yr$Q_D),
            Q_SRS = trend_test(yr$WaterYear, yr$Q_S),
            dQpct = trend_test(yr$WaterYear, yr$dQpct)))

## ---- 6. Closure residual by year (= -dS/dt, not a bug) ------------------------
close_yr <- m |>
  group_by(WaterYear) |>
  summarise(clD = mean(close_D), clS = mean(close_S), .groups = "drop")
cat("\n--- closure residual varies by year ---\n")
cat("sd of annual closure: DRS =", round(sd(close_yr$clD), 2), " SRS =", round(sd(close_yr$clS), 2), "\n")
print(close_yr |> mutate(across(where(is.numeric), \(x) round(x, 1))))
write.csv(close_yr, file.path(out, "closure_by_year.csv"), row.names = FALSE)

## ---- 7. Zero-flow / catchment-memory diagnostics ------------------------------
cat("\n--- dry-year / memory behaviour ---\n")
cat("cell-years with zero baseflow: DRS =", sum(m$SumBase_D <= 0), " SRS =", sum(m$SumBase_S <= 0), "\n")

n_years <- n_distinct(m$WaterYear)
cell_zero <- m |>
  group_by(lat, lon) |>
  summarise(zD = sum(SumBase_D <= 0), zS = sum(SumBase_S <= 0), .groups = "drop")
cat("cells that NEVER have baseflow (of", n_years, "yrs): DRS =", sum(cell_zero$zD == n_years),
    " SRS =", sum(cell_zero$zS == n_years), "\n")

## lag-1 autocorrelation of annual Q and corr with previous year's P, per cell then averaged
memory_stats <- m |>
  arrange(lat, lon, WaterYear) |>
  group_by(lat, lon) |>
  summarise(
    ac1_D  = if (sd(Qmm_D) > 0) cor(Qmm_D[-1], Qmm_D[-n()]) else NA_real_,
    ac1_S  = if (sd(Qmm_S) > 0) cor(Qmm_S[-1], Qmm_S[-n()]) else NA_real_,
    lagP_D = if (sd(Qmm_D) > 0) cor(Qmm_D[-1], SumPrecip_D[-n()]) else NA_real_,
    lagP_S = if (sd(Qmm_S) > 0) cor(Qmm_S[-1], SumPrecip_D[-n()]) else NA_real_,
    .groups = "drop")
cat("mean lag-1 autocorr of Q : DRS =", round(mean(memory_stats$ac1_D, na.rm = TRUE), 3),
    " SRS =", round(mean(memory_stats$ac1_S, na.rm = TRUE), 3), "\n")
cat("mean corr(Q_t, P_t-1)    : DRS =", round(mean(memory_stats$lagP_D, na.rm = TRUE), 3),
    " SRS =", round(mean(memory_stats$lagP_S, na.rm = TRUE), 3), "\n")

sink()

## ---- 8. Figure ----------------------------------------------------------------
base_theme <- theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 11))

p1 <- ggplot(yr, aes(WaterYear)) +
  geom_ribbon(aes(ymin = Q_S, ymax = Q_D), fill = "#332288", alpha = 0.15) +
  geom_line(aes(y = Q_D, colour = "DRS (TRH-RSP)")) + geom_point(aes(y = Q_D, colour = "DRS (TRH-RSP)")) +
  geom_line(aes(y = Q_S, colour = "SRS (T-RSP)"))   + geom_point(aes(y = Q_S, colour = "SRS (T-RSP)")) +
  scale_colour_manual(values = c("DRS (TRH-RSP)" = "#332288", "SRS (T-RSP)" = "#888888")) +
  base_theme + labs(y = "Annual streamflow (mm)", x = NULL, colour = NULL,
                    title = "Every year, DRS is wetter")

p2 <- ggplot(yr, aes(P, dQ, colour = WaterYear)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", linewidth = 0.6) +
  scale_colour_viridis_c() + base_theme +
  labs(x = "Annual precipitation (mm)", y = "Q difference SRS - DRS (mm)",
       title = "Absolute gap widens in wet years", colour = "Water\nyear")

p3 <- ggplot(yr, aes(WaterYear, dQpct, fill = dQpct < -15)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#CC6677", `FALSE` = "#332288"), guide = "none") +
  geom_hline(yintercept = mean(yr$dQpct), linetype = "dashed") +
  base_theme + labs(y = "Q difference (%)", x = NULL, title = "Relative loss is persistent, every year")

p4 <- ggplot(close_yr, aes(WaterYear, clD - 100, fill = clD > 100)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#CC6677", `FALSE` = "#44AA99"), guide = "none") +
  geom_hline(yintercept = 0) +
  base_theme + labs(y = "Closure residual (% of P)", x = NULL,
                    title = "Residual = storage change, oscillates about zero")

save_fig((p1 | p2) / (p3 | p4), file.path(out, "rsp_interannual.png"), 13, 9, dpi = 300)

}  # end RUN_INTERANNUAL_CHECK
message("Done.")
