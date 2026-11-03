
# Climate change impacts on snowmelt contributions to streamflow and implications on ecosystem services.

## Authors

Analysis and R code: Supriya Savalkar.
Code cleanup and documentation assisted by Claude (Anthropic).

# Rain–snow partitioning and ecosystem services in the Columbia River Basin

Data-processing and plotting code for *Savalkar, Yourek, Brady & Rajagopalan* – how the choice of
rain–snow partitioning changes the snowmelt fraction of streamflow (fQ,snow), the water balance,
and flood, hydropower and minimum-flow ecosystem services under climate change.

| Name in the paper | Meaning | Code in the data files |
|---|---|---|
| **TRH-RSP** | partitioning from temperature **+ relative humidity** | `DRS` |
| **T-RSP** | partitioning from **temperature only** | `SRS` |

Differences are always **T-RSP minus TRH-RSP**. Every figure is clipped to the Columbia River Basin (CRB).

> This repository contains the **analysis after the hydrologic models**. It does **not** contain the
> VIC-CropSyst, routing, bias-correction or reservoir-model code, and it does **not** contain data

## Workflow – how to reproduce the results

Three scripts, run in this order. 

```
 STEP 0  (outside this repo)   VIC-CropSyst run TWICE, once per partitioning method (DRS and SRS),
         │                     historical + every future GCM x RCP  ->  daily flux file per grid cell
         │
         ├─ STEP 1  scripts/01_fqsnow.R
         │            Part A  fQ,snow per grid cell, both methods    
         │            Part B  Figs 3, 4, 5, A2, A3, A4 + CRB-level fQ,snow, clipped to the CRB
         │
         ├─ STEP 2  scripts/02_water_balance.R             (mass-balance check, both methods)
         │            Part A  collate the flux files to annual water balances 
         │            Part B  Table A1 and Figure A1, CRB only
         │            Part C  optional annual-resolution bias check (all cells)
         │
         └─ STEP 3  (outside this repo, TRH-RSP / DRS ONLY)
                      routing  ->  bias correction  ->  reservoir model
                      ->  daily dam flows + monthly firm energy
                        └─ STEP 4  scripts/03_ecosystem_services_TRHRSP.R
                                     hydropower (Figs 7, A5), flood at The Dalles (Fig 6),
                                     minimum-flow heatmaps (Fig 8, A6-A13) and all their tables
```

0. **Run VIC-CropSyst with both partitioning methods.** Same forcing, two partitioning schemes (TRH-RSP = `DRS`,
   T-RSP = `SRS`), for the historical period and for each future GCM × RCP. Keep the daily output per
   grid cell as tab-delimited files named `flux_<lat>_<lon>` (32 columns, no header; column order in
   `scripts/01_fqsnow.R`), one folder per method.
1. **fQ,snow** – `01_fqsnow.R` turns each cell's daily flux file into a water-year fQ,snow and averages it over
   the period (Part A), then aggregates the cells to HUC6 sub-basins (runoff-weighted), clips everything to
   the CRB and draws the figures (Part B).
2. **Mass-balance check** – `02_water_balance.R` collapses the raw flux files to annual water balances for both
   methods (Part A), computes Table A1 and Fig A1 for the CRB (Part B: differences in evapotranspiration,
   runoff and baseflow, with precipitation identical in both methods) and optionally repeats the check
   year by year (Part C, not CRB-clipped).
3. **Only TRH-RSP** continues to the water-resources side: its runoff is routed to the dam locations,
   bias-corrected, and run through the reservoir model. T-RSP is not carried through this chain, so the
   ecosystem-service results are TRH-RSP only.
4. **Ecosystem services** – `03_ecosystem_services_TRHRSP.R` reads the reservoir-model output and produces the
   hydropower, flood (The Dalles) and minimum-flow results.

## Setup

1. R ≥ 4.1 (native pipe `|>` and `\(x)` are used).
2. Packages: `dplyr tidyr purrr stringr ggplot2 sf terra patchwork scales ggtext broom` (+ `data.table`
   recommended for speed). Scripts stop with an install hint if one is missing. Consider `renv::init()`.
3. `cp config/config.example.R config/config.R`, then set `data_root`, `huc6_shp`, `crb_shp`, `out_root`.
4. Open `rsp-columbia.Rproj` (or `setwd()` to the repository folder). Run a script with
   `source("scripts/01_fqsnow.R")`, or everything with `source("run_all.R")`.
5. The calculation parts are slow (one file read per cell, twice per method). Run them once; afterwards set
   `RUN_CALCULATION <- FALSE` (script 01) / `RUN_COLLATE <- FALSE` (script 02) to redraw figures only.
   Before a first full run, test with `MAX_CELLS <- 50`.
6. Script 01 future runs are off by default (`RUN_FUTURE <- FALSE`); fill in `FUTURE_RUNS`,
   `future_flux_dir()` and `FUTURE_PERIODS` first.

## Scripts

| Script | Paper | Reads | Writes |
|---|---|---|---|
| `01_fqsnow.R` | 3.1, 3.2; Figs 3, 4, 5, A2, A3, A4 | raw flux files (A); HUC6 + CRB shapefiles (B) | `gridlevel_hist_{DRS,SRS}.rds`, `<GCM>_<RCP>_<DRS\|SRS>.rds`; figures, `crb_level_fqsnow.csv`, basin tables |
| `02_water_balance.R` | §2.9, Table A1, Fig A1 | raw flux files (A); CRB shapefile (B) | `Waterbalance_{DRS,SRS}*.rds`; Table A1 CSVs, Fig A1, diagnostics; optional annual-check files |
| `03_ecosystem_services_TRHRSP.R` | 3.3.1–3.3.3; Figs 6, 7, A5, 8, A6–A13 | reservoir-model master files, Table A2 | figures + CSVs |
