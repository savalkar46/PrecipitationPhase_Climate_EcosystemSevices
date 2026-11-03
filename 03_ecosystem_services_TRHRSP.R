## =============================================================================
## Ecosystem-service figures and data  --  TRH-RSP only   (ONE self-contained script)
## Author:      Supriya Savalkar
## Created:     2025-10-09
## Last updated: 2026-10-05
## Code cleanup and documentation assisted by Claude (Anthropic)
##
## "DRS" in the data files = TRH-RSP (temperature + relative humidity)
##
## Produces (all in OUT_DIR):
##   1. HYDROPOWER   TRHRSP_hydropower_percent_deficit.png        (Fig. 7)
##                   TRHRSP_hydropower_achieved_firm_energy.png   (Fig. A5)
##                   TRHRSP_hydropower_monthly_summary.csv
##   2. FLOOD        TRHRSP_flooding_TheDalles.png                (Fig. 6)
##                   TRHRSP_flood_stats_TheDalles.csv
##                   TRHRSP_flood_yearly_exceedance_TheDalles.csv
##   3. MIN FLOWS    TRHRSP_minflow_<DAM>.png   (nine dams, monthly heatmaps)
##                   TRHRSP_minflow_monthly_targets.csv   (target per dam x month)
##                   TRHRSP_minflow_yearly.csv            (flow, target, flag per year)
##                   TRHRSP_minflow_monthly_metrics.csv   (frequency + shortfall; Table 1)
##
## Minimum flows are evaluated MONTHLY: the month's mean flow is compared with ONE
## monthly target. For months with several windows in Table A2 (Apr, Jun, Jul, Oct)
## that target is the average of the window targets (see TARGET_AVERAGE).
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

## ---- 0. SETTINGS (paths come from config/config.R via R/setup.R) --------------
source("R/setup.R")
OUT_DIR  <- out_dir("ecosystem_services_TRHRSP")
F_HYDRO  <- paths$master_mop       # firm energy, monthly (output of the reservoir model, TRH-RSP)
F_DAMS   <- paths$master_dams      # daily dam flows       (output of the reservoir model, TRH-RSP)
F_MINFLW <- paths$min_flow_rules   # Table A2 (data/Minimum_Flows_Mainstem_CRB.csv)

PARTITION        <- "DRS"                 # DRS = TRH-RSP
EXCLUDED_GCMS    <- c("CanESM2", "IPSL-CM5A-LR", "CNRM-CM5")
TARGET_AVERAGE   <- "day_weighted"        # or "simple" (plain mean of the window targets)
USE_ACTUAL_YEARS <- TRUE                  # frequency = failing years / years present (30);
                                          # FALSE divides by 31 as in the old slide code
ACTION_CFS <- 4.5e5                       # The Dalles action stage
DAMAGE_CFS <- 6.0e5                       # The Dalles damage stage

out <- function(f) file.path(OUT_DIR, f)

## ---- 1. Shared helpers ---------------------------------------------------------
tf_levels <- c("Historical", "2040s", "2080s")
tf_cols   <- c("Historical" = "grey35", "2040s" = "#4A90C8", "2080s" = "#08306B")

rcp_label <- function(s) case_when(s == "Future45" ~ "RCP4.5", s == "Future85" ~ "RCP8.5", TRUE ~ s)

## Historical is drawn in BOTH RCP panels so each reads Historical | 2040s | 2080s
duplicate_historical <- function(df) {
  h <- filter(df, Timeframe == "Historical")
  bind_rows(filter(df, Timeframe != "Historical"),
            mutate(h, RCP = "RCP4.5"), mutate(h, RCP = "RCP8.5")) |>
    mutate(RCP = factor(RCP, levels = c("RCP4.5", "RCP8.5")),
           Timeframe = factor(Timeframe, levels = tf_levels))
}

## Timeframes. closed_start = TRUE keeps the first year (31-yr windows, hydropower);
## FALSE drops it (30-yr windows, daily flows), exactly as in the original code.
add_timeframe <- function(df, closed_start) {
  ge <- if (closed_start) `>=` else `>`
  bind_rows(
    df |> filter(Scenario == "Historical", ge(year, 1975), year <= 2005) |> mutate(Timeframe = "Historical"),
    df |> filter(ge(year, 2025), year <= 2055) |> mutate(Timeframe = "2040s"),
    df |> filter(ge(year, 2064), year <= 2094) |> mutate(Timeframe = "2080s"))
}

theme_paper <- theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),
        axis.text = element_text(colour = "black", face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 18),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 14),
        legend.key.size = unit(1, "cm"), legend.box = "horizontal")

save_png <- function(p, file, w, h) { ggsave(out(file), p, width = w, height = h, dpi = 300); message("Saved ", file) }

## =============================================================================
## 2. HYDROPOWER  (Section 3.3.2; Figs. 7 and A5)
## =============================================================================
firm_ref <- tibble::tribble(
  ~Month,       ~Firm_MWh_week, ~Firm_thousand_GWh_week,
  "January",    11382919.76,    11.38291976,
  "February",    8471194.688,    8.471194688,
  "March",       7447029.632,    7.447029632,
  "April",       7228699.872,    7.228699872,
  "May",        10499676.64,    10.49967664,
  "June",        8312409.408,    8.312409408,
  "July",        8217138.24,     8.21713824,
  "August",      9533071.248,    9.533071248,
  "September",   6962734.528,    6.962734528,
  "October",     8137745.6,      8.1377456,
  "November",    7851932.096,    7.851932096,
  "December",    9146032.128,    9.146032128) |>
  mutate(month = match(Month, month.name))

hydro <- readRDS(F_HYDRO) |>
  separate(Model_Folder, into = c("Model", "Scenario"), sep = "_", remove = FALSE) |>
  filter(!Model %in% EXCLUDED_GCMS, Category == PARTITION) |>
  add_timeframe(closed_start = TRUE) |>
  group_by(month, year, Model, Scenario, Timeframe) |>
  summarise(Firmdeficit_total = sum(FirmEnergy), .groups = "drop") |>     # sum over weeks
  group_by(month, year, Scenario, Timeframe) |>
  summarise(Firmdeficit_mean_models = mean(Firmdeficit_total), .groups = "drop") |>  # mean over GCMs
  left_join(firm_ref, by = "month") |>
  mutate(Percent_deficit = round(Firmdeficit_mean_models / Firm_MWh_week * 100, 2),
         RCP = rcp_label(Scenario)) |>
  duplicate_historical() |>
  mutate(month_f = factor(month, levels = 1:12, labels = month.abb))

## Fig. 7: bars = mean over years, whiskers = +/-1 SD over years
hp_sum <- hydro |>
  group_by(RCP, Timeframe, month_f) |>
  summarise(mean_deficit_pct = mean(Percent_deficit), sd_deficit_pct = sd(Percent_deficit),
            .groups = "drop")

p_hp1 <- ggplot(hp_sum, aes(month_f, mean_deficit_pct, fill = Timeframe)) +
  geom_col(position = position_dodge(0.85), width = 0.8, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_deficit_pct - sd_deficit_pct, 0),
                    ymax = mean_deficit_pct + sd_deficit_pct),
                position = position_dodge(0.85), width = 0.3, linewidth = 0.3) +
  facet_wrap(~RCP) +
  scale_fill_manual(values = tf_cols) +
  labs(x = "Month", y = "Percent shortfall of firm hydropower production (%)") +
  theme_paper
save_png(p_hp1, "TRHRSP_hydropower_percent_deficit.png", 12, 7.5)

## Fig. A5: achieved firm energy = monthly target - mean shortfall
hp_ach <- hydro |>
  mutate(Achieved = Firm_thousand_GWh_week - Firmdeficit_mean_models / 1e6) |>
  group_by(RCP, Timeframe, month_f, Firm_thousand_GWh_week) |>
  summarise(Achieved = mean(Achieved), .groups = "drop")

p_hp2 <- ggplot(hp_ach, aes(month_f, Achieved, fill = Timeframe)) +
  geom_col(position = position_dodge(0.85), width = 0.8, colour = "black", linewidth = 0.2) +
  geom_line(aes(y = Firm_thousand_GWh_week, group = 1), colour = "purple", linewidth = 0.7) +
  geom_point(aes(y = Firm_thousand_GWh_week), colour = "purple", shape = 17, size = 1.8) +
  facet_wrap(~RCP) +
  scale_fill_manual(values = tf_cols) +
  labs(x = "Month", y = "Thousand GWh / week") +
  theme_paper
save_png(p_hp2, "TRHRSP_hydropower_achieved_firm_energy.png", 12, 7.5)

write.csv(left_join(hp_sum, select(hp_ach, RCP, Timeframe, month_f, Firm_thousand_GWh_week, Achieved),
                    by = c("RCP", "Timeframe", "month_f")),
          out("TRHRSP_hydropower_monthly_summary.csv"), row.names = FALSE)

## =============================================================================
## 3. DAILY DAM FLOWS (shared by flood and minimum flows)
## =============================================================================
DAMS <- c("WANAP", "WELLS", "CHIEF", "DALLE", "JDAYY", "MCNAR", "PRIRA", "RISLA", "ROCKY")

flows <- readRDS(F_DAMS) |>
  separate(Model_Folder, into = c("Model", "Scenario"), sep = "_", remove = FALSE) |>
  filter(!Model %in% EXCLUDED_GCMS, Scenario != "baseline", Category == PARTITION) |>
  select(date, month, week, day, year, Model, Scenario, Category, all_of(DAMS)) |>
  add_timeframe(closed_start = FALSE) |>
  pivot_longer(all_of(DAMS), names_to = "Dams", values_to = "Flow") |>
  mutate(Flow_cfs = Flow / 13.88)

## =============================================================================
## 4. FLOOD RISK AT THE DALLES  (Section 3.3.1; Fig. 6)
## =============================================================================
dalles <- flows |>
  filter(Dams == "DALLE") |>
  group_by(year, month, Model, Scenario, Timeframe) |>
  summarise(meanflow = mean(Flow_cfs), .groups = "drop") |>
  mutate(RCP = rcp_label(Scenario)) |>
  duplicate_historical()

stages <- data.frame(stage = c("Action stage", "Damage stage"), y = c(ACTION_CFS, DAMAGE_CFS))

p_flood <- ggplot(dalles, aes(factor(month, levels = 1:12, labels = month.abb), meanflow, fill = Timeframe)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, position = position_dodge2(preserve = "single")) +
  geom_hline(data = stages, aes(yintercept = y, colour = stage, linetype = stage)) +
  facet_wrap(~RCP, nrow = 1) +
  scale_fill_manual(values = tf_cols, name = NULL) +
  scale_colour_manual(values = c("Action stage" = "red", "Damage stage" = "blue"), name = NULL) +
  scale_linetype_manual(values = c("Action stage" = "dashed", "Damage stage" = "dotted"), name = NULL) +
  labs(x = "Month", y = "Mean flow (cfs)") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.key.width = unit(1.5, "cm"))
save_png(p_flood, "TRHRSP_flooding_TheDalles.png", 12, 7.5)

## Numbers for the text. Unit = one GCM-year of monthly mean flow.
## Historical rows are identical in the two RCP blocks (duplicated for plotting).
flood_stats <- dalles |>
  group_by(RCP, Timeframe, month) |>
  summarise(n_model_years = n(),
            median_cfs = median(meanflow),
            q25_cfs = quantile(meanflow, 0.25), q75_cfs = quantile(meanflow, 0.75),
            p95_cfs = quantile(meanflow, 0.95), max_cfs = max(meanflow),
            pct_ge_action = 100 * mean(meanflow >= ACTION_CFS),
            pct_ge_damage = 100 * mean(meanflow >= DAMAGE_CFS), .groups = "drop") |>
  mutate(month_name = month.abb[month])

flood_yearly <- dalles |>       # a YEAR counts if ANY GCM exceeds / if the ensemble median exceeds
  group_by(RCP, Timeframe, month, year) |>
  summarise(any_ge_action = any(meanflow >= ACTION_CFS), any_ge_damage = any(meanflow >= DAMAGE_CFS),
            ens_median_cfs = median(meanflow), .groups = "drop") |>
  group_by(RCP, Timeframe, month) |>
  summarise(pct_years_any_gcm_action = 100 * mean(any_ge_action),
            pct_years_any_gcm_damage = 100 * mean(any_ge_damage),
            pct_years_ens_median_action = 100 * mean(ens_median_cfs >= ACTION_CFS), .groups = "drop") |>
  mutate(month_name = month.abb[month])

write.csv(flood_stats,  out("TRHRSP_flood_stats_TheDalles.csv"), row.names = FALSE)
write.csv(flood_yearly, out("TRHRSP_flood_yearly_exceedance_TheDalles.csv"), row.names = FALSE)

## Printed summaries used in the 3.3.1 paragraph
print(flood_stats |> filter(month %in% 2:7) |> select(RCP, Timeframe, month_name, median_cfs) |>
        mutate(median_cfs = round(median_cfs, -3)) |>
        pivot_wider(names_from = month_name, values_from = median_cfs) |> arrange(RCP, Timeframe),
      width = Inf)
print(flood_stats |> filter(pct_ge_action > 0 | pct_ge_damage > 0) |> group_by(RCP, Timeframe) |>
        summarise(months_with_action = paste(month_name[pct_ge_action > 0], collapse = ", "),
                  months_with_damage = paste(month_name[pct_ge_damage > 0], collapse = ", "),
                  .groups = "drop"), width = Inf)

## =============================================================================
## 5. MINIMUM FLOWS  (Section 3.3.3; Table 1; Figs. 8 and A6-A13)  -- MONTHLY
## =============================================================================
rules <- read.csv(F_MINFLW) |>
  pivot_longer(cols = 4:12, names_to = "Dams", values_to = "Lowflow") |>   # columns 4-12 = nine dams
  mutate(Lowflow = as.numeric(gsub(",", "", Lowflow)),
         n_days  = day_end - day_start + 1)

## 5a. One target per dam x month. Only Apr, Jun, Jul, Oct have several windows.
targets <- rules |>
  group_by(Dams, month) |>
  summarise(n_windows = n(),
            target_dayweighted = weighted.mean(Lowflow, n_days),   # each window counts by its days
            target_simple      = mean(Lowflow),                    # each window counts equally
            .groups = "drop") |>
  mutate(Lowflow = if (TARGET_AVERAGE == "day_weighted") target_dayweighted else target_simple)
write.csv(targets, out("TRHRSP_minflow_monthly_targets.csv"), row.names = FALSE)

## 5b. Mean flow of each year-month (days inside the rule windows = the whole month)
##     NOTE: grouping omits Model, so each year's flow is the mean over GCMs and days,
##     as in the original code. Add Model to group_by()/n_distinct() to judge each GCM separately.
yearly <- flows |>
  inner_join(distinct(rules, month, day_start, day_end), by = "month") |>   # many-to-many warning here is expected
  filter(day >= day_start, day <= day_end) |>
  group_by(year, month, Scenario, Timeframe, Dams) |>
  summarise(meanflow = mean(Flow_cfs), .groups = "drop") |>
  left_join(select(targets, Dams, month, Lowflow), by = c("Dams", "month")) |>
  mutate(flag = as.integer(meanflow < Lowflow),
         pct_short = (Lowflow - meanflow) / Lowflow * 100)
write.csv(yearly, out("TRHRSP_minflow_yearly.csv"), row.names = FALSE)

## 5b-2. HISTORICAL CHECK (for review with Matt): the same test at four levels
##   unit   = pooled over GCMs (one value per year)   vs   each GCM separately (GCM-years)
##   period = whole month (day-weighted target)       vs   each Table A2 window (biweekly target)
hist_flows <- flows |> filter(Timeframe == "Historical") |>
  inner_join(rules, by = c("Dams", "month"), relationship = "many-to-many") |>   # adds day_start, day_end, Lowflow (window)
  filter(day >= day_start, day <= day_end)
mon_tgt <- select(targets, Dams, month, target = Lowflow)

fail_stats <- function(df, level, keys)
  df |> mutate(fail = meanflow < target, short = (target - meanflow) / target * 100) |>
    group_by(across(all_of(keys))) |>
    summarise(n_units = n(), target_cfs = mean(target), mean_flow_cfs = mean(meanflow),
              worst_flow_cfs = min(meanflow), worst_margin_pct = min(-short),   # < 0 = below target
              n_fail = sum(fail), pct_fail = 100 * mean(fail),
              mean_shortfall_pct = ifelse(n_fail > 0, mean(short[fail]), 0), .groups = "drop") |>
    mutate(level = level)

win_keys <- c("Dams", "month", "day_start", "day_end")
hist_check <- bind_rows(
  hist_flows |> group_by(year, month, Dams) |> summarise(meanflow = mean(Flow_cfs), .groups = "drop") |>
    left_join(mon_tgt, by = c("Dams", "month")) |>
    fail_stats("1 month | GCMs pooled (used in paper)", c("Dams", "month")),
  hist_flows |> group_by(year, month, Dams, day_start, day_end, target = Lowflow) |>
    summarise(meanflow = mean(Flow_cfs), .groups = "drop") |>
    fail_stats("2 biweekly window | GCMs pooled (old method)", win_keys),
  hist_flows |> group_by(year, month, Model, Dams) |> summarise(meanflow = mean(Flow_cfs), .groups = "drop") |>
    left_join(mon_tgt, by = c("Dams", "month")) |>
    fail_stats("3 month | each GCM separately", c("Dams", "month")),
  hist_flows |> group_by(year, month, Model, Dams, day_start, day_end, target = Lowflow) |>
    summarise(meanflow = mean(Flow_cfs), .groups = "drop") |>
    fail_stats("4 biweekly window | each GCM separately", win_keys)
) |> select(level, Dams, month, day_start, day_end, everything()) |> arrange(level, Dams, month, day_start)
write.csv(hist_check, out("TRHRSP_minflow_historical_check.csv"), row.names = FALSE)
print(hist_check |> group_by(level) |>
        summarise(units = sum(n_units), failures = sum(n_fail), pct = round(100 * failures / units, 2),
                  worst_margin_pct = round(min(worst_margin_pct), 1)), width = Inf)

## 5c. Frequency (% of years failing) and magnitude (mean % shortfall in failing years)
row_levels <- c("Historical", "2040s RCP4.5", "2040s RCP8.5", "2080s RCP4.5", "2080s RCP8.5")  # bottom -> top
mf <- yearly |>
  group_by(Dams, month, Scenario, Timeframe) |>
  summarise(n_years = n_distinct(year), n_flag = sum(flag),
            mean_short = ifelse(n_flag > 0, mean(pct_short[flag == 1]), 0), .groups = "drop") |>
  mutate(freq = n_flag / ifelse(USE_ACTUAL_YEARS, n_years, 31) * 100,
         row_label = factor(ifelse(Timeframe == "Historical", "Historical",
                                   paste(Timeframe, rcp_label(Scenario))), levels = row_levels),
         month_f = factor(month.abb[month], levels = month.abb))

write.csv(mf |> select(Dams, month, row_label, n_years, n_flag, freq, mean_short) |>
            arrange(Dams, row_label, month),
          out("TRHRSP_minflow_monthly_metrics.csv"), row.names = FALSE)

## 5d. Heatmaps, one file per dam (frequency on top, shortfall below). Zero cells are grey.
heat_cols <- c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000")

heat <- function(df, val, title, lim_name) {
  df <- mutate(df, is_zero = .data[[val]] == 0)
  ggplot(df, aes(month_f, row_label)) +
    geom_tile(data = filter(df, is_zero), fill = "grey96", colour = "grey80", linewidth = 0.3) +
    geom_tile(data = filter(df, !is_zero), aes(fill = .data[[val]]), colour = "grey80", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is_zero, "0", round(.data[[val]])), colour = is_zero),
              size = 5, fontface = "bold") +
    scale_colour_manual(values = c("TRUE" = "grey50", "FALSE" = "black"), guide = "none") +
    scale_fill_gradientn(colours = heat_cols, name = lim_name, limits = c(0, NA)) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +   # keep all months/rows in order
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "grey60", fill = NA),
          axis.text = element_text(face = "bold", colour = "black", size = 14),
          plot.title = element_text(face = "bold", size = 16),
          legend.title = element_text(face = "bold", size = 13),
          legend.text = element_text(size = 12))
}

for (d in unique(mf$Dams)) {
  sub <- filter(mf, Dams == d)
  p <- heat(sub, "freq", paste0(d, " - Frequency (% years not meeting minimum flow)"), "Frequency (%)") /
       heat(sub, "mean_short", paste0(d, " - Mean shortfall when minimum flow not met"), "Mean shortfall (%)")
  save_png(p, paste0("TRHRSP_minflow_", d, ".png"), 12, 9)
}

message("Done. Figures and CSVs in: ", OUT_DIR)
