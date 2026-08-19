

# Global ------------------------------------------------------------------

# reading KE-DHS data & trying vaccine timeliness on Measles
pacman::p_load(posterior, tidybayes, rstanarm, marginaleffects, data.table, brms,
               purrr, dplyr, haven, ggplot2, janitor, lubridate, stringr , survival,
               ggsurvfit, icenReg, sf, kableExtra, tidyr, viridisLite, autoReg, flexsurv,
               survey, lubridate, here)
mvs <- naniar::miss_var_summary
source('../R/autoReg-modifier.R')

# Data --------------------------------------------------------------------

# DHS Geo data - spatial join — same CRS etc.
g1 <- st_read("../data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST)
shp <- readRDS('../data/shp/gadm/gadm41_NGA_2_pk.rds') |> terra::unwrap() |>
  st_as_sf() |> select(GID_2)
g2 <- st_join(g1, shp[, "GID_2"], left = TRUE) |> st_drop_geometry()

# ldata - vax-data | cdata - covariates
ldata <- readRDS('../data/processed/vaxdata-components.rds')
cdata <- readRDS('../data/processed/dhs-covariates.rds')
master_data <- readRDS("../data/processed/master-survey-dataset.rds")

# processed temperature data - at pixel level (buffered)
cds_geoloc <- arrow::read_parquet('../data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

# processed data on the heat index
hi_data <- arrow::read_parquet('../data/processed/heatindex-processed.parquet') |>
  mutate(cluster = as.character(cluster))
setDT(hi_data)

# Table 1 -----------------------------------------------------------------

pacman::p_load(gtsummary, dplyr)

# data
antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1')
usenames <- c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')

ldata <- readRDS('data/processed/vaxdata-components.rds')
cdata <- readRDS('data/processed/dhs-covariates.rds')
use_ldata <- ldata[antigen]

for (i in 1:length(use_ldata)) use_ldata[[i]] <- use_ldata[[i]] |> mutate(antigen = usenames[i])

d1 <- merge(bind_rows(use_ldata) |> mutate(across(c(wt, caseid, bidx), as.character))
            , cdata |> mutate(wt = as.character(wt)) |> select(-num_anc_visits),
            by = c('caseid', 'bidx', 'admin', 'cluster', 'wt', 'residence'), all.x = T)


d2 <- d1 |>
  transmute(

    Antigen = factor(antigen, labels = usenames, levels = usenames),

    Residence = factor(
      residence,
      levels = c("urban", "rural"),
      labels = c("Urban", "Rural")
    ),

    `Birth order` = birth_order,

    `Place of delivery` = factor(
      place_delivery,
      levels = c("home", "institution"),
      labels = c("Home", "Health facility")
    ),

    `Child's gender` = factor(
      child_gender,
      levels = c("male", "female"),
      labels = c("Male", "Female")
    ),

    `Travel time to the\nnearest health facility` = factor(
      time_to_hf,
      levels = c("<30 mins", "31-60 mins", "1-2 hrs", "2+ hrs"),
      labels = c(
        "<30 minutes",
        "31–60 minutes",
        "1–2 hours",
        ">2 hours"
      )
    ),

    `Wealth index` = factor(
      wealth_index,
      levels = c("poorest", "poorer", "middle", "richer", "richest"),
      labels = c("Poorest", "Poorer", "Middle", "Richer", "Richest")
    ),

    `Highest level of\neducation attained by mother` = factor(
      meduc,
      levels = c("no education", "primary", "secondary", "higher"),
      labels = c(
        "No formal education",
        "Primary school",
        "Secondary school",
        "Higher education"
      )
    ),

    `Mother's age group\nat birth` = case_when(
      mother_age_birth <= 19 ~ "<=19 yrs",
      mother_age_birth >= 20 & mother_age_birth <= 24 ~ "20-24 yrs",
      mother_age_birth >= 25 & mother_age_birth <= 29 ~ "25-29 yrs",
      mother_age_birth >= 30 & mother_age_birth <= 34 ~ "30-34 yrs",
      mother_age_birth >= 35 & mother_age_birth <= 39 ~ "35-39 yrs",
      mother_age_birth >= 40 & mother_age_birth <= 44 ~ "40-44 yrs",
      mother_age_birth >= 45 ~ "45+ yrs",
      TRUE ~ NA_character_
    ),
    `Mother's age group\nat birth` = factor(
      `Mother's age group\nat birth`,
      levels = c("<=19 yrs", "20-24 yrs","25-29 yrs","30-34 yrs",
                 "35-39 yrs","40-44 yrs","45+ yrs"),
      labels = c("≤19 years","20–24 years","25–29 years","30–34 years",
                 "35–39 years","40–44 years","45 years or older")
    ),

    `Geopolitical zone` = factor(geozone,
                                 levels = c("north west", "north east", "north central", "south east",
                                            "south south", "south west"),
                                 labels = c("North West", "North East", "North Central", "South East",
                                            "South South", "South West"))
  )

d2 |>
  tbl_summary(
    by = Antigen,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      `Birth order` ~ "{mean} ({sd})"
    ),
    missing = "no",
    label = list(
      Residence ~ "Residence",
      `Birth order` ~ "Birth order",
      `Place of delivery` ~ "Place of delivery",
      `Child's gender` ~ "Child's gender",
      `Travel time to the\nnearest health facility` ~ "Travel time to nearest health facility",
      `Wealth index` ~ "Wealth index",
      `Highest level of\neducation attained by mother` ~ "Maternal education",
      `Mother's age group\nat birth` ~ "Mother's age at birth",
      `Geopolitical zone` ~ "Geopolitical zone"
    )
  ) |>
  bold_labels() |>
  modify_header(label ~ "**Characteristic**") |>
  modify_spanning_header(all_stat_cols() ~ "**Antigen**")




# Heat Index plot ---------------------------------------------------------

# shapefile
admin2_res <- bind_cols(
  data.frame(st_centroid(g1) |> st_coordinates()),
  cluster = g1$cluster
) |>
  setNames(c('lon', 'lat', 'cluster')) |>
  merge(hi_data |>
          mutate(
            heatclass = case_when(
              heatindex < 26.67                       ~ "Tolerable",
              heatindex >= 26.67 & heatindex < 32.22  ~ "Caution",
              heatindex >= 32.22 & heatindex < 39.44  ~ "Extreme Caution",
              heatindex >= 39.44 & heatindex < 51.67  ~ "Danger",
              heatindex >= 51.67                      ~ "Extreme Danger",
              TRUE                                    ~ NA_character_
            ),
            heatclass = factor(
              heatclass,
              levels = c("Tolerable", "Caution", "Extreme Caution",
                         "Danger", "Extreme Danger"
              ),
              labels = c("Tolerable", "Caution", "Extreme Caution",
                         "Danger", "Extreme Danger"
              ),
              ordered = T
            )
          ),
        by = 'cluster', all = T)

admin2_res <- admin2_res |>
  mutate(
    cluster = factor(cluster,
                     levels = admin2_res |> distinct(cluster, lat) |> arrange(lat) |> pull(cluster)
    )
  )

# heat index classification
p <- ggplot(admin2_res) +
  geom_tile(aes(x = date, y = cluster, fill = heatclass)) +
  # Swapped continuous gradient for explicit NWS categorical palette
  scale_fill_manual(
    values = c(
      "Tolerable"       = "#a6d96a",  # Soft NWS Light Green
      "Caution"         = "#ffff00",  # NWS Yellow
      "Extreme Caution" = "#ffaa00",  # NWS Amber / Dark Yellow
      "Danger"          = "#ff6600",  # NWS Orange
      "Extreme Danger"  = "#cc0000"   # NWS Red
    ),
    name = "Heat Index\nCategory",
    drop = FALSE # Keeps all levels in the legend even if missing in data subset
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = levels(admin2_res$cluster)[c(1, length(levels(admin2_res$cluster)))],
    labels = c("Southern\nNigeria", "Northern\nNigeria"),
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, title = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.title.y = element_blank(),

    # Updated legend adjustments for discrete blocks
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.8, "cm"),

    axis.text.x = element_text(size = 10, colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),

    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  paste0("output/img/heatmaps/", "heat-index-class.png"),
  p,
  width = 12,
  height = 9,
  dpi = 1000
)



# Heat Index exposure -----------------------------------------------------

# data
antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1')
usenames <- c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')
use_ldata <- ldata[antigen];

for (i in 1:length(use_ldata)) use_ldata[[i]] <- use_ldata[[i]] |> mutate(antigen = usenames[i])

vax_data <- bind_rows(use_ldata)
use_window_min <- ifelse(vax_data$antigen == 'BCG', 0, 28)
setDT(vax_data)

# the 7-daywindow bound
vax_data[, `:=`(
  start_dt = due_date - use_window_min,
  end_dt   = due_date + 28,
  cluster  = as.character(cluster)
)]

heatwave <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                      .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
                      by = .EACHI] |>
  setNames(c('cluster', 'start_dt', 'end_dt', 'heatwave')) |>
  distinct()
heatindex <- hi_data[vax_data,
                     on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                     .(heatindex = mean(heatindex, na.rm = TRUE)),
                     by = .EACHI] |>
  setnames(c('cluster', 'start_dt', 'end_dt', 'heatindex')) |>
  unique()

# merging to have the heatwave and temperature variable
vax_data <- merge(
  vax_data, heatwave,
  by = c('cluster', 'start_dt', 'end_dt'),
  all.x = T
) |>
  merge(heatindex, by = c('cluster', 'start_dt', 'end_dt'), all.x = T) |>
  mutate(
    heatwave = ifelse(heatwave == 0, 'absent', 'present'),

    heatclass = case_when(
      heatindex < 26.67                       ~ "Tolerable",
      heatindex >= 26.67 & heatindex < 32.22  ~ "Caution",
      heatindex >= 32.22 & heatindex < 39.44  ~ "Extreme Caution",
      heatindex >= 39.44 & heatindex < 51.67  ~ "Danger",
      heatindex >= 51.67                      ~ "Extreme Danger",
      TRUE                                    ~ NA_character_
    ),
    heatclass = factor(
      heatclass,
      levels = c("Tolerable", "Caution", "Extreme Caution",
                 "Danger", "Extreme Danger"
      ),
      labels = c("Tolerable", "Caution", "Extreme Caution",
                 "Danger", "Extreme Danger"
      ),
      ordered = T
    )
    )


vax_data |>
  select(`Heatwave` = heatwave, `Heat Index` = heatclass, antigen)





# Kaplan meier analysis ---------------------------------------------------

antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1')
usenames <- c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')
antigen_times <- pmax(1, c(0, c(6, 10, 14)*7, 6*30.4, 9*30.4) - 7)
duedates <- antigen_times + 28 + c(0, rep(7, 5))

plt <- list(); maxtime <- cov_duedate <- cov_3y <- median_aav <- vector()

for (i in 1:length(antigen)) {

  message('Currently: ', usenames[i])

  use_window_min <- ifelse(antigen[i] == 'bcg', 0, 7)
  vax_data <- ldata[[antigen[i]]]; setDT(vax_data)

  # the 7-daywindow bound
  vax_data[, `:=`(
    start_dt = due_date - use_window_min,
    end_dt   = due_date + 7,
    cluster  = as.character(cluster)
  )]
  vax_data <- vax_data |>
    mutate(
      wt = as.numeric(wt),
      delay = as.numeric(vaxx_date - due_date),
      delayclass = ifelse(delay > 28, 1, 0)
    )

  vax_data <- vax_data |> filter(delay >= -7 | is.na(delay))
  age_days <- as.numeric(vax_data$interview_date - vax_data$birth_date)
  age_id_omit <- which(age_days < (antigen_times[i]))
  {if (length(age_id_omit) == 1) {vax_data <- vax_data[-age_id_omit, ]}}


  vax_data_plot <- vax_data |>
    mutate(

      # 1. First, create the exact Surv-compatible time windows
      time_start = case_when(
        outcome_event == 0  ~ event_time,   # Exact event: starts at vaccine day
        outcome_event == -1 ~ antigen_times[i],            # Interval censored: starts at vaccine due date (7 days before schedule) #* was initially 0
        outcome_event == 1  ~ censor_time,  # Right censored: starts at interview day
        TRUE ~ NA_real_
      ),

      time_end = case_when(
        outcome_event == 0  ~ event_time,   # Exact event: ends at vaccine day
        outcome_event == -1 ~ censor_time,  # Interval censored: ends at interview day
        outcome_event == 1  ~ NA_real_,     # Right censored: open-ended upper bound
        TRUE ~ NA_real_
      ),

      # 2. Now map your outcome_event to R's Surv(..., type="interval") standards
      # 0 = right censored, 1 = exact event, 2 = left censored, 3 = interval censored
      status = case_when(
        outcome_event == 0  ~ 1,  # Exact event
        outcome_event == -1 ~ 3,  # Interval censored was 2=Left censored
        outcome_event == 1  ~ 0,  # Right censored
        TRUE ~ NA_real_
      )
    )

  {if (i == 3) vax_data_plot <- vax_data_plot |> mutate(time_end = ifelse(time_start > time_end, time_start + 1, time_end))}

  # BCG checks to modify time points for people receiving vax before b.day
  {if (antigen[i] == 'bcg') {
    vax_data_plot$time_start <- ifelse(vax_data_plot$time_start <= 0, 1, vax_data_plot$time_start)
    vax_data_plot$time_end <- ifelse(vax_data_plot$time_end <= 0, 1, vax_data_plot$time_end)
  }}

  survmodel <- survfit(Surv(time_start, time_end, status, type = "interval") ~ 1,
                      data = vax_data_plot)
  maxtime[i] <- max((summary(survmodel))$time)
  fit_summary_28d <- summary(survmodel, times = duedates[i])
  fit_summary_3y  <- summary(survmodel, times = maxtime[i])
  cov_duedate[i] <- round((1 - fit_summary_28d$surv) * 100, 2)
  cov_3y[i] <- round((1 - fit_summary_3y$surv) * 100, 2)
  median_aav[i] <- summary(survmodel)$table["median"] |> as.numeric()

  # plt[[i]] <- survmodel |>
  #   ggsurvfit(type = "risk") +
  #   add_confidence_interval() +
  #   add_censor_mark(shape = '+') +
  #   scale_ggsurvfit(
  #     x_scales = list(breaks = seq(0, 1095, 100)),
  #     y_scales = list(labels = scales::percent, limits = c(0, 1))
  #   ) +
  #   add_risktable(
  #     risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
  #                         "{format(round(n.event, 0), nsmall = 0)}"),
  #     stats_label = c("At risk",
  #                     "Events")
  #   ) +
  #   theme_ggsurvfit_KMunicate() +
  #   labs(
  #     title = usenames[i],
  #     x = "Age (in days)",
  #     y = "Cumulative Vaccination Coverage"
  #   )
  # ggsave(filename = paste0('../output/img/', antigen[i], "/overall.png"), plot = plt[[i]], height = 8, width = 8, dpi = 1000)


}

data.frame(
  `Vaccine` = usenames,
  `Median age at\nvaccination (days)` = median_aav,
  `Vaccination coverage\nwithin 28 days of due date` = paste0(cov_duedate, '%'),
  `Vaccination coverage\nat 3 years` = paste0(cov_3y, '%'),
  check.names = F
)

# plotting them all

ggsave(filename = paste0('../output/img/all-coverage.png'),
       plot = ((plt[[1]] | plt[[2]] | plt[[3]]) / (plt[[4]] | plt[[5]] | plt[[6]])),
       height = 6, width = 12, dpi = 1000)


# Heat and Timeliness -----------------------------------------------------

tbl <- list()
for (i in 1:length(antigen)) {
  f <- readRDS(paste0('../output/models/', antigen[i], '-file.rds'))

  tbl[[i]] <- lapply(f, \(x) {
    confint(x) |>
      data.frame() |>
      mutate(Est = coef(x), Var = names(coef(x))) |>
      setNames(c('lwr', 'upr', 'est', 'var'))
  }) |>
    bind_rows() |>
    filter(str_detect(var, 'heat')) |>
    mutate(vaccine = usenames[i])
  rownames(tbl[[i]]) <- 1:nrow(tbl[[i]])
}

tbl |>
  bind_rows() |>
  mutate(
    across(c(upr, lwr, est), exp),
    vaccine = factor(vaccine, labels = usenames, levels = usenames),
    model = ifelse(str_detect(var, 'wave'), 'Heatwave', 'Heat Index'),
    model = factor(model, labels = c('Heatwave', 'Heat Index'), levels = c('Heatwave', 'Heat Index')),
    var = case_when(
      str_detect(var, 'wave') ~ 'Heatwave',
      str_detect(var, 'Extreme') ~ 'Extreme Caution',
      T ~ 'Caution'
    )
  ) |>

  ggplot() +

  geom_point(aes(y = vaccine, x = est, col = var),
             shape = 15, size = 2,
             position = position_dodge2(width = 0.3, preserve = "total")) +
  geom_errorbar(aes(y = vaccine, xmin = lwr, xmax = upr, col = var),
                width = .3, position = position_dodge2(width = 0.3, preserve = "total")) +
  geom_vline(aes(xintercept = 1), col = 'black', linetype = 'dashed') +

  labs(
    x = 'Event Time Ratios (ETR)',
    y = NULL,
    col = 'Category\n(Ref: Tolerable / No Heatwave)'
  ) +
  facet_wrap(~model, nrow = 1, scales = "free_x") +

  scale_colour_manual(
    values = c(
      "Heatwave"        = "black",
      "Caution"         = "#fe9929",
      "Extreme Caution" = "#d95f02"
    ),
    breaks = c("Caution", "Extreme Caution")
  ) +
  theme_bw(base_family = 'Times New Roman') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(colour = "black", size = 10),
    axis.ticks = element_line(colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave('../output/img/heat-coef-plot.png',
       dpi = 1e3, height = 6, width = 12)


# avg marginal cum-vax-cov ------------------------------------------------

# main function for computing the standardizatiom
standardize_curve <- function(x, scenario_value = NULL, heat_var = "heatindex",
                              days, w) {
  cf_data <- x$model
  if (!is.null(scenario_value)) cf_data[[heat_var]] <- scenario_value

  lp     <- predict(x, newdata = cf_data, type = "lp")  # Xb per child
  shape  <- 1 / x$scale                                  # Weibull shape (check for strata!)
  scale_i <- exp(lp)                                     # Weibull scale per child

  sapply(days, function(t) {
    S_i <- pweibull(t, shape = shape, scale = scale_i)  # P(T <= t | X_i, scenario)
    weighted.mean(S_i, w = w, na.rm = TRUE) * 100
  })
}

# curvs for any given number of days
curves <- map_dfr(c("Tolerable", "Caution", "Extreme Caution"), function(s) {
  data.frame(
    Scenario = s,
    Day = 1:1095,
    Coverage = standardize_curve(x, scenario_value = s, days = 1:1095, w = w)
  )
})

ggplot(curves) + geom_line(aes(x = Day, y = Coverage, col =Scenario))

# coverage as at timeliness thresholds
mvxd <- avxd <- list()
antigen_times <- pmax(1, c(0, c(6, 10, 14)*7, 6*30.4, 9*30.4))

for (i in 1:length(antigen)) {

  # load up the model
  f <- readRDS(paste0('../output/models/', antigen[i], '-file.rds'))
  f_hi <- f$fit_hi; f_hw <- f$fit_hw
  w_hi <- weights(f_hi$survey.design);  w_hw <- weights(f_hw$survey.design)

  target_milestones <- (7 * (1:4)) + antigen_times[i]

  # heat index: marginals
  m_vx_hi <- map_dfr(c("Tolerable", "Caution", "Extreme Caution"), function(s) {
    data.frame(Scenario = s, Target_Day = target_milestones,
               Coverage = standardize_curve(f_hi, scenario_value = s,
                                            days = target_milestones, w = w_hi))
  }) |>
    mutate(Exposure = 'Heat Index')

  # heatwaves: marginals
  m_vx_hw <- map_dfr(c('No heatwave', 'Heatwave'), function(s) {
    data.frame(Scenario = s, Target_Day = target_milestones,
               Coverage = standardize_curve(f_hw, scenario_value = s,
                                            days = target_milestones, w = w_hw))
  }) |>
    mutate(Exposure = 'Heatwave')

  mvxd[[i]] <- bind_rows(m_vx_hi, m_vx_hw) |> mutate(Antigen = usenames[i])

  # heatwave and heat index: actual data
  a_vx_hi <- data.frame(Target_Day = 1:1094,
                        Coverage = standardize_curve(f_hi, scenario_value = NULL,
                                                     days = 1:1094, w = w_hi)) |>
    mutate(Exposure = 'Heat Index')

  # heatwaves: marginals
  a_vx_hw <- data.frame(Target_Day = 1:1094,
                        Coverage = standardize_curve(f_hw, scenario_value = NULL,
                                                     days = 1:1094, w = w_hw)) |>
    mutate(Exposure = 'Heatwave')

  avxd[[i]] <- bind_rows(a_vx_hi, a_vx_hw) |> mutate(Antigen = usenames[i])

}

# Marginal: Model-standardised predicted vaccination coverage by timeliness threshold under different heat-index scenarios
p1 <- bind_rows(mvxd) |>
  filter(Exposure == 'Heat Index') |>
  mutate(
    Scenario = factor(Scenario,
                      levels = c("Tolerable", "Caution", "Extreme Caution"),
                      labels = c("Tolerable", "Caution", "Extreme Caution")),
    DueDate = case_when(
      Antigen == 'BCG' ~ antigen_times[1],
      Antigen == 'Pentavalent 1' ~ antigen_times[2],
      Antigen == 'Pentavalent 2' ~ antigen_times[3],
      Antigen == 'Pentavalent 3' ~ antigen_times[4],
      Antigen == 'Vitamin A' ~ antigen_times[5],
      T ~ antigen_times[6]
    ),
    Target_Day = Target_Day - DueDate,
    Antigen = factor(Antigen, levels = usenames, labels = usenames),
    Percent_Label = paste0(round(Coverage, 0), '%'),
    Coverage = Coverage / 100
  ) |>


  ggplot(aes(x = factor(Target_Day), y = Scenario,)) +
  geom_tile(aes(fill = Coverage)) +
  geom_text(aes(label = Percent_Label), color = "black", size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", labels = scales::percent) +
  labs(subtitle = 'Heat Index exposure', x = 'Timeliness threshold\n(days past recommended vccination date)', y = NULL) +
  facet_wrap(~Antigen) +
  theme_bw(base_family = 'Times New Roman')

p2 <- bind_rows(mvxd) |>
  filter(Exposure != 'Heat Index') |>
  mutate(
    Scenario = factor(Scenario,
                      levels = c('No heatwave', 'Heatwave'),
                      labels = c('No heatwave', 'Heatwave')),
    DueDate = case_when(
      Antigen == 'BCG' ~ antigen_times[1],
      Antigen == 'Pentavalent 1' ~ antigen_times[2],
      Antigen == 'Pentavalent 2' ~ antigen_times[3],
      Antigen == 'Pentavalent 3' ~ antigen_times[4],
      Antigen == 'Vitamin A' ~ antigen_times[5],
      T ~ antigen_times[6]
    ),
    Target_Day = Target_Day - DueDate,
    Antigen = factor(Antigen, levels = usenames, labels = usenames),
    Percent_Label = paste0(round(Coverage, 0), '%'),
    Coverage = Coverage / 100
  ) |>


  ggplot(aes(x = factor(Target_Day), y = Scenario,)) +
  geom_tile(aes(fill = Coverage)) +
  geom_text(aes(label = Percent_Label), color = "black", size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", labels = scales::percent) +
  labs(subtitle = 'Heatwave exposure', x = 'Timeliness threshold\n(days past recommended vccination date)', y = NULL) +
  facet_wrap(~Antigen) +
  theme_bw(base_family = 'Times New Roman')

patchwork::wrap_plots(p1, p2, nrow = 2, heights = c(2, 1))
ggsave('../output/img/marginal-preds.png', height = 8, width = 12, dpi = 1e3)


# Performance comparison --------------------------------------------------

avxd |>
  bind_rows() |>
  mutate(
    Antigen = factor(Antigen, levels = usenames, labels = usenames)
  ) |>

  ggplot() +
  geom_step(aes(x = Target_Day, y = Coverage, col = Antigen)) +
  geom_vline(aes(xintercept = 28), col = 'black', lty = 'dashed') +
  labs(x = 'Age (in days)', x = 'Vaccine coverage') +
  facet_wrap(~Exposure) +
  theme_bw(base_family = 'Times New Roman')

