
# reading KE-DHS data & trying vaccine timeliness on Measles
pacman::p_load(posterior, tidybayes, rstanarm, marginaleffects, data.table, brms,
               ggpubr, dplyr, haven, ggplot2, janitor, lubridate, stringr , survival,
               ggsurvfit, icenReg, sf)
mvs <- naniar::miss_var_summary

# Data --------------------------------------------------------------------

# DHS Geo data - spatial join — same CRS etc.
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST)
shp <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> terra::unwrap() |>
  st_as_sf() |> select(GID_2)
g2 <- st_join(g1, shp[, "GID_2"], left = TRUE) |> st_drop_geometry()

# ldata - vax-data | cdata - covariates
ldata <- readRDS('data/processed/vaxdata-components.rds')
cdata <- readRDS('data/processed/dhs-covariates.rds')

# processed temperature data - at pixel level (buffered)
cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .9
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

# IN CASE WE NEED IT: processed temperature data - at areal level (zonal aggregation)
# cds_areal <- arrow::read_parquet('data/processed/admin2-processed.parquet')
# cds_areal$heatwave <- cds_areal$p >= .9
#
# cds_areal <- merge(g2, cds_areal, by = 'GID_2', all.x = T)
# cds_areal <- cds_areal |> mutate(cluster = as.character(cluster))
# setDT(cds_areal); setkey(cds_areal, cluster, date)


# Merging -----------------------------------------------------------------

vax_data <- ldata[['penta1']]; setDT(vax_data)
vax_used_folder <- 'output/img/penta1/'

# the 14-day window bound
vax_data[, `:=`(
  start_dt = due_date - 28, #* change to 7
  end_dt   = due_date + 28,
  cluster  = as.character(cluster)
)]

# find all rows in cds_geoloc where cluster matches AND date is between start/end
# and sum the 'heatwave' column for @ child
results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
               .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
               by = .EACHI]

vax_data$heatwave <- ifelse(results$heatwave_sum == 0, 'absent', 'present')
(table(vax_data$heatwave))
(table(is.na(vax_data$vaxx_date), vax_data$heatwave))

# merging with covariates
vax_data <- merge(vax_data |> mutate(across(c(wt, caseid, bidx), as.character))
                  , cdata |> mutate(wt = as.character(wt)),
                  by = c('caseid', 'bidx', 'admin', 'cluster', 'wt', 'residence'), all.x = T)

vax_data <- vax_data |> filter(!is.na(num_anc_visits)) |>
  mutate(wt = as.numeric(wt),
         # across(where(is.factor), as.character)
         )

# 0 - event | -1 - left censored | 1 - right censored

# Data exploration --------------------------------------------------------

glimpse(vax_data)

vax_data_plot <- vax_data |>
  mutate(
    outcome_event = case_when(
      outcome_event == 0 ~ 1,
      outcome_event == 1 ~ 0,
      outcome_event == -1 ~ 3 # 2 is left censored, and 3 is interval
    ),

    time_start = case_when(
      outcome_event == 2 ~ (30.4*9)-7,            # Left: happened between birth and now / penta 3 always given as at 14*7 days / vit a: (6*30)-7
      outcome_event == 0 ~ time_outcome,    # Right: started at interview, ends never
      outcome_event == 1 ~ time_outcome     # Exact: happened at this time
    ),

    time_end = case_when(
      outcome_event == 2 ~ time_outcome,    # Left: upper bound is interview
      outcome_event == 0 ~ NA_real_,        # Right: no upper bound
      outcome_event == 1 ~ time_outcome     # Exact: same as start
    )
  ) |>
  rename(
    "Heatwave" = heatwave,
    "BirthOrder" = birth_order,
    "ANCVisits" = num_anc_visits,
    "PlaceDelivery" = place_delivery,
    "ChildGender" = child_gender,
    "MotherOccupation" = mother_occupation,
    "TimeHealthFacility" = time_to_hf,
    "WealthIndex" = wealth_index,
    "MotherEducation" = meduc,
    "MotherAge" = mother_age_group
  )

# kaplan meier plots based on: heatwave, birth order, anc, delivery, gender, occupation, time to hf, wealth, education, mother age
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ 1,
                   data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  labs(
    title = "Cumulative Vaccination Coverage",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )
ggsave(filename = paste0(vax_used_folder, "overall.png"), height = 8, width = 8, dpi = 1000)

# heatwaves
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ Heatwave,
                   data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Heatwave",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "heatwave.png"), height = 8, width = 8, dpi = 1000)

# Birth Order
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ ANCVisits,
                   data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By ANC Visits",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "ANC Visits.png"), height = 10, width = 10, dpi = 1000)

# place of delivery
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ PlaceDelivery,
                   data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Place of Delivery",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Place of delivery.png"), height = 10, width = 10, dpi = 1000)

# Child Gender
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ ChildGender,
                   data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Child's Gender",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Child Gender.png"), height = 10, width = 10, dpi = 1000)

# Mother Occupation
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ MotherOccupation,
                   data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  # add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Mother's Occupation",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Mother Occupation.png"), height = 10, width = 10, dpi = 1000)


# Time to Health Facility
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ TimeHealthFacility,
        data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  # add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Time to Health Facility",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Time to Health Facility.png"), height = 10, width = 10, dpi = 1000)


# Wealth Index
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ WealthIndex,
        data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  # add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Wealth Index",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Wealth Index.png"), height = 10, width = 10, dpi = 1000)


# Mother Education
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ MotherEducation,
        data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Mother's Eduction Attainment",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Mother Education.png"), height = 10, width = 10, dpi = 1000)


# Mother Education
survfit(Surv(time_start, time_end, outcome_event, type = "interval") ~ MotherAge,
        data = vax_data_plot) |>
  ggsurvfit(type = "risk") +
  # add_confidence_interval() +
  add_censor_mark(shape = '+') +
  scale_ggsurvfit(
    x_scales = list(breaks = seq(0, 1050, 100)),
    y_scales = list(labels = scales::percent, limits = c(0, 1))
  ) +
  add_risktable(
    risktable_stats = c("{format(round(n.risk, 0), nsmall = 0)}",
                        "{format(round(n.event, 0), nsmall = 0)}"),
    stats_label = c("At risk",
                    "Events")
  ) +
  theme_ggsurvfit_KMunicate() +
  theme(legend.position = 'bottom') +
  labs(
    title = "Cumulative Vaccination Coverage By Mother's Age",
    x = "Time since birth (days)",
    y = "Probability of Being Vaccinated"
  )

ggsave(filename = paste0(vax_used_folder, "Mother Age.png"), height = 10, width = 10, dpi = 1000)




# Model -------------------------------------------------------------------

b1 <- brm(
  time_outcome | weights(wt) + cens(outcome_event) ~

    heatwave + residence + admin +

    birth_order + num_anc_visits + place_delivery + child_gender +
    mother_occupation + time_to_hf + wealth_index + meduc + mother_age_birth,
    #(1 | caseid),

  family = weibull(),
  data = vax_data, # |> mutate(time_outcome = ifelse(time_outcome == 0, 1, time_outcome)),

  chains = 4,
  iter = 2000,
  warmup = 1000
)

summary(b2)
conditional_effects(b2, "heatwave")
bayesplot::pp_check(b1)

# conditional effects for timeliness
# type response - ignores individual residual randomness | prediction - incorporates it
# rvar format for posterior package
pred <- avg_predictions(b1, newdata = vax_data,
                        re_formula = NULL, type = 'response', by = "heatwave", wts = 'wt')
draws <- get_draws(pred, "rvar")

quantile2(draws$rvar, c(0.025, 0.5, 0.975)) # posterior quantiles
E(draws$rvar) # expected value of posterior dist
Pr(draws$rvar <= ((14+2)*7)) # posterior mass below 9 months + 2 week buffer

# Posterior predictions and quantiles
pp <- posterior_predict(b1, newdata = vax_data, re_formula = NULL)
tvax_mean <- colMeans(pp)
aggregate(tvax_mean ~ vax_data$heatwave, FUN = mean)

timely_draws <- pp <= ((14+2)*7)
timely_prob <- colMeans(timely_draws)
aggregate(timely_prob ~ vax_data$heatwave, FUN = mean)


predictions(
  b2,
  newdata = datagrid(heatwave = c('absent', 'present'),
                     caseid = unique),
  by = "heatwave",
  re_formula = NULL
)
