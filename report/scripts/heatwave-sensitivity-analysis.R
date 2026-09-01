

# globals -----------------------------------------------------------------

pacman::p_load(posterior, tidybayes, rstanarm, marginaleffects, data.table, brms,
               purrr, dplyr, haven, ggplot2, janitor, lubridate, stringr , survival,
               ggsurvfit, icenReg, sf, kableExtra, tidyr, viridisLite, autoReg, flexsurv,
               survey, here, terra, zoo, mev, gamlss, gamlss.dist, gamlssx, progressr)
mvs <- naniar::miss_var_summary
handlers(global = TRUE)

# vaxx data
# ldata - vax-data | cdata - covariates
ldata <- readRDS('../data/processed/vaxdata-components.rds')
cdata <- readRDS('../data/processed/dhs-covariates.rds')
master_data <- readRDS("../data/processed/master-survey-dataset.rds")

# the grid
grd <- expand.grid(
  rolling = c(3, 5, 7),
  quantile = 4 * c(1, 2, 3),
  window = c(3, 7, 14)
)

# the GMST covariate
gmst <- read.csv('../data/GMST/tas_ERA5.csv')[-c(1:30), 2:3] |>
  setNames(c('year', 'gmst')) |>
  filter(year %in% 1940:2025) |>
  mutate(across(everything(), as.numeric))

# the cluster geo codes
g1 <- st_read("../data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)
bf_gps <- g1 |> st_transform(3857) |> st_buffer(dist = 10000) |> st_transform(4326)

shp <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> terra::unwrap() |>
  st_as_sf() |> select(GID_2)
g2 <- st_join(g1, shp[, "GID_2"], left = TRUE) |> st_drop_geometry()

full_dates <- seq(from=ymd('1940-01-01'), to = ymd('2024-12-31'), by = 1)

ls_full <- list.files('../data/cds/', pattern = "*.nc") |> sort()
ls_full <- paste0("../data/cds/", ls_full)

# full-stack 1940-2019
t1 <- rast(ls_full)
t1 <- t1 - 273.15

t2_master <- terra::extract(t1, vect(bf_gps), fun = mean, na.rm = TRUE, ID = FALSE) |>
  setNames(paste0('dt.', full_dates)) |>
  mutate(
    cluster = bf_gps |> pull(cluster),
    residence = bf_gps |> pull(residence)
  ) |>
  pivot_longer(-c(cluster, residence), names_to = "layer", values_to = "temp") |>
  mutate(layer = ymd(substr(layer, 4, 13))) |>
  group_by(cluster) |>
  group_split()

# grid & loop -------------------------------------------------------------

full_data <- data.frame()

for (g in 1:nrow(grd)) {

  rolling_window <- grd[g, 'rolling'] - 1
  use_quantile <- (grd[g, 'quantile'] - 1) / grd[g, 'quantile']
  exposure_window <- grd[g, 'window']

  message('Iter: ', g, ' | Rolling window: ', rolling_window + 1, ' | Quantile: ', use_quantile, ' | Exposure window: ', exposure_window)

  # Heatwave construction ---------------------------------------------------
  {
    # computation of cluster TX5X index
    t2 <- t2_master |>
      lapply(
        FUN = \(dt) {
          dt |>
            arrange(layer) |>
            mutate(
              year = year(layer),
              tx5x = slider::slide_dbl(temp, mean, .before = rolling_window, .complete = TRUE),
              tx5x = na.locf(tx5x, fromLast = T)
            ) |>
            select(cluster, residence, year, date = layer, tx5x)
        }
      ) |>
      bind_rows()

    # construction of the training and testing set
    train <- t2 |>
      filter(year %in% 1940:2019) |>
      group_by(cluster, residence, year) |>
      reframe(mx = max(tx5x))
    test <- t2 |> filter(year %in% 2020:2024)
    rm(t2)

    d3 <- train |> merge(gmst, by = 'year', all.x = T)
    split_data <- d3 |> group_by(cluster, residence) |> group_split()

    with_progress({
      p <- progressor(steps = length(split_data))

      for (i in 1:length(split_data)) {

        grp_data <- split_data[[i]]
        tmp_model <- fitGEV(mx ~ gmst, data = grp_data, n.cyc = 500)
        tmp_data <- data.frame(
          cluster = unique(grp_data$cluster), residence = unique(grp_data$residence),
          icpt  = as.numeric(tmp_model$mu.coefficients[1]),
          coeftime  = as.numeric(tmp_model$mu.coefficients[2]),
          scale = exp(as.numeric(tmp_model$sigma.coefficients)),
          shape = as.numeric(tmp_model$nu.coefficients)
        )

        if (i == 1) nonstationary_est <- tmp_data
        else nonstationary_est <- bind_rows(nonstationary_est, tmp_data)

        p(message = paste("Doing step", i))
      }

    })


    # merging the model fit with the test data
    s1 <- merge(
      merge(test, gmst, by = 'year', all.x = T),
      nonstationary_est, by = c('cluster', 'residence'), all.x = T)

    # computing Pr(Y <= y) and return period. Also heatwave = p > 0.90
    s2 <- s1 |>
      mutate(
        # the eta/location is b0 + b1*gmst
        loc = icpt + (coeftime * gmst),

        p = mev::pgev(q = tx5x, loc = loc, scale = scale, shape = shape, lower.tail = T),
        return = 1 / (1 - p),

        # 1-in-10 year return level threshold (90th percentile of annual max distribution)
        thresh = mev::qgev(p = rep(use_quantile, nrow(s1)), loc = loc, scale = scale, shape = shape),
        heatwave = tx5x >= thresh
      ) |>
      select(cluster, year, date, tx5x, icpt, coeftime, loc, scale, shape, p, return, thresh, heatwave)
  }


  # vax data modelling ------------------------------------------------------

  cds_geoloc <- s2 |> mutate(cluster = as.character(cluster))
  setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

  antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1')
  antigen_times <- pmax(1, c(0, c(6, 10, 14)*7, 6*30.4, 9*30.4) - 7)

  for (i in 1:length(antigen)) {

    use_window_min <- ifelse(antigen[i] == 'bcg', 0, exposure_window)
    vax_data <- ldata[[antigen[i]]]; setDT(vax_data)

    # the exposure_window-day window bound
    vax_data[, `:=`(
      start_dt = due_date - use_window_min,
      end_dt   = due_date + exposure_window,
      cluster  = as.character(cluster)
    )]

    # find all rows in cds_geoloc where cluster matches AND date is between start/end
    # and sum the 'heatwave' column for @ child
    results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                          .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
                          by = .EACHI] |>
      setNames(c('cluster', 'start_dt', 'end_dt', 'heatwave')) |>
      distinct()

    # merging to have the heatwave and temperature variable
    vax_data <- merge(
      vax_data, results,
      by = c('cluster', 'start_dt', 'end_dt'),
      all.x = T
    ) |>
      mutate(heatwave = ifelse(heatwave == 0, 'absent', 'present'))

    # merging with covariates
    vax_data <- merge(vax_data |> mutate(across(c(wt, caseid, bidx), as.character))
                      , cdata |> mutate(wt = as.character(wt)),
                      by = c('caseid', 'bidx', 'admin', 'cluster', 'wt', 'residence'), all.x = T)

    vax_data <- vax_data |>
      mutate(
        wt = as.numeric(wt),
        delay = as.numeric(vaxx_date - due_date),
        delayclass = ifelse(delay > 28, 1, 0)
      )

    # omitting people vaxxed way before due dates (-7 is max before)
    # omitting people who have not met vaccine age requirements: 7 days before vaxx
    vax_data <- vax_data |> filter(delay >= -7 | is.na(delay))
    age_days <- as.numeric(vax_data$interview_date - vax_data$birth_date)
    age_id_omit <- which(age_days < (antigen_times[i]))
    {if (length(age_id_omit) == 1) {vax_data <- vax_data[-age_id_omit, ]}}
    min_age <- min(as.numeric(vax_data$interview_date - vax_data$birth_date))

    model_data <- vax_data |>
      data.frame() |>
      transmute(
        caseid, bidx,

        wt, v021, v022, geozone,

        event_time, outcome_event, censor_time,

        heatwave = factor(
          heatwave,
          levels = c("absent", "present"),
          labels = c("No heatwave", "Heatwave")
        ),

        residence = factor(
          residence,
          levels = c("urban", "rural"),
          labels = c("Urban", "Rural")
        ),

        birth_order,

        place_delivery = factor(
          place_delivery,
          levels = c("home", "institution"),
          labels = c("Home", "Health facility")
        ),

        child_gender = factor(
          child_gender,
          levels = c("male", "female"),
          labels = c("Male", "Female")
        ),

        time_to_hf = factor(
          time_to_hf,
          levels = c("<30 mins", "31-60 mins", "1-2 hrs", "2+ hrs"),
          labels = c(
            "<30 minutes",
            "31–60 minutes",
            "1–2 hours",
            ">2 hours"
          )
        ),

        wealth_index = factor(
          wealth_index,
          levels = c("poorest", "poorer", "middle", "richer", "richest"),
          labels = c("Poorest", "Poorer", "Middle", "Richer", "Richest")
        ),

        meduc = factor(
          meduc,
          levels = c("no education", "primary", "secondary", "higher"),
          labels = c(
            "No formal education",
            "Primary school",
            "Secondary school",
            "Higher education"
          )
        ),

        mother_age_group = case_when(
          mother_age_birth <= 19 ~ "<=19 yrs",
          mother_age_birth >= 20 & mother_age_birth <= 24 ~ "20-24 yrs",
          mother_age_birth >= 25 & mother_age_birth <= 29 ~ "25-29 yrs",
          mother_age_birth >= 30 & mother_age_birth <= 34 ~ "30-34 yrs",
          mother_age_birth >= 35 & mother_age_birth <= 39 ~ "35-39 yrs",
          mother_age_birth >= 40 & mother_age_birth <= 44 ~ "40-44 yrs",
          mother_age_birth >= 45 ~ "45+ yrs",
          TRUE ~ NA_character_
        ),
        mother_age_group = factor(
          mother_age_group,
          levels = c("<=19 yrs", "20-24 yrs","25-29 yrs","30-34 yrs",
                     "35-39 yrs","40-44 yrs","45+ yrs"),
          labels = c("≤19 years","20–24 years","25–29 years","30–34 years",
                     "35–39 years","40–44 years","45 years or older")
        ),

        geozone = factor(geozone,
                         levels = c("north west", "north east", "north central", "south east",
                                    "south south", "south west"),
                         labels = c("North West", "North East", "North Central", "South East",
                                    "South South", "South West"))
      )
    stopifnot(sum(is.na(model_data$heatindex)) == 0)

    {
      model_data$geozone           <- setLabel(model_data$geozone, "Geographic region")
      model_data$heatwave          <- setLabel(model_data$heatwave, "Heatwave")
      model_data$residence         <- setLabel(model_data$residence, "Place of residence")
      model_data$birth_order       <- setLabel(model_data$birth_order, "Birth order")
      model_data$place_delivery    <- setLabel(model_data$place_delivery, "Place of delivery")
      model_data$child_gender      <- setLabel(model_data$child_gender, "Child's sex")
      model_data$time_to_hf        <- setLabel(model_data$time_to_hf, "Travel time to health facility")
      model_data$wealth_index      <- setLabel(model_data$wealth_index, "Household wealth")
      model_data$meduc             <- setLabel(model_data$meduc, "Mother's education")
      model_data$mother_age_group  <- setLabel(model_data$mother_age_group, "Mother's age")
    }

    model_data <- model_data |>
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

    # a stop check if time_end < time_start: only managed to catch one case in i=3
    {if (i == 3) model_data <- model_data |> mutate(time_end = ifelse(time_start > time_end, time_start + 1, time_end))}
    # stopifnot(model_data |> filter(status == 3) |> mutate(check = time_end < time_start) |> pull(check) |> sum() == 0)
    # table(model_data$time_start[model_data$time_start < 0]);summary(model_data$time_end);summary(model_data$time_start)

    # BCG checks to modify time points for people receiving vax before b.day
    {if (antigen[i] == 'bcg') {
      model_data$time_start <- ifelse(model_data$time_start <= 0, 1, model_data$time_start)
      model_data$time_end <- ifelse(model_data$time_end <= 0, 1, model_data$time_end)
    }}

    # creating the full design
    mstd <- merge(master_data,
                  model_data |>
                    mutate(across(c(bidx, v021, v022), as.character),
                           selector = 'Include') |>
                    select(-wt),
                  by = c('caseid', 'bidx', 'v021', 'v022'), all.x = T) |>
      mutate(selector = ifelse(is.na(selector), 'Not included', selector))

    full_surv_design <- svydesign(
      id = ~v021,          # Primary Sampling Unit / Cluster
      strata = ~v022,      # Stratification variable/022
      weights = ~wt,     # DHS weight (remember to divide v005 by 1,000,000 first)
      data = mstd,
      nest = TRUE
    )
    surv_design <- subset(full_surv_design, selector == 'Include')
    options(survey.lonely.psu = 'adjust')

    fit_hw <- svysurvreg(
      Surv(time_start, time_end, status, type = "interval") ~
        heatwave + residence + birth_order + place_delivery + child_gender +
        time_to_hf + wealth_index + meduc + mother_age_group + geozone,
      dist = 'weibull',
      design = surv_design,
      data = mstd
    )

    # saveRDS(list(fit_hw = fit_hw, fit_hi = fit_hi),
    #         paste0('output/models/', antigen[i], '-file.rds'))

    confint(fit_hw)[2, 1]

    full_data <- bind_rows(
      full_data,
      data.frame(
        antigen = antigen[i],
        prop_exposed = as.numeric(table(model_data$heatwave)[2]/nrow(model_data)),
        aft_est = as.numeric(coef(fit_hw)[2]) |> exp(),
        aft_lwr = exp(confint(fit_hw)[2, 1]),
        aft_upr = exp(confint(fit_hw)[2, 2])
      ) |>
        crossing(grd[g, ])
    )



  }

}

saveRDS(full_data, '../data/processed/heatwave-sensitivity-analysis.rds')


# Plotting results --------------------------------------------------------

full_data <- readRDS('../data/processed/heatwave-sensitivity-analysis.rds')
{
  full_data <- readRDS('../data/processed/heatwave-sensitivity-analysis.rds') |>
    mutate(rolling = paste0('Rolling window (days): ', rolling),
           quantile = paste0('Return level (years): ', quantile),
           window = paste0('Exposure window (days): ', window))
  glimpse(full_data)

  true <- full_data |>
    mutate(rownum = row_number()) |>
    filter(rolling == 'Rolling window (days): 5' &
             quantile == 'Return level (years): 4' &
             window == 'Exposure window (days): 7')
  sens <- full_data #[-true$rownum, ]


  true <- tidyr::crossing(
    true |> select(-quantile, -rownum),
    quantile = unique(sens$quantile)
  ) |>
    mutate(antigen = factor(antigen, levels = c("bcg", "penta1", "penta2", "penta3", "vita1", "mcv1"), labels = c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')))

  order_exposure_window <- c('Exposure window (days): 3', 'Exposure window (days): 7', 'Exposure window (days): 14')
  order_return_level <- c('Return level (years): 4', 'Return level (years): 8', 'Return level (years): 12')
  order_rolling_window <- c('Rolling window (days): 3', 'Rolling window (days): 5' ,'Rolling window (days): 7')

  sens <- sens |>
    mutate(
      antigen = factor(antigen, levels = c("bcg", "penta1", "penta2", "penta3", "vita1", "mcv1"), labels = c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')),
      rolling = factor(rolling, levels = order_rolling_window, labels = order_rolling_window),
      quantile = factor(quantile, levels = order_return_level, labels = order_return_level),
      window = factor(window, levels = order_exposure_window, labels = order_exposure_window)
    )
}


# plot 1: using 4 year return level only
ggplot(sens |> #mutate(across(starts_with('aft'), \(x) ifelse(x > 7, NA, x))) |>
         filter(quantile == 'Return level (years): 4'),
       aes(x = aft_est, y = window, shape = rolling)) + # , linetype = rolling

  # Primary model CI
  geom_rect(
    data = true,
    aes(xmin = aft_lwr, xmax = aft_upr, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +

  # Null
  geom_vline(xintercept = 1, linetype = "dashed") +

  # Sensitivity estimates
  geom_errorbarh(
    aes(xmin = aft_lwr, xmax = aft_upr),
    position = position_dodge(width = 0.6),
    height = 0.15
  ) +

  geom_point(
    position = position_dodge(width = 0.6),
    size = 2.5
  ) +

  facet_wrap( ~ antigen, scales = 'free_x') +

  labs(
    x = "Event-time ratio (AFT coefficient)",
    y = NULL, #"Exposure window (days)",
    shape = "Rolling window (days)",
    linetype = "Rolling window (days)"
  ) +

  theme_bw() +
  theme(
    axis.line = element_blank(),
    axis.line.y = element_blank(),
    legend.position = 'bottom'
  )
ggsave('../output/img/sensitivity/aft-est.png', height = 7, width = 9, dpi = 300)

# table
d1 <- sens |> filter(quantile == 'Return level (years): 4') |>
  mutate(
    across(starts_with('aft'), \(x) round(x, 2)),
    coef = paste0(aft_est, '(', aft_lwr, ' - ', aft_upr, ')'),
    # prop_exposed = paste0(round(prop_exposed * 100, 2), '%')
  ) |>
  select(Antigen = antigen, `Proportion exposed` = prop_exposed,
         `Rolling window (days)` = rolling,
         `Return level (years)` = quantile,
         `Exposure window (days)` = window,
         `AFT Estimate (95% CI)` = coef)

# proportion exposed
ggplot(
  d1,
  aes(
    x = `Exposure window (days)`,
    y = Antigen,
    fill = `Proportion exposed`
  )
) +
  geom_tile() +
  geom_text(
    aes(label = paste0(round(`Proportion exposed` * 100, 2), "%")),
    size = 3.5,
    col = 'white'
  ) +
  facet_wrap(~ `Rolling window (days)`) +
  scale_fill_viridis_c(
    name = "Proportion exposed",
    limits = c(0.001273, 0.010723),
    breaks = c(0.001, 0.003, 0.005, 0.007, 0.009, 0.011),
    labels = percent_format(accuracy = 0.1),
    guide = guide_colorbar(
      barwidth = 12,
      barheight = 0.6
    )
  ) +
  scale_x_discrete(
    labels = \(x) gsub("Exposure window \\(days\\): ", "", x)
  ) +
  labs(
    x = "Exposure window (days)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    panel.grid = element_blank(),
    strip.background = element_rect()
  )
ggsave('../output/img/sensitivity/prop-exposed.png', height = 6, width = 12, dpi = 300)


# plot 2
ggplot(sens |> mutate(across(starts_with('aft'), \(x) ifelse(x > 7, NA, x))),
       aes(x = aft_est, y = window, shape = rolling)) + # , linetype = rolling

  # Primary model CI
  geom_rect(
    data = true,
    aes(xmin = aft_lwr, xmax = aft_upr, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +

  # Null
  geom_vline(xintercept = 1, linetype = "dashed") +

  # Sensitivity estimates
  geom_errorbarh(
    aes(xmin = aft_lwr, xmax = aft_upr),
    position = position_dodge(width = 0.6),
    height = 0.15
  ) +

  geom_point(
    position = position_dodge(width = 0.6),
    size = 2.5
  ) +

  facet_grid(antigen ~ quantile, scales = 'free') +

  labs(
    x = "Event-time ratio (AFT coefficient)",
    y = NULL, #"Exposure window (days)",
    shape = "Rolling window (days)",
    linetype = "Rolling window (days)"
  ) +

  theme_bw() +
  theme(
    axis.line = element_blank(),
    axis.line.y = element_blank(),
    legend.position = 'bottom'
  )


