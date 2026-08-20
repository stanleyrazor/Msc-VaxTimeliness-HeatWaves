

# libs | funcs ------------------------------------------------------------

pacman::p_load(dplyr, tidyr, lubridate, terra, stringr, zoo, purrr, mev,
               gamlss, gamlss.dist, gamlssx, ggplot2, sf)

# plotting theme
theme_use <- theme_bw(base_family = 'Times New Roman') +
  theme(panel.grid = element_line(0))

# function for missingness analysis
mvs <- naniar::miss_var_summary

# Notes -------------------------------------------------------------------

# GMST obtained from first facet of: https://www.jkclimate.fr/Dashboard2025/global_mean_temperature.html
# 2020-12-06 is earliest birth day possible in DHS. Thus, model fitting should only use data before this period.
# use 1940:2019 for model fitting

# DHS Geo-code level analysis  --------------------------------------

# the GMST data
gmst <- read.csv('data/GMST/tas_ERA5.csv')[-c(1:30), 2:3] |>
  setNames(c('year', 'gmst')) |>
  filter(year %in% 1940:2025) |>
  mutate(across(everything(), as.numeric))

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)
bf_gps <- g1 |> st_transform(3857) |> st_buffer(dist = 10000) |> st_transform(4326)

full_dates <- seq(from=ymd('1940-01-01'), to = ymd('2024-12-31'), by = 1)

{
  ls_full <- list.files('data/cds/', pattern = "*.nc") |> sort()
  ls_full <- paste0("data/cds/", ls_full)

  # full-stack 1940-2019
  t1 <- rast(ls_full)
  t1 <- t1 - 273.15

  # computation of cluster TX5X index
  t2 <- terra::extract(t1, vect(bf_gps), fun = mean, na.rm = TRUE, ID = FALSE) |>
    setNames(paste0('dt.', full_dates)) |>
    mutate(
      cluster = bf_gps |> pull(cluster),
      residence = bf_gps |> pull(residence)
    ) |>
    pivot_longer(-c(cluster, residence), names_to = "layer", values_to = "temp") |>
    mutate(layer = ymd(substr(layer, 4, 13))) |>
    group_by(cluster) |>
    group_split() |>
    lapply(
      FUN = \(dt) {
        dt |>
          arrange(layer) |>
          mutate(
            year = year(layer),
            tx5x = slider::slide_dbl(temp, mean, .before = 4, .complete = TRUE),
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
  }

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
      # heatwave = p >= .90

      # 1-in-10 year return level threshold (90th percentile of annual max distribution)
      thresh_10yr = mev::qgev(p = rep(0.90, nrow(s1)), loc = loc, scale = scale, shape = shape),
      heatwave = tx5x >= thresh_10yr
    ) |>
    select(cluster, year, date, tx5x, icpt, coeftime, loc, scale, shape, p, return, thresh_10yr, heatwave)

  arrow::write_parquet(x = s2, sink = 'data/processed/cluster-processed.parquet')
}


# Admin-2 level analysis ------------------------------------------------------

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)

g2 <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> st_as_sf()

{
  ls_full <- list.files('data/cds/', pattern = "*.nc") |> sort()
  ls_full <- paste0("data/cds/", ls_full)

  # full-stack 1940-2019
  t1 <- rast(ls_full)
  t1 <- t1 - 273.15

  # computation of cluster TX5X index
  t2 <- terra::extract(t1, vect(g2 |> select(GID_2, geometry)), ID = FALSE, fun = 'mean', bind = TRUE) |>
    as.data.frame() |>
    setNames(c('GID_2', paste0('dt.', full_dates))) |>
    pivot_longer(-GID_2, names_to = "layer", values_to = "temp") |>
    mutate(layer = ymd(substr(layer, 4, 13))) |>
    group_by(GID_2) |>
    group_split() |>
    lapply(
      FUN = \(dt) {
        dt |>
          arrange(layer) |>
          mutate(
            year = year(layer),
            tx5x = slider::slide_dbl(temp, mean, .before = 4, .complete = TRUE),
            tx5x = na.locf(tx5x, fromLast = T)
          ) |>
          select(GID_2, year, date = layer, tx5x)
      }
    ) |>
    bind_rows()

  # construction of the training and testing set
  train <- t2 |>
    filter(year %in% 1940:2019) |>
    group_by(GID_2, year) |>
    reframe(mx = max(tx5x))
  test <- t2 |> filter(year %in% 2020:2024)
  rm(t2)

  d3 <- train |> merge(gmst, by = 'year', all.x = T)
  split_data <- d3 |> group_by(GID_2) |> group_split()

  for (i in 1:length(split_data)) {

    grp_data <- split_data[[i]]
    tmp_model <- fitGEV(mx ~ gmst, data = grp_data, n.cyc = 500)
    tmp_data <- data.frame(
      GID_2 = unique(grp_data$GID_2),
      icpt  = as.numeric(tmp_model$mu.coefficients[1]),
      coeftime  = as.numeric(tmp_model$mu.coefficients[2]),
      scale = exp(as.numeric(tmp_model$sigma.coefficients)),
      shape = as.numeric(tmp_model$nu.coefficients)
    )

    if (i == 1) nonstationary_est <- tmp_data
    else nonstationary_est <- bind_rows(nonstationary_est, tmp_data)
  }

  # merging the model fit with the test data
  s1 <- merge(
    merge(test, gmst, by = 'year', all.x = T),
    nonstationary_est, by = c('GID_2'), all.x = T)

  # computing Pr(Y <= y) and return period. Also heatwave = p > 0.90
  s2 <- s1 |>
    mutate(
      # the eta/location is b0 + b1*gmst
      loc = icpt + (coeftime * gmst),

      p = mev::pgev(q = tx5x, loc = loc, scale = scale, shape = shape, lower.tail = T),
      return = 1 / (1 - p),
      # heatwave = p >= .90

      # 1-in-10 year return level threshold (90th percentile of annual max distribution)
      thresh_10yr = mev::qgev(p = rep(0.90, nrow(s1)), loc = loc, scale = scale, shape = shape),
      heatwave = tx5x >= thresh_10yr
    ) |>
    select(GID_2, year, date, tx5x, icpt, coeftime, loc, scale, shape, p, return, thresh_10yr, heatwave)

  arrow::write_parquet(x = s2, sink = 'data/processed/admin2-processed.parquet')
}

