

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

# Pixel-level analysis ----------------------------------------------------

# Test: only test years: 2020-2024
{
  ls_test <- c("m2_temp_2020.nc", "m2_temp_2021.nc", "m2_temp_2022.nc", "m2_temp_2023.nc", "m2_temp_2024.nc")
  ls_test <- paste0('data/cds/', ls_test)

  # full-test-stack 2020-2024
  t1 <- rast(ls_test)
  t1 <- t1 - 273.15

  # Tx5x creation
  t2 <- roll(t1, n = 5, fun = mean, type = 'to', na.rm = TRUE)
  all_dates <- as.Date("2020-01-01") + 0:(nlyr(t2) - 1)
  terra::time(t2) <- all_dates
  test <- as.data.frame(t2, xy = TRUE, wide = FALSE, time = TRUE) |>
    select(x, y, time, values)

  # clear memory
  rm(t1, t2); gc()
}

# Train: all years 1940 to 2019 values
{
  ls_full <- list.files('data/cds/', pattern = "*.nc") |> sort()
  ls_full <- paste0("data/cds/", ls_full)[-(81:85)]

  # full-stack 1940-2019
  t1 <- rast(ls_full)
  t1 <- t1 - 273.15

  # Tx5x creation
  t2 <- roll(t1, n = 5, fun = mean, type = 'to', na.rm = TRUE)
  all_dates <- as.Date("1940-01-01") + 0:(nlyr(t2) - 1)
  terra::time(t2) <- all_dates

  # annual maxima
  t3 <- tapp(t2, index = "years", fun = max)
  t4 <- as.data.frame(t3, xy = T, wide = FALSE)
  head(t4)

  rm(t1, t2, t3); gc()
}

# sample plotting max for sample pixel: 2.5 | 3.5
t4 |> filter(x == 2.5 & y == 3.5) |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max)) +
  theme_use

# plotting all on a plot:
t4 |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max,
                group = interaction(x, y)),
            lwd = 0.05, col = 'black') +
  labs(title = 'Trend of annual maximum Tx5x across 2,385 pixels_full',
       x = 'Year', y = 'Maximum Tx5x') +
  theme_use


# Time series GEV on pixels -------------------------------------

t5 <- t4 |> mutate(year = as.numeric(year) - 2019)

# x <- t5 |> filter(x == 2.5 & y == 3.5)
# f1 <- fit.gev(x$max); f1
# fitGEV(max ~ year, data = x)

split_data <- t5 |> group_by(x, y) |> group_split()
for (i in 1:length(split_data)) {

  grp_data <- split_data[[i]]
  tmp_model <- fitGEV(max ~ year, data = grp_data, n.cyc = 500)
  tmp_data <- data.frame(
    x = unique(grp_data$x), y = unique(grp_data$y),
    icpt  = as.numeric(tmp_model$mu.coefficients[1]),
    time  = as.numeric(tmp_model$mu.coefficients[2]),
    scale = exp(as.numeric(tmp_model$sigma.coefficients)),
    shape = as.numeric(tmp_model$nu.coefficients)
  )

  if (i == 1) nonstationary_est <- tmp_data
  else nonstationary_est <- bind_rows(nonstationary_est, tmp_data)

}

# merging the stationary fit with the test data
s1 <- merge(test |> mutate(year = year(date) - 2019),
            nonstationary_est, by = c('x', 'y'), all.x = T)

# computing Pr(Y <= y) and return period. Also heatwave = p > 0.95
s2 <- s1 |>
  mutate(
    # the eta/location is b0 + b1*year
    loc = icpt + (time * year),

    p = mev::pgev(q = values, loc = loc, scale = scale, shape = shape, lower.tail = T),
    return = 1 / (1 - p),
    heatwave = p >= .95
  ) |>
  select(x, y, year, date, values, icpt, coeftime = time, loc, scale, shape, p, return, heatwave)

# plotting date with highest number of heatwaves
ggplot(s2 |> filter(date == '2024-03-30'),
       aes(x = x, y = y)) +
  geom_tile(aes(fill = heatwave)) +
  theme_use



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

