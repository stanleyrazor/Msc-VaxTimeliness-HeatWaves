

# libs | funcs ------------------------------------------------------------

pacman::p_load(dplyr, tidyr, lubridate, terra, stringr, zoo, purrr, mev,
               gamlss, gamlss.dist, gamlssx)

# plotting theme
theme_use <- theme_bw(base_family = 'Times New Roman') +
  theme(panel.grid = element_line(0))

mvs <- naniar::miss_var_summary
# CDF of GEV
# pgev <- \(y, xi, tau, eta, sigma = 1) {
#   exp(-(1 + xi * sqrt(tau * sigma) * (y - eta))^(-1 / xi))
# }

# Notes -------------------------------------------------------------------

# numpixel*daily for one year gives: 872,910.
# for 85 years it is 74,197,350 rows (not possible to operate this on my machine)
# for simplification we compute tx5x on annual - ignore chaining at year ends

# 2020-12-06 is earliest birth day possible in DHS. Thus, model fitting should only use data before this period.
# use 1940:2019 for model fitting


# Data --------------------------------------------------------------------

# only test years: 2020-2024
ls_test <- c("m2_temp_2020.nc", "m2_temp_2021.nc", "m2_temp_2022.nc", "m2_temp_2023.nc", "m2_temp_2024.nc")
d0 <- map(ls_test, \(x) {

  yr <- str_split(x, '_', n = 3) |> lapply(FUN = last) |> unlist() |> str_split('[.]') |> lapply(FUN = first) |> unlist()

  # load and process
  t1 <- rast(paste0('data/cds/', x))
  t1 <- t1 - 273.15

  # Using terra::roll for 5-day moving window averages
  t2 <- roll(t1, n = 5, fun = mean, type = 'to', na.rm = T)

  as.data.frame(t2, xy = T, wide = F) |>
    mutate(date = ymd(paste0(yr, '-01-01')) + as.numeric(str_remove_all(layer, 't2m_valid_time='))) |>
    select(x, y, date, values)
})
test <- bind_rows(d0)

# all years max values
ls_full <- list.files('data/cds/', pattern = "*.nc")
d1 <- map(ls_full, \(x) {

  yr <- str_split(x, '_', n = 3) |> lapply(FUN = last) |> unlist() |> str_split('[.]') |> lapply(FUN = first) |> unlist()

  # load and process
  t1 <- rast(paste0('data/cds/', x))
  t1 <- t1 - 273.15

  # Using terra::roll for 5-day moving window averages
  t2 <- roll(t1, n = 5, fun = mean, type = 'to', na.rm = T)

  t3 <- as.data.frame(t2, xy = T, wide = F) |>
    mutate(date = ymd(paste0(yr, '-01-01')) + as.numeric(str_remove_all(layer, 't2m_valid_time='))) |>
    select(x, y, date, values)

  tmax <- t3 |>
    group_by(x, y) |>
    reframe(max = max(values)) |>
    mutate(year = yr)

  tmax
})

d2 <- bind_rows(d1)

# plotting max for sample pixel: 2.5 | 3.5
d2 |> filter(x == 2.5 & y == 3.5) |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max)) +
  theme_use

# plotting all on a plot:
d2 |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max,
                group = interaction(x, y)),
            lwd = 0.05, col = 'black') +
  labs(title = 'Trend of annual maximum Tx5x across 2,385 pixels_full',
       x = 'Year', y = 'Maximum Tx5x') +
  theme_use



# fitting Stationary GEV per location -------------------------------------

# data to use for model fitting:
d3 <- d2 |> filter(!year %in% as.character(2020:2024))

# x <- d3 |> filter(x == 2.5 & y == 3.5)
# f1 <- fit.gev(x$max); f1

stationary_est <- d3 |>
  group_by(x, y) |>
  reframe(
    loc = fit.gev(max)$estimate[1],
    scale = fit.gev(max)$estimate[2],
    shape = fit.gev(max)$estimate[3]
  )

# plot to understand:
# plot(seq(from = 25, to = 40, length = 1e3),
#      dgev(seq(from = 25, to = 40, length = 1e3),
#           loc = 33.80914, scale = 1.1128443, shape = -0.4762559),
#      type = 'l')
# abline(v = 36.34209, col = 'red')

# merging the stationary fit with the test data
s1 <- merge(test, stationary_est, by = c('x', 'y'), all.x = T)

# computing Pr(Y <= y) and return period. Also heatwave = p > 0.95
s2 <- s1 |>
  mutate(
    p = mev::pgev(q = values, loc = loc, scale = scale, shape = shape, lower.tail = T),
    return = 1 / (1 - p),
    heatwave = p >= .95
  )

# plotting date with highest number of heatwaves
ggplot(s2 |> filter(date == '2024-03-31'),
       aes(x = x, y = y)) +
  geom_tile(aes(fill = heatwave)) +
  theme_use

# Fitting non-stationary GEV ----------------------------------------------

d3 <- d2 |> filter(!year %in% as.character(2020:2024)) |>
  mutate(year = as.numeric(year) - 1940)

# x <- d3 |> filter(x == 2.5 & y == 3.5)
# f1 <- fit.gev(x$max); f1
# fitGEV(max ~ year, data = x)

split_data <- d3 |> group_by(x, y) |> group_split()
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
s1 <- merge(test |> mutate(year = year(date) - 1940),
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


# DHS Geo-code for first stage triangulations -----------------------------

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)

# displacement is up to 5km for rural points (with 1% of rural cluster being
# randomly displaced up to 10km) and up to 2km for urban
# (all urban points have the same 2km randomization applied)


# only test years: 2020-2024
ls_test <- c("m2_temp_2020.nc", "m2_temp_2021.nc", "m2_temp_2022.nc", "m2_temp_2023.nc", "m2_temp_2024.nc")
d0 <- map(ls_test, \(x) {

  yr <- str_split(x, '_', n = 3) |> lapply(FUN = last) |> unlist() |> str_split('[.]') |> lapply(FUN = first) |> unlist()

  # load and process
  t1 <- rast(paste0('data/cds/', x))
  t1 <- t1 - 273.15

  # Using terra::roll for 5-day moving window averages
  t2 <- roll(t1, n = 5, fun = mean, type = 'to', na.rm = T)

  # 2km urban, 5km rural. since temperature raster is so coarse (0.25 degrees ~ 27.75 km),
  # we simply extract the the pixel value - method simple.
  t3 <- terra::extract(
    t2, bf_gps,          # All clusters included
    method = "simple", weights = FALSE,
    na.rm = TRUE,    # Ignore missing values
    bind = TRUE      # Recombine with original `bf_gps` data
  )

  # 1380 * 365
  t4 <- as.data.frame(t3) |> pivot_longer(-c(cluster, residence)) |>
    mutate(date = ymd(paste0(yr, '-01-01')) + as.numeric(str_remove_all(name, 't2m_valid_time.'))) |>
    select(cluster, residence, date, value)
  t4
})
test <- bind_rows(d0)


# all years max values
ls_full <- list.files('data/cds/', pattern = "*.nc")
bf_gps <- sf::st_buffer(g1, dist = ifelse(g1$residence == 'U', 2/111, 5/111))

# extracting max per cluster by year
d1 <- map(ls_full, \(x) {

  yr <- str_split(x, '_', n = 3) |> lapply(FUN = last) |> unlist() |> str_split('[.]') |> lapply(FUN = first) |> unlist()

  # load and process - convert from kelvin to celsius
  t1 <- rast(paste0('data/cds/', x))
  t1 <- t1 - 273.15

  # Using terra::roll for 5-day moving window averages
  t2 <- roll(t1, n = 5, fun = mean, type = 'to', na.rm = T)

  # 2km urban, 5km rural. since temperature raster is so coarse (0.25 degrees ~ 27.75 km),
  # we simply extract the the pixel value - method simple.
  t3 <- terra::extract(
    t2, bf_gps,          # All clusters included
    method = "simple", weights = FALSE,
    na.rm = TRUE,    # Ignore missing values
    bind = TRUE      # Recombine with original `bf_gps` data
  )

  # 1380 * 365
  t4 <- as.data.frame(t3) |> pivot_longer(-c(cluster, residence)) |>
    mutate(date = ymd(paste0(yr, '-01-01')) + as.numeric(str_remove_all(name, 't2m_valid_time.'))) |>
    select(cluster, residence, date, value)

  tmax <- t4 |>
    group_by(cluster, residence) |>
    reframe(max = max(value)) |>
    mutate(year = yr)

  tmax
})

d2 <- bind_rows(d1)

# plotting max for sample cluster: 1000
d2 |> filter(cluster == 1000) |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max)) +
  theme_use

# plotting all clusters on a plot:
d2 |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max,
                group = interaction(cluster, residence)),
            lwd = 0.05, col = 'black') +
  labs(title = 'Trend of annual maximum Tx5x across 1,380 clusters covered',
       x = 'Year', y = 'Annual Maximum Tx5x') +
  theme_use


# Time-series GEV on clusters ---------------------------------------------

d3 <- d2 |> filter(!year %in% as.character(2020:2024)) |>
  mutate(year = as.numeric(year) - 1940)

split_data <- d3 |> group_by(cluster, residence) |> group_split()
for (i in 1:length(split_data)) {

  grp_data <- split_data[[i]]
  tmp_model <- fitGEV(max ~ year, data = grp_data, n.cyc = 500)
  tmp_data <- data.frame(
    cluster = unique(grp_data$cluster), residence = unique(grp_data$residence),
    icpt  = as.numeric(tmp_model$mu.coefficients[1]),
    time  = as.numeric(tmp_model$mu.coefficients[2]),
    scale = exp(as.numeric(tmp_model$sigma.coefficients)),
    shape = as.numeric(tmp_model$nu.coefficients)
  )

  if (i == 1) nonstationary_est <- tmp_data
  else nonstationary_est <- bind_rows(nonstationary_est, tmp_data)

}

# merging the stationary fit with the test data
s1 <- merge(test |> mutate(year = year(date) - 1940),
            nonstationary_est, by = c('cluster', 'residence'), all.x = T)

# computing Pr(Y <= y) and return period. Also heatwave = p > 0.95
s2 <- s1 |>
  mutate(
    # the eta/location is b0 + b1*year
    loc = icpt + (time * year),

    p = mev::pgev(q = value, loc = loc, scale = scale, shape = shape, lower.tail = T),
    return = 1 / (1 - p),
    heatwave = p >= .95
  ) |>
  select(cluster, year, date, value, icpt, coeftime = time, loc, scale, shape, p, return, heatwave)

arrow::write_parquet(x = s2, sink = 'data/processed/cluster-processed.parquet')



