

# libs | funcs ------------------------------------------------------------

pacman::p_load(dplyr, tidyr, lubridate, terra, stringr, zoo, purrr, mev,
               gamlss, gamlss.dist, gamlssx)

# plotting theme
theme_use <- theme_bw(base_family = 'Times New Roman') +
  theme(panel.grid = element_line(0))

# function for missingness analysis
mvs <- naniar::miss_var_summary

# Notes -------------------------------------------------------------------

# numpixel*daily for one year gives: 872,910.
# for 85 years it is 74,197,350 rows (not possible to operate this on my machine)
# for simplification we compute tx5x on annual - ignore chaining at year ends

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


# fitting Stationary GEV per pixel -------------------------------------

# x <- t4 |> filter(x == 2.5 & y == 3.5)
# f1 <- fit.gev(x$max); f1

stationary_est <- t4 |>
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

# Fitting non-stationary GEV per pixel -------------------------------------

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


# DHS Geo-code for first stage triangulation ------------------------------

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)

bf_gps <- sf::st_buffer(g1, dist = ifelse(g1$residence == 'U', 2/111, 5/111))

# displacement is up to 5km for rural points (with 1% of rural cluster being
# randomly displaced up to 10km) and up to 2km for urban
# (all urban points have the same 2km randomization applied)

# only test years: 2020-2024
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

  # 1380 * 365
  t3 <- terra::extract(t2, vect(bf_gps), ID = FALSE, method = 'simple') |>
    as.data.frame() |>
    mutate(cluster = bf_gps$cluster) |>
    pivot_longer(-cluster, names_to = "layer", values_to = "tx5x") |>
    group_by(cluster) |>
    mutate(date = all_dates) |> # confirmed
    ungroup()

  # 3. Final Join with cluster info
  test <- bf_gps |> st_drop_geometry() |> left_join(t3, by = "cluster") |>
    select(cluster, residence, time = date, tx5x)

  # clear memory
  rm(t1, t2, t3); gc()
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

  # 1380 * 365
  train <- terra::extract(t3, vect(bf_gps), ID = FALSE, bind = TRUE, method = 'simple') |>
    as.data.frame() |>
    pivot_longer(-c(cluster, residence), names_to = "year", values_to = "max") |>
    mutate(year = str_remove_all(year, 'y_'))

  # clear memory
  rm(t1, t2, t3); gc()
}

# plotting max for sample cluster: 1000
train |> filter(cluster == 1000) |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max)) +
  theme_use

# plotting all clusters on a plot:
train |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max,
                group = interaction(cluster, residence)),
            lwd = 0.05, col = 'black') +
  labs(title = 'Trend of annual maximum Tx5x across 1,380 clusters covered',
       x = 'Year', y = 'Annual Maximum Tx5x') +
  theme_use

# Time-series GEV on clusters ---------------------------------------------

d3 <- train |> mutate(year = as.numeric(year) - 2019)

split_data <- d3 |> group_by(cluster, residence) |> group_split()
for (i in 1:length(split_data)) {

  grp_data <- split_data[[i]]
  tmp_model <- fitGEV(max ~ year, data = grp_data, n.cyc = 500)
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
s1 <- merge(test |> mutate(year = year(time) - 2019),
            nonstationary_est, by = c('cluster', 'residence'), all.x = T)

# computing Pr(Y <= y) and return period. Also heatwave = p > 0.95
s2 <- s1 |>
  mutate(
    # the eta/location is b0 + b1*year
    loc = icpt + (coeftime * year),

    p = mev::pgev(q = tx5x, loc = loc, scale = scale, shape = shape, lower.tail = T),
    return = 1 / (1 - p),
    heatwave = p >= .95
  ) |>
  select(cluster, year, date = time, tx5x, icpt, coeftime, loc, scale, shape, p, return, heatwave)

arrow::write_parquet(x = s2, sink = 'data/processed/cluster-processed.parquet')


# Admin-2 processing ------------------------------------------------------

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)

g2 <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> st_as_sf()

# displacement is up to 5km for rural points (with 1% of rural cluster being
# randomly displaced up to 10km) and up to 2km for urban
# (all urban points have the same 2km randomization applied)

# only test years: 2020-2024
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

  # 1380 * 365
  t3 <- terra::extract(t2, vect(g2 |> select(GID_2, geometry)), ID = FALSE, fun = 'mean', bind = TRUE) |>
    as.data.frame() |>
    pivot_longer(-GID_2, names_to = "layer", values_to = "tx5x") |>
    group_by(GID_2) |>
    mutate(date = all_dates) |> # confirmed
    ungroup()

  # 3. Final Join with cluster info
  test <- g2 |> select(GID_2, NAME_2) |>
    st_drop_geometry() |>
    left_join(t3, by = "GID_2") |>
    select(GID_2, NAME_2, time = date, tx5x)

  # clear memory
  rm(t1, t2, t3); gc()
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

  # 1380 * 365
  train <- terra::extract(t3, vect(g2), ID = TRUE, bind = FALSE, fun = 'mean') |>
    as.data.frame() |>
    pivot_longer(-ID, names_to = "year", values_to = "max") |>
    mutate(
      year = as.numeric(str_remove_all(year, 'y_')),
      GID_2 = g2$GID_2[ID]
    )

  # clear memory
  rm(t1, t2, t3); gc()
}

# plotting max for sample cluster: 1000
train |> filter(GID_2 == 'NGA.1.1_1') |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max)) +
  theme_use

# plotting all clusters on a plot:
train |>
  ggplot() +
  geom_line(aes(x = as.numeric(year), y = max,
                group = GID_2),
            lwd = 0.05, col = 'black') +
  labs(title = 'Trend of annual maximum Tx5x across 775 Admin-2 regions',
       x = 'Year', y = 'Annual Maximum Tx5x') +
  theme_use

# Time-series GEV on clusters ---------------------------------------------

d3 <- train |> mutate(year = as.numeric(year) - 2019)

split_data <- d3 |> group_by(GID_2) |> group_split()
for (i in 1:length(split_data)) {

  grp_data <- split_data[[i]]
  tmp_model <- fitGEV(max ~ year, data = grp_data, n.cyc = 500)
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
s1 <- merge(test |> mutate(year = year(time) - 2019),
            nonstationary_est, by = c('GID_2'), all.x = T)

# computing Pr(Y <= y) and return period. Also heatwave = p > 0.95
s2 <- s1 |>
  mutate(
    # the eta/location is b0 + b1*year
    loc = icpt + (coeftime * year),

    p = mev::pgev(q = tx5x, loc = loc, scale = scale, shape = shape, lower.tail = T),
    return = 1 / (1 - p),
    heatwave = p >= .95
  ) |>
  select(GID_2, year, date = time, tx5x, icpt, coeftime, loc, scale, shape, p, return, heatwave)

arrow::write_parquet(x = s2, sink = 'data/processed/admin2-processed.parquet')



