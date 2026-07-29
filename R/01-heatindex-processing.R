


# libs | funcs ------------------------------------------------------------

pacman::p_load(dplyr, tidyr, lubridate, terra, stringr, zoo, purrr, mev,
               gamlss, gamlss.dist, gamlssx, ggplot2, sf, weathermetrics)

# function for missingness analysis
mvs <- naniar::miss_var_summary
formals(table)$useNA = 'always'


# Data --------------------------------------------------------------------

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)
bf_gps <- g1 |> st_transform(3857) |> st_buffer(dist = 10000) |> st_transform(4326)

full_dates <- seq(from=ymd('2020-01-01'), to = ymd('2024-12-31'), by = 1)



# Analysis ----------------------------------------------------------------

# dewpoint temp
{
  ls_full <- list.files('data/cds/2m_dewpoint_temperature/', pattern = "*.nc") |> sort()
  ls_full <- paste0("data/cds/2m_dewpoint_temperature/", ls_full)

  # full-stack 1940-2019
  d1 <- rast(ls_full)
  d1 <- d1 - 273.15

  # computation of cluster TX5X index
  d2 <- terra::extract(d1, vect(bf_gps), fun = mean, na.rm = TRUE, ID = FALSE) |>
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
          select(cluster, residence, year, date = layer, dewpoint = tx5x)
      }
    ) |>
    bind_rows()

}

# 2m air temp
{
  ls_full <- list.files('data/cds/', pattern = "*.nc") |> sort()
  ls_full <- paste0("data/cds/", ls_full) |> tail(5)

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
          select(cluster, residence, year, date = layer, airtemp = tx5x)
      }
    ) |>
    bind_rows()

}

# merging
head(d2); head(t2)
hi_data <- merge(t2, d2, by = c('cluster', 'residence', 'year', 'date')) |>
  mutate(
    heatindex = heat.index(t = airtemp, dp = dewpoint, temperature.metric = 'celsius', output.metric = 'celsius')
  )

arrow::write_parquet(x = hi_data, sink = 'data/processed/heatindex-processed.parquet')
rm(hi_data, t1, t2, d1, d2)

