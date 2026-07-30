


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
    mutate(date = ymd(substr(layer, 4, 13)),
           year = year(date)) |>
    select(cluster, residence, year, date, dewpoint = temp)

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
    mutate(date = ymd(substr(layer, 4, 13)),
           year = year(date)) |>
    select(cluster, residence, year, date, airtemp = temp)
}

# merging
head(d2); head(t2)
hi_data <- merge(t2, d2, by = c('cluster', 'residence', 'year', 'date')) |>
  mutate(
    heatindex = heat.index(t = airtemp, dp = dewpoint,
                           temperature.metric = 'celsius', output.metric = 'celsius',
                           round = Inf)
  )

arrow::write_parquet(x = hi_data, sink = 'data/processed/heatindex-processed.parquet')


# Plotting ----------------------------------------------------------------

# obtain the cluster lon-lat, and rank by lat (northern-southern)
plt_data <- hi_data |>
  mutate(month = month.abb[month(date)],
         month = factor(month, levels = month.abb, labels = month.abb)) |>
  group_by(year, month, cluster, residence) |>
  reframe(heatindex = mean(heatindex)) |>

  merge(
    data.frame(st_centroid(g1) |>
                 st_coordinates()) |>
      setNames(c('lon', 'lat')) |>
      bind_cols(cluster = g1$cluster),
    by = 'cluster',
    all.x = T
  )

yrs <- sort(unique(plt_data$year))
for(i in seq_along(yrs)){

  tmp <- merge(g1, filter(plt_data, year == yrs[i]), by = "cluster")

  p <- ggplot() +
    geom_sf(data = shp, colour = "grey70", fill = NA, linewidth = 0.1) +
    geom_sf(data = tmp, aes(colour = heatindex)) +
    facet_wrap(~month, nrow = 3) +
    scale_color_viridis_b(
      option = "inferno", # "inferno" or "magma" are perfect for heat data
      direction = 1,      # 1 puts darker/cooler colors at 15 and bright yellow at 39
      breaks = seq(15, 39, length.out = 10),
      guide = guide_coloursteps(
        show.limits = TRUE,
        even.steps = TRUE
      ),
      name = expression(Heat-Index~"("*degree*C*")")
    ) +
    labs(title = yrs[i]) +
    theme_void(base_family = "Times New Roman") +
    theme(
      legend.position = "right",

      legend.key.height = unit(15, "mm"),
      legend.title = element_text(margin = margin(b = 10)),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ## optionally save
  ggsave(
    paste0("output/img/heatmaps/heat index/", yrs[i], ".png"),
    p,
    width = 12,
    height = 9,
    dpi = 1000
  )
}

# Cleanup -----------------------------------------------------------------


rm(hi_data, t1, t2, d1, d2, mvs, g1, bf_gps, full_dates, ls_full, plt_data)

