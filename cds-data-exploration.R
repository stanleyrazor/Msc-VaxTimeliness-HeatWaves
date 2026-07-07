

# libs | funcs ------------------------------------------------------------

pacman::p_load(dplyr, tidyr, lubridate, terra, stringr, zoo, purrr, mev, tidyr,
               sf, ggplot2, gamlss, gamlss.dist, gamlssx, tidyterra, gganimate)

# plotting theme
theme_use <- theme_bw(base_family = 'Times New Roman') +
  theme(panel.grid = element_line(0))

# function for missingness analysis
mvs <- naniar::miss_var_summary

# Notes -------------------------------------------------------------------

# 2020-12-06 is earliest birth day possible in DHS. Thus, model fitting should only use data before this period.
# use 1940:2019 for model fitting


# Data --------------------------------------------------------------------

# the cluster geo codes
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST, residence = URBAN_RURA, geometry)

g2 <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> st_as_sf()

# Pixel-level analysis ----------------------------------------------------

use_adm <- g2[which.max(st_area(g2)),]
use_clust <-  st_filter(g1, use_adm)

bf_gps <- sf::st_buffer(use_clust, dist = ifelse(use_clust$residence == 'U', 5/111, 5/111))


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

  # CROP AND MASK HERE:
  sel_t1 <- crop(t3, use_adm, touches = T, snap = 'out', mask = T)

  use_clust_buffer <- use_clust %>%
    st_transform(32631) %>%
    st_buffer(5e3) %>% # 5000 meters = 5km
    st_transform(4326)

  # 2. Build the Plot
  ggplot() +
    # Plot the Raster tiles
    geom_spatraster(data = sel_t1[[80]]) +

    # Add the Admin boundary (no fill so we see the raster)
    geom_sf(data = use_adm, fill = NA, color = "white", size = 0.8) +

    # Add the 5km buffers
    geom_sf(data = use_clust_buffer, fill = "red", alpha = 0.2, color = "red") +

    # Add the cluster points
    geom_sf(data = use_clust, col = 'red') +

    # Styling
    scale_fill_viridis_c(name = "Temp (°C)", na.value = NA) +
    theme_void() +
    labs(
      title = "2m Temperature in Baruten (year 2019)",
      subtitle = "5km radius buffers shown around clusters\nBaruten is the largest Admin 2 region by area",
      x = "Longitude",
      y = "Latitude"
    )
  ggsave("output/img/baruten-img.png", width = 8, height = 8, dpi = 2000)

  # extracting to see if i get what I think
  terra::extract(sel_t1[[2]], vect(bf_gps), ID = FALSE, method = 'simple') |>
    as.data.frame() |>
    mutate(cluster = bf_gps$cluster)

  terra::extract(
    sel_t1[[2]],
    vect(use_clust_buffer),
    fun = mean,
    exact = TRUE,
    na.rm = TRUE
  ) |>
    as.data.frame() |>
    mutate(cluster = use_clust_buffer$cluster)


  # Animation ---------------------------------------------------------------

  df <- as.data.frame(sel_t1, xy = TRUE) |>
    pivot_longer(-c(x, y), names_to = 'year', values_to = 'maxtemp') |>
    mutate(year = as.numeric(gsub("y_", "", year)))

  p <- ggplot() +
    geom_raster(
      data = df,
      aes(x = x, y = y, fill = maxtemp)
    ) +

    geom_sf(
      data = use_adm,
      fill = NA,
      color = "white",
      linewidth = 0.8
    ) +

    geom_sf(
      data = use_clust_buffer,
      fill = "red",
      alpha = 0.2,
      color = "red",
      linewidth = 0.4
    ) +

    geom_sf(
      data = use_clust,
      color = "red",
      size = 2
    ) +

    scale_fill_viridis_c(
      name = "Temp (°C)",
      na.value = NA
    ) +

    coord_sf(expand = FALSE) +

    theme_void(base_size = 18) +

    theme(
      plot.title = element_text(
        size = 28,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 14,
        hjust = 0.5
      ),
      plot.caption = element_text(
        size = 11
      ),
      legend.position = "right"
    ) +

    labs(
      title = "{closest_state}",
      caption = "2m temp surface in Baruten LGA. In red are 5km DHS cluster buffers."
    ) +

    transition_states(
      year,
      transition_length = 1,
      state_length = 1
    ) +

    ease_aes("linear")

  animate(
    p,
    nframes = 160,
    fps = 4,
    width = 2400,
    height = 1800,
    res = 300,
    renderer = gifski_renderer("output/img/baruten_temperature_animation.gif")
  )




  t4 <- as.data.frame(t3, xy = T, wide = FALSE)
  head(t4)

  rm(t1, t2, t3); gc()
}
