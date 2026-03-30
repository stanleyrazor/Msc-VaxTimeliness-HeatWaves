

# libs | dir --------------------------------------------------------------

pacman::p_load(dplyr, tidyr, lubridate, terra, stringr, zoo, purrr, gamlss, gamlss.dist)
theme_use <- theme_bw(base_family = 'Times New Roman') +
  theme(panel.grid = element_line(0))


# Notes -------------------------------------------------------------------

# numpixel*daily for one year gives: 872,910.
# for 85 years it is 74,197,350 rows (not possible to operate this on my machine)
# for simplification we compute tx5x on annual - ignore chaining at year ends

# 2020-12-06 is earliest birth day possible in DHS. Thus, model fitting should only use data before this period.
# use 1940:2019 for model fitting


# Data --------------------------------------------------------------------

ls <- list.files('data/cds/', pattern = "*.nc")
d1 <- map(ls, \(x) {

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
  labs(title = 'Trend of annual maximum Tx5x across 2,385 pixels',
       x = 'Year', y = 'Maximum Tx5x') +
  theme_use

# fitting the GEV on the maxima's per location

x <- data.frame(x = rexp(100000, .2)) |> mutate(mx = -x)
b1 <- brm(x ~ 1,
          data = x, iter = 4e3, chains = 4, cores = 4,
          family = gen_extreme_value())

g1 <- gamlss(mx ~ 1, data = x, family = RGE())




