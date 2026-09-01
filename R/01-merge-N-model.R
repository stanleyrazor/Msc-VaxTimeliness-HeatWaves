
# reading KE-DHS data & trying vaccine timeliness on Measles
pacman::p_load(posterior, tidybayes, rstanarm, marginaleffects, data.table, brms,
               purrr, dplyr, haven, ggplot2, janitor, lubridate, stringr , survival,
               ggsurvfit, icenReg, sf, kableExtra, tidyr, viridisLite, autoReg, flexsurv,
               survey, lubridate, here, scales)
mvs <- naniar::miss_var_summary
source('R/autoReg-modifier.R')

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
master_data <- readRDS("data/processed/master-survey-dataset.rds")

# processed temperature data - at pixel level (buffered)
cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

# processed data on the heat index
hi_data <- arrow::read_parquet('data/processed/heatindex-processed.parquet') |>
  mutate(cluster = as.character(cluster))
setDT(hi_data)

# IN CASE WE NEED IT: processed temperature data - at areal level (zonal aggregation)
# cds_areal <- arrow::read_parquet('data/processed/admin2-processed.parquet')
# cds_areal$heatwave <- cds_areal$p >= .9
#
# cds_areal <- merge(g2, cds_areal, by = 'GID_2', all.x = T)
# cds_areal <- cds_areal |> mutate(cluster = as.character(cluster))
# setDT(cds_areal); setkey(cds_areal, cluster, date)

# Tabling heatwave stats --------------------------------------------------

g <- expand.grid(
  antigen = c('bcg', 'penta1', 'penta2', 'penta3', 'mcv1'),
  return = c(.75, .8, .85, .9, .95, .99),
  final_duedate = c(7, 14, 21, 28),
  full_data_heatwavenum = NA,
  vaxcard_heatwavenum = NA
)

for (i in 1:nrow(g)) {
  antigen <- g[i, 1] |> as.character()
  p <- g[i, 2]

  cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
  cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
  cds_geoloc$heatwave <- cds_geoloc$p >= p
  setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

  vax_data <- ldata[[antigen]]; setDT(vax_data)
  max <- g[i, 3]
  min <- ifelse(antigen == 'bcg', 0, max)

  # the 7-daywindow bound
  vax_data[, `:=`(
    start_dt = due_date - min, #* change to 7
    end_dt   = due_date + max,
    cluster  = as.character(cluster)
  )]

  # find all rows in cds_geoloc where cluster matches AND date is between start/end
  # and sum the 'heatwave' column for @ child
  results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                        .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
                        by = .EACHI]

  vax_data$heatwave <- ifelse(results$heatwave_sum == 0, 0, 1)
  g[i, 4] <- (sum(vax_data$heatwave))
  g[i, 5] <- sum(vax_data$heatwave[!is.na(vax_data$vaxx_date)])

  message(paste(g[i, ], collapse = ' | '))
}

tab <- g |>
  mutate(
    antigen = factor(antigen, levels = c('bcg', 'penta1', 'penta2', 'penta3', 'mcv1'),
                     labels = c('BCG', 'Penta-1', 'Penta-2', 'Penta-3', 'MCV1')),
    return = paste0('Quantile: ', return),
  ) |>
  pivot_wider(
    id_cols = c(antigen, final_duedate),
    names_from = return,
    values_from = c(full_data_heatwavenum, vaxcard_heatwavenum),
    names_vary = "slowest"
  )

kbl(
  tab,
  booktabs = TRUE,
  escape = FALSE,
  align = "ll" %+% strrep("c", ncol(tab) - 2),
  col.names = c(
    "Antigen",
    "Exposure window",
    rep(c("Full data", "VaxCard only"),
        (ncol(tab) - 2) / 2)
  ),
  caption = "Number of observations classified as heatwaves under different return periods."
) |>
  add_header_above(
    c(
      " " = 2,
      setNames(rep(2, length(unique(g$return))),
               paste0('Quantile: ', sort(unique(g$return))))
    )
  ) |>
  kable_classic(full_width = FALSE, html_font = "Times New Roman") |>
  row_spec(0, bold = TRUE)




# Maps on heatwaves on clusters (space-time) ------------------------------

# processed temperature data - at pixel level (buffered)
cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

results <- cds_geoloc |> select(cluster, date, heatwave) |>
  data.frame() |>
  mutate(year = year(date),
         month = month(date),
         month = factor(month, levels = 1:12, labels = month.abb)) |>
  group_by(cluster, year, month) |>
  reframe(hdays = sum(heatwave)) |>
  mutate(
    hdays = ifelse(hdays > 5, '>5', hdays),
    hdays = factor(
      hdays,
      levels = c("0", "1", "2", "3", "4", "5", ">5"),
      ordered = TRUE
    )
  )

cols <- c(
  "0"   = "black",
  setNames(
    viridis(6, option = "D", begin = 0.15, end = 1),
    c("1", "2", "3", "4", "5", ">5"))
)
yrs <- sort(unique(results$year))
for(i in seq_along(yrs)){

  tmp <- merge(g1, filter(results, year == yrs[i]), by = "cluster")

  p <- ggplot() +
    geom_sf(data = shp,
            colour = "grey70",
            fill = NA,
            linewidth = 0.1) +
    geom_sf(data = tmp,
            aes(colour = hdays,
                size = hdays)) +
    facet_wrap(~month, nrow = 3) +
    scale_size_manual(
      values = c(
        "0" = .3,
        "1" = 2,
        "2" = 2,
        "3" = 2,
        "4" = 2,
        "5" = 2,
        ">5" = 2
      ),
      guide = "none"
    ) +
    scale_colour_manual(
      values = cols,
      drop = FALSE,
      name = "Heatwave\ndays"
    ) +
    labs(title = yrs[i]) +
    theme_void(base_family = "Times New Roman") +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(
    paste0("output/img/heatmaps/heat wave/", yrs[i], ".png"),
    p,
    width = 12,
    height = 9,
    dpi = 300
  )
}





# Plotting temperatures and heatwave data ---------------------------------


# processed temperature data - at pixel level (buffered)
cds_geoloc <- arrow::read_parquet('data/processed/admin2-processed.parquet')
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc)

# shapefile
admin2_res <- bind_cols(
  shp |> st_drop_geometry(),
  data.frame(st_centroid(shp) |> st_coordinates()) |> setNames(c('lon', 'lat'))
) |>
  merge(cds_geoloc, by = 'GID_2', all.y = T)

admin2_res <- admin2_res |>
  mutate(
    GID_2 = factor(GID_2,
                   levels = admin2_res |> distinct(GID_2, lat) |> arrange(lat) |> pull(GID_2)
    )
  )

# TX5X plot
mu <- mean(admin2_res$tx5x)
p <- ggplot(admin2_res) +
  geom_tile(aes(x = date, y = GID_2, fill = tx5x)) +
  scale_fill_gradient2(
    low = "#1a9850",
    mid = "white",
    high = "#d73027",
    midpoint = mu,
    breaks = seq(15, 37.5, length.out = 10),
    guide = guide_coloursteps(
      show.limits = TRUE,
      even.steps = TRUE
    ),
    name = expression(TX5X~"("*degree*C*")")
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = levels(admin2_res$GID_2)[c(1, length(levels(admin2_res$GID_2)))],
    labels = c("Southern\nNigeria", "Northern\nNigeria"),
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, title = NULL) + # title = 'The 5-day rolling average temperature'
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),

    legend.key.height = unit(1.4, "cm"),
    legend.key.width = unit(0.8, "cm"),

    axis.text.x = element_text(size = 10, colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),

    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  paste0("output/img/heatmaps/", "tx5x.png"),
  p,
  width = 12,
  height = 9,
  dpi = 300
)

# Threshold plot
admin2_res <- admin2_res |>
  mutate(
    thresh = mev::qgev(
      p = rep(.8, n()),
      loc = loc,
      scale = scale,
      shape = shape
    ),
    GID_2 = as.character(GID_2),
    GID_2 = factor(
      GID_2,
      levels = admin2_res |>
        distinct(GID_2, lat) |>
        arrange(lat) |>
        pull(GID_2)
    ),
    year = factor(year, levels = 2020:2024)
  )
mu <- mean(admin2_res$thresh)


p <- ggplot(admin2_res) +
  geom_tile(aes(x = GID_2, y = year, fill = thresh)) +
  scale_fill_gradient2(
    low = "#1a9850",
    mid = "white",
    high = "#d73027",
    midpoint = mu,
    breaks = seq(floor(min(admin2_res$thresh)) |> round(2),
                 ceiling(max(admin2_res$thresh)) |> round(2),
                 length.out = 10),
    guide = guide_coloursteps(
      show.limits = TRUE,
      even.steps = TRUE
    ),
    name = expression(Threshold~"("*degree*C*")")
  ) +
  scale_x_discrete(
    breaks = levels(admin2_res$GID_2)[c(1, length(levels(admin2_res$GID_2)))],
    labels = c("Southern Nigeria", "Northern Nigeria"),
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL #"Temperature thresholds for a 1-in-4 year event"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.line.x = element_blank(),

    axis.text.y = element_text(size = 11, colour = "black"),
    axis.ticks.y = element_line(colour = "black"),

    legend.key.height = unit(1.4, "cm"),
    legend.key.width = unit(0.8, "cm"),

    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  paste0("output/img/heatmaps/heat wave/", "temp-thresh-4.png"),
  p,
  width = 12,
  height = 4,
  dpi = 300
)

# plotting a line plot:
p <- ggplot(admin2_res |> select(lat, thresh, year) |> distinct(),
            aes(x = lat, y = thresh, col = year)) +
  geom_point(alpha = .3) +
  geom_smooth(method = 'loess', se = F) +
  scale_x_continuous(
    breaks = c(min(admin2_res$lat), max(admin2_res$lat)),
    labels = c("Southern Nigeria", "Northern Nigeria"),
    expand = expansion(mult = c(.02, .02))
  ) +
  labs(
    x = "",
    y = 'Threshold (degree celsius)',
    title = NULL # "Temperature thresholds for a 1-in-4 year event"
  ) +
  theme_bw(base_family = 'Times New Roman') +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
  )

ggsave(
  paste0("output/img/heatmaps/heat wave/", "temp-thresh-4(2).png"),
  p,
  width = 12,
  height = 5,
  dpi = 300
)

# heatwave classification over time, and by location
admin2_res <- admin2_res |>
  mutate(
    GID_2 = as.character(GID_2),
    GID_2 = factor(GID_2,
                   levels = admin2_res |> distinct(GID_2, lat) |> arrange(lat) |> pull(GID_2)
    ),
    heatwave = ifelse(p >= 0.8, 'Heatwave', 'Non-Heatwave')
  )

p <- ggplot(admin2_res) +
  geom_tile(aes(x = date, y = GID_2, fill = heatwave)) +
  scale_fill_manual(
    values = c(
      "Non-Heatwave" = "grey92",
      "Heatwave" = "#d73027"
    )
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = levels(admin2_res$GID_2)[c(1, length(levels(admin2_res$GID_2)))],
    labels = c("Southern\nNigeria", "Northern\nNigeria"),
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, title = NULL, #'The daily heatwave classification',
       fill = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),

    axis.text.x = element_text(size = 10, colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),

    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  paste0("output/img/heatmaps/heat wave/", "heatwave-classification.png"),
  p,
  width = 12,
  height = 9,
  dpi = 300
)


# another way, monthly
tmp <- admin2_res |>
  mutate(
    GID_2 = as.character(GID_2),
    GID_2 = factor(GID_2,
                   levels = admin2_res |> distinct(GID_2, lat) |> arrange(lat) |> pull(GID_2)
    ),
    heatwave = ifelse(p >= 0.8, 1, 0),
    month = month(date)
  ) |>
  group_by(GID_2, year, month) |>
  reframe(heatwave = sum(heatwave)) |>
  mutate(
    heatwave = ifelse(heatwave > 5, '>5', heatwave),
    heatwave = factor(
      heatwave,
      levels = c("0", "1", "2", "3", "4", "5", ">5"),
      ordered = TRUE
    ),
    date = ym(paste0(year, '-', month))
  )

cols <- c(
  "0"   = "gray",
  setNames(
    viridis(6, option = "D", begin = 0.15, end = 1),
    c("1", "2", "3", "4", "5", ">5"))
)

p <- ggplot(tmp) +
  geom_tile(aes(x = date, y = GID_2, fill = heatwave), col = NA, linewidth = 0) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = cols,
    drop = FALSE,
    name = "Number of\nheatwave days"
  ) +
  labs(x = NULL, title = NULL, #'The daily heatwave classification',
       fill = NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    axis.text.x = element_text(size = 10, colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),

    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  paste0("output/img/heatmaps/heat wave/", "heatwave-classification(agg).png"),
  p,
  width = 12,
  height = 9,
  dpi = 300
)

# plotting map of num days heatwave exposed
cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc)

clusterloc <- merge(g1,
                    cds_geoloc |>
                      group_by(year, cluster) |>
                      reframe(heatwave = sum(heatwave))
                    , by = 'cluster', all = T) |>
  mutate(
    class = case_when(
      heatwave == 0                ~ "0",
      heatwave >= 1 & heatwave <= 4 ~ "1-4",
      heatwave >= 5 & heatwave <= 9 ~ "5-9",
      heatwave >= 10 & heatwave <= 14 ~ "10-14",
      heatwave >= 15 & heatwave <= 19 ~ "15-19",
      heatwave >= 20               ~ ">=20",
      TRUE                         ~ NA_character_ # Safely handles any unexpected NAs
    ),
    # Convert to an ordered factor so R knows the correct sequence
    class = factor(
      class,
      levels = c("0", "1-4", "5-9", "10-14", "15-19", ">=20"),
      ordered = TRUE
    )
  )

cols <- c(
  "0"   = "black",
  setNames(
    viridis(6, option = "D", begin = 0.15, end = 1),
    c("0", "1-4", "5-9", "10-14", "15-19", ">=20"))
)

p <- ggplot() +
  geom_sf(data = shp,
          colour = "grey70",
          fill = NA,
          linewidth = 0.1) +
  geom_sf(data = clusterloc,
          aes(colour = class,
              size = class)) +
  facet_wrap(~year, nrow = 2) +
  scale_size_manual(
    values = c(
      "0" = .3,
      "1-4" = 2,
      "5-9" = 2,
      "10-14" = 2,
      "15-19" = 2,
      ">=20" = 2
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = cols,
    drop = FALSE,
    name = "Heatwave frequency\n(days)"
  ) +
  theme_void(base_family = "Times New Roman") +
  theme(
    plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt"),

    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
  # coord_sf(expand = FALSE)

ggsave(
  paste0("output/img/heatmaps/heat wave/heatwave-map.png"),
  p,
  width = 12,
  height = 9,
  dpi = 300
)




# Heat index exploration --------------------------------------------------

# shapefile
admin2_res <- bind_cols(
  data.frame(st_centroid(g1) |> st_coordinates()),
  cluster = g1$cluster
) |>
  setNames(c('lon', 'lat', 'cluster')) |>
  merge(hi_data |>
          mutate(
            heatclass = case_when(
              heatindex < 26.67                       ~ "Tolerable",
              heatindex >= 26.67 & heatindex < 32.22  ~ "Caution",
              heatindex >= 32.22 & heatindex < 39.44  ~ "Extreme Caution",
              heatindex >= 39.44 & heatindex < 51.67  ~ "Danger",
              heatindex >= 51.67                      ~ "Extreme Danger",
              TRUE                                    ~ NA_character_
            ),
            heatclass = factor(
              heatclass,
              levels = c("Tolerable", "Caution", "Extreme Caution",
                         "Danger", "Extreme Danger"
              ),
              labels = c("Tolerable", "Caution", "Extreme Caution",
                         "Danger", "Extreme Danger"
              ),
              ordered = T
            )
          ),
        by = 'cluster', all = T)

admin2_res <- admin2_res |>
  mutate(
    cluster = factor(cluster,
                   levels = admin2_res |> distinct(cluster, lat) |> arrange(lat) |> pull(cluster)
    )
  )

# heatindex plot
mu <- mean(admin2_res$heatindex)
p <- ggplot(admin2_res) +
  geom_tile(aes(x = date, y = cluster, fill = heatindex)) +
  scale_fill_gradient2(
    low = "#1a9850",
    mid = "white",
    high = "#d73027",
    midpoint = mu,
    breaks = seq(15, 37.5, length.out = 10),
    guide = guide_coloursteps(
      show.limits = TRUE,
      even.steps = TRUE
    ),
    name = expression(Heat-Index~"("*degree*C*")")
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = levels(admin2_res$cluster)[c(1, length(levels(admin2_res$cluster)))],
    labels = c("Southern\nNigeria", "Northern\nNigeria"),
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, title = 'The Heat Index') +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),

    legend.key.height = unit(1.4, "cm"),
    legend.key.width = unit(0.8, "cm"),

    axis.text.x = element_text(size = 10, colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),

    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  paste0("output/img/heatmaps/", "heat-index-num.png"),
  p,
  width = 12,
  height = 9,
  dpi = 300
)


# heat index classification
p <- ggplot(admin2_res) +
  geom_tile(aes(x = date, y = cluster, fill = heatclass)) +
  # Swapped continuous gradient for explicit NWS categorical palette
  scale_fill_manual(
    values = c(
      "Tolerable"       = "#a6d96a",  # Soft NWS Light Green
      "Caution"         = "#ffff00",  # NWS Yellow
      "Extreme Caution" = "#ffaa00",  # NWS Amber / Dark Yellow
      "Danger"          = "#ff6600",  # NWS Orange
      "Extreme Danger"  = "#cc0000"   # NWS Red
    ),
    name = "Heat Index\nCategory",
    drop = FALSE # Keeps all levels in the legend even if missing in data subset
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = levels(admin2_res$cluster)[c(1, length(levels(admin2_res$cluster)))],
    labels = c("Southern\nNigeria", "Northern\nNigeria"),
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, title = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.title.y = element_blank(),

    # Updated legend adjustments for discrete blocks
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.8, "cm"),

    axis.text.x = element_text(size = 10, colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),

    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  paste0("output/img/heatmaps/", "heat-index-class.png"),
  p,
  width = 12,
  height = 9,
  dpi = 300
)




# Kaplan meier analysis ---------------------------------------------------

cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1')
usenames <- c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')
antigen_times <- pmax(1, c(0, c(6, 10, 14)*7, 6*30.4, 9*30.4, 15*30.4) - 7)
plt <- list()

for (i in 1:length(antigen)) {

  use_window_min <- ifelse(antigen[i] == 'bcg', 0, 7)
  vax_data <- ldata[[antigen[i]]]; setDT(vax_data)

  # the 7-daywindow bound
  vax_data[, `:=`(
    start_dt = due_date - use_window_min,
    end_dt   = due_date + 7,
    cluster  = as.character(cluster)
  )]

  # find all rows in cds_geoloc where cluster matches AND date is between start/end
  # and sum the 'heatwave' column for @ child
  results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                        .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
                        by = .EACHI] |>
    setNames(c('cluster', 'start_dt', 'end_dt', 'heatwave')) |>
    distinct()
  tx5x <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                     .(tx5x = mean(tx5x, na.rm = TRUE)),
                     by = .EACHI] |>
    setNames(c('cluster', 'start_dt', 'end_dt', 'tx5x')) |>
    distinct()
  heatindex <- hi_data[vax_data,
                       on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                       .(heatindex = mean(heatindex, na.rm = TRUE)),
                       by = .EACHI] |>
    setnames(c('cluster', 'start_dt', 'end_dt', 'heatindex')) |>
    unique()
  # summary(heatindex$heatindex)

  # merging to have the heatwave and temperature variable
  vax_data <- merge(
    vax_data, results,
    by = c('cluster', 'start_dt', 'end_dt'),
    all.x = T
  ) |>
    merge(heatindex, by = c('cluster', 'start_dt', 'end_dt'), all.x = T) |>
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

      heatindex = case_when(
        heatindex < 26.67                       ~ "Tolerable",
        heatindex >= 26.67 & heatindex < 32.22  ~ "Caution",
        heatindex >= 32.22 & heatindex < 39.44  ~ "Extreme Caution",
        # heatindex >= 39.44 & heatindex < 51.67  ~ "Danger",
        # heatindex >= 51.67                      ~ "Extreme Danger",
        TRUE                                    ~ NA_character_
      ),
      heatindex = factor(
        heatindex,
        levels = c("Tolerable", "Caution", "Extreme Caution"
                   # "Danger", "Extreme Danger"
        ),
        labels = c("Tolerable", "Caution", "Extreme Caution"
                   # "Danger", "Extreme Danger"
        )
      ),

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
    model_data$heatindex         <- setLabel(model_data$heatindex, "Heat Index")
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
  # mstd <- merge(master_data,
  #               model_data |>
  #                 mutate(across(c(bidx, v021, v022), as.character),
  #                        selector = 'Include') |>
  #                 select(-wt),
  #               by = c('caseid', 'bidx', 'v021', 'v022'), all.x = T) |>
  #   mutate(selector = ifelse(is.na(selector), 'Not included', selector))
  #
  # full_surv_design <- svydesign(
  #   id = ~v021,          # Primary Sampling Unit / Cluster
  #   strata = ~v022,      # Stratification variable/022
  #   weights = ~wt,     # DHS weight (remember to divide v005 by 1,000,000 first)
  #   data = mstd,
  #   nest = TRUE
  # )
  # surv_design <- subset(full_surv_design, selector == 'Include')
  # options(survey.lonely.psu = 'adjust')
  #
  sv_hw <- survfit(Surv(time_start, time_end, status, type = "interval") ~ heatwave,
                   data = model_data, weights = wt)
  sv_hi <- survfit(Surv(time_start, time_end, status, type = "interval") ~ heatindex,
                   data = model_data, weights = wt)

  km_df_hw <- tidy(sv_hw) |> subset(time > 0 & estimate > 0 & estimate < 1)
  km_df_hi <- tidy(sv_hi) |> subset(time > 0 & estimate > 0 & estimate < 1)

  # 3. Create the survey-weighted log-log plot
  hw_plt <- ggplot(km_df_hw, aes(x = log(time), y = log(-log(estimate)), color = strata)) +
    geom_step(linewidth = 1) +
    theme_bw() +
    labs(
      col = NULL,
      title = NULL,
      x = "Log of Time (log(t))",
      y = "Log-Log Cumulative Hazard\n(log(-log(S(t))))"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  hi_plt <- ggplot(km_df_hi, aes(x = log(time), y = log(-log(estimate)), color = strata)) +
    geom_step(linewidth = 1) +
    theme_minimal() +
    labs(
      col = NULL,
      title = NULL,
      x = "Log of Time (log(t))",
      y = "Log-Log Cumulative Hazard\n(log(-log(S(t))))"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  plt[[i]] <- patchwork::wrap_plots(hw_plt, hi_plt, nrow = 1) + ggtitle(usenames[i])

}

wrap_plots(plt, nrow = 6)
ggsave('output/img/sensitivity/log-log-weibull-plots.png',
       height = 20, width = 14, dpi = 300)

# Models ------------------------------------------------------------------

cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .75
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1', 'mcv2')
antigen_times <- pmax(1, c(0, c(6, 10, 14)*7, 6*30.4, 9*30.4, 15*30.4) - 7)

for (i in 1:length(antigen)) {

  use_window_min <- ifelse(antigen[i] == 'bcg', 0, 7)
  vax_data <- ldata[[antigen[i]]]; setDT(vax_data)

  # the 7-daywindow bound
  vax_data[, `:=`(
    start_dt = due_date - use_window_min,
    end_dt   = due_date + 7,
    cluster  = as.character(cluster)
  )]

  # find all rows in cds_geoloc where cluster matches AND date is between start/end
  # and sum the 'heatwave' column for @ child
  results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                        .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
                        by = .EACHI] |>
    setNames(c('cluster', 'start_dt', 'end_dt', 'heatwave')) |>
    distinct()
  tx5x <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                     .(tx5x = mean(tx5x, na.rm = TRUE)),
                     by = .EACHI] |>
    setNames(c('cluster', 'start_dt', 'end_dt', 'tx5x')) |>
    distinct()
  heatindex <- hi_data[vax_data,
                       on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                       .(heatindex = mean(heatindex, na.rm = TRUE)),
                       by = .EACHI] |>
    setnames(c('cluster', 'start_dt', 'end_dt', 'heatindex')) |>
    unique()
  # summary(heatindex$heatindex)

  # merging to have the heatwave and temperature variable
  vax_data <- merge(
    vax_data, results,
    by = c('cluster', 'start_dt', 'end_dt'),
    all.x = T
  ) |>
    merge(heatindex, by = c('cluster', 'start_dt', 'end_dt'), all.x = T) |>
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

      heatindex = case_when(
        heatindex < 26.67                       ~ "Tolerable",
        heatindex >= 26.67 & heatindex < 32.22  ~ "Caution",
        heatindex >= 32.22 & heatindex < 39.44  ~ "Extreme Caution",
        # heatindex >= 39.44 & heatindex < 51.67  ~ "Danger",
        # heatindex >= 51.67                      ~ "Extreme Danger",
        TRUE                                    ~ NA_character_
      ),
      heatindex = factor(
        heatindex,
        levels = c("Tolerable", "Caution", "Extreme Caution"
                   # "Danger", "Extreme Danger"
        ),
        labels = c("Tolerable", "Caution", "Extreme Caution"
                   # "Danger", "Extreme Danger"
        )
      ),

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
    model_data$heatindex         <- setLabel(model_data$heatindex, "Heat Index")
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

  # fit <- survreg(
  #   Surv(time_start, time_end, status, type = "interval") ~
  #     heatwave + residence + birth_order + place_delivery + child_gender +
  #     time_to_hf + wealth_index + meduc + mother_age_group,
  #   weights = wt,
  #   dist = 'weibull',
  #   data = model_data,
  #   robust = TRUE,
  # )

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
  # surv_design <- svydesign(
  #   id = ~v021,          # Primary Sampling Unit / Cluster
  #   strata = ~v022,      # Stratification variable/022
  #   weights = ~wt,     # DHS weight (remember to divide v005 by 1,000,000 first)
  #   data = model_data,
  #   nest = TRUE
  # )
  options(survey.lonely.psu = 'adjust')

  fit_hw <- svysurvreg(
    Surv(time_start, time_end, status, type = "interval") ~
      heatwave + residence + birth_order + place_delivery + child_gender +
      time_to_hf + wealth_index + meduc + mother_age_group + geozone,
    dist = 'weibull',
    design = surv_design,
    data = mstd
  )

  fit_hi <- svysurvreg(
    Surv(time_start, time_end, status, type = "interval") ~
      heatindex + residence + birth_order + place_delivery + child_gender +
      time_to_hf + wealth_index + meduc + mother_age_group + geozone,
    dist = 'weibull',
    design = surv_design,
    data = mstd
  )

  saveRDS(list(fit_hw = fit_hw, fit_hi = fit_hi),
          paste0('output/models/', antigen[i], '-file.rds'))

  result <- autoReg(fit_hw, uni = F)
  flextable::save_as_image(
    x = result %>% myft(),
    path = paste0('output/img/', antigen[i], '/HW-model-table.png'),
    res = 500
  )

  result <- autoReg(fit_hi, uni = F)
  flextable::save_as_image(
    x = result %>% myft(),
    path = paste0('output/img/', antigen[i], '/HI-model-table.png'),
    res = 500
  )

  p <- mod_modelPlot(fit_hw)
  png(
    filename = paste0("output/img/", antigen[i], "/HW-model-plot.png"),
    width = 12,
    height = 9,
    units = "in",
    res = 500
  )
  print(p)
  dev.off()

  p <- mod_modelPlot(fit_hi)
  png(
    filename = paste0("output/img/", antigen[i], "/HI-model-plot.png"),
    width = 12,
    height = 9,
    units = "in",
    res = 500
  )
  print(p)
  dev.off()

}



# Cluster locations -------------------------------------------------------

{

  g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
    select(cluster = DHSCLUST, residence = URBAN_RURA) |>
    mutate(residence = ifelse(residence == 'U', 'Urban', 'Rural'))

  shp <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> terra::unwrap() |>
    st_as_sf() |> select(GID_2)

  p <- ggplot() +
    geom_sf(data = shp,
            colour = "grey70",
            fill = NA,
            linewidth = 0.1) +
    geom_sf(data = g1, aes(colour = residence)) +
    labs(color = 'Residence') +
    theme_void(base_family = "Times New Roman") +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ## optionally save
  ggsave(
    paste0("output/img/cluster-locations.png"),
    p,
    width = 7,
    height = 7,
    dpi = 300
  )
}
