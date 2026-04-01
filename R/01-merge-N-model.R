
# reading KE-DHS data & trying vaccine timeliness on Measles
pacman::p_load(posterior, tidybayes, rstanarm, marginaleffects, data.table, brms,
               ggpubr, dplyr, haven, ggplot2, janitor, lubridate, stringr)
mvs <- naniar::miss_var_summary

# Data --------------------------------------------------------------------

# DHS Geo data - spatial join — same CRS etc.
g1 <- st_read("data/dhs/NG_2024_DHS_03262026_919_211396/NGGE8AFL/NGGE8AFL.shp") |>
  select(cluster = DHSCLUST)
shp <- readRDS('data/shp/gadm/gadm41_NGA_2_pk.rds') |> st_as_sf() |> select(GID_2)
g2 <- st_join(g1, shp[, "GID_2"], left = TRUE) |> st_drop_geometry()


# vax-data: saved preciously
ldata <- readRDS('data/processed/vaxdata-components.rds')

# processed temperature data - at pixel level (buffered)
cds_geoloc <- arrow::read_parquet('data/processed/cluster-processed.parquet')
cds_geoloc <- cds_geoloc |> mutate(cluster = as.character(cluster))
cds_geoloc$heatwave <- cds_geoloc$p >= .9
setDT(cds_geoloc); setkey(cds_geoloc, cluster, date)

# processed temperature data - at areal level (zonal aggregation)
cds_areal <- arrow::read_parquet('data/processed/admin2-processed.parquet')
cds_areal$heatwave <- cds_areal$p >= .95

cds_areal <- merge(g2, cds_areal, by = 'GID_2', all.x = T)
cds_areal <- cds_areal |> mutate(cluster = as.character(cluster))
setDT(cds_areal); setkey(cds_areal, cluster, date)

# -------------------------------------------------------------------------

# was penta3
bcg <- ldata[['mcv2']]; setDT(bcg)

# the 14-day window bound
bcg[, `:=`(
  start_dt = due_date - 7,
  end_dt   = due_date + 7,
  cluster  = as.character(cluster)
)]

# find all rows in cds_geoloc where cluster matches AND date is between start/end
# and sum the 'heatwave' column for @ child
results <- cds_areal[bcg, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
               .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
               by = .EACHI]

bcg$heatwave <- ifelse(results$heatwave_sum == 0, 'absent', 'present')
(table(bcg$heatwave))

bcg <- as.data.frame(bcg) |>
  mutate(meduc = as.character(meduc) |> as.factor())

# Model -------------------------------------------------------------------

b3 <- brm(time_outcome | weights(wt) + cens(outcome_event) ~
            heatwave + residence + wealth + sex + bord + delivery + meduc + hhsize,
          family = weibull(),
          data = bcg |> mutate(time_outcome = ifelse(time_outcome == 0, .01, time_outcome)),

          chains = 4,
          iter = 1000,
          warmup = 500
)

summary(b3)
conditional_effects(b3, "heatwave")

# conditional effects for timeliness
pred <- avg_predictions(b3, by = "heatwave", wts = 'wt')
draws <- get_draws(pred, "rvar") # rvar format for posterior package

E(draws$rvar) # expected value of posterior dist
quantile2(draws$rvar, c(0.05, 0.5, 0.95)) # posterior quantiles
Pr(draws$rvar <= (16*7)) # posterior mass below 16*7

# Marginal effects --------------------------------------------------------

library(marginaleffects)

# define a function that returns timeliness
fn_timely <- function(vec) as.integer(vec <= (16*7))

# Average predicted timeliness
a1 <- avg_predictions(b3, variables = "heatwave", type = "response", wts = 'wt'); a1
a2 <- avg_predictions(b3, variables = "heatwave", type = "response",
                      transform = fn_timely, wts = 'wt'); a2

threshold <- (16 * 7)
# 112 days = 16 weeks
threshold <- 112

a2 <- avg_predictions(
  b3, variables = "heatwave",
  comparison = function(x, ...) mean(x <= threshold),
  wts = "wt"
)
a2

# Manual computation of marginal effects: averaging over.
abs <- bcg |> mutate(heatwave = 'absent')
prs <- bcg |> mutate(heatwave = 'present')

abspred <- posterior_epred(b3, newdata = abs) |> colMeans() |> as.numeric()
prspred <- posterior_epred(b3, newdata = prs) |> colMeans() |> as.numeric()

weighted.mean(abspred, bcg$wt)
weighted.mean(prspred, bcg$wt)

mean(abspred)
mean(prspred)


# 1. Use epred to get the 4000 x N matrix
p_matrix <- posterior_epred(b3, newdata = prs)

# 2. Calculate the weighted mean for EACH of the 4000 draws
draw_means <- apply(p_matrix, 1, \(row) weighted.mean(row, w = bcg$wt))

# 3. Take the mean of those 4000 results
mean(draw_means)
# This should now match 1095 much more closely





