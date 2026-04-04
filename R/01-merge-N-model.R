
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
cds_areal$heatwave <- cds_areal$p >= .9

cds_areal <- merge(g2, cds_areal, by = 'GID_2', all.x = T)
cds_areal <- cds_areal |> mutate(cluster = as.character(cluster))
setDT(cds_areal); setkey(cds_areal, cluster, date)

# -------------------------------------------------------------------------

# was penta3
vax_data <- ldata[['penta3']]; setDT(vax_data)

# the 14-day window bound
vax_data[, `:=`(
  start_dt = due_date - 7,
  end_dt   = due_date + 7,
  cluster  = as.character(cluster)
)]

# find all rows in cds_geoloc where cluster matches AND date is between start/end
# and sum the 'heatwave' column for @ child
# results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
#                .(heatwave_sum = sum(heatwave, na.rm = TRUE)),
#                by = .EACHI]
results <- cds_geoloc[vax_data, on = .(cluster = cluster, date >= start_dt, date <= end_dt),
                      .(p = mean(p, na.rm = TRUE)),
                      by = .EACHI]

# vax_data$heatwave <- ifelse(results$heatwave_sum == 0, 'absent', 'present')
vax_data$heatwave <- results$p
(table(vax_data$heatwave))

# vax_data <- as.data.frame(vax_data) |>
#   mutate(meduc = as.character(meduc) |> as.factor())

# Model -------------------------------------------------------------------

b3 <- brm(
  time_outcome | weights(wt) + cens(outcome_event) ~ heatwave + residence + (1 | caseid),

  family = weibull(),
  data = vax_data |> mutate(time_outcome = ifelse(time_outcome == 0, 1, time_outcome)),

  chains = 4,
  iter = 2000,
  warmup = 1000
)

summary(b3)
conditional_effects(b3, "heatwave")
bayesplot::pp_check(b3)

# conditional effects for timeliness
# type response - ignores individual residual randomness | prediction - incorporates it
# rvar format for posterior package
pred <- avg_predictions(b3, newdata = vax_data,
                        re_formula = NULL, type = 'response', by = "heatwave", wts = 'wt')
draws <- get_draws(pred, "rvar")

quantile2(draws$rvar, c(0.025, 0.5, 0.975)) # posterior quantiles
E(draws$rvar) # expected value of posterior dist
Pr(draws$rvar <= ((4 * 7))) # posterior mass below 9 months + 2 week buffer

# Posterior predictions and quantiles
pp <- posterior_predict(b3, newdata = vax_data, re_formula = NULL)
tvax_mean <- colMeans(pp)
aggregate(tvax_mean ~ vax_data$heatwave, FUN = mean)

timely_draws <- pp <= 28
timely_prob <- colMeans(timely_draws)
aggregate(timely_prob ~ vax_data$heatwave, FUN = mean)


predictions(
  b3,
  newdata = datagrid(heatwave = c('absent', 'present'),
                     caseid = unique),
  by = "heatwave",
  re_formula = NULL
)

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
x <- vax_data |> mutate(pred = colMeans(posterior_epred(b3))) |> select(time_outcome, pred, everything())
abs <- vax_data |> mutate(heatwave = 'absent')
prs <- vax_data |> mutate(heatwave = 'present')

abspred <- posterior_epred(b3, newdata = abs) |> colMeans() |> as.numeric()
prspred <- posterior_epred(b3, newdata = prs) |> colMeans() |> as.numeric()

weighted.mean(abspred, vax_data$wt)
weighted.mean(prspred, vax_data$wt)

mean(abspred)
mean(prspred)


# 1. Use epred to get the 4000 x N matrix
p_matrix <- posterior_epred(b3, newdata = prs)

# 2. Calculate the weighted mean for EACH of the 4000 draws
draw_means <- apply(p_matrix, 1, \(row) weighted.mean(row, w = vax_data$wt))

# 3. Take the mean of those 4000 results
mean(draw_means)
# This should now match 1095 much more closely





