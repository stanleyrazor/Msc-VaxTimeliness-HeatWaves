
# reading KE-DHS data & trying vaccine timeliness on Measles
pacman::p_load(posterior, tidybayes, rstanarm, brms, ggpubr, dplyr, haven, ggplot2, janitor, lubridate, stringr)
mvs <- naniar::miss_var_summary

# Data --------------------------------------------------------------------

load("~/Documents/GitHub/KDHS Model Updates/data/KEBR8AFL_2022.RData")
d1 <- birth; rm(birth)

d2 <- d1 |>
  select(caseid, v001, v005, v024, v025, b18,
         v190, m15, m14, bord, b4, b8, v106, v136, v130, v131, v501, v717,
         h9, h9d, h9m, h9y) |>

  # omitting "don't know" and "inconsistent" - 1 obs
  mutate(across(starts_with('h9'), \(x) ifelse(x %in% 97:98, NA, x))) |>

  # omitting NA on h9
  filter(!is.na(h9d)) |>

  # computing birth date
  mutate(bday = as.Date(b18 - 1, origin = "1900-01-01"),
         vday = make_date(h9y, h9m, h9d),
         due_date = (bday %m+% months(9)), # + weeks(4),
         delay = as.integer(pmax(0, as.numeric(vday - due_date)))) |>

  # dealing with covariates
  mutate(v005 = v005 / 1e6,
         m15 = ifelse(m15 %in% c(10,11,12,96), 'home', 'facility'),
         m14 = ifelse(m14 == 98, NA, m14),
         across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)) |>

  select(caseid, wt = v005, county = v024, residence = v025, bday, vday, due_date, delay,
         wealth = v190, sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
         meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
         curr_marital = v501, occupation = v717) |>
  select(-occupation, -curr_marital, -ethnicity, -religion) |>
  na.omit()
mvs(d2)


m1 <- glm(delay ~ residence + wealth + sex + bord + anc + delivery +
      meduc + hhsize,
    weights = wt,
    family = poisson(),
    data = d2)

summary(m1)

ggdensity(simulate(m1, 10) |> mutate(across(everything(), \(x) ifelse(x > 28, 1, 0))) |>
            rowMeans())


# Trying out a cox model --------------------------------------------------

d3 <- d1 %>%

  filter(b19 < 36, b5 == 1, h9 != 8) %>%

  mutate(
    eligible = (b19 > 9),

    iday = as.Date(v008a - 1, origin = "1900-01-01"),
    bday = as.Date(b18 - 1, origin = "1900-01-01"),
    vday = make_date(h9y, h9m, h9d),

    event_time = as.integer(as.numeric(vday - bday)),   # NA if missing vday
    censor_time = as.integer(as.numeric(iday - bday)),

    time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
    time_outcome = pmax(0, time_outcome),

    # descriptive event label
    event_label = case_when(
      h9 == 1           ~ 0,   # card with date
      h9 %in% c(2,3)    ~ -1,    # reported / card but no date
      h9 == 0           ~ 1     # no vaxx
    )
  ) |>

  # dealing with covariates
  mutate(v005 = v005 / 1e6,
         m15 = ifelse(m15 %in% c(10,11,12,96), 'home', 'facility'),
         m14 = ifelse(m14 == 98, NA, m14),
         across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)) |>

  select(caseid, wt = v005, county = v024, residence = v025, eligible, time_outcome, event_label,
         wealth = v190, sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
         meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
         curr_marital = v501, occupation = v717) |>
  select(-occupation, -curr_marital, -ethnicity, -religion) |>
  na.omit()
mvs(d3)

b1 <- brm(time_outcome | weights(wt) + cens(event_label) ~
            eligible + residence + wealth + sex + bord + anc + delivery + meduc + hhsize,
          family = cox(),
          data = d3,

          chains = 4,
          iter = 2000,
          warmup = 1000
)

# caseid for mother, eligible for whether eligible to receive the vax
b2 <- brm(time_outcome | weights(wt) + cens(event_label) ~
            eligible + residence + wealth + sex + bord + anc + delivery + meduc + hhsize,
          family = weibull(),
          data = d3 |> mutate(time_outcome = ifelse(time_outcome == 0, 1, time_outcome)),

          chains = 4,
          iter = 1000,
          warmup = 500
)

b3 <- brm(time_outcome | weights(wt) + cens(event_label) ~
            eligible + residence + wealth + sex + bord + anc + delivery + meduc + hhsize +
            (1 | caseid),
          family = weibull(),
          data = d3 |> mutate(time_outcome = ifelse(time_outcome == 0, 1, time_outcome)),

          chains = 4,
          iter = 1000,
          warmup = 500
)

d3 |> mutate(time_outcome = ifelse(time_outcome == 0, 1, time_outcome)) |>
  mutate(pred = predict(b3),
         timeliness = timeliness)


fit <- posterior_predict(b3)
timeliness_posterior <- fit |>
  t() |>
  data.frame() |>
  mutate(across(everything(), \(x) x < 365.25))

timeliness <- rowMeans(timeliness_posterior)
hist(timeliness, breaks=60)


d4 <- d3 |> mutate(time_outcome = ifelse(time_outcome == 0, 1, time_outcome)) |>
  mutate(pred = predict(b3)[, 1],
         timeliness = timeliness) |>
  select(caseid, county, residence, eligible, time_outcome, event_label,
         pred, timeliness, everything())


ggplot(d4) +
  geom_point(aes(x = pred, y = time_outcome, col = eligible), pch = '*') +
  geom_abline(aes(intercept = 0, slope = 1), col = 'red', lty = 'dashed') +
  geom_hline(aes(yintercept = 365.25), col = 'black', lwd = .4) +
  geom_vline(aes(xintercept = 365.25), col = 'black', lwd = .4) +
  facet_wrap(~event_label) +
  theme_bw() + theme(panel.grid = element_line(0))



conditional_effects(b2, effects = 'anc')
conditional_effects(b2, effects = 'anc', method = 'posterior_epred',
                    spaghetti = T, ndraws=10, mean=T, rug=T,
                    theme = 'theme_classic')

