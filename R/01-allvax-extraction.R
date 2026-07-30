


# function for processing one & multiple vaccines -------------------------

make_vax_one <- function(d1, vaccine, identifier, period, timeframe = c("weeks", "months")) {
  timeframe <- match.arg(timeframe)
  period <- as.numeric(period)

  #* age_min <- if (timeframe == "weeks") floor((period * 7) / 30.4375) else period
  age_min <- ifelse(timeframe == 'weeks', floor(period * 7), floor(period * 30.4375))

  due_expr <- if (timeframe == "weeks") {
    expr(birth_date %m+% weeks(!!period))
  } else {
    expr(birth_date %m+% months(!!period))
  }

  ycol <- paste0(identifier, "y")
  mcol <- paste0(identifier, "m")
  dcol <- paste0(identifier, "d")

  d1 |>

    # computing ages, in order to threshold by age
    mutate(
      interview_date = as.Date(v008a - 1, origin = "1900-01-01"),
      birth_date = as.Date(b18 - 1, origin = "1900-01-01"),
      age_days = as.integer(interview_date - birth_date)
    ) |>

    #* filter(b19 <= 36, b19 >= age_min, b5 == 1, .data[[identifier]] != 8) |>
    filter(b19 <= 36, age_days >= age_min, b5 == 1, .data[[identifier]] != 8) |>

    mutate(
      across(starts_with(identifier), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x)),
      across(c(v024, v025), as_factor),
      v005 = v005 / 1e6,

      vaxx_date = make_date(.data[[ycol]], .data[[mcol]], .data[[dcol]]),
      due_date = !!due_expr,

      event_time = as.integer(as.numeric(vaxx_date - birth_date)),
      censor_time = as.integer(as.numeric(interview_date - birth_date)),

      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # 0 - event | -1 - left censored | 1 - right censored
      outcome_event = case_when(
        .data[[identifier]] == 1 ~ 0L,
        .data[[identifier]] %in% c(2, 3) ~ -1L,
        .data[[identifier]] == 0 ~ 1L,
        TRUE ~ NA_integer_
      ),

      outcome_event = case_when(
        is.na(vaxx_date) & .data[[identifier]] == 1 ~ -1L,
        TRUE ~ outcome_event
      )
    ) |>
    select(
      caseid, bidx, wt = v005, cluster = v001, residence = v025, admin = v024,
      interview_date, birth_date, vaxx_date, due_date,
      event_time, censor_time, time_outcome, outcome_event
    )

}

make_vax_all <- function(d1, spec) {
  out <- vector("list", nrow(spec))

  for (i in seq_len(nrow(spec))) {
    out[[i]] <- make_vax_one(
      d1 = d1,
      vaccine = spec$vaccine[i],
      identifier = spec$identifier[i],
      period = spec$period[i],
      timeframe = spec$timeframe[i]
    )
  }

  names(out) <- spec$vaccine
  out
}


# Data --------------------------------------------------------------------

spec <- data.frame(
  vaccine     = c(
    "bcg", "opv0", "hepb0",
    "penta1", "pneumo1", "opv1",
    "penta2", "pneumo2", "opv2",
    "penta3", "pneumo3", "opv3",
    "vita1", "mcv1", "mcv2"
  ),
  identifier  = c(
    "h2", "h0", "h50",
    "h51", "h54", "h4",
    "h52", "h55", "h6",
    "h53", "h56", "h8",
    "h33", "h9", "h9a"
  ),
  period      = c(
    0, 0, 0,     # week 0
    6, 6, 6,     # week 6
    10, 10, 10,  # week 10
    14, 14, 14,  # week 14
    6,           # 6 months
    9,           # 9 months
    15           # 15 months
  ),
  timeframe   = c(
    "weeks","weeks","weeks",
    "weeks","weeks","weeks",
    "weeks","weeks","weeks",
    "weeks","weeks","weeks",
    "months","months","months"
  ),
  stringsAsFactors = FALSE
)

# Results -----------------------------------------------------------------

# d1 - birth recode from dhs
# spec - the spec file
d1 <- read_dta("data/dhs/NG_2024_DHS_03262026_919_211396/NGBR8BDT/NGBR8BFL.dta")
vax_data <- make_vax_all(d1, spec)

saveRDS(vax_data, "data/processed/vaxdata-components.rds")

# creating the master dataset with only necessqry columns for survey design
d2 <- d1 |>
  select(v021, v022, caseid, bidx, wt = v005) |>
  as_factor() |>
  mutate(
    across(c(bidx, v021, v022), as.character),
    wt = as.numeric(wt) / 1e6)
saveRDS(d2, "data/processed/master-survey-dataset.rds")

# compatibility -----------------------------------------------------------

# vax_dfs |>
#   lapply(FUN = nrow) |>
#   unlist()
#
# ldata |>
#   lapply(FUN = nrow) |>
#   unlist()
#
# # pr(time < 1e3)
# prop_less_1e3 <- \(x) mean(x$time_outcome < 1e3)
# vax_dfs |>
#   lapply(FUN = prop_less_1e3) |>
#   unlist()
#
# ldata |>
#   lapply(FUN = prop_less_1e3) |>
#   unlist()
