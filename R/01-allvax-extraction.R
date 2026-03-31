
# Birth vaccines ----------------------------------------------------------

# BCG - h2 | OPV - h0 | Hep.B - h50
d3 <- d1 |>

  filter(b19 <= 36 & b5 == 1) |>

  select(caseid, b18, h0, h0d, h0m, h0y, h2, h2d, h2m, h2y, h50, h50d, h50m, h50y) |>

  mutate(
    bday = as.Date(b18 - 1, origin = "1900-01-01"),
    opv = make_date(h0y, h0m, h0d),
    bcg = make_date(h2y, h2m, h2d),
    hep = make_date(h50y, h50m, h50d)
  ) |>

  select(h0, h2, h50, bday, opv, bcg, hep)


{
  bcg <- d1 |>

    # 14 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((0*7)/30.4375) & b5 == 1 & h2 != 8) |>
    mutate(across(starts_with('h2'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h2y, h2m, h2d),

      due_date = (bday %m+% weeks(0)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h2 == 1 ~ 0,  # card with date
        h2 %in% c(2, 3) ~ -1, # reported / card but no date
        h2 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h2d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h2 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h2,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)

}

{
  opv0 <- d1 |>

    # 14 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((0*7)/30.4375) & b5 == 1 & h0 != 8) |>
    mutate(across(starts_with('h0'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h0y, h0m, h0d),

      due_date = (bday %m+% weeks(0)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h0 == 1 ~ 0,  # card with date
        h0 %in% c(2, 3) ~ -1, # reported / card but no date
        h0 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h0d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h0 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h0,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)


}

{
  hep <- d1 |>

    # 14 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((0*7)/30.4375) & b5 == 1 & h50 != 8) |>
    mutate(across(starts_with('h50'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h50y, h50m, h50d),

      due_date = (bday %m+% weeks(0)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h50 == 1 ~ 0,  # card with date
        h50 %in% c(2, 3) ~ -1, # reported / card but no date
        h50 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h50d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h50 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h50,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)


}


# 6-week doses ------------------------------------------------------------

# Penta1 - h51 | Pneumo1 - h54 | Rota1 - h57 | IPV1 - h60 | opv1 - h4
d3 <- d1 |>

  filter(b19 <= 36 & b19 >= floor((6*7)/30.4375) & b5 == 1) |>

  select(caseid, b18,
         h4, h4d, h4m, h4y, h51, h51d, h51m, h51y, h60, h60d, h60m, h60y,
         h57, h57d, h57m, h57y, h54, h54d, h54m, h54y) |>

  mutate(
    birth = as.Date(b18 - 1, origin = "1900-01-01"),
    opv1 = make_date(h4y, h4m, h4d),
    penta1 = make_date(h51y, h51m, h51d),
    ipv1 = make_date(h60y, h60m, h60d),
    rota1 = make_date(h57y, h57m, h57d),
    pneumo1 = make_date(h54y, h54m, h54d)
  ) |>

  select(caseid, h4, h51, h60, h57, h54, birth ,opv1 ,penta1 ,ipv1, rota1, pneumo1)


{
  opv1 <- d1 |>

    # 6 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((6*7)/30.4375) & b5 == 1 & h4 != 8) |>
    mutate(across(starts_with('h4'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h4y, h4m, h4d),

      due_date = (bday %m+% weeks(6)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h4 == 1 ~ 0,  # card with date
        h4 %in% c(2, 3) ~ -1, # reported / card but no date
        h4 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h4d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h4 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h4,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



  }

{
  penta1 <- d1 |>

    # 6 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((6*7)/30.4375) & b5 == 1 & h51 != 8) |>
    mutate(across(starts_with('h51'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h51y, h51m, h51d),

      due_date = (bday %m+% weeks(6)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h51 == 1 ~ 0,  # card with date
        h51 %in% c(2, 3) ~ -1, # reported / card but no date
        h51 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h51d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h51 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h51,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



  }

{
  pneumo1 <- d1 |>

    # 6 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((6*7)/30.4375) & b5 == 1 & h54 != 8) |>
    mutate(across(starts_with('h54'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h54y, h54m, h54d),

      due_date = (bday %m+% weeks(6)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h54 == 1 ~ 0,  # card with date
        h54 %in% c(2, 3) ~ -1, # reported / card but no date
        h54 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h54d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h54 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h54,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



}


# 10-week doses -----------------------------------------------------------

# Penta2 - h52 | Pneumo2 - h55 | Rota2 - h57 | OPV2 - h6
d3 <- d1 |>

  filter(b19 <= 36 & b19 >= floor((10*7)/30.4375) & b5 == 1) |>

  select(caseid, b18,
         h52, h52d, h52m, h52y,
         h55, h55d, h55m, h55y,
         h58, h58d, h58m, h58y,
          h6, h6d, h6m, h6y) |>

  mutate(
    birth = as.Date(b18 - 1, origin = "1900-01-01"),
    penta2 = make_date(h52y, h52m, h52d),
    pneumo2 = make_date(h55y, h55m, h55d),
    rota2 = make_date(h58y, h58m, h58d),
    opv2 = make_date(h6y, h6m, h6d)
  ) |>

  select(caseid, h52, h55, h58, h6, birth, penta2, pneumo2, rota2, opv2)


{
  penta2 <- d1 |>

    # 10 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((10*7)/30.4375) & b5 == 1 & h52 != 8) |>
    mutate(across(starts_with('h52'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h52y, h52m, h52d),

      due_date = (bday %m+% weeks(10)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h52 == 1 ~ 0,  # card with date
        h52 %in% c(2, 3) ~ -1, # reported / card but no date
        h52 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h52d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h52 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h52,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



  }

{
  pneumo2 <- d1 |>

    # 10 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((10*7)/30.4375) & b5 == 1 & h55 != 8) |>
    mutate(across(starts_with('h55'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h55y, h55m, h55d),

      due_date = (bday %m+% weeks(10)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h55 == 1 ~ 0,  # card with date
        h55 %in% c(2, 3) ~ -1, # reported / card but no date
        h55 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h55d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h55 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h55,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



}

{
  opv2 <- d1 |>

    # 10 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((10*7)/30.4375) & b5 == 1 & h6 != 8) |>
    mutate(across(starts_with('h6'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h6y, h6m, h6d),

      due_date = (bday %m+% weeks(10)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h6 == 1 ~ 0,  # card with date
        h6 %in% c(2, 3) ~ -1, # reported / card but no date
        h6 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h6d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h6 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h6,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



}

# 14-week doses -----------------------------------------------------------

# Penta3 - h53 | Pneumo3 - h56 | Rota3 - h59 | OPV3 - h8 | ipv2 - sipv2
d3 <- d1 |>

  filter(b19 <= 36 & b19 >= floor((14*7)/30.4375) & b5 == 1) |>

  select(caseid, b18,
         h53, h53d, h53m, h53y,
         h56, h56d, h56m, h56y,
         h59, h59d, h59m, h59y,
         h8, h8d, h8m, h8y,
         sipv2, sipv2d, sipv2m, sipv2y) |>

  mutate(
    birth = as.Date(b18 - 1, origin = "1900-01-01"),
    penta3 = make_date(h53y, h53m, h53d),
    pneumo3 = make_date(h56y, h56m, h56d),
    rota3 = make_date(h59y, h59m, h59d),
    opv3 = make_date(h8y, h8m, h8d),
    ipv2 = make_date(sipv2y, sipv2m, sipv2d)
  ) |>

  select(caseid, h53 ,h56 ,h59 ,h8 ,sipv2, birth, penta3, pneumo3, rota3, opv3, ipv2)

{
  penta3 <- d1 |>

    # 14 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((14*7)/30.4375) & b5 == 1 & h53 != 8) |>
    mutate(across(starts_with('h53'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h53y, h53m, h53d),

      due_date = (bday %m+% weeks(14)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h53 == 1 ~ 0,  # card with date
        h53 %in% c(2, 3) ~ -1, # reported / card but no date
        h53 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h53d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h53 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h53,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



  }

{
  pneumo3 <- d1 |>

    # 14 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((14*7)/30.4375) & b5 == 1 & h56 != 8) |>
    mutate(across(starts_with('h56'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h56y, h56m, h56d),

      due_date = (bday %m+% weeks(14)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h56 == 1 ~ 0,  # card with date
        h56 %in% c(2, 3) ~ -1, # reported / card but no date
        h56 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h56d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h56 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h56,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



}

{
  opv3 <- d1 |>

    # 14 weeks eligibility.
    filter(b19 <= 36 & b19 >= floor((14*7)/30.4375) & b5 == 1 & h8 != 8) |>
    mutate(across(starts_with('h8'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h8y, h8m, h8d),

      due_date = (bday %m+% weeks(14)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h8 == 1 ~ 0,  # card with date
        h8 %in% c(2, 3) ~ -1, # reported / card but no date
        h8 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h8d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h8 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h8,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



}

# 6-month doses -----------------------------------------------------------

# vit. A - h33
d3 <- d1 |>

  filter(b19 <= 36 & b19 >= 6 & b5 == 1) |>

  select(caseid, b18,
         h33, h33d, h33m, h33y) |>

  mutate(
    birth = as.Date(b18 - 1, origin = "1900-01-01"),
    vit_a = make_date(h33y, h33m, h33d),
  ) |>

  select(caseid, h33, birth, vit_a)

{
  vita1 <- d1 |>

    # 6 months eligibility.
    filter(b19 <= 36 & b19 >= 6 & b5 == 1 & h33 != 8) |>
    mutate(across(starts_with('h33'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h33y, h33m, h33d),

      due_date = (bday %m+% months(6)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h33 == 1 ~ 0,  # card with date
        h33 %in% c(2, 3) ~ -1, # reported / card but no date
        h33 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h33d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h33 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h33,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



  }

# 9-month doses -----------------------------------------------------------

# mcv-1 - h9
d3 <- d1 |>

  filter(b19 <= 36 & b19 >= 9 & b5 == 1) |>

  select(caseid, b18,
         h9, h9d, h9m, h9y) |>

  mutate(
    birth = as.Date(b18 - 1, origin = "1900-01-01"),
    mcv1 = make_date(h9y, h9m, h9d),
  ) |>

  select(caseid, h9, birth, mcv1)

{
  mcv1 <- d1 |>

    # 9 months eligibility.
    filter(b19 <= 36 & b19 >= 9 & b5 == 1 & h9 != 8) |>
    mutate(across(starts_with('h9'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h9y, h9m, h9d),

      due_date = (bday %m+% months(9)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h9 == 1 ~ 0,  # card with date
        h9 %in% c(2, 3) ~ -1, # reported / card but no date
        h9 == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h9d/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h9 == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h9,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)



  }


# 15-month doses ----------------------------------------------------------

# MCV2 - h9a
d3 <- d1 |>
  # b19 >= floor((6*7)/30.4375)
  filter(b19 <= 36 & b19 >= 15 & b5 == 1) |>

  select(caseid, b18,
         h9a, h9ad, h9am, h9ay) |>

  mutate(
    birth = as.Date(b18 - 1, origin = "1900-01-01"),
    mcv2 = make_date(h9ay, h9am, h9ad),
  ) |>

  select(caseid, h9a, birth, mcv2)

{
  mcv2 <- d1 |>

    # 15 months eligibility.
    filter(b19 <= 36 & b19 >= 15 & b5 == 1 & h9a != 8) |>
    mutate(across(starts_with('h9a'), \(x) ifelse(x %in% c(9997:9998, 97:98), NA, x))) |>

    mutate(

      iday = as.Date(v008a - 1, origin = "1900-01-01"),
      bday = as.Date(b18 - 1, origin = "1900-01-01"),
      vday = make_date(h9ay, h9am, h9ad),

      due_date = (bday %m+% months(15)),
      event_time = as.integer(as.numeric(vday - bday)),

      # NA if missing vday
      censor_time = as.integer(as.numeric(iday - bday)),
      time_outcome = ifelse(!is.na(event_time), event_time, censor_time),
      time_outcome = pmax(0, time_outcome),

      # descriptive event label
      event_label = case_when(
        h9a == 1 ~ 0,  # card with date
        h9a %in% c(2, 3) ~ -1, # reported / card but no date
        h9a == 0 ~ 1,   # no vaxx
        T ~ NA
      ), # Don't know (fallback)

      # some people have NA in their h9ad/m/y so force them to be censored.
      # so if vday is NA, event_label is left/right - depends
      event_label = case_when(
        is.na(vday) & h9a == 1 ~ -1, # missing either d/m/y, so make them left censored
        T ~ event_label)
    ) |>

    # dealing with covariates
    mutate(
      v005 = v005 / 1e6,
      m15 = ifelse(m15 %in% c(10, 11, 12, 96), 'home', 'facility'),
      m14 = ifelse(m14 == 98, NA, m14),

      across(c(v717, b4, v501, v131, v130, v106, v190, v024, v025), as_factor)
    ) |>
    select(
      iday, bday, due_date,
      caseid, wt = v005,h9a,vday,
      county = v024, residence = v025, time_outcome, event_label, wealth = v190,
      sex = b4, bord, anc = m14, delivery = m15, curr_age = b8,
      meduc = v106, hhsize = v136, religion = v130, ethnicity = v131,
      curr_marital = v501, occupation = v717
    ) |>
    select(-occupation, -curr_marital, -ethnicity, -religion)

}





# Saving the datasets to memory -------------------------------------------

ldata <- list(
  bcg = bcg, opv0 = opv0, hep = hep,
  penta1 = penta1, pneumo1 = pneumo1, opv1 = opv1,
  penta2 = penta2, pneumo2 = pneumo2, opv2 = opv2,
  penta3 = penta3, pneumo3 = pneumo3, opv3 = opv3,
  mcv1 = mcv1, mcv2 = mcv2, vita1 = vita1
)

saveRDS(ldata, "data/processed/vaxdata-components.rds")

