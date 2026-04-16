

# covariates:
# - sampling dependent: case id, bidx, cluster, admin location, residence type sampling weights
# - place of delivery
# - wealth index
# - mother's education
# - mother's occupation
# - mother's age
# - child's sex
# - Number of ANC Visits
# - Distance to nearest facilities
# - Number of children (birth order of child) / Household size
# - vaccine hesitancy questions (if any)

# Data
d1 <- read_dta("data/dhs/NG_2024_DHS_03262026_919_211396/NGBR8BDT/NGBR8BFL.dta")


d2 <- d1 |>

  # omitting over 36 months old children & alive kids
  filter(b19 <= 36 & b5 == 1) |>

  select(caseid, bidx, wt = v005, cluster = v001, residence = v025, admin = v024,
         birth_order = bord, num_anc_visits = m14, place_delivery = m15,
         child_gender = b4, mother_occupation = v717, time_to_hf = v483a,
         wealth_index = v190, meduc = v106, m_current_age = v012, child_age = b19

  ) |>

  mutate(

    # Sampling weights:
    wt = wt / 1e6,

    # mother's age at birth: current age - child current age
    mother_age_birth = m_current_age - floor(child_age / 12),

    # facility
    place_delivery = ifelse(place_delivery %in% c(10, 11, 12, 90, 96), 'home', 'institution'),

    # changing everything to factor
    across(everything(), as_factor),

    # ANC visits categorization
    num_anc_visits = as.character(num_anc_visits),
    num_anc_visits = case_when(
      num_anc_visits == 'no antenatal visits' ~ '0 visits',
      num_anc_visits %in% as.character(1:4) ~ "1-4 visits",
      num_anc_visits %in% as.character(5:8) ~ "5-8 visits",
      num_anc_visits %in% as.character(9:20) ~ "9+ visits",
      num_anc_visits == "don't know" ~ "Unknown visits",
      T ~ NA
    ),

    # travel time to nearest facility
    time_to_hf = case_when(
      time_to_hf %in% as.character(0:30) ~ "<30 mins",
      time_to_hf %in% as.character(31:60) ~ "31-60 mins",
      time_to_hf %in% as.character(61:120) ~ "1-2 hrs",
      time_to_hf %in% c('600+', as.character(121:600)) ~ "2+ hrs",
      T ~ NA
    ),

    # numeric variables
    across(c(birth_order, mother_age_birth, child_age), \(x) x |> as.character() |> as.numeric())


  ) |>

  select(-m_current_age)

glimpse(d2)

saveRDS(d2, "data/processed/dhs-covariates.rds")
