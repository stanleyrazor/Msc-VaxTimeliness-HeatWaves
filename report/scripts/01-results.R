

# Table 1 -----------------------------------------------------------------

pacman::p_load(gtsummary, dplyr)

# data
antigen <- c('bcg', 'penta1', 'penta2', 'penta3', 'vita1', 'mcv1')
usenames <- c('BCG', 'Pentavalent 1', 'Pentavalent 2', 'Pentavalent 3', 'Vitamin A', 'MCV1')

ldata <- readRDS('data/processed/vaxdata-components.rds')
cdata <- readRDS('data/processed/dhs-covariates.rds')
use_ldata <- ldata[antigen]

for (i in 1:length(use_ldata)) use_ldata[[i]] <- use_ldata[[i]] |> mutate(antigen = usenames[i])

d1 <- merge(bind_rows(use_ldata) |> mutate(across(c(wt, caseid, bidx), as.character))
            , cdata |> mutate(wt = as.character(wt)) |> select(-num_anc_visits),
            by = c('caseid', 'bidx', 'admin', 'cluster', 'wt', 'residence'), all.x = T)


d2 <- d1 |>
  transmute(

    Antigen = factor(antigen, labels = usenames, levels = usenames),

    Residence = factor(
      residence,
      levels = c("urban", "rural"),
      labels = c("Urban", "Rural")
    ),

    `Birth order` = birth_order,

    `Place of delivery` = factor(
      place_delivery,
      levels = c("home", "institution"),
      labels = c("Home", "Health facility")
    ),

    `Child's gender` = factor(
      child_gender,
      levels = c("male", "female"),
      labels = c("Male", "Female")
    ),

    `Travel time to the\nnearest health facility` = factor(
      time_to_hf,
      levels = c("<30 mins", "31-60 mins", "1-2 hrs", "2+ hrs"),
      labels = c(
        "<30 minutes",
        "31–60 minutes",
        "1–2 hours",
        ">2 hours"
      )
    ),

    `Wealth index` = factor(
      wealth_index,
      levels = c("poorest", "poorer", "middle", "richer", "richest"),
      labels = c("Poorest", "Poorer", "Middle", "Richer", "Richest")
    ),

    `Highest level of\neducation attained by mother` = factor(
      meduc,
      levels = c("no education", "primary", "secondary", "higher"),
      labels = c(
        "No formal education",
        "Primary school",
        "Secondary school",
        "Higher education"
      )
    ),

    `Mother's age group\nat birth` = case_when(
      mother_age_birth <= 19 ~ "<=19 yrs",
      mother_age_birth >= 20 & mother_age_birth <= 24 ~ "20-24 yrs",
      mother_age_birth >= 25 & mother_age_birth <= 29 ~ "25-29 yrs",
      mother_age_birth >= 30 & mother_age_birth <= 34 ~ "30-34 yrs",
      mother_age_birth >= 35 & mother_age_birth <= 39 ~ "35-39 yrs",
      mother_age_birth >= 40 & mother_age_birth <= 44 ~ "40-44 yrs",
      mother_age_birth >= 45 ~ "45+ yrs",
      TRUE ~ NA_character_
    ),
    `Mother's age group\nat birth` = factor(
      `Mother's age group\nat birth`,
      levels = c("<=19 yrs", "20-24 yrs","25-29 yrs","30-34 yrs",
                 "35-39 yrs","40-44 yrs","45+ yrs"),
      labels = c("≤19 years","20–24 years","25–29 years","30–34 years",
                 "35–39 years","40–44 years","45 years or older")
    ),

    `Geopolitical zone` = factor(geozone,
                                 levels = c("north west", "north east", "north central", "south east",
                                            "south south", "south west"),
                                 labels = c("North West", "North East", "North Central", "South East",
                                            "South South", "South West"))
  )

d2 |>
  tbl_summary(
    by = Antigen,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      `Birth order` ~ "{mean} ({sd})"
    ),
    missing = "no",
    label = list(
      Residence ~ "Residence",
      `Birth order` ~ "Birth order",
      `Place of delivery` ~ "Place of delivery",
      `Child's gender` ~ "Child's gender",
      `Travel time to the\nnearest health facility` ~ "Travel time to nearest health facility",
      `Wealth index` ~ "Wealth index",
      `Highest level of\neducation attained by mother` ~ "Maternal education",
      `Mother's age group\nat birth` ~ "Mother's age at birth",
      `Geopolitical zone` ~ "Geopolitical zone"
    )
  ) |>
  bold_labels() |>
  modify_header(label ~ "**Characteristic**") |>
  modify_spanning_header(all_stat_cols() ~ "**Antigen**")



