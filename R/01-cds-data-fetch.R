
require(ecmwfr)


# Simple request for a single year ----------------------------------------

request <- list(
  dataset_short_name = "derived-era5-single-levels-daily-statistics",
  product_type = "reanalysis",
  variable = "2m_temperature",
  year = "2024",
  month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
          "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
          "23", "24", "25", "26", "27", "28", "29", "30", "31"),
  daily_statistic = "daily_mean",
  time_zone = "utc+00:00",
  frequency = "6_hourly",
  area = c(14.5, 2.5, 3.5, 15.5),
  target = '2m_temperature_2024.nc'
)

Sys.sleep(60*2)
file <- wf_request(
  request  = request,
  transfer = TRUE,
  path     = "/data/cds/"
)



# Multiple year requests --------------------------------------------------

dynamic_request <- wf_archetype(
  request = list(
    dataset_short_name = "derived-era5-single-levels-daily-statistics",
    product_type = "reanalysis",
    variable = "2m_temperature",
    year = "1940",
    month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
    day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
            "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
            "23", "24", "25", "26", "27", "28", "29", "30", "31"),
    daily_statistic = "daily_mean",
    time_zone = "utc+00:00",
    frequency = "6_hourly",
    area = c(14.5, 2.5, 3.5, 15.5),
    target = 'tmp_target'
  ),
  dynamic_fields = c("year", "target"))

# setting the day value
batch_request <- list(
  dynamic_request(year = "1970", target = "m2_temp_1970.nc"),
  dynamic_request(year = "1976", target = "m2_temp_1976.nc"),
  dynamic_request(year = "1977", target = "m2_temp_1977.nc"),
  dynamic_request(year = "1978", target = "m2_temp_1978.nc"),
  dynamic_request(year = "1979", target = "m2_temp_1979.nc"),
  dynamic_request(year = "1980", target = "m2_temp_1980.nc"),
  dynamic_request(year = "1997", target = "m2_temp_1997.nc")
)

# cant exceed 20 workers
files <- wf_request_batch(
  batch_request,
  workers = 20,
  path = "data/cds/"
)


# Dew point data fetching -------------------------------------------------

dynamic_request <- wf_archetype(
  request = list(
    dataset_short_name = "derived-era5-single-levels-daily-statistics",
    product_type = "reanalysis",
    variable = "2m_dewpoint_temperature",
    year = "1940",
    month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
    day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
            "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
            "23", "24", "25", "26", "27", "28", "29", "30", "31"),
    daily_statistic = "daily_mean",
    time_zone = "utc+00:00",
    frequency = "6_hourly",
    area = c(14.5, 2.5, 3.5, 15.5),
    target = 'tmp_target'
  ),
  dynamic_fields = c("year", "target"))

# setting the day value
batch_request <- list(
  dynamic_request(year = "2020", target = "m2_dewpoint_temp_2006.nc"),
  dynamic_request(year = "2021", target = "m2_dewpoint_temp_2007.nc"),
  dynamic_request(year = "2022", target = "m2_dewpoint_temp_2008.nc"),
  dynamic_request(year = "2023", target = "m2_dewpoint_temp_2009.nc"),
  dynamic_request(year = "2024", target = "m2_dewpoint_temp_2010.nc")
)

# cant exceed 20 workers
files <- wf_request_batch(
  batch_request,
  workers = 20,
  path = "data/cds/2m_dewpoint_temperature/"
)
