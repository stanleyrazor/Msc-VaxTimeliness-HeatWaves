
pacman::p_load(terra, geodata)
# x <- rast("~/Downloads/data.grib")
x <- rast("~/Downloads/7b282c777e1ea581d4dec1ae61957aea.nc")
x

# 1. Define the day you want
target_day <- as.Date("2020-06-01")

# 2. Subset layers where the date matches
# x[[...]] is the standard way to subset layers in terra
day_raster <- x[[ as.Date(time(x)) == target_day ]]
day_raster <- project(day_raster, "EPSG:4326")
day_raster <- day_raster[[1]]
# 3. Plot the 24 hourly layers
plot(day_raster)


# GADM Map of Nigeria Admin 0: --------------------------------------------

adm0 <- gadm('Nigeria', level = 0, path = 'data/shp/')
plot(adm0)


day_nigeria_crop <- crop(day_raster, adm0)
day_nigeria_final <- mask(day_nigeria_crop, adm0)

plot(day_raster, main="Temperature in Nigeria (EPSG:4326)")
points(adm0, pch = '.')



# -------------------------------------------------------------------------

# from: ERA5 post-processed daily statistics on single levels from 1940 to present

x <- rast("~/Downloads/7b282c777e1ea581d4dec1ae61957aea.nc")
x <- rast("~/Downloads/4a0c726aff1a3a3e6974730809b76687.nc")
