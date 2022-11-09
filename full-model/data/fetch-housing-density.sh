#!/bin/zsh

# Download data from internet
# wget -O data/raw/us_pbg00.zip http://silvis.forest.wisc.edu/WebData/pbg00_old/gis/us_pbg00.zip
# 7z e data/raw/us_pbg00.zip -odata/raw/us_pbg00.gdb/

# ogr2ogr -f 'ESRI Shapefile' data/processed/us_pbg data/raw/us_pbg00.gdb

# rasterize for each year of interest 1980 - 2030 (4km grid):
SHP="data/processed/us_pbg/us_pbg00_Oct_8_2007.shp"
gdal_rasterize -a HDEN80 -of GTiff -tr 4000 4000 $SHP data/processed/hden80.tif
gdal_rasterize -a HDEN90 -of GTiff -tr 4000 4000 $SHP data/processed/hden90.tif
gdal_rasterize -a HDEN00 -of GTiff -tr 4000 4000 $SHP data/processed/hden00.tif
gdal_rasterize -a HDEN10 -of GTiff -tr 4000 4000 $SHP data/processed/hden10.tif
gdal_rasterize -a HDEN20 -of GTiff -tr 4000 4000 $SHP data/processed/hden20.tif
gdal_rasterize -a HDEN30 -of GTiff -tr 4000 4000 $SHP data/processed/hden30.tif