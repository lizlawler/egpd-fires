#!/bin/zsh

# Download data from internet
wget -O data/raw/conus_blk20.zip http://silvis.forest.wisc.edu/GeoData/block_change_2020/zip/fgdb/CONUS_blk20_Census_change_1990_2020_HU_PLA_fgdb.zip
7z e data/raw/conus_blk20.zip -odata/raw/conus_blk20_v2/

# ogr2ogr -f 'ESRI Shapefile' data/processed/conus_blk20 data/raw/conus_blk20.gdb

# rasterize for each year of interest 1990 - 2020 (4km grid):
# SHP="data/processed/conus_blk20.gpkg"
# gdal_rasterize -l buffered -at -a HUDEN1990 -of GTiff -tr 4000 4000 $SHP data/processed/huden1990.tif
# gdal_rasterize -l buffered -at -a HUDEN2000 -of GTiff -tr 4000 4000 $SHP data/processed/huden2000.tif
# gdal_rasterize -l buffered -at -a HUDEN2010 -of GTiff -tr 4000 4000 $SHP data/processed/huden2010.tif
# gdal_rasterize -l buffered -at -a HUDEN2020 -of GTiff -tr 4000 4000 $SHP data/processed/huden2020.tif