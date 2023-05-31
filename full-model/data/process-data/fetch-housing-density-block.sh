#!/bin/zsh

# Download data from internet
wget -O raw/conus_wui_blk20.zip https://geoserver.silvis.forest.wisc.edu/geodata/wui_change_2020/zip/fgdb/CONUS_WUI_block_1990_2020_change_fgdb.zip
7z e raw/conus_wui_blk20.zip -oraw/conus_wui_blk20.gdb/

# ogr2ogr -f 'ESRI Shapefile' data/processed/conus_blk20 data/raw/conus_blk20.gdb

# rasterize for each year of interest 1990 - 2020 (4km grid):
# SHP="data/processed/conus_blk20.gpkg"
# gdal_rasterize -l buffered -at -a HUDEN1990 -of GTiff -tr 4000 4000 $SHP data/processed/huden1990.tif
# gdal_rasterize -l buffered -at -a HUDEN2000 -of GTiff -tr 4000 4000 $SHP data/processed/huden2000.tif
# gdal_rasterize -l buffered -at -a HUDEN2010 -of GTiff -tr 4000 4000 $SHP data/processed/huden2010.tif
# gdal_rasterize -l buffered -at -a HUDEN2020 -of GTiff -tr 4000 4000 $SHP data/processed/huden2020.tif