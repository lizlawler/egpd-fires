#!/bin/zsh

# Check boundaries of contiguous US to filter fires from MTBS
# These can change slightly over time, but likely not enough to necessitate downloading most recent 
# tiger files and proceeding with the below process
# However, this process runs quickly (5-10sec total), so if you want a sanity check, run again

# Download data from internet and unzip
wget -O data/raw/tl_2021_us_state.zip https://www2.census.gov/geo/tiger/TIGER2021/STATE/tl_2021_us_state.zip
7z e data/raw/tl_2021_us_state.zip -odata/raw/tl_2021_us_state/

# Create sequence of state FIPS codes to filter by
fips=($(seq -s "," -f "'%02g'" -t "'56'" 1 55))

# Convert to shapefile and limit to lower 48 states + DC 
cd data/raw
ogr2ogr -f 'ESRI Shapefile' -sql "SELECT * FROM tl_2021_us_state WHERE STATEFP IN (${fips:-1}) AND STATEFP NOT IN ('02','15')" lower48 tl_2021_us_state

# Look at info for lower 48; "extent" piece at the top gives long and lat boundaries 
ogrinfo -so lower48/tl_2021_us_state.shp tl_2021_us_state