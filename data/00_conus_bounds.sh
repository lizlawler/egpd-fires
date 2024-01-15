#!/bin/zsh

##############################################################################
##  This shell script checks the boundaries of the CONUS to filter          ##
##  fires from the MTBS dataset. These coordinates can change slightly      ##
##  over time, but likely not enough to necessitate downloading the most    ##
##  recent tiger files and proceeding with the below process. However, the  ##
##  code runs quickly (5-10sec), so can easily be run again.                ##
##############################################################################

# Download data from internet and unzip
wget -O raw/tl_2021_us_state.zip https://www2.census.gov/geo/tiger/TIGER2021/STATE/tl_2021_us_state.zip
7z e raw/tl_2021_us_state.zip -oraw/tl_2021_us_state/

# Create sequence of state FIPS codes to filter by
fips=($(seq -s "," -f "'%02g'" -t "'56'" 1 55))

# Convert to shapefile and limit to lower 48 states & DC 
cd raw
ogr2ogr -f 'ESRI Shapefile' -sql "SELECT * FROM tl_2021_us_state WHERE STATEFP IN (${fips:-1}) AND STATEFP NOT IN ('02','15')" lower48 tl_2021_us_state

# Look at info for lower 48; "extent" information at the top gives longitude and latitude boundaries 
ogrinfo -so lower48/tl_2021_us_state.shp tl_2021_us_state