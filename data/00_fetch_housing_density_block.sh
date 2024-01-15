#!/bin/zsh

##############################################################################
##  This shell script downloads and unzips the CONUS housing density        ##  
## data from the SILVIS lab.                                                ##
##############################################################################

wget -O raw/conus_blk20.zip \
https://geoserver.silvis.forest.wisc.edu/geodata/block-change-2020/zip/fgdb/CONUS_block20_change_1990_2020_PLA4_fgdb.zip
7z e raw/conus_blk20.zip -oraw/conus_blk20.gdb/
