##############################################################################
##  This script contains the source code for a function that                ## 
##  downloads EPA ecoregions and extracts the shapefile.                    ##
##############################################################################

load_ecoregions <- function() {
  if (!dir.exists("data/raw/us_eco_l3/")) {
    download.file("https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip",
                  destfile = "data/raw/us_eco_l3.zip")
    unzip("data/raw/us_eco_l3.zip",
          exdir = "data/raw/us_eco_l3/")
  }
  ecoregion_shp <- terra::vect(x = './data/raw/us_eco_l3', layer = 'us_eco_l3') 
  ecoregion_shp
}
