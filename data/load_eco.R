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
  # fix names for Chihuahuan Deserts (L3 ecoregion) 
  ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                    'Chihuahuan Deserts',
                                    ecoregion_shp$NA_L3NAME)
  # fix names for Upper Gila Mountains (L2 ecoregion) 
  ecoregion_shp$NA_L2NAME <- ifelse(ecoregion_shp$NA_L2NAME == 'UPPER GILA MOUNTAINS (?)',
                                    'UPPER GILA MOUNTAINS',
                                    ecoregion_shp$NA_L2NAME)
  ecoregion_shp
}
