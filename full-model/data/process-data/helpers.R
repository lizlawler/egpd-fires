
# Helper functions for fire extremes project ------------------------------

load_ecoregions <- function() {
  if (!dir.exists("full-model/data/raw/us_eco_l3/")) {
    download.file("https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip",
                  destfile = "full-model/data/raw/us_eco_l3.zip")
    unzip("full-model/data/raw/us_eco_l3.zip",
          exdir = "full-model/data/raw/us_eco_l3/")
  }
  ecoregion_shp <- vect(x = './full-model/data/raw/us_eco_l3', layer = 'us_eco_l3') 
  ecoregion_shp
}
