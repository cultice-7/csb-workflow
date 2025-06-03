libraries <- c("tidyverse",
               "tidytext",
               "tictoc",
               "Hmisc",
               "data.table",
               "here",
               "future", 
               "future.apply",
               "furrr",
               "pryr",
               "tigris",
               "tidycensus",
               "sf",
               "rvest",
               "R.utils",
               "lubridate",
               "matrixStats",
               "terra",
               "tigris",
               "haven",
               "rvest",
               "jsonlite",
               "httr")
pacman::p_load(char = libraries)

#####
### Set up FIPS code
#####

dt.fips      <- fips_codes                           %>%
  as.data.table()                                    %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]    %>%
  .[state_code %in% c("17", "18", "26", "39", "55")]

#####
### Yearly intersection with NLCD and HUC8-Filtered for Watersheds
#####

f.huc_nlcd <- function(){
  v.years  <- c("2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019", "2021")
  v.watersheds <- c("Lower Maumee", "Sugar", "Upper Fox", "Macoupin", "Maple")
  sf.huc <- st_read(here("Data", "USGS", "WBD", "USGS_HUC8.gpkg")) %>%
    filter(name %in% v.watersheds) %>%
    .[,c("name", "geom")]
  
  f.nlcd_year <- function(i.year,
                          i.sf.huc = sf.huc){
    ### Load GL NLCD
    rast.nlcd = terra::rast(here("Data", "NLCD", "Raw", paste0("NLCD_GreatLakes_", i.year, ".tif")))
    
    ### Obtain levels and reclassify
    vals = levels(rast.nlcd)[[1]]
    
    ### convert crs in huc file and create sv
    sv.huc <- st_transform(i.sf.huc,
                           crs = terra::crs(rast.nlcd)) %>%
      terra::vect(.)
    
    ### Load data from huc file
    n              <- dim(sv.huc)[1]
    dt.sv.huc   <- copy(sv.huc) %>%
      as.data.table()           %>%
      .[,ID     := seq.int(n)]
    
    ### Intersection
    dt.extract    <- terra::extract(rast.nlcd, sv.huc, 
                                    touches = FALSE,
                                    cells   = TRUE)  %>%
      as.data.table(key = c("ID", "NLCD Land Cover Class")) %>%
      merge(dt.sv.huc, 
            all.x = TRUE,
            by = "ID") %>%
      .[,!c("ID")] %>%
      .[,Area_m2 := .N * 900, by = c("NLCD Land Cover Class", "name")] %>%
      .[,c("Area_m2", "NLCD Land Cover Class", "name")] %>%
      .[,Year := i.year] %>%
      merge(vals, 
            by = "NLCD Land Cover Class",
            all.x = TRUE) %>%
      setnames(c("name", "value"),
               c("huc8_name", "NLCD_value")) %>%
      unique()
    
    ### Return raw pixel intersections
    return(dt.extract)
  }
  dt.full <- map(v.years, f.nlcd_year) %>% rbindlist() %>%
    dcast(huc8_name + Year ~ `NLCD Land Cover Class`, value.var = c("Area_m2")) %>%
    .[,Unit := "m^2"]
  fwrite(dt.full,
         here("Data", "NLCD", "NLCD_HUC8_Aggregated.csv"))
}
f.huc_nlcd()