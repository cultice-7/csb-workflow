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
               "httr",
               "sf",
               "terra")
pacman::p_load(char = libraries)

#####
### Preparing gSSURGO data
#####

### Merge gSSURGO into MSA files

f_msamerge <- function(msa.i){
  cnty.sf <- tigris::counties(year = 2021)       %>%
    filter(CBSAFP == "16980")                    %>%
    mutate(FIPS_cnty = str_c(STATEFP, COUNTYFP)) %>%
    st_transform(., crs = "ESRI:102008")
  
  st.v <- unique(cnty.sf$STATEFP)
  stabbv.v <- as.data.table(fips_codes) %>%
    .[state_code %in% st.v]             %>%
    .[,c("state")]                      %>%
    unique()
    
  files.v  <- list.files(here("Data", "SSURGO"),
                         pattern = paste0(stabbv.v$state, collapse = "|")) %>%
    .[str_detect(., ".tif$")]
  ssurgo.l <- sprc(map(here("Data", "SSURGO", files.v), terra::rast))
  ssurgo.r <- merge(ssurgo.l)
  
  ### Merge together raw data (e.g. mapunits and walks)
  
  
}
### Grab mapunit datasets
f_mapunitmerge <- function(st.i){
  ### Note: Load layers of interest from 
  # ssurgo.layers <- st_layers(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")))$name
  
  layer1.dt     <- st_read(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")),
                           layer = "Valu1") %>%
    as.data.table()                         %>%
    .[,c("nccpi3all", "mukey")]
  
  component.dt  <- st_read(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")),
                           layer = "component") %>%
    as.data.table()                             %>%
    .[,c("cokey", "comppct_r", "mukey")]
  
  horizon.dt <- st_read(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")),
                        layer = "chorizon") %>%
    as.data.table()                         %>%
    .[,c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact", "cokey", "chkey")] %>%
    .[,kwfact := as.numeric(kwfact)]                                               %>%
    .[,(c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact")) := lapply(.SD, mean, na.rm = TRUE),
      by = c("cokey"),
      .SDcols = c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact")]          %>%
    .[,!c("chkey")]                                                                %>%
    unique()
    
  ssurgo_walk.dt <- copy(horizon.dt) %>%
    merge(copy(component.dt),
          by    = "cokey",
          all.x = TRUE)              %>%
    merge(copy(layer1.dt),
          by    = "mukey",
          all.x = TRUE)              %>%
    .[!is.na(comppct_r)]             %>%
    .[,comppct_r := comppct_r/100]   %>%
    .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact"))) := lapply(.SD, function(x) comppct_r*x),
      .SDcols = c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact")] %>%
    .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact"))) := round(.SD/sum(comppct_r),2), 
      by = c("mukey"),
      .SDcols = str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact"))]
}