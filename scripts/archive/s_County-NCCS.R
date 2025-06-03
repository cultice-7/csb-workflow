################################################################################
### INFEWS/NIFA: NCCS-County Intersections
### Description: This script extracts daily weather data for years 2000-2050 from 
###              the NCCS CMIP6 ACCESS-CM2 collection. The scripts then process 
###              the resulting data to create a county average for each day
###              which is used to generate a set of county-level predictors for
###              our simulations
################################################################################

libraries <- c("tidyverse",
               "tictoc",
               "tibble",
               "data.table",
               "terra", # Raster
               "sf",    # Vector (shapefiles)
               "here",  
               "future", 
               "future.apply",
               "raster", # Raster
               "furrr", 
               "pryr",
               "rmapshaper",
               "tigris")
pacman::p_load(char = libraries)

### Point to weather data and capture CSB
weat.loc  <- "K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\NCCS"
nccs.crs  <- terra::rast(paste0(weat.loc, "\\sfcWind\\sfcWind_day_ACCESS-CM2_ssp245_r1i1p1f1_gn_2015.nc")) %>%
  terra::crs(.)

### Load county shapefile
cnty.sf   <- tigris::counties(state = NULL,
                              cb    = TRUE,
                              year  = NULL)            %>%
  filter(STATEFP %in% c("17", "18", "26", "39", "55")) %>%
  st_transform(., crs = nccs.crs)                      %>%
  mutate(FIPS_st   = STATEFP,
         FIPS_cnty = str_c(STATEFP, COUNTYFP),
         county_name = NAME,
         state_name  = STATE_NAME)                     %>%
  .[,c("FIPS_st", "FIPS_cnty", "county_name", "state_name",
       "geometry")]

f_county_nccs_int <- function(yr.i,
                              data.i,
                              scenario.i,
                              cnty.sf.i  = cnty.sf,
                              weat.loc.i = weat.loc){
  
  ### Step 1: Process county vector
  cnty.sv.i <- copy(cnty.sf.i) %>%
    vect()
    
  cnty_sv.dt.i <- as.data.table(cnty.sf.i) %>%
    .[,ID     := .I]                       %>%
    .[,c("ID", "FIPS_cnty")]
    
  ### Step 2: Load rasters and filter the desired raster  
  rast.loc <- paste0(weat.loc.i,
                     "\\",
                     data.i)
  rast.v   <- list.files(rast.loc,
                         pattern = ".nc$") %>%
    as.data.table()                         %>%
    .[,Year := str_extract(.,
                           pattern = "2[:digit:]{3}")] %>%
    .[,Scenario := str_extract(.,
                               pattern = "historical|ssp245")] %>%
    setorder(Year)                                             %>%
    .[Scenario == scenario.i]                                 %>%
    .[Year == yr.i]
  
  ### Step 3: Load the world shapefile and change extent to reflect western hem
  ext.r    <- ext(180, 360, 0, 90)
  rast.r   <- terra::rast(str_c(rast.loc, "\\", rast.v$.)) %>%
    terra::crop(., ext.r)
  set.ext(rast.r, c(-180, 0, 0, 90))
  
  ### Step 4: Extract weather for each county
  rast.ext <- terra::extract(rast.r, cnty.sv.i,
                             ID = TRUE)            %>%
    as.data.table()                                %>%
    merge.data.table(cnty_sv.dt.i,
                     by = "ID")                    %>%
    .[,Weather := data.i]                          %>%
    .[,!("ID")]                                    %>%
    .[,Year := yr.i]
  
  ### Step 5: Replace w/ average values
  vars.v <- names(rast.ext)[!(names(rast.ext) %in% c("FIPS_cnty", "Weather", "Year"))]
  rast.ext <- rast.ext %>%
    .[,(vars.v) := lapply(.SD, mean, na.rm = TRUE), by = "FIPS_cnty",
      .SDcols = vars.v] %>%
    unique()
  
  f_createdate <- function(x){
    num <- str_extract_all(x,
                           "[:digit:]")[[1]] %>%
      paste0(., collapse = "")               %>%
      as.numeric() - 1
    
    date.start <- lubridate::mdy(paste0("01-01-", yr.i))
    date.new   <- str_c("Date_", str_remove_all(as.character(date.start + num),
                                                "-"))
    return(date.new)
  }
    colA = names(rast.ext)[str_detect(names(rast.ext),data.i)]
    colB = map(colA, f_createdate) %>% unlist()
    rast.ext <- setnames(rast.ext,
                         old = colA,
                         new = colB)        %>%
      melt(., measure = list(colB), 
           value.name     = c("Value"),
           variable.name  = c("Date"))      %>%
      .[,Date := str_remove(Date, "Date_")] %>%
      setorder(Date, FIPS_cnty)
    fwrite(rast.ext,
           paste0("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\County\\", "County-NCCS_", scenario.i, "_", data.i, "_", yr.i, ".csv"))
    return("done")
}

### Run county intersections for historical data
inp1 <- rep(2000:2014, each = 6)
inp2 <- rep(c("huss", "pr", "sfcWind", "tas", "tasmax", "tasmin"), times = length(2000:2014))
inp3 <- rep("historical", times = length(inp1))
cl          <- plan(multisession, workers = 3)
demint.l    <- future_pmap(list(inp1, inp2, inp3),
                           f_county_nccs_int,
                           .options = furrr_options(scheduling = FALSE))
future:::ClusterRegistry("stop")

### Run county intersections for scenarios of interest
inp1 <- rep(2015:2050, each = 6)
inp2 <- rep(c("huss", "pr", "sfcWind", "tas", "tasmax", "tasmin"), times = length(2015:2050))
inp3 <- rep("ssp245", times = length(inp1))
cl          <- plan(multisession, workers = 3)
demint.l    <- future_pmap(list(inp1, inp2, inp3),
                           f_county_nccs_int,
                           .options = furrr_options(scheduling = FALSE))
future:::ClusterRegistry("stop")

################################################################################
### Create variables of interest for each county for use in BMP modeling
################################################################################
cnty_weat_files.dt <- list.files("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\County\\",
                                 full.names = TRUE)                                      %>%
  as.data.table()                                                                        %>%
  .[,Year     := as.numeric(str_extract(., "[:digit:]{4}"))]                             %>%
  .[,Scenario := str_extract(.,
                             "historical|ssp245")]
cnty_weat.dt <- lapply(cnty_weat_files.dt$., fread) %>% rbindlist()

cnty_rain.dt <- copy(cnty_weat.dt)  %>%
  .[Weather == "pr"]                %>%
  .[, Value := Value*60*60*24*(1/25.4)]  %>% # Conversion of mm/s to inches per day
  .[, Date  := lubridate::ymd(as.character(Date))] %>%
  .[, Year  := lubridate::year(Date)]              %>%
  .[, Month := lubridate::month(Date)]             %>%
  .[, Day   := lubridate::day(Date)]

cnty_rain_avg.dt <- copy(cnty_rain.dt) %>%
  .[Year <= 2020]                      %>%
  .[Month >= 4 & Month <= 10, GrowingSeason_RainTotal := sum(Value), by = c("Year", "FIPS_cnty")] %>%
  .[,c("FIPS_cnty", "Year", "GrowingSeason_RainTotal")]                                           %>%
  .[!is.na(GrowingSeason_RainTotal)]                                                              %>%
  unique()                                                                                        %>%
  .[,GrowingSeason_RainAvg := mean(GrowingSeason_RainTotal, na.rm = TRUE), by = c("FIPS_cnty")]   %>%
  .[,!c("Year", "GrowingSeason_RainTotal")]                                                               %>%
  unique()

cnty_rain_shock.dt <- copy(cnty_rain.dt)                                                          %>%
  .[Year > 2020]                                                                                  %>%
  .[Month >= 4 & Month <= 10, GrowingSeason_RainTotal := sum(Value), by = c("Year", "FIPS_cnty")] %>%
  .[,c("FIPS_cnty", "Year", "GrowingSeason_RainTotal")]                                           %>%
  .[!is.na(GrowingSeason_RainTotal)]                                                              %>%
  unique()                                                                                        %>%
  merge(cnty_rain_avg.dt,
        by = "FIPS_cnty")                                                                         %>%
  .[,GrowingSeason_RainfallShock_in]
  

# cnty_files.v <- list.files("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\County\\") %>%
#   as.data.table()                                                                  %>%
#   .[,Year := as.numeric(str_extract(., "[:digit:]{4}"))]                           %>%
#   .[Year >= 2015]                                                                  %>%
#   .[,newname := str_replace(.,
#                             "historical",
#                             "ssp245")]
# file.rename(paste0("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\County\\", cnty_files.v$.),
#             paste0("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\County\\", cnty_files.v$newname))  
