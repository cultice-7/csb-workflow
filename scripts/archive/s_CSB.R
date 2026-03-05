######################################################################
##### CSB - Create 
##### Description: Process CSB Boundaries
#####   - 
##### Date: 9/12/2023
######################################################################

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

################################################################################
### Note: Create Grid to intersect with fields
###   - Due to the complexity of the field files, we don't want to conduct the 
###     intersections between the CSB boundaries and other spatial data. Instead
###     we create a grid of 30m x 30m centroids which we first assign to given
###     field polygons (for 2022) and then intersect the relevant centroids with
###     soil, weather, and terrain data
###   - Possibly, these centroids can act as stand-ins for "fields", which we then
###     predict the probability of them changing over time using the land use change
###     model
################################################################################
f_gridcentroids <- function(i){
  
  #####
  ### Step 1: Load county shapefile
  #####
  
  cnty.sf <- tigris::counties(year = 2021, cb = TRUE)    %>%
    filter(STATEFP %in% c("17", "18", "26", "39", "55")) %>%
    mutate(FIPS_cnty  = str_c(STATEFP, COUNTYFP))         %>%
    mutate(FIPS_st    = STATEFP)                          %>%
    mutate(CountyName = NAME)                             %>%
    .[,c("FIPS_cnty", "FIPS_st", "CountyName", "geometry")] %>%
    st_transform(., crs = "ESRI:102008")
  
  #####
  ### Step 2: Generate grid centroid based on county bounding box
  #####
  cnty_bb.sf <- st_bbox(cnty.sf)
  
  # Create mesh of centroids; create coords variables
  grid.sf <- st_make_grid(cnty_bb.sf,
                          cellsize = 300,
                          what     = "centers") %>%
    st_as_sf()
  grid.x <- unlist(purrr::transpose(grid.sf$x)[[1]])
  grid.y <- unlist(purrr::transpose(grid.sf$x)[[2]])
  grid.sf <- grid.sf %>%
    mutate(lon_proj = grid.x,
           lat_proj = grid.y) %>%
    mutate(GridID = str_c("G_", 1:dim(.)[1]))
  
  # Create focal mesh of gridcells
  focal.sf <- st_make_grid(cnty_bb.sf,
                           cellsize = 3000,
                           what     = "polygons") %>%
    st_as_sf()                                    %>%
    mutate(FocalID = str_c("F_", 1:dim(.)[1]))
  
  #####
  ### Step 3: Basic intersections between centroids and county sf; create lat/long
  #####
  grid.sf <- st_intersection(grid.sf, cnty.sf) %>%
    st_intersection(.,focal.sf)                %>%
    st_transform(., 4326)
  grid.x <- unlist(purrr::transpose(grid.sf$x)[[1]])
  grid.y <- unlist(purrr::transpose(grid.sf$x)[[2]])
  
  grid.sf <- grid.sf     %>%
    mutate(lon = grid.x,
           lat = grid.y) %>%
    st_transform(., crs = "ESRI:102008")
  st_write(grid.sf,
           here("Data", "Gridded", "GL_Gridded_300m-3km.gpkg"),
           append=FALSE)
  return(print("Gridded Saved"))
}
map(1, f_gridcentroids)

################################################################################
### Creating CSB Datasets
################################################################################

f_csb_stsave <- function(statefips.i){
  
  #####
  ### Step 1: Filter states for csb boundaries
  #####
  q.v <- paste0("SELECT * FROM nationalGIS where STATEFIPS = '", 
                as.character(statefips.i),
                "'")
  
  csb.sf <- st_read(here("Data", "USDA", "CSB2022", "CSB1522.gdb"),
                    layer = "nationalGIS",
                    query = q.v)
  
  #####
  ### Step 2: Filter for each county
  #####
  fips.dt <- tigris::fips_codes                        %>%
    as.data.table()                                    %>%
    .[state_code == statefips.i]                       %>%
    .[,FIPS_cnty := str_c(state_code, county_code)]
  
  f_csb_cntysave <- function(FIPS_cnty.i,
                             csb.sf.i = csb.sf){
    
    csb.sf.i <- csb.sf.i %>%
      filter(CNTYFIPS == str_sub(FIPS_cnty.i,
                                 start = 3))
      
    #####
    ### Step 2: Filter for outlier fields wrt shape; create centroids and save latlon
    #####
    if (dim(csb.sf.i)[1] != 0){
      csb.sf.i <- csb.sf.i                                          %>%
        mutate(PerAr_Ratio = round(Shape_Length/Shape_Area, 6))     %>%
        filter(PerAr_Ratio <= quantile(PerAr_Ratio, 0.95))          %>%
        st_transform(., crs = "ESRI:102008")
      
      csb_cen.sf <- st_centroid(csb.sf.i) %>%
        st_transform(., 4326)
      cen.x <- unlist(purrr::transpose(csb_cen.sf$Shape)[[1]])
      cen.y <- unlist(purrr::transpose(csb_cen.sf$Shape)[[2]])
      
      csb_cen.sf <- csb_cen.sf               %>%
        mutate(lon = cen.x,
               lat = cen.y)                  %>%
        st_transform(., crs = "ESRI:102008")
      
      #####
      ### Step 3: Write shapefiles of full field and centroid for all counties
      #####
      st_write(csb.sf.i,
               here("Data", "USDA", "CSB", paste0("CSB_", FIPS_cnty.i, "_2022.gpkg")),
               append = FALSE)
      st_write(csb_cen.sf,
               here("Data", "USDA", "CSB", paste0("CSB_", FIPS_cnty.i, "_2022_centroid.gpkg")),
               append = FALSE)
      return(print(paste0("County ", FIPS_cnty.i, " is present and saved")))
    }else{
      return(print(paste0("County ", FIPS_cnty.i, " doesn't exist in the data")))
    }
  }
  cl          <- plan(multisession, workers = 4)
  future_lapply(fips.dt$FIPS_cnty, 
                f_csb_cntysave,
                future.scheduling = FALSE)
  future:::ClusterRegistry("stop")
  
  return(print(paste0("State ", statefips.i, " completed")))
}
map(c("18", "26", "39", "55"),
    f_csb_stsave)

################################################################################
### Intersecting CSB w/ other data: SSURGO and DEM
################################################################################

cnty.dt <- fips_codes %>%
  as.data.table()     %>%
  .[state_code %in% c("17", "18", "26", "39", "55")] %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]

fips.v_error <- cnty.dt$FIPS_cnty[which(cnty.dt$FIPS_cnty=="55099"):length(cnty.dt$FIPS_cnty)]

f_csb_spatial <- function(cntyfips.i){
  
  #####
  ### Step 0: Set up crosswalks and other items
  #####
  stname.dt <- fips_codes        %>%
    as.data.table()              %>%
    .[,c("state", "state_code")] %>%
    unique()
  stname.v  <- stname.dt[state_code == str_sub(cntyfips.i,
                                               end = 2)]$state
  
  #####
  ### Step 1: Load csb centroids file; rename; recode rotations
  #####
  
  csb_cen.sf <- st_read(here("Data", "USDA", "CSB", paste0("CSB_", cntyfips.i, "_2022_centroid.gpkg"))) %>%
    rename(FIPS_st = STATEFIPS, 
           `FIPS_st-ASD` = STATEASD,
           USDA_ASD = ASD,
           Name_cnty = CNTY,
           FIPS_cnty = CNTYFIPS)                            %>%
    dplyr::select(!c("INSIDE_X", "INSIDE_Y", "Shape_Leng"))

  #####
  ### Step 2: Terrain intersections (Slope and aspect intersections)
  #####
  dem_tags.dt <- fread(here("Data", "DEM", "Raw", "DEM_tags.csv"))
  files.v    <- list.files(here("Data", "DEM", "Raw"))
  filesize.v <- file.size(here("Data", "DEM", "Raw", files.v))
  dem_dl.dt  <- data.table("SaveName" = files.v,
                           "FileSize" = filesize.v) %>%
    merge(dem_tags.dt,
          by = "SaveName")                          %>%
    .[FileSize > 10000]                             %>%
    .[Date == MostRecent]                           %>%
    .[FIPS_cnty == cntyfips.i]                      %>%
    .[,s_SaveName := str_c("s_", SaveName)]         %>%
    .[,a_SaveName := str_c("a_", SaveName)]
  
  f_demint   <- function(i,
                         dem_dl.dt.i  = dem_dl.dt,
                         csb_cen.sf.i = csb_cen.sf){
    elev.r   <- terra::rast(here("Data", "DEM", "Raw", unique(dem_dl.dt.i[i]$SaveName))) %>%
      terra::project(st_crs(csb_cen.sf.i)$proj4string)
    slope.r  <- terra::rast(here("Data", "DEM", "Raw", unique(dem_dl.dt.i[i]$s_SaveName))) %>%
      terra::project(st_crs(csb_cen.sf.i)$proj4string)
    asp.r    <- terra::rast(here("Data", "DEM", "Raw", unique(dem_dl.dt.i[i]$a_SaveName))) %>%
      terra::project(st_crs(csb_cen.sf.i)$proj4string)
    
    csb_cen.sf.i <- csb_cen.sf.i %>%
      mutate(ID = 1:dim(.)[1])
    csb_cen.sv   <- copy(csb_cen.sf.i) %>%
      terra::vect(.)
    st_geometry(csb_cen.sf.i) <- NULL
    csb_cen.dt.i <- copy(csb_cen.sf.i) %>%
      as.data.table()

    ### Step 2: Intersect Buffers w/ NLCD
    elev.dt <- terra::extract(elev.r, csb_cen.sv,
                               ID      = TRUE)   %>%
      as.data.table()                            %>%
      setnames(c("ID", "Elev"))                  %>%
      .[!is.nan(Elev)]                           %>%
      merge(csb_cen.dt.i,
            by = "ID")                           %>%
      .[,Elev   := round(Elev, digits = 1)]      %>%
      .[,c("CSBID", "Elev")]                     %>%
      unique()                                   %>%
      .[,TileDEM := dem_dl.dt.i[i]$AreaID]
    
    slope.dt <- terra::extract(slope.r, csb_cen.sv,
                               ID      = TRUE)   %>%
      as.data.table()                          %>%
      .[!is.nan(slope)]                        %>%
      merge(csb_cen.dt.i,
            by = "ID")                         %>%
      .[,slope   := round(slope, digits = 1)]  %>%
      .[,c("CSBID", "slope")]                  %>%
      unique()                                 %>%
      .[,TileDEM := dem_dl.dt.i[i]$AreaID]
    
    asp.dt   <- terra::extract(asp.r, csb_cen.sv,
                               ID      = TRUE) %>%
      as.data.table()                          %>%
      .[!is.nan(aspect)]                       %>%
      merge(csb_cen.dt.i,
            by = "ID")                         %>%
      .[,aspect   := round(aspect, digits = 1)]%>%
      .[,c("CSBID", "aspect")]                 %>%
      unique()                                 %>%
      .[,TileDEM := dem_dl.dt.i[i]$AreaID]
      return(list(elev.dt, slope.dt, asp.dt))
  }
  cl          <- plan(multisession, workers = 4)
  demint.l    <- future_lapply(1:dim(dem_dl.dt)[1], f_demint,
                               future.seed = TRUE)
  future:::ClusterRegistry("stop")
  
  ### Pull slope and aspect from the results
  f_dempull   <- function(i,
                          var){
    dt <- demint.l[[i]][[var]]
  }
  elevint.dt  <- map2(1:length(demint.l),
                      rep(1, times = length(demint.l)),
                      f_dempull) %>% rbindlist() %>%
    .[, Elev := round(mean(Elev), 1), by = "CSBID"] %>%
    .[,!c("TileDEM")]                                   %>%
    unique()
  
  slpint.dt   <- map2(1:length(demint.l),
                      rep(2, times = length(demint.l)), 
                      f_dempull) %>% rbindlist() %>%
    .[, slope := round(mean(slope), 1), by = "CSBID"] %>%
    .[,!c("TileDEM")]                                 %>%
    unique()
  
  aspint.dt   <- map2(1:length(demint.l),
                      rep(3, times = length(demint.l)),
                      f_dempull) %>% rbindlist() %>%
    .[, aspect := round(mean(aspect), 1), by = "CSBID"] %>%
    .[,!c("TileDEM")]                                   %>%
    unique()
  
  
  #####
  ### Step 2: gSSURGO Intersection
  #####
  ssurgo.v <- list.files(here("Data", "SSURGO"),
                         pattern = stname.v) %>%
    .[str_detect(., ".tif$")]
  ssurgo.r <- terra::rast(here("Data", "SSURGO", ssurgo.v))
  csb_cen_ssurgo.sf <- copy(csb_cen.sf) %>%
    st_transform(terra::crs(ssurgo.r))
  
  csb_cen_ssurgo.sv <- copy(csb_cen_ssurgo.sf) %>%
    vect(.)
  
  n          <- dim(csb_cen_ssurgo.sv)[1]
  csb_cen_ssurgo_sv.dt <- csb_cen_ssurgo.sv@ptr$getDF() %>%
    as.data.table()                                     %>%
    .[,ID     := seq.int(n)]                            %>%
    .[,c("ID", "CSBID")]
  
  ssurgo.ext <- terra::extract(ssurgo.r, csb_cen_ssurgo.sv,
                               ID      = TRUE)          %>%
    as.data.table(key = c("ID"))                        %>%
    merge(csb_cen_ssurgo_sv.dt,
          by = "ID")                                    %>%
    .[,!c("ID")]
  
  #####
  ### Step 3: Merge intersections and save
  #####
  
  csb_cen.sf <- csb_cen.sf %>%
    left_join(elevint.dt,
              by = "CSBID") %>%
    left_join(slpint.dt,
              by = "CSBID") %>%
    left_join(aspint.dt,
              by = "CSBID") %>%
    left_join(ssurgo.ext,
              by = "CSBID")
  st_write(csb_cen.sf,
           here("Data", "USDA", "CSB", paste0("CSB_", cntyfips.i, "_2022_centroid_DEM-SSURGO.gpkg")),
           append = FALSE)
  return(print(paste("Intersections for", cntyfips.i, "Complete")))
}
map(fips.v_error, f_csb_spatial)

################################################################################
### Note: Weather intersections
###   - Using PRISM right now; tmin, tmax, and precip
###   - Storing weather data for each CSB as a separate file; process as needed
################################################################################
cnty.dt <- fips_codes %>%
  as.data.table()     %>%
  .[state_code %in% c("17", "18", "26", "39", "55")] %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]

prism.crs <- terra::rast("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\PRISM\\Raw\\2000\\PRISM_ppt_stable_4kmD2_20000101_bil.bil") %>%
  terra::crs(.)
nccs.crs  <- terra::rast("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\NCCS\\sfcWind\\sfcWind_day_ACCESS-CM2_ssp245_r1i1p1f1_gn_2015.nc") %>%
  terra::crs(.)

f_csb_prism <- function(cntyfips.i,
                        yr.i,
                        prism.crs.i = prism.crs){
  
  csb_cen.sf <- st_read(here("Data", "USDA", "CSB", 
                             paste0("CSB_", cntyfips.i, "_2022_centroid_DEM-SSURGO.gpkg"))) %>%
    st_transform(., prism.crs.i)
  csb_cen.sv <- copy(csb_cen.sf) %>%
    vect()
  
  n          <- dim(csb_cen.sv)[1]
  csb_cen_sv.dt <- csb_cen.sv@ptr$getDF() %>%
    as.data.table()                       %>%
    .[,ID     := seq.int(n)]              %>%
    .[,c("ID", "CSBID")]
  
  ### For each year and each data type, set up extract function; run multiproc
  f_extract_prism <- function(data.i,
                              yr.i2           = yr.i,
                              csb_cen.sv.i    = csb_cen.sv,
                              csb_cen_sv.dt.i = csb_cen_sv.dt){
    rast.loc <- paste0("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\PRISM\\Raw\\",
                       yr.i2)
    rast.v   <- list.files(rast.loc,
                           pattern = ".bil$") %>%
      .[str_detect(.,data.i)]                 %>%
      as.data.table()                         %>%
      .[,Date := str_sub(.,
                         start = str_locate(.,
                                            "_20")[,1] +1,
                         end   = str_locate(.,
                                            "_20")[,1] +8)] %>%
      setorder(Date)
    
    rast.r   <- terra::rast(str_c(rast.loc, "\\", rast.v$.))
    rast.ext <- terra::extract(rast.r, csb_cen.sv.i,
                               ID = TRUE)            %>%
      setnames(c("ID", str_c("Date_", rast.v$Date))) %>%
      as.data.table()                                %>%
      merge.data.table(csb_cen_sv.dt.i,
                       by = "ID")                    %>%
      .[,Weather := data.i]                          %>%
      .[,!("ID")]
    return(rast.ext)
  }
  weather.v     <- c("ppt", "tmax", "tmin",
                     "tmean", "vpdmax", "vpdmin")
  
  cl          <- plan(multisession, workers = 6)
  demint.l    <- future_lapply(weather.v, 
                               f_extract_prism,
                               future.scheduling = FALSE,
                               future.seed = TRUE) %>% rbindlist()
  fwrite(demint.l,)
  future:::ClusterRegistry("stop")
}

f_csb_nccs <-  function(cntyfips.i,
                        nccs.crs.i = nccs.crs){
  csb.loc    <- "K:\\CFAES\\AEDE\\Arcdata\\BrianCultice\\INFEWS_LandUse\\Data\\USDA\\CSB"
  csb_cen.sf <- st_read(paste0(csb.loc, "\\", 
                               paste0("CSB_", cntyfips.i, "_2022_centroid_DEM-SSURGO.gpkg"))) %>%
    st_transform(., nccs.crs.i)
  datatype.v <- c("huss", "pr", "sfcWind", "tas", "tasmax", "tasmin")
  ### For each year and each data type, set up extract function; run multiproc
  f_extract_prism <- function(data.i,
                              yr.i,
                              cntyfips.i2     = cntyfips.i,
                              csb_cen.sf.i    = csb_cen.sf){
    csb_cen.sv.i <- copy(csb_cen.sf.i) %>%
      vect()
    n               <-
    csb_cen_sv.dt.i <- as.data.table(csb_cen.sf.i) %>%
      .[,ID     := .I]                             %>%
      .[,c("ID", "CSBID")]
    
    
    rast.loc <- paste0("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\NCCS\\",
                       data.i)
    rast.v   <- list.files(rast.loc,
                           pattern = ".nc$") %>%
      as.data.table()                         %>%
      .[,Year := str_extract(.,
                             pattern = "2[:digit:]{3}")] %>%
      setorder(Year)                                     %>%
      .[Year == yr.i]
    
    ext.r    <- ext(180, 360, 0, 90)
    rast.r   <- terra::rast(str_c(rast.loc, "\\", rast.v$.)) %>%
      terra::crop(., ext.r)
    set.ext(rast.r, c(-180, 0, 0, 90))
    rast.ext <- terra::extract(rast.r, csb_cen.sv.i,
                               ID = TRUE)            %>%
      as.data.table()                                %>%
      merge.data.table(csb_cen_sv.dt.i,
                       by = "ID")                    %>%
      .[,Weather := data.i]                          %>%
      .[,!("ID")]                                    %>%
      .[,Year := yr.i]
    
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
                         new = colB)
    fwrite(rast.ext,
           here("Data", "CSB", paste0("CSB-NCCS_", cntyfips.i2, "_", data.i, "_", yr.i, ".csv")))
    return("done")
  }
  inp1 <- rep(datatype.v, times = length(2015:2050))
  inp2 <- rep(2015:2050, each = length(datatype.v))
  cl          <- plan(multisession, workers = 3)
  demint.l    <- future_map2(inp1, inp2,
                             f_extract_prism,
                             .options = furrr_options(scheduling = FALSE))
  future:::ClusterRegistry("stop")
  
  return(print(cntyfips.i))
}
map(cnty.dt$FIPS_cnty[-1], f_csb_nccs)
