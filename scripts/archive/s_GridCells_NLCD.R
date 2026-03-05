################################################################################
### Land Use Model
### Creating Grid-cell Dataset
################################################################################

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

### External data file locations
prism.loc <- "K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM"
fips.dt      <- fips_codes                           %>%
  as.data.table()                                    %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]    %>%
  .[state_code %in% c("17", "18", "26", "39", "55")]

#####
### Step 1: Creation of grid for GL region
#####
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
  
  # Create grid as polygon for specific intersections
  grid_poly.sf <- st_make_grid(cnty_bb.sf,
                               cellsize = 300,
                               what     = "polygons") %>%
    st_as_sf()                                        %>%
    st_set_crs(.,"ESRI:102008")                       %>%
    mutate(GridID = str_c("G_", 1:dim(.)[1]))

  # Create mesh of centroids; create coords variables
  grid.sf <- st_make_grid(cnty_bb.sf,
                          cellsize = 300,
                          what     = "centers") %>%
    st_as_sf()                                  %>%
    st_set_crs(.,"ESRI:102008") 
  
  grid.x <- unlist(purrr::transpose(grid.sf$x)[[1]])
  grid.y <- unlist(purrr::transpose(grid.sf$x)[[2]])
  grid.sf <- grid.sf          %>%
    mutate(lon_proj = grid.x,
           lat_proj = grid.y) %>%
    mutate(GridID = str_c("G_", 1:dim(.)[1]))
  
  # Create focal mesh of gridcells
  focal.sf <- st_make_grid(cnty_bb.sf,
                           cellsize = 3000,
                           what     = "polygons") %>%
    st_as_sf()                                    %>%
    mutate(FocalID = str_c("F_", 1:dim(.)[1]))
  
  focal_cen.sf <- st_centroid(focal.sf)
  f_grid_proj.x <- unlist(purrr::transpose(focal_cen.sf$x)[[1]])
  f_grid_proj.y <- unlist(purrr::transpose(focal_cen.sf$x)[[2]])
  
  focal_cen.sf <- focal_cen.sf %>%
    st_transform(., 4326)
  f_grid.x <- unlist(purrr::transpose(focal_cen.sf$x)[[1]])
  f_grid.y <- unlist(purrr::transpose(focal_cen.sf$x)[[2]])
  
  focal.sf <- focal.sf %>%
    mutate(centroid_lon_proj = f_grid_proj.x,
           centroid_lat_proj = f_grid_proj.y,
           centroid_lon      = f_grid.x,
           centroid_lat      = f_grid.y)
  
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
  
  ### Filter grid_poly for centroids in counties
  
  grid.dt <- copy(grid.sf)
  st_geometry(grid.dt) <- NULL
  
  grid_poly.sf <- grid_poly.sf %>%
    filter(GridID %in% grid.sf$GridID) %>%
    left_join(grid.dt,
              by = "GridID")
  st_write(grid_poly.sf,
           here("Data", "Gridded", "GL_Gridded_300m-3km_polygons.gpkg"),
           append = FALSE)
  
  ### Check for focal zones w/ ANY 300m gridcells contained
  focal.sf <- focal.sf %>%
    filter(FocalID %in% grid.sf$FocalID)
  st_write(focal.sf,
           here("Data", "Gridded", "GL_Gridded_300m-3km_Focal.gpkg"))
  return(print("Gridded Saved"))
}
map(1, f_gridcentroids)

#####
### Step 2: Intersection w/ NLCD
#####
f_grid.makecounty <- function(i){
  
  ### Step 1: Load full polygon file
  grid_poly.sf <- st_read(here("Data", "Gridded", "GL_Gridded_300m-3km_polygons.gpkg"))
  fips.dt      <- fips_codes                           %>%
    as.data.table()                                    %>%
    .[,FIPS_cnty := str_c(state_code, county_code)]    %>%
    .[state_code %in% c("17", "18", "26", "39", "55")]
  
  f_grid.filcounty <- function(FIPS_cnty.i,
                               grid_poly.sf.i = grid_poly.sf){
    ### Filter the county of interest from the main file
    grid_poly.sf.i <- grid_poly.sf.i   %>%
      filter(FIPS_cnty == FIPS_cnty.i)
    st_write(grid_poly.sf.i,
             here("Data", "Gridded", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, "_polygons.gpkg")),
             append = FALSE)
    return(paste0("Filter ", FIPS_cnty.i, " completed."))
  }
  map(fips.dt$FIPS_cnty, f_grid.filcounty)
}
map(1, f_grid.makecounty)

f_grid.NLCD <- function(i){
  
  ### Step 1: Save CRS for NLCD data
  nlcd.crs <- terra::rast(here("Data", "NLCD", "Raw", "NLCD_Ohio_2008.tif")) %>%
    terra::crs(.)
  
  fips.dt      <- fips_codes                           %>%
    as.data.table()                                    %>%
    .[,FIPS_cnty := str_c(state_code, county_code)]    %>%
    .[state_code %in% c("17", "18", "26", "39", "55")]
  
  
  ### Step 2: Filter out a county of interest
  f_grid.NLCD.cnty <- function(FIPS_cnty.i,
                               nlcd.crs.i     = nlcd.crs){
    cnty.sv <- tigris::counties(cb = TRUE) %>%
      filter(GEOID == FIPS_cnty.i)         %>%
      st_transform(.,
                   crs = nlcd.crs.i)       %>%
      terra::vect(.)
    
    
    grid_poly.sf <- st_read(here("Data", "Gridded", 
                                 paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, "_polygons.gpkg"))) %>%
      st_transform(.,
                   crs = nlcd.crs.i)
    
    ### Note: Get the fips code. This crosswalks different naming structures for states
    #       because of some poor naming choices
    fips    <- c("MI", "OH", "IL", "IN", "WI",
                 "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
                 "26", "39", "17", "18", "55")
    fipsmat <- matrix(data = fips,
                      nrow = 5,
                      ncol = 3,
                      byrow = FALSE)
    fips.abbv <- fipsmat[fipsmat[,3]==str_sub(FIPS_cnty.i, end = 2),1]
    fips.n    <- fipsmat[fipsmat[,3]==str_sub(FIPS_cnty.i, end = 2),2]
    
    f_grid.NLCD.int <- function(yr.i,
                                FIPS_cnty.i2      = FIPS_cnty.i,
                                stname.i          = fips.n,
                                grid_poly.sf.cnty = grid_poly.sf,
                                nlcd.crs.i2       = nlcd.crs.i){
      
      rast.r   <- terra::rast(here("Data", "NLCD", "Raw", paste0("NLCD_", stname.i, "_", yr.i, ".tif")))
      rast_crop.r <- terra::crop(rast.r,
                                 cnty.sv)
      
      rast.dt     <- data.frame("lon" = xFromCell(rast_crop.r, 1:ncell(rast_crop.r)),
                                "lat" = yFromCell(rast_crop.r, 1:ncell(rast_crop.r)),
                                "Value" = extract(rast_crop.r, xyFromCell(rast_crop.r, 1:ncell(rast_crop.r)))) %>%
        st_as_sf(.,
                 coords = c("lon", "lat"),
                 crs    = nlcd.crs.i2)
        
      grid.int.sf <- st_intersection(grid_poly.sf.cnty, rast.dt)      %>%
        as.data.table()                                               %>%
        .[,Year := yr.i]                                              %>%
        .[,c("GridID", "FIPS_cnty", "Year", "NLCD.Land.Cover.Class")] %>%
        .[,Count := .N, by = c("GridID", "NLCD.Land.Cover.Class")]    %>%
        unique()
      fwrite(grid.int.sf,
             here("Data", "Gridded", "NLCD", paste0("NLCD-GL_Gridded_300m-3km_", FIPS_cnty.i2, "-", yr.i, ".csv")))
      return(paste0("Intersection for ", yr.i, " and ", FIPS_cnty.i2))
    }
    yr.v <- c("2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019", "2021")
    map(yr.v,
        f_grid.NLCD.int)
  }
  cl          <- plan(multisession, workers = 4)
  out.dist.dt <- future_mapply(f_grid.NLCD.cnty, 
                               fips.dt$FIPS_cnty,
                               future.seed = TRUE,
                               future.scheduling = FALSE)
  future:::ClusterRegistry("stop")
}
map(1, f_grid.NLCD)

#####
### Soils and slope/elev
#####
f_grid.makecounty_cen <- function(i){
  ### Step 1: Load full polygon file
  grid_cen.sf <- st_read(here("Data", "Gridded", "GL_Gridded_300m-3km.gpkg"))
  fips.dt      <- fips_codes                           %>%
    as.data.table()                                    %>%
    .[,FIPS_cnty := str_c(state_code, county_code)]    %>%
    .[state_code %in% c("17", "18", "26", "39", "55")]
  
  f_grid.filcounty <- function(FIPS_cnty.i,
                               grid_cen.sf.i = grid_cen.sf){
    ### Filter the county of interest from the main file
    grid_cen.sf.i <- grid_cen.sf.i   %>%
      filter(FIPS_cnty == FIPS_cnty.i)
    st_write(grid_cen.sf.i,
             here("Data", "Gridded", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, ".gpkg")),
             append = FALSE)
    return(paste0("Filter ", FIPS_cnty.i, " completed."))
  }
  map(fips.dt$FIPS_cnty, f_grid.filcounty)
}
map(1, f_grid.makecounty_cen)

f_grid.ssurgo_dem <- function(FIPS_cnty.i){
  
  #####
  ### Step 0: Set up crosswalks and other items
  #####
  stname.dt <- fips_codes        %>%
    as.data.table()              %>%
    .[,c("state", "state_code")] %>%
    unique()
  stname.v  <- stname.dt[state_code == str_sub(FIPS_cnty.i,
                                               end = 2)]$state
  
  ### Step 1: Load grid centroids
  grid_cen.sf <- st_read(here("Data", "Gridded", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, ".gpkg")))
  
  #####
  ### Step 2: Load DEM files (slope, aspect, elev)
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
    .[FIPS_cnty == FIPS_cnty.i]                     %>%
    .[,s_SaveName := str_c("s_", SaveName)]         %>%
    .[,a_SaveName := str_c("a_", SaveName)]
  
  f_demint   <- function(i,
                         dem_dl.dt.i  = dem_dl.dt,
                         csb_cen.sf.i = grid_cen.sf){
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
      .[,c("GridID", "Elev")]                     %>%
      unique()                                   %>%
      .[,TileDEM := dem_dl.dt.i[i]$AreaID]
    
    slope.dt <- terra::extract(slope.r, csb_cen.sv,
                               ID      = TRUE)   %>%
      as.data.table()                          %>%
      .[!is.nan(slope)]                        %>%
      merge(csb_cen.dt.i,
            by = "ID")                         %>%
      .[,slope   := round(slope, digits = 1)]  %>%
      .[,c("GridID", "slope")]                  %>%
      unique()                                 %>%
      .[,TileDEM := dem_dl.dt.i[i]$AreaID]
    
    asp.dt   <- terra::extract(asp.r, csb_cen.sv,
                               ID      = TRUE) %>%
      as.data.table()                          %>%
      .[!is.nan(aspect)]                       %>%
      merge(csb_cen.dt.i,
            by = "ID")                         %>%
      .[,aspect   := round(aspect, digits = 1)]%>%
      .[,c("GridID", "aspect")]                 %>%
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
    .[, Elev := round(mean(Elev), 1), by = "GridID"] %>%
    .[,!c("TileDEM")]                                %>%
    unique()
  
  slpint.dt   <- map2(1:length(demint.l),
                      rep(2, times = length(demint.l)), 
                      f_dempull) %>% rbindlist() %>%
    .[, slope := round(mean(slope), 1), by = "GridID"] %>%
    .[,!c("TileDEM")]                                 %>%
    unique()
  
  aspint.dt   <- map2(1:length(demint.l),
                      rep(3, times = length(demint.l)),
                      f_dempull) %>% rbindlist() %>%
    .[, aspect := round(mean(aspect), 1), by = "GridID"] %>%
    .[,!c("TileDEM")]                                    %>%
    unique()
  
  #####
  ### Step 3: gSSURGO Intersection
  #####
  ssurgo.v <- list.files(here("Data", "SSURGO"),
                         pattern = stname.v) %>%
    .[str_detect(., ".tif$")]
  ssurgo.r <- terra::rast(here("Data", "SSURGO", ssurgo.v))
  grid_cen_ssurgo.sf <- copy(grid_cen.sf) %>%
    st_transform(terra::crs(ssurgo.r))
  
  grid_cen_ssurgo.sv <- copy(grid_cen_ssurgo.sf) %>%
    vect(.)
  
  n          <- dim(grid_cen_ssurgo.sv)[1]
  grid_cen_ssurgo_sv.dt <- grid_cen_ssurgo.sv@ptr$getDF() %>%
    as.data.table()                                     %>%
    .[,ID     := seq.int(n)]                            %>%
    .[,c("ID", "GridID")]
  
  ssurgo.ext <- terra::extract(ssurgo.r, grid_cen_ssurgo.sv,
                               ID      = TRUE)          %>%
    as.data.table(key = c("ID"))                        %>%
    merge(grid_cen_ssurgo_sv.dt,
          by = "ID")                                    %>%
    .[,!c("ID")]
  
  #####
  ### Step 4: Merge intersections and save
  #####
  
  grid_cen.sf <- grid_cen.sf %>%
    left_join(elevint.dt,
              by = "GridID") %>%
    left_join(slpint.dt,
              by = "GridID") %>%
    left_join(aspint.dt,
              by = "GridID") %>%
    left_join(ssurgo.ext,
              by = "GridID") %>%
    as.data.table()         %>%
    .[,c("GridID", "slope", "aspect", "Elev", "MUKEY")]
  
  fwrite(grid_cen.sf,
         here("Data", "Gridded", "SSURGO-DEM", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, "SSURGO-DEM_.csv")))
  return(print(paste("Intersections for", FIPS_cnty.i, "Complete")))
}
map(fips.dt$FIPS_cnty, f_grid.ssurgo_dem)

#####
### Create SSURGO data for merge
#####

f_mapunitmerge <- function(st.i){
  ### Note: Load layers of interest from 
  # ssurgo.layers <- st_layers(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")))$name
  
  layer1.dt     <- st_read(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")),
                           layer = "Valu1")  %>%
    as.data.table()                          %>%
    .[,c("nccpi3all", "pwsl1pomu", "pctearthmc", "rootznemc", "rootznaws", "droughty", "mukey")]
  
  component.dt  <- st_read(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")),
                           layer = "component") %>%
    as.data.table()                             %>%
    .[,c("cokey", "comppct_r", "mukey")]
  
  horizon.dt <- st_read(here("Data", "SSURGO", paste0("gSSURGO_", st.i, ".gdb")),
                        layer = "chorizon") %>%
    as.data.table()                         %>%
    .[,c("sandtotal_r", "silttotal_r", "claytotal_r", "hzdept_r", "cokey", "chkey")] %>%
    setorder(cokey, hzdept_r)                                                        %>%
    .[,HorizonCount := seq(.N), by = "cokey"]                                        %>%
    .[HorizonCount <= 3]                                                             %>%
    .[,(c("sandtotal_r", "silttotal_r", "claytotal_r")) := lapply(.SD, mean, na.rm = TRUE),
      by = c("cokey"),
      .SDcols = c("sandtotal_r", "silttotal_r", "claytotal_r")]          %>%
    .[,!c("chkey", "HorizonCount")]                                      %>%
    unique()
  
  ssurgo_walk.dt <- copy(horizon.dt) %>%
    merge(copy(component.dt),
          by    = "cokey",
          all.x = TRUE)              %>%
    merge(copy(layer1.dt),
          by    = "mukey",
          all.x = TRUE)              %>%
    .[!is.na(sandtotal_r) & !is.na(silttotal_r) & !is.na(claytotal_r)] %>%
    .[sandtotal_r + silttotal_r + claytotal_r == 100]                  %>%
    .[!is.na(comppct_r)]                                               %>%
    .[,comppct_r := comppct_r/100]                                     %>%
    .[,comppct_r := comppct_r/sum(comppct_r), by = "mukey"]            %>%
    .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r"))) := (lapply(.SD, function(x) comppct_r*x)),
      .SDcols = c("sandtotal_r", "silttotal_r", "claytotal_r")] %>%
    .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r"))) := (lapply(.SD, function(x) sum(x, na.rm = TRUE))), 
      by = c("mukey"),
      .SDcols = str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r"))] %>%
    .[,c(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r")), 
         "nccpi3all", "pwsl1pomu", "pctearthmc", "rootznemc", "rootznaws", "droughty", "mukey")] %>% 
    unique() %>%
    setnames("mukey", "MUKEY") %>%
    .[,MUKEY := as.numeric(MUKEY)] %>%
    .[,MajorityClay := ifelse(w_claytotal_r > w_sandtotal_r & w_claytotal_r > w_silttotal_r,
                              1,
                              0)]  %>%
    .[,MajoritySilt := ifelse(w_silttotal_r > w_sandtotal_r & w_silttotal_r > w_claytotal_r,
                              1,
                              0)]  %>%
    .[,MajoritySand := ifelse(w_sandtotal_r > w_silttotal_r & w_sandtotal_r > w_claytotal_r,
                              1,
                              0)]
  return(ssurgo_walk.dt)
}
ssurgo_mu.dt <- map(unique(fips.dt$state), f_mapunitmerge) %>% rbindlist() %>%
  unique()

#####
### Other intersections of potential use
#####

### Protected Areas
f_PAD <- function(i,
                  cnty.v = fips.dt$FIPS_cnty){
  
  ### Protected lands (e.g. federal lands, state lands, etc)
  pad.sf <- st_read(here("Data", "USGS", "PAD", "PADUS3_0_Region_3_GeoPackage", "PADUS3_0Region3.gpkg"),
                    layer = "PADUS3_0Combined_Region3") %>%
    filter(Category %in% c("Fee", "Designation"))       %>%
    filter(Own_Type %in% c("FED", "STAT", "LOC", "DIST")) %>%
    .[,c("Own_Type", "Own_Name", "d_GAP_Sts", "FeatClass", "Des_Tp", "Mang_Type", "Shape")]
  
  
  
  ### Load each county centroid; intersect to confirm whether protected (Yes/no)
  f_PADCounty <- function(FIPS_cnty.i,
                          pad.sf.i = pad.sf){
    grid_cen.sf <- st_read(here("Data", "Gridded", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, ".gpkg"))) %>%
      .[,c("FocalID", "GridID", "FIPS_cnty", "geom")] %>%
      st_transform(., st_crs(pad.sf.i))
    
    grid_pad_int.sf <- st_intersection(grid_cen.sf, pad.sf.i) %>%
      as.data.table()                                         %>%
      .[,is_pad := 1]                                         %>%
      .[,is_pad_Fed := ifelse(Own_Type == "FED",
                              1,
                              0)]                             %>%
      .[,is_pad_St  := if_else(Own_Type == "STAT",
                               1,
                               0)]                            %>%
      .[,is_pad_Loc := if_else(Own_Type == "LOC",
                               1,
                               0)]                            %>%
      .[,c("GridID", "is_pad", "is_pad_Fed", "is_pad_St", "is_pad_Loc")]
    
    grid_cen.sf <- as.data.table(grid_cen.sf) %>%
      merge(grid_pad_int.sf,
            by = "GridID",
            all.x = TRUE)
    fwrite(grid_cen.sf,
           here("Data", "Gridded", "PAD", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, "PAD.csv")))
  }
  map(cnty.v, f_PADCounty)
}
map(1, f_PAD)

### Roadways (2010 roadways due to strange data in 2000 census)
f_county_roads <- function(FIPS_cnty.i){
  if (!file.exists(here("Data", "Gridded", "Roads", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, "-Roads2010.gpkg")))){
    ### Step 1: Load Roads
    roads.sf <- tigris::primary_secondary_roads(state = str_sub(FIPS_cnty.i, end = 2),
                                                year  = 2011)
    proads.sf <- copy(roads.sf) %>%
      filter(MTFCC == "S1100")
    sroads.sf <- copy(roads.sf) %>%
      filter(MTFCC == "S1200")
    
    ### Step 2: Load centroids
    grid.sf  <- st_read(here("Data", "Gridded", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, ".gpkg"))) %>%
      st_transform(., crs = st_crs(proads.sf))
    
    ### Step 3: calculate distances via functions
    f_finddist <- function(i,
                           roads.sf.i = roads.sf){
      nearest.index_p <- st_nearest_feature(grid.sf[i,], roads.sf.i[roads.sf.i$MTFCC == "S1100",])
      nearest.index_s <- st_nearest_feature(grid.sf[i,], roads.sf.i[roads.sf.i$MTFCC == "S1200",])
      dist_p          <- st_distance(grid.sf[i,], roads.sf.i[nearest.index_p,])
      dist_s          <- st_distance(grid.sf[i,], roads.sf.i[nearest.index_s,])
      return(list(dist_p, dist_s))
    }
    dist.l <- map(1:dim(grid.sf)[1], f_finddist) %>%
      purrr::transpose(.)
    
    ### Step 4: Create variables for point data
    out.dt <- data.table("GridID" = grid.sf$GridID,
                         "PRoad_dist_km" = as.numeric(dist.l[[1]])/1000,
                         "SRoad_dist_km" = as.numeric(dist.l[[2]])/1000)
    fwrite(out.dt,
           file = here("Data", "Gridded", "Roads", paste0("GL_Gridded_300m-3km_", FIPS_cnty.i, "-Roads2010.csv")))
  }
}
cl          <- plan(multisession, workers = 4)
future_lapply(fips.dt$FIPS_cnty, 
              f_county_roads,
              future.seed = TRUE,
              future.scheduling = FALSE)
future:::ClusterRegistry("stop")

################################################################################
### Focal Areas Intersections
################################################################################

### i to j focal distance matrix (maybe network distance later??)


### Weather
f_focal_PRISM <- function(i){
  ### Step 1: Load Focal Areas
  fgrid.sf <- st_read(here("Data", "Gridded", "GL_Gridded_300m-3km_Focal.gpkg"))
  
  ### Step 2: Program and run function to load PRISM and intersect (note that PRISM is daily)
  f_focal_PRISM_yr <- function(yr.i){
    prism.loc <- "K:/CFAES/AEDE/Arcdata/Data_PRISM/Data/PRISM/Raw"
    ### Load Rasters from File
    rast.r <- terra::rast(paste0(prism.loc, "/", yr.i))
  }
}


