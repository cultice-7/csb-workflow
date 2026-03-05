################################################################################
### Urban Growth Grid Cell Dataset Construction
################################################################################
terraOptions(tempdir="K:\\CFAES\\AEDE\\Arcdata\\BrianCultice\\Temp")
###
# TO DO
# - Primary road distance
# - Cumulative Land Consumption for Each Year
# - Process the slope data (grab cutoffs eg 15 degrees)
# - Keep majority census tract; save separate file of actual shares
###

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
### Load County and MSA files
#####

cbsa.dt <- read_dta(here("Data", "all_pop_centroids.dta")) %>%
  as.data.table()                                          %>%
  .[name != "Honolulu"]
f_gridcreate <- function(msa.i,
                         gridsize.i,
                         cellsize_tol.i){
  
  ##############################################################################
  ### Step 1: Pull counties of interest for the msa
  ##############################################################################
  
  cnty.sf <- tigris::counties(year = 2021)       %>%
    filter(CBSAFP == msa.i)                      %>%
    mutate(FIPS_cnty = str_c(STATEFP, COUNTYFP)) %>%
    st_transform(., crs = "ESRI:102008")
  
  st.v    <- unique(cnty.sf$STATEFP)
  f_trct  <- function(st.i,
                      yr.i){
    trct.sf.o <- tigris::tracts(state = st.i, year = yr.i) 
    stvar.o   <- names(trct.sf.o)[str_detect(names(trct.sf.o), "STATE")][1]
    cntyvar.o <- names(trct.sf.o)[str_detect(names(trct.sf.o), "COUNTY")][1]
    landvar.o <- names(trct.sf.o)[str_detect(names(trct.sf.o), "LAND")][1]
    
    trct.sf.o <- trct.sf.o                                 %>%
      filter(get(landvar.o) != 0)                          %>%
      mutate(FIPS_cnty = str_c(STATEFP, COUNTYFP))         %>%
      filter(FIPS_cnty %in% cnty.sf$FIPS_cnty)
  }
  f_trct90 <- function(st.i,
                       yr.i){
    trct.sf.o <- tigris::tracts(state = st.i, year = yr.i, cb=TRUE)  %>%
      mutate(FIPS_cnty = str_c(ST, CO))                     %>%
      mutate(GEOID     = str_c(ST, CO, TRACTBASE, TRACTSUF))%>%
      filter(FIPS_cnty %in% cnty.sf$FIPS_cnty)
    return(trct.sf.o)
  }
  tryCatch({
    trct20.sf <- map2(st.v, rep(2021, times = length(st.v)), f_trct) %>% do.call(rbind, .) %>%
      st_transform(., crs = "ESRI:102008")
    },
    error = function(cond) {
      trct20.sf <- map2(st.v, rep(2022, times = length(st.v)), f_trct) %>% do.call(rbind, .) %>%
        st_transform(., crs = "ESRI:102008")
    },
    warning = function(cond) {
      trct20.sf <- map2(st.v, rep(2022, times = length(st.v)), f_trct) %>% do.call(rbind, .) %>%
        st_transform(., crs = "ESRI:102008")
    }
  )
  # trct20.sf <- map2(st.v, rep(2021, times = length(st.v)), f_trct) %>% do.call(rbind, .) %>%
  #   st_transform(., crs = "ESRI:102008")
  
  trct10.sf <- map2(st.v, rep(2011, times = length(st.v)), f_trct) %>% do.call(rbind, .) %>%
    st_transform(., crs = "ESRI:102008")
  
  trct00.sf <- map2(st.v, rep(2000, times = length(st.v)), f_trct) %>% do.call(rbind, .) %>%
    st_transform(., crs = "ESRI:102008")                                                 %>%
    rename(GEOID = CTIDFP00)
  
  trct90.sf <- map2(st.v, rep(1990, times = length(st.v)), f_trct90) %>% do.call(rbind, .) %>%
    st_transform(., crs = "ESRI:102008")
  
  
  # Grab bounding box from cnty.dt
  trct_bb.sf <- st_bbox(trct20.sf)
  
  ##############################################################################
  ### Step 2: Create Gridcells; generate gridcell id; identify city center gridcell; distance calcs to city center/other tracts
  ##############################################################################
  
  # Create gridcells w/ counties in MSA.
  grid.sf    <- st_make_grid(trct_bb.sf, cellsize = c(gridsize.i,gridsize.i)) %>%
    st_as_sf()
  
  # Intersect gridcells with tracts; keep only those gridcells that intersect at all
  int.l <- st_intersects(grid.sf, trct20.sf) %>%
    lapply(., length)                        %>%
    unlist()                                 %>%
    as.data.table()                          %>%
    .[,keep.v := ifelse(. > 0,
                        TRUE,
                        FALSE)]
  
  grid.sf <- grid.sf %>%
    .[int.l$keep.v,] %>%
    mutate(CellID = str_c(msa.i, "-", rownames(.)))
  
  # Identify city center (gridcell) via Nic's population weighted approach
  cc.sf <- read_dta(here("Data", "all_pop_centroids.dta")) %>%
    filter(msa_id == msa.i)                                             %>%
    st_as_sf(., coords = c("xcoord", "ycoord"), crs = 4326)             %>%
    st_transform(., crs = st_crs(grid.sf))
  cc_gc.sf <- st_intersection(cc.sf, grid.sf)
  if (dim(cc_gc.sf)[1] == 0){
    cc_gc.sf <- grid.sf[st_nearest_feature(cc.sf, grid.sf),]
  }
  grid.sf <- grid.sf %>%
    mutate(PopWeight_cc = ifelse(CellID == cc_gc.sf$CellID,
                                 1,
                                 0))
  
  # Identify city center via LODES data on workplaces
  files.v <- list.files(here("Data", "LODES", "rolled_up_wac_files"), full.names = TRUE)
  wac.dt  <- lapply(files.v, fread) %>% rbindlist() %>%
    .[,WFH_sh := round(rowSums(.SD),3)/100,
      .SDcols   = patterns("trct")]                 %>%
    .[,WFH_jobs := round(WFH_sh*C000, 3)]           %>%
    .[,FIPS_trt := str_pad(h_geocode, width = 11, side = "left", pad = "0")] %>%
    .[FIPS_trt %in% trct10.sf$GEOID]                                           %>%
    .[C000 == max(C000)]
  
  wac_cc.sf <- copy(trct10.sf)   %>%
    mutate(WAC_cc = ifelse(GEOID == wac.dt$FIPS_trt,
                           1,
                           0)) %>%
    filter(WAC_cc == 1)        %>%
    st_centroid(.)
  wac_cc.sf <- st_intersection(grid.sf, wac_cc.sf)
  grid.sf <- grid.sf %>%
    mutate(WAC_cc = ifelse(CellID == wac_cc.sf$CellID,
                           1,
                           0))
  ##############################################################################
  ### Cell distance to city centers
  ##############################################################################
  grid_cen.sf <- copy(grid.sf) %>%
    st_centroid(.)
  
  cc_dist_pop.m <- as.numeric(st_distance(grid_cen.sf, grid_cen.sf[grid_cen.sf$PopWeight_cc == 1,])/1000)
  cc_dist_wac.m <- as.numeric(st_distance(grid_cen.sf, grid_cen.sf[grid_cen.sf$WAC_cc == 1,])/1000)
  
  grid.sf <- grid.sf %>%
    mutate(DistCC_Pop_km = cc_dist_pop.m,
           DistCC_Wac_km = cc_dist_wac.m)
  
  ##############################################################################
  ### Cell Distance to Primary Roadways
  ##############################################################################
  
  ### Step 1: Load all road networks
  proad90.sf <- st_read(here("Data", "Census", "Roads", paste0("USCensus1990_PrimaryRoads_", msa.i, ".gpkg"))) %>%
    st_union()
  proad00.sf <- st_read(here("Data", "Census", "Roads", paste0("USCensus2000_PrimaryRoads_", msa.i, ".gpkg"))) %>%
    st_union()
  
  proad10.sf <- primary_roads(2011)      %>%
    st_transform(., crs = "ESRI:102008") %>%
    st_intersection(cnty.sf)             %>%
    st_union()
  
  proad20.sf <- primary_roads(2020)      %>%
    st_transform(., crs = "ESRI:102008") %>%
    st_intersection(cnty.sf)             %>%
    st_union()
  
  
  ### Step 2: Calculate distances to each centroid
  r_dist_90.m <- as.numeric(st_distance(grid_cen.sf, proad90.sf))/1000
  r_dist_00.m <- as.numeric(st_distance(grid_cen.sf, proad00.sf))/1000
  r_dist_10.m <- as.numeric(st_distance(grid_cen.sf, proad10.sf))/1000
  r_dist_20.m <- as.numeric(st_distance(grid_cen.sf, proad20.sf))/1000
  
  grid_road.dt <- data.table("CellID" = grid_cen.sf$CellID) %>%
    .[,PrimaryRoad90_dist_km := r_dist_90.m]                %>%
    .[,PrimaryRoad00_dist_km := r_dist_00.m]                %>%
    .[,PrimaryRoad10_dist_km := r_dist_10.m]                %>%
    .[,PrimaryRoad20_dist_km := r_dist_20.m]
  
  ##############################################################################
  ### Tract to Cell Intersections
  ##############################################################################
  
  f_tractint <- function(yr.i){
    yr2.i <- str_sub(yr.i, start = 3)
    grid_trt.sf <- st_intersection(grid.sf, get(paste0("trct", yr2.i, ".sf"))) %>%
      mutate(area_int = as.numeric(st_area(.)))                                %>%
      group_by(CellID)                                                         %>%
      mutate(totarea_int = sum(area_int))                                      %>%
      ungroup()                                                                %>%
      mutate(gridshare_trt = round(area_int/totarea_int,3))
    st_geometry(grid_trt.sf) <- NULL
    grid_trt.dt <- as.data.table(grid_trt.sf)     %>%
      setorder(CellID, gridshare_trt)             %>%
      .[,trtintID := seq(.N), by = "CellID"]      %>%
      setnames(c("GEOID","gridshare_trt"),
               c(str_c("FIPS_trt", yr.i), str_c("gridshare_trt", yr.i)))          %>%
      dcast(CellID + totarea_int ~ trtintID, value.var = c(str_c("FIPS_trt", yr.i), str_c("gridshare_trt", yr.i))) %>%
      setnames("totarea_int",
               str_c("totarea_int", yr.i))
  }
  grid_trt.out.dt <- map(c("2020", "2010", "2000", "1990"), f_tractint)
  grid_trt.m.dt   <- grid_trt.out.dt[[1]] %>%
    merge(grid_trt.out.dt[[2]],
          by = "CellID") %>%
    merge(grid_trt.out.dt[[3]],
          by = "CellID") %>%
    merge(grid_trt.out.dt[[4]],
          by = "CellID")
  fwrite(grid_trt.m.dt,
         here("Data", "Gridded", paste0("GridData", gridsize.i, "m_TractIntersection_CBSA",msa.i, ".csv")))
  
  ##############################################################################
  ### Step 3: Create dataset of assessor data w/ merged grid cells
  ##############################################################################
  
  f_createcons <- function(fips.i,
                           grid.sf.i = grid.sf){
    asmt.dt <- tryCatch(
      {
        asmt.dt <- fread(here("Data", "Raw", "ZAsmt_Raw", paste0("ZAsmt_RRYBFilter_", as.numeric(fips.i), "-220429.csv"))) %>%
          .[str_detect(PropertyLandUseStndCode, "RR|RI")]                   %>%
          .[!is.na(YearBuilt)]                                              %>%
          .[!is.na(PropertyAddressLongitude)]                               %>%
          .[!is.na(PropertyAddressLatitude)]                                %>%
          .[!is.na(LotSizeAcres) | !is.na(LotSizeSquareFeet)]
      },
      error= function(cond) {
        asmt.dt <- data.table("Fail" = vector(mode = "character"))
        print(paste0("County ", fips.i, " Zillow Data Missing."))
        return(asmt.dt)
      }
    )
    
    if (dim(asmt.dt)[1] > 0){
      
      f_precons <- function(i,
                            grid.sf.i2 = grid.sf.i,
                            asmt.dt.i2 = asmt.dt){
        
        ### Note: Filter for year built before the first year and create a geometry
        asmt.dt.i2 <- asmt.dt.i2[YearBuilt < 1990]
        asmt.sf.i2 <- st_as_sf(asmt.dt.i2,
                               coords = c("PropertyAddressLongitude", "PropertyAddressLatitude"),
                               crs    = 4326) %>%
          st_transform(st_crs(grid.sf.i2))
        
        ### Intersections with grid cells
        asmt_grid.sf <- st_intersection(asmt.sf.i2,
                                        grid.sf.i2)
        st_geometry(asmt_grid.sf) <- NULL
        asmt_grid.dt <- as.data.table(asmt_grid.sf)    %>%
          .[,LandCons := sum(LotSizeAcres), by = "CellID"] %>%
          .[,c("LandCons", "CellID")]                      %>%
          unique()
        
        grid.out.dt <- data.table("CellID" = grid.sf.i2$CellID) %>%
          merge(asmt_grid.dt,
                by = "CellID",
                all.x = TRUE)                                   %>%
          .[,Year := "Pre1990"]                                 %>%
          .[is.na(LandCons), LandCons := 0]
        return(grid.out.dt)
      }
      out.l1 <- map(1, f_precons)[[1]]
      
      f_yearcons <- function(yr.i,
                             grid.sf.i2 = grid.sf.i,
                             asmt.dt.i2 = asmt.dt){
        
        ### Note: Filter for year built before the first year and create a geometry
        asmt.dt.i2 <- asmt.dt.i2[YearBuilt == yr.i]
        asmt.sf.i2 <- st_as_sf(asmt.dt.i2,
                               coords = c("PropertyAddressLongitude", "PropertyAddressLatitude"),
                               crs    = 4326) %>%
          st_transform(st_crs(grid.sf.i2))
        
        ### Intersections with grid cells
        asmt_grid.sf <- st_intersection(asmt.sf.i2,
                                        grid.sf.i2)
        st_geometry(asmt_grid.sf) <- NULL
        asmt_grid.dt <- as.data.table(asmt_grid.sf)    %>%
          .[,LandCons := sum(LotSizeAcres), by = "CellID"] %>%
          .[,c("LandCons", "CellID")]                      %>%
          unique()
        
        grid.out.dt <- data.table("CellID" = grid.sf.i2$CellID) %>%
          merge(asmt_grid.dt,
                by = "CellID",
                all.x = TRUE)                                   %>%
          .[,Year := as.character(yr.i)]                        %>%
          .[is.na(LandCons), LandCons := 0]
        return(grid.out.dt)
      }
      out.l2 <- map(1990:2020, f_yearcons) %>% rbindlist() %>%
        rbind(out.l1)                                      %>%
        .[,FIPS_cnty := fips.i]
      return(out.l2)
    }else{
      out.l2 <- data.table("CellID" = NA_character_,
                           "LandCons" = NA_integer_,
                           "LandCons_bi" = NA_integer_,
                           "Year"     = NA_character_,
                           "FIPS_cnty" = fips.i)
      return(out.l2)
    }
  }
  
  cl     <- plan(multisession, workers = 4)
  cnty.v <- cnty.sf$GEOID
  grid_cons.dt <- future_lapply(cnty.v, f_createcons,
                                future.seed = TRUE) %>% rbindlist(fill = TRUE, use.names = TRUE)
  
  ### CHECK FOR MISSING COUNTIES
  countycheck.dt <- copy(grid_cons.dt) %>%
    .[, count              := .N, by = "FIPS_cnty"] %>%
    .[, NonMissingGeocodes := ifelse(count > 1,
                                     1,
                                     0)]                         %>%
    .[,c("FIPS_cnty", "NonMissingGeocodes")]        %>%
    unique()
  
  grid_cons.dt <- grid_cons.dt                                             %>%
    .[!(FIPS_cnty %in% countycheck.dt[NonMissingGeocodes == 0]$FIPS_cnty)] %>%
    .[,!c("FIPS_cnty")]                                                    %>%
    .[,LandCons    := sum(LandCons), by = c("CellID", "Year")]             %>%
    unique()                                                               %>%
    .[,LandCons_bi := ifelse(LandCons > 0,
                             1,
                             0)]
  
  grid_pre90.dt <- copy(grid_cons.dt)[Year == "Pre1990"] %>%
    .[,!c("LandCons_bi")]                                %>%
    setnames("LandCons", "LandCons_Pre1990")
  grid_cons.dt  <- grid_cons.dt[Year != "Pre1990"] %>%
    merge(grid_pre90.dt[,!c("Year")],
          by = "CellID",
          all.x = TRUE)                            %>%
    dcast(CellID + LandCons_Pre1990 ~ Year, value.var = c("LandCons_bi", "LandCons"))
  future:::ClusterRegistry("stop")
  
  ### WITH COUNTYCHECK; REMOVE ANY GRIDCELL IN MISSING COUNTIES
  cntygrab.v       <- names(grid_trt.m.dt) %>%
    .[str_detect(., "FIPS_trt2020")]
  n_cntygrab.v     <- str_c("FIPS_cnty2020_", seq(1:length(cntygrab.v)))
  grid.nocounty.dt <- copy(grid_trt.m.dt) %>%
    .[,(n_cntygrab.v) := (lapply(.SD, function(x) str_sub(x, end = 5))),
      .SDcols = cntygrab.v] %>%
    melt(measure    = patterns("^FIPS_trt2020", "^gridshare_trt2020", "^FIPS_cnty2020"),
         value.name = c("FIPS_trt2020", "gridshare_trt2020", "FIPS_cnty2020")) %>%
    .[,c("CellID", "FIPS_trt2020", "gridshare_trt2020", "FIPS_cnty2020", "totarea_int2020")] %>%
    .[!is.na(gridshare_trt2020)]                                                             %>%
    .[,missing_tag := ifelse((FIPS_cnty2020 %in% countycheck.dt[NonMissingGeocodes == 0]$FIPS_cnty),
                             1,
                             0)] %>%
    .[,totarea_int2020_zillowadj := sum((1-missing_tag)*gridshare_trt2020*totarea_int2020),
      by = "CellID"]                             %>%
    .[,c("CellID", "totarea_int2020_zillowadj")] %>%
    unique()
  
  grid_trt.m.dt <- grid_trt.m.dt %>%
    merge(grid.nocounty.dt,
          by = "CellID")
  
  ### Note: Create cumulative sums by Year
  names.v <- names(grid_cons.dt) %>%
    .[!str_detect(.,"LandCons_bi_")]
  grid_cons_cum.dt <- copy(grid_cons.dt) %>%
    .[,..names.v]                        %>%
    melt(measure = patterns("LandCons_"),
         value.name = "LandAcres",
         variable.name = "Year")         %>%
    .[Year == "LandCons_Pre1990", Year := "1989"] %>%
    .[,Year := as.numeric(str_extract(Year, "[:digit:]{4}"))] %>%
    setorder(Year)                                            %>%
    .[,LandAcres_Cumulative := cumsum(LandAcres), by = "CellID"] %>%
    .[,LandAcres_bi         := ifelse(LandAcres == 0,
                                      0,
                                      1)]                        %>%
    .[,LandAcres_Cumulative_bi := ifelse(LandAcres_Cumulative == 0,
                                         0,
                                         1)]
  
  ##############################################################################
  ### Developable Space (Rule out non-overlaps, space w/ water)
  ##############################################################################
  f_waterload <- function(cnty.i){
    water.sf.o <- tigris::area_water(state = str_sub(cnty.i, end = 2), county = str_sub(cnty.i, start = 3))
    return(water.sf.o)
  }
  water.sf <- map(cnty.sf$FIPS_cnty, f_waterload) %>% bind_rows() %>%
    st_transform(., crs = "ESRI:102008")
  
  grid_water.dt <- st_intersection(grid.sf, water.sf) %>%
    mutate(water_int = as.numeric(st_area(.)))        %>%
    group_by(CellID)                                  %>%
    mutate(totwater_int = sum(water_int)) 
  st_geometry(grid_water.dt) <- NULL
  grid_water.dt <- grid_water.dt %>%
    .[,c("CellID", "totwater_int")] %>%
    unique()
  
  #####
  ### Slope
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
    .[,s_SaveName := str_c("s_", SaveName)]
  
  
  ### Bounding Box for MSA
  cnty_bb.dt <- copy(cnty.sf) %>%
    st_transform(., 4326)     %>%
    st_bbox(.)                %>%
    round()
  cnty_bb.dt[c(1:2)] <- cnty_bb.dt[c(1:2)] - 2
  cnty_bb.dt[c(3:4)] <- cnty_bb.dt[c(3:4)] + 2
  
  inpx <- seq(cnty_bb.dt[1], cnty_bb.dt[3]-1)
  inpy <- seq(cnty_bb.dt[2] + 1, cnty_bb.dt[4])
  inp.g <- expand.grid(inpx, inpy)
  f_slopeint <- function(i,
                         inp.g.i      = inp.g,
                         dem_dl.dt.i  = dem_dl.dt,
                         grid.sf.i    = grid.sf){
    ### former inputs
    minX <- inp.g[[i,1]]
    maxY <- inp.g[[i,2]]
    
    tile.v   <- paste0("n", maxY, "w", str_pad(as.character(abs(minX)), width = 3, side = "left", pad = "0"))
    slope.dt <- copy(dem_dl.dt.i) %>%
      .[AreaID == tile.v]         %>%
      .[,c("s_SaveName")]         %>%
      unique()
    if (dim(slope.dt)[1] != 0){
      slope.r  <- terra::rast(here("Data", "DEM", "Raw", unique(slope.dt$s_SaveName))) %>%
        terra::project(st_crs(grid.sf.i)$proj4string)
      
      grid.sf.i  <- grid.sf.i %>%
        mutate(ID = 1:dim(.)[1])
      grid_ll.sv <- copy(grid.sf.i) %>%
        terra::vect(.)
      st_geometry(grid.sf.i) <- NULL
      grid.dt.i <- copy(grid.sf.i) %>%
        as.data.table()
      ### Step 2: Intersect Buffers w/ NLCD
      ext.dt <- terra::extract(slope.r, grid_ll.sv, 
                               weights = TRUE,
                               ID      = TRUE)   %>%
        as.data.table()                          %>%
        .[!is.nan(slope)]                        %>%
        merge(grid.dt.i,
              by = "ID")                         %>%
        .[,slope   := round(slope, digits = 1)]  %>%
        .[,AreaInt_m := 900*weight]                %>%
        .[,AreaSum_m := sum(AreaInt_m), by = c("CellID", "slope")] %>%
        .[,c("AreaSum_m", "CellID", "slope")]                      %>%
        unique()                                                   %>%
        .[,TileDEM := tile.v]
      return(ext.dt)
    }else{
      ext.dt <- data.table("TileDEM" = tile.v)
      return(ext.dt)
    }
  }
  cl     <- plan(multisession, workers = 4)
  slopeint.dt <- future_lapply(1:dim(inp.g)[1], f_slopeint,
                               future.seed = TRUE) %>% rbindlist(fill = TRUE) %>%
    .[,AreaSlope_m2 := sum(AreaSum_m), by = c("CellID", "slope")]                    %>%
    .[,c("CellID", "slope", "AreaSlope_m2")]                                         %>%
    unique()
  future:::ClusterRegistry("stop")
  
  fwrite(slopeint.dt,
         here("Data", "Gridded", paste0("GridData", gridsize.i, "m_SlopeDist_CBSA",msa.i, ".csv")))
  
  ### Note: Create slope variables (% over 15, % over 30, mean slope, max slope)
  slopeint.dt <- slopeint.dt                                                 %>%
    .[,share := AreaSlope_m2/sum(AreaSlope_m2, na.rm = TRUE), by = "CellID"] %>%
    .[,MeanSlope_deg := sum(slope*share, na.rm = TRUE), by = "CellID"]       %>%
    .[,MaxSlope_deg  := max(slope, na.rm = TRUE), by = "CellID"]             %>%
    .[,Slope_gt15deg  := round(sum(ifelse(slope > 15,
                                          share,
                                          0))/sum(share),4), by = "CellID"]  %>%
    .[,Slope_gt30deg  := round(sum(ifelse(slope > 30,
                                          share,
                                          0))/sum(share),4), by = "CellID"]  %>%
    .[,c("CellID", "MeanSlope_deg", "MaxSlope_deg", "Slope_gt15deg", "Slope_gt30deg")] %>%
    unique()
  
  #####
  ### Soil Type
  #####
  
  ### Step 1: Assemble map units (and underlying raw data)
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
      .[!is.na(sandtotal_r) & !is.na(silttotal_r) & !is.na(claytotal_r)] %>%
      .[sandtotal_r + silttotal_r + claytotal_r == 100]                  %>%
      .[!is.na(comppct_r)]                                               %>%
      .[,comppct_r := comppct_r/100]                                     %>%
      .[,comppct_r := comppct_r/sum(comppct_r), by = "mukey"]            %>%
      .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact", "nccpi3all"))) := (lapply(.SD, function(x) comppct_r*x)),
        .SDcols = c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact", "nccpi3all")] %>%
      .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact", "nccpi3all"))) := (lapply(.SD, function(x) sum(x, na.rm = TRUE))), 
        by = c("mukey"),
        .SDcols = str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact", "nccpi3all"))] %>%
      .[,c(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "kwfact", "nccpi3all")), "mukey")] %>% 
      unique() %>%
      setnames("mukey", "MUKEY") %>%
      .[,MUKEY := as.numeric(MUKEY)]
    return(ssurgo_walk.dt)
  }
  stabbv.v <- as.data.table(fips_codes) %>%
    .[state_code %in% unique(st.v)]     %>%
    .[,c("state")]                      %>%
    unique()
  ssurgo_mu.dt <- map(stabbv.v$state, f_mapunitmerge) %>% rbindlist() %>%
    unique()
  
  ### Step 2: Load and merge rasters for MSA
  f_msamerge <- function(msa.i2,
                         stabbv.v2 = stabbv.v$state,
                         crs.v     = st_crs(grid.sf)$proj4string){
    files.v  <- list.files(here("Data", "SSURGO"),
                           pattern = paste0(stabbv.v2, collapse = "|")) %>%
      .[str_detect(., ".tif$")]
    ssurgo.l <- sprc(map(here("Data", "SSURGO", files.v), terra::rast))
    ssurgo.out <- terra::merge(ssurgo.l) %>%
      terra::project(crs.v)
    return(ssurgo.out)
  }
  ssurgo.r <- map(msa.i, f_msamerge)[[1]]
  
  ### Step 3: Extract map units from gSSURGO by grid cells
  grid.sv <- copy(grid.sf) %>%
    vect(.)
  
  n          <- dim(grid.sv)[1]
  grid.sv.dt <- grid.sv@ptr$getDF() %>%
    as.data.table()                %>%
    .[,ID     := seq.int(n)]       %>%
    .[,c("ID", "CellID")]
  
  ssurgo.ext <- terra::extract(ssurgo.r, grid.sv,
                               weights = TRUE,
                               ID      = TRUE)                             %>%
    as.data.table(key = c("ID"))                                           %>%
    .[,AreaInt_m := 900*weight]                                            %>%
    .[,AreaSumMU_m := sum(AreaInt_m, na.rm = TRUE), by = c("ID", "MUKEY")] %>%
    .[,c("AreaSumMU_m", "ID", "MUKEY")]                      %>%
    unique()                                               %>%
    merge(grid.sv.dt,
          by = "ID",
          all.x = TRUE)                                    %>%
    .[,!c("ID")]                                           %>%
    .[,AreaTot_m   := sum(AreaSumMU_m, na.rm = TRUE), by = c("CellID")]  %>%
    .[,mu_share    := AreaSumMU_m/AreaTot_m]                             %>%
    .[,MUKEY := as.double(as.character(MUKEY))]                          %>%
    merge(ssurgo_mu.dt,
          by    = "MUKEY",
          all.x = TRUE)                                        %>%
    .[!is.na(MUKEY)]                                           %>%
    .[!is.na(w_sandtotal_r)]                                   %>%
    .[,mu_share    := mu_share/sum(mu_share, na.rm = TRUE), by = "CellID"] %>%
    .[,(c("gw_sandtotal_r", "gw_silttotal_r", "gw_claytotal_r", "gw_kwfact", "gw_nccpi3all")) :=
        (lapply(.SD, function(x) sum(x*mu_share, na.rm = TRUE))),
      by = c("CellID"),
      .SDcols = c("w_sandtotal_r", "w_silttotal_r", "w_claytotal_r", "w_kwfact", "w_nccpi3all")] %>%
    .[,c("CellID","gw_sandtotal_r", "gw_silttotal_r", "gw_claytotal_r", "gw_kwfact", "gw_nccpi3all")] %>%
    unique()
  
  ##############################################################################
  ### Create Final Single Dataset (single row item per cell)
  ##############################################################################
  
  ### Note: Just take max tract for current purposes
  grid_trt.m.dt <- grid_trt.m.dt %>%
    .[,c("CellID", "totarea_int2020", "totarea_int2020_zillowadj",
         "FIPS_trt2020_1", "FIPS_trt2010_1", "FIPS_trt2000_1", "FIPS_trt1990_1")]
  
  grid.dt       <- copy(grid.sf)
  st_geometry(grid.dt) <- NULL
  grid.dt       <- grid.dt %>%
    as.data.table()
  
  
  grid_f.dt <- copy(grid_cons_cum.dt)        %>%
    merge(grid.dt,
          by = "CellID",
          all.x = TRUE)                      %>%
    merge(grid_road.dt,
          by = "CellID",
          all.x = TRUE)                      %>%
    merge(grid_water.dt[,c("CellID", "totwater_int")],
          by = "CellID",
          all.x = TRUE)                      %>%
    merge(grid_trt.m.dt,
          by = "CellID",
          all.x = TRUE)                      %>%
    merge(slopeint.dt,
          by = "CellID",
          all.x = TRUE)                      %>%
    merge(ssurgo.ext,
          by = "CellID",
          all.x = TRUE)                      %>%
    .[totarea_int2020_zillowadj >= cellsize_tol.i*gridsize.i^2]
  fwrite(grid_f.dt,
         here("Data", "Gridded", paste0("GridData", gridsize.i, "m_CBSA", msa.i, ".csv")))
  
  grid.sf <- grid.sf[,c("x", "CellID")] %>%
    left_join(grid_trt.m.dt,
              by = "CellID")            %>%
    filter(!is.na(totarea_int2020_zillowadj))
  st_write(grid.sf,
           here("Data", "Gridded", paste0("Grid", gridsize.i, "m_CBSA", msa.i, ".gpkg")),
           append = FALSE)
  
  ### Clear Temp File
  tempfiles.v <- list.files(path = "K:\\CFAES\\AEDE\\Arcdata\\BrianCultice\\Temp",
                            full.names = TRUE)
  file.remove(tempfiles.v)
}
terraOptions(tempdir="K:\\CFAES\\AEDE\\Arcdata\\BrianCultice\\Temp")

### From Fresno (e.g)
index.v <- which(cbsa.dt$msa_id == "12060")
inp.v   <- as.character(cbsa.dt$msa_id)[index.v:length(as.character(cbsa.dt$msa_id))]
pmap(list(inp.v,
          rep(1000, length(inp.v)),
          rep(1, length(inp.v))),
     f_gridcreate)
