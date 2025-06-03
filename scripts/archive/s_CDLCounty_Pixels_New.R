######################################################################
##### CDL - Pixel by Pixel transitions
##### Description: Identify county level land use ineractions
##### Date: 4/20/22
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
               "rgdal",  # Projection schemes
               "furrr", 
               "pryr",
               "rmapshaper",
               "tigris")
pacman::p_load(char = libraries)

### Step 1: Cut the raw CDL rasters into states. This isn't essential for this to run, 
#         but it is essential if you want/need to run these scripts in parallel.
#         For my purposes, I'm doing this for these 5 states. You can make this a vector
#         of all states as needed.
states.v <- c("Ohio", "Indiana", "Illinois", "Michigan", "Wisconsin")
inp1     <- 2008:2018
CDL_CropSave <- function(year.i,
                         state.v.i = states.v){
  
  ### Note: Identify the location of the raw CDL raster data in the folder you moved them:
  #       for me, that is in the ./Data/Raw/CDL folder
  rasters  <- list.files(here("Data", "CDL", "Raw")) %>%
    as.data.table()                                  %>%
    .[str_detect(., ".img$")]                        %>%
    .[str_detect(., as.character(year.i))]
  
  ### Note: Use that location identified above for year year_v and load in the 
  #       full raster data
  cdl.i     <- terra::rast(here("Data", "CDL", "Raw", rasters[[1,1]]))
  cdl.crs.i <- terra::crs(cdl.i)
  
  ### Note: Load in a US state shapefile; I use this to filter a specific state,
  #       then intersect that specific state w/ the raster to create a state specific
  #       version of the CDL data. Be sure the CRS is the same for both geospatial data
  state.loc   <- here("Data", "Census", "US_state_2018.shp")
  state.sf.i  <- st_read(state.loc) %>%
    st_transform(., cdl.crs.i) # read it in
  
  f_statecrop <- function(state.i2,
                          year.i2     = year.i,
                          state.sf.i2 = state.sf.i,
                          cdl.i2      = cdl.i){
    state.sf.i2 <- state.sf.i2[state.sf.i2$NAME == state.i2,]
    state.sv.i2 <- terra::vect(state.sf.i2)
    
    ### Note: Cut out land use from stack (crop first to make this go more quickly)
    cdl_crop.i2  <- terra::crop(cdl.i2, state.sv.i2)
    cdl_crop.i2[is.nan(cdl_crop.i2)] <- 0
    
    ### Note: Save the cropped file
    tempname <- here("Data", "CDL", "Raw", paste0("CDL_", state.i2, "_", year.i2, ".tif"))
    terra::writeRaster(cdl_crop.i2, tempname, 
                       overwrite = TRUE)
    return(paste0("Crop for ", state.i2, year.i2, " complete."))
  }
  map(state.v.i, f_statecrop)
}
map(inp1, CDL_CropSave)


### Step 2: Intersect Counties w/ CDL. This grabs all the pixels that fall within a given
#         county. For my purposes, I'm doing this for these 5 states. You can make this a vector
#         of all states as needed.
st.v <- rep(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"),
            times = length(inp1))
yr.v <- rep(inp1,
            each = 5)
f_CDLCountyCellwise <- function(state.i, year.i){
  ### Note: Get the fips code. This crosswalks different naming structures for states
  #       because of some poor naming choices
  fips    <- c("MI", "OH", "IL", "IN", "WI",
               "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
               "26", "39", "17", "18", "55")
  fipsmat <- matrix(data = fips,
                    nrow = 5,
                    ncol = 3,
                    byrow = FALSE)
  fips   <- fipsmat[fipsmat[,2]==state.i,1]
  fips.n <- fipsmat[fipsmat[,2]==state.i,3]
  ### Note: Load data and change projection schemes
  county.sf.i  <- st_read(here("Data", "Census", paste0("US_county_2018.shp"))) %>%
    filter(STATEFP == fips.n)                                                   %>%
    .[,c("GEOID", "geometry")]
  
  f_cnty <- function(cnty.i2, 
                     county.sf.i2 = county.sf.i, 
                     state.i2     = state.i, 
                     year.i2     = year.i){
    
    ### Note: Load in raster file
    cdl.i2  <- terra::rast(here("Data", "CDL", "Raw", paste0("CDL_", state.i2, "_", year.i2, ".tif")))
    
    ### Note: Confirm that coordinate reference systems for the two files im intersecting are the same
    county.sf.i2 <- st_transform(county.sf.i2, terra::crs(cdl.i2))
    county.sf.i2 <- county.sf.i2 %>%
      filter(GEOID == cnty.i2)
    
    ### Note: TUrn SF object into a terra object
    county.sv.i2      <- terra::vect(county.sf.i2)
    
    ### Note: Create county data to merge
    n              <- dim(county.sv.i2)[1]
    county.dt.i2   <- county.sv.i2@ptr$getDF() %>%
      as.data.table()                          %>%
      .[,ID     := seq.int(n)]
    
    ### Note: Read in CDL raster
    extract    <- terra::extract(cdl.i2, county.sv.i2, 
                                 touches = TRUE,
                                 cells   = TRUE,
                                 xy      = TRUE)  %>%
      as.data.table(key = c("ID", "Class_Names")) %>%
      merge(county.dt.i2, 
            all.x = TRUE,
            by = "ID")                            %>%
      .[,!c("ID")]
    fwrite(extract, 
           here("Data", "CDL", "Raw", paste0("CDL_County_Cellwise_Raw_", cnty.i2, "_", year.i2, ".csv")))
    return(paste0("extract for ", cnty.i2, year.i2, " complete"))
  }
  cnty.v <- county.sf.i$GEOID
  map(cnty.v, f_cnty)
  return(paste0(state.i, year.i, " complete"))
}
#map2(st.v, yr.v, f_CDLCountyCellwise)
cl          <- plan(multisession, workers = 4)
out.dist.dt <- future_mapply(f_CDLCountyCellwise, 
                             st.v,
                             yr.v,
                             future.seed = TRUE,
                             future.scheduling = FALSE)
future:::ClusterRegistry("stop")

### Step 3: Merge together the county-year cdl pixel files into a single county file
f_cntymerge <- function(st.i){
  ### Note: Get the fips code
  fips    <- c("MI", "OH", "IL", "IN", "WI",
               "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
               "26", "39", "17", "18", "55")
  fipsmat <- matrix(data = fips,
                    nrow = 5,
                    ncol = 3,
                    byrow = FALSE)
  fips   <- fipsmat[fipsmat[,2]==st.i,1]
  fips.n <- fipsmat[fipsmat[,2]==st.i,3]
  
  ### Note: Load data and change projection schemes
  fips.dt <- fips_codes %>%
    as.data.table()     %>%
    .[,FIPS_cnty := str_c(state_code, county_code)] %>%
    .[state_name %in% fipsmat[fipsmat[,2]==st.i,2]]
  cnty.v <- fips.dt$FIPS_cnty
  
  ### Note: Load in and merge together datasets wide by cell ID
  f_merge <- function(cnty.i2, 
                      st.i2 = st.i){
    
    f_bind <- function(yr.i3, 
                       cnty.i3 = cnty.i2, 
                       st.i3 = st.i2){
      dt <- fread(here("Data", "CDL", "Raw", paste0("CDL_County_Cellwise_Raw_", cnty.i3, "_", yr.i3, ".csv"))) %>%
        .[,Year := ..yr.i3]
      return(dt)
    }
    yr.v.i2 <- 2008:2018
    out.dt  <- map(yr.v.i2, f_bind) %>% rbindlist(use.names = TRUE) %>%
      .[,Year := str_c("LU_", Year)]                                 %>%
      dcast(cell + x + y + GEOID ~ Year, 
            value.var = c("Class_Names"))
    fwrite(out.dt,
           here("Data", "CDL", paste0("CDL_County_Cellwise_Raw_", cnty.i2, ".csv")))
  }
  map(cnty.v, f_merge)
}
st.v <- c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin")
#map(st.v, f_cntymerge)

cl          <- plan(multisession, workers = 4)
out.dist.dt <- future_lapply(st.v,
                             f_cntymerge,
                             future.seed = TRUE,
                             future.scheduling = FALSE)
future:::ClusterRegistry("stop")

### Step 3b: Land Use Categories to Aggregates
f_cdlcats <- function(yr.i){
  
  ### Note: Identify the location of the raw CDL raster data in the folder you moved them:
  #       for me, that is in the ./Data/Raw/CDL folder
  rasters  <- list.files(here("Data", "CDL", "Raw")) %>%
    as.data.table()                                  %>%
    .[str_detect(., ".img$")]                        %>%
    .[str_detect(., as.character(yr.i))]
  
  cdl.r <- terra::rast(here("Data", "CDL", "Raw", rasters[[1,1]]))
  cdl.v <- levels(cdl.r)
  cdl.dt <- data.table("Value" = 1:length(cdl.v[[1]]),
                       "Class_Names" = cdl.v[[1]]) %>%
    .[,Year := yr.i]
  return(cdl.dt)
}

cdl.out.dt <- map(inp1, f_cdlcats) %>% rbindlist() %>%
  dcast(Value + Class_Names ~ Year, value.var = "Class_Names") %>%
  setnames(c("Value", "Class_Names", str_c("LU_", inp1)))
fwrite(cdl.out.dt,
       here("Data", "CDL", paste0("CDL_Values_", inp1[1], "-", inp1[length(inp1)], ".csv")))

### Step 4: Create land use transition shares (no transistions for CDL necessary right now)
f_trans        <- function(st.i){
  ### Note: Get the fips code
  ### Note: Get the fips code. This crosswalks different naming structures for states
  #       because of some poor naming choices
  fips    <- c("MI", "OH", "IL", "IN", "WI",
               "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
               "26", "39", "17", "18", "55")
  fipsmat <- matrix(data = fips,
                    nrow = 5,
                    ncol = 3,
                    byrow = FALSE)
  fips   <- fipsmat[fipsmat[,2]==st.i,1]
  fips.n <- fipsmat[fipsmat[,2]==st.i,3]
  
  ### Note: Load data and change projection schemes
  fips.dt <- fips_codes %>%
    as.data.table()     %>%
    .[,FIPS_cnty := str_c(state_code, county_code)] %>%
    .[state_name %in% fipsmat[fipsmat[,2]==st.i,2]]
  
  cdl_tags.dt <- fread(here("Data", "CDL", "CDL_DRFEWSTags_2008-2018.csv")) %>%
    .[,LU_DRFEWS_CCasPasture := ifelse(LU_DRFEWS_CoverCropCheck == "CC",
                                       "Pasture",
                                       LU_DRFEWS)]
  
  ### Note: Load in each county file, summarize w/ cdl bins to create county by county transitions 
  cnty.v <- fips.dt$FIPS_cnty
  f_transmat <- function(cnty.i2, 
                         st.i2       = st.i,
                         cdl.cat.i2  = cdl_tags.dt){
    
    ### Note: Read in Cellwise Raw Data and convert to categories; save for gridcell dataset
    
    dt <- fread(here("Data", "CDL", paste0("CDL_County_Cellwise_Raw_", cnty.i2, ".csv")))
    for (yr.i2 in 2008:2018){
      dt <- dt %>%
        .[cdl.cat.i2, on = paste0("LU_", yr.i2), paste0("LU_", yr.i2, "_DRFEWS") := i.LU_DRFEWS_CCasPasture]
    }
    fwrite(dt,
           here("Data", "CDL", paste0("CDL_County_Cellwise_", cnty.i2, ".csv")))
    
    ### Note: Create a "shares" dataset for each county w/ the associated shares of land use (and acres)
    luwater.v <- names(dt)[str_detect(names(dt), "LU_")]
    share.dt <- copy(dt) %>%
      .[, if_water := apply(dt[,..luwater.v], 1, function(x) any(str_detect(x,"Water$")))]  %>%
      .[if_water == FALSE]                                                                  %>%
      .[,!c("if_water")]                                                                    %>%
      melt(measure = patterns("LU_"), 
           value.name = ("LU"), 
           variable.name = "Year")       %>%
      .[,N_Pixels := .N, by = c("Year")] %>%
      .[,N_LU   := .N, by = c("Year", "LU")]                                %>%
      .[,c("GEOID", "Year", "N_Pixels", "LU", "N_LU")]                      %>%
      unique()
    fwrite(share.dt,
           here("Data", "CDL", paste0("CDL_CountyShares_DRFEWS_", cnty.i2, ".csv")))
  }
  cl          <- plan(multisession, workers = 4)
  out.dist.dt <- future_mapply(f_transmat, 
                               cnty.v,
                               future.seed = TRUE,
                               future.scheduling = FALSE)
  future:::ClusterRegistry("stop")
}
map(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"), f_trans)

### Step 5a: Pull together shares into a single file
shares.dt <- map(here("Data", "CDL", paste0("CDL_CountyShares_DRFEWS_", fips.dt$FIPS_cnty, ".csv")), fread) %>% rbindlist()
fwrite(shares.dt,
       here("Data", "CDL", "CDL_CountyShares_DRFEWS.csv"))
                
                