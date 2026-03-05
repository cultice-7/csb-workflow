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
states <- c("Ohio", "Indiana", "Illinois", "Michigan", "Wisconsin")
inp1   <- as.character(rep(2008:2018, each = 5))
inp2   <- rep(states, times = length(2008:2018))
CDL_CropSave <- function(year_v, state_v){
  
  ### Note: Identify the location of the raw CDL raster data in the folder you moved them:
    #       for me, that is in the ./Data/CDL folder
  rasters  <- list.files(here("Data", "CDL")) %>%
    as.data.table()                   %>%
    .[str_detect(., ".img$")]         %>%
    .[!str_detect(., ".vat")]         %>%
    .[str_detect(., year_v)]          %>%
    .[str_detect(., "30m")]
  
  ### Note: Use that location identified above for year year_v and load in the 
    #       full raster data
  CDL     <- terra::rast(here("Data", "CDL", rasters[[1,1]]))
  CDL_crs <- terra::crs(CDL)
  
  ### Note: Load in a US state shapefile; I use this to filter a specific state,
    #       then intersect that specific state w/ the raster to create a state specific
    #       version of the CDL data. Be sure the CRS is the same for both geospatial data
  state_loc <- here("Data", "Census", "US_state_2018.shp")
  state_sf  <- st_read(state_loc) # read it in
  state_sf2 <- state_sf[state_sf$NAME == state_v
                        ,c("NAME")] # filter it
  state_sf3  <- st_transform(state_sf2, crs(CDL)) # transform the crs
  state_sv   <- terra::vect(state_sf3) # create a "terra" version to use with the raster data package
  rm(state_sf, state_sf2, state_sf3) 
  
  ### Note: Cut out land use from stack (crop first to make this go more quickly)
  CDL_t  <- terra::crop(CDL, state_sv)
  CDL_t[is.nan(CDL_t)] <- 0
  
  ### Note: Save the cropped file
  tempname <- here("CDL Data", paste0("CDL_", state_v, "_", year_v, ".tif"))
  terra::writeRaster(CDL_t, tempname, 
                     overwrite = TRUE)
}
map2(inp1, inp2, CDL_CropSave)


### Step 2: Intersect Counties w/ CDL. This grabs all the pixels that fall within a given
  #         county. For my purposes, I'm doing this for these 5 states. You can make this a vector
  #         of all states as needed.
st.v <- rep(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"),
            times = length(2008:2018))
yr.v <- rep(2008:2018,
            each = 5)
f_CDLCountyCellwise <- function(st, yr){
  ### Note: Get the fips code. This crosswalks different naming structures for states
    #       because of some poor naming choices
  fips    <- c("MI", "OH", "IL", "IN", "WI",
               "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
               "26", "39", "17", "18", "55")
  fipsmat <- matrix(data = fips,
                    nrow = 5,
                    ncol = 3,
                    byrow = FALSE)
  fips <- fipsmat[fipsmat[,2]==st,1]
  
  ### Note: Load data and change projection schemes
  county_sf  <- st_read(here("Data", "Census", paste0(fips, "_County_NU.shp")))
  
  f_cnty <- function(cnty, county_sf.inp = county_sf, st.inp = st, yr.inp = yr){
    
    ### Note: Load in raster file
    CDL            <- terra::rast(here("Data", "CDL", paste0("CDL_", st.inp, "_", yr.inp, ".tif")))
    
    ### Note: Confirm that coordinate reference systems for the two files im intersecting are the same
    county_sf2 <- st_transform(county_sf.inp, terra::crs(CDL))
    county_sf2 <- county_sf2 %>%
      filter(CountyID == cnty)
    
    ### Note: TUrn SF object into a terra object
    county_sv      <- terra::vect(county_sf2)
    
    ### Note: Create county data to merge
    n           <- dim(county_sv)[1]
    county_dt   <- county_sv@ptr$getDF() %>%
      as.data.table()                    %>%
      .[,ID     := seq.int(n)]
      
      ### Note: Read in CDL raster
    extract    <- terra::extract(CDL, county_sv, 
                                 touches = TRUE,
                                 cells   = TRUE,
                                 xy      = TRUE)  %>%
      as.data.table(key = c("ID", "Layer_1"))     %>%
      merge(county_dt, 
            all.x = TRUE,
            by = "ID")                            %>%
      .[,!c("ID")]
    fwrite(extract, 
           here("Data", "CDL", paste0("CountyLU_Cellwise_Raw_", st.inp, "_", cnty, "_", yr.inp, ".csv")))
  }
  cnty.v <- county_sf$CountyID
  map(cnty.v, f_cnty)
}
map2(st.v, yr.v, f_CDLCountyCellwise)

### Step 3: Merge together the county-year cdl pixel files into a single county file
f_cntymerge <- function(st){
  ### Note: Get the fips code
  fips    <- c("MI", "OH", "IL", "IN", "WI",
               "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
               "26", "39", "17", "18", "55")
  fipsmat <- matrix(data = fips,
                    nrow = 5,
                    ncol = 3,
                    byrow = FALSE)
  fips <- fipsmat[fipsmat[,2]==st,1]
  
  ### Note: Load data and change projection schemes
  county_sf  <- st_read(here("Data", "Census", paste0(fips, "_County_NU.shp")))
  cnty.v     <- county_sf$CountyID
  
  ### Note: Load in and merge together datasets wide by cell ID
  f_merge <- function(cnty, st.inp = st){
    
    f_bind <- function(yr, cnty.inp = cnty, st.inp2 = st.inp){
      dt <- fread(here("Data", "CDL", paste0("CountyLU_Cellwise_Raw_", st.inp2, "_", cnty.inp, "_", yr, ".csv"))) %>%
        .[,Year := ..yr]
      return(dt)
    }
    out.dt <- map(2008:2018, f_bind) %>% rbindlist(use.names = TRUE) %>%
      .[,Year := str_c("LU_", Year)]                                 %>%
      dcast(cell + x + y + StateID + CountyID + FIPS + CountyName + CSAID + CBSAID ~ Year, 
            value.var = c("Layer_1"))
    fwrite(out.dt,
           here("Data", "CDL", paste0("CountyLU_Cellwise_", st.inp, "_", cnty, ".csv")))
  }
  map(cnty.v, f_merge)
}
map(st.v, f_cntymerge)

### Step 4: Create land use transition matrices for states. I'm not sure if you need to do
  #         this part or not. 
f_transitions <- function(st){
  ### Note: Get the fips code
  fips    <- c("MI", "OH", "IL", "IN", "WI",
               "Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin",
               "26", "39", "17", "18", "55")
  fipsmat <- matrix(data = fips,
                    nrow = 5,
                    ncol = 3,
                    byrow = FALSE)
  fips <- fipsmat[fipsmat[,2]==st,1]
  
  ### Note: Load data and change projection schemes
  county_sf  <- st_read(here("Data", "Census", paste0(fips, "_County_NU.shp")))
  
  ### Note: Load Transition Tags
  bins.c <- fread(here("Data", "CDL", "cdl_tags_cropland.csv"))
  bins.p <- fread(here("Data", "CDL", "cdl_tags_pasture.csv"))
  bins.u <- fread(here("Data", "CDL", "cdl_tags_urban.csv"))
  bins.f <- fread(here("Data", "CDL", "cdl_tags_forest.csv"))
  bins.w <- fread(here("Data", "CDL", "cdl_tags_wetlands.csv"))
  bins   <- rbindlist(list(bins.c, bins.p, bins.u, bins.f, bins.w)) %>%
    .[,!c("Label")]
  fwrite(bins, here("Data", "CDL", "cdl_tags.csv"))
  
  ### Note: Load in each county file, summarize w/ cdl bins to create county by county transitions 
  cnty.v <- county_sf$CountyID
  f_transmat <- function(cnty, st.inp = st, bins.inp = bins){
    
    ### Note: Read in Cellwise Raw Data
    dt <- fread(here("Data", "CDL", paste0("CountyLU_Cellwise_", st.inp, "_", cnty, ".csv")))
    
    cntyid <- dt %>%
      .[,c("StateID", "CountyID", "FIPS", "CSAID", "CBSAID")] %>%
      unique()
    
    dt <- dt %>%
      .[,N_Pixels := dim(.)[1]]                                                               %>%
      .[.(LU_2008 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2008", LU_2008 := i.to] %>%
      .[.(LU_2009 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2009", LU_2009 := i.to] %>%
      .[.(LU_2010 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2010", LU_2010 := i.to] %>%
      .[.(LU_2011 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2011", LU_2011 := i.to] %>%
      .[.(LU_2012 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2012", LU_2012 := i.to] %>%
      .[.(LU_2013 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2013", LU_2013 := i.to] %>%
      .[.(LU_2014 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2014", LU_2014 := i.to] %>%
      .[.(LU_2015 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2015", LU_2015 := i.to] %>%
      .[.(LU_2016 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2016", LU_2016 := i.to] %>%
      .[.(LU_2017 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2017", LU_2017 := i.to] %>%
      .[.(LU_2018 = bins.inp$Layer_1, to = bins.inp$LU_Ind), on = "LU_2018", LU_2018 := i.to] %>%
      .[,':=' (LUC_0809 = str_c(LU_2008, LU_2009),
               LUC_0910 = str_c(LU_2009, LU_2010),
               LUC_1011 = str_c(LU_2010, LU_2011),
               LUC_1112 = str_c(LU_2011, LU_2012),
               LUC_1213 = str_c(LU_2012, LU_2013),
               LUC_1314 = str_c(LU_2013, LU_2014),
               LUC_1415 = str_c(LU_2014, LU_2015),
               LUC_1516 = str_c(LU_2015, LU_2016),
               LUC_1617 = str_c(LU_2016, LU_2017),
               LUC_1718 = str_c(LU_2017, LU_2018),
               LUC_0813 = str_c(LU_2008, LU_2013),
               LUC_1318 = str_c(LU_2013, LU_2018))]          %>%
      .[,LUC_0809_N := .N/N_Pixels, by = "LUC_0809"]         %>%
      .[,LUC_0910_N := .N/N_Pixels, by = "LUC_0910"]         %>%
      .[,LUC_1011_N := .N/N_Pixels, by = "LUC_1011"]         %>%
      .[,LUC_1112_N := .N/N_Pixels, by = "LUC_1112"]         %>%
      .[,LUC_1213_N := .N/N_Pixels, by = "LUC_1213"]         %>%
      .[,LUC_1314_N := .N/N_Pixels, by = "LUC_1314"]         %>%
      .[,LUC_1415_N := .N/N_Pixels, by = "LUC_1415"]         %>%
      .[,LUC_1516_N := .N/N_Pixels, by = "LUC_1516"]         %>%
      .[,LUC_1617_N := .N/N_Pixels, by = "LUC_1617"]         %>%
      .[,LUC_1718_N := .N/N_Pixels, by = "LUC_1718"]         %>%
      .[,LUC_0813_N := .N/N_Pixels, by = "LUC_0813"]         %>%
      .[,LUC_1318_N := .N/N_Pixels, by = "LUC_1318"]
    
    f_pull <- function(yr, dt.inp = dt){
      vars <- c(paste0("LUC_", yr), paste0("LUC_", yr, "_N"))
      dt.out <- dt.inp %>%
        .[,..vars]     %>%
        .[,Year_trans := ..yr] %>% 
        setnames(c(paste0("LUC_", yr), paste0("LUC_", yr, "_N")),
                 c("Transition", "Share"))      %>%
        unique()
      return(dt.out)
    }
    out.l <- map(c("0809", "0910", "1011", "1112", "1213",
                   "1314", "1415", "1516", "1617", "1718",
                   "0813", "1318"),
                 f_pull)             %>%
      rbindlist()                    %>%
      .[,StateID  := cntyid$StateID] %>%
      .[,CountyID := ..cnty]         %>%
      .[,FIPS     := cntyid$FIPS]    %>%
      .[,CSAID    := cntyid$CSAID]   %>%
      .[,CBSAID   := cntyid$CBSAID]  %>%
      .[str_length(Transition)==2]
    return(out.l)
  }
  out.l2 <- map(cnty.v, f_transmat) %>%
    rbindlist()                     %>%
    fwrite(here("Data", "CDL", paste0("CountyLU_Cellwise_Trans_", st, ".csv")))
}
map(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"), f_transitions)


### Note: merge labels for clarity
f_labels <- function(st){
  cdl <- fread(here("Data", "CDL", "cdl_tags.csv")) %>%
    .[,c("LU", "LU_Ind")]                           %>%
    unique()
  dt  <- fread(here("Data", "CDL", paste0("CountyLU_Cellwise_Trans_", st, ".csv"))) %>%
    .[,':=' (trans1 = as.numeric(str_sub(as.character(Transition), start = 1, end = 1)),
             trans2 = as.numeric(str_sub(as.character(Transition), start = 2)))] %>%
    merge(cdl,
          by.x = "trans1",
          by.y = "LU_Ind",
          all.x = TRUE) %>%
    merge(cdl,
          by.x = "trans2",
          by.y = "LU_Ind",
          all.x = TRUE) %>%
    setnames(c("LU.x", "LU.y"),
             c("LU_Year1", "LU_Year2")) %>%
    .[,Transition_Names := str_c(LU_Year1, "-", LU_Year2)]
  fwrite(dt, here("Data", "CDL", paste0("CountyLU_Cellwise_Trans_edit_", st, ".csv")))
}
map(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"), f_labels)
