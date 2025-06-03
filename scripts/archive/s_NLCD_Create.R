######################################################################
##### NLCD - Create 
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
               "furrr", 
               "pryr",
               "rmapshaper",
               "tigris")
pacman::p_load(char = libraries)

### Step 1: Crop raw NLCD rasters into smaller, regional pieces for running intersections
v.region <- c("Ohio", "Indiana", "Illinois", "Michigan", "Wisconsin")
v.years  <- c("2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019", "2021")
NLCD_CropSave <- function(i.year,
                          i.region = v.region,
                          i.regionname = "GreatLakes"){
  
  ### Note: Identify the location of the raw NLCD raster data in the folder you moved them:
  #       for me, that is in the ./Data/Raw/NLCD folder
  rasters  <- list.files(here("Data", "NLCD", "Raw")) %>%
    as.data.table()                   %>%
    .[str_detect(., ".img$")]         %>%
    .[str_detect(., i.year)]
  
  ### Note: Use that location identified above for year year_v and load in the 
  #       full raster data
  rast.nlcd <- terra::rast(here("Data", "NLCD", "Raw", rasters[[1,1]]))
  crs.nlcd  <- terra::crs(rast.nlcd)
  
  ### Note: Load in a US state shapefile; I use this to filter a specific state,
  #       then intersect that specific state w/ the raster to create a state specific
  #       version of the CDL data. Be sure the CRS is the same for both geospatial data
  state.loc   <- here("Data", "Census", "US_state_2018.shp")
  sf.region  <- st_read(state.loc)     %>%
    .[,c("NAME", "GEOID", "geometry")] %>%
    filter(NAME %in% i.region)         %>%
    st_transform(., crs.nlcd)          %>%
    st_union()
  
  sv.region <- terra::vect(sf.region)
    
  ### Note: Cut out land use from stack (crop first to make this go more quickly)
  rast.nlcd_crop  <- terra::crop(rast.nlcd, sv.region)
  rast.nlcd_crop[is.nan(rast.nlcd_crop)] <- 0
    
  ### Note: Save the cropped file
  tempname <- here("Data", "NLCD", "Raw", paste0("NLCD_", i.regionname, "_", i.year, ".tif"))
  terra::writeRaster(rast.nlcd_crop, tempname, 
                     overwrite = TRUE)
}
map(v.years, NLCD_CropSave)


### Step 2: Intersect Counties w/ CDL. This grabs all the pixels that fall within a given
#         county. For my purposes, I'm doing this for these 5 states. You can make this a vector
#         of all states as needed.
st.v <- rep(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"),
            times = length(inp1))
yr.v <- rep(inp1,
            each = 5)
f_NLCDCountyCellwise <- function(state.i, year.i){
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
    filter(STATEFP == fips.n) 
  
  f_cnty <- function(cnty.i2, 
                     county.sf.i2 = county.sf.i, 
                     state.i2     = state.i, 
                     year.i2     = year.i){
    
    ### Note: Load in raster file
    NLCD.i2  <- terra::rast(here("Data", "NLCD", "Raw", paste0("NLCD_", state.i2, "_", year.i2, ".tif")))
    
    ### Note: Confirm that coordinate reference systems for the two files im intersecting are the same
    county.sf.i2 <- st_transform(county.sf.i2, terra::crs(NLCD.i2))
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
    extract    <- terra::extract(NLCD.i2, county.sv.i2, 
                                 touches = TRUE,
                                 cells   = TRUE,
                                 xy      = TRUE)  %>%
      as.data.table(key = c("ID", "NLCD Land Cover Class"))   %>%
      merge(county.dt.i2, 
            all.x = TRUE,
            by = "ID")                            %>%
      .[,!c("ID")]
    fwrite(extract, 
           here("Data", "NLCD", "Raw", paste0("NLCD_County_Cellwise_Raw_", cnty.i2, "_", year.i2, ".csv")))
    return(paste0("extract for ", cnty.i2, year.i2, " complete"))
  }
  cnty.v <- county.sf.i$GEOID
  map(cnty.v, f_cnty)
  return(paste0(state.i, year.i, " complete"))
}
#map2(st.v, yr.v, f_NLCDCountyCellwise)
cl          <- plan(multisession, workers = 4)
out.dist.dt <- future_mapply(f_NLCDCountyCellwise, 
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
      dt <- fread(here("Data", "NLCD", "Raw", paste0("NLCD_County_Cellwise_Raw_", cnty.i3, "_", yr.i3, ".csv"))) %>%
        .[,Year := ..yr.i3]
      return(dt)
    }
    yr.v.i2 <- c("2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019")
    out.dt  <- map(yr.v.i2, f_bind) %>% rbindlist(use.names = TRUE) %>%
      .[,Year := str_c("LU_", Year)]                                 %>%
      dcast(cell + x + y + GEOID + CSAFP + CBSAFP ~ Year, 
            value.var = c("NLCD Land Cover Class"))
    fwrite(out.dt,
           here("Data", "CDL", paste0("NLCD_County_Cellwise_Raw_", cnty.i2, ".csv")))
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
f_nlcdcats <- function(yr.i){
  nlcd.r <- terra::rast(here("Data", "NLCD", "Raw", paste0("nlcd_", yr.i, "_land_cover_l48_20210604.img")))
  nlcd.v <- levels(nlcd.r)
  nlcd.dt <- data.table("Value" = 1:length(nlcd.v[[1]]),
                        "NLCD_Cat" = nlcd.v[[1]]) %>%
    .[,Year := yr.i]
  return(nlcd.dt)
}

# Note: New Classifications based on our aggregate categories
inp1        <- c("2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019")
lu.v        <- c("Other", "Cropland", "Forest", "Developed", 
                 "Developed", "Developed", "Developed", "Wetlands", "Forest", 
                 "Hay/Pasture", "Grasslands", "Forest", "Water", 
                 "Other", "Grasslands", "Other", "Wetlands")
lu2.v       <- c("Other", "Cropland", "Forest", "Developed", 
                 "Developed", "Developed", "Developed", "Wetlands", 
                 "Forest", "Pasture", "Pasture", "Forest", 
                 "Water", "Other", "Pasture", "Other", "Wetlands")

nlcd.out.dt <- map(inp1, f_nlcdcats) %>% rbindlist() %>%
  .[NLCD_Cat != ""]                  %>%
  dcast(NLCD_Cat ~ Year, value.var = "Value") %>%
  .[,Summary_Cat := ..lu.v]                   %>%
  .[,DRFEWS_Cat  := ..lu2.v]                  %>%
  .[,c("NLCD_Cat", "2001", "Summary_Cat", "DRFEWS_Cat")] %>%
  setnames("2001", "Value")

### Step 4: Create land use transition matrices for states.
f_trans_drfews <- function(st.i,
                           nlcd.cat.i = nlcd.out.dt){
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
  
  ### Note: Load in each county file, summarize w/ cdl bins to create county by county transitions 
  cnty.v <- fips.dt$FIPS_cnty
  f_transmat <- function(cnty.i2, 
                         st.i2       = st.i,
                         nlcd.cat.i2 = nlcd.cat.i){
    
    ### Note: Read in Cellwise Raw Data and convert to categories; save for gridcell dataset
    dt <- fread(here("Data", "CDL", paste0("NLCD_County_Cellwise_Raw_", cnty.i2, ".csv")))               %>%
      .[.(LU_2001 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2001", LU_2001 := i.to] %>%
      .[.(LU_2004 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2004", LU_2004 := i.to] %>%
      .[.(LU_2006 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2006", LU_2006 := i.to] %>%
      .[.(LU_2008 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2008", LU_2008 := i.to] %>%
      .[.(LU_2011 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2011", LU_2011 := i.to] %>%
      .[.(LU_2013 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2013", LU_2013 := i.to] %>%
      .[.(LU_2016 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2016", LU_2016 := i.to] %>%
      .[.(LU_2019 = nlcd.cat.i2$NLCD_Cat, to = nlcd.cat.i2$DRFEWS_Cat), on = "LU_2019", LU_2019 := i.to] %>%
      .[LU_2001 != "Water" & LU_2004 != "Water" & LU_2006 != "Water" & LU_2008 != "Water" & LU_2011 != "Water" & LU_2013 != "Water" & LU_2016 != "Water" & LU_2019 != "Water"]
    fwrite(dt,
           here("Data", "NLCD", paste0("NLCD_County_Cellwise_", cnty.i2, ".csv")))
    
    ### Note: Create a "shares" dataset for each county w/ the associated shares of land use (and acres)
    share.dt <- dt                       %>%
      melt(measure = patterns("LU_"), 
           value.name = ("LU"), 
           variable.name = "Year")       %>%
      .[,N_Pixels := .N, by = c("Year")] %>%
      .[,SumVal   := .N, by = c("Year", "LU")]                              %>%
      .[,c("GEOID", "CSAFP", "CBSAFP", "Year", "N_Pixels", "LU", "SumVal")] %>%
      unique()
    fwrite(share.dt,
           here("Data", "NLCD", paste0("NLCD_CountyShares_DRFEWS_", cnty.i2, ".csv")))
    
    ### Note: Transitions Dataset (e.g. i to j transitions)
    f_transfunc <- function(index.i3){
      yr1 <- yr.v.i3[index.i3]
      yr2 <- yr.v.i3[index.i3+1]
      
      var.i3 <- str_c(dt[,get(paste0("LU_", yr1))], "-", dt[,get(paste0("LU_", yr2))])
      dt <- dt %>%
        .[,temp := var.i3] %>%
        setnames("temp",
                 paste0("LUC_", yr1, "-", yr2))
      return("Done")
    }
    yr.v.i3 <- c("2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019")
    map(1:(length(yr.v.i3)-1), f_transfunc)
    
    dt <- dt                              %>%
      .[, grep("LU_", names(dt)) := NULL] %>%
      melt(measure = patterns("LUC"), 
           value.name = ("LUC"), 
           variable.name = "Year")        %>%
      .[!str_detect(LUC, "Other")]        %>%
      .[,N_Pixels := .N, by = "Year"]     %>%
      .[,SumVal   := .N, by = c("Year", "LUC")]                              %>%
      .[,c("GEOID", "CSAFP", "CBSAFP", "Year", "N_Pixels", "LUC", "SumVal")] %>%
      unique()
    fwrite(dt,
           here("Data", "NLCD", paste0("NLCD_CountyLUC_DRFEWS_", cnty.i2, ".csv")))
    return(paste0("Transitions for ", cnty.i2, " done."))
  }
  cl          <- plan(multisession, workers = 4)
  out.dist.dt <- future_mapply(f_transmat, 
                               cnty.v,
                               future.seed = TRUE,
                               future.scheduling = FALSE)
  future:::ClusterRegistry("stop")
}
map(c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin"), f_trans_drfews)

### Step 5: Pull Together Transitions into Single File
fips.dt <- fips_codes %>%
  as.data.table()     %>%
  .[state_name %in% c("Michigan", "Ohio", "Illinois", "Indiana", "Wisconsin")] %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]
trans.dt <- map(here("Data", "NLCD", paste0("NLCD_CountyLUC_DRFEWS_", fips.dt$FIPS_cnty, ".csv")), fread) %>% rbindlist()
fwrite(trans.dt,
       here("Data", "NLCD", "NLCD_CountyLUC_DRFEWS.csv"))

### Step 5b: Pull together shares into a single file
shares.dt <- map(here("Data", "NLCD", paste0("NLCD_CountyShares_DRFEWS_", fips.dt$FIPS_cnty, ".csv")), fread) %>% rbindlist()
fwrite(shares.dt,
       here("Data", "NLCD", "NLCD_CountyShares_DRFEWS.csv"))
                
                