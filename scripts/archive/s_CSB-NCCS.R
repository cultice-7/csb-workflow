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

#####
### Create nccs-csb intersection function that is flexible to scenario of interest
#####
nccs.crs  <- terra::rast("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\NCCS\\sfcWind\\sfcWind_day_ACCESS-CM2_ssp245_r1i1p1f1_gn_2015.nc") %>%
  terra::crs(.)
fips.dt <- fips_codes                           %>%
  as.data.table()                               %>%
  .[state %in% c("IL", "IN", "OH", "MI", "WI")] %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]

f_csb_nccs <-  function(cntyfips.i,
                        scenario.i,
                        nccs.crs.i = nccs.crs){
  csb.loc    <- "K:\\CFAES\\AEDE\\Arcdata\\BrianCultice\\INFEWS_LandUse\\Data\\USDA\\CSB"
  csb_cen.sf <- st_read(paste0(csb.loc, "\\", 
                               paste0("CSB_", cntyfips.i, "_2022_centroid_DEM-SSURGO.gpkg"))) %>%
    st_transform(., nccs.crs.i)
  datatype.v <- c("hurs", "rlds", "rsds")
  ### For each year and each data type, set up extract function; run multiproc
  f_extract_prism <- function(data.i,
                              yr.i,
                              scenario.i2     = scenario.i,
                              cntyfips.i2     = cntyfips.i,
                              csb_cen.sf.i    = csb_cen.sf){
    csb_cen.sv.i <- copy(csb_cen.sf.i) %>%
      vect()
    
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
      .[,Scenario := str_extract(.,
                                 pattern = "historical|ssp245")] %>%
      setorder(Year)                                             %>%
      .[Scenario == scenario.i2]                                 %>%
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
           paste0("K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\CSB\\", "CSB-NCCS_", cntyfips.i2, "_", scenario.i2, "_", data.i, "_", yr.i, ".csv"))
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
map2(fips.dt$FIPS_cnty,
     rep("ssp245", times = length(fips.dt$FIPS_cnty)), f_csb_nccs)