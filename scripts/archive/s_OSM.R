################################################################################
##### OSM Download and Shapefile Generation
##### 4-4-23
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
               "r5r",
               "osmextract")
pacman::p_load(char = libraries)

#####
### Load County and MSA files
#####

cnty.dt              <- tigris::counties(year = 2011)
st_geometry(cnty.dt) <- NULL
cnty.dt              <- cnty.dt %>%
  as.data.table()               %>%
  .[,FIPS_cnty := str_c(STATEFP, COUNTYFP, sep = "")]

cnty.sf <- tigris::counties(year = 2011) %>%
  st_transform(., crs = 3857)

msa.dt              <- tigris::core_based_statistical_areas(year = 2011)
st_geometry(msa.dt) <- NULL
msa.dt              <- msa.dt %>%
  as.data.table()
msa.v <- c("16980", "35620", "17460", "31100", "33460", "12580", "37980", "41180", "38300", "45300")

### Step 2: Download all pbf files for all states
# f_dlosm <- function(st.i){
#   details <- oe_match(st.i, match_by = "iso3166_2")
#   
#   oe_download(
#     file_url = details$url, 
#     file_size = details$file_size,
#     provider = "test",
#     download_directory = here("Data", "OSM")
#   )
#   return(paste0("Download for ", st.i, "complete. Big baller doing big shit"))
# }
# file.v <- oe_match_pattern("US", match_by = "iso3166_2")$geofabrik
# map(file.v, f_dlosm)
# map(list.files(here("Data", "OSM"), 
#                pattern = ".pbf",
#                full.names = TRUE), 
#     oe_vectortranslate)
# 
# ### Step 2a: Pull out the shapefiles (e.g. land use map?) from OSM
# map(list.files(here("Data", "OSM"),
#                pattern = ".pbf",
#                full.names = TRUE),
#     layer = "multipolygons",
#     oe_vectortranslate, 
# )

f_osmpark <- function(msa.i,
                      cnty.dt.i = cnty.dt,
                      cnty.sf.i = cnty.sf){
  st.v <- cnty.dt.i %>%
    .[CBSAFP %in% msa.i]
  st.v    <- unique(st.v$STATEFP)
  st.dt.i <- fips_codes     %>%
    as.data.table()         %>%
    .[state_code %in% st.v] %>%
    .[, state_name := str_to_lower(str_replace_all(state_name, " ", "-"))]
  
  ### Note: Land use/leisure definitions to keep for park
  keep.v <- c("park", "nature_reserve", "nature_reserve;park", "playground", "dog_park", "recreation_ground")
  f_loadosm <- function(stname.i2,
                        keep.v.i = keep.v){
    osm.sf.i2 <- st_read(here("Data", "OSM", paste0("test_", stname.i2, "-latest.gpkg")), layer = "multipolygons") %>%
      filter(leisure %in% keep.v.i)
  }
  osm.sf.i <- map(unique(st.dt.i$state_name), f_loadosm) %>% rbindlist() %>% st_as_sf() %>%
    st_transform(., 3857)                                                               %>%
    st_buffer(., dist = 0)
  
  ### Note: filter for msa's, intersect, save in msa specific folder
  cnty.sf.i <- cnty.sf.i    %>%
    filter(CBSAFP == msa.i) %>%
    st_union()
  
  osm.sf.int.i <- st_intersection(osm.sf.i, cnty.sf.i)
  
  st_write(osm.sf.int.i, here("Data", "OSM", msa.i, paste0("OSM_Parks_", msa.i, ".gpkg")),
           append = FALSE)
  
  ### Note: Return flag
  return(paste0("Intersection for", msa.i, " Complete"))
}
map(msa.v, f_osmpark)

f_parks_OSM <- function(msa.i,
                        dldate.i = "220429"){
  ### Note: Load sales data for msa; create geometry
  sales.sf <- fread(here("Data", paste0("Trans_Tract_MSA", msa.i, "-", dldate.i, ".csv"))) %>%
    .[,tempid := .I]                                                                       %>%
    st_as_sf(., coords = c("lon", "lat"), crs = st_crs("ESRI:102003"))                     %>%
    st_transform(., crs = 5070)
  
  ### Note: Load OSM Parks data
  osm.sf <- st_read(here("Data", "OSM", msa.i, paste0("OSM_Parks_", msa.i, ".gpkg")))      %>%
    st_transform(., crs = 5070) %>%
    .[!is.na(.$name),]          %>%
    mutate(parkarea_acres = as.numeric(st_area(.))*0.000247105)
  
  osm.sf_halfAcre <- copy(osm.sf) %>%
    filter(parkarea_acres > 0.5)
  
  ### Note: Load in county shapefile and filter for msa to calculate land available
  if (msa.i == "31100"){
    msa.i_fix <- "31080"
  }else{
    msa.i_fix <- msa.i
  }
  cnty.sf <- st_read(here("Data", "Census", "US_county_2018.shp"))      %>%
    filter(CBSAFP == msa.i_fix)                                         %>%
    st_transform(., crs = 5070)
  
  nearest.v          <- st_nearest_feature(sales.sf, osm.sf)
  nearest.v_halfAcre <- st_nearest_feature(sales.sf, osm.sf_halfAcre)
  f_dist <- function(row.i,
                     index.i           = nearest.v,
                     index_ha.i        = nearest.v_halfAcre,
                     osm.sf.i          = osm.sf,
                     osm.sf_halfAcre.i = osm.sf_halfAcre,
                     sales.sf.i = sales.sf){
    dist.i    <- as.numeric(st_distance(sales.sf.i[row.i,], osm.sf.i[index.i[row.i],]))
    dist_ha.i <- as.numeric(st_distance(sales.sf.i[row.i,], osm.sf_halfAcre.i[index_ha.i[row.i],]))
    size.i    <- osm.sf.i$parkarea_acres[index.i[row.i]]
    size_ha.i <- osm.sf_halfAcre.i$parkarea_acres[index_ha.i[row.i]]
    type.i    <- osm.sf.i$leisure[index.i[row.i]]
    type_ha.i <- osm.sf_halfAcre.i$leisure[index_ha.i[row.i]]
    dt        <- data.table("tempid"    = sales.sf.i[row.i,]$tempid,
                            "Nearest_m" = dist.i,
                            "Nearest_hAcre_m" = dist_ha.i,
                            "Nearest_ParkSize_Acre" = size.i,
                            "Nearest_ParkSize_hAcre_Acre" = size_ha.i,
                            "Nearest_Type" = type.i,
                            "Nearest_Type_hAcre" = type_ha.i)
    return(dt)
  }
  cl          <- plan(multisession, workers = 8)
  out.dist.dt <- future_lapply(1:dim(sales.sf)[1], f_dist,
                               future.seed = TRUE) %>%
    rbindlist()
  future:::ClusterRegistry("stop")
  
  #1:dim(sales.sf)[1]
  ### Note: Buffer zones calculation (surrounding green space)
  f_buffers <- function(buffer,
                        osm.sf.i   = osm.sf,
                        sales.sf.i = sales.sf,
                        cnty.sf.i  = cnty.sf){
    
    ### Step 1: Create Buffers for Sales
    sales_b.sf.i <- st_buffer(sales.sf.i, dist = buffer)
    
    ### Step 2: Intersect Buffers w/ County to Identify total Land
    l_area.sf.i  <- st_intersection(sales_b.sf.i, cnty.sf.i) %>%
      mutate(Area_land = st_area(.))
    st_geometry(l_area.sf.i) <- NULL
    l_area.sf.i  <- l_area.sf.i                       %>%
      as.data.table()                                 %>%
      .[, Area_land := as.numeric(sum(Area_land)), by = "tempid"] %>%
      .[,c("tempid", "Area_land")]                                %>%
      unique()
    
    ### Step 3: Intersect Buffers w/ OSM for Park Land
    osm.sf.i      <- st_union(osm.sf.i)
    osm_area.sf.i <- st_intersection(sales_b.sf.i, osm.sf.i) %>%
      mutate(Area = as.numeric(st_area(.)))
    
    ### Step 4: Calculate Share of Park/Greenspace for each sale
    st_geometry(osm_area.sf.i) <- NULL
    osm_area.sf.i <- osm_area.sf.i %>%
      as.data.table()              %>%
      .[, Area_t := as.numeric(sum(Area)), by = "tempid"] %>%
      .[, c("tempid", "Area_t")]                          %>%
      unique()
    
    ### Step 5: Merge to counties and change missing to zero (eg no parks )
    l_area.sf.i <- copy(l_area.sf.i)                %>%
      merge(osm_area.sf.i,
            by    = "tempid",
            all.x = TRUE)                     %>%
      .[is.na(Area_t), Area_t := 0]                       %>%
      .[,parkshare   := Area_t/as.numeric(Area_land)]     %>%
      .[,buffer_size := ..buffer]                         %>%
      .[,.(tempid, parkshare, buffer_size)]               %>%
      setorder(tempid)                                    %>%
      .[parkshare > 1, parkshare := 1]
    return(l_area.sf.i)
  }
  buff.v <- c(100, 200, 400, 800, 1600)
  cl     <- plan(multisession, workers = 3)
  out.buff.dt <- future_lapply(buff.v, f_buffers,
                               future.seed = TRUE,
                               future.scheduling = FALSE) %>% rbindlist(fill = TRUE) %>%
    dcast(tempid ~ buffer_size, value.var = c("parkshare"))
  future:::ClusterRegistry("stop")
  
  ### Note restore coords and cbind park share and distance to sales
  out.full.dt <- copy(out.buff.dt) %>%
    setnames(c("tempid", str_c("Buffer_",  c(100, 200, 400, 800, 1600)))) %>%
    merge(out.dist.dt,
          by = "tempid")                                                  %>%
    .[,FIPS_cbsa := msa.i]
  
  # sales.sf <- sales.sf %>%
  #   st_transform(., 4326)
  # lon.v <- unlist(purrr::transpose(sales.sf$geometry)[[1]])
  # lat.v <- unlist(purrr::transpose(sales.sf$geometry)[[2]])
  # 
  # st_geometry(sales.sf) <- NULL
  # sales.sf <- sales.sf     %>%
  #   as.data.table()        %>%
  #   .[,':=' (lon = ..lon.v,
  #            lat = ..lat.v,
  #            ParkDist_m = out.dist.dt$NearestPark_m)] %>%
  #   merge(out.buff.dt,
  #         by = "tempid")
  
  fwrite(out.full.dt, 
         here("Data", "OSM", paste0("Trans_OSMDistance_CBSA", msa.i, ".csv")))
  return(paste0("int for ", msa.i, " complete"))
}
o_raw.dt <- map(msa.v[-1], f_parks_OSM)

### Merge together all files into single distance table
files.v <- list.files(here("Data", "OSM"),
                      pattern = "OSMDistance_CBSA",
                      full.names = TRUE)
full.dt <- lapply(files.v, fread) %>% rbindlist()
fwrite(full.dt,
       here("Data", "OSM", "OSMDistance_CBSA_allMSA.csv"))

#############################################################
### Addendum: Dist function just for nearest park with size
###           and type merge
#############################################################

f_parks_OSM_dist <- function(msa.i,
                             dldate.i = "220429"){
  ### Note: Load sales data for msa; create geometry
  sales.sf <- fread(here("Data", paste0("Trans_Tract_MSA", msa.i, "-", dldate.i, ".csv"))) %>%
    .[,tempid := .I]                                                                       %>%
    st_as_sf(., coords = c("lon", "lat"), crs = st_crs("ESRI:102003"))                     %>%
    st_transform(., crs = 5070)
  
  ### Note: Load OSM Parks data
  osm.sf <- st_read(here("Data", "OSM", msa.i, paste0("OSM_Parks_", msa.i, ".gpkg")))      %>%
    st_transform(., crs = 5070) %>%
    .[!is.na(.$name),]          %>%
    mutate(parkarea_acres = as.numeric(st_area(.))*0.000247105)
  
  osm.sf_halfAcre <- copy(osm.sf) %>%
    filter(parkarea_acres > 0.5)
  
  ### Note: Load in county shapefile and filter for msa to calculate land available
  if (msa.i == "31100"){
    msa.i_fix <- "31080"
  }else{
    msa.i_fix <- msa.i
  }
  cnty.sf <- st_read(here("Data", "Census", "US_county_2018.shp"))      %>%
    filter(CBSAFP == msa.i_fix)                                         %>%
    st_transform(., crs = 5070)
  
  nearest.v          <- st_nearest_feature(sales.sf, osm.sf)
  nearest.v_halfAcre <- st_nearest_feature(sales.sf, osm.sf_halfAcre)
  f_dist <- function(row.i,
                     index.i           = nearest.v,
                     index_ha.i        = nearest.v_halfAcre,
                     osm.sf.i          = osm.sf,
                     osm.sf_halfAcre.i = osm.sf_halfAcre,
                     sales.sf.i = sales.sf){
    dist.i    <- as.numeric(st_distance(sales.sf.i[row.i,], osm.sf.i[index.i[row.i],]))
    dist_ha.i <- as.numeric(st_distance(sales.sf.i[row.i,], osm.sf_halfAcre.i[index_ha.i[row.i],]))
    size.i    <- osm.sf.i$parkarea_acres[index.i[row.i]]
    size_ha.i <- osm.sf_halfAcre.i$parkarea_acres[index_ha.i[row.i]]
    type.i    <- osm.sf.i$leisure[index.i[row.i]]
    type_ha.i <- osm.sf_halfAcre.i$leisure[index_ha.i[row.i]]
    dt        <- data.table("tempid"    = sales.sf.i[row.i,]$tempid,
                            "Nearest_m" = dist.i,
                            "Nearest_hAcre_m" = dist_ha.i,
                            "Nearest_ParkSize_Acre" = size.i,
                            "Nearest_ParkSize_hAcre_Acre" = size_ha.i,
                            "Nearest_Type" = type.i,
                            "Nearest_Type_hAcre" = type_ha.i)
    return(dt)
  }
  cl          <- plan(multisession, workers = 8)
  out.dist.dt <- future_lapply(1:dim(sales.sf)[1], f_dist,
                               future.seed = TRUE) %>%
    rbindlist()                                    %>%
    .[,FIPS_cbsa := msa.i]
  future:::ClusterRegistry("stop")
  
  fwrite(out.dist.dt, 
         here("Data", "OSM", paste0("Trans_OSMDistance_ParkSize_CBSA", msa.i, ".csv")))
  return(paste0("int for ", msa.i, " complete"))
}
o_raw.dt <- map(msa.v, f_parks_OSM_dist)

files.v <- list.files(here("Data", "OSM"),
                      pattern = "OSMDistance_ParkSize_CBSA",
                      full.names = TRUE)
full.dt <- lapply(files.v, fread) %>% rbindlist()
fwrite(full.dt,
       here("Data", "OSM", "OSMDistance_ParkSize_CBSA_allMSA.csv"))
