################################################################################
### Building distributions of fields by ag district (and county) for soil chars
################################################################################

### Note: Packages of Interest
libraries <- c("tidyverse",
               "tictoc",
               "tibble",
               "data.table",
               "terra",
               "sf",
               "here",
               "future", 
               "future.apply",
               "raster",
               "furrr", 
               "pryr",
               "reticulate",
               "RSelenium",
               "rvest",
               "tigris")
pacman::p_load(char = libraries)

# set future max
options(future.globals.maxSize= 891289600)

#####
### Step 1: Load all csb's intersected with MUKEY's; bind into dataset
#####


f_load_csb <- function(FIPS_cnty.i){
  csb.sf <- st_read(here("Data", "USDA", "CSB", paste0("CSB_", FIPS_cnty.i, "_2022_centroid_DEM-SSURGO.gpkg")))
  st_geometry(csb.sf) <- NULL
  return(csb.sf)
}
fips.dt <- fips_codes %>%
  as.data.table()     %>%
  .[state %in% c("IL", "IN", "MI", "OH", "WI")] %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]
csb_full.dt <- map(fips.dt$FIPS_cnty, f_load_csb) %>% rbindlist() %>%
  .[,MUKEY := as.numeric(MUKEY)]


#####
### Step 2: Using gSSURGO, assemble the variables needed
#####
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
    .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "nccpi3all"))) := (lapply(.SD, function(x) comppct_r*x)),
      .SDcols = c("sandtotal_r", "silttotal_r", "claytotal_r", "nccpi3all")] %>%
    .[,(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "nccpi3all"))) := (lapply(.SD, function(x) sum(x, na.rm = TRUE))), 
      by = c("mukey"),
      .SDcols = str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "nccpi3all"))] %>%
    .[,c(str_c("w_", c("sandtotal_r", "silttotal_r", "claytotal_r", "nccpi3all")), "mukey")] %>% 
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
### Step 3: Merge ssurgo to CSB data; calculate necessary variables for integrated model
#####
csb_full.dt <- merge(csb_full.dt,
                     ssurgo_mu.dt,
                     by = "MUKEY",
                     all.x = TRUE) %>%
  .[,MedFieldSize_County := median(CSBACRES, na.rm = TRUE), by = "FIPS_cnty"]   %>%
  .[,MedFieldSize_AgDist := median(CSBACRES, na.rm = TRUE), by = "FIPS_st.ASD"]
fwrite(csb_full.dt,
       here("Data", "USDA", "CSB-SSURGO_FieldsSample_15-22.csv"))

#####
### Step 4: HUC Intersection
#####

v.huc2 = c("04", "05", "07")
f.huc <- function(i.huc){
  sf.huc <- st_read(here("Data", "USGS", "WBD", paste0("WBD_", i.huc, "_HU2_GPKG.gpkg")),
                    layer = "WBDHU8") %>%
    st_transform(3857)
  return(sf.huc)
}
l.sf.huc <- map(v.huc2, f.huc) %>% bind_rows(.) %>%
  .[,c("huc8", "name", "shape")]
st_write(l.sf.huc,
         here("Data", "USGS", "WBD", "USGS_HUC8.gpkg"),
         layer = "HUC8",
         delete_dsn = TRUE)

#####
### Create geometry from csv and intersect
#####

### Load csb; clean up from previous bad merges
dt.csb = fread(here("Data", "USDA", "CSB-SSURGO_FieldsSample_15-22.csv")) %>%
  .[, FIPS_cnty := str_c(FIPS_st, str_pad(FIPS_cnty, pad = "0", width = 3, side = "left"))] %>%
  setnames(c("w_sandtotal_r.x", "w_silttotal_r.x", "w_claytotal_r.x", "w_nccpi3all.x", 
             "MajorityClay.x", "MajoritySilt.x", "MajoritySand.x"),
           c("w_sandtotal_r", "w_silttotal_r", "w_claytotal_r", "w_nccpi3all",
             "MajorityClay", "MajoritySilt", "MajoritySand")) %>%
  .[,c("w_sandtotal_r.y", "w_silttotal_r.y", "w_claytotal_r.y", "w_nccpi3all.y", 
       "MajorityClay.y", "MajoritySilt.y", "MajoritySand.y") := NULL] %>%
  .[,CSBID := as.character(CSBID)]

f.huc_intersection <- function(i.fips,
                               i.dt.csb = dt.csb,
                               i.sf.huc = l.sf.huc){
  
  ### Step 1: Filter CSB and create geometry; keep only csb id and geom
  i.dt.csb = i.dt.csb[FIPS_cnty %in% i.fips] %>%
    st_as_sf(.,
             coords = c("lon", "lat"),
             crs = 4326) %>%
    st_transform(3857) %>%
    .[,c("CSBID", "geometry")]
  
  ### Step 2: Intersect w/ huc
  i.dt.csb = st_intersection(i.dt.csb, i.sf.huc)
  st_geometry(i.dt.csb) = NULL
  i.dt.csb = as.data.table(i.dt.csb) %>%
    setnames(c("name", "huc8"),
             c("huc8_name", "huc8_code"))
  print(paste0("Run complete for ", i.fips))
  return(i.dt.csb)
}
plan(multisession, workers = 3)
l.dt.int    <- future_map(unique(dt.csb$FIPS_cnty),
                          f.huc_intersection) %>% rbindlist()
future:::ClusterRegistry("stop")

### Merge intersections to dt.csb and write file
fwrite(l.dt.int,
       here("Data", "USDA", "CSB-HUC8_15-22.csv"))