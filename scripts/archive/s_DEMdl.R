################################################################################
##### Downloading National Map Data: Direct Download/API calls
################################################################################

### Note: Packages of Interest
libraries <- c("tidyverse",
               "tictoc",
               "tibble",
               "data.table",
               "sf",
               "terra",
               "here",
               "lubridate",
               "rvest",
               "jsonlite",
               "httr",
               "RSelenium",
               "haven")
pacman::p_load(char = libraries)

#####
### Working with the NM API to download data
#####

# Step 1: Create bounding box (in degree coordinates) of counties in MSA
fips.dt <- tigris::fips_codes                        %>%
  as.data.table()                                    %>%
  .[state_code %in% c("17", "18", "26", "39", "55")] %>%
  .[,FIPS_cnty := str_c(state_code, county_code)]

cnty.sf <- tigris::counties(year = 2021)       %>%
  mutate(FIPS_cnty = str_c(STATEFP, COUNTYFP)) %>%
  st_transform(., 
               crs = 4326)
  

f_slopedl_gl <- function(cnty.i,
                         cnty.sf.i = cnty.sf){
  
  # Step 1: Create bounding box (in degree coordinates) of counties in MSA
  cnty.sf.i <- cnty.sf.i %>%
    filter(FIPS_cnty == cnty.i)
  
  cnty_bb.dt <- copy(cnty.sf.i) %>%
    st_bbox(.)
  cnty_bb.dt[1] <- floor(cnty_bb.dt[1])
  cnty_bb.dt[2] <- floor(cnty_bb.dt[2])
  cnty_bb.dt[3] <- ceiling(cnty_bb.dt[3])
  cnty_bb.dt[4] <- ceiling(cnty_bb.dt[4])
  # cnty_bb.dt[c(1:2)] <- cnty_bb.dt[c(1:2)] - 1
  # cnty_bb.dt[c(3:4)] <- cnty_bb.dt[c(3:4)] + 1
  
  # Step 2: Use bounding box in API call to DEM
  base.dom  <- "https://tnmaccess.nationalmap.gov/api/v1/products?datasets=National%20Elevation%20Dataset%20(NED)%201%20arc-second"
  bbox.dom  <- paste0(cnty_bb.dt, collapse = ",")
  api.dom   <- str_c(base.dom, "&bbox=", bbox.dom, "&max=400")
  
  api.test  <- 500
  iter      <- 1
  out.dt_t  <- NULL
  while(api.test != 200 & is.null(out.dt_t) & iter < 10){
    tryCatch(test.json <- GET(api.dom),
             error = function(e) iter <- 100)
    api.test <- test.json$status_code
    iter     <- iter + 1
    
    if (api.test == 200){
      out.json  <- fromJSON(rawToChar(test.json$content))
      out.dt_t  <- dim(out.json$items)[1]
      Sys.sleep(10)
    }else{
      Sys.sleep(30)
    }
  }
  if (!is.null(out.dt_t)){
    out.json  <- fromJSON(rawToChar(test.json$content))
    out.dt    <- out.json$items                                         %>%
      as.data.table()                                                   %>% 
      .[,Date   := (str_sub(publicationDate, end = 10))]                %>%
      .[,AreaID := str_extract(title, "n[:digit:]{2}w[:digit:]{3}")]    %>%
      .[str_detect(format, "GeoTIFF")]                                  %>%
      .[,MostRecent := max(lubridate::ymd(Date), na.rm = TRUE), by = AreaID]   %>%
      .[Date == MostRecent]                                                    %>%
      .[,SaveName := str_c("USGS_", AreaID, "_", str_remove_all(Date, "[:punct:]"), ".tif")] 
    
    # Step 3: Download URLs
    f_dl_dem <- function(x){
      if (!file.exists(here("Data", "DEM", "Raw", out.dt[x]$SaveName))){
        tryCatch(download.file(out.dt[x]$downloadURL,
                               here("Data", "DEM", "Raw", out.dt[x]$SaveName),
                               method = "curl"),
                 error = function(e) print(paste0("skip ", x)))
        return("yes")
      }else{
        return("Already Exists")
      }
    }  
    map(1:dim(out.dt)[1], f_dl_dem)  
    
    out_filetags.dt <- copy(out.dt)                     %>%
      .[,c("SaveName", "AreaID", "Date", "MostRecent")] %>%
      .[,FIPS_cnty := cnty.i]
    fwrite(out_filetags.dt,
           here("Data", "DEM", "DL-Tags", paste0("DL-Tags_", cnty.i, ".csv")))
    print(cnty.i)
    Sys.sleep(5)
    return(out_filetags.dt)
  }else{
    out_filetags.dt <- data.table("SaveName" = "",
                                  "AreaID"   = "",
                                  "Date"     = "",
                                  "MostRecent" = "",
                                  "FIPS_cnty"  = cnty.i)
    fwrite(out_filetags.dt,
           here("Data", "DEM", "DL-Tags", paste0("DL-Tags_", cnty.i, ".csv")))
    Sys.sleep(5)
    return(out_filetags.dt)
  }
}
dem_tags.dt <- map(fips.dt$FIPS_cnty, f_slopedl_gl)

### Find failed downloads in list; rerun
fips_rr.v <- lapply(list.files(here("Data", "DEM", "DL-Tags"),
                               full.names = TRUE),
                    fread, colClasses = 'character') %>% rbindlist() %>%
  .[SaveName != ""]

dem_tags_m.dt <- map(fips.dt$FIPS_cnty[!(fips.dt$FIPS_cnty %in% fips_rr.v$FIPS_cnty)], f_slopedl_gl)


### After rerun, bind together tags and save; be sure to check for missing
dem_tags_full.dt <- lapply(list.files(here("Data", "DEM", "DL-Tags"),
                                      full.names = TRUE),
                           fread, colClasses = 'character') %>% rbindlist() %>%
  .[SaveName != ""]
fwrite(dem_tags_full.dt,
       here("Data", "DEM", "Raw", "DEM_tags.csv"))

### Clean up downloads to remove missing tiles
# dem_tags.dt <- fread(here("Data", "DEM", "Raw", "DEM_tags.csv"))
files.v    <- list.files(here("Data", "DEM", "Raw"))
filesize.v <- file.size(here("Data", "DEM", "Raw", files.v))
dem_dl.dt  <- data.table("SaveName" = files.v,
                         "FileSize" = filesize.v) %>%
  merge(dem_tags_full.dt,
        by = "SaveName")                          %>%
  .[FileSize > 10000]                             %>%
  .[Date == MostRecent]                           %>%
  .[,s_SaveName := str_c("s_", SaveName)]         %>%
  .[,a_SaveName := str_c("a_", SaveName)]

f_slopecreate <- function(file.i,
                          save.i,
                          dem_dl.dt.i = dem_dl.dt){
  temp.r  <- terra::rast(here("Data", "DEM", "Raw", file.i))
  if (!file.exists(here("Data", "DEM", "Raw", save.i))){
    temp.r  <- terra::rast(here("Data", "DEM", "Raw", file.i))
    slope.r <- terra::terrain(temp.r,
                              v = "slope")
    terra::writeRaster(slope.r,
                       here("Data", "DEM", "Raw", save.i))
    return(paste0(save.i, " created and saved"))
  }else{
    return(paste0(save.i, " already exists"))
  }
}
inp1 <- unique(dem_dl.dt$SaveName)
inp2 <- unique(dem_dl.dt$s_SaveName)
map2(inp1, inp2, f_slopecreate)

f_aspcreate <- function(file.i,
                        save.i,
                        dem_dl.dt.i = dem_dl.dt){
  temp.r  <- terra::rast(here("Data", "DEM", "Raw", file.i))
  if (!file.exists(here("Data", "DEM", "Raw", save.i))){
    temp.r  <- terra::rast(here("Data", "DEM", "Raw", file.i))
    slope.r <- terra::terrain(temp.r,
                              v = "aspect")
    terra::writeRaster(slope.r,
                       here("Data", "DEM", "Raw", save.i))
    return(paste0(save.i, " created and saved"))
  }else{
    return(paste0(save.i, " already exists"))
  }
}
inp1 <- unique(dem_dl.dt$SaveName)
inp2 <- unique(dem_dl.dt$a_SaveName)
map2(inp1, inp2, f_aspcreate)



# f_slopemerge <- function(msa.i,
#                          dem_dl.dt.i = dem_dl.dt){
#   dl.dt    <- copy(dem_dl.dt.i)[FIPS_cbsa == msa.i]
#   rc.l     <- vector(mode   = "list",
#                      length = dim(dl.dt)[1])
#   for (i in 1:dim(dl.dt)[1]){
#     temp.r <- terra::rast(here("Data", "DEM", "Raw", dl.dt$s_SaveName[i]))
#     rc.l[[i]] <- temp.r
#     rm(temp.r)
#   }
#   
#   out.r <- rc.l[[1]]
#   for (i in 2:dim(dl.dt)[1]){
#     out.r <- merge(out.r, rc.l[[i]])
#   }
#   
#   terra::writeRaster(out.r,
#                      here("Data", 'DEM', paste0("Slope_1ArcDegree_CBSA", msa.i, ".tif")))
# }
# map(unique(dem_dl.dt$FIPS_cbsa), f_slopemerge)
