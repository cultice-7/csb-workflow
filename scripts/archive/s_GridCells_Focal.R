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
               "httr",
               "Matrix",
               "arrow")
pacman::p_load(char = libraries)

################################################################################
### Focal Areas Intersections
################################################################################

grid.loc <- "K:/CFAES/AEDE/Arcdata/BrianCultice/INFEWS_LandUse/"

### i to j focal distance matrix (maybe network distance later??)
f_distmat <- function(self,
                      grid.loc.i = grid.loc,
                      searchdist.i = 50000){
  
  ### Step 1: Load Focal Data
  fgrid.sf     <- st_read(paste0(grid.loc, "Data/Gridded/GL_Gridded_300m-3km_Focal.gpkg"))
  fgrid_cen.sf <- st_centroid(fgrid.sf)  
  fgrid_buf.sf <- st_buffer(fgrid_cen.sf,
                            dist = searchdist.i)
  
  idwalk.dt <- data.table("FocalID" = fgrid_cen.sf$FocalID) %>%
    setorder(FocalID)                                       %>%
    .[,DistMatrixID := .I]
  
  ### Step 2: Calculate distance matrix
  f_dist <- function(i,
                     fgrid_cen.sf.i = fgrid_cen.sf,
                     fgrid_buf.sf.i = fgrid_buf.sf,
                     idwalk.dt.i    = idwalk.dt){
    
    dist_clear.sf <- st_intersection(fgrid_buf.sf.i[i,], fgrid_cen.sf.i[,c("FocalID", "geom")]) %>%
      mutate(Dist_ij_km = round(as.numeric(st_distance(., fgrid_cen.sf.i[i,]))/1000, 0))        %>%
      as.data.table()                                                                           %>%
      .[,c("FocalID", "FocalID.1", "Dist_ij_km")]                                               %>%
      setnames(c("FocalID", "FocalID.1"),
               c("FocalID.i", "FocalID.j"))
    return(dist_clear.sf)
  }
  cl          <- plan(multisession, workers = 4)
  dist.l <- future_lapply(1:dim(fgrid_cen.sf)[1], 
                          f_dist,
                          future.seed = TRUE,
                          future.scheduling = FALSE)
  future:::ClusterRegistry("stop")
  
  ### Step 3: Bind tables and save as parquet file
  dist.l <- rbindlist(dist.l)
  write_parquet(dist.l,
                paste0(grid.loc.i, "/Data/Gridded/Focal/Distances/GL_Gridded_300m-3km_FocalDist.parquet"))
}

### Weather
f_focal_PRISM <- function(i){
  ### Step 1: Load Focal Areas
  fgrid.sf <- st_read(here("Data", "Gridded", "GL_Gridded_300m-3km_Focal.gpkg"))
  
  ### Step 2: Program and run function to load PRISM and intersect (note that PRISM is daily)
  f_focal_PRISM_yr <- function(yr.i){
    prism.loc <- "K:/CFAES/AEDE/Arcdata/Data_PRISM/Data/PRISM/Raw/Data/Gridded\Focal\Distances"
    ### Load Rasters from File
    rast.r <- terra::rast(paste0(prism.loc, "/", yr.i))
  }
}