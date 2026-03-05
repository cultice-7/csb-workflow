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
               "rgdal",
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
               "terra",
               "exactextractr")
pacman::p_load(char = libraries)

dirloc.v   <- "K:\\CFAES\\AEDE\\Arcdata\\Data_PRISM\\Data\\NCCS"
datatype.v  <- c("huss", "pr", "sfcWind", "tas", "tasmax", "tasmin")
datatype.v2 <- c("hurs", "rlds", "rsds")
ssp.v      <- rep("historical", times = length(datatype.v2))
f_dl_nccs <- function(datatype.i,
                      ssp.i){
  
  dir.v <- paste0(dirloc.v, "\\", datatype.i)
  if(!dir.exists(dir.v)){
    dir.create(dir.v)
  }
  f_findfile <- function(x){
    slash.dt   <- str_locate_all(x,
                               "/")[[1]]
    n         <- dim(slash.dt)[1]
    slashloc  <- slash.dt[n,1]
    
    filename  <- str_sub(x,
                         start = slashloc + 1)
    return(filename)
  }
  
  ### Load all the various links from the catalog site; adding in trycatch for dl errors
  dl.link <- paste0("https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6/ACCESS-CM2/",
                    ssp.i, 
                    "/r1i1p1f1/",
                    datatype.i,
                    "/catalog.html")
  success.index <- 0
  iter          <- 1
  while (success.index == 0 & iter < 10){
    tryCatch(
      {
        nccs.html <- read_html(dl.link) %>%
          html_elements("table")        %>%
          html_elements("tr")           %>%
          html_elements("td")           %>%
          html_elements("a")            %>%
          html_attrs()                  %>%
          unlist()                      %>%
          str_remove(.,
                     pattern = "catalog\\.html\\?dataset=")
        
        success.index <- 1
      },
      error = function(e){
        iter = iter + 1
        Sys.sleep(15)
      }
    )
    if (length(nccs.html) == 0){
      iter = iter + 1
      Sys.sleep(15)
    }else{
      html.dt <- data.table("dl_link" = nccs.html) %>%
        .[,dl_link_full := paste0("https://ds.nccs.nasa.gov/thredds/fileServer/", dl_link)]%>%
        .[,Year         := str_extract(dl_link_full,
                                       "[:digit:]{4}")]                                   %>%
        .[,Data         := datatype.i]                                                     %>%
        .[,filename     := unlist(lapply(dl_link_full,
                                         f_findfile))]                                     %>%
        .[,saveloc      := str_c(dir.v, "\\", filename)]                                   %>%
        .[, version_i   := ifelse(str_detect(dl_link_full,
                                             "_v"),
                                  1,
                                  0)]                                                      %>%
        .[,max_ver      := max(version_i), by = "Year"]                                    %>%
        .[max_ver == version_i]
      
      success.index <- 1
    }
  }
  
  ### Download files
  f_dl_each <- function(yr.i){
    dt <- html.dt[Year == yr.i]
    
    success.index = 0
    iter          = 1
    if (!file.exists(dt$saveloc)){
      while (success.index == 0 & iter < 10){
        tryCatch(
          {
            download.file(url      = dt$dl_link_full,
                          destfile = dt$saveloc,
                          mode = "wb")
            success.index = 1
          },
          error = function(e){
            iter = iter + 1
            Sys.sleep(15)
          },
          warning = function(w){
            print(paste0("Issue w/ download ", dt$dl_link_full, ": incomplete dl"))
            success.index = 0
            Sys.sleep(15)
          }
        )
      }
    } else{
      success.index = 1
    }
    track.dt <- data.table("Year"    = yr.i,
                           "Success" = success.index)
    return(track.dt)
  }
  track_dl.dt <- map(2000:2014, f_dl_each) %>% rbindlist() %>%
    .[, data := datatype.i]
  fwrite(track_dl.dt,
         paste0(dirloc.v, "\\", paste0("dl_tracking_", ssp.i, "-", datatype.i, ".csv")))
  return(print(datatype.i))
}
map2(datatype.v2, ssp.v, f_dl_nccs)




### Post Check (Some files incomplete dl; check to rerun first time; using 1% error
datatype.v <- c("huss", "pr", "sfcWind", "tas", "tasmax", "tasmin")
f_dl_check <- function(datatype.i){
  files.dt <- data.table("filename" = list.files(path = paste0(dirloc.v, "\\", datatype.i),
                                                 full.names = TRUE)) %>%
    .[,filesize := unlist(lapply(filename, file.size))]              %>%
    .[,dif      := abs(filesize - max(filesize))/max(filesize)]      %>%
    .[,ind      := ifelse(dif > 0.01,
                          1,
                          0)]                                        %>%
    .[ind == 1]
  return(files.dt)
}
dl_check.dt <- map(datatype.v, f_dl_check) %>% rbindlist()
file.remove(dl_check.dt$filename)
