################################################################################
### Downloading and Processing PRISM Climate Data for Multiple Purposes
### Description
###   - Loads and saves raw rasters
###   - Generates summary stats for locales (e.g. avg temp, growing season, etc)
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
               "sf",
               "terra",
               "exactextractr")
pacman::p_load(char = libraries)

################################################################################
### Loading Raw Data
################################################################################

weather.v     <- c("ppt", "tmax", "tmin",
                   "tmean", "vpdmax", "vpdmin")
inp1          <- rep(1990:1999,
                     each = length(weather.v))
inp2          <- rep(weather.v,
                     times = length(1990:1999))

f_prism_daily <- function(yr.i,
                          data.i){
  if (!dir.exists(here("Data", "PRISM", "Raw", yr.i))){
    dir.create(here("Data", "PRISM", "Raw", yr.i))
  }
  
  url.v  <- "https://ftp.prism.oregonstate.edu/daily"
  url.dl <- paste0(url.v, "/", data.i, "/", yr.i)
  html.dt <- read_html(url.dl) %>%
    html_elements("table")     %>%
    html_table()               %>%
    .[[1]]                     %>%
    as.data.table()            %>%
    .[str_detect(Name, "PRISM")] %>%
    .[,Date_loc := str_locate(Name, yr.i)[,1]] %>%
    .[,Date := str_sub(Name, 
                       start = Date_loc,
                       end = Date_loc + 7)]    %>%
    .[,httploc  := paste0(..url.dl, "/", Name)] %>%
    .[,dlloc    := here("Data", "PRISM", "Raw", ..yr.i, Name)] %>%
    .[,unziploc := here("Data", "PRISM", "Raw", ..yr.i,
                        str_remove(Name, ".zip"))]
  
   
  ### Download rasters
  f_dl <- function(i){
    if (!file.exists(html.dt$dlloc[i])){
      download.file(url      = html.dt$httploc[i],
                    destfile = html.dt$dlloc[i])
      Sys.sleep(1)
      return(print(html.dt$Date[i]))
    }
  }
  map(1:dim(html.dt)[1], f_dl)
  
  ### Unzip file
  f_unzip <- function(i){
    unzip(zipfile = html.dt$dlloc[i],
          exdir   = here("Data", "PRISM", "Raw", yr.i))
    print(html.dt$Date[i])
  }
  map(1:dim(html.dt)[1], f_unzip)
}
map2(as.character(inp1), inp2, f_prism_daily)
