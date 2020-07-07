# utilities.R
# Author: Margaret Swift
# Contact: margaret.swift@duke.edu
# Created: 11/17/19
# Last Updated: 4/6/20

# Loads libraries and package files for GJAM Time Series.
# Holds all universal Utility functions, sources, libraries, 
# and constants for KNP Project

################################################################################
# Libraries
################################################################################
message("---------------------------------------------------------------------")

message('loading packages...')
if (!require("pacman")) install.packages("pacman")
packages <- c('ggplot2','gjam','gridExtra','downloader','MASS','RcppArmadillo','Rcpp','sf','tidyverse')
suppressPackageStartupMessages({ pacman::p_load(char=packages) })
cat('   ')
message(paste(packages, collapse=", "))
message('   packages loaded.')

message('loading local clark files...')
clark <- '../clarkFiles/'
files <- list.files(clark)
cat('   ')
message(paste(files, collapse=", "))
for (f in files) {
  path <- file.path(clark, f)
  if (grepl('.R', toupper(f))) source(path)
  else Rcpp::sourceCpp(path)
}
create('../output')
message('   files loaded.')

message('sourcing from github...')
d <- "https://github.com/jimclarkatduke/gjam/blob/master/gjamTimeFunctions.R?raw=True"
source_url(d, prompt=F)
message('   github files downloaded.')

message('loading demo data...')
f <- '../data/demoData.RData'
load(f)
message(paste0('   ', f))
message('   data loaded.')

dir.create('../output', showWarnings = FALSE)

################################################################################
# Functions
################################################################################

message('loading functions...')
message('   filePaths(), knpPlot(), plotSpecies()')
filePaths <- function(base, files) {
  for (i in 1:length(files)) files[i] <- file.path(base, files[i])
  files
}
knpPlot <- function(fill="transparent", color="black", grid=NULL, 
              shp_name='../data/borderFiles/knpBorder.shp',
              p4s="+proj=longlat +zone=36 +ellps=WGS84 +datum=WGS84 +no_defs",
              popdata=NULL, years=1989, sps=1) {
  
  # Read in shapefile and plot KNP border
  g <- st_read(shp_name) %>% st_transform(p4s)
  p <- ggplot() + geom_sf(data=g$geometry, color=color, fill=fill) + 
    scale_x_continuous(breaks=c(31, 31.5, 32)) +
    ggtitle('Kruger National Park')
  
  # If there's a grid file, add it to the map.
  if (!is.null(grid)) {
    p <- p + geom_sf(data=grid$., fill='transparent')
    
    # If there're species count data included, add those to the map too.
    # First we clean up the input a bit to make things more robust.
    if (!is.null(popdata) && length(popdata[,sps])) {
      # handle bad input
      if (!any(grepl("_", rownames(popdata)))) { 
        message('incorrect species data format. please use original ydata.') 
      } else {
        plots <- list()
        count <- 0
        for (i in 1:length(sps)) {
          for (j in 1:length(years)) {
            count <- count + 1
            plots[[count]] <- plotSpecies(sps[i], years[j], popdata, grid, p)
          }
        }
        do.call(grid.arrange, plots)
        return('plotting done')
      }
    }
  }
  p
}

plotSpecies <- function(sp, year, popdata, grid, p) {
  rows <- rownames(popdata)
  if (class(sp)=='numeric') sp <- colnames(popdata)[sp] #get sp name if inx
  inx.p <- which(grepl(paste0("_", year, "$"), rows))
  
  if (!length(inx.p)) { message(year, " is not a valid year!")
  } else {
    pops <- 0
    # Match population values to appropriate grid cells
    ids <- as.numeric(gsub("grid_|_[0-9]{4}$", '', rows))
    for (i in inx.p) {
      pop <- popdata[i,sp]
      g.inx <- which(grid$grid_id == ids[i])
      pops[g.inx] <- pop
    }
    # Plot populations in a grid
    q <- p + geom_sf(data=grid$., aes(fill=pops), color='transparent') + 
      ggtitle(paste(sp, 'pop. in', year)) + 
      scale_fill_gradient(low = "#CDE5FC",high = "#89011D") +
      labs(fill = "population") 
  }
  return(q)
}


message('   functions loaded.')
message("---------------------------------------------------------------------")

# EOF
