# gjamTimeExample
example dataset and script for gjamTime.

The purpose of this file is to run a Generalized Joint Attribute Model on 
herbivore count data from Kruger National Park in South Africa, from 1989 to
1996. Herbivore counts were collected during Kruger's annual Ecological 
Aerial Survey, conducted via airplane by four observers over a series of 
lateral transect lines across the park.

Count data and covariates were aggregated into approximately 400 grid squares 
spanning the area of the park. Although effort was constant throughout the 
survey, edge grid cells were assigned 'effort' proportional to the area 
contained within the park.

Covariates (aggregated across each grid cell) include: 

## Percentage of:
  - clay content in soil
  - basalt content in underlying geologies
  
## Centroid distance from:
  - water sources (rivers)
  
## Mean over all years, and annual anomalies from the full mean, for:
  - grass biomass
  - rainfall
  - minimum temperature

The main run file is found in ./scripts/gjamTime.R

For more about GJAM, open the vignette: vignette('gjamVignette') or go to "https://cran.r-project.org/web/packages/gjam/vignettes/gjamVignette.html". For more about the Time Series version, see './information/GJAMTime.html'.
