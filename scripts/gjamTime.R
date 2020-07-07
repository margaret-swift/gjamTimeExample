# gjamTime.R
# Author: Margaret Swift
# Contact: margaret.swift@duke.edu
# Created: 11/17/19
# Last Updated: 4/6/20

# This file runs as a tutorial for GJAM Time Series using data collected at
#    Kruger National Park in South Africa. Response data are herbivore count
#    data from 1989 to 1996; covariates include precipitation, minimum temps,
#    soil and geology types, distance from water sources, and grass biomass.
# For more about GJAM, open the vignette: vignette('gjamVignette')
# For more about the Time Series version, see 'information/GJAMTime.html'.

################################################################################
# SETUP
################################################################################

# Clean and set working directory; source utilities file (loads functions, data,
# libraries, and GJAM files.)
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()[['path']]))
source('utilities.R')

################################################################################
# Read in data & plot grids
################################################################################

# saves ydata from further manipulations
species.data <- ydata 
# plots grid overlaid on KNP map
knpPlot(grid=knp.grid) 
# you can also plot species populations by year.
knpPlot(grid=knp.grid, popdata=species.data, sp='sp1', year=1991)

################################################################################
# Fill Missing Times
################################################################################

# Fill in missing data for time series & add initialized time-0 row.
tmp   <- gjamFillMissingTimes(xdata,ydata,edata, missingEffort=0.1, FILLMEANS=T,
                              groupCol='grid', timeCol='year',
                              groupVars=c('grid'), typeNames='DA')

# Reset variables with the missing time fills.
xdata <- tmp$xdata
edata <- tmp$edata
ydata <- tmp$ydata
effort<- list(columns=1:ncol(ydata), values=edata)
tlist <- tmp$timeList
S <- ncol(ydata)

################################################################################
# SET PRIORS
################################################################################

#------------------------------------------------------------------------------
# Beta Priors: Density-Independent Environment Effects ( B * x * x' )
#------------------------------------------------------------------------------

# We are trying to estimate the amount of em-/immigration into/out of each grid
# that is affected by each covariate. 'Lo' of 0 means we think it will have a 
# positive effect. We leave 'hi' to be figured out by the data, except for 'dist
# from rivers', which we expect to have a negative effect on immigration; and 
# 'anom.grass', which we expect to have a pretty strong effect (5% increase)
# on immigration.

# FORMULAE
formulaB <- as.formula(~ clay + anom.grass + dist.rivers + anom.ppt + veld.sav)

# PRIORS
priorB <- list(lo = list(dist.rivers=-0.3, anom.grass=-0.3), 
               hi = list(dist.rivers=0.3, anom.grass=0.3) )

#------------------------------------------------------------------------------
# Rho Priors: Density-Dependent Environment-Species Effects ( R * v * v' )
#------------------------------------------------------------------------------

# We are trying to estimate the amount of population growth that is affected by 
# each covariate. We only include covariates in the rho.form if we think they 
# have an effect. 'Lo' of 0 means we think it will have a  positive effect. 'hi' 
# of 0 is a negative effect. For this example, I think that 'distance from 
# rivers' and 'dry-season rainfall anomaly from the mean' will have a positive 
# effect on species growth.
formulaR <- as.formula(~ dist.rivers + anom.ppt)
priorR  <- list(lo = list(dist.rivers=0, anom.ppt=0),
                hi = list(dist.rivers=2, anom.ppt=2) )

#------------------------------------------------------------------------------
# Alpha Priors: Density-Independent Species Interactions ( A * w * w' )
#------------------------------------------------------------------------------

# Try different combinations of interaction signs. For this demo, let's say I 
# know that sp3 and sp4 often cooperate, and I make their interaction positive.
signA <- matrix(-1, S, S)
signA[3,4] <- signA[4,3] <- 1
colnames(signA) <- rownames(signA) <- colnames(ydata)

################################################################################
# MODEL PARAMETERS
################################################################################

priorList <- list( formulaBeta = formulaB, formulaRho = formulaR, 
                   betaPrior = priorB, rhoPrior = priorR, alphaSign = signA)
priors    <- gjamTimePrior( xdata, ydata, edata, priorList )
timeList  <- mergeList( tlist, priors )

# Try setting the ng to 2000 and burnin to 500 for a longer run.
modelList <- list( typeNames = 'DA', ng = 500, burnin = 50,  
                  timeList = timeList, effort = effort ) 

################################################################################
# RUN THE MODEL
################################################################################

output <- .gjam(formulaB, xdata, ydata, modelList, verbose=T)
save(output, file=file.path('../data', 'outGJAMTime.RData'))

################################################################################
# PLOT OUTPUT
################################################################################

# Set "SAVEPLOTS=T" to save the plots to the data folder.
plotPars <- list(GRIDPLOTS=T, CLUSTERPLOTS=T, PLOTALLY=T, SAVEPLOTS=F, 
                 outFolder='../data/gjamTimePlots')
.gjamPlot(output, plotPars)

# you may note that the contribution of immigration/emigration is small; this is
#   probably due to most of the park being fenced at this point in time.
# another thing to notice is that anom.grass doesn't seem to affect most of the 
#   species in this demo. try running this again with thi variable removed from 
#   the analysis.


# EOF
