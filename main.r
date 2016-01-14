# MGI 22, Joris Wever, Erwin van den Berg
# January 2016
# Exercise 8 of the Geo-Scripting Course from WUR
# Objective: create linear model object and predict vcf values for the Gewata area

# libraries
library(raster)

source('R/importFromUrl.R')

importFromUrl("https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/archive/gh-pages.zip")

#load the raster layers from the data folder
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB1.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB2.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB3.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB4.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB5.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/GewataB7.rda")
load("data/AdvancedRasterAnalysis-gh-pages/data/vcfGewata.rda")

# remove unrealistic values from the vcf (tree cover) file
vcfGewata[vcfGewata > 100] <- NA

# put all bands in a raster brick
alldata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(alldata) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")

# plot the bands in pairs to show correlation
pairs(alldata)

# highly correlated bands
pairs(brick(vcfGewata, GewataB3))
pairs(brick(vcfGewata, GewataB7))

# extract all data to a data.frame
df <- as.data.frame(getValues(alldata))

##### step 1: construct a linear model
modelLM <- lm(formula = VCF ~ band1 + band2 + band3 + band4 + band5 + band7, data = df) 
summary(modelLM)
# In the summary a column named ‘Pr(>|t|)’ is given. Small values mean that the probability that 
# this correlation results from coincidence is really small, so the bands with the small p-values 
# are important. So they are all important, except band 7. 

##### step 2: predict tree cover
# construct new brick layer without vcf_Gewata
gewata_novcf <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7)
names(gewata_novcf) <- c("band1", "band2", "band3", "band4", "band5", "band7")

# predict vcf values for the Gewata area based on Landsat bands
predLC <- predict(gewata_novcf, model=modelLM, na.rm=TRUE)
predLC[predLC < 0] <- NA

# plot predicted tree cover
plot(predLC, main="Predicted tree coverage", zlim=c(0,100))
# plot actual tree cover
plot(vcfGewata, main="Actual tree coverage", zlim=c(0,100))

# Calculate RMSE(Root Mean Squarred Error)
sqdif <- (vcfGewata-predLC)^2
sqdifdf <- as.data.frame(sqdif)
RMSE <- sqrt(mean(sqdifdf$layer, na.rm=TRUE))
RMSE

##### step 3: calculate difference between predicted and actual tree coverage per class

# prepare training polygons
load("data/AdvancedRasterAnalysis-gh-pages/data/trainingPoly.rda")
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
trainingPoly@data$Code
ndvi <- overlay(GewataB4, GewataB3, fun=function(x,y){(x-y)/(x+y)})
classes <- rasterize(trainingPoly, ndvi, field='Code')

# calculate RMSE for each class
averagediff<-zonal(sqdif, classes)
RMSEClasses<-sqrt(averagediff[,2])
names(RMSEClasses) <- c("cropland", "forest", "wetland")

# RMSE for each class
RMSEClasses

# The class wetland shows the biggest difference between predicted and actual tree cover

