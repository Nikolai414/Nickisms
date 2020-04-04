library(zoo)
library(tseries)
library(forecast)
library(xts)

# since this data were cleaned (extensively!!!) to an annual trend, I have already "decomposed" the data and removed any potential seasonal components
cod_1902.2018 <- read.csv("C:/Users/nikol/Desktop/Class/2020Spring/Quantitative-Ecology/Historic_Data/Data Sets Used for modelling/cod_1902.2017.csv")
#Time Series Analysis
{
land.mt <- ts(cod_1902.2018$Land_mt, start =1, frequency = 1)
par(mfrow=c(1,1), mai=c(0.25,0.8,0.1,0.1))
plot( land.mt, typ="l", 
      ylab= "Weight (mt)", main = "Atlantic cod commerical landings")                             
lines(tsclean(land.mt), col = "blue")

max(cod_1902.2018$Land_mt)
#Assumption 1: test for stationarity
##p-value < 0.05 indicates the Time Series is stationary
adf.test(land.mt)

#Developing autocorrelation plots- lag max was determined by the length of A.cod lifespan (25 years). 
par(mfrow=c(2,1))
acf(land.mt, lag.max= 25, main="ACF landings")
pacf(land.mt, lag.max = 25, main = "PACF landings")

# Fitting the ARIMA model
arima.land.mt1 <- auto.arima(land.mt, trace =T)

tsdisplay(residuals(arima.land.mt1), lag.max=25, main = "ARIMA residuals for cod landings") # Works great! everything falls within the bounds meaning there are no significant correlations present. 

AIC(arima.land.mt1)
  #2144.029

par(mfrow=c(1,1))
plot(land.mt, typ="l"); lines(fitted(arima.land.mt1), col = "blue")


#Test for independence:
checkresiduals(arima.land.mt1, lag = 20)
#p-value is greater than 0.05 so this is independent

#forecast metric tonnage 10 years into the future. 
par(mfrow=c(2,1))
plot(forecast(arima.land.mt1, h = 10))
plot(forecast(arima.land.mt1, h = 25))
#plot(forecast(arima.land.mt1, h = 50))

# Creating a time series object
temp <- ts(cod_1902.2018$Temp.avg.C, start =1, frequency = 1)

par(mfrow=c(1,1), mai = c(0.3, 0.8, 0.1,0.1))
plot(temp, typ = "l", ylab="Temperature (C)", xlab = "")

# remove outliers
plot(temp, typ = "l", ylab="Temperature (C)", xlab = "")
lines(tsclean(temp),col = "orange")

temp <- tsclean(temp)

# Again, not going to "decompose" since this is on an annual basis and there are not periods visible here. 

#test for stationarity
adf.test(temp)
#p-value is less than 0.05

# look for significant lags
ccf(diff(temp), land.mt, na.action = na.pass, lag.max = 30, plot = T)
# no significant correlation in this graph!

#fit an arima
arima.land.mt2_T <- auto.arima(land.mt, xreg=c(diff(temp),0), trace = TRUE)

# Check index for model comparison:
AIC(arima.land.mt1, arima.land.mt2_T)
# The first model was better but the difference is fairly small. 

# consider extreme temperatures (especially since eggs are only viable between [0.5 -10.3])
temp.i <- temp
temp.i[temp.i > 10.3] <- 0
temp.i[temp.i >= 30] <- 1

arima.land.mt3_Ti <- auto.arima(land.mt, xreg = temp.i, trace = T)

AIC(arima.land.mt1, arima.land.mt3_Ti)
# barely made any differenece, it looks like even extreme temperatures are not playing as much of a role on eggs. 

checkresiduals(arima.land.mt3_Ti, lag = 20)
# still above 0.05 so it's independent. 

par(mfrow = c(1,1))
plot(land.mt, typ = "l"); lines(fitted(arima.land.mt3_Ti), col = "red")

# Repeat with Salinity

sal <- ts(cod_1902.2018$Sal.avg.PSS, start =1, frequency = 1)
par(mfrow=c(1,1), mai = c(0.3, 0.8, 0.1,0.1))
plot(sal, typ = "l", ylab="Salinity (PSS)", xlab = "")

# remove outliers
plot(sal, typ = "l", ylab="Salinity (PSS)", xlab = "")
lines(tsclean(sal),col = "orange")

sal <- tsclean(sal)

# Again, not going to "decompose" since this is on an annual basis and there are not periods visible here. 

#test for stationarity
adf.test(diff(sal))
#p-value is less than 0.05

# look for significant lags
ccf(diff(sal), land.mt, na.action = na.pass, lag.max = 30, plot = T)
# no significant correlation in this graph!

#fit an arima
arima.land.mt4_S <- auto.arima(land.mt, xreg=c(diff(sal),0), trace = TRUE)

# Check index for model comparison:
AIC(arima.land.mt1, arima.land.mt4_S)
# The first model was better but the difference is fairly small. 

# consider extreme temperatures (especially since eggs are only viable between [0.5 -10.3])
sal.i <- 
sal.i[sal.i < 10.3] <- 0
sal.i[sal.i >= 25] <- 1

arima.land.mt3_Ti <- auto.arima(land.mt, xreg = sal.i, trace = T)

AIC(arima.land.mt1, arima.land.mt3_Ti)
# barely made any differenece, it looks like even extreme temperatures are not playing as much of a role on eggs. 

checkresiduals(arima.land.mt4_S, lag = 20)
# still above 0.05 so it's independent. 

par(mfrow = c(1,1))
plot(land.mt, typ = "l"); lines(fitted(arima.land.mt3_Ti), col = "purple")


}

#Species Distribution Model
Cod_Dist <- read.csv("C:/Users/nikol/Desktop/Class/2020Spring/Quantitative-Ecology/Historic_Data/Data Sets Used for modelling/Cod_Dist.csv")

library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library("marmap")
library("sdmpredictors")
# Check the data to make sure it loaded correctly
summary(Cod_Dist)


max.lat <- ceiling(max(Cod_Dist$Lat))
min.lat <- floor(min(Cod_Dist$Lat))
max.lon <- ceiling(max(Cod_Dist$Long))
min.lon <- floor(min(Cod_Dist$Long))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))



Cod.dist.map <- getNOAA.bathy(lon1 = min.lon-1, lon2 = max.lon+1,
                        lat1 = min.lat-1, lat2 = max.lat+1, resolution = 1)

colorramp<- colorRampPalette(c("purple", "darkblue", "lightblue", "white"))

par(mfrow=c(1,1), mai = c(0.9, 0.8, 0.3,0.3))
plot(Cod.dist.map, 
     lwd = 0.1,
     bpal =  colorramp(50),
     image = TRUE)
scaleBathy(Cod.dist.map, deg = 2, x = "topleft", inset = 5)
# Add the points for individual observation
points(x = Cod_Dist$Long, 
       y = Cod_Dist$Lat, 
       col = "orange", 
       pch = 20, 
       cex = 0.75)

# And draw a little box around the graph
box()



###Building a model and visualizing results


datasets<- list_datasets(terrestrial= FALSE, marine = TRUE)
layers <- list_layers(datasets)
#Temperature Bio-ORACLE	BO2_tempmean_bdmean	Sea water temperature (mean at mean depth)	Mean sea water temperature at the bottom at mean bottom d
#Salinity Bio-ORACLE	BO2_salinitymean_bdmean	Sea water salinity (mean at mean depth)
{load_layers(layers[layers$name %in% c("Sea water temperature (mean at mean depth)", 
                                      "Sea water salinity (mean at mean depth)") & 
                      layers$dataset_code == "Bio-ORACLE",], datadir = ("C:/Users/nikol/Desktop/Class/2020Spring/Quantitative-Ecology/Historic_Data/update as of 1Apr/"))
}
layercodes <- c("BO2_tempmean_bdmean", "BO2_salinitymean_bdmean")

env<- load_layers(layercodes, equalarea= F, datadir= ("C:/Users/nikol/Desktop/Class/2020Spring/Quantitative-Ecology/Historic_Data/update as of 1Apr/"))


# Build species distribution model
Lat.Long <- data.frame(Cod_Dist$Long, Cod_Dist$Lat)
sdm.data <- crop(x = env, y = geographic.extent)

bc.model <- bioclim(x = sdm.data, p = Lat.Long)

# Predict presence from model
predict.presence <- dismo::predict(object = bc.model, x = env, ext = geographic.extent)


par(mfrow=c(1,1), mai = c(0.9, 0.8, 0.3,0.3))
plot(Cod.dist.map, 
     lwd = 0.1,
     bpal =  colorramp(50),
     image = TRUE)
scaleBathy(Cod.dist.map, deg = 2, x = "topleft", inset = 5)
# Add the points for individual observation
points(x = Cod_Dist$Long, 
       y = Cod_Dist$Lat, 
       col = "orange", 
       pch = 20, 
       cex = 0.75)
box()

# Add model probabilities
library("maptools")
data(wrld_simpl)

plot(wrld_simpl, xlim= c(min.lon-1, max.lon+1),
     ylim= c(min.lat-1,max.lat+1), axes = T, col = "Green")

colorramp.2 <- colorRampPalette(c("darkred", "red", "purple", "pink", "white"))
plot(predict.presence, add = TRUE, bpal = colorramp.2(5))
points(x = Cod_Dist$Long, 
       y = Cod_Dist$Lat, 
       col = "purple", 
       pch = 20, 
       cex = 0.75)
box()

#psuedo absence points were generated with another program.
Pseudo.Absence <- read.csv("C:/Users/nikol/Desktop/Class/2020Spring/Quantitative-Ecology/Historic_Data/Data Sets Used for modelling/Pseudo.Absence.csv")


par(mfrow=c(1,1), mai = c(0.9, 0.8, 0.3,0.3))
plot(Cod.dist.map, 
     lwd = 0.1,
     bpal =  colorramp(50),
     image = TRUE)
scaleBathy(Cod.dist.map, deg = 2, x = "topleft", inset = 5, main = "Presence and Pseudo-absence Points")
# Add the points for individual observation
box() 
  
points(x = Pseudo.Absence$Long,
       y = Pseudo.Absence$Lat,
       col = "grey", pch = 1, cex = 0.75) 

points(x = Cod_Dist$Long, 
       y = Cod_Dist$Lat, 
       col = "orange", 
       pch = 20, 
       cex = 0.75)


testing.group <- 1
group.presence <- kfold(x = Lat.Long, k = 5)
head(group.presence)
table(group.presence)


presence.train <- Lat.Long[group.presence != testing.group, ]
presence.test <- Lat.Long[group.presence == testing.group, ]

# Repeat the process for pseudo-absence points
group.background <- kfold(x = Pseudo.Absence, k = 5)
background.train <- Pseudo.Absence[group.background != testing.group, ]
background.test <- Pseudo.Absence[group.background == testing.group, ]


# Build a model using training data
bc.model <- bioclim(x = env, p = presence.train)

# Predict presence from model (same as previously, but with the update model)
predict.presence <- dismo::predict(object = bc.model, 
                                   x = env, 
                                   ext = geographic.extent)


# Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test,   # The presence testing data
                    a = background.test, # The absence testing data
                    model = bc.model,    # The model we are evaluating
                    x = env)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "prevalence")



#Where the 
plot(Cod.dist.map, 
     lwd = 0.1,
     bpal =  colorramp(50), main = "Georges Bank, Atlantic Ocean", 
     image = TRUE)
scaleBathy(Cod.dist.map, deg = 2, x = "topleft", inset = 7)

box()

plot(predict.presence > bc.threshold, add = TRUE, legend = F, col = c(NA, "darkred"))
points(x = Cod_Dist$Long, 
       y = Cod_Dist$Lat, 
       col = "green", 
       pch = 20, 
       cex =0.25)
plot(Cod.dist.map, add = TRUE, border = "chocolate")
box()
