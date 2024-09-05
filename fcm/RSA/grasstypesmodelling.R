# Install and load necessary packages
install.packages("geodata")
install.packages("Recocrop")
install.packages("tmap")
library(geodata)
library(Recocrop)
library(tmap)

# Define working directory
work_dir <- "D:/OneDrive - CGIAR/CIAT/PastureModelling_Turkana_Case/Coding/
Turkana_County_Pasture_Modelling"
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# Download climate data for Kenya
rain <- geodata::worldclim_country("Kenya", var="prec", path=work_dir)
tavg <- geodata::worldclim_country("Kenya", var="tavg", path=work_dir)
adm <- geodata::gadm("Kenya", level=1, path=work_dir)

# Extract data for Turkana County
turkana_adm <- adm[adm$NAME_1 == "Turkana", ]

# Mask the climate data to Turkana County
rain_turkana <- mask(rain, turkana_adm)
tavg_turkana <- mask(tavg, turkana_adm)

# Create an Ecocrop model for "Cenchrus ciliaris"
crop <- Recocrop::ecocropPars("Cenchrus ciliaris L.")
m <- Recocrop::ecocrop(crop)

# Use the model to make predictions for Turkana County
plant <- predict(m, prec=rain, tavg=tavg)
# Check the raw prediction values
print(summary(plant))

# Process predictions to suitability index
p <- classify(plant > 0, cbind(0,NA)) * 1:12
pm <- median(p, na.rm=TRUE)
hv <- pm + m$duration
hv <- ifel(hv > 12, hv - 12, hv)

# Normalize suitability values to 0-100 range
hv <- (pm / 12) * 100

# Plot the suitability map
plot(hv)
lines(turkana_adm)

# Create an Ecocrop model for "Chloris gayana"
crop <- Recocrop::ecocropPars("Chloris gayana")
m <- Recocrop::ecocrop(crop)

# Use the model to make predictions for Turkana County
plant <- predict(m, prec=rain, tavg=tavg)
# Check the raw prediction values
print(summary(plant))

# Process predictions to suitability index
p <- classify(plant > 0, cbind(0,NA)) * 1:12
pm <- median(p, na.rm=TRUE)
hv <- pm + m$duration
hv <- ifel(hv > 12, hv - 12, hv)

# Normalize suitability values to 0-100 range
hv <- (pm / 12) * 100

# Plot the suitability map
plot(hv)
lines(turkana_adm)

# Create an Ecocrop model for "Eragrostis superba Peyr."
crop <- Recocrop::ecocropPars("Eragrostis superba Peyr.")
m <- Recocrop::ecocrop(crop)

# Use the model to make predictions for Turkana County
plant <- predict(m, prec=rain, tavg=tavg)
# Check the raw prediction values
print(summary(plant))

# Process predictions to suitability index
p <- classify(plant > 0, cbind(0,NA)) * 1:12
pm <- median(p, na.rm=TRUE)
hv <- pm + m$duration
hv <- ifel(hv > 12, hv - 12, hv)

# Normalize suitability values to 0-100 range
hv <- (pm / 12) * 100

# Plot the suitability map
plot(hv)
lines(turkana_adm)

