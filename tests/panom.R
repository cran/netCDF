library(netCDF)

precip<-open.netCDF("panom.nc")
precip
summary(precip)

panom<-read.netCDF(precip)
is.na(panom$PANOM)<- panom$PANOM==attr(panom$PANOM,"missing_value")

image(panom[[1]],panom[[2]],t(panom$PANOM),xlab="longitude",
      ylab="latitude",main="precipitation anomaly") 

close(precip)
q("no")

