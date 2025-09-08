
setwd("~/Dropbox/Home/Projects/Tornado")
load("~/Dropbox/Home/Projects/Tornado/UStornado2022.rda")

library( fields)
library( viridis)

x<- cbind( UStornado$lon, UStornado$lat)
UStornado$x<- x
ind<- UStornado$Fscale >=3 & UStornado$x[,1]>-102 & 
  UStornado$x[,2]>=29&
  UStornado$x[,1]< -83 & !is.na( UStornado$Fscale)
subsetTornado<- UStornado[ind,]
plot( subsetTornado$x)
US( add=TRUE)

temp<- discretize.image(subsetTornado$x, m=60, n=60 )
loc<- make.surface.grid( temp$grid)
y<- c( temp$hist)
ctab<- tim.colors(max(y))
# set zero counts to grey 
ctab<- c("grey95",ctab)
imagePlot( as.surface( loc,y), col=ctab )
US( add=TRUE, col="magenta")

save( loc, y, file="tornadoCountsSpace.rda")






