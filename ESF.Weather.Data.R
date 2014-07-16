syr<-read.table("SyracuseNYNOAAWeatherHistory.csv",sep=",",header=TRUE)
library(dplyr)

str(syr)
View(syr)
syr$DATE<-as.character(syr$DATE)
syr$DATE<-as.Date(syr$DATE,"%Y%m%d")


#looks like SUNY ESF is the closest station with all the needed data
syr.esf<-filter(syr,STATION=="GHCND:USC00308386")
View(syr.esf)
str(syr.esf)

syr.esf$Month<-c(11,12,1:12,1:10)
syr.esf$Year<-c(11,11,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13)
View(syr.esf)
str(syr.esf)
# Put temperatures in degrees celcius
syr.esf$MMXT<-syr.esf$MMXT/10
syr.esf$MMNT<-syr.esf$MMNT/10
syr.esf$MNTM<-syr.esf$MNTM/10

#Put precipitation in mm instead of 10th of mm
syr.esf$TPCP<-syr.esf$TPCP/10

with(syr.esf, plot(DATE, MMXT))
with(syr.esf, plot(DATE, MMNT))
with(syr.esf, plot(DATE, MNTM))
syr.esf$Month<-as.character(syr.esf$Month)

layout(rbind(1,2), heights=c(15,1))
with(syr.esf, plot(Month,MNTM,type="n", xlim=c(3,10),xlab="Month",ylab="Degrees C",main="Mean Monthly Temperature"))
with(syr.esf[5:12,], lines(Month,MNTM,col="red")) 
with(syr.esf[5:12,], points(Month,MNTM,col="red",pch=19))
with(syr.esf[17:24,], lines(Month, MNTM, col="blue"))
with(syr.esf[17:24,], points(Month, MNTM, col="blue", pch=19))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c(2012,2013),pch=19,
       col=c("red","blue"),ncol=2,bty ="n")

dev.off()
str(syr.esf)
?layout
lay<-layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(15,1,15,1))
layout.show(lay)

layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(15,1))
with(syr.esf, plot(Month,TPCP,type="n", xlim=c(3,10),xlab="Month",ylab="mm",main="Total Monthly Precipitation"))
with(syr.esf[5:12,], lines(Month,TPCP,col="red")) 
with(syr.esf[5:12,], points(Month,TPCP,col="red",pch=19))
with(syr.esf[17:24,], lines(Month, TPCP, col="blue"))
with(syr.esf[17:24,], points(Month, TPCP, col="blue", pch=19))
with(syr.esf, plot(Month,MNTM,type="n", xlim=c(3,10),xlab="Month",ylab="Degrees C",main="Mean Monthly Temperature"))
with(syr.esf[5:12,], lines(Month,MNTM,col="red")) 
with(syr.esf[5:12,], points(Month,MNTM,col="red",pch=19))
with(syr.esf[17:24,], lines(Month, MNTM, col="blue"))
with(syr.esf[17:24,], points(Month, MNTM, col="blue", pch=19))
par(mar=c(0,0,0,0))
plot.new()
legend('center','groups',c(2012,2013),pch=19,
       col=c("red","blue"),ncol=2,bty ="n")



dev.off()


