library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
library("Hmisc", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

alldata<-read.table("All_Integrated_Output.csv",sep=",",header=TRUE)
alldata$date<-as.Date(alldata$date, "%Y-%m-%d")

###Adjustments### ----
# Remove October observations
alldata<-alldata[alldata$Month!=10,]
alldata<-alldata[alldata$Month!=9,]
alldata$logcyano<-log(1+alldata$cyano, 10)

#Setting up subsets to plot -----
rl.2012<-filter(alldata,Year==2012, lake=="RL")
rl.2013<-filter(alldata,Year==2013, lake=="RL")
gl.2012<-filter(alldata, Year==2012, lake=="GL")
gl.2013<-filter(alldata, Year==2013, lake=="GL")
all.2012<-filter(alldata,Year==2012)
all.2013<-filter(alldata,Year==2013)
rl.all<-filter(alldata,lake=="RL")
gl.all<-filter(alldata,lake=="GL")

alldata$bin1.cyanos<-log((alldata$bin1.cyanos+1), 10)
alldata$bin2.cyanos<-log((alldata$bin2.cyanos+1), 10)
alldata$bin3.cyanos<-log((alldata$bin3.cyanos+1), 10)
alldata$bin4.cyanos<-log((alldata$bin4.cyanos+1), 10)

#log of cyanophage abundance

alldata$log.rlcp1<-log((alldata$rlcp1+1),10)
alldata$log.rlcp2a<-log((alldata$rlcp2a+1),10)
alldata$log.rlcp4<-log((alldata$rlcp4+1),10)

dev.off()



layout(matrix(c(1,2,3,4,5,6,7,8,9,10),5,2,byrow=TRUE),heights=c(15,15,15,15,3))
par(mar=c(2,2,1,1))
with(all.2012, plot(date,bin2.cyanos,type="n",ylim=c(0,5),xlab="Date",ylab="log(cells/ml)",main="Cyanobacteria GL 2012"))
with(gl.2012,lines(date,bin1.cyanos, col="forestgreen"))
with(gl.2012, points(date,bin1.cyanos, pch=19, col="forestgreen"))
with(gl.2012,lines(date,bin2.cyanos, col="purple"))
with(gl.2012, points(date,bin2.cyanos, pch=19, col="purple"))
with(gl.2012,lines(date,bin3.cyanos, col="orange"))
with(gl.2012, points(date,bin3.cyanos, pch=19, col="orange"))
with(gl.2012,lines(date,bin4.cyanos, col="blue"))
with(gl.2012, points(date,bin4.cyanos, pch=19, col="blue"))
#gl 2012
with(all.2012, plot(date,log.rlcp1,type="n", xlab="Date",ylim=c(0,4),ylab="log(copies/ml)",main="Cyanophage GL 2012"))
with(gl.2012,lines(date,log.rlcp1, col="navy"))
with(gl.2012, points(date,log.rlcp1, pch=18, col="navy"))
with(gl.2012,lines(date,log.rlcp2a, col="magenta"))
with(gl.2012, points(date,log.rlcp2a, pch=18, col="magenta"))
with(gl.2012,lines(date,log.rlcp4, col="turquoise"))
with(gl.2012, points(date,log.rlcp4, pch=18, col="turquoise"))
#RL 2012
with(all.2012, plot(date,bin3.cyanos,type="n",xlab="Date",ylim=c(0,5),ylab="log(cells/ml)",main="Cyanobacteria RL 2012"))
with(rl.2012, lines(date,bin1.cyanos,col="forestgreen"))
with(rl.2012, points(date,bin1.cyanos,pch=19, col="forestgreen"))
with(rl.2012, lines(date,bin2.cyanos,col="purple"))
with(rl.2012, points(date,bin2.cyanos,pch=19, col="purple"))
with(rl.2012, lines(date,bin3.cyanos,col="orange"))
with(rl.2012, points(date,bin3.cyanos,pch=19, col="orange"))
with(rl.2012, lines(date,bin4.cyanos,col="blue"))
with(rl.2012, points(date,bin4.cyanos,pch=19, col="blue"))

with(all.2012, plot(date,log.rlcp1,type="n", xlab="Date",ylim=c(0,4),ylab="log(copies/ml)",main="Cyanophage RL 2012"))
with(rl.2012,lines(date,log.rlcp1, col="navy"))
with(rl.2012, points(date,log.rlcp1, pch=18, col="navy"))
with(rl.2012,lines(date,log.rlcp2a, col="magenta"))
with(rl.2012, points(date,log.rlcp2a, pch=18, col="magenta"))
with(rl.2012,lines(date,log.rlcp4, col="turquoise"))
with(rl.2012, points(date,log.rlcp4, pch=18, col="turquoise"))
#GL 2013
with(all.2013, plot(date,bin1.cyanos,type="n", ylim=c(0,5), xlab="Date",ylab="log(cells/ml)",main="Cyanobacteria GL 2013"))
with(gl.2013,lines(date,bin1.cyanos, col="forestgreen"))
with(gl.2013, points(date,bin1.cyanos, pch=19, col="forestgreen"))
with(gl.2013,lines(date,bin2.cyanos, col="purple"))
with(gl.2013, points(date,bin2.cyanos, pch=19, col="purple"))
with(gl.2013,lines(date,bin3.cyanos, col="orange"))
with(gl.2013, points(date,bin3.cyanos, pch=19, col="orange"))
with(gl.2013,lines(date,bin4.cyanos, col="blue"))
with(gl.2013, points(date,bin4.cyanos, pch=19, col="blue"))

with(all.2013, plot(date,log.rlcp2a,type="n", ylim=c(0,4),xlab="Date",ylab="log(copies/ml)",main="Cyanophage GL 2013"))
with(gl.2013,lines(date,log.rlcp1, col="navy"))
with(gl.2013, points(date,log.rlcp1, pch=18, col="navy"))
with(gl.2013,lines(date,log.rlcp2a, col="magenta"))
with(gl.2013, points(date,log.rlcp2a, pch=18, col="magenta"))
with(gl.2013,lines(date,log.rlcp4, col="turquoise"))
with(gl.2013, points(date,log.rlcp4, pch=18, col="turquoise"))

with(all.2013, plot(date,bin1.cyanos,type="n", ylim=c(0,5), xlab="Date",ylab="log(cells/ml)",main="Cyanobacteria RL 2013"))
with(rl.2013, lines(date,bin1.cyanos,col="forestgreen"))
with(rl.2013, points(date,bin1.cyanos,pch=19, col="forestgreen"))
with(rl.2013, lines(date,bin2.cyanos,col="purple"))
with(rl.2013, points(date,bin2.cyanos,pch=19, col="purple"))
with(rl.2013, lines(date,bin3.cyanos,col="orange"))
with(rl.2013, points(date,bin3.cyanos,pch=19, col="orange"))
with(rl.2013, lines(date,bin4.cyanos,col="blue"))
with(rl.2013, points(date,bin4.cyanos,pch=19, col="blue"))


with(all.2013, plot(date,log.rlcp2a,type="n", ylim=c(0,4), xlab="Date",ylab="log(copies/ml)",main="Cyanophage RL 2013"))
with(rl.2013,lines(date,log.rlcp1, col="navy"))
with(rl.2013, points(date,log.rlcp1, pch=18, col="navy"))
with(rl.2013,lines(date,log.rlcp2a, col="magenta"))
with(rl.2013, points(date,log.rlcp2a, pch=18, col="magenta"))
with(rl.2013,lines(date,log.rlcp4, col="turquoise"))
with(rl.2013, points(date,log.rlcp4, pch=18, col="turquoise"))

par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("Type1","Type2","Type3","Type4"),pch=19,
       col=c("forestgreen","purple","orange","blue"),ncol=4,bty ="n")
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("rlcp1","rlcp2a","rlcp4"),pch=18,
       col=c("navy","magenta","turquoise"),ncol=3,bty ="n")

dev.off()





