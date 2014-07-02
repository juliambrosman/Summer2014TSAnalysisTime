###Load output table from "All Integrated" script:
library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

alldata<-read.table("All_Integrated_Output.csv",sep=",",header=TRUE)
alldata$date<-as.Date(alldata$date, "%Y-%m-%d")

###Adjustments### ----
# Remove October observations
alldata<-alldata[alldata$Month!=10,]

#log of cyano abundance:
alldata$bin1.cyanos<-log((alldata$bin1.cyanos+1), 10)
alldata$bin2.cyanos<-log((alldata$bin2.cyanos+1), 10)
alldata$bin3.cyanos<-log((alldata$bin3.cyanos+1), 10)
alldata$bin4.cyanos<-log((alldata$bin4.cyanos+1), 10)

#log of cyanophage abundance

alldata$log.rlcp1<-log((alldata$rlcp1+1),10)
alldata$log.rlcp2a<-log((alldata$rlcp2a+1),10)
alldata$log.rlcp4<-log((alldata$rlcp4+1),10)


#Setting up subsets to plot -----
rl.2012<-filter(alldata,Year==2012, lake=="RL")
rl.2013<-filter(alldata,Year==2013, lake=="RL")
gl.2012<-filter(alldata, Year==2012, lake=="GL")
gl.2013<-filter(alldata, Year==2013, lake=="GL")
all.2012<-filter(alldata,Year==2012)
all.2013<-filter(alldata,Year==2013)
palette(c("forestgreen","purple"))


##playing with the data##----
#boxplots: -----
boxplot(alldata$bac~alldata$Year+alldata$lake, main="bac")
boxplot(alldata$bac~alldata$Year+alldata$lake, main="bac")
boxplot(alldata$cyano~alldata$Year+alldata$lake, main="cyano")
boxplot(alldata$chla~alldata$Year+alldata$lake, main="chla")
boxplot(alldata$nitrate~alldata$Year+alldata$lake, main="nitrate")
boxplot(alldata$srp~alldata$Year+alldata$lake, main="srp")
boxplot(alldata$temp~alldata$Year+alldata$lake, main="temp")


##time series plots of counts/ml per year----
#Library with errbars function:
library("Hmisc", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

#vlp 2012
with(all.2012, errbar(date,vlp,
                      yplus=vlp+vlp.sd,yminus=vlp-vlp.sd,
                      col=lake,pch=19, xlab="Date",ylab="VLP/ml"))
title(main="Virus Like Particle Abundance 2012")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))

#vlp 2013
with(all.2013, errbar(date,vlp,
                      yplus=vlp+vlp.sd,yminus=vlp-vlp.sd,
                      col=lake,pch=19, xlab="Date",ylab="VLP/ml"))
title(main="Virus Like Particles 2013")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))
#bac 2012
with(all.2012, errbar(date,bac,
                      yplus=bac+bac.sd,yminus=bac-bac.sd,
                      col=lake,pch=19, xlab="date",ylab="bac/ml"))
title(main="Bacterial Abundance 2012")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))

#bac 2013
with(all.2013, errbar(date,bac,
                      yplus=bac+bac.sd,yminus=bac-bac.sd,
                      col=lake,pch=19, main="bac 2013", xlab="date",ylab="bac/ml"))
title(main="Bacterial Abundance 2013")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))

#cyano 2012
with(all.2012, errbar(date,cyano,
                      yplus=cyano+cyano.sd,yminus=cyano-cyano.sd,
                      col=lake,pch=19, main="cyano 2012", xlab="date",ylab="cyano/ml"))
title(main="Picocyanobacteria 2012")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))


#cyano 2013
with(all.2013, errbar(date,cyano,
                      yplus=cyano+cyano.sd,yminus=cyano-cyano.sd,
                      col=lake,pch=19, main="cyano 2013", xlab="date",ylab="cyano/ml"))
title(main="Picocyanobacteria 2013")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))

#plotting cyanophage genotype abundance over time----
library("Hmisc", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
#rlcp1 2013
with(all.2013, errbar(date,rlcp1,
                      yplus=rlcp1+rlcp1.sd,yminus=rlcp1-rlcp1.sd,
                      col=lake,pch=19, main="rlcp1 2013", xlab="date",ylab="copies/ml"))
title(main="Cyanophagetype rlcp1 2013")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))


#rlcp1 2012
with(all.2012, errbar(date,rlcp1,
                      yplus=rlcp1+rlcp1.sd,yminus=rlcp1-rlcp1.sd,
                      col=lake,pch=19, main="rlcp1 2013", xlab="date",ylab="copies/ml"))
title(main="Cyanophage type rlcp1 2012")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))


#rlcp2a 2013
with(all.2013, errbar(date,rlcp2a,
                      yplus=rlcp2a+rlcp2a.sd,yminus=rlcp2a-rlcp2a.sd,
                      col=lake,pch=19, main="rlcp2a 2013", xlab="date",ylab="copies/ml"))
title(main="Cyanophage type rlcp2a 2013")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))


#rlcp2a 2012
with(all.2012, errbar(date,rlcp2a,
                      yplus=rlcp2a+rlcp2a.sd,yminus=rlcp2a-rlcp2a.sd,
                      col=lake,pch=19, main="rlcp2a 2013", xlab="date",ylab="rlcp2a/ml"))
title(main="Cyanophage type rlcp2a 2012")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))


#rlcp4 2013
with(all.2013, errbar(date,rlcp4,
                      yplus=rlcp4+rlcp4.sd,yminus=rlcp4-rlcp4.sd,
                      col=lake,pch=19, main="rlcp4 2013", xlab="date",ylab="copies/ml"))
title(main="Cyanophage type rlcp4 2013")
legend("topright",c("GL","RL"),pch=19,col=c("forestgreen","purple"))


#rlcp4 2012
with(all.2012, errbar(date,rlcp4,
                      yplus=rlcp4+rlcp4.sd,yminus=rlcp4-rlcp4.sd,
                      col=lake,pch=19, main="rlcp4 2013", xlab="date",ylab="copies/ml"))
title(main="Cyanophage type rlcp4 2012")

##plotting metadata----

##Nitrate
#nitrate 2012
with(all.2012, plot(date,nitrate,
                      col=lake,pch=19, main="nitrate", xlab="date",ylab="uM Nitrate"))
#nitrate 2013
with(all.2013, plot(date,nitrate,
                    col=lake,pch=19, main="nitrate", xlab="date",ylab="uM Nitrate"))
##Nitrite
#nitrite 2012
with(all.2012, plot(date,nitrite,
                    col=lake,pch=19, main="nitrite", xlab="date",ylab="uM nitrite"))
#nitrite 2013
with(all.2013, plot(date,nitrite,
                    col=lake,pch=19, main="nitrite", xlab="date",ylab="uM nitrite"))
##SRP
#srp 2012
with(all.2012, plot(date,srp,
                    col=lake,pch=19, main="srp", xlab="date",ylab="uM srp"))
#srp 2013
with(all.2013, plot(date,srp,
                    col=lake,pch=19, main="srp", xlab="date",ylab="uM srp"))
##ChlA
#chla 2013
with(all.2013, errbar(date,chla,
                      yplus=chla+chla.sd,yminus=chla-chla.sd,
                      col=lake,pch=19, main="chla 2013", xlab="date",ylab="chla/ml"))
title(main="chla 2013")
#chla 2012
with(all.2012, errbar(date,chla,
                      yplus=chla+chla.sd,yminus=chla-chla.sd,
                      col=lake,pch=19, main="chla 2013", xlab="date",ylab="chla/ml"))
title(main="chla 2012")

##Plotting Cyanobacterial Genotype Abundance##----

#bin1.cyanos 2012
with(all.2012, plot(date,bin1.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 1 2012"))
with(gl.2012,lines(date,bin1.cyanos, col="forestgreen"))
with(rl.2012, lines(date,bin1.cyanos,col="purple"))
with(gl.2012, points(date,bin1.cyanos, pch=19, col="forestgreen"))
with(rl.2012, points(date,bin1.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))


#bin1.cyanos 2013
with(all.2013, plot(date,bin1.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 1 2013"))
with(gl.2013,lines(date,bin1.cyanos, col="forestgreen"))
with(rl.2013, lines(date,bin1.cyanos,col="purple"))
with(gl.2013, points(date,bin1.cyanos, pch=19, col="forestgreen"))
with(rl.2013, points(date,bin1.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))

#bin2.cyanos 2012
with(all.2012, plot(date,bin2.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 2 2012"))
with(gl.2012,lines(date,bin2.cyanos, col="forestgreen"))
with(rl.2012, lines(date,bin2.cyanos,col="purple"))
with(gl.2012, points(date,bin2.cyanos, pch=19, col="forestgreen"))
with(rl.2012, points(date,bin2.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))


#bin2.cyanos 2013
with(all.2013, plot(date,bin2.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 2 2013"))
with(gl.2013,lines(date,bin2.cyanos, col="forestgreen"))
with(rl.2013, lines(date,bin2.cyanos,col="purple"))
with(gl.2013, points(date,bin2.cyanos, pch=19, col="forestgreen"))
with(rl.2013, points(date,bin2.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))

#bin3.cyanos 2012
with(all.2012, plot(date,bin3.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 3 2012"))
with(gl.2012,lines(date,bin3.cyanos, col="forestgreen"))
with(rl.2012, lines(date,bin3.cyanos,col="purple"))
with(gl.2012, points(date,bin3.cyanos, pch=19, col="forestgreen"))
with(rl.2012, points(date,bin3.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))


#bin3.cyanos 2013
with(all.2013, plot(date,bin3.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 3 2013"))
with(gl.2013,lines(date,bin3.cyanos, col="forestgreen"))
with(rl.2013, lines(date,bin3.cyanos,col="purple"))
with(gl.2013, points(date,bin3.cyanos, pch=19, col="forestgreen"))
with(rl.2013, points(date,bin3.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))

#bin4.cyanos 2012
with(all.2012, plot(date,bin4.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 4 2012"))
with(gl.2012,lines(date,bin4.cyanos, col="forestgreen"))
with(rl.2012, lines(date,bin4.cyanos,col="purple"))
with(gl.2012, points(date,bin4.cyanos, pch=19, col="forestgreen"))
with(rl.2012, points(date,bin4.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))

#bin4.cyanos 2013
with(all.2013, plot(date,bin4.cyanos,type="n", xlab="Date",ylab="cells/ml",main="Cyano Type 4 2013"))
with(gl.2013,lines(date,bin4.cyanos, col="forestgreen"))
with(rl.2013, lines(date,bin4.cyanos,col="purple"))
with(gl.2013, points(date,bin4.cyanos, pch=19, col="forestgreen"))
with(rl.2013, points(date,bin4.cyanos,pch=19, col="purple"))
legend("topright",c("GL","RL"),pch=19, col=c("forestgreen","purple"))

##Plotting Cyanobacterial Genotype Abundance Per Lake Per Year----
?jpeg

#GL 2012

layout(rbind(1,2), heights=c(15,1))
with(all.2012, plot(date,bin3.cyanos,type="n", xlab="Date",ylab="log(cells/ml)",main="Cyanobacterial Genotype Abundance GL 2012"))
with(gl.2012,lines(date,bin1.cyanos, col="forestgreen"))
with(gl.2012, points(date,bin1.cyanos, pch=19, col="forestgreen"))
with(gl.2012,lines(date,bin2.cyanos, col="purple"))
with(gl.2012, points(date,bin2.cyanos, pch=19, col="purple"))
with(gl.2012,lines(date,bin3.cyanos, col="orange"))
with(gl.2012, points(date,bin3.cyanos, pch=19, col="orange"))
with(gl.2012,lines(date,bin4.cyanos, col="blue"))
with(gl.2012, points(date,bin4.cyanos, pch=19, col="blue"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("Type1","Type2","Type3","Type4"),pch=19,
       col=c("forestgreen","purple","orange","blue"),ncol=4,bty ="n")
dev.off()


#RL 2012
layout(rbind(1,2), heights=c(15,1))
with(all.2012, plot(date,bin3.cyanos,type="n", xlab="Date",ylab="log(cells/ml)",main="Cyanobacterial Genotype Abundance RL 2012"))
with(rl.2012, lines(date,bin1.cyanos,col="forestgreen"))
with(rl.2012, points(date,bin1.cyanos,pch=19, col="forestgreen"))
with(rl.2012, lines(date,bin2.cyanos,col="purple"))
with(rl.2012, points(date,bin2.cyanos,pch=19, col="purple"))
with(rl.2012, lines(date,bin3.cyanos,col="orange"))
with(rl.2012, points(date,bin3.cyanos,pch=19, col="orange"))
with(rl.2012, lines(date,bin4.cyanos,col="blue"))
with(rl.2012, points(date,bin4.cyanos,pch=19, col="blue"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("Type1","Type2","Type3","Type4"),pch=19,
       col=c("forestgreen","purple","orange","blue"),ncol=4,bty ="n")
dev.off()

#GL 2013
layout(rbind(1,2), heights=c(7,1))
with(all.2013, plot(date,bin1.cyanos,type="n", xlab="Date",ylab="log(cells/ml)",main="Cyanobacterial Genotype Abundance GL 2013"))
with(gl.2013,lines(date,bin1.cyanos, col="forestgreen"))
with(gl.2013, points(date,bin1.cyanos, pch=19, col="forestgreen"))
with(gl.2013,lines(date,bin2.cyanos, col="purple"))
with(gl.2013, points(date,bin2.cyanos, pch=19, col="purple"))
with(gl.2013,lines(date,bin3.cyanos, col="orange"))
with(gl.2013, points(date,bin3.cyanos, pch=19, col="orange"))
with(gl.2013,lines(date,bin4.cyanos, col="blue"))
with(gl.2013, points(date,bin4.cyanos, pch=19, col="blue"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("Type1","Type2","Type3","Type4"),pch=19,
       col=c("forestgreen","purple","orange","blue"),ncol=4,bty ="n")
dev.off()


#RL 2013
layout(rbind(1,2), heights=c(15,1))
with(all.2013, plot(date,bin1.cyanos,type="n", xlab="Date",ylab="log(cells/ml)",main="Cyanobacterial Abundance RL 2013"))
with(rl.2013, lines(date,bin1.cyanos,col="forestgreen"))
with(rl.2013, points(date,bin1.cyanos,pch=19, col="forestgreen"))
with(rl.2013, lines(date,bin2.cyanos,col="purple"))
with(rl.2013, points(date,bin2.cyanos,pch=19, col="purple"))
with(rl.2013, lines(date,bin3.cyanos,col="orange"))
with(rl.2013, points(date,bin3.cyanos,pch=19, col="orange"))
with(rl.2013, lines(date,bin4.cyanos,col="blue"))
with(rl.2013, points(date,bin4.cyanos,pch=19, col="blue"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("Type1","Type2","Type3","Type4"),pch=19,
       col=c("forestgreen","purple","orange","blue"),ncol=4,bty ="n")
dev.off()


##Plotting Bar Graph of all Cyano Genotypes per Lake:----

##pull out just cyano representation:
cyano.bins<-alldata[,c(1,2,21:24)]

##Make individual DF per bin and assign bin value in bin column
cyano.b1<-cyano.bins[,c(1,2,3)]
View(cyano.b1)
cyano.b1$bin<-1
colnames(cyano.b1)[3]<-"log.cyano.ml"

cyano.b2<-cyano.bins[,c(1,2,4)]
cyano.b2$bin<-2
colnames(cyano.b2)[3]<-"log.cyano.ml"

cyano.b3<-cyano.bins[,c(1,2,5)]
cyano.b3$bin<-3
colnames(cyano.b3)[3]<-"log.cyano.ml"

cyano.b4<-cyano.bins[,c(1,2,6)]
cyano.b4$bin<-4
colnames(cyano.b4)[3]<-"log.cyano.ml"

#bind rows for all the results into one data matrix
cyano.combined<-rbind(cyano.b1, cyano.b2, cyano.b3, cyano.b4)

library(ggplot2)
library(plyr)
library(reshape2)

#Extract 2012 Values
#make date a searchable feature:
cyano.combined$date2 <- as.POSIXlt(cyano.combined$date)
#Extract 2012 values from Green Lake
cyano.combined.GL.2012<-cyano.combined[cyano.combined$date2$year==112 & cyano.combined$lake=="GL",]

#sort by date
cyano.combined.GL.2012<-cyano.combined.GL.2012[order(cyano.combined.GL.2012$date),]

#make date into a factor instead of a date
cyano.combined.GL.2012$date<-as.factor(cyano.combined.GL.2012$date)

#make bin variable into a factor instead of being interpreted as a number... should've assigned a through e
cyano.combined.GL.2012$bin<-as.factor(cyano.combined.GL.2012$bin)

#plot it.  still not perfect, but it's somthing...
ggplot(cyano.combined.GL.2012,aes(date,log.cyano.ml,fill=bin))+
  geom_bar(stat="identity",position="dodge")

## Plotting graph of log cyanophage types per year per lake----

#GL 2012
layout(rbind(1,2), heights=c(15,1))
with(all.2012, plot(date,log.rlcp1,type="n", xlab="Date",ylab="log(copies/ml)",main="Cyanophage Genotype Abundance GL 2012"))
with(gl.2012,lines(date,log.rlcp1, col="lightblue"))
with(gl.2012, points(date,log.rlcp1, pch=18, col="lightblue"))
with(gl.2012,lines(date,log.rlcp2a, col="magenta"))
with(gl.2012, points(date,log.rlcp2a, pch=18, col="magenta"))
with(gl.2012,lines(date,log.rlcp4, col="aquamarine"))
with(gl.2012, points(date,log.rlcp4, pch=18, col="aquamarine"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("rlcp1","rlcp2a","rlcp4"),pch=18,
       col=c("lightblue","magenta","aquamarine"),ncol=3,bty ="n")

dev.off()

#RL 2012
layout(rbind(1,2), heights=c(15,1))
with(all.2012, plot(date,log.rlcp1,type="n", xlab="Date",ylab="log(copies/ml)",main="Cyanophage Genotype Abundance RL 2012"))
with(rl.2012,lines(date,log.rlcp1, col="lightblue"))
with(rl.2012, points(date,log.rlcp1, pch=18, col="lightblue"))
with(rl.2012,lines(date,log.rlcp2a, col="magenta"))
with(rl.2012, points(date,log.rlcp2a, pch=18, col="magenta"))
with(rl.2012,lines(date,log.rlcp4, col="aquamarine"))
with(rl.2012, points(date,log.rlcp4, pch=18, col="aquamarine"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("rlcp1","rlcp2a","rlcp4"),pch=18,
       col=c("lightblue","magenta","aquamarine"),ncol=3,bty ="n")

dev.off()

#GL 2013
layout(rbind(1,2), heights=c(15,1))
with(all.2013, plot(date,log.rlcp2a,type="n", xlab="Date",ylab="log(copies/ml)",main="Cyanophage Genotype Abundance GL 2013"))
with(gl.2013,lines(date,log.rlcp1, col="lightblue"))
with(gl.2013, points(date,log.rlcp1, pch=18, col="lightblue"))
with(gl.2013,lines(date,log.rlcp2a, col="magenta"))
with(gl.2013, points(date,log.rlcp2a, pch=18, col="magenta"))
with(gl.2013,lines(date,log.rlcp4, col="aquamarine"))
with(gl.2013, points(date,log.rlcp4, pch=18, col="aquamarine"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("rlcp1","rlcp2a","rlcp4"),pch=18,
       col=c("lightblue","magenta","aquamarine"),ncol=3,bty ="n")

dev.off()

#RL 2013
layout(rbind(1,2), heights=c(15,1))
with(all.2013, plot(date,log.rlcp2a,type="n", xlab="Date",ylab="log(copies/ml)",main="Cyanophage Genotype Abundance RL 2013"))
with(rl.2013,lines(date,log.rlcp1, col="lightblue"))
with(rl.2013, points(date,log.rlcp1, pch=18, col="lightblue"))
with(rl.2013,lines(date,log.rlcp2a, col="magenta"))
with(rl.2013, points(date,log.rlcp2a, pch=18, col="magenta"))
with(rl.2013,lines(date,log.rlcp4, col="aquamarine"))
with(rl.2013, points(date,log.rlcp4, pch=18, col="aquamarine"))
par(mar=c(0, 0, 0, 0))
plot.new()
legend('center','groups',c("rlcp1","rlcp2a","rlcp4"),pch=18,
       col=c("lightblue","magenta","aquamarine"),ncol=3,bty ="n")

dev.off()


#Comparing Data Vectors... playing with stats----
with(alldata, plot(bac, vlp, col=lake,pch=19))
alldata$logbac<-with(alldata, log(bac, 10))
alldata$logvlp<-with(alldata, log(vlp,10))
with(alldata, plot(logbac, logvlp, col=lake, pch=19))
#looks like a relationship between viruses and bacteria... not surprising, let's apply a linear model:
lm<-lm(logbac~logvlp,data=alldata)
summary(lm)
hist(lm$residuals)
stres<-rstandard(lm)
pred<-predict(lm)
plot(pred,stres)

alldata$vbr<-alldata$vlp/alldata$bac
with(alldata, plot(date, vbr, col=lake, pch=19))
with(alldata, plot(cyano, vlp,col=lake,pch=19))
alldata$logcyano<-log(alldata$cyano)
alldata$logvlp<-log(alldata$vlp)
with(alldata, plot(logcyano, vlp, col=lake, pch=19))

with(alldata, plot(logrlcp1,logbin2, col=lake,pch=19))
alldata$logbin2<-log(alldata$bin2.cyanos, 10)

with(alldata,plot(logrlcp1,logbin3, col=lake,pch=19))
alldata$logrlcp1<-log(alldata$rlcp1)
alldata$logbin3<-log(alldata$bin3.cyanos)
lm1<-lm(logbin3~logrlcp1,data=alldata)

##to add a legend to any of these...----
legend("topright",c("GL","RL"),pch=19,col=c("green","blue"))

##Statistics....----
num.data.all<-alldata[,c("temp","do","salinity","conductivity","pH","rlcp1","rlcp2a","rlcp4","bac","vlp","cyano","nitrate","nitrite","srp","chla","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bin5.cyanos")]
num.mat<-as.matrix(num.data.all)

hist(alldata$nitrate)
hist(alldata$srp)
hist(alldata$lake~alldata$nitrate)
table(alldata$lake,alldata$nitrate)
table(alldata$nitrate,alldata$lake)
boxplot(alldata$nitrate~alldata$lake)
boxplot(alldata$bac~alldata$lake+alldata$Year)
boxplot(alldata$vlp~alldata$lake+alldata$Year)


vlp.anova<-aov(alldata$vlp~alldata$lake+alldata$Year)
bac.anova<-aov(alldata$bac~alldata$lake+alldata$Year)
vlp.anova
bac.anova
