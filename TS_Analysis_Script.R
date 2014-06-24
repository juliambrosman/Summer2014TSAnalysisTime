##First Look at ARISA Data##----
arisa.data=read.table("20140613_ARISA_Updated.csv", sep=",",header=TRUE)
arisa.data$Size<-round(arisa.data$Size, digits=0)
plot(table(arisa.data$Size), xlab="fragment length (bp)", ylab="# Samples with Peak")

#take all measurements greater than 1000 bp to look more closely at data
arisa.data.big<-subset(arisa.data, Size>=1000)

#same wiht 1100bp
arisa.data.big2<-subset(arisa.data, Size>=1100)
plot(table(arisa.data.big2$Size), xlab="fragment length (bp)", ylab="# Samples with Peak")
big2table<-table(arisa.data.big2$Size)


#Looked at this table/plot and determined five different groups of cyanobacterial peaks:
#1108bp to 1125bp-->group 1  Sequence Group 1129 to 1137
#1130bp to 1142bp-->group 2  Sequence Group 1146 to 1148
#1211bp to 1240bp-->group 3  Sequence group with no sequenced reps
#1245bp to 1304bp-->group 4  Sequence Group 1300 to 1310
#1310bp to 1336bp-->group 5  Sequence Group 1350-1360bp

#So now to pull out only peaks within these different groups:
arisa.data.big2$size.bin[arisa.data.big2$Size>=1108 & arisa.data.big2$Size<=1125]=1
arisa.data.big2$size.bin[arisa.data.big2$Size>=1130 & arisa.data.big2$Size<=1142]=2
arisa.data.big2$size.bin[arisa.data.big2$Size>=1211 & arisa.data.big2$Size<=1240]=3
arisa.data.big2$size.bin[arisa.data.big2$Size>=1245 & arisa.data.big2$Size<=1304]=4
arisa.data.big2$size.bin[arisa.data.big2$Size>=1310 & arisa.data.big2$Size<=1336]=5

#remove rows with 'na' in size.bin column
arisa.binned<-arisa.data.big2[!(is.na(arisa.data.big2$size.bin)),]

#load dplyr package for data manipulation:----
library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

#Trying to make a Summary Data Matrix

summary<-group_by(arisa.binned,date,lake)

total.area<-summarise(summary,count=n(),area=sum(Area.in.BP))

#separate the different bin sizes 
bin1<-filter(arisa.binned, size.bin==1)
bin2<-filter(arisa.binned, size.bin==2)
bin3<-filter(arisa.binned, size.bin==3)
bin4<-filter(arisa.binned, size.bin==4)
bin5<-filter(arisa.binned, size.bin==5)

summary(bin1$Size)
summary(bin2$Size)
summary(bin3$Size)
summary(bin4$Size)
summary(bin5$Size)

#group each bin by date and lake
bin1grouped<-group_by(bin1,date,lake)
bin2grouped<-group_by(bin2,date,lake)
bin3grouped<-group_by(bin3,date,lake)
bin4grouped<-group_by(bin4,date,lake)
bin5grouped<-group_by(bin5,date,lake)

#calculate total area for each bin at each time point and location
bin1.area<-summarise(bin1grouped,bin1.area=sum(Area.in.BP))
bin2.area<-summarise(bin2grouped,bin2.area=sum(Area.in.BP))
bin3.area<-summarise(bin3grouped,bin3.area=sum(Area.in.BP))
bin4.area<-summarise(bin4grouped,bin4.area=sum(Area.in.BP))
bin5.area<-summarise(bin5grouped,bin5.area=sum(Area.in.BP))

#now merging all the documents into one for additional calculations...
mergetrial<-merge(total.area, bin1.area, by=c("date","lake"), all.x=TRUE)
mergetrial<-merge(mergetrial, bin2.area, by=c("date","lake"), all.x=TRUE)
mergetrial<-merge(mergetrial, bin3.area, by=c("date","lake"), all.x=TRUE)
mergetrial<-merge(mergetrial, bin4.area, by=c("date","lake"), all.x=TRUE)
mergetrial<-merge(mergetrial, bin5.area, by=c("date","lake"), all.x=TRUE)

#Create column of fraction of representation of each of the peaks
mergetrial$bin1.rep<-(mergetrial$bin1.area/mergetrial$area)
mergetrial$bin2.rep<-(mergetrial$bin2.area/mergetrial$area)
mergetrial$bin3.rep<-(mergetrial$bin3.area/mergetrial$area)
mergetrial$bin4.rep<-(mergetrial$bin4.area/mergetrial$area)
mergetrial$bin5.rep<-(mergetrial$bin5.area/mergetrial$area)

Relative.arisa<-mergetrial[,c(1:2,10:14)]

##Now bring in all TS Metadata##----
alldata=read.table("AllTSMetaData.csv", sep=",",header=TRUE)

#remove empty columns
alldata$X<-NULL
alldata$X.1<-NULL
alldata$X.2<-NULL
#Remove empty rows
alldata<-alldata[1:51,]

cyano.counts<-alldata[,c(1,5,21,22)]

cyano.counts.arisa<-merge(cyano.counts, Relative.arisa, by=c("date","lake"), all.x=TRUE)

cyano.counts.arisa$bin1.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin1.rep
cyano.counts.arisa$bin2.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin2.rep
cyano.counts.arisa$bin3.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin3.rep
cyano.counts.arisa$bin4.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin4.rep
cyano.counts.arisa$bin5.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin5.rep

AbundCyano<-cyano.counts.arisa[,c(1,2,10:14)]

#Merge cyano bin counts with the rest of the meta data 
alldata<-merge(alldata, AbundCyano, by=c("date","lake"), all.x=TRUE)

##Cyano Abundance Data is now combined with all other data in the data table, next small fixes## ----

alldata$date <- as.Date(alldata$date, "%m/%d/%y") 
alldata$lake<-droplevels(alldata$lake)
# Remove October observations
alldata<-alldata[alldata$Month!=10,]
str(alldata)
#log of cyano abundance:
alldata$bin1.cyanos<-log(alldata$bin1.cyanos, 10)
alldata$bin2.cyanos<-log(alldata$bin2.cyanos, 10)
alldata$bin3.cyanos<-log(alldata$bin3.cyanos, 10)
alldata$bin4.cyanos<-log(alldata$bin4.cyanos, 10)
alldata$bin5.cyanos<-log(alldata$bin5.cyanos, 10)

##playing with the data##----
#boxplots: -----
boxplot(alldata$bac~alldata$Year+alldata$lake, main="bac")
boxplot(alldata$bac~alldata$Year+alldata$lake, main="bac")
boxplot(alldata$cyano~alldata$Year+alldata$lake, main="cyano")
boxplot(alldata$chla~alldata$Year+alldata$lake, main="chla")
boxplot(alldata$nitrate~alldata$Year+alldata$lake, main="nitrate")
boxplot(alldata$srp~alldata$Year+alldata$lake, main="srp")
boxplot(alldata$temp~alldata$Year+alldata$lake, main="temp")

#Setting up subsets to plot -----
rl.2012<-filter(alldata,Year==2012, lake=="RL")
rl.2013<-filter(alldata,Year==2013, lake=="RL")
gl.2012<-filter(alldata, Year==2012, lake=="GL")
gl.2013<-filter(alldata, Year==2013, lake=="GL")
all.2012<-filter(alldata,Year==2012)
all.2013<-filter(alldata,Year==2013)
palette(c("green","blue"))

##time series plots of counts/ml per year----
#Library with errbars function:
library("Hmisc", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

#vlp 2012
with(all.2012, errbar(date,vlp,
                      yplus=vlp+vlp.sd,yminus=vlp-vlp.sd,
                      col=lake,pch=19, main="vlp 2012", xlab="date",ylab="vlp/ml"))
title(main="vlp 2012")
?plot
#vlp 2013
with(all.2013, errbar(date,vlp,
                      yplus=vlp+vlp.sd,yminus=vlp-vlp.sd,
                      col=lake,pch=19, xlab="date",ylab="vlp/ml"))
title(main="vlp 2013")
#bac 2012
with(all.2012, errbar(date,bac,
                      yplus=bac+bac.sd,yminus=bac-bac.sd,
                      col=lake,pch=19, xlab="date",ylab="bac/ml"))
title(main="bac 2012")
#bac 2013
with(all.2013, errbar(date,bac,
                      yplus=bac+bac.sd,yminus=bac-bac.sd,
                      col=lake,pch=19, main="bac 2013", xlab="date",ylab="bac/ml"))
title(main="bac 2013")
#cyano 2012
with(all.2012, errbar(date,cyano,
                      yplus=cyano+cyano.sd,yminus=cyano-cyano.sd,
                      col=lake,pch=19, main="cyano 2012", xlab="date",ylab="cyano/ml"))
title(main="cyano 2012")

#cyano 2013
with(all.2013, errbar(date,cyano,
                      yplus=cyano+cyano.sd,yminus=cyano-cyano.sd,
                      col=lake,pch=19, main="cyano 2013", xlab="date",ylab="cyano/ml"))
title(main="cyano 2013")

#plotting cyanophage genotype abundance over time----

#rlcp1 2013
with(all.2013, errbar(date,rlcp1,
                      yplus=rlcp1+rlcp1.sd,yminus=rlcp1-rlcp1.sd,
                      col=lake,pch=19, main="rlcp1 2013", xlab="date",ylab="rlcp1/ml"))
title(main="rlcp1 2013")
#rlcp1 2012
with(all.2012, errbar(date,rlcp1,
                      yplus=rlcp1+rlcp1.sd,yminus=rlcp1-rlcp1.sd,
                      col=lake,pch=19, main="rlcp1 2013", xlab="date",ylab="rlcp1/ml"))
title(main="rlcp1 2012")

#rlcp2a 2013
with(all.2013, errbar(date,rlcp2a,
                      yplus=rlcp2a+rlcp2a.sd,yminus=rlcp2a-rlcp2a.sd,
                      col=lake,pch=19, main="rlcp2a 2013", xlab="date",ylab="rlcp2a/ml"))
title(main="rlcp2a 2013")
#rlcp2a 2012
with(all.2012, errbar(date,rlcp2a,
                      yplus=rlcp2a+rlcp2a.sd,yminus=rlcp2a-rlcp2a.sd,
                      col=lake,pch=19, main="rlcp2a 2013", xlab="date",ylab="rlcp2a/ml"))
title(main="rlcp2a 2012")

#rlcp4 2013
with(all.2013, errbar(date,rlcp4,
                      yplus=rlcp4+rlcp4.sd,yminus=rlcp4-rlcp4.sd,
                      col=lake,pch=19, main="rlcp4 2013", xlab="date",ylab="rlcp4/ml"))
title(main="rlcp4 2013")
#rlcp4 2012
with(all.2012, errbar(date,rlcp4,
                      yplus=rlcp4+rlcp4.sd,yminus=rlcp4-rlcp4.sd,
                      col=lake,pch=19, main="rlcp4 2013", xlab="date",ylab="rlcp4/ml"))
title(main="rlcp4 2012")

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
with(all.2012, plot(date,bin1.cyanos,
                    col=lake,pch=19, main="bin1.cyanos 2012", xlab="date",ylab="log(cells/ml)"))
#bin1.cyanos 2013
with(all.2013, plot(date,bin1.cyanos,
                    col=lake,pch=19, main="bin1.cyanos 2013", xlab="date",ylab="log(cells/ml)"))
#bin2.cyanos 2012
with(all.2012, plot(date,bin2.cyanos,
                    col=lake,pch=19, main="bin2.cyanos 2012", xlab="date",ylab="log(cells/ml)"))
#bin2.cyanos 2013
with(all.2013, plot(date,bin2.cyanos,
                    col=lake,pch=19, main="bin2.cyanos 2013", xlab="date",ylab="log(cells/ml)"))

#bin3.cyanos 2012
with(all.2012, plot(date,bin3.cyanos,
                    col=lake,pch=19, main="bin3.cyanos 2012", xlab="date",ylab="log(cells/ml)"))
#bin3.cyanos 2013
with(all.2013, plot(date,bin3.cyanos,
                    col=lake,pch=19, main="bin3.cyanos 2013", xlab="date",ylab="log(cells/ml)"))

#bin4.cyanos 2012
with(all.2012, plot(date,bin4.cyanos,
                    col=lake,pch=19, main="bin4.cyanos 2012", xlab="date",ylab="log(cells/ml)"))
#bin4.cyanos 2013
with(all.2013, plot(date,bin4.cyanos,
                    col=lake,pch=19, main="bin4.cyanos 2013", xlab="date",ylab="log(cells/ml)"))
#bin5.cyanos 2012
with(all.2012, plot(date,bin5.cyanos,
                    col=lake,pch=19, main="bin5.cyanos 2012", xlab="date",ylab="log(cells/ml)"))
#bin5.cyanos 2013
with(all.2013, plot(date,bin5.cyanos,
                    col=lake,pch=19, main="bin5.cyanos 2013", xlab="date",ylab="log(cells/ml)"))

##Plotting Bar Graph of all Cyano Genotypes per Lake:

##pull out just cyano representation:
cyano.bins<-alldata[,c(1,2,29:33)]

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

cyano.b5<-cyano.bins[,c(1,2,7)]
cyano.b5$bin<-5
colnames(cyano.b5)[3]<-"log.cyano.ml"

#bind rows for all the results into one data matrix
cyano.combined<-rbind(cyano.b1, cyano.b2, cyano.b3, cyano.b4, cyano.b5)
View(cyano.combined)

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

#Comparing Data Vectors----
with(alldata, plot(bac, vlp, col=lake,pch=19))
with(alldata, plot(chla,vlp,col=lake,pch=19))
with(alldata, plot(bac,chla, col=lake,pch=19))
with(alldata,plot(temp,vlp, col=lake,pch=19))

##to adda legend to any of these...----
legend("topright",c("GL","RL"),pch=19,col=c("green","blue"))

##Statistics....----
num.data.all<-alldata[,c("temp","do","salinity","conductivity","pH","rlcp1","rlcp2a","rlcp4","bac","vlp","cyano","nitrate","nitrite","srp","chla","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bin5.cyanos")]
num.mat<-as.matrix(num.data.all)

hist(alldata$nitrate)
hist(alldata$srp)
hist(alldata$nitrate~alldata$lake)
table(alldata$lake,alldata$nitrate)
table(alldata$nitrate,alldata$lake)
boxplot(alldata$nitrate~alldata$lake)
boxplot(alldata$bac~alldata$lake+alldata$Year)

vlp.anova<-aov(alldata$vlp~alldata$lake+alldata$Year)
bac.anova<-aov(alldata$bac~alldata$lake+alldata$Year)
vlp.anova
bac.anova
