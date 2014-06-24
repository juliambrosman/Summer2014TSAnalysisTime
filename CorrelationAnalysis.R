##Integrate ARISA Data into Data Matrix:----
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

str(alldata)

###Correlation Analyses Attempts----
#extract biologial variables
biovar<-alldata[,c("rlcp1", "rlcp2a","rlcp4","bac","vlp","cyano","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bin5.cyanos")]

#remove NA's from counts data
biovar<-biovar[-which(is.na(c(biovar$vlp, biovar$bac, biovar$cyano))),]

#Change NA's in cyano counts to zeros
biovar[is.na(biovar)]<-0

#correlation analysis of all biological variables:
biomat<-as.matrix(biovar)
cor1<-cor(biomat, y=NULL, use="everything",method="spearman")
View(cor1)

##Individual Lakes/Time Correlations:
#remove null counts measurements
cor.data<-alldata[-which(is.na(c(alldata$vlp, alldata$bac, alldata$cyano))),]

#change "NA" cyanotype measurements to 0
cor.data$bin1.cyanos[is.na(cor.data$bin1.cyanos)]<-0
cor.data$bin2.cyanos[is.na(cor.data$bin2.cyanos)]<-0
cor.data$bin3.cyanos[is.na(cor.data$bin3.cyanos)]<-0
cor.data$bin4.cyanos[is.na(cor.data$bin4.cyanos)]<-0
cor.data$bin5.cyanos[is.na(cor.data$bin5.cyanos)]<-0

#set correlation sets to look at
cor.rl.2012<-filter(cor.data,Year==2012, lake=="RL")
cor.rl.2013<-filter(cor.data,Year==2013, lake=="RL")
cor.gl.2012<-filter(cor.data, Year==2012, lake=="GL")
cor.gl.2013<-filter(cor.data, Year==2013, lake=="GL")
cor.all.2012<-filter(cor.data,Year==2012)
cor.all.2013<-filter(cor.data,Year==2013)
cor.all.gl<-filter(cor.data, lake=="GL")
cor.all.rl<-filter(cor.data, lake=="RL")


cor.2012.matrix<-as.matrix(cor.all.2012[,c("do","salinity","conductivity","pH","rlcp1","rlcp2a","rlcp4",
                       "bac","vlp","cyano","nitrate","nitrite","srp","chla",
                       "bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bin5.cyanos")])
cor.2012all<-cor(cor.2012.matrix, y=NULL, use="complete.obs",method="spearman")
str(cor.2012all)

#comparing between lakes:----
cor.gl<-cor.all.gl[,c("date","do","salinity","conductivity","pH","rlcp1","rlcp2a","rlcp4",
                      "bac","vlp","cyano","nitrate","nitrite","srp","chla",
                      "bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bin5.cyanos")]
cor.rl<-cor.all.rl[,c("date","do","salinity","conductivity","pH","rlcp1","rlcp2a","rlcp4",
                      "bac","vlp","cyano","nitrate","nitrite","srp","chla",
                      "bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bin5.cyanos")]
#change col names to distinguish lake measurements:
colnames(cor.gl)[2:20]<-paste("gl",colnames(cor.gl)[2:20],sep="_")
colnames(cor.rl)[2:20]<-paste("rl",colnames(cor.rl)[2:20],sep="_")

#merge rl and gl values
sbs.rl.gl<-merge(cor.rl, cor.gl, by="date",all.x=TRUE)
head(sbs.rl.gl)
sbs.mat<-as.matrix(sbs.rl.gl[,2:39])

sbscor<-cor(sbs.mat, y=NULL, use="everything",method="spearman")
View(sbscor)
