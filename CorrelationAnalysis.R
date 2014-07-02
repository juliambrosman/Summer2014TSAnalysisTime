#Load Data from Integration script: 

alldata<-read.table("All_Integrated_Output.csv",sep=",",header=TRUE)
alldata$date<-as.Date(alldata$date, "%Y-%m-%d")
alldata$date<-as.POSIXlt(alldata$date)
alldata$Year<-as.factor(alldata$Year)

#Load dplyr package:
library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

#Gather all biological data into one matrix to work from for biological correlations:
allbiodata<-alldata[,c("date","lake","Year","chla","bin1.rep","bin2.rep","bin3.rep","bin4.rep","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bac","vlp","cyano","rlcp1","rlcp2a","rlcp4")]

#group data by Lake and Year:
allbiodata<-group_by(allbiodata, lake, Year)

# create dataframe of just 2012 data:
allbiodata2012<-filter(allbiodata, Year=="2012")

# Extract data into a matrix: **using cyano-normalized cyanotype data for first go
allbiomatrix2012<-as.matrix(allbiodata2012[,c("chla","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bac","vlp","cyano","rlcp1","rlcp2a","rlcp4")])

# Run correlation:
cor1<-cor(x=allbiomatrix2012, use="na.or.complete",method="spearman")
View(cor1)
# Looks to be a relationship between rlcp1 and cyanotype 2 and 3

# Now look at just individual lakes
allbiodata2012.gl<-filter(allbiodata, Year=="2012", lake=="GL")
allbiomatrix2012.gl<-as.matrix(allbiodata2012.gl[,c("chla","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bac","vlp","cyano","rlcp1","rlcp2a","rlcp4")])
cor.gl.2012<-cor(x=allbiomatrix2012.gl, use="na.or.complete",method="spearman")

allbiodata2012.rl<-filter(allbiodata, Year=="2012", lake=="RL")
allbiomatrix2012.rl<-as.matrix(allbiodata2012.rl[,c("chla","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bac","vlp","cyano","rlcp1","rlcp2a","rlcp4")])
cor.rl.2012<-cor(x=allbiomatrix2012.rl, use="na.or.complete",method="spearman")
View(cor.rl.2012)

# OK, now looking to skew 2012 virus data by one week to see how that affects results:
#sort by date:
allbiodata2012.gl<-allbiodata2012.gl[order(allbiodata2012.gl$date),]
allbiodata2012.gl.vir<-allbiodata2012.gl[,c("vlp","rlcp1","rlcp2a","rlcp4")]
allbiodata2012.gl.bac<-allbiodata2012.gl[,c("bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bac","cyano","chla")]
allbiodata2012.gl.vir$numbered<-c(1:17)
allbiodata2012.gl.bac$numbered<-c(2:18)

allbiodata.gl.virskewed<-merge(allbiodata2012.gl.vir,allbiodata2012.gl.bac,by="numbered")
View(allbiodata.gl.virskewed)
gl2012.virskewed.mat<-as.matrix(allbiodata.gl.virskewed[,2:12])

cor.gl2012.virskewed<-cor(gl2012.virskewed.mat, use="na.or.complete",method="spearman")
View(cor.gl2012.virskewed)
