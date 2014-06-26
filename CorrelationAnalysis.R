#First Run "Cyanotype integration" script

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
