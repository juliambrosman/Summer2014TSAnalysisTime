##playing with vegan to figure out how to visualize and analyze my time series data
#going to work with my own data to make some funky plots:
library(vegan)
library(MASS)

#load data, add a second date column
alldata<-read.table("All_Integrated_Output.csv",sep=",",header=TRUE)
alldata$date<-as.Date(alldata$date, "%Y-%m-%d")
alldata$date2<-as.POSIXlt(alldata$date)
alldata$Year<-as.factor(alldata$Year)
alldata$unique<-paste(alldata$lake,alldata$date,sep="_")
rownames(alldata)<-alldata$unique

#separate out different variables:
#Gather all biological data into one matrix to work from for biological correlations:
allbiodata<-alldata[,c("date","lake","Year","chla","bin1.rep","bin2.rep","bin3.rep","bin4.rep","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","bac","vlp","cyano","rlcp1","rlcp2a","rlcp4")]
allphysdata<-alldata[,c("date","lake","Year","temp","do","salinity","conductivity","pH","nitrate","nitrite","srp")]
allcyanos<-alldata[,c("date","lake","Year","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos")]
allcyphage<-alldata[,c("date","lake","Year","rlcp1","rlcp2a","rlcp4")]
allcounts<-alldata[,c("date","lake","Year","bac","vlp","cyano")]
allcy<-alldata[c("date","lake","Month","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","rlcp1","rlcp2a","rlcp4")]
allbiomat<-as.matrix(allbiodata[,9:18])
str(alldata)
#comparing lakes----
allcy.tran<-as.matrix(t(allcy[,3:9]))
allcy.tran.cor<-cor(allcy.tran, method="spearman")
View(allcy.tran.cor)
library(corrplot)
corrplot(allcy.tran.cor,method="circle")

#play with mds plots of just the cyanobacterial and cyanophage genotypes:----

allcym<-as.matrix(allcy[,3:9])
allcym.dist<-vegdist(allcym,na.rm=TRUE)
allcy.mds0<-isoMDS(allcym.dist)
ordiplot(allcy.mds0, type="t")
allcy.mds1<-metaMDS(allcym.dist, trace=FALSE)
plot(allcy.mds1,type="t")
pro2<-procrustes(allcy.mds0, allcy.mds1)
plot(pro2)
#now PCA:
?rda
allcy.pca<-rda(X=allcym, scale=TRUE)
plot(allcy.pca, scale=3)
allcy.dca<-decorana(allcym)
plot(allcy.dca)
pd<-as.matrix(allphysdata[2:11])
cca1<-cca(allcym~pd$temp+alldata$lake)
cca1.plot<- plot(cca1, choices=c(1,2))

allcy.ccora<-CCorA(pd, allcym)

allbio.dist<-vegdist(allbiomat, na.rm=TRUE)
allbio.mds1<-isoMDS(allbio.dist)
stressplot(allbio.mds1,allbio.dist)
ordiplot(allbio.mds1, type="t")

#make big ordination plot with all values considered, must remove "na"
allbiomat<-na.omit(allbiomat)
allbio.mds<-metaMDS(allbiomat,trace=FALSE)
plot(allbio.mds, type="t")
abdis<-wisconsin(sqrt(allbiomat))
abdist<-vegdist(abdis)
allbio.mds1<-isoMDS(abdist,trace=0)
pro1<-procrustes(allbio.mds,allbio.mds1)
pro1
plot(pro1)
plot(pro1, kind=2)

allbiopca<-rda(allbiomat, scale=TRUE)
allbiopca
plot(allbiopca, scaling=3)
biplot(allbiopca, scaling=3)

allbio.ca<-cca(allbiomat)
allbio.ca
plot(allbio.ca, scaling=1)

allbio.dca<-decorana(allbiomat)
allbio.dca
plot(allbio.dca, display="sites")

##Trying DFA analysis based on Fuhrman et al 2006 ----
fit <- lda(G ~ x1 + x2 + x3, data=mydata,
           na.action="na.omit", CV=TRUE)

# SO, need to pull out all data by, say, month.  Added "Month" to allcy
allcy<-alldata[,c("Month","bin1.cyanos","bin2.cyanos","bin3.cyanos","bin4.cyanos","rlcp1","rlcp2a","rlcp4")]
head(allcy)
str(allcy)
allcy$Month<-as.factor(allcy$Month)
fit<-lda(Month ~ rlcp1+rlcp2a+rlcp4, data=allcy, na.action="na.omit",CV=TRUE)
fit
?lda
ct <- table(allcy$Month, fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))
plot(fit)
plot(fit, dimen=1, type="both")
