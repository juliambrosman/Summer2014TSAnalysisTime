##import qPCR data
qpcr.data=read.table("TimeSeriesqPCR_data.csv", sep=",",header=TRUE)
str(qpcr.data)

##read data data as a date:
qpcr.data$date <- as.Date(qpcr.data$date, "%m/%d/%y")
View(qpcr.data)

#calculate copies per sample and adjust sd
qpcr.data$rlcp1.smpl<-(qpcr.data$rlcp1.rxn.copies*qpcr.data$elute.vol)/qpcr.data$filt.vol
qpcr.data$rlcp1.smpl.sd<-(qpcr.data$rlcp1.rxn.stdev*qpcr.data$elute.vol)/qpcr.data$filt.vol
qpcr.data$rlcp2a.smpl<-(qpcr.data$rlcp2a.rxn.copies*qpcr.data$elute.vol)/qpcr.data$filt.vol
qpcr.data$rlcp2a.smpl.sd<-(qpcr.data$rlcp2a.rxn.stdev*qpcr.data$elute.vol)/qpcr.data$filt.vol
qpcr.data$rlcp4.smpl<-(qpcr.data$rlcp4.rxn.copies*qpcr.data$elute.vol)/qpcr.data$filt.vol
qpcr.data$rlcp4.smpl.sd<-(qpcr.data$rlcp4.rxn.stdev*qpcr.data$elute.vol)/qpcr.data$filt.vol

##Load dplyr to start grouping data
library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

#group data by sampling location and time
quant.groups<-group_by(qpcr.data,date,lake)


df$df <- df$n-1
pooledSD <- sqrt( sum(df$sd^2 * df$df) / sum(df$df) )
per.locale<-summarise(quant.groups,count=n(), rlcp1=mean(rlcp1.smpl),
                      rlcp1.sd=sqrt(sum(rlcp1.smpl.sd^2*3)/((count-1)*3)))
View(per.locale)
          

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
