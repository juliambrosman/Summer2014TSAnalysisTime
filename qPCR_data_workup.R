##import qPCR data
qpcr.data=read.table("TimeSeriesqPCR_data.csv", sep=",",header=TRUE)

##read data data as a date:
qpcr.data$date <- as.Date(qpcr.data$date, "%m/%d/%y")

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

#Calculating mean abundance per timepoint and location and pooled variance according to a calculation I found on the internet here:
#http://stackoverflow.com/questions/16974389/how-to-calculate-a-pooled-standard-deviation-in-r
#Example pooled standard deviation from which I'm working:
#df$df <- df$n-1
#pooledSD <- sqrt( sum(df$sd^2 * df$df) / sum(df$df) )
per.locale<-summarise(quant.groups,count=n(), rlcp1=mean(rlcp1.smpl),
                      rlcp1.sd=sqrt(sum((rlcp1.smpl.sd^2)*(3-1))/((count*3)-1)), 
                      rlcp2a=mean(rlcp2a.smpl),
                      rlcp2a.sd=sqrt(sum((rlcp2a.smpl.sd^2)*(3-1))/((count*3)-1)),
                      rlcp4=mean(rlcp4.smpl),
                      rlcp4.sd=sqrt(sum((rlcp4.smpl.sd^2)*(3-1))/((count*3)-1)))
View(per.locale)

#per.local now holds all the data
          

