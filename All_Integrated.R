### Working with bacterial and viral counts first:----
library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")

#load data:
bv.counts<-read.table("bac-vir_counts.csv", sep=",",header=TRUE)
bv.counts$date<-as.Date(bv.counts$date, "%m/%d/%y")

#define viewing area term:
va=976553.2
#define magnification:
m=100
#add a column for cells/ml:
bv.counts$cells.ml<-((bv.counts$average.per.box*va)/(bv.counts$filt.vol))*m/100
#add column for standard deviation adjustment:
bv.counts$cells.ml.sd<-((bv.counts$stdev.per.box*va)/(bv.counts$filt.vol))*m/100

#separate out bacterial and viral counts:
b.counts<-filter(bv.counts,cat=="b")
v.counts<-filter(bv.counts, cat=="v")


#now let's group the data based on date and lake:
bac.groups<-group_by(b.counts,date,lake)
vir.groups<-group_by(v.counts,date,lake)

#and summarize the averages and pooled standard deviations:
vir.counts<-summarise(vir.groups,count=n(), vlp=mean(cells.ml),
                      vlp.sd=sqrt(sum((cells.ml.sd^2)*(10-1))/((count*10)-1)))
bac.counts<-summarise(bac.groups,count=n(), bac=mean(cells.ml),
                      bac.sd=sqrt(sum((cells.ml.sd^2)*(10-1))/((count*10)-1)))
#remove count column:
vir.counts$count<-NULL
bac.counts$count<-NULL

### Now let's load the cyano data and do the same (ish) thing:
c.counts<-read.table("cyano_counts.csv", sep=",",header=TRUE)
c.counts$date<-as.Date(c.counts$date,"%m/%d/%y")

#because magnification was at 40x for cyanobacterial counts:
mc=40

#calculate cells/ml
c.counts$cells.ml<-((c.counts$cyano.boxcount*va)/(c.counts$filt.vol))*mc/100
c.counts$cells.ml.sd<-((c.counts$cyano.boxcount.sd*va)/(c.counts$filt.vol))*mc/100

#group the data based on lake and date:
cyanocount.groups<-group_by(c.counts,date,lake)

#and summarize mean and pooled standard deviation:
cyano.counts<-summarise(cyanocount.groups,count=n(), cyano=mean(cells.ml),
                      cyano.sd=sqrt(sum((cells.ml.sd^2)*(n.counts-1))/(sum(n.counts-1))))
#remove count column:
cyano.counts$count<-NULL

#now there are three files that have mean counts and pooled standard deviation data:
#cyano.counts
#bac.counts
#vir.counts

#Time to merge the files:
all.cell.counts<-merge(bac.counts, cyano.counts, by=c("date","lake"),all.x=TRUE, all.y=TRUE)
all.cell.counts<-merge(all.cell.counts,vir.counts,by=c("date","lake"),all.x=TRUE, all.y=TRUE)

#all.cell.counts now holds all the data

##Now load qPCR data----
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

#group data by sampling location and time
quant.groups<-group_by(qpcr.data,date,lake)

#Calculating mean abundance per timepoint and location and pooled variance according to a calculation I found on the internet here:
#http://stackoverflow.com/questions/16974389/how-to-calculate-a-pooled-standard-deviation-in-r
#Example pooled standard deviation from which I'm working:
#df$df <- df$n-1
#pooledSD <- sqrt( sum(df$sd^2 * df$df) / sum(df$df) )
per.locale<-summarise(quant.groups,count=n(), rlcp1=mean(rlcp1.smpl, na.rm=TRUE),
                      rlcp1.sd=sqrt(sum((rlcp1.smpl.sd^2)*(3-1), na.rm=TRUE)/((count*3)-1)), 
                      rlcp2a=mean(rlcp2a.smpl, na.rm=TRUE),
                      rlcp2a.sd=sqrt(sum((rlcp2a.smpl.sd^2)*(3-1), na.rm=TRUE)/((count*3)-1)),
                      rlcp4=mean(rlcp4.smpl, na.rm=TRUE),
                      rlcp4.sd=sqrt(sum((rlcp4.smpl.sd^2)*(3-1),na.rm=TRUE)/((count*3)-1)))

#per.local now holds all the data

#Now merge qPCR data with counts data:
counts.qpcr.data<-merge(all.cell.counts,per.locale,by=c("date","lake"),all.x=TRUE,all.y=TRUE)

##Incorporating ARISA Data##----
arisa.data=read.table("20140613_ARISA_Updated.csv", sep=",",header=TRUE)
arisa.data$Size<-round(arisa.data$Size, digits=0)
arisa.data$date<-as.Date(arisa.data$date, "%m/%d/%y")

# Only looking at peaks greater than 1100bp
arisa.data.big2<-subset(arisa.data, Size>=1100)
big2table<-table(arisa.data.big2$Size)

#Looked at this table/plot and determined five different groups of cyanobacterial peaks:
#1108bp to 1125bp-->group 1  Sequence Group 1129 to 1137
#1130bp to 1142bp-->group 2  Sequence Group 1146 to 1148
#1240bp to 1304bp-->group 3  Sequence Group 1300 to 1310
#1310bp to 1336bp-->group 4  Sequence Group 1350-1360bp

#So now to pull out only peaks within these different groups:
arisa.data.big2$size.bin[arisa.data.big2$Size>=1108 & arisa.data.big2$Size<=1125]=1
arisa.data.big2$size.bin[arisa.data.big2$Size>=1130 & arisa.data.big2$Size<=1142]=2
arisa.data.big2$size.bin[arisa.data.big2$Size>=1240 & arisa.data.big2$Size<=1304]=3
arisa.data.big2$size.bin[arisa.data.big2$Size>=1310 & arisa.data.big2$Size<=1336]=4

#remove rows with 'na' in size.bin column
arisa.binned<-arisa.data.big2[!(is.na(arisa.data.big2$size.bin)),]


#Make a grouped Summary Data Matrix

summary<-group_by(arisa.binned,date,lake)

total.area<-summarise(summary,count=n(),area=sum(Area.in.BP))

#separate the different bin sizes 
bin1<-filter(arisa.binned, size.bin==1)
bin2<-filter(arisa.binned, size.bin==2)
bin3<-filter(arisa.binned, size.bin==3)
bin4<-filter(arisa.binned, size.bin==4)

#group each bin by date and lake
bin1grouped<-group_by(bin1,date,lake)
bin2grouped<-group_by(bin2,date,lake)
bin3grouped<-group_by(bin3,date,lake)
bin4grouped<-group_by(bin4,date,lake)

#calculate total area for each bin at each time point and location
bin1.area<-summarise(bin1grouped,bin1.area=sum(Area.in.BP))
bin2.area<-summarise(bin2grouped,bin2.area=sum(Area.in.BP))
bin3.area<-summarise(bin3grouped,bin3.area=sum(Area.in.BP))
bin4.area<-summarise(bin4grouped,bin4.area=sum(Area.in.BP))

#now merging all the documents into one for additional calculations...
mergetrial<-merge(total.area, bin1.area, by=c("date","lake"), all.x=TRUE, all.y=TRUE)
mergetrial<-merge(mergetrial, bin2.area, by=c("date","lake"), all.x=TRUE, all.y=TRUE)
mergetrial<-merge(mergetrial, bin3.area, by=c("date","lake"), all.x=TRUE, all.y=TRUE)
mergetrial<-merge(mergetrial, bin4.area, by=c("date","lake"), all.x=TRUE, all.y=TRUE)

#Create column of fraction of representation of each of the peaks
mergetrial$bin1.rep<-(mergetrial$bin1.area/mergetrial$area)
mergetrial$bin2.rep<-(mergetrial$bin2.area/mergetrial$area)
mergetrial$bin3.rep<-(mergetrial$bin3.area/mergetrial$area)
mergetrial$bin4.rep<-(mergetrial$bin4.area/mergetrial$area)

Relative.arisa<-mergetrial[,c(1:2,9:12)]

cyano.counts.arisa<-merge(cyano.counts, Relative.arisa, by=c("date","lake"), all.x=TRUE, all.y=TRUE)

cyano.counts.arisa$bin1.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin1.rep
cyano.counts.arisa$bin2.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin2.rep
cyano.counts.arisa$bin3.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin3.rep
cyano.counts.arisa$bin4.cyanos<-cyano.counts.arisa$cyano*cyano.counts.arisa$bin4.rep

AbundCyano<-cyano.counts.arisa[,c(1,2,9:12)]

arisa.rel.norm<-merge(Relative.arisa,AbundCyano,by=c("date","lake"), all.x=TRUE, all.y=TRUE)
arisa.rel.norm[is.na(arisa.rel.norm)]<-0

#Merge the relative arisa and cyano-count normalized arisa data with other data:
counts.qpcr.arisa<-merge(arisa.rel.norm,counts.qpcr.data, by=c("date","lake"), all.x=TRUE, all.y=TRUE)

##Integrate physical data ----
#Now bring in physical data:
phys.data=read.table("physical_ts_data.csv", sep=",",header=TRUE)
phys.data$date<-as.Date(phys.data$date, "%m/%d/%y")

# Merge physical data with the larger table:
all.data<-merge(phys.data, counts.qpcr.arisa, by=c("date","lake"),all.x=TRUE, all.y=TRUE)
#adjust
all.data$count<-NULL
#biovar[is.na(biovar)]<-0
all.data[is.na(all.data[,17:24])]<-0

##Write table to file----
write.table(all.data, file="All_Integrated_Output.csv", sep=",")

