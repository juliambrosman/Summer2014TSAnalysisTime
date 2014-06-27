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
vir.counts<-summarise(vir.groups,count=n(), vir=mean(cells.ml),
                      vir.sd=sqrt(sum((cells.ml.sd^2)*(10-1))/((count*10)-1)))
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
View(all.cell.counts)

#all.cell.counts now holds all the data

