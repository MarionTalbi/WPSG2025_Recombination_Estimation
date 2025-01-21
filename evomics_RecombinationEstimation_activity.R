#################################################
########          pyrho Activity.        ########
######    Marion Talbi - 10th of January   ######
#################################################
setwd("~/workshop_materials/21_recombination_estimation")

#################################################
#Load the data
maps_Ben_bp15<-read.table("./pyrho/output_files/LS420024.2_benthic_bp15_pyrho.optimize",col.names = c("start","end","rate"))
maps_Lit_bp15<-read.table("./pyrho/output_files/LS420024.2_littoral_bp15_pyrho.optimize",col.names = c("start","end","rate"))
maps_Ben_bp30<-read.table("./pyrho/output_files/LS420024.2_benthic_bp30_pyrho.optimize",col.names = c("start","end","rate"))
maps_Lit_bp30<-read.table("./pyrho/output_files/LS420024.2_littoral_bp30_pyrho.optimize",col.names = c("start","end","rate"))

#################################################
#####1.  To convert r in cM / Mb

maps_Ben_bp15$"cM/Mb"<-(maps_Ben_bp15$rate/0.01)*10^6
maps_Lit_bp15$"cM/Mb"<-(maps_Lit_bp15$rate/0.01)*10^6

#Quick look at the recombination landscape
plot(y=maps_Ben_bp15$`cM/Mb`,x=(maps_Ben_bp15$start+maps_Ben_bp15$end)/2,type="l",lwd=1.5,col="steelblue4",ylim=c(0,(y_lim_max/0.01)*10^6),main="benthic - bp 15",ylab="r")

#################################################
##### To obtain the recombination map length
right_snp=5806
left_snp=2539
rec_rate_in_cM=0.455158616243
map_length=0.00457837
new_map_length=map_length+((right_snp-left_snp)*rec_rate_in_cM/10^6)
new_map_length

#################################################
#####2.  To plot recombination landscapes
par(mfrow=c(2,1))
y_lim_max=1e-06 #You can change this value

##Benthic
plot(y=maps_Ben_bp15$rate,x=(maps_Ben_bp15$start+maps_Ben_bp15$end)/2,type="l",lwd=1.5,col="steelblue4",xlab="Genomic position",ylim=c(0,y_lim_max),main="benthic - bp 15",ylab="r")
plot(y=maps_Ben_bp30$rate,x=(maps_Ben_bp30$start+maps_Ben_bp30$end)/2,type="l",lwd=1.5,col="steelblue4",xlab="Genomic position",ylim=c(0,y_lim_max),main="benthic - bp 30",ylab="r")

##Littoral
plot(y=maps_Lit_bp15$rate,x=(maps_Lit_bp15$start+maps_Lit_bp15$end)/2,type="l",lwd=1.5,col="gold",xlab="Genomic position",ylim=c(0,y_lim_max),main="littoral - bp 15",ylab="r")
plot(y=maps_Lit_bp30$rate,x=(maps_Lit_bp30$start+maps_Lit_bp30$end)/2,type="l",lwd=1.5,col="gold",xlab="Genomic position",ylim=c(0,y_lim_max),main="littoral - bp 30",ylab="r")

##Both pops
plot(y=maps_Ben_bp15$rate,x=(maps_Ben_bp15$start+maps_Ben_bp15$end)/2,type="l",lwd=1.5,col="steelblue4",xlab="Genomic position",ylim=c(0,y_lim_max),main="benthic - bp 15",ylab="r")
plot(y=maps_Lit_bp15$rate,x=(maps_Lit_bp15$start+maps_Lit_bp15$end)/2,type="l",lwd=1.5,col="gold",xlab="Genomic position",ylim=c(0,y_lim_max),main="littoral - bp 15",ylab="r")

#################################################
#####3.  Extra activity - Cumulative Curve
ecdf_Lit_bp15<-maps_Lit_bp15
ecdf_Lit_bp15$Rate_pond<-maps_Lit_bp15$rate*(maps_Lit_bp15$end-maps_Lit_bp15$start)
ecdf_Lit_bp15<-ecdf_Lit_bp15[order(ecdf_Lit_bp15$rate,decreasing = T),] ##We order the rate by decreasing order here
ecdf_Lit_bp15$cum_sum<-cumsum(ecdf_Lit_bp15$end-ecdf_Lit_bp15$start)
ecdf_Lit_bp15$cum_sum_pond<-cumsum(ecdf_Lit_bp15$Rate_pond)

ecdf_Ben_bp15<-maps_Ben_bp15
ecdf_Ben_bp15$Rate_pond<-maps_Ben_bp15$rate*(maps_Ben_bp15$end-maps_Ben_bp15$start)
ecdf_Ben_bp15<-ecdf_Ben_bp15[order(ecdf_Ben_bp15$rate,decreasing = T),]
ecdf_Ben_bp15$cum_sum<-cumsum(ecdf_Ben_bp15$end-ecdf_Ben_bp15$start)
ecdf_Ben_bp15$cum_sum_pond<-cumsum(ecdf_Ben_bp15$Rate_pond)

SubData_ecdf_Lit_bp15<-ecdf_Lit_bp15[seq(1,length(ecdf_Lit_bp15$cum_sum),1000),]
SubData_ecdf_Ben_bp15<-ecdf_Ben_bp15[seq(1,length(ecdf_Ben_bp15$cum_sum),1000),]

par(mfrow=c(1,1),mai=c(1,1,0.5,0.5))
plot(SubData_ecdf_Lit_bp15$cum_sum/max(SubData_ecdf_Lit_bp15$cum_sum),SubData_ecdf_Lit_bp15$cum_sum_pond/max(SubData_ecdf_Lit_bp15$cum_sum_pond,na.rm=T),
     type='l',col="gold",lwd=2,lty=1,cex=2,xlab="Proportion of the Genome",ylab="Proportion of recombination event",main="Cumulative curve")
lines(SubData_ecdf_Ben_bp15$cum_sum/max(SubData_ecdf_Ben_bp15$cum_sum),SubData_ecdf_Ben_bp15$cum_sum_pond/max(SubData_ecdf_Ben_bp15$cum_sum_pond,na.rm=T),col="steelblue",type='l',lwd=2,lty=1,cex=2)

