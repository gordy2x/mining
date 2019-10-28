
rain_moist=read.csv("all_probes_2012_2017.csv", header=TRUE)
library("ggplot2")
library(cowplot)
library(dplyr)
library(coxme)
library(rms)
library(gridExtra)
library(tidyr)

#define probe 
rain_moist$Probe=with(rain_moist,paste0(Swamp,Probe_no,Probe_depth,veg_type,Impact))
probes=unique(rain_moist$Probe)

for(probe in probes){
  
  #which rows of data correspond to this probe
  rel.rows=rain_moist$Probe==probe
  #extract just those rows
  dat.loc=rain_moist[rel.rows,]
  #when is the rain greater than 25.4
  rain_event=dat.loc$rainfall_daily_total>30.5
  
  #this code extracts the start of a rain event
  events=diff(rain_event)
  day_before_rain=which(events==1)
  
  #if this probe has no rain events, write a data frame with one row
  if(length(day_before_rain)==0){
    event=data.frame(counted=FALSE,reason="No rain events",rel=NA,
                     date=NA,
                     dbr=NA,ldr=NA,dbnr=NA,
                     Swamp=unique(dat.loc$Swamp),Probe_no=unique(dat.loc$Probe_no),
                     Probe_depth=unique(dat.loc$Probe_depth),
                     veg_type=unique(dat.loc$veg_type),
                     Impact=unique(dat.loc$Impact),
                     soil_moist_day_before=NA,days.above.50pcSM=NA,
                     days.above.25pcSM=NA,days.above.75pcSM=NA,
                     da50C=FALSE,da25C=FALSE,da75C=FALSE,
                     n.days.rain=NA,n.days.no.rain=NA,average.days.no.rain=NA,total.rain.volume=NA)

    write.csv(event,paste0("probes/",probe,".csv"))
    next
  }
  
  #otherwise create dataset of this probe
  event=data.frame(counted=TRUE,reason=NA,rel=NA,
                   date=dat.loc$date[day_before_rain],
                   dbr=day_before_rain,ldr=NA,dbnr=NA,
                   Swamp=unique(dat.loc$Swamp),Probe_no=unique(dat.loc$Probe_no),
                   Probe_depth=unique(dat.loc$Probe_depth),
                   veg_type=unique(dat.loc$veg_type),
                   Impact=unique(dat.loc$Impact),
                   soil_moist_day_before=NA,days.above.50pcSM=NA,
                   days.above.25pcSM=NA,days.above.75pcSM=NA,
                   da50C=TRUE,da25C=TRUE,da75C=TRUE,
                   n.days.rain=NA,n.days.no.rain=NA,average.days.no.rain=NA,total.rain.volume=NA)
  event$rel=which(rel.rows==T)[1]
  
  #now look inside this probe
  #for each rain event
  for(i in 1:length(day_before_rain)){
    #cycle though data
    dbr=day_before_rain[i]
    dontstop=T
    k=dbr
    while(dontstop){
      k=k+1
      #if missing rain data found, stop, and record not to count this event
      if(is.na(events[k])){
        event$counted[i]=FALSE
        event$reason[i]="Missing rain data"
        dontstop=F
        #if we reach the end of the data, stop
      }else if(k==length(events)){
        event$counted[i]=FALSE
        event$reason[i]="End of data"
        dontstop=F
      #if events is -1, this means this is the last day of rain
      }else if(events[k]==(-1)){
        event$ldr[i]=k #extract last day of rain
      #if events is 1 then this is the day before the next rain
      }else if(events[k]==1){
        event$dbnr[i]=k #extract day before new rain
        dontstop=F
      }
    }
    
  }
  
  #save
  write.csv(event,paste0("probes/",probe,".csv"))
}


#for all those probes datasets saved above
for(probe in probes){
  #read in data
  rel.rows=rain_moist$Probe==probe
  dat.loc=rain_moist[rel.rows,]
  event=read.csv(paste0("probes/",probe,".csv"),stringsAsFactors=FALSE)
  
  for(i in 1:dim(event)[1]){
    #if event is counted
    if(event$counted[i]==TRUE){
      #extract soil moisture the day before rain
      event$soil_moist_day_before[i]=dat.loc$soil_moisture[event$dbr[i]]
      between_rain=(event$ldr[i]+1):(event$dbnr[i])
      
      #extract moisture between rain events
      moisture_loc=dat.loc$soil_moisture[between_rain]
      
      #if any missing data, record as not counted event
      if(any(is.na(moisture_loc))){
        event$counted[i]=FALSE
        event$reason[i]="Missing soil moisture"
        next
      }
      #count number of days above 50% , 25%  and 75%  
      mg=moisture_loc>50
      if(sum(mg)==length(mg)){
        event$days.above.50pcSM[i]=sum(mg)
        event$da50C[i]=FALSE #false is alive
      } else{
        event$days.above.50pcSM[i]=min(which(mg==FALSE)-1)
      }
      
      mg=moisture_loc>25
      if(sum(mg)==length(mg)){
        event$days.above.25pcSM[i]=sum(mg)
        event$da25C[i]=FALSE
      } else{
        event$days.above.25pcSM[i]=min(which(mg==FALSE)-1)
      } 
      
      mg=moisture_loc>75
      if(sum(mg)==length(mg)){
        event$days.above.75pcSM[i]=sum(mg)
        event$da75C[i]=FALSE
      } else{
        event$days.above.75pcSM[i]=min(which(mg==FALSE)-1)
      }
      
      #record number of days of rain, no rain and rain volume
      event$n.days.rain[i]=event$ldr[i]-event$dbr[i]
      event$n.days.no.rain[i]=event$dbnr[i]-event$ldr[i]
      event$total.rain.volume[i]=sum(dat.loc$rainfall_daily_total[(event$dbr[i]+1):event$ldr[i]])
    }
    event$average.days.no.rain=mean(event$n.days.no.rain,na.rm=T)
    
  }
  write.csv(event,paste0("probes/",probe,".csv"))
}


#combine into one dataset
all_probes=event[0,]
for(probe in probes){
event=read.csv(paste0("probes/",probe,".csv"),stringsAsFactors=FALSE)
all_probes=rbind(all_probes,event)
}



#add NARCLIM data
NARCLIM=read.csv("NARCLIM.csv",stringsAsFactors=F)
colnames(NARCLIM)[3]="Probe_no"
NARCLIM$Swamp[NARCLIM$Swamp=="Leech"]="leech"
NARCLIM$Swamp[NARCLIM$Swamp=="1a"]="den_1a"
NARCLIM$Swamp[NARCLIM$Swamp=="5"]="den_5"
NARCLIM$Swamp=factor(NARCLIM$Swamp)
probes=left_join(all_probes,NARCLIM)

#add swamp and catchment area data
area=read.csv("area.csv")
swamp_area=area[area$Type=="Swamp",-2]
colnames(swamp_area)[1]="swamp_area"
catchment_area=area[area$Type=="Catchment",-2]
colnames(catchment_area)[1]="catchment_area"
probes=left_join(probes,swamp_area)
probes=left_join(probes,catchment_area)

#add dist from escarpment
Scarpdist=read.csv("Probe_Scarp_Dist.csv",stringsAsFactors = FALSE)
Scarpdist=dplyr::select(Scarpdist,Swamp,Probechann,scarpdist)
colnames(Scarpdist)[2]="Probe_no"
Scarpdist$Swamp[Scarpdist$Swamp=="Leech"]="leech"
Scarpdist$Swamp[Scarpdist$Swamp=="1a"]="den_1a"
Scarpdist$Swamp[Scarpdist$Swamp=="5"]="den_5"
unique(Scarpdist$Swamp)
unique(probes$Swamp)

probes=left_join(probes,Scarpdist)

probes=probes[,-c(1:2)]
write.csv(probes,paste0("probes/all_probes.csv"))




