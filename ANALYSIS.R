packages <-c('sf','terra','ggplot2','leastcostpath','lwgeom','dplyr','tidyr',
             'ggbeeswarm','hrbrthemes','ggridges','ggpubr','terrainr')
for(p in packages) if(p %in% rownames(installed.packages()) == F) { install.packages(p) }
for(p in packages) suppressPackageStartupMessages(library(p,quietly=T,character.only=T))



# set working directory

setwd("")
theme_set(theme_bw())
sf_use_s2(TRUE)

############### STEP 0: IMPORT DATA & GIVE STRUCTURE ###########################
study_area<- data.frame(lat=c(35.877, 35.877, 35.982,35.982),
                        long=c(-108.759, -108.543,-108.759,-108.543),
                        box=c("example","example","example","example"))%>%
  st_as_sf(coords=c("long","lat"),remove=F,crs=4326)%>%
  st_bbox()%>%
  st_as_sfc()

dem_tiles <- get_tiles(study_area,
                 services=c("elevation"),
                 resolution= 1)

dem <- terra::rast(dem_tiles[["elevation"]][[1]])%>%
  aggregate(fact=10)


sites <- data.frame(lat=c(35.92,35.98,35.93,35.95),
                    long=c(-108.74,-108.68,-108.72,-108.7),
                    type = c("origin","destination","conducive site","conducive site"))%>%
  st_as_sf(coords=c("long","lat"),remove=F,crs=4326)
origin <- subset(sites,type=="origin")
destination <- subset(sites,type=="destination")


plot(dem)
points(sites)


###### STEP 1: BUILD ORIGINAL COST SURFACES ####################################

## Build time-based cost surface
time <-create_slope_cs(x=dem,cost_function="tobler",neighbours=8)

## Build energy-based cost surface
energy <- create_slope_cs(x=dem,cost_function="herzog",neighbours=8)

## Caculate two lcps
time_lcp <- create_lcp(x=time,origin=origin,destination=destination)
energy_lcp <- create_lcp(x=energy,origin=origin,destination=destination)


###### STEP 2: BUILD COMBINED COST SURFACE #####################################
# combine cost raster by standardizing values for each 
# reverting one or the other so that higher costs (in energy and time expense)
# are both represented by higher values 
# and multiplying rasters together

# standardize time cs
time_raster <- rasterise(time)
mm <- minmax(time_raster)    
time_raster_stand <- (time_raster - mm[1,]) / (mm[2,] - mm[1,])
time_raster_stand_cs <- create_cs(time_raster_stand,neighbours=4, dem=NULL,max_slope=NULL, exaggeration=F)

# ensure original time-based lcp and standardized time-based lcp are similar 
time_stand_lcp <- create_lcp(x=time_raster_stand_cs,origin=origin,destination=destination)

plot(dem)
points(sites)
lines(time_lcp,col="white",lwd=2.5)
lines(time_stand_lcp,col='dodgerblue3',lwd=2.5)


# standardize energy cs
energy_raster <- rasterise(energy)
mm <- minmax(energy_raster)
energy_raster_stand <- (energy_raster - mm[1,]) / (mm[2,] - mm[1,])

# ensure original time-based lcp and standardized time-based lcp are similar 
energy_raster_stand_cs <- create_cs(energy_raster_stand,neighbours=4, dem=NULL,max_slope=NULL, exaggeration=F)
energy_stand_lcp <- create_lcp(x=energy_raster_stand_cs,origin=origin,destination=destination)

plot(dem)
points(sites)
lines(energy_lcp,col="white",lwd=2.5)
lines(energy_stand_lcp,col='brown3',lwd=2.5)




## to average together we need to inverse values of energy based raster
energy_raster_stand <- (energy_raster - mm[1,]) / (mm[1,] - mm[2,])+1

# average together
tande <- (energy_raster_stand+time_raster_stand)/2

# rebuild into conductance surface
te <- create_cs(tande,neighbours=4, dem=NULL,max_slope=NULL, exaggeration=F)


# conduct another lcp to see how it changes 
te_lcp <- create_lcp(x=te,origin=origin,destination=destination)


plot(dem)
points(sites)
lines(time_stand_lcp,col='dodgerblue3',lwd=2.5)
lines(energy_stand_lcp,col='brown3',lwd=2.5)
lines(te_lcp,col="darkorchid3",lwd=2.5)






############### STEP 3: INCORPORATE GRAVITY ####################################
## First write some functions just to make things easier
##buffer function simplified
buffer_fun <- function(d){
  st_buffer(sites,d)
}

sites_cond <- subset(sites,type=="conducive site")

## run a loop that builds nested polygons for each feature type in site list
sitetypes <- as.data.frame(unique(sites_cond$type))

for(i in 1:nrow(sitetypes)){
  
  name_sitetype <- sitetypes[i,1]
  site_subset<-subset(sites,type==name_sitetype)
  
  buff <-st_buffer(site_subset,100)%>%
    st_union
  
  first_name <- paste(name_sitetype,"Buff_diff0")
  st_write(buff, paste("./output/github example/",first_name,".shp",sep=""))

  for (j in 1:9){
  
  # paste new buffer name
  buff_name <- paste("Buff",j,sep="")
  buff_name2 <- paste("Buff",j+1,sep="")
  
  #build buffer
  buff <-st_buffer(site_subset,20*(j*5))%>%
    st_union
  
  buff2 <- st_buffer(site_subset,20*((j+1)*5))%>%
    st_union
  
  # assign buffer new buffer name
  assign(buff_name,buff)
  assign(buff_name2,buff2)
  
  # new difference b/w buffer name
  buff_diff_name <- paste("Buff_diff",j,sep="")
  
  # calculate symettrical difference b/w buffers
  buff_diff <- st_sym_difference(buff,buff2)
  
  # assign buffer diff new buffer diff name
  #assign(buff_diff_name,buff_diff)
  second_name <- paste(name_sitetype,buff_diff_name)
  
  st_write(buff_diff, paste("./output/github example/",second_name,".shp",sep=""))
}
}



## READ THOSE BACK IN
for(i in 0:9){
  name <- paste("cond_site_",i,sep="")
  
  test <- sf::read_sf(paste("./output/github example/conducive site Buff_diff",i,".shp",sep=""))%>%
    st_transform(4326)
  
  assign(name,test)
}

plot(cond_site_1,col="red")


plot(dem)
points(sites)
plot(cond_site_0,add=T,col='#66ffcc')
plot(cond_site_1,add=T,col='#68fbce')
plot(cond_site_2,add=T,col='#6af6d0')
plot(cond_site_3,add=T,col='#6decd3')
plot(cond_site_4,add=T,col='#73d9d9')
plot(cond_site_5,add=T,col='#80b3e6')
plot(cond_site_6,add=T,col='#8d8df3')
plot(cond_site_7,add=T,col='#8d8df3')
plot(cond_site_8,add=T,col='#937af9')
plot(cond_site_9,add=T,col='#9670fc')

### incorporate areas as conduits 
update <- update_values(te,cond_site_0, FUN=function(j){replace(x=j,values=j/.9)})%>%
  update_values(.,cond_site_1,FUN=function(j){replace(x=j,values=j/.91)})%>%
  update_values(.,cond_site_2,FUN=function(j){replace(x=j,values=j/.92)})%>%
  update_values(.,cond_site_3,FUN=function(j){replace(x=j,values=j/.93)})%>%
  update_values(.,cond_site_4,FUN=function(j){replace(x=j,values=j/.94)})%>%
  update_values(.,cond_site_5,FUN=function(j){replace(x=j,values=j/.95)})%>%
  update_values(.,cond_site_6,FUN=function(j){replace(x=j,values=j/.96)})%>%
  update_values(.,cond_site_7,FUN=function(j){replace(x=j,values=j/.97)})%>%
  update_values(.,cond_site_8,FUN=function(j){replace(x=j,values=j/.98)})%>%
  update_values(.,cond_site_9,FUN=function(j){replace(x=j,values=j/.99)})


# conduct another lcp to see how it changes 
cond_lcp <- create_lcp(x=update,origin=origin,destination=destination)

plot(dem)
points(sites)
lines(te_lcp,col="darkorchid3",lwd=2.5)
lines(cond_lcp,col="white",lwd=2.5)










### RE-CREATE PLOTS FROM 
master <- read.csv("https://raw.githubusercontent.com/sfield2/Conductance-Centric-LCA/main/DATA/master_results.csv")

master_sub <- master%>%
  mutate(Average_Conductance=((Conductance.Value.Great.House+Conductance.Value.Earthen.Mound+
                                 Conductance.Value.Great.Kiva+Conductance.Value.Herradura)/4))

## things were combined to share difference between two
master_sub1 <- master_sub[,c(1,3,5,10:14)]%>%
  gather("Buffer","Buffer_Value",-Running.Test.Number,-Conductance.Value.Great.House,
         -Conductance.Value.Great.Kiva,-Conductance.Value.Earthen.Mound,
         -Conductance.Value.Herradura,-Average_Conductance)

master_sub2 <- master_sub[,c(1,4,6,10:14)]%>%
  gather("PDI","PDI_Value",-Running.Test.Number,-Conductance.Value.Great.House,
         -Conductance.Value.Great.Kiva,-Conductance.Value.Earthen.Mound,
         -Conductance.Value.Herradura,-Average_Conductance)

master_all<-merge(master_sub2,master_sub1)%>%
  mutate_at(c("Buffer_Value"), ~ (scale(.)%>% as.vector))%>%
  mutate_at(c("PDI_Value"), ~ (scale(.)%>% as.vector))%>%
  mutate(Average_Conductance_norm = Average_Conductance)%>%
  mutate_at(c("Average_Conductance_norm"), ~ (scale(.)%>%as.vector))%>%
  mutate(PDI_Value_rev = (PDI_Value * -1) )


## cond master
master_sub_conductance <- master_all%>%
  subset(.,PDI == "Conductance.Centric.PDI")%>%
  subset(.,Buffer == "Conductance.Centric.Buffer")%>%
  mutate(AC_Factor = ntile(Average_Conductance_norm,n=4))

master_cond <- master%>%
  subset(.,select=c("Conductance.Value.Great.House","Conductance.Value.Great.Kiva","Conductance.Value.Earthen.Mound",
                    "Conductance.Value.Herradura"))%>% 
  gather("Buffer_Type","Buffer_Value")

master_cond$Buffer_Type <-factor(master_cond$Buffer_Type,levels = c("Conductance.Value.Herradura", "Conductance.Value.Earthen.Mound",
                                                                    "Conductance.Value.Great.Kiva","Conductance.Value.Great.House"))

master_cond <- master_cond%>%
  mutate(Buffer_Value_2 = Buffer_Value *2)%>%
  mutate(Buffer_Value_3 = Buffer_Value *3)%>%
  mutate(Buffer_Value_4 = Buffer_Value *4)%>%
  mutate(Buffer_Value_5 = Buffer_Value *5)%>%
  mutate(Buffer_Value_6 = Buffer_Value *6)%>%
  mutate(Buffer_Value_7 = Buffer_Value *7)%>%
  mutate(Buffer_Value_8 = Buffer_Value *8)%>%
  mutate(Buffer_Value_9 = Buffer_Value *9)%>%
  mutate(Buffer_Value_10 = Buffer_Value *10)




# Figure 5
ggplot()+
  geom_jitter(data=master_cond,aes(Buffer_Value,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_2,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_3,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_4,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_5,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_6,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_7,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_8,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_9,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  geom_jitter(data=master_cond,aes(Buffer_Value_10,as.factor(Buffer_Type),color=as.factor(Buffer_Type)),size=5,alpha=0.05)+
  
  stat_summary(data=master_cond,aes(Buffer_Value_5,as.factor(Buffer_Type)),fun="mean",geom="crossbar",color="white")+
  
  scale_color_manual(values=c("#3A5B9466","#3680A466","#3DA5AB66","#63C8B199"))+
  scale_x_continuous(limit=c(0,0.06),breaks=seq(0,0.06,by=0.02),labels = c("0","2%","4%","6%"))+
  scale_y_discrete(labels=c("Herradura","Earthen Mound","Great Kiva","Great House"))+
  labs(title="",
       x="Increase in Conductance",
       y="")+
  theme(plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=24, face="bold"),
        axis.text = element_text(size=24),
        panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(.5, "cm"),
        axis.ticks.y=element_blank(),
        legend.position= "none")



### FIGURE 7
master_all$PDI <-factor(master_all$PDI,levels = c("Cost.Centric.PDI","Conductance.Centric.PDI"))
master_all$Buffer <-factor(master_all$Buffer,levels = c("Cost.Centric.Buffer","Conductance.Centric.Buffer"))


p1<-ggplot(data=master_all,aes(PDI_Value_rev,PDI,fill=PDI))+
  geom_density_ridges(quantile_lines=TRUE, quantiles = 2,col="white")+
  scale_x_continuous(limit=c(-1.5,1.5),breaks=seq(-1.4,1.4,by=2.8),labels=c("Less Similar","More Similar"))+
  scale_fill_manual(values=c("#BE3A7699","#3DA5AB99"))+
  
  annotate("text", x = -1.25, y = 2.3, label = "Conductance Centric",col="#3680a4",size=8,fontface="bold")+
  annotate("text", x = -1.4, y = 1.3, label = "Cost Centric",col="#882881",size=8,fontface="bold")+

  labs(title="Normalized PDI Index",
       x="",
       y="No. of Trials")+

  theme(plot.title = element_text(size=26,face="bold"),
        axis.title= element_blank(),
        axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(.5, "cm"),
        axis.ticks.y=element_blank(),
        legend.position= "none")


p2<-ggplot(data=master_all,aes(Buffer_Value,Buffer,fill=Buffer))+
  geom_density_ridges(quantile_lines=TRUE, quantiles = 2,col="white",scale=3)+
  scale_x_continuous(limit=c(-2,4),breaks=seq(-1.8,3.8,by=5.6),labels=c("Less Similar","More Similar"))+
  scale_fill_manual(values=c("#BE3A7699","#3DA5AB99"))+
  labs(title="Normalized Buffer Index",
       x="Similarity With Known Route",
       y="No. of Trials")+
  
  theme(plot.title = element_text(size=26,face="bold"),
        axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=24),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(.5, "cm"),
        axis.ticks.y=element_blank(),
        legend.position= "none")


ggarrange(p1,p2, ncol = 1, nrow = 2)





# Figure 8
p3<-ggplot()+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*2,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*3,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*4,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*5,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*6,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*7,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*8,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*9,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.House*10,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*2,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*3,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*4,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*5,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*6,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*7,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*8,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*9,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Great.Kiva*10,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*2,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*3,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*4,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*5,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*6,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*7,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*8,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*9,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Earthen.Mound*10,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*2,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*3,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*4,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*5,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*6,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*7,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*8,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*9,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+
  geom_jitter(data=master_sub_conductance,aes(Conductance.Value.Herradura*10,as.factor(AC_Factor),color=as.factor(AC_Factor)),alpha=0.1)+ 
  
  scale_x_continuous(limit=c(0,0.045),breaks=seq(0,0.04,by=0.02),labels = c("0","2%","4%"))+
  scale_y_discrete(limits=1:4,labels=c("Lowest
Quartile","","","Highest
Quartile"))+
  scale_color_manual(values=c("#3A5B9466","#3680A466","#3DA5AB66","#63C8B199"))+
  labs(title="Conductance Values",
       x="Increase in Conductance",
       y="Average Conductance")+
  theme(plot.title = element_text(size=26,face="bold"),
        axis.title = element_text(size=24, face="bold"),
        axis.text = element_text(size=24),
        panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(.5, "cm"),
        axis.ticks.y=element_blank(),
        legend.position= "none")


p4<-ggplot(data=master_sub_conductance,aes(PDI_Value_rev,as.factor(AC_Factor),fill=as.factor(AC_Factor)))+
  geom_density_ridges(quantile_lines=TRUE, quantiles = 2,col="white")+
  scale_x_continuous(limit=c(-0.6,1.2),breaks=seq(-0.5,1.1,by=1.6),labels=c("Less Similar","More Similar"))+
  scale_y_discrete(limits=1:5)+
  # scale_y_discrete(limits=1:6,breaks=1:5,labels=c("low","medlow","med","medhigh","high"))+
  scale_fill_manual(values=c("#3A5B9466","#3680A466","#3DA5AB66","#63C8B199"))+
  labs(title="Normalized PDI Index",
       x="Similarity with Known Route",
       y="")+
  theme(plot.title = element_text(size=26,face="bold"),
        axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=24),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(.5, "cm"),
        axis.ticks.y=element_blank(),
        legend.position= "none")

p5<-ggplot(data=master_sub_conductance,aes(Buffer_Value,as.factor(AC_Factor),fill=as.factor(AC_Factor)))+
  geom_density_ridges(quantile_lines=TRUE, quantiles = 2,col="white")+
  scale_x_continuous(limit=c(-2.5,4),breaks=seq(-2.14,3.64,by=5.78),labels=c("Less Similar","More Similar"))+
  scale_y_discrete(limits=1:5)+
  scale_fill_manual(values=c("#3A5B9466","#3680A466","#3DA5AB66","#63C8B199"))+
  labs(title="Normalized Buffer Index",
       x="Similarity with Known Route",
       y="")+
  theme(legend.position= "none")+
  theme(plot.title = element_text(size=26,face="bold"),
        axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=24),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(.5, "cm"),
        axis.ticks.y=element_blank(),
        legend.position= "none")



ggarrange(p3,p4,p5, ncol = 3, nrow = 1,widths = c(0.8,1.3,1.3))

















