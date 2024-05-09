packages <-c('sf','terra','ggplot2','leastcostpath','lwgeom','dplyr','tidyr',
             'ggbeeswarm','hrbrthemes','ggridges','ggpubr','terrainr','tidyterra')
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
  aggregate(fact=50)


sites <- data.frame(lat=c(35.92,35.98,35.93,35.95),
                    long=c(-108.716,-108.68,-108.72,-108.7),
                    type = c("origin","destination","conducive site","conducive site"))%>%
  st_as_sf(coords=c("long","lat"),remove=F,crs=4326)

 
origin <- subset(sites,type=="origin")
destination <- subset(sites,type=="destination")

road <- data.frame(lat=c(35.92,35.98),
                   long=c(-108.716,-108.68),
                   type=c("road","road"))%>%
  st_as_sf(coords=c("long","lat"),remove=F,crs=4326)%>%
  group_by(type)%>%
  dplyr::summarize(do_union=F)%>%
  st_cast("LINESTRING")

ggplot()+
  geom_spatraster(data=dem)+
  geom_sf(data=sites,aes(col=type))+
  geom_sf(data=road,col="white")


############### STEP 1: IMPORT FUNCTIONS #######################################
### Function for creating buffers around sites
# @ param conducive_sites: sf object of site locations (class "sf"  "data.frame")
# @ param buffer_distance: distance of each buffer in meters (class "integer")
# @ param buffer_number: number of consecutive buffers (class "integer")
# @ param output_location: location of output folder (class "char")

build_buffers <- function(conducive_sites,buffer_distance,buffer_number,output_location){
  
  sitetypes <- as.data.frame(unique(conducive_sites$type))
  
  for(i in 1:nrow(sitetypes)){
    name_sitetype <- sitetypes[i,1]
    site_subset<-subset(sites,type==name_sitetype)
    
    for(j in 1:buffer_number){
      
      if (j == 1){
        first_buffer <-st_buffer(site_subset,buffer_distance)%>%
          st_union
        first_name <- paste(name_sitetype,"buffer 1")
        st_write(first_buffer, paste(output_location,first_name,".shp",sep=""),delete_dsn=T)
      } else {
        prior_buffer <- st_buffer(site_subset,(buffer_distance + (buffer_distance * (j-2))))%>%
          st_union
        next_buffer <- st_buffer(site_subset,(buffer_distance + (buffer_distance * (j-1))))%>%
          st_union
        buffer_difference <- st_sym_difference(prior_buffer,next_buffer)
        next_name <- paste(name_sitetype,"buffer",j,sep=" ")
        st_write(buffer_difference, paste(output_location,next_name,".shp",sep=""),delete_dsn=T)
      }
    }
  }
}


### Function for combining two or more cost surfaces
# @ param elev: dem surface (class "spatraster")
# @ param cost_function_1: name of cost algorithm used in creation of first surface (class "integer")
# @ param cost_function_2: name of cost algorithm used in creation of second surface (class "integer")
# @ param neigh: number of neighbours (class "integer")

# code suggestions courtesy of anonymous reviewer
combine_slope_cs <- function(elev,cost_function_1,cost_function_2,neigh){
  first_surface <- create_slope_cs(x=elev,cost_function = cost_function_1,neighbours=neigh)
  second_surface <- create_slope_cs(x=elev,cost_function = cost_function_2,neighbours=neigh)
  
  mm <- quantile(first_surface$conductanceMatrix@x, c(0, 1))
  first_surface$conductanceMatrix@x <- (first_surface$conductanceMatrix@x - mm[1]) / (mm[2] - mm[1])
  
  mm <- quantile(second_surface$conductanceMatrix@x, c(0, 1))
  second_surface$conductanceMatrix@x <- (second_surface$conductanceMatrix@x - mm[1]) / (mm[2] - mm[1])
  
  combined_surface <- first_surface
  combined_surface$conductanceMatrix@x <- (first_surface$conductanceMatrix@x + second_surface$conductanceMatrix@x) / 2
  
  return(combined_surface)
}


############### STEP 2: BUILD COMBINED COST SURFACE ############################
combined_surface <- combine_slope_cs(elev=dem,cost_function_1="tobler",cost_function_2="herzog",neigh = 8)
logistic_lcp <- create_lcp(x=combined_surface,origin=origin,destination=destination)

plot(combined_surface)


############### STEP 3: INCORPORATE GRAVITY ####################################
conducive_sites <- subset(sites,type=="conducive site")

# build buffers and read them back in
build_buffers(conducive_sites,500,5,"../output/github example/")

for(i in 1:5){
  name <- paste("cond_site_buff_",i,sep="")
  
  test <- sf::read_sf(paste("../output/github example/conducive site buffer ",i,".shp",sep=""))%>%
    st_transform(4326)
  
  assign(name,test)
}

# incorporate areas as conduits 
update <- update_values(combined_surface,cond_site_buff_1, FUN=function(j){replace(x=j,values=j/.85)})%>%
  update_values(.,cond_site_buff_2,FUN=function(j){replace(x=j,values=j/.87)})%>%
  update_values(.,cond_site_buff_3,FUN=function(j){replace(x=j,values=j/.89)})%>%
  update_values(.,cond_site_buff_4,FUN=function(j){replace(x=j,values=j/.91)})%>%
  update_values(.,cond_site_buff_5,FUN=function(j){replace(x=j,values=j/.93)})

plot(update)

############### STEP 4: Create Updated LCP #####################################
social_lcp <- create_lcp(x=update,origin=origin,destination=destination)

ggplot()+
  geom_spatraster(data=dem)+
  geom_sf(data=sites,aes(col=type))+
  geom_sf(data=road,col="white")+
  geom_sf(data=logistic_lcp,col="darkslategray2")+
  geom_sf(data=social_lcp,col="darkolivegreen2")


############### STEP 5: Example of a simulation ################################
iterations_master <- 1 # specify number of new dems to create
iterations_per <- 25   # specify number of stochastic lcps to produce for each dem

for(j in 1:iterations_master){
  ## build df for tabular outcomes
  results <- as.data.frame(matrix(NA,(iterations_per),8))
  colnames(results) <- c("Test Number",
                         "Cost Centric Buffer",
                         "Cost Centric PDI",
                         "Conductance Centric Buffer",
                         "Conductance Centric PDI",
                         "DEM error",
                         "Stochasticticity error",
                         "Conductance Value Site Type 1")
  
  # ensure TRUE to reduce computational costs
  sf_use_s2(TRUE)
  
  # step 5.1: incorporate RSME into dem (Lewis 2020)
  n <- abs(rnorm(1, mean=0,sd=0.82))
  r_err <- add_dem_error(x=dem,rmse=n,type="u")
  
  # step 5.2: build cost surfaces
  combined_surface <- combine_slope_cs(elev=r_err,cost_function_1="tobler",cost_function_2="herzog",neigh = 8)
  
  for (i in 1:iterations_per){
    # step 5.3: incorporate global stochasticity into each cs (Pinto and Keitt 2009)
    m <- abs(rnorm(1, mean=0.15,sd= 0.05)) # mean and sd changed in this example to produce different outputs
                                            # case study in Field et al. (2024) used: bs(rnorm(1, mean=0.005,sd= 0.025))
    combined_surface <- add_global_stochasticity(combined_surface,percent_quantile = m)
    
    # step 5.3: create lcp and assign to environment with name based on iterations_per
    lcp <- create_lcp(x=combined_surface,origin=origin,destination=destination)
    lcp_name <- paste("lcp_",i,sep="")
    assign(lcp_name,lcp)
    
    # ensure this is FALSE to save some computational cost
    sf_use_s2(FALSE)
    
    # step 5.4: run validation results on lcp and known route
    buff_val <- buffer_validation(lcp,road,dist=c(0.001))
    pdi_val <- as.data.frame(PDI_validation(lcp,road))
    
    results[i,2]<- buff_val$similarity
    results[i,3]<- pdi_val$pdi
    
    # step 5.5: create a random conductance value and use to increase conductance around sites within buffer regions 
    cv <- abs(rnorm(1, mean=0.02,sd= 0.05)) # mean and sd changed in this example to produce different outputs
                                              # case study in Field et al. (2024) used: bs(rnorm(1, mean=0.003,sd= 0.001)) 
    update_surface <- update_values(combined_surface,cond_site_buff_1, FUN=function(j){replace(x=j,values=j/(1-(cv*10)))})%>%
      update_values(.,cond_site_buff_2,FUN=function(j){replace(x=j,values=j/(1-(cv*9)))})%>%
      update_values(.,cond_site_buff_3,FUN=function(j){replace(x=j,values=j/(1-(cv*8)))})%>%
      update_values(.,cond_site_buff_4,FUN=function(j){replace(x=j,values=j/(1-(cv*7)))})%>%
      update_values(.,cond_site_buff_5,FUN=function(j){replace(x=j,values=j/(1-(cv*6)))})
    
    # step 5.6: create new lcp and assign to environment with name based on iterations_per
    lcp_update <- create_lcp(x=update_surface,origin=origin,destination=destination)
    lcp_update_name <- paste("lcp_update_",i,sep="")
    assign(lcp_update_name,lcp_update)
    
    # step 5.7: run validation results on update lcp and known route
    buff_val <- buffer_validation(lcp_update,road,dist=c(0.001))
    pdi_val <- as.data.frame(PDI_validation(lcp_update,road))
    
    # step 5.8: record other important parameters and metrics
    results[i,1]<- i    # iteration number 
    results[i,6] <- n   # RMSE error
    results[i,7] <- m   # Stochasticity
    results[i,8] <- cv  # conductance value 

  }
  
  # step 5.9: export lcps and tabular data
  logistic_lcp <- dplyr::bind_rows(list(lcp_1,lcp_2,lcp_3,lcp_4,lcp_5,lcp_6,lcp_7,lcp_8,lcp_9,lcp_10,
                                        lcp_11,lcp_12,lcp_13,lcp_14,lcp_15,lcp_16,lcp_17,lcp_18,lcp_19,lcp_20,
                                        lcp_21,lcp_22,lcp_23,lcp_24,lcp_25))
  log_lcp <- st_union(logistic_lcp)
  st_write(log_lcp,paste("../output/github example/logistic_paths_",j,".shp",sep=""),append=F)
  
  social_lcp <- dplyr::bind_rows(list(lcp_update_1,lcp_update_2,lcp_update_3,lcp_update_4,lcp_update_5,lcp_update_6,lcp_update_7,lcp_update_8,lcp_update_9,lcp_update_10,
                                       lcp_update_11,lcp_update_12,lcp_update_13,lcp_update_14,lcp_update_15,lcp_update_16,lcp_update_17,lcp_update_18,lcp_update_19,lcp_update_20,
                                       lcp_update_21,lcp_update_22,lcp_update_23,lcp_update_24,lcp_update_25))
  soc_lcp <- st_union(social_lcp)
  st_write(soc_lcp,paste("../output/github example/social_paths_",j,".shp",sep=""),append=F)
  
  write.csv(results,paste("../output/github example/tabular_results_",j,".csv",sep=""))
}


ggplot()+
  geom_spatraster(data=dem)+
  geom_sf(data=sites,aes(col=type))+
  geom_sf(data=road,col="white")+
  geom_sf(data=logistic_lcp,col="darkslategray2",lwd=2)+
  geom_sf(data=social_lcp,col="darkolivegreen2")










### RE-CREATE PLOTS FROM 
master <- read.csv("https://raw.githubusercontent.com/sfield2/Conductance-Centric-LCA/main/DATA/master_results.csv")

master_sub <- master%>%
  mutate(Average_Conductance=((Conductance.Value.Great.House+Conductance.Value.Earthen.Mound+
                                 Conductance.Value.Great.Kiva+Conductance.Value.Herradura)/4))

## things were combined to share difference between two
master_sub1 <- master_sub%>%
  select(c("Running.Test.Number","Cost.Centric.Buffer","Conductance.Centric.Buffer",
           "Conductance.Value.Great.House","Conductance.Value.Great.Kiva",
           "Conductance.Value.Earthen.Mound","Conductance.Value.Herradura","Average_Conductance"))%>%
  gather("Buffer","Buffer_Value",-Running.Test.Number,-Conductance.Value.Great.House,
         -Conductance.Value.Great.Kiva,-Conductance.Value.Earthen.Mound,
         -Conductance.Value.Herradura,-Average_Conductance)


master_sub2 <- master_sub%>%
  select(c("Running.Test.Number","Cost.Centric.PDI","Conductance.Centric.PDI",
           "Conductance.Value.Great.House","Conductance.Value.Great.Kiva",
           "Conductance.Value.Earthen.Mound","Conductance.Value.Herradura","Average_Conductance"))%>%
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
mean(master_cond$Buffer_Value,na.rm=T)

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

# REMOVE WHEN FULL VALUES 
master_cond <- subset(master_cond, Buffer_Value >0)
master_all <- subset(master_all,Conductance.Value.Great.House >0)



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
ggsave("../FIGURES/Figure 5_update.png",dpi=600,width = 12000, height = 7000,units=c("px"))




### FIGURE 6 INSET
master_all$PDI <-factor(master_all$PDI,levels = c("Cost.Centric.PDI","Conductance.Centric.PDI"))
master_all$Buffer <-factor(master_all$Buffer,levels = c("Cost.Centric.Buffer","Conductance.Centric.Buffer"))


p1<-ggplot(data=master_all,aes(PDI_Value_rev,PDI,fill=PDI))+
  geom_density_ridges(quantile_lines=TRUE, quantiles = 2,col="white")+
  scale_x_continuous(limit=c(-1.5,1.5),breaks=seq(-1.4,1.4,by=2.8),labels=c("Less Similar","More Similar"))+
  scale_fill_manual(values=c("#BE3A7699","#3DA5AB99"))+
  
  annotate("text", x = -1.25, y = 2.5, label = "Social Model",col="#3680a4",size=8,fontface="bold")+
  annotate("text", x = -1.25, y = 1.5, label = "Logistic Model",col="#882881",size=8,fontface="bold")+

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


plot1<-ggarrange(p1,p2, ncol = 1, nrow = 2)
plot1
ggsave("../FIGURES/Figure 6 inset_update.png",dpi=600,width = 9000, height = 5000,units=c("px"))




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



plot2<-ggarrange(p3,p4,p5, ncol = 3, nrow = 1,widths = c(0.8,1.3,1.3))

ggsave("../FIGURES/Figure 8_update.png",dpi=600,width = 18000, height = 6000,units=c("px"))
