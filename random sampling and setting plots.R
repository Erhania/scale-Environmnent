setwd() #set the direction

#randomly sample 500 location across the globe
library(raster);library(dismo) #use a global data as background raster
m1<-raster("wc2.0_bio_30s_01.tif") 
randomplot<-randomPoints(m1,500) #to get longitude and latitude

#plotting biomes involved by these locations across the globe
library(plotbiomes);library(ggplot2)
#Whittaker biomes (defined by annual temperature and precipitation)
m1<-raster("wc2.0_bio_30s_01.tif") #annual mean temperature (MAT) from worldclim
m2<-raster("wc2.0_bio_30s_12.tif") #annual mean precipitation (MAP)from worldclim
m<-stack(m1,m2)
x<-extract(m,randomplot) #extract MAT and MAP using coordinates of 500 locations 
x<-as.data.frame(x);names(x)<-c("MAT","MAP");x$MAP<-x$MAP/10  #convert to same unit in Whittaker biomes figure
whittaker_base_plot() +
  theme(legend.position = c(0.2, 0.75),
        panel.background = element_blank(),
        panel.grid.major = element_line(gray(0.7)),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size=12),
        axis.title = element_text(size = 14))+
  geom_point(data = randomplot,aes(x=MAT,y=MAP),
             size=2,shape=21,color="gray95",fill="black",stroke=1,alpha=0.5)

#WWFb biomes:
biome<-raster("WWFbiome.tif")
x<-extract(biome,randomplot)
x<-cbind(randomplot,x)
table(x[,3]) #frequence of biomes in which 500 locations occur

####
#NOTE: the sampling step could be repeated for many times until that all the locations seem to cover most of biomes in earth 
####

#set plots with different spaital extents from 1° × 1° to 10° × 10°
#save data sets in the "plotsdata" file (provided as supplementary data in our study)
randomplot
l<-list.files(file.path(getwd(),"plantdata")) #plant distribution data files downloaded from GBIF and manipulated in R (to reduce sample bias) before modeling
for (i in 1:length(nfam)) {
  mdat<-read.csv(file.path(getwd(),"plantdata",l[i])) #read plant data for a family at one time (including family, genus and species name, geographical coordinates and corresponding environmental variables)
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+0.5)&mdat$lon>(randomplot[t,1]-0.5))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+0.5)&mda$lat>(randomplot[t,2]-0.5))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("//plotsdata//",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+1)&mdat$lon>(randomplot[t,1]-1))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+1)&mda$lat>(randomplot[t,2]-1))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("//plotsdata//",t,"_",i,".csv",sep=""),row.names = F)
    } 
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+1.5)&mdat$lon>(randomplot[t,1]-1.5))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+1.5)&mda$lat>(randomplot[t,2]-1.5))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("//plotsdata//",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+2)&mdat$lon>(randomplot[t,1]-2))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+2)&mda$lat>(randomplot[t,2]-2))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+2.5)&mdat$lon>(randomplot[t,1]-2.5))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+2.5)&mda$lat>(randomplot[t,2]-2.5))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+3)&mdat$lon>(randomplot[t,1]-3))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+3)&mda$lat>(randomplot[t,2]-3))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+3.5)&mdat$lon>(randomplot[t,1]-3.5))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+3.5)&mda$lat>(randomplot[t,2]-3.5))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+4)&mdat$lon>(randomplot[t,1]-4))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+4)&mda$lat>(randomplot[t,2]-4))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+4.5)&mdat$lon>(randomplot[t,1]-4.5))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+4.5)&mda$lat>(randomplot[t,2]-4.5))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
  for (t in 1:nrow(randomplot)) {
    mda<-subset(mdat,mdat$lon<(randomplot[t,1]+5)&mdat$lon>(randomplot[t,1]-5))
    mda<-subset(mda,mda$lat<(randomplot[t,2]+5)&mda$lat>(randomplot[t,2]-5))
    if(nrow(mda)>1){
      write.csv(mda,file = paste("",t,"_",i,".csv",sep=""),row.names = F)
    }
  }
}







