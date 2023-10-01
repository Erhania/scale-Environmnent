LocalDir<-setwd(...)
#1,random forest
#running random forest (RF)
library(randomForest)
N0<-c("family","genus","species") #running RF for family,genus and species taxa respectively
for (ii in 1:length(N0)) {
  l<-list.files(file.path(LocalDir,N0[ii]))
  nnl<-strsplit(l,split = ".csv",fixed = T);nnl<-unlist(lapply(nnl,head,1))
  for (i in 1:length(l)){
    mdat<-read.csv(file.path(LocalDir.N0[ii],l[i]))
    mdat<-mdat[,c("bio_pc1","bio_pc2","bio_pc3","soil_pc1","soil_pc2","soil_pc3","soil_pc4",
                  "soil_pc5","GHI","elevation","slope","aspect","presence")]
    n<-unique(mdat$presence);n<-subset(n,!n==0) #get taxa names in a plot at specific spatial and taxonomic scale
    n1<-data.frame()
    for(j in 1:length(n)){
      mda<-subset(mdat,mdat$presence%in%c(n[j],"0")) #0represent absence from background data
      #if the sample size of data set are too large, kmeans can be used for shorten the operation time
      model<-randomForest(x=mda[,-13], y=as.factor(as.character(mda[,13])), ntree=500,importance=T,proximity=T)
      im<-importance(model);im<-as.data.frame(im)
      x1<-im$MeanDecreaseGini
      x1<-t(x1);n1<-rbind(n1,x1)
      save(file=paste("RFresults/",N0[ii],"/",nnl[i],"_",j,".rdata",sep=""),model) #save models for every taxon respectively
    }
    names(n1)<-c("bio_pc1","bio_pc2","bio_pc3","soil_pc1","soil_pc2","soil_pc3","soil_pc4","soil_pc5",
                 "GHI","elevation","slope","aspect")
    rownames(n1)<-n
    write.csv(n,file = paste("RFresults/",N0[ii],"/G_",l[i],sep=""),row.names = T) #save MDG result seperately for every data set
  }
}

#synthesize MDG value from all data sets
randomplot500<-read.csv("E:\\erhan\\scalerf\\500select.csv")
A<-data.frame()
for (ii in 1:length(N0)) {
  l<-list.files(file.path(LocalDir,"RFresults",,N0[ii]))
  n<-strsplit(l,split = "A_");n<-unlist(lapply(n,tail,1));n<-strsplit(n,split = "_");nspatial<-unlist(lapply(n,head,1))
  tmp<-unlist(lapply(n,tail,1));tmp<-strsplit(tmp,split = ".rds");np<-unlist(lapply(tmp,head,1))
  for (i in 1:length(nn)) {
    mdat<-read.csv(file.path(LocalDir,N0[ii],l[i]))
    mdat$space<-nspatial[i];mdat$plot<-np[i];mdat$taxa<-N0[ii]
    A<-rbind(A,mdat)
  }
}
write.csv(A,"MDG.csv",row.names = F)


#standardize MDG
a<-read.csv("MDG.csv")
y<-apply(a[,1:12], 1, min);a[,1:12]<-a[,1:12]-y
y<-rowSums(a[,1:12]);a[,1:12]<-a[,1:12]/y
a$taxaname<-rownames(a)
write.csv(a,"standardMDG.csv",row.names = F)



#2, mixed effect model
#sum of MDG of variables in a group (climate/soil/topography)
a<-read.csv("standardMDG.csv")
#climate
d<-subset(a,a$variable%in%c("bio_pc1","bio_pc2","bio_pc3","GHI"))
n<-unique(d[,c("plot","space","taxa")]);dd<-vector()
for (i in 1:nrow(n)) {
  d1<-subset(d,d$plot==n$plot[i]&d$space==n$space[i]&d$taxa==n$taxa[i])
  dd[i]<-sum(d1$value)
}
x<-cbind(n,dd);names(x)[4]<-"value"
write.csv(x,"clim_rf.csv",row.names = F)
#soil
d<-subset(a,a$variable%in%c("soil_pc1","soil_pc2","soil_pc3","soil_pc4","soil_pc5"))
n<-unique(d[,c("plot","space","taxa")]);dd<-vector()
for (i in 1:nrow(n)) {
  d1<-subset(d,d$plot==n$plot[i]&d$space==n$space[i]&d$taxa==n$taxa[i])
  dd[i]<-sum(d1$value)
}
x<-cbind(n,dd);names(x)[4]<-"value"
write.csv(x,"soil_rf.csv",row.names = F)
#topography
d<-subset(a,a$variable%in%c("elevation","slope","aspect"))
n<-unique(d[,c("plot","space","taxa")]);dd<-vector()
for (i in 1:nrow(n)) {
  d1<-subset(d,d$plot==n$plot[i]&d$space==n$space[i]&d$taxa==n$taxa[i])
  dd[i]<-sum(d1$value)
}
x<-cbind(n,dd);names(x)[4]<-"value"
write.csv(x,"topo_rf.csv",row.names = F)

#running mixed effect model
library(lme4);library(reshape2);library(MuMIn)
d<-read.csv("standMDG_climate")
d$plot<-as.factor(d$plot)
d$space<-as.factor(d$space)
d$point_taxa = paste(d$plot, d$taxa, sep="_")
d$point_scale = paste(d$plot, d$space, sep="_")
d$point_ts = paste(d$point_taxa, d$space, sep="_")

m1<-lmer(value ~taxa
         + space
         + taxa:space
         +(1|plot)
         +(1|point_scale)
         +(1|point_taxa),data=d,
         control=lmerControl(check.nobs.vs.nlev = "ignore",
                             check.nobs.vs.rankZ = "ignore",
                             check.nobs.vs.nRE="ignore"))
#check model assumptions
res <- resid(m1)
plot(fitted(m1), res);abline(0,0)
qqnorm(res);qqline(res)
plot(density(res))
#anova and R2 for the mixed effect models
anova(m1)
library(MuMIn)
r.squaredGLMM(m1)

##repeat this analyses for climate, soil and topography, as well as for 12 variables respectively.
