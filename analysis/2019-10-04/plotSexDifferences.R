#test sex differences in NF

require(synapser)
synLogin()
require(tidyverse)

path.data<-read.csv(synGet("syn20684776")$path,check.names=F)
gene.ex.data<-synTableQuery("select distinct specimenID,tumorType,diagnosis,sex,isCellLine from syn20449214")$asDataFrame()


gene.ex.data$specimenID=sapply(gene.ex.data$specimenID,function(x) gsub('Plexiform ','p',gsub('Neurofibroma','NF',gsub('Malignant Peripheral Nerve Sheath Tumor','MPNST',x))))
gene.ex.data$sex=tolower(gene.ex.data$sex)
colnames(path.data)[1]='pathName'
path.tidied=path.data%>%gather(key='specimenID',value='score',-pathName)

tab<-path.tidied%>%left_join(gene.ex.data,by='specimenID')


tab<-subset(tab,diagnosis=='Neurofibromatosis 1')
#cl=which(tab$isCellLine=='TRUE')
#if(length(cl)>0)
#  tab=tab[-cl,]
#tab$tumorType[which(tab$tumorType%in%c("Ependymoma","Ganglioglioma"))]<-'Other'

#tab$tumorType[which(tab$tumorType=="Malignant peripheral nerve sheath tumor")]<-"Malignant Peripheral Nerve Sheath Tumor"


##now what do we see on a tissue level? 
res<-tab%>%
  spread(key=sex,value=score)%>%
  group_by(pathName)%>%
  mutate(pval=t.test(female,male)$p.value)%>%
  select(pathName,pval)%>%distinct()%>%
  ungroup()%>%
  mutate(correctedP=p.adjust(pval))

sigs.all<-subset(res,correctedP<0.05)
View(sigs.all)

##now what do we see on a tissue level? 
res.c<-tab%>%
  spread(key=sex,value=score)%>%
  group_by(pathName,tumorType)%>%
  mutate(pval=t.test(female,male)$p.value)%>%
  select(pathName,pval,tumorType)%>%distinct()%>%
  ungroup()%>%
  group_by(tumorType)%>%
  mutate(correctedP=p.adjust(pval))

sigs<-subset(res.c,correctedP<0.05)
View(sigs)
for(ct in unique(sigs$pathName)){
  sigs.t=subset(sigs,pathName==ct)
#  for(tu in unique(sigs$tumorType)){
  tab.t=subset(tab,pathName%in%sigs.t$pathName)%>%subset(pathName==ct)
    
  #  tab.p<-subset(tab.t,method==meth)%>%subset(cell_type%in%(sigs.t$cell_type))
 p<-ggplot(tab.t,palette='jco')+geom_boxplot(aes(x=pathName,fill=sex,y=score))+facet_grid(.~tumorType)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
 # p<-ggboxplot(subset(tab.t,cell_type=='Neutrophil'),x='tumorType',y='score',color='sex',palette='jco')+stat_compare_means(method='t.test')+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  ggsave(paste0(gsub(' ','',gsub('/','',ct)),'SigSexdiffs.png'))
  
}


for(meth in unique(sigs$pathName)){
  sigs.t=subset(sigs,pathName==meth)
  for(tu in sigs.t$tumorType){
  #  for(tu in unique(sigs$tumorType)){
  tab.t=subset(tab,tumorType==tu)%>%
      subset(pathName==meth)
    
  #  tab.p<-subset(tab.t,method==meth)%>%subset(cell_type%in%(sigs.t$cell_type))
  p<-ggplot(tab.t,palette='jco')+geom_boxplot(aes(x=tumorType,fill=sex,y=score))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_y_log10()+ggtitle(paste(meth,'scores'))
  # p<-ggboxplot(subset(tab.t,cell_type=='Neutrophil'),x='tumorType',y='score',color='sex',palette='jco')+stat_compare_means(method='t.test')+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  ggsave(paste0(meth,gsub(' ','',gsub('/','',tu)),'diffs.png'))
  }
}
  

