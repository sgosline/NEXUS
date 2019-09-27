#test sex differences in NF

require(synapser)
synLogin()
require(tidyverse)

tab<-synapser::synTableQuery("select * from syn20710536")$asDataFrame()

tab<-subset(tab,diagnosis=='Neurofibromatosis 1')

tab$tumorType[which(tab$tumorType%in%c("Ependymoma","Ganglioglioma"))]<-'Other'

tab$tumorType[which(tab$tumorType=="Malignant peripheral nerve sheath tumor")]<-"Malignant Peripheral Nerve Sheath Tumor"


##which samples are significantly different
res<-tab%>%spread(key=sex,value=score)%>%
  group_by(method,cell_type)%>%
    mutate(pval=t.test(female,male)$p.value)%>%
  select(method,cell_type,pval)%>%distinct()%>%
  group_by(method)%>%
  mutate(correctedP=p.adjust(pval))


sigs<-subset(res,pval<0.05)

for(meth in unique(tab$method)){
  tab.p<-subset(tab,method==meth)%>%subset(cell_type%in%(sigs$cell_type))

  p<-ggboxplot(tab.p,x='sex',y='score',facet.by='cell_type',color='sex',palette='jco')+stat_compare_means(method='t.test')
  print(p)
  ggsave(paste0(meth,'diffs.png'))
}

##now what do we see on a tissue level? 
res.c<-tab%>%subset(tumorType!='Other')%>%
  spread(key=sex,value=score)%>%
  group_by(method,tumorType,cell_type)%>%
  mutate(pval=t.test(female,male)$p.value)%>%
  select(method,cell_type,pval,tumorType)%>%distinct()%>%
  group_by(method)%>%
  mutate(correctedP=p.adjust(pval))

sigs<-subset(res.c,pval<0.05)

for(meth in unique(sigs$method)){
  sigs.t=subset(sigs,method==meth)
#  for(tu in unique(sigs$tumorType)){
  tab.t=subset(tab,tumorType%in%sigs.t$tumorType)%>%subset(method==meth)%>%subset(cell_type%in%sigs.t$cell_type)
    
  #  tab.p<-subset(tab.t,method==meth)%>%subset(cell_type%in%(sigs.t$cell_type))
 p<-ggplot(tab.t,palette='jco')+geom_boxplot(aes(x=cell_type,fill=sex,y=score))+facet_grid(.~tumorType)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
 # p<-ggboxplot(subset(tab.t,cell_type=='Neutrophil'),x='tumorType',y='score',color='sex',palette='jco')+stat_compare_means(method='t.test')+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  ggsave(paste0(meth,'TTdiffs.png'))
  
}


for(meth in unique(sigs$method)){
  sigs.t=subset(sigs,method==meth)
  for(tu in sigs.t$tumorType){
  #  for(tu in unique(sigs$tumorType)){
  tab.t=subset(tab,tumorType==tu)%>%
      subset(method==meth)%>%
        subset(cell_type%in%sigs.t$cell_type)
  #  tab.p<-subset(tab.t,method==meth)%>%subset(cell_type%in%(sigs.t$cell_type))
  p<-ggplot(tab.t,palette='jco')+geom_boxplot(aes(x=cell_type,fill=sex,y=score))+facet_grid(.~tumorType)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # p<-ggboxplot(subset(tab.t,cell_type=='Neutrophil'),x='tumorType',y='score',color='sex',palette='jco')+stat_compare_means(method='t.test')+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  ggsave(paste0(meth,tu,'diffs.png'))
  }
}
  

