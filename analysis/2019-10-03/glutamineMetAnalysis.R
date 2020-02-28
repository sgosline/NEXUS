library(synapser)
synLogin()

require(tidyverse)

path.data<-read.csv(synGet("syn20684776")$path,check.names=F)

glut.enzymes=c("PPAR","PFAS","GMPS","CPS-II","CTPS","NADSYN","GLS1","GLUL","ASNS","GFAT")

wgs.vars=synTableQuery("SELECT * FROM syn20551862 where \"Hugo_Symbol\" in ('TP53','EED','SUZ12','CDKN2A','NF1')")$asDataFrame()
exome.vars=synTableQuery("SELECT * FROM syn20554939 where \"Hugo_Symbol\" in ('TP53','EED','SUZ12','CDKN2A','NF1')")$asDataFrame()

gene.ex.data<-synTableQuery("select * from syn20449214")$asDataFrame()

##bracket various pathways on mutations
all.vars<-rbind(select(wgs.vars,'Hugo_Symbol','Protein_position','sex','specimenID','IMPACT','species','disease','tumorType','studyId','study'),select(exome.vars,'Hugo_Symbol','Protein_position','sex','specimenID','IMPACT','species','disease','tumorType','studyId','study'))


colnames(path.data)[1]='pathName'
path.tidied=path.data%>%gather(key='specimenID',value='score',-pathName)

##get specimens
specs.in.studies=subset(gene.ex.data,!is.na(study))%>%
    subset(study%in%all.vars$study)%>%
    select(specimenID)%>%distinct()

var.data=all.vars%>%left_join(path.tidied,by='specimenID')
genes=all.vars%>%
    mutate(mutated=ifelse(is.na(IMPACT),'WT','Mutated'))%>%
    spread(key=Hugo_Symbol,value='mutated',fill='WT')%>%
    select(specimenID,sex,disease,tumorType,CDKN2A,EED,NF1,SUZ12,TP53)

path.with.genes<-path.tidied%>%
    subset(specimenID%in%specs.in.studies$specimenID)%>%
  left_join(genes,by='specimenID')


muts=c('TP53','EED','SUZ12','CDKN2A','NF1')
library(ggpubr)
glut.path=path.with.genes[grep('GLUTAM',path.with.genes$pathName),]
pvals=data.frame()
for(mgene in muts){
  gp=select(glut.path,pathName,score,gene=mgene)
  gp$gene[is.na(gp$gene)]<-'WT'
  if(length(unique(gp$gene))==1)
    next
  cm=compare_means(score~gene,gp,group.by="pathName",method='t.test')
p<-ggplot(gp)+
  geom_boxplot(aes(x=gene,y=score,fill=pathName))+
  scale_fill_brewer(palette="Dark2")
  ggtitle(paste(mgene,'status in expression pathways'))
  ggsave(paste0(mgene,'statusInGlutaminPathways.png'))
  pvals=rbind(pvals,cbind(cm,gene=mgene))
}

glut.path=path.with.genes[grep('PURIN',path.with.genes$pathName),]

for(mgene in muts){
  gp=select(glut.path,pathName,score,gene=mgene)
  gp$gene[is.na(gp$gene)]<-'WT'
  if(length(unique(gp$gene))==1)
    next
  cm=compare_means(score~gene,gp,group.by="pathName",method='t.test')
  
  p<-ggplot(gp)+
    geom_boxplot(aes(x=pathName,y=score,fill=gene),position='dodge')+
    scale_fill_brewer(palette="Dark2")+
    theme(text=element_text(size=8),axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste(mgene,'status in expression pathways'))
  ggsave(paste0(mgene,'statusInPurinePathways.png'))
  pvals=rbind(pvals,cbind(cm,gene=mgene))
  
}
glut.path=path.with.genes[grep('PYRIM',path.with.genes$pathName),]

for(mgene in muts){
  gp=select(glut.path,pathName,score,gene=mgene)
  gp$gene[is.na(gp$gene)]<-'WT'
  if(length(unique(gp$gene))==1)
    next
  cm=compare_means(score~gene,gp,group.by="pathName",method='t.test')
  
  p<-ggplot(gp)+
    geom_boxplot(aes(x=pathName,y=score,fill=gene),position='dodge')+
    scale_fill_brewer(palette="Dark2")+
    theme(text=element_text(size=8),axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste(mgene,'status in expression pathways'))
  ggsave(paste0(mgene,'statusInPyrimidinePathways.png'))
  pvals=rbind(pvals,cbind(cm,gene=mgene))
  
}
