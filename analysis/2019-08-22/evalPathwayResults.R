#eval pathway results
require(synapser)
synLogin()

#get paths
gsva<-read.csv(synGet('syn20684776')$path,check.names=F,row.names=1)

#get samples
tab<-synTableQuery('select distinct specimenID,individualID,isCellLine,diagnosis,tumorType,transplantationType from syn20449214')$asDataFrame()

tab$sampleType<-apply(tab,1,function(x){
  if(!is.na(x[['transplantationType']]) && x[['transplantationType']]=='xenograft')
    return('xenograft')
  else if (!is.na(x[['isCellLine']])&&x[['isCellLine']]==TRUE)
    return('cell line')
  else
    return('tissue')})

tab$tumorType[is.na(tab$tumorType)]<-'None'
tab$diagnosis[is.na(tab$diagnosis)]<-'None'
tab$isCellLine[is.na(tab$isCellLine)]<-'FALSE'
tab$isCellLine<-toupper(tab$isCelline)

tab$specimenID<-sapply(tab$specimenID,function(x)
  gsub('Neurofibroma','NF',
    gsub("Plexiform ",'p',
      gsub('Cutaneous ',"c",
        gsub('Malignant Peripheral Nerve Sheath Tumor','MPNST',x)))))

#tidy GSVA pathways
require(tidyverse)
gsva$pathway=rownames(gsva)
gsva.tid<-gsva%>%
    gather(key=specimenID,value=score,-pathway)%>%
  left_join(tab,by='specimenID')

gsva.tis=subset(gsva.tid,sampleType=='tissue')

ggplot(gsva.tis[grep('COMPLEMENT',gsva.tis$pathway),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('complement.png')

ggplot(gsva.tis[grep('_ALK_',gsva.tis$pathway),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('alk.png')
ggplot(gsva.tis[intersect(grep('FIBROBLAST',gsva.tis$pathway),grep('ASSOCIATED',gsva.tis$pathway)),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('fib.png')
ggplot(gsva.tis[grep('VEGF',gsva.tis$pathway),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('vegf.png')

ggplot(gsva.tis[grep('FGF',gsva.tis$pathway),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('fgf.png')

gsva.mat<-gsva.tis%>%
    select(specimenID,pathway,score)%>%
    spread(key=pathway,value=score)%>%
    setRownames(specimenID)

#load up limma
#do design matrix with tumorType



