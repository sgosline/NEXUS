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

ggplot(gsva.tid[grep('COMPLEMENT',gsva.tid$pathway),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gsva.tid[grep('TARGETS',gsva.tid$pathway),])+geom_boxplot(aes(x=pathway,y=score,fill=tumorType))+theme(axis.text.x = element_text(angle = 45, hjust = 1))



