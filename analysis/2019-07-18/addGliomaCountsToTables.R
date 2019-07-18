###
require(synapser)
synLogin()
require(tidyverse)
#this gets the gzipped counts files
gliomanf1=synTableQuery("SELECT id,specimenID,individualID,species,sex,age,tumorType,isCellLine,study,diagnosis,dataType,consortium,studyName,fundingAgency,resourceType,nf1Genotype,nf2Genotype FROM syn11614207 WHERE ( ( \"assay\" = 'rnaSeq' ) )")$asDataFrame()

###COlumbia Glioma
gli.genes=do.call('rbind',lapply(gliomanf1$id,function(x){
  ent=synGet(x)
  f=ent$path
  tab<-read.table(gzfile(f),header=F)
  colnames(tab)<-c('ensembl','Symbol','Counts')
  tab<-tab%>%group_by(Symbol)%>%summarize(totalCounts=sum(Counts))
  tab$zScore=(tab$totalCounts-mean(tab$totalCounts,na.rm=T))/sd(tab$totalCounts,na.rm=T)
  
  data.frame(tab,
      used=rep(x,nrow(tab)),
      path=rep(basename(f),nrow(tab)),
      parent=rep(ent$properties$parentId,nrow(tab)))
        
}))

full.tab<-gli.genes%>%left_join(rename(gliomanf1,used='id'),by='used')

pub.table='syn20449214'
synStore(Table(pub.table,full.tab))
priv.table='syn20370978'
synStore(Table(priv.table,full.tab))

