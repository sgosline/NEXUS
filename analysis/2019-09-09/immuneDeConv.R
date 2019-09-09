###reprocess RNA-seq

require(tidyverse)

#######eventually break these out into separate files!!!
require(synapser)
synLogin()

##get table query
tab<-synTableQuery('SELECT * FROM syn20449214 WHERE ( ( "isCellLine" IS NULL ) OR "isCellLine" = \'FALSE\' ) AND ("transplantationType" is NULL OR "transplantationType" = \'\')')$asDataFrame()

tab$tumorType[is.na(tab$tumorType)]<-'None'
tab$diagnosis[is.na(tab$diagnosis)]<-'None'
tab$sex<-tolower(tab$sex)


prefix='allPublicTissueData'



runImmuneDeconv<-function(tab,method){
  #run MCP counter

  mat<-reshape2::acast(tab,Symbol~specimenID,value.var='totalCounts',fun.aggregate=mean,na.rm=T)
  nas<-which(apply(mat,1,function(x) any(is.na(x))))
  if(length(nas)>0)
    mat<-mat[-nas,]

    library(immunedeconv)  
  res<-deconvolute(mat,method)
  
  df<-dplyr::select(tab,c(study,sex,tumorType,specimenID))%>%unique()
  rownames(df)<-df$specimenID

  #save as heatmap with metadata
  mtab<-res%>%select(-cell_type)%>%as.data.frame()
  rownames(mtab)<-res$cell_type

  library(pheatmap)
  pheatmap(log2(mtab+0.01),annotation_col=select(df,-specimenID),
      cellheight = 10,cellwidth=10,
      file=paste0(prefix,'_',method,'Preds.pdf'),
    labels_col=rep(" ",ncol(mtab)))
  
  ##now tidy up data to table
  td<-tidyr::gather(res,key="specimenID",value="score",-cell_type )%>%
    left_join(df,by='specimenID')
  td$method=method
  return(td)

}

#store
this.script='https://raw.githubusercontent.com/sgosline/NEXUS/master/analysis/2019-09-09/immuneDeConv.R'

set_cibersort_binary('../../bin/CIBERSORT.R')
set_cibersort_mat('../../bin/LM22.txt')
synapse_table='syn20710536'
for(m in c('cibersort')){
  res<-runImmuneDeconv(tab,m)
  synapser::synStore(synapser::Table(synapse_table,res),used='syn20449214',executed=this.script)
}


parentid='syn20710537'
for(fi in list.files('.')[grep('Preds.pdf',list.files('.'))])
  synapser::synStore(synapser::File(fi,parentId=parentid,annotations=list(resourceType='analysis',isMultiSpecimen='TRUE',isMultiIndividual='TRUE')),used='syn20449214',executed=this.script)
