##compare tumor immune scores to metaviper predictions

require(synapser)
synLogin()
require(tidyverse)

deconv_scores='syn20710536'
metaviper_scores='syn20503291'

#get immune predictions
dtab<-synapser::synTableQuery(paste('select * from',deconv_scores))$asDataFrame()%>%
  rename(immScore='score')

##get metaviper scores
mtab<-read.csv(synapser::synGet(metaviper_scores)$path,sep='\t')%>%
  rename(specimenID='sample')%>%
  rename(protScore='counts')


dtab<-subset(dtab,method!='xcell')
combined=dtab%>%select(c(cell_type,method,specimenID,immScore))%>%
  left_join(mtab,by='specimenID')

corVals=combined%>%group_by(cell_type,gene,method)%>%summarize(corVal=cor(immScore,protScore,use='pairwise.complete.obs'))

##now how do we bracket them?
##plot correlation distributions by cell type and method. 
require(ggplot2)
p<-ggplot(corVals)+geom_boxplot(aes(x=cell_type,y=corVal,fill=method))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("Correlation of metaviper proteins with predicted cell type")
print(p)
ggsave('metaviper_cellType_correlations.pdf')
corthresh=0.65

##now filter to the cell types with correlated proteins
cor_cell_types=subset(corVals,corVal>corthresh)%>%ungroup()%>%
  select(cell_type,method)%>%unique()
print(paste('we found',nrow(cor_cell_types),'cell types with some protein correlation greater than',corthresh))


apply(cor_cell_types,1,function(x){
  ct=x[['cell_type']]
  m=x[['method']]

  #for each gene and cell type
  genes=subset(corVals,cell_type==ct)%>%
        subset(corVal>corthresh)%>%
    subset(method==m)%>%arrange(desc(corVal))%>%
      ungroup()

    if(nrow(genes)>12){
    new.corthresh=format(genes$corVal[15],digits=3)
    genes=genes[1:12,]
  }else{
    new.corthresh=corthresh
  }

  scores=subset(combined,gene%in%genes$gene)%>%subset(cell_type==ct)%>%subset(method==m)

  p2<- ggplot(scores)+
      geom_point(aes(x=immScore,y=protScore,
          col=gene,shape=conditions))+
    scale_x_log10()+
      ggtitle(paste(m,'predictions of',ct,'correlation >',new.corthresh))
  print(p2)
  ggsave(paste0(m,'predictions of',gsub(" ","",gsub("/","",ct)),'cor',new.corthresh,'.pdf'))
})

this.script='https://raw.githubusercontent.com/sgosline/NEXUS/master/analysis/2019-09-26/metaviperDeconvCorrelations.R'
parentid='syn20710537'
for(fi in list.files('.')[grep('tions',list.files('.'))])
  synapser::synStore(synapser::File(fi,parentId=parentid,annotations=list(resourceType='analysis',isMultiSpecimen='TRUE',isMultiIndividual='TRUE')),used=c(deconv_scores,metaviper_scores),executed=this.script)


