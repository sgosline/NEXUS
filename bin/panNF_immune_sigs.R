
library(synapser)
synLogin()

getMatrixFromTab<-function(tab){
  
}

#' getGeneList grabs the gene list from teh synapse table
#' @export 
getGeneList <- function(method='cibersort'){
  geneListTable <- 'syn12211688'
  # require(synapser,quietly=T)
  require(synapser)
#  synapse <- import("synapseclient")
#  syn <- synapse$Synapse()
  synLogin() 
  
  require(dplyr)
  
  tab <-synTableQuery(paste('select * from',geneListTable))$asDataFrame()%>%dplyr::select(Gene=`Gene Name`,Cell=`Cell Type`,Source,Operator)
  
  ##first make into a list of lists
  tab.list<-lapply(split(tab,tab$Source),function(x) lapply(split(x,x$Cell),function(y) y$Gene))
  
  if(method%in%(unique(tab$Source)))
    tab.list<-tab.list[[method]]
  
  #print(tab.list)
  
  return(tab.list)
  
}

runGsvaOnMat<-function(mat){
  gene.lists<-c('mcpcounter','cibersort','LyonsEtAl','Wallet','SchelkerEtAl')
  suppressMessages(library(GSVA))

  for(g in gene.lists){
    g.list=getGeneList(g)
    g.res<-gsva(as.matrix(mat),g.list,method='ssgsea',rnaseq=TRUE,verbose=FALSE)
    colnames(g.res)<-rownames(mat)
    
    pheatmap(t(g.res),cellwidth=10,clustering_distance_rows='correlation',clustering_distance_cols='correlation',clustering_method='ward.D2',annotation_row=cell.annotations[colnames(g.res),],show_rownames=F))
cat('\n\n') 
  }
}
