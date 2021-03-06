###reprocess RNA-seq

require(tidyverse)

sampleCounts<-function(tab,prefix){
  counts<-dplyr::select(tab,c('specimenID','individualID','nf1Genotype','nf2Genotype','tumorType','diagnosis','isCellLine','transplantationType'))%>%
    group_by(.dots=c('diagnosis','tumorType','nf1Genotype','nf2Genotype','isCellLine','transplantationType'))%>%
    summarize(individuals=n_distinct(individualID),samples=n_distinct(specimenID))
  counts
}

# PCA plot
plotPCA<-function(dds,prefix=''){
  #tab, xlim=0,ylim=0,scale=FALSE,prefix=''){
  require(ggplot2)
  vsd <- vst(dds, blind = FALSE)
  
  p <- DESeq2::plotPCA(vsd, intgroup = c("individualID", "tumorType", 
    "sampleType"), returnData = TRUE)
  
  percentVar <- round(100 * attr(p, "percentVar"))
  pdf(paste(prefix,'RNASeqDataPCA.pdf',sep=''),width='960',height='960')
  
  p2<-ggplot(p, aes(PC1, PC2, color = tumorType, shape = sampleType, 
    label = vsd$individualID)) + geom_point(size = 3) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    coord_fixed() + ggrepel::geom_text_repel()
  
  print(p2)
  dev.off()
  
}

plotCounts<-function(tab,prefix=''){
  require(ggplot)
  pdf(paste0(prefix,'rnaSeqCountsBoxplot.pdf'),width='960',height='960')
  p<-ggplot(tab)+geom_boxplot(aes(x=specimenID,y=zScore,fill=tumorType))+scale_y_log10()+theme(axis.text.x = element_text(angle = 90)) 
  print(p)
  dev.off()
}




runGSVA<-function(genes.with.meta,prefix=''){
  
  ###now do gsva
  library(GSVA)
  library(GSVAdata)
  library(pheatmap)
  mat<-reshape2::acast(genes.with.meta,Symbol~specimenID,value.var='zScore',fun.aggregate=mean)
  missing<-which(apply(mat,1,function(x) any(is.na(x))))
  mat<-mat[-missing,]
  data("c2BroadSets")
  
  library(org.Hs.eg.db)
  
  map<-AnnotationDbi::select(org.Hs.eg.db,columns=c("SYMBOL",'ENTREZID'),keys=keys(org.Hs.eg.db,'ENTREZID'),multiVals=unique(tab$Symbol))
  
  entrez<-map$ENTREZID[match(rownames(mat),map$SYMBOL)]
  mat<-mat[which(!is.na(entrez)),]
  rownames(mat)<-entrez[!is.na(entrez)]
  res=gsva(mat,method='ssgsea',gset.idx.list=c2BroadSets)
  library(pheatmap)
  vars<-apply(res,1,var)
  annotes=genes.with.meta%>%dplyr::select(specimenID,sex,tumorType,diagnosis,isCellLine,transplantationType)%>%unique()
  rownames(annotes)<-annotes$specimenID
  
  pheatmap(res[names(sort(vars,decreasing=T)[1:50]),],labels_col=rep("",ncol(res)),fontsize_row = 8,clustering_method = 'ward.D2',annotation_col = dplyr::select(annotes,-specimenID),cellwidth = 10,cellheight = 10,filename=paste0(prefix,'top50GsvaPaths.pdf'))
  
  res
}

tumorTypeDiffEx<-function(tab){
   library(org.Hs.eg.db)
  require(clusterProfiler)
  require(DESeq2)
  
  samps<-dplyr::select(tab,tumorType,specimenID,sampleType,individualID)%>%distinct()
  rownames(samps)<-samps$specimenID
  samps<-samps%>%dplyr::select(tumorType,sampleType)
  mat<-reshape2::acast(tab,Symbol~specimenID,value.var='totalCounts',fun.aggregate=sum)
  
  
  map<-AnnotationDbi::select(org.Hs.eg.db,columns=c("SYMBOL",'ENTREZID'),keys=keys(org.Hs.eg.db,'ENTREZID'),multiVals=unique(tab$Symbol))
  
  OrgDb <- org.Hs.eg.db  # can also be other organisms
  full.dat<-NULL
  for(con in unique(tab$tumorType)){
    samps$isTumor<-sapply(samps$tumorType,function(x) x==con)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(mat[,rownames(samps)]), 
      colData =samps,design = ~isTumor)
    #copied from Xengie's markdwon
    ### filter out reads with low counts
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    
    dds <- DESeq2::DESeq(dds)
    
    res <- results(dds, alpha = 0.05)
  
  genes<-map$ENTREZ[match(rownames(subset(res,padj<0.05)),map$SYMBOL)]
  genes <- na.omit(genes)
  
  ego <- clusterProfiler::enrichGO(gene = genes, OrgDb = OrgDb, 
    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05, readable = TRUE)
  
  kk <- clusterProfiler::enrichKEGG(gene = genes, organism = "hsa", 
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  combined<-rbind(data.frame(paths=rep('GO',nrow(ego)),ego),
      data.frame(paths=rep("KEGG",nrow(kk)),kk))
  message(paste(dim(combined),collapse=','))
  combined$condition=rep(con,nrow(combined))
  message(paste(dim(combined),collapse=','))
  combined$numGenes=rep(as.character(length(genes)),nrow(combined))
  full.dat<-rbind(full.dat,combined)
  }
  return(full.dat)

}




getDDS<-function(tab){
  require(DESeq2)
  detach('package:synapser', unload=TRUE)
  unloadNamespace('PythonEmbedInR')
  samps<-dplyr::select(tab,tumorType,specimenID,sampleType,individualID)%>%distinct()
  rownames(samps)<-samps$specimenID
  mat<-reshape2::acast(tab,Symbol~specimenID,value.var='totalCounts',fun.aggregate=sum)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(mat[,samps$specimenID]), 
    colData =samps,design = ~tumorType)
  #copied from Xengie's markdwon
  ### filter out reads with low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  dds
}


#######eventually break these out into separate files!!!
require(synapser)
synLogin()

##get table query
tab<-synTableQuery('select * from syn20449214 where ( ( "studyName" = \'A Nerve Sheath Tumor Bank from Patients with NF1\' ) )')$asDataFrame()

tab$sampleType<-apply(tab,1,function(x){
  if(!is.na(x[['transplantationType']]) && x[['transplantationType']]=='xenograft')
    return('xenograft')
  else if (!is.na(x[['isCellLine']])&&x[['isCellLine']]==TRUE)
    return('cell line')
  else
    return('tissue')})

tab$specimenID<-sapply(tab$specimenID,function(x)
  gsub('Neurofibroma','NF',
    gsub("Plexiform ",'p',
      gsub('Cutaneous ',"c",
        gsub('Malignant Peripheral Nerve Sheath Tumor','MPNST',x)))))

prefix='jhuBiobank'
#get sample counds
countsTab<-sampleCounts(tab,prefix)
write.csv(countsTab,paste0(prefix,'countsTab.csv'))


#plotCounts
counts<-plotCounts(tab,prefix)

#doGSVA
#gs<-runGSVA(tab,prefix)
#write.csv(gs,paste0(prefix,'GSVAPathwayScores.csv'))

#doDiffEx
diffex<-tumorTypeDiffEx(tab)
write.csv(diffex,paste0(prefix,'EnrichedGenesPaths.csv'))

dds<-getDDS(tab)

#doPCA
res<-plotPCA(dds,prefix)

#store
this.script='https://raw.githubusercontent.com/sgosline/NEXUS/master/analysis/2019-08-07/rnaSeqStandardAnalysis.R'

synapser::synLogin()
parentid='syn20683372'
for(fi in list.files('.'))
  synapser::synStore(synapser::File(fi,parentId=parentid,annotations=list(resourceType='analysis',diagnosis='Neurofibromatosis 1',isMultiSpecimen='TRUE',isMultiIndividual='TRUE')),used='syn20449214',executed=this.script)
