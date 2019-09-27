require(synapser)
synLogin()

#get drug data
source("../../bin/plotDrugsAcrossCells.R")


comps=findVariableDrugs()
res=lapply(subset(comps,response_type=='AUC_Simpson')$std_name[1:20],plotDoseResponseCurve)
res=lapply(subset(comps,response_type=='AUC_Simpson')$std_name[1:20],plotDrugByCellAndTumor)

res=lapply(subset(comps,response_type=='IC50_rel')$std_name[1:20],plotDoseResponseCurve)
res=lapply(subset(comps,response_type=='IC50_rel')$std_name[1:20],plotDrugByCellAndTumor)
plotDrugByCellAndTumor('SIROLIMUS')
plotDrugByCellAndTumor('IMATINIB')
plotDrugByCellAndTumor('SELUMETINIB')
