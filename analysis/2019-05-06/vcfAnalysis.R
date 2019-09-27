##vcf analysis
library(vcfR)
library(synapser)

synLogin()
vcf.with.genotype<-synGet('syn18079703')$path
vcf <- read.vcfR(vcf.with.genotype, verbose = FALSE)
#ref<- ape::read.dna(synGet('syn18082228')$path,format='fasta')
#gff <- read.table(gff_file, sep="\t", quote="")
#chrom <- create.chromR(name="Allsamples", vcf=vcf, seq=ref)
gt <- extract.gt(vcf, element = 'GT', as.numeric = TRUE)

##need to map variants to genes. 
