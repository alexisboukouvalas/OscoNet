
GenesDescr =read.table(genedet,sep='\t',header=TRUE, quote="")
mygeneUniverse = as.vector(GenesDescr[,1]) 
 myselectedGenes=sample(mygeneUniverse,10)
 
 
 library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("refseq_dna", "go_molecular_function_id",
           "go_molecular_function__dm_name_1006"), filters = "refseq_dna",
           values = c("NM_030621"), mart = mart)
######serve il 
sampleGOdata <- new("topGOdata",
                      description = "Simple session", ontology = "BP",
                      allGenes = mygeneUniverse, geneSel = myselectedGenes,
                      nodeSize = 10,
                      annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")

