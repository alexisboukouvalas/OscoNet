
library(STRINGdb)
 string_db <- STRINGdb$new( version="10", species=10090,score_threshold=0, input_directory="" )
 data(diff_exp_example1)
 example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = FALSE)
 hits <- example1_mapped$STRING_id[1:80]

enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE ) 
enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE ) 
head(enrichmentGO, n=7)

eh <- string_db$enrichment_heatmap( list( hits[1:100], hits[101:200]), list("list1","list2"), title="My Lists" )



##############

filecomm='../Communities/Networks02/W45fastgreedy.txt'
thresh=8


 string_db <- STRINGdb$new( version="10", species=10090,score_threshold=0, input_directory="" )

AAA=read.delim(filecomm)
AAA_mapped=string_db$map( AAA, "mgi_symbol", removeUnmapped=FALSE)
#backgroundV <- AAA_mapped$STRING_id
#string_db$set_background(backgroundV)
uc=unique(AAA_mapped$Community.Id)
nuc=length(uc)
universe=as.vector(AAA_mapped$mgi_symbol)
gene_sets=list()
names=NULL
sizes=NULL

for (i in c(1:nuc)){
	dum=which(AAA_mapped$Community.Id%in%uc[i])
	ld=length(dum)
	if(ld>=thresh){
		print(i)
		group= AAA_mapped$STRING_id[dum]
		gene_sets[length(gene_sets)+1]=list(group)
		names=c(names,paste("comm",uc[i],sep=''))
		sizes=c(sizes,ld)
		}
		}
		names(gene_sets)=names
		ngs=length(names)
		#################set background
#backgroundV <- example1_mapped$STRING_id[1:2000] # as an example, we use the first 2000 genes > string_db$set_background(backgroundV) You can also set the background when you instantiate the STRINGdb #object: > #string_db <- STRINGdb$new( score_threshold=0, backgroundV = backgroundV ) I
		GoRes=list()
		KeggRes=list()
		for (i in c(1:ngs)){
			enrichmentGO <- string_db$get_enrichment( gene_sets[[i]], category = "Process", methodMT = "fdr", iea = TRUE ) 	
			enrichmentKEGG <- string_db$get_enrichment( gene_sets[[i]], category = "KEGG", methodMT = "fdr", iea = TRUE ) 
			GoRes[length(GoRes)+1]=list(enrichmentGO)
			KeggRes[length(KeggRes)+1]=list(enrichmentKEGG)

		}
	
#eh <- string_db$enrichment_heatmap(gene_sets[1:2],list(names[1:2]),title="netwa")