library('org.Mm.eg.db')
##########preso da sclvm on GITHUB
getEnsembl <- function(term, species = 'mMus'){
if(!(species %in%c('mMus','Hs'))){stop("'species' needs to be either 'mMus' or 'Hs'")}

if(species=='mMus'){
if(require(org.Mm.eg.db)){
xxGO <- AnnotationDbi::as.list(org.Mm.egGO2EG)
#org.Mm.egENSEMBL
x <- org.Mm.egSYMBOL}else{
stop("Install org.Mm.eg.db package for retrieving gene lists from GO")
}
}else{
if(require(org.Hs.eg.db)){
xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG)
#org.Hs.egENSEMBL for ensamble
x <-org.Hs.egSYMBOL}else{
stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
}
}
cell_cycleEG <-unlist(xxGO[term])
#get ENSEMBLE ids

mapped_genes <- mappedkeys(x)
xxE <- as.list(x[mapped_genes])
ens_ids_cc<-unlist(xxE[cell_cycleEG]) 

ens_ids_cc
}

######Esempio d'uso:
######any GO term! Questo e' quello per CC!
# symbol_ids_cc= getEnsembl('GO:0007049')


# write.table(ens_ids_cc,'../OscopeData/CCrelatedsymbol.txt',sep='\t',quote=FALSE)
