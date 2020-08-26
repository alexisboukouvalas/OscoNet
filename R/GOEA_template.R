require("categoryCompare")
require("Homo.sapiens")
require("KEGG.db")

do.GOEA = function(cls,univ) {
  univ = subset(select(Homo.sapiens, keys=univ, columns=c("ENTREZID"), keytype="SYMBOL"),!is.na(ENTREZID))
  geneLists = vector(mode = "list",length = length(cls))
  names(geneLists) = names(cls)
  for (i in 1:length(geneLists)) {
    ann <- subset(select(Homo.sapiens, keys=cls[[i]], columns=c("ENTREZID"), keytype="SYMBOL"),!is.na(ENTREZID))
    geneLists[[i]] <- list(genes=unique(ann$ENTREZID), universe=unique(univ$ENTREZID), annotation='org.Hs.eg.db')
  }
  
  geneLists <- new("ccGeneList", geneLists, ccType=c("BP","CC","MF","KEGG"), pvalueCutoff=1, fdr=100 )
  print("Enrichment ...")
  enrichLists <- ccEnrich(geneLists)
  return(enrichLists)
}

get.table = function(enrichLists,proc,i,th=0.05,th_fdr = 1) {
  t = summary(enrichLists[[proc]][[i]])
  t$FDR = p.adjust(t$Pvalue,method = "fdr")
  t = t[,-c(4,5)]
  if (nrow(subset(t,Pvalue < th & FDR < th_fdr))==0) {
    t = subset(t,Pvalue < min(th,0.01))
  } else {
    t = subset(t,Pvalue < th & FDR < th_fdr)
  }
  
  if (proc == "KEGG") {
    xx <- as.list(KEGGPATHID2NAME)
    t$desc = unlist(xx[t[,1]])
  } else {
    desc = select(Homo.sapiens, keys=t$ID, columns=c("TERM"), keytype="GO")
    t$desc = desc$TERM[match(t$ID,desc$GO)]
  }
  
  t$groupID = names(enrichLists[["BP"]])[i]
  return(t)
}



get.df = function(enrichLists,th=0.05,th_fdr = 1) {
  n = length(enrichLists[["BP"]])
  
  df.BP = NULL
  df.CC = NULL
  df.MF = NULL
  df.KEGG = NULL
  for (i in 1:n) {
    df.BP = rbind(df.BP,get.table(enrichLists,"BP",i,th,th_fdr))
    df.CC = rbind(df.CC,get.table(enrichLists,"CC",i,th,th_fdr))
    df.MF = rbind(df.MF,get.table(enrichLists,"MF",i,th,th_fdr))
    df.KEGG = rbind(df.KEGG,get.table(enrichLists,"KEGG",i,th,th_fdr))
  }
  
  l = list(df.BP,df.CC,df.MF,df.KEGG)
  names(l) = c("BP","CC","MF","KEGG")
  return(l)
}


## EXAMPLE
universe <- c("TP53","MYC","TFEB","DDIT3","APC","DDIT4","XBP1","MCLN1","MLH1")
gene_set1 <- c("DDIT3","XBP1","TP53")
gene_set2 <- c("MYC","TP53")

gene_sets <- list(gene_set1,gene_set2)
names(gene_sets) <- c("set1","set2")

enrichLists <- do.GOEA(gene_sets,universe)
results <- get.df(enrichLists,th = 2,th_fdr = 2)

names(results)

head(results[["BP"]])

head(subset(results[["KEGG"]],groupID=="set2"))




