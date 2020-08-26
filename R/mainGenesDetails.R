#####you need to install Biomart if you dont have it and load the library
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#source("http://bioconductor.org/biocLite.R")
#biocLite("GO.db")
library(GO.db)
######################
library(biomaRt)
#####################Read the file name containing the genenames:
#genesdir='../GenesLabels/'
genesdir='../OscopeData/GenesLabels/'
genesout='../OscopeData/GenesLabels/Mapping/'
lf=dir(genesdir,pattern=".txt")
###we can make up a loop here
# i=1
# print(c('reading file: ',lf[i]))
# filename=paste(genesdir,lf[i],sep='')
# out='SVZa.txt'
# outnm='SVZaNotMapped.txt'
# i=2
# print(c('reading file: ',lf[i]))
# filename=paste(genesdir,lf[i],sep='')
# out='SVZq.txt'
# outnm='SVZqNotMapped.txt'
# #filename='../GenesLabels/WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s123res_C73G9200_geneNames.txt'
for (i in c(1:4)){
	
print(c('reading file: ',lf[i]))
filename=paste(genesdir,lf[i],sep='')
if (length(grep('s123',filename))>0){
	out='Water123.txt'
	outnm='Water123NotMapped.txt'}
	if (length(grep('s45',filename))>0){
	out='Water45.txt'
	outnm='Water45NotMapped.txt'}
	if (length(grep('aNSC',filename))>0){
	out='SVZa.txt'
	outnm='SVZaNotMapped.txt'}
	if (length(grep('qNSC',filename))>0){
	out='SVZq.txt'
	outnm='SVZqNotMapped.txt'}
	
#
genes=read.csv(filename, header = FALSE, sep = ",")
genes=as.matrix(genes)
genes=as.vector(genes)

ng=length(genes)

###choose among available datasets
#ensembl=useMart("ensembl")
#listDatasets(ensembl)
      
#ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")


                   

#vv=unlist(strsplit(res_tot$Row.names[],'[.]'))
#ind=seq(1,length(vv),2)
#genes=vv[ind]



##########Connect to biomart to retrieve the available informations



mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))


 
 



if (i==1|i==2) {
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
"mgi_symbol","description"),values=genes,mart= mart)} else {
	G_list <- getBM(filters= "mgi_symbol", attributes= c("ensembl_gene_id",
"mgi_symbol","description"),values=genes,mart= mart)}	
gooo=G_list[,3]
G_list$Term=Term(gooo)

     
detail=as.matrix(G_list)
x=genes
y=detail[,1]
Nonmapped=x[is.na(match(x,y))]

# aaa=match(genes,G_list$ensembl_gene_id)
# inda=aaa[!(is.na(aaa))]
# #indb=aaa[(is.na(aaa))]
# M=matrix(data=NA,ng,3)
# M[inda,1]=G_list$ensembl_gene_id
# M[inda,2]=G_list$mgi_symbol
# M[inda,3]=G_list$description

# rownames(M)=genes
# dum=merge(M,res_tot,by.x="row.names",by.y="Row.names",all.x=TRUE)
fout=paste(genesout,out,sep='')
foutnm=paste(genesout,outnm,sep='')


write.table(detail,fout,sep="\t",row.names=FALSE,quote=FALSE)
write.table(Nonmapped,foutnm,sep="\t",row.names=FALSE,quote=FALSE)

}
