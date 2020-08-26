##############################################################################
#This code will merge a selected list of genes in "nmeth.3549-S8.txt", with the data in your input data file "GSE64016_H1andFUCCI_normalized_EC.csv". The output is stored in "mergedH1.csv"
###############################################################################
H1Oscope=read.delim("../casestudy/H1/nmeth.3549-S8.txt",quote="",sep='\t')
H1Oscope=as.matrix(H1Oscope)

XX=H1Oscope
print('Subselected number of genes:')
print(length(XX))

#######Read DATA
data=read.csv("../casestudy/H1/GSE64016_H1andFUCCI_normalized_EC.csv")
A=merge(XX,data,by.x=1,by.y='X')
####selecting H1 experiment data
H1=A[,c(2:214)]
row.names(H1)=A[,1]
###Saving the output in mergedH1.csv
write.csv(H1,file='../casestudy/H1/mergedH1.csv')
