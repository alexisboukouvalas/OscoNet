
# Dircom='../Communities/Networks02/'
# file='W45fastgreedy.txt'
# ont='BP'
# Nlimit=2
# plim=0.05
queryGo<-function(Dircom,file,ont,Nlimit,plim){
######
# 
#Input Dircom: Directory base for the communities annotated 
#			   file is the name of the community annotate file
#			   ont is the term to show enrichment (BP, CC, MF)
#              Nlimit is the min num of elements present in GO/Kegg term per community
#              Plim is the threshold for the pvalues (we should chnge to FDR that we compute but we are less stringent now)
#Output list: selected, list containing a summary for each community top GO (pval<0.05, count>2)
#             Max, the best GO (pval<0.05, max count)
######

filecomm=paste(Dircom,file, sep='')
#outGo is the file output of the mainGOana.R that performs a GO and Kegg Enrichment 
outGo=paste(Dircom,'GO/Go_',file,sep='')
A=read.delim(outGo,header=TRUE,sep='\t')
iBP=which(A$Ont=='BP')
Abp=A[iBP,]
ipcom=grep('P.',colnames(A))
nc=length(ipcom)
selected=list()
Max=matrix(NA,nrow=nc,ncol=2)
for( i in c(1:nc)){
	dumC=Abp[,c(c(1:4),c((ipcom[i]-nc),ipcom[i]))]
	pval=dumC[,6]
	dumC$FDR=p.adjust(pval,'BH')
	isig=which(dumC[,6]<plim)
	sub=dumC[isig,]
	iN=which(sub[,5]>=Nlimit)
	dumsub=sub[iN,]
	selected[length(selected)+1]=list(dumsub)
	im=which.max(dumsub[,5])
	Max[i,]=as.matrix(dumsub[im,c(1:2)])
	rownames(Max)=colnames(A)[ipcom-nc]
	}
	myout=list(selected,Max)
	return(myout)
	}