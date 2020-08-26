#dir('../OscopeData/dataAdjacencyMatrixe02/',pattern='.csv')
rm(list=ls())
inputdir="../Communities/Networks02"
outputdir="../Communities/Networks02/Overlapping"

# namef1="W45fastgreedy.txt"
# namef2="W123fastgreedy.txt"
# nameOutAll="W123_W45_fastgreedyI.txt"
# nameOutSig="W123_W45_fastgreedySigI.txt"

# namef1="W45walktrap.txt"
# namef2="W123walktrap.txt"
# nameOutAll="W123_W45_walktrapI.txt"
# nameOutSig="W123_W45_walktrapSigI.txt"

# namef1="W45infomap.txt"
# namef2="W123infomap.txt"
# nameOutAll="W123_W45_infomapI.txt"
# nameOutSig="W123_W45_infomapSigI.txt"

namef1="Sawalktrap.txt"
namef2="Sqwalktrap.txt"
nameOutAll="Sa_Sq_walktrapI.txt"
nameOutSig="Sa_Sq_walktrapSigI.txt"



file1=paste(inputdir,namef1,sep="/")#"../Communities/Networks02/W45fastgreedy.txt"
file2=paste(inputdir,namef2,sep="/")#"../Communities/Networks02/W123fastgreedy.txt"
outfileAll=paste(outputdir, nameOutAll,sep="/")
outfileSig=paste(outputdir, nameOutSig,sep="/")

library(igraph)
	M1=read.delim(file1,header=TRUE,sep="\t",quote="")
	comm1=M1[,2]
	ucomm1=sort(unique(comm1))
	nuc1=length(ucomm1)
	g1all=M1[,1]
	M2=read.delim(file2,header=TRUE,sep="\t",quote="")
	comm2=M2[,2]
	ucomm2=sort(unique(comm2))
	nuc2=length(ucomm2)
	g2all=M2[,1]
	uniong1g2=union(g1all,g2all)
	intersectg1g2=intersect(g1all,g2all)
	N=length(uniong1g2)
	Ni=length(intersectg1g2)
	nc1all=NULL
	indexAll=matrix(NA,nuc1*nuc2,6)
	count=0
	for (h in c(1:nuc1)){
		c1=M1[which(comm1==ucomm1[h]),1]
		c1i=intersect(c1,intersectg1g2)
        	nc1i=length(c1i)
        nc1=length(c1)
        for (k in c(1:nuc2)){
        	count=count+1
        	#perform hypergeometric test on c1h^c2k
        	c2=M2[which(comm2==ucomm2[k]),1]
        	c2i=intersect(c2,intersectg1g2)
        	nc2=length(c2)
        	nc2i= length(c2i)
        	commonc1c2=intersect(c1,c2)
        	ncommonc1c2=length(commonc1c2)
        	commonc1c2i=intersect(c1i,c2i)
        	ncommonc1c2i=length(commonc1c2i)
        	sprintf('h=%d k=%d intersez=%d',h,k, ncommonc1c2)
        ######tutte le palline sono le uniong1g2, le biache sono le c1 e l'estrazione e' lunga c2 e il numero di bianche osservate e' c1^c2
        	#pval=phyper(ncommonc1c2,nc1,N-nc1,nc2)
        	pval=phyper(ncommonc1c2i,nc1i,Ni-nc1i,nc2i,lower.tail=FALSE)

        #	indexAll[count,]=c(ucomm1[h],ucomm2[k],nc1,nc2,ncommonc1c2,pval)
         	indexAll[count,]=c(ucomm1[h],ucomm2[k],nc1i,nc2i,ncommonc1c2i,1)
        if(ncommonc1c2i>0){
        	indexAll[count,]=c(ucomm1[h],ucomm2[k],nc1i,nc2i,ncommonc1c2i,pval)}

        	


        	}
       }
        

					
#}

vetp=indexAll[,6]
sortp=sort(vetp,index.return=TRUE)
ntot=length(vetp)
pvaladjvec=vetp*ntot/sortp$ix
pvaladjvec[pvaladjvec>=1]=1
summaryAll=cbind(indexAll,pvaladjvec)
colnames(summaryAll)=c(file1,file2,"size c1","size c2","size common","pval","adj_pval")
summarySig=summaryAll[which(summaryAll[,"adj_pval"]<0.05),]
write.table(summaryAll,file=outfileAll,sep="\t",quote=FALSE)
write.table(summarySig,file=outfileSig,sep="\t",quote=FALSE)


