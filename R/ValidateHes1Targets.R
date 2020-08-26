###########validate Hes1 targets
#netfile='../Communities/Networks02/W45infomap.txt'
#N=read.table(netfile,sep='\t',header=TRUE,quote="")
#pimW45=merge(A,N,by.x=1, by.y=1)

#######################
#This piece of code gives in output also which are the communities you want to look for the targets
############
target= 'HES1'#'CC'#'HES1'
dirlist=c('../Communities/Networks02/','../Communities/NetworksPiDiv9/','../Communities/FirstData/')

#communityDir='../Communities/Networks02/'
#communityDir='../Communities/NetworksPiDiv9/'
#communityDir='../Communities/FirstData'

if (target=='CC'){
	#targetfile='../OscopeData/CCrelated.txt'
    targetfile='../OscopeData/CCrelatedsymbol.txt'


}
if (target=='HES1'){
	targetfile='../OscopeData/Hes1Targets.txt'

}
for (h in c(1:3)){
	
	communityDir=dirlist[h]
	print('##############')
	print(communityDir)
    print('##############')

    A=read.delim(targetfile)
    netfiles=dir(communityDir,pattern='.txt')
    for (j in c(1:length(netfiles))){
    	dumf=netfiles[j]
    	print(dumf)
    	outfile=paste(communityDir,target,'/',dumf,sep='')
    	##########scelta annotation file
    	if (length(grep('s123',dumf))>0){
    		   genedet='../OscopeData/GenesLabels/Mapping/Water123.txt'}
    		   if (length(grep('s45',dumf))>0){
    		   	 genedet='../OscopeData/GenesLabels/Mapping/Water45.txt'}
    		   	 if (length(grep('aNSC',filename))>0){
    		   	 	genedet='../OscopeData/GenesLabels/Mapping/SVZa.txt'}
    		   	 	if (length(grep('qNSC',filename))>0){
    		   	 		genedet='../OscopeData/GenesLabels/Mapping/SVZq.txt'}
	
	netfile=paste(communityDir,netfiles[j],sep='')
	print(netfile)
	N=read.table(netfile,sep='\t',header=TRUE,quote="")
	common=merge(A,N,by.x=1,by.y='mgi_symbol')
	print('total genes')
    print(length(N[,1]))
	print('overlapping with targets')
	print(length(common[,1]))
	write.table(common,file=outfile,sep='\t',quote=FALSE,row.names=FALSE)
	}
	}

#######################
#This piece of code gives in output only if there is an intersection of the target list with the network
############

library(igraph)

target='HES1'#'CC'

expdirList=c('../OscopeData/dataAdjacencyMatrixe02/','../OscopeData/dataAdjacencyMatrixePiDiv9/','../OscopeData/FirstData/')
if (target=='HES1'){
	outputList=c('../OscopeData/Hes1OverlapNet02.txt','../OscopeData/Hes1OverlapNetPiDiv9.txt','../OscopeData/Hes1OverlapNetNoFilter.txt')
}
if (target=='CC'){
	outputList=c('../OscopeData/CCOverlapNet02.txt','../OscopeData/CCOverlapNetPiDiv9.txt','../OscopeData/CCOverlapNetNoFilter.txt')
}

for (h in c(1:3)){
	expdir=expdirList[h]
	output=outputList[h]
#expdir='../OscopeData/dataAdjacencyMatrixe02/'
#output='../OscopeData/Hes1OverlapNet02.txt'
#output='../OscopeData/CCOverlapNet02.txt'

#expdir='../OscopeData/dataAdjacencyMatrixePiDiv9/'
#output='../OscopeData/Hes1OverlapNetPiDiv9.txt'
#output='../OscopeData/CCOverlapNetPiDiv9.txt'

#expdir='../OscopeData/FirstData/'
#output='../OscopeData/Hes1OverlapNetNoFilter.txt'
#output='../OscopeData/CCOverlapNetNoFilter.txt'


if (target=='CC'){
	#targetfile='../OscopeData/CCrelated.txt'
    targetfile='../OscopeData/CCrelatedsymbol.txt'


}
if (target=='HES1'){
	targetfile='../OscopeData/Hes1Targets.txt'

}

A=read.delim(targetfile)



filelist=dir(expdir,pattern='.csv')
lf=length(filelist)
M=matrix('NA',lf,3)
for (i in c(1:lf)){
	
	dumf=filelist[i]

	
	filename=paste(expdir,filelist[i],sep='')
	print(filename)
	
	##########scelta annotation file
	if (length(grep('s123',dumf))>0){
	   genedet='../OscopeData/GenesLabels/Mapping/Water123.txt'
	   data=read.csv(file=filename, header= TRUE,row.names=1)
	   G=graph_from_data_frame(data,directed=FALSE)
       g1name=get.vertex.attribute(G)$name
       details=read.table(genedet,sep="\t",header=TRUE,quote="")
       merg1=merge(g1name,details,by.x=1,by.y='mgi_symbol',all.x=TRUE)
       colnames(merg1)=c('mgi_symbol','ensembl_gene_id','description')

	}
	if (length(grep('s45',dumf))>0){
       genedet='../OscopeData/GenesLabels/Mapping/Water45.txt'
       data=read.csv(file=filename, header= TRUE,row.names=1)
	   G=graph_from_data_frame(data,directed=FALSE)
       g1name=get.vertex.attribute(G)$name
       details=read.table(genedet,sep="\t",header=TRUE,quote="")
       merg1=merge(g1name,details,by.x=1,by.y='mgi_symbol',all.x=TRUE)
       colnames(merg1)=c('mgi_symbol','ensembl_gene_id','description')
	}
	if (length(grep('aNSC',filename))>0){
       genedet='../OscopeData/GenesLabels/Mapping/SVZa.txt'
       data=read.csv(file=filename, header= TRUE,row.names=1)
	   G=graph_from_data_frame(data,directed=FALSE)
       g1name=get.vertex.attribute(G)$name
       details=read.table(genedet,sep="\t",header=TRUE,quote="")
       merg1=merge(g1name,details,by.x=1,by.y='ensembl_gene_id',all.x=TRUE)
       colnames(merg1)=colnames(details)
	}
	if (length(grep('qNSC',filename))>0){
	   genedet='../OscopeData/GenesLabels/Mapping/SVZq.txt'
	   data=read.csv(file=filename, header= TRUE,row.names=1)
	   G=graph_from_data_frame(data,directed=FALSE)
       g1name=get.vertex.attribute(G)$name
       details=read.table(genedet,sep="\t",header=TRUE,quote="")
       merg1=merge(g1name,details,by.x=1,by.y='ensembl_gene_id',all.x=TRUE)
       colnames(merg1)=colnames(details)
    }
	
	
	
	
	
# # 	data=read.csv(file=filename, header= TRUE,row.names=1)
	# G=graph_from_data_frame(data,directed=FALSE)
    # g1name=get.vertex.attribute(G)$name
    
    # details=read.table(genedet,sep="\t",header=TRUE,quote="")
    # merg1=merge(g1name,details,by.x=0,by.y='ensembl_gene_id',all.x=TRUE)


    
    
    print('total genes')
    lg1=length(g1name)
    print(lg1)
    print('overlapping with targets')
    common=intersect(A[[1]],merg1$mgi_symbol)
	ncommon=length(common)

    
    # if (target=='CC'){
    	    # common=intersect(A[[1]],merg1$ensembl_gene_id)
    	    # ncommon=length(common)
    # }
   # if (target=='HES1'){
   	  # common=intersect(A[[1]],merg1$mgi_symbol)
	  # ncommon=length(common)

    # }

    print(ncommon)
    print(common)
    M[i,1]=lg1
    M[i,2]=ncommon
    aaa=capture.output(cat(c(common)))
    if (    length(aaa)>0){
    	    M[i,3]=aaa

    }
    #MM=data.frame(M)
    #MM$targetsName[1]=aaa

}
colnames(M)=c('co-osclillating','#targets',"target Names")
rownames(M)=filelist
write.table(M,output,sep='\t',quote=FALSE)
}