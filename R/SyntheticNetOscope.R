SyntheticNet<-function(filename){
	
	library(mclust)

  library(igraph)
  print(filename)
  A=read.csv(filename, header = TRUE, sep = ",", quote = "")
  A=as.matrix(A)
  A=A[,c(2:dim(A)[2])]
  B=matrix(0,dim(A)[1],dim(A)[2])
  B[which(A=="True")]=1

colnames(B)=c(1:dim(B)[2])
rownames(B)=colnames(B)
#This is now a 01 network matrix, and we can build a graph
VV=which(B==1,arr.ind = TRUE)


G=graph_from_data_frame(VV,directed=FALSE)
G=igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE,
  edge.attr.comb = igraph_opt("edge.attr.comb"))

is_simple(G)

wt=walktrap.community(G)
fg=fastgreedy.community(G)
im=infomap.community(G)
lp=label.propagation.community(G)

True=c(1*rep(1,60),2*rep(1,60),3*rep(1,60),0*rep(1,820))
names(True)=c(1:1000)
################################################################################
randindex=c(0,1,3)
names(randindex)=c('WALKTRAP','FASTGREEDY','INFOMAP')

print('1-WALKTRAP')
Pred=membership(wt)
iless=which(sizes(wt)<9)
zeroclass=NULL
if(is.atomic(iless)==FALSE|length(iless)>=1){
for(i in c(1:length(iless))){
	zeroclass=c(zeroclass,wt[[iless[i]]])
    }
  }
zeroclass=as.integer(zeroclass)
connected=as.integer(names(Pred))
labelPred=rep(0,1,1000)
labelPred[connected]=Pred
labelPred[zeroclass]=0
print(table(labelPred,True))
randindex[1]=adjustedRandIndex(labelPred,True)
print(randindex[1])


print('2-FASTGREEDY')
Pred=membership(fg)
iless=which(sizes(fg)<9)
zeroclass=NULL
if(is.atomic(iless)==FALSE|length(iless)>=1){
for(i in c(1:length(iless))){
	zeroclass=c(zeroclass,fg[[iless[i]]])
    }
   }
zeroclass=as.integer(zeroclass)
connected=as.integer(names(Pred))
labelPred=rep(0,1,1000)
labelPred[connected]=Pred
labelPred[zeroclass]=0
print(table(labelPred,True))
randindex[2]=adjustedRandIndex(labelPred,True)
print(randindex[2])


print('3-INFOMAP')
Pred=membership(im)
iless=which(sizes(im)<9)
zeroclass=NULL
if(is.atomic(iless)==FALSE|length(iless)>=1){
for(i in c(1:length(iless))){
	zeroclass=c(zeroclass,im[[iless[i]]])
    }
  }
zeroclass=as.integer(zeroclass)
connected=as.integer(names(Pred))
labelPred=rep(0,1,1000)
labelPred[connected]=Pred
labelPred[zeroclass]=0
print(table(labelPred,True))
randindex[3]=adjustedRandIndex(labelPred,True)
print(randindex[3])

return(randindex)
	
}