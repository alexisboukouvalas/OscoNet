
filename1='../OscopeData/dataAdjacencyMatrixe02/WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s123res_C73G9200_edgeNetwork.csv'
filename2='../OscopeData/dataAdjacencyMatrixePiDiv9/WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s123res_C73G9200_edgeNetwork.csv'
filename3='../OscopeData/FirstData/WaterfallbootstrapSequential_CPUDPSF_TF_2000Boot_data_s123res_C73G9200_edgeNetwork.csv'
data1=read.csv(file=filename1, header= TRUE,row.names=1)
data2=read.csv(file=filename2, header= TRUE,row.names=1)
data3=read.csv(file=filename3, header= TRUE,row.names=1)


G1=graph_from_data_frame(data1,directed=FALSE)
G2=graph_from_data_frame(data2,directed=FALSE)
G3=graph_from_data_frame(data3,directed=FALSE)

g1name=get.vertex.attribute(G1)$name
g2name=get.vertex.attribute(G2)$name
g3name=get.vertex.attribute(G3)$name

g1g2=intersect(g1name, g2name)
g1g3=intersect(g1name, g3name)
g2g3=intersect(g2name, g3name)
g1g2g3=intersect(g1name,g2g3)
n12=length(g1g2)
n13=length(g1g3)
n23=length(g2g3)
n123=length(g1g2g3)

print(n12)
print(n13)
print(n23)
print(n123)
n1=length(g1name)
n2=length(g2name)
n3=length(g3name)

print(n1)
print(n2)
print(n3)

grid.newpage()
draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12 = n12, n23 =23, n13 = n13, 
    n123 = n123, category = c("quiescent (a)", "quiescent (b)", "quiescent(c)"), lty = "blank", 
    fill = c("skyblue", "pink1", "mediumorchid"))