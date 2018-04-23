load("ZerosReplaceANDlogTransformedDatasetTaxLevelL3.RData")
load("vectorOfComb.RData")

library(energy)
library(zCompositions)
#### hierarchical clustering

EHCres=list()
mat<- matrix(NaN, nrow = 4, ncol=4, byrow = TRUE, dimnames = list(c("sq","bl","gbm","czm") ,c("tss","cenLR","addLR","ilr")))
#rownames(mat)<-vectorOfComb

for(combi in vectorOfComb ){


subsetIdx.1=grep("\\.S1\\.1$",rownames(combiObj[[combi]]))

subsetIdx.2=grep("\\.S2\\.1$",rownames(combiObj[[combi]]))

subsetIdx.3=grep("\\.S3\\.1$",rownames(combiObj[[combi]]))

subsetIdx.4=grep("\\.S1\\.2$",rownames(combiObj[[combi]]))

subsetIdx.5=grep("\\.S2\\.2$",rownames(combiObj[[combi]]))

subsetIdx.6=grep("\\.S3\\.2$",rownames(combiObj[[combi]]))


df2prcomp=rbind(

combiObj[[combi]][subsetIdx.1,],
combiObj[[combi]][subsetIdx.2,],
combiObj[[combi]][subsetIdx.3,],
combiObj[[combi]][subsetIdx.4,],
combiObj[[combi]][subsetIdx.5,],
combiObj[[combi]][subsetIdx.6,]

)

subsetIdx=c(subsetIdx.1,subsetIdx.2,subsetIdx.3,subsetIdx.4,subsetIdx.5,subsetIdx.6)

EHCres[[combi]][[1]]<-energy.hclust(dist(combiObj[[combi]][subsetIdx,]))

# the ground truth : 16 individuals each associated with 3 different sample collection method and 2 isolation kits; 16*3*2=96 

g <- setNames(c(1:96), rownames(combiObj[[combi]])[subsetIdx]) 
for(i in 1:16){
#cat(i,"\n")
togrep=substr(names(cutree(EHCres[[combi]][[1]], k=16))[i],2,4)
mysubset=grep(togrep,names(cutree(EHCres[[combi]][[1]], k=16)))
g[mysubset]<-i
}

par(cex=0.7, mar=c(15, 8, 4, 1))
plot(EHCres[[combi]][[1]],hang=-1)
EHCres[[combi]][[2]]<-table(cutree(EHCres[[combi]][[1]], k=16) == g)

cat(combi,"\n")

show(EHCres[[combi]][[2]])
}
