library(ModelMetrics)
library(randomForest)
library(genefilter)
library(hgu133plus2frmavecs)
library(hgu133plus2.db)
library(hgu133a.db)
library(hgu133afrmavecs)
library(annotate)

Brain = read.csv("Brain_GSE50161.csv",check.names=F,row.names = T)
rownames(Brain) = Brain[,1]
Brain = Brain[,-1]

### use probe with largest IQR to represent the gene

exp = na.omit(t(Brain))[-1,]

arrayIQR<-apply(exp,1,IQR)
probe<-rownames(exp)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,"hgu133plus2.db")

exp2<-exp[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu133plus2.db")

rownames(exp2)<-geneSymbol
write.csv(exp2, "Data/Brain_new.csv")

## Leukemia

Leukemia = read.csv("Leukemia_GSE28497.csv",check.names=F,row.names = 1)

### use probe with largest IQR to represent the gene

exp = t(Leukemia[,-1])

arrayIQR<-apply(exp,1,IQR)
probe<-rownames(exp)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,"hgu133a.db")

exp2<-exp[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu133a.db")

rownames(exp2)<-geneSymbol

Leukemia_new = data.frame(type = Leukemia$type, t(exp2))
write.csv(Leukemia_new, "Data/Leukemia_new.csv")

