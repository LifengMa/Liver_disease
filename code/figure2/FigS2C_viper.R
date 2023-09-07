## Figure S2C
setwd("/media/ggj/ggjlab2/hezuo/gwh/cancer/Scissor/viper/")
library(viper)
library(bcellViper)
library(mixtools)
library(Seurat)

load("./rawdata.RData")


Scissor<-readRDS("pse20.rds")
pbmc<-CreateSeuratObject(Scissor)
pbmc<-NormalizeData(pbmc)
Scissor_matrix_normalize<-GetAssayData(pbmc@assays$RNA,slot="data")
Scissor_matrix_normalize<-as.data.frame(as.matrix(Scissor_matrix_normalize))
gene<-as.data.frame(rownames(Scissor_matrix_normalize))
colnames(gene)<-"gene"
Scissor_matrix_normalize<-cbind(gene,Scissor_matrix_normalize)
colnames(Scissor_matrix_normalize) <- gsub("\\+","_pos",colnames(Scissor_matrix_normalize))
colnames(Scissor_matrix_normalize) <- gsub("\\-","_neg",colnames(Scissor_matrix_normalize))
write.table(Scissor_matrix_normalize,quote = F,file="./Scissor_matrix_normalize.txt",col.names = T,row.names = F,sep="\t")

#####################linux########################
#step1
## calculate threshold with a fixed seed
cd /media/ggj/ggjlab2/hezuo/gwh/cancer/Scissor/viper/
java -Xmx5G -jar /home/ggj/Documents/ARACNe-AP/dist/aracne.jar -e Scissor_matrix_normalize.txt -o outputFolder --tfs tfs.txt --pvalue 1E-8 --seed 1 --calculateThreshold 

#run ARACNe on a single bootstrap
java -Xmx5G -jar /home/ggj/Documents/ARACNe-AP/dist/aracne.jar -e Scissor_matrix_normalize.txt -o outputFolder --tfs tfs.txt --pvalue 1E-8 --seed 1 

#step2 viper
library(viper)
library(bcellViper)
library(mixtools)
Scissor<-read.csv("./Scissor_matrix_normalize.txt",header = T,sep='\t')
rownames(Scissor)<-Scissor[,1]
Scissor<-Scissor[,-1]
Scissor[1:5,1:5]
Scissor<-as.matrix(Scissor)
######rm network.txt header
regul<-aracne2regulon("./outputFolder/bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt",Scissor,format = c("3col"))
print(regul)

pheno<-read.table("pse20.pheotype",sep=",",header=T)
head(pheno)
pos<-as.character(pheno[pheno$Celltype%in%'Scissor+',]$Sample_ID)
neg<-as.character(pheno[pheno$Celltype%in%'Scissor-',]$Sample_ID)
pos <- gsub("\\+","_pos",pos)
neg <- gsub("\\-","_neg",neg)
matrix2<-Scissor[,pos]
matrix1<-Scissor[,neg]
dim(matrix1)
dim(matrix2)

signature<-rowTtest(matrix1,matrix2)
signature<-(qnorm(signature$p.value/2,lower.tail = FALSE)*sign(signature$statistic))[,1]

nullmodel<-ttestNull(matrix2,matrix1,per=1000,repos=TRUE,verbose=FALSE)

mrs<-msviper(signature,regul,nullmodel,verbose=FALSE)
summary(mrs,40)
plot(mrs,cex=0.7,30)
pdf("./viper.pdf",height = 15,width = 10)
plot(mrs,30)
dev.off()
save(mrs,signature,regul,nullmodel,matrix1,matrix2,Scissor,file="./viper.RData")

