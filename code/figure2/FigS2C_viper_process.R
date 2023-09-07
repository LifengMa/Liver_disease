setwd("/media/ggj/ggjlab2/hezuo/gwh/cancer/Scissor/")
library(Seurat)
load("./Scissor_HCC_survival_0.25_result.RData")

Idents(HCC) <- HCC$scissor
DimPlot(HCC,reduction="tsne",label=T)

Scissor<-subset(HCC,idents=c(1,2))


table(Scissor$scissor)
Scissor$scissor <- gsub("1","Scissor+",Scissor$scissor)
Scissor$scissor <- gsub("2","Scissor-",Scissor$scissor)
Idents(Scissor) <- Scissor$scissor
Scissor_pos<-subset(Scissor,idents=c("Scissor+"))
Scissor_neg<-subset(Scissor,idents=c("Scissor-"))
DimPlot(Scissor_pos)
dir.create("viper")
setwd("./viper/")
save(Scissor,Scissor_pos,Scissor_neg,file = "rawdata.RData")

###pseducell
DimPlot(Scissor,label=T)
Inter.id<-as.data.frame(Scissor@active.ident)
Inter.id$Cell_id<-rownames(Inter.id)
colnames(Inter.id)<-cbind("Celltype","Cell_id")
Inter<-as.data.frame(as.matrix(GetAssayData(Scissor@assays$RNA,slot = "counts")))
#rownames(Inter)<-Gene.name[,1]
#colnames(Inter)<-Cell.name[,1]
#Inter.id<-cbind(as.character(id$Cell_id),as.character(id$Group))
#Inter.id<-as.data.frame(Inter.id)
#colnames(Inter.id)<-cbind("Cell_id","Celltype")
#rownames(Inter.id)<-Inter.id$Cell_id
nn<-intersect(colnames(Inter),Inter.id$Cell_id)
Inter1<-Inter[,as.character(nn)]
Inter<-Inter1
Inter.id1<-Inter.id[as.character(nn),]
Inter.id<-Inter.id1
raw.10e6<-t(t(Inter)/colSums(Inter))*1000000
Inter<-raw.10e6
Inter<-as.data.frame(Inter)


pseudocell.size = 20 ## 10 test
new_ids_list = list()
for (i in 1:length(levels(Inter.id$Celltype))) {
  cluster_id = levels(Inter.id$Celltype)[i]
  cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)		
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)	
  new_ids_list[[i]] <- pseudo_ids		
}

new_ids <- unlist(new_ids_list)
new_ids <- as.data.frame(new_ids)
new_ids_length <- table(new_ids)

new_colnames <- rownames(new_ids)  ###add
all.data<-Inter[,as.character(new_colnames)] ###add
all.data <- t(all.data)###add

new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                    list(name=new_ids[,1]),FUN=mean)
rownames(new.data)<-new.data$name
new.data<-new.data[,-1]

new_ids_length<-as.matrix(new_ids_length)##
short<-which(new_ids_length<10)##
new_good_ids<-as.matrix(new_ids_length[-short,])##
new_good_ids<-as.matrix(new_ids_length)
result<-t(new.data)[,rownames(new_good_ids)]



#Human------------------------------------------------------------
#colnames(result)<-paste("CN",colnames(result),sep="_")
rownames(result)<-rownames(Inter)
saveRDS(result,file="pse20.rds") ###
cc<-gsub("[_]Cell.*$","",colnames(result))
new.phe<-cbind(colnames(result),'Scissor',cc)
colnames(new.phe)<-c("Sample_ID","Study_ID","Celltype")
write.table(new.phe,file="pse20.pheotype",quote=F,row.names=F,sep=",") ###

posname<-paste0("Scissor+_Cell",c(0:131))
negname<-paste0("Scissor-_Cell",c(0:131))
pos_matrix<-result[,colnames(result)%in%posname]
neg_matrix<-result[,colnames(result)%in%negname]
write.table(pos_matrix,quote = F,file="pos_matrix.txt",col.names = T,row.names = F,sep="\t")
write.table(neg_matrix,quote = F,file="neg_matrix.txt",col.names = T,row.names = F,sep="\t")
