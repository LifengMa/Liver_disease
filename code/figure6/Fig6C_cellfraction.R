setwd("/media/ggj/ggjlab2/hezuo/gwh/endothelial/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("./endothelial_harmony.RData")
DimPlot(pbmc.harmony,label=T)
# cell fraction -----------------------------------------------------------
##ggplot2
#pbmc.harmony$type <- factor(pbmc.harmony$type,levels=c("CA","Fibrosis","AD","Normal"))
unique(pbmc.harmony$Batch2)
final<-data.frame(Idents(pbmc.harmony),pbmc.harmony$Batch2)
colnames(final)<-c("cluster","state")
#final$state <- gsub("S1|S2","",final$state)
unique(final$state)
final$state <- as.factor(final$state)
#table(pbmc.harmony$type)

cluster1<-final[final$cluster==0,]
table(cluster1$state)
clustername<-data.frame(table(cluster1$cluster))
statename<-data.frame(table(cluster1$state))
clusternum<-nrow(data.frame(table(final$cluster)))
state<-data.frame(matrix(nrow = nrow(statename),ncol = 0),row.names = statename$Var1)
for (i in 1: clusternum ) {
  temp<-final[final$cluster==(i-1),]
  state0<-data.frame(table(temp$state))
  state<-cbind(state,state0[,2])
}
colnames(state) <- paste0("cluster",c(0:(length(colnames(state))-1)))

#state$num <- apply(state,1,sum)
state <- apply(state,1,function(x){x/sum(x)})
state <- as.data.frame(t(state))
state$patient <- rownames(state)

state.use <- melt(state)
colnames(state.use) <- c("donor","cluster","fraction")
state.use$type <- state.use$donor
state.use$type <- gsub("[0-9]","",state.use$type)
state.use$type <- gsub("ICC","CA",state.use$type)
state.use$type <- gsub("S","",state.use$type)
state.use$type <- gsub("Liver","Normal",state.use$type)
state.use$type <- gsub("Fib","Cirhosis",state.use$type)
state.use$type <- gsub("CAAD","ADJ",state.use$type)
unique(state.use$type)
state.use$fraction <- state.use$fraction*100
source("/media/ggj/ggjlab2/hezuo/gwh/errorbar.R")
colnames(state.use)
head(state.use)
#df2 <- data_summary(state.use,varname = "fraction",groupnames = c("type","cluster"))

df2 <- state.use
library(ggplot2)
df2$type <- factor(df2$type,levels=c("CA","Cirhosis","ADJ","Normal"))
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(df2$cluster))))

ggplot(df2,aes(x=type,y=fraction,fill=cluster))+geom_bar(stat = "summary",fun="mean",position = position_dodge())+
  stat_summary(fun.data = 'mean_se',geom="errorbar",width=0.15,position = position_dodge(0.9))+
  #stat_compare_means(label="p.signif",method="wilcox.test",ref.group=".all.")+
  facet_grid(.~cluster)+ylab("Celluar fraction(%)")+scale_fill_manual(values=col_flg)+
  theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+NoLegend()


###Normal distribution test
library(dplyr)
shapiro.test(state.use$fraction)

#homogeneity test for variance
bartlett.test(state.use$fraction,state.use$type)

library(ggpubr)
library(gridExtra)
fig <- list()

for(i in 1:length(unique(df2$cluster))){
  temp <- df2[df2$cluster==unique(df2$cluster)[i],]
  #temp2 <- state.use[state.use$cluster==unique(df2$cluster)[i],]
  if(i==1){
    p <-  ggplot(temp,aes(x=type,y=fraction,fill=cluster))+geom_bar(stat = "summary",fun="mean",position = position_dodge())+
      stat_summary(fun.data = 'mean_se',geom="errorbar",width=0.15,position = position_dodge(0.9))+
      stat_compare_means(label="p.signif",method="wilcox.test",ref.group = 'Normal',hide.ns = TRUE)+
      #geom_hline(yintercept = mean(temp$fraction),linetype=2)+
      ylab("Celluar fraction(%)")+scale_fill_manual(values=col_flg[i])+
      ggtitle(unique(df2$cluster)[i])+xlab("")+
      theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+NoLegend()
  }
  else{
    p <-  ggplot(temp,aes(x=type,y=fraction,fill=cluster))+geom_bar(stat = "summary",fun="mean",position = position_dodge())+
      stat_summary(fun.data = 'mean_se',geom="errorbar",width=0.15,position = position_dodge(0.9))+
      stat_compare_means(label="p.signif",method="wilcox.test",ref.group="Normal",hide.ns = TRUE)+
      #geom_hline(yintercept = mean(temp$fraction),linetype=2)+
      ylab("")+scale_fill_manual(values=col_flg[i])+
      ggtitle(unique(df2$cluster)[i])+xlab("")+
      theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+NoLegend()
  }
  
  fig[[i]] <- p
}

fig[['nrow']] <- 1
fig[['ncol']] <- 9

pdf("./cellfraction.pdf",width = 28,height = 4)
do.call('grid.arrange', fig)
dev.off()
#facet_grid(.~cluster)
#compare_means(data=temp2, fraction~type,label="p.signif",method="wilcox.test",ref.group = ".all.")


