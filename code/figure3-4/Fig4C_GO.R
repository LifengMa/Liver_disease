#Fig4C
setwd("/media/ggj/ggjlab2/hezuo/gwh/myeloid/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("./myeloid_harmony_1129.RData")


# GO_cluster --------------------------------------------------------------

library(gprofiler2)
library(openxlsx)

#library(BiocManager)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

cluster <- pbmc.harmony.markers[pbmc.harmony.markers$cluster==9,]
cluster.use <- cluster[cluster$p_val_adj<0.05&cluster$avg_log2FC>0.25,]
symbol=as.character(cluster.use$gene)
eg = bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
head(id)
length(gene)
length(id)

##GO
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
View(ego@result)
##FABP5+
use<- ego@result[c("GO:0050900","GO:0006869","GO:0001666","GO:0045765","GO:0019882"),]
##Kupffer
use<- ego@result[c("GO:0048002","GO:0002697","GO:0042110","GO:0019882","GO:0001819"),]
##CLU+
use<- ego@result[c("GO:0046686","GO:0006979","GO:0071276","GO:0097193","GO:1903706"),]
##APOE+
use<- ego@result[c("GO:0002443","GO:0019882","GO:0002399","GO:0002377","GO:0050863"),]

head(use)
use$p.adjust_log <- -log10(use$p.adjust)


#use<-use[order(aaa$p.value ),]
#dotplot(ego)
#library(topGO)
#plotGOgraph(ego)
library(ggpubr)
ggdotchart(use,x="Description",y="p.adjust_log",add="segments",color = "#F7DC68",add.params = list(color="lightgray",size=1),
           dot.size = 5,sorting = "descending",rotate = T,
           ggtheme = theme_pubclean())+
           ylab("-log10(p.adjust)")+xlab("GO term")+ggtitle("Macrophage FABP5+")
ggsave("./myeloid/GO_fabp5.pdf",width = 6,height = 4)

ggdotchart(use,x="Description",y="p.adjust_log",add="segments",color = "#2E9599",add.params = list(color="lightgray",size=1),
           dot.size = 5,sorting = "descending",rotate = T,
           ggtheme = theme_pubclean())+
  ylab("-log10(p.adjust)")+xlab("GO term")+ggtitle("Kupffer cell")
ggsave("./myeloid/GO_kupffer.pdf",width = 6,height = 4)

ggdotchart(use,x="Description",y="p.adjust_log",add="segments",color = "#F46C3F",add.params = list(color="lightgray",size=1),
           dot.size = 5,sorting = "descending",rotate = T,
           ggtheme = theme_pubclean())+
  ylab("-log10(p.adjust)")+xlab("GO term")+ggtitle("Macrophage CLU+")
ggsave("./myeloid/GO_clu.pdf",width = 6,height = 4)

ggdotchart(use,x="Description",y="p.adjust_log",add="segments",color = "#A7226F",add.params = list(color="lightgray",size=1),
           dot.size = 5,sorting = "descending",rotate = T,
           ggtheme = theme_pubclean())+
  ylab("-log10(p.adjust)")+xlab("GO term")+ggtitle("Macrophage")
ggsave("./myeloid/GO_macrophage.pdf",width = 7,height = 4)
# ggplot(data=use)+
#   geom_bar(aes(x=reorder(Description,c(Count)),y=Count, fill=-log10(p.adjust)), stat='identity') + 
#   coord_flip() +
#   scale_fill_gradient(expression(-log10(p.adjust)),low="#81c5f4",high = "firebrick3") +
#   xlab("") +
#   ylab("Gene count") +
#   scale_y_continuous(expand=c(0, 0))+
#   theme(
#     axis.text.x=element_text(color="black",size=rel(0.8)),
#     axis.text.y=element_text(color="black", size=rel(0.8)),
#     axis.title.x = element_text(color="black", size=rel(0.8)),
#     legend.text=element_text(color="black",size=rel(0.7)),
#     legend.title = element_text(color="black",size=rel(0.7))
#   )

