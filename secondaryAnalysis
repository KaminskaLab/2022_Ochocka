library(tidyverse)
library(HDCytoData)
library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(scater)
library(plyr)

setwd("../../cyTOF/merged/out/clustering_B1plusB2")

fcsFiles <- list.files(path = ".", pattern = ".fcs$", full = TRUE, ignore.case = TRUE)

fcsFiles <- fcsFiles[which(substr(fcsFiles,5,9) %in% lol)] 


sceFiles <- mclapply(fcsFiles, function(x) {
  fs <- read.flowSet(x,
    transformation = FALSE, truncate_max_range = FALSE,
    column.pattern = "Y89Di|Sm147Di|Nd146Di|Nd150Di|Nd143Di|Tm169Di|Pr141Di|Gd160Di|Bi209Di|Nd144Di|Lu175Di|Eu151Di|Bi209Di|Dy163Di|Dy164Di|Er167Di|Er168Di|Er170Di|Eu151Di|Eu153Di|Gd155Di|Gd156Di|Gd158Di|Gd160Di|Ho165Di|Lu175Di|Nd142Di|Nd143Di|Nd144Di|Nd145Di|Nd146Di|Nd150Di|Pr141Di|Sm147Di|Sm149Di|Sm152Di|Tb159Di|Tm169Di|Y89Di|Yb171Di|Yb172Di|Yb173Di|Yb174Di|Yb176Di", invert.pattern = FALSE
  )


3>1.5
  names <- rownames(fs@phenoData)
  for (i in names) {
    fs@frames[[i]]@exprs <- fs@frames[[i]]@exprs[
      which(fs@frames[[i]]@exprs[, 27] > 0),
    ]
  }
  names <- rownames(fs@phenoData) 
  for (i in names) {
    fs@frames[[i]]@exprs <- fs@frames[[i]]@exprs[
      which(fs@frames[[i]]@exprs[, 1] > 1.5),
    ]
  }
  names <- rownames(fs@phenoData)
  for (i in names) {
    fs@frames[[i]]@exprs <- fs@frames[[i]]@exprs[
      which(fs@frames[[i]]@exprs[, 6] < 1.5),
    ]
  }

  names <- rownames(fs@phenoData)
  for (i in names) {
    fs@frames[[i]]@exprs <- fs@frames[[i]]@exprs[sample
    (nrow(fs@frames[[i]]@exprs), 80000), ]
  }

  write.flowSet(fs)
})

fcsFiles <- list.files(path = ".", pattern = ".*fcs$", full = TRUE, ignore.case = TRUE)


# nie wczytujemy pr141(cd49d i cd5 bo nie sa w obu batchach)  i Er170 - CD3 (bo usunelismy wysokie wartosci, a potem zostaly 0-0,3 i przeskalowane na 0-1 bardzo mocno wplywaja na klastrowaniae, a nie powinny)
fs <- read.flowSet(fcsFiles,
                   transformation = FALSE, truncate_max_range = FALSE,
                   column.pattern = "Er166Di|Y89Di|Sm147Di|Nd146Di|Nd150Di|Nd143Di|Tm169Di|Gd160Di|Bi209Di|Nd144Di|Lu175Di|Eu151Di|Bi209Di|Dy163Di|Dy164Di|Er167Di|Er168Di|Er170Di|Eu151Di|Eu153Di|Gd155Di|Gd156Di|Gd158Di|Gd160Di|Ho165Di|Lu175Di|Nd143Di|Nd144Di|Nd145Di|Nd146Di|Nd150Di|Sm147Di|Sm149Di|Sm152Di|Tb159Di|Tm169Di|Y89Di|Yb171Di|Yb172Di|Yb173Di|Yb174Di|Yb176Di", invert.pattern = FALSE
)


panel <- pData(parameters(fs[[1]]))
panel[4,2] <- "166Er_CADM1"

md <- read.csv2("../../../clinical_data.csv", na.strings = c("", "N/A"), sep = "\t")

md <- md[, c(1, 2, 3, 5, 6)]

md <- md[which(md$ID %in% substr(fcsFiles,13,17)),] #wybieramy tylko pacjentow GBM

md$tumor_type <- md$Diagnosis

md[md$Diagnosis %in% c(
     "Endometrium carcinoma BrM", "Squamous cell carcinoma BrM",
     "Melanoma BrM", "NSCLC BrM"
), "tumor_type"] <- "BrM"
md[md$Diagnosis %in% c("Glioblastoma", "Glioblastoma recurrent"), "tumor_type"] <- "GBM"
md[md$Diagnosis %in% c(
     "Anaplastic astrocytoma", "Anaplastic ependymoma", "Anaplastic oligodendroglioma",
     "Diffuse astrocytoma"
), "tumor_type"] <- "Glioma IDHwt"

md$Sex <- factor(md$Sex, levels = c("m", "f"))
md$file_name <- substr(fcsFiles, 3, 30)


sce <- prepData(fs, panel, md,
                cofactor = 5,
                transform = TRUE,
                panel_cols = list(channel = "name", antigen = "desc"),
                md_cols = list(
                     file = "file_name", id = "ID",
                     factors = c("tumor_type", "Sex", "Age")
                )
)

# zmieniamy wyniki >99percentyla na 99percentyl
for (i in colnames(macExpAsinh)) {
     q99 <- quantile(macExpAsinh[, i], probs = 0.99)
     macExpAsinh[(which(macExpAsinh[, i] > q99)), i] <- q99
}

macExpAsinh <- t(macExpAsinh)
colnames(macExpAsinh) <- NULL
dim(macExpAsinh) == dim(counts(sce))
colnames(macExpAsinh) == colnames(counts(sce))
rownames(macExpAsinh) == rownames(counts(sce))
assay(sce, "exprs") <- macExpAsinh

#### normalizacja 0-1 #####
macExp <- assay(sce, i = "exprs")
macExp <- t(macExp)

normalize <- function(x) {
     return((x - min(x)) / (max(x) - min(x)))
}
for (x in 1:length(colnames(macExp))) {
     macExp[, x] <- normalize(macExp[, x])
}
macExpAsinh <- macExp
colnames(macExpAsinh) <- names(sce)


#### koniec:normalizacja 0-1 #####


set.seed(1234)
sceCluster <- cluster(x = sce, features = rownames(sce), maxK = 20, seed = 1234)
#
set.seed(1)
sceClusterUMAP <- runDR(sceCluster, "UMAP", cells = 10000, features = rownames(sceCluster))
 
plotDR(sceClusterUMAP, "UMAP", color_by = "meta12", scale = FALSE)

markers <- rownames(sceClusterUMAP)


exp_plots <- mclapply(markers, function(i){
  plotDR(sceClusterUMAP, "UMAP", color_by = i, scale = FALSE)  +
    ggtitle(i)+
    theme_classic()+
    theme(  title = element_text(size=40),
            axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title = element_text(hjust=0), legend.title = element_blank(),
            legend.position = "right", legend.text=element_text(size=14),
            axis.title.x=element_text(size=16),
            axis.title.y=element_text(size=16))
},mc.cores=31)

library(ggpubr)
plotName <- paste("sceClusterUMAP",".png",sep="")
print(plotName)
png(plotName,
    width=4800, height=2700)
print(ggarrange(plotlist = exp_plots, ncol =7, nrow=5))
dev.off()

macExpr <- t(assay(sceClusterUMAP, i="exprs"))
metadata <- sceClusterUMAP@colData
redDim <- reducedDim(sceClusterUMAP)
metadata <- cbind(metadata,redDim)
Tmac <- cbind(metadata,macExpr)
colnames(Tmac)[6] <- "UMAP1"
colnames(Tmac)[7] <- "UMAP2"


library(Rphenograph)
library(gridExtra)
library(ggpubr)

macExp <- Tmac[8:38]
macExp <- as.data.frame(macExp)
colnames(macExp) <- substr(colnames(macExp),2,50)
phenograph_out <- Rphenograph(macExp, k=70)


modularity(phenograph_out[[2]])
membership(phenograph_out[[2]])

macExp$phenograph_cluster <- factor(membership(phenograph_out[[2]]))
macExp <- cbind(Tmac$sample_id,Tmac$Sex,Tmac$tumor_type,Tmac$UMAP1,Tmac$UMAP2,macExp)
colnames(macExp)[[1]] <- "sample_id"
colnames(macExp)[[2]] <- "Sex"
colnames(macExp)[[3]] <- "tumor_type"
colnames(macExp)[[4]] <- "UMAP1"
colnames(macExp)[[5]] <- "UMAP2"


saveRDS(macExp,"myeloMacExpKmeans_minust141Pr_170EU.rds")


umapClusters <- ggplot(macExp, aes(x=UMAP1, y=UMAP2, col=phenograph_cluster)) + geom_point(size = 0.1) +
  ggtitle("kMeans full 1 batch")+
  theme_classic()+
  guides(color = guide_legend(override.aes = list(size=20)))+
  theme(  title = element_text(size=40),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_text(hjust=0), legend.title = element_blank(),
          legend.position = "right", legend.text=element_text(size=30),
          #legend.key.size=unit(30,"line"),
          axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16))

png("kMeansClusters.png",
    width=1500, height=1500)
print(umapClusters)
dev.off()




##################### SZUKAMY ROZNIC POMIEDZY M/F ########################


library(tidyverse)
library(plyr)
library(ggplot2)
library(ggpubr)
library(viridis)

#CD56 - NK
#CD5  - T-cells
setwd("~/cyTOF/merged/out/clustering_B1plusB2")
macExp <- readRDS("B1plusB2ExprsMatrixesKmeans_minust141Pr.rds")

data <- macExp
metadata<-md

#patient occuring in both panels ZH683

#Glioblastoma patients

# ID Sex              Diagnosis Myeloid.focused.panel.CyTOF
# 7  ZH678   m           Glioblastoma                           1
# 9  ZH720   m           Glioblastoma                           1
# 10 ZH736   f           Glioblastoma                           1
# 11 ZH746   m           Glioblastoma                           2
# 14 ZH761   m           Glioblastoma                           2
# 15 ZH784   m Glioblastoma recurrent                           2
# 16 ZH791   m           Glioblastoma                           2
# 17 ZH794   f           Glioblastoma                           2
# 18 ZH802   f           Glioblastoma                           2
# 19 ZH810   f           Glioblastoma                           2
# 20 ZH816   f           Glioblastoma                           2
# 21 ZH818   m           Glioblastoma                           2
# 23 ZH813   f           Glioblastoma                           2
# 24 ZH961   f Glioblastoma recurrent                           1




clustersUMAP <- ggplot(full, aes(x=UMAP1, y=UMAP2, color=phenograph_cluster))+
     geom_jitter(size=0.7, alpha=0.5)+
     geom_text(aes(label=phenograph_cluster), size=4, alpha=0.5, check_overlap=T, color="black")+
     #facet_grid(.~day)+
     guides(colour = guide_legend(override.aes = list(size=7)))

png(file = "clustersUMAP_minus141.png", width = 1500, height = 1500)
clustersUMAP
dev.off()



#expresion plots
theme_feature<-theme_classic()+
     theme(axis.text = element_blank(), axis.ticks = element_blank(),
           axis.title = element_text(hjust=0), legend.title = element_blank(),
           legend.position = "right", legend.text=element_text(size=14),
           title = element_text(size=40), axis.title.x=element_text(size=16),
           axis.title.y=element_text(size=16))
YlOrRd<-colorRampPalette(c("Yellow","Orange", "Red"))

genes<-colnames(full)

graphs<-list()
graphs <- mclapply(7:37, function(i) {
     max <- round(max(full[, genes[i]]), digits = 1)
     min <- 0 # ceiling(min(panel_99[,genes[i]]))
     plot <- ggplot(full, aes(x = UMAP1, y = UMAP2)) +
          geom_jitter(size = 1, alpha = 0.5, aes_string(color = as.name(genes[i]))) +
          labs(title = genes[i]) +
          scale_color_gradientn(colors = (c(
               viridis(50),
               YlOrRd(50)
          ))) +
          # breaks=c(min,max))+
          # , labels=c("min", "max"))+
          xlab("UMAP_1") +
          ylab("UMAP_2") +
          theme_feature +
          theme(legend.text = element_text(size = 20), legend.position = "right")
     # facet_grid(sex~day)
     return(plot)
}, mc.cores = 31)
png(file = "expression_minus141.png", width = 4000, height = 2000)
ggarrange(plotlist = graphs, ncol = 8, nrow = 4, common.legend = F)
dev.off()


#annotate clusters
#ZH683_1
anno<-data.frame(        "cluster"=c(   5,   11,    8, 1, 2, 14, 15, 9,  10, 3, 13),
                         "cell_annotation"=c("NK", "NK", "Mo",  rep("Mphi", 5),rep("MG", 3)))

#ZH683_2
anno<-data.frame(         "cluster"=c(   3, 9, 16,   10, 11, 13, 5,4,14,    8, 6, 2),
                          "cell_annotation"=c(rep("NK",3), "Mo", rep("Mphi", 5), rep("MG",3)))

a<-length(unique(patient$phenograph_cluster))
anno<-rbind(anno,
            data.frame("cluster"=c(1:a)[!c(1:a) %in% anno$cluster],
                       "cell_annotation"=rep("UN", length(c(1:a)[!c(1:a) %in% anno$cluster]))))

patient$cell_annotation<-anno[match(patient$phenograph_cluster, anno$cluster),2]
table(patient$cell_annotation)

#veryfy annotation on a plot
png(paste("plots/Cell_annotation", patient_id), width=500, height=500)
ggplot(patient, aes(x=UMAP1, y=UMAP2, color=cell_annotation))+
  geom_jitter(size=0.7, alpha=0.5)+
  geom_text(aes(label=phenograph_cluster), size=4, alpha=0.5, check_overlap=T, color="black")+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  theme_bw()
dev.off()

#save new patient list with annotated clusters
#patients_anno<-list()

a<-list(patient)
names(a)<-patient_id
patients_anno<-c(patients_anno, a)

saveRDS(patients_anno, "data/patients_anno.RDS")



#____________________________________________________________________________________
#append cell type distribution

group<-c("MG", "Mo", "Mphi")
patient<-patient[patient$cell_annotation %in% group,]
#cell_type_distribution<-list()
b<-table(patient$cell_annotation) %>% as.vector() %>% list
names(b)<-patient_id
cell_type_distribution<-c(cell_type_distribution, b)

saveRDS(cell_type_distribution, "data/cell_type_distribution.RDS")

colnames(cell_type_distribution)<-c("patient_id", group)


#____________________________________________________________________________________
#extract MHC level

#MHCII<-list()
c<-c(aggregate(patient[,"150Nd_HLA-DR"], list(patient$cell_annotation), mean)$x,
     mean(patient[,"150Nd_HLA-DR"]))%>% list()
names(c)<-patient_id
MHCII<-c(MHCII, c)

saveRDS(MHCII, "data/MHCII.RDS")


colnames<-c(group, "total")



mikroglej <- c(11,2,1)
monocyty <- c(8,23) 
makrofagi <- c(10,3,9,4,26)


females <- nmacExp[which(nmacExp$Sex == "f"),]
males <- nmacExp[which(nmacExp$Sex == "m"),]

malesMikro <- males[which(males$phenograph_cluster == 11 |males$phenograph_cluster == 2 |males$phenograph_cluster == 1),]
malesMakro <- males[which(males$phenograph_cluster == 10 |males$phenograph_cluster == 3 |males$phenograph_cluster == 9|males$phenograph_cluster == 26|males$phenograph_cluster == 4),]
malesMono <- males[which(males$phenograph_cluster == 8 |males$phenograph_cluster == 23),]

valuesPie <- c(dim(malesMikro)[[1]],dim(malesMakro)[[1]],dim(malesMono)[[1]])
labelsPie <- c("Mikroglej","Makrofagi","Monocyty")
pct <- round(valuesPie/sum(valuesPie)*100)
labelsPie <- paste(labelsPie,pct,"%")

pie(valuesPie, labels=labelsPie, main = "Male")




femalesMikro <- females[which(females$phenograph_cluster == 11 |females$phenograph_cluster == 2 |females$phenograph_cluster == 1),]
femalesMakro <- females[which(females$phenograph_cluster == 10 |females$phenograph_cluster == 3 |females$phenograph_cluster == 9|females$phenograph_cluster == 26|females$phenograph_cluster == 4),]
femalesMono <- females[which(females$phenograph_cluster == 8 |females$phenograph_cluster == 23),]


valuesPie <- c(dim(femalesMikro)[[1]],dim(femalesMakro)[[1]],dim(femalesMono)[[1]])
labelsPie <- c("Mikroglej","Makrofagi","Monocyty")
pct <- round(valuesPie/sum(valuesPie)*100)
labelsPie <- paste(labelsPie,pct,"%")

pie(valuesPie, labels=labelsPie, main = "Female")




#annotate clusters variant II

cell_types<-list(
     NK=c(2, 12, 14, 17, 25), #CD56
     Granulocytes=c(15,19), #CD66b
     cDC=c(8), #CD1c
     TAM=c(1,4,6, 11,22,23), #CCR2-, CD49d, CD45, CD68, CD14
     Mo=c(3,20,21), #CCR2+, CD141
     BAMs= c(5,7), #CD206, CD169
     "Plasma cells"=c(18),
     UN=c(9,10,13,16,20,24))

#annotate clusters variant III
cell_types<-list(
     NK=c(2, 12, 14, 17, 25), #CD56
     Granulocytes=c(15,19), #CD66b
     cDC=c(8), #CD1c
     "TAM HLA-DR low"=c(1, 4, 11),
     "TAM HLA-DR high"=c(6,22,23), #CCR2-, CD49d, CD45, CD68, CD14
     Mo=c(3,20,21), #CCR2+, CD141
     BAMs= c(5,7), #CD206, CD169
     "Plasma cells"=c(18),
     UN=c(9,10,13,16,20,24))



male <- full[which(full$Sex=="m"),]
female <- full[which(full$Sex=="f"),]





















