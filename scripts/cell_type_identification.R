library(Matrix)
library(tidyverse)
library(biomaRt)
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)

data_norm<-samplesIntegratedStandard@assays$integrated@data
data_cord<-samplesIntegratedStandard@reductions$umap@cell.embeddings
meta<-samplesIntegratedStandard@meta.data
prot_int<-samplesIntegratedStandard@assays$integrated_ADT@data
prot_raw<-samplesIntegratedStandard@assays$ADT@data

#generate ensembl annotations
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)

anno<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
            filters = "ensembl_gene_id", values = data_norm@Dimnames[1] , mart = ensembl)


#------------------select genes of interests-------------------------#
#Immunemarker panel
markers<-read.csv("data/supp_table1.csv", sep=",")[,1] 

#top expressed genes
top50<-getTopMarkersGnames(samplesIntegratedStandardMarkers_1.1, 50)
top_macro<-readRDS("markersList_res.1.1_integratedStandard_allMoMphi.RDS") %>%  rbindlist()

#Interferon pathway genes
IFN<-read.csv("data/INTERFEROM_allTypes_clear.csv", sep="\t")
IFN2<- anno[anno$external_gene_name %like% 'Ifi|Isg|Ifn|Irf',2]

#Dendritc cell markers
DC<-c("CD45_prot", "CD11b_prot", "Itgax", "MHCII_prot", "Cd8a", "Cd8b", "Cd4", 
      "Itgae", "Cd207", "Epcam", "Cd24a", "Btla","Kit", "Dpp4", "Xcr1",
      "Cd36", "Cst3", "Clec9a", "Cadm1", "Ly75", "Cx3cr1", "Cd209b", "Cd209d", "Cd209e", "Adgre1",
      "Sirpa", "Fcgr1", "Ly6C_prot", "Cd83", "Cd80")
markers<-c(markers, top50$Gene.name, IFN$Gene.Name, IFN2, DC) %>% unique

#extract genes of interest to count matrix
counts<-data_norm[markers_present$ensembl_gene_id,] %>% as.matrix() %>% as.data.frame() 

#convert ensembl id to external gene name
rownames(counts)<-anno[match(rownames(counts), anno$ensembl_gene_id), "external_gene_name"] 
counts<-t(counts) %>% as.data.frame()

#Extract protein counts
proteins1<-c( "CD11b-prot", "CD45-prot", "Ly6C-prot", "Gr-1-prot", "PD-L1-prot",
              "CD49d-prot", "MHCII-prot", "CD74-prot","Tmem119-prot", "CD52-prot")
prot<-prot_int[proteins1,]  %>% as.matrix() %>% t() %>% as.data.frame() 


#prepare data frame compatible with ggplot required format
panel<-data.frame(data_cord)
panel[,c(3:9)]<-meta[rownames(panel),c("day", "sex", "condition", "HTO_classification", "replicate", 
                                       "integrated_snn_res.0.5", "integrated_snn_res.1.1" )]
panel[,c(10:20)]<-prot[rownames(panel),]
panel[,c(21:(20+length(markers_present$ensembl_gene_id)))]<-counts[rownames(panel),]

#annotate cell types
cell_type<-read.csv("data/cell_types_1_1.csv", row.names = 1)
a<-panel[,c(1:9)]
a[,c(10:12)]<-cell_type[match(a$integrated_snn_res.1.1, cell_type$Cluster), c(2:4)]


panel<-cbind(a,panel[c(10:length(panel))])

colnames(panel)<-gsub('\\-', '_', colnames(panel))

write.csv(panel, "data/panel.csv")
panel<-read.csv("data/panel.csv", row.names = 1)

#trim extreme values to 99th percentile that would otherwise skew coloring scale
percentile_99<- apply(panel[,c(13:length(colnames(panel)))], 2, quantile, probs=0.99)
percentile_01<- apply(panel[,c(13:length(colnames(panel)))], 2, quantile, probs=0.01)
panel_99<-panel[,c(13:length(colnames(panel)))]
for (i in 1:length(colnames(panel_99))){
  panel_99[panel_99[,i] > percentile_99[i], i]<- percentile_99[i]
  panel_99[panel_99[,i] < percentile_01[i], i]<- percentile_01[i]
}
panel_99<-cbind(panel[,1:12], panel_99)

write.csv(panel_99, "data/panel_99.csv")
panel_99<- read.csv("panel_99.csv", row.names = 1)

#set universal plot settings
theme_feature<-theme_classic()+
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(hjust=0), legend.title = element_blank(),
        legend.position = "right", legend.text=element_text(size=14),
        title = element_text(size=40), axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))
YlOrRd<-colorRampPalette(c("Yellow","Orange", "Red"))

#mix up cells to oveid overplotting by subsequent groups

panel_99<-panel_99[order(sample(rownames(panel_99))),]


#Print UMAP plots

# Ccell clusters
panel_99$integrated_snn_res.1.1<-factor(panel_99$integrated_snn_res.1.1, levels=c(0:22)) 

pdf(file="UMAP_res1.1.pdf", width=12, height=9)
ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2, color=integrated_snn_res.1.1))+
  geom_jitter(size=0.7, alpha=0.5)+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  theme_feature
dev.off()

# Replicates
pdf(file="UMAP_replicate.pdf", width=11, height=9)
ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2, color=replicate))+
  geom_jitter(size=0.6, alpha=0.7)+
  scale_color_manual(values=c( "green1", "skyblue2","gold1"))+
  theme_feature+
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

# control vs tumor
panel_99$day<-factor(panel_99$day)

pdf(file="plots/UMAP_day.pdf", width=12, height=9)
ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2, color=day))+
  geom_jitter(size=0.5, alpha=0.8)+
  scale_color_manual(values=c("grey", "#786BC0", "#F44686"))+
  theme_feature+
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

#sex
pdf(file="pdf plots for pub/UMAP_sex.pdf", width=10, height=9)
ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2, color=sex))+
  geom_jitter(size=0.5, alpha=0.8)+
  scale_color_manual(values=c("gold",  "royalblue"))+
  theme_feature+
  guides(colour = guide_legend(override.aes = list(size=7)))

dev.off()


#print cell type UMAP

#adjust colors


panel_99$Cell_annotation3<-factor(panel_99$Cell_annotation3)
panel_99[panel_99$Cell_annotation3=="DCs", "Cell_annotation3"]<-"cDCs"
panel_99$Cell_annotation3<-factor(panel_99$Cell_annotation3, levels=c("Hom-MG_1", "Hom-MG_2", "Hom-MG_3",
                                                                      "Act-MG_1", "Act-MG_2","Act-MG_3","Act-MG_4",
                                                                       "ncMo",  "Mo", "Int", "Mphi_1", "Mphi_2","Mo-DCs", "cDCs", "BAM","NK cells", 
                                                                       "T-cells", "CD24a+", 
                                                                       "Ncam1+","UN"))

color_code<-read.csv("data/color_code3.csv", row.names = 1)
color_code<-color_code[order(match(color_code$Cell_annotation3, levels(panel_99$Cell_annotation3) )), ]

pdf(file="pdf_pub2//cell_type_3_v3_legend.pdf", width=9, height=11)
ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2, color=Cell_annotation3))+
  geom_jitter(size=0.6, alpha=0.5)+
  scale_color_manual(values= color_code$color)+
  theme_feature+
  #facet_grid(.~day)+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  theme(legend.position = "bottom")

dev.off()


#print top expressed genes per cluster
clusters<-panel_99$integrated_snn_res.1 %>% unique %>% as.character
#for (e in clusters){
e<-0
  filename<-paste0("plots/Top_clust_",e,".jpeg")
  print(filename)
  genes<-top50_1[top50_1$cluster==e, "Gene.name2"][c(1:20),] %>% as_vector()
  graphs<-list()
  for (i in seq_along(genes)){
    max<-round(max(panel_99[,genes[i]]), digits=1)
    min<-0 #ceiling(min(panel_99[,genes[i]]))
    plot<-ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2))+
      geom_jitter(size=1, alpha=0.5, aes_string(color=as.name(genes[i])))+
      labs(title=genes[i])+
      scale_color_gradientn(colors= (c(viridis(50),
                                     YlOrRd(50))))+
      #breaks=c(min,max))+
      #, labels=c("min", "max"))+
      xlab("UMAP_1")+
      ylab("UMAP_2")+
      theme_feature+
      theme(legend.text = element_text(size=20), legend.position = "right")
    graphs[[i]]<-plot
    print(paste(max,min))
  } 
  jpeg(file=filename, width=2500, height=2000)
  ggarrange(plotlist=graphs, ncol=5, nrow=4, common.legend = F)
  dev.off() 
  print("done")
#}



#cell heatmaps
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(scales)
library(colorspace)
library(stats)
library(dendsort)
library(data.table)

MMmarkers<-c("Tmem119", "P2ry12", "Selplg", "Fcrls", "Gpr34", "Olfml3",  "Sparc", "H2.Aa", "H2.Ab1",
              "Ptprc","Cd52", "Lgals3", "Itga4", "Tgfbi", "Ly6c2", "Ly6i", "Cd74")

annotation<-panel_99[panel_99$Cell_annotation1 %in% c("MG", "Mo/Mphi"),]
annotation$Cell_annotation3<-factor(annotation$Cell_annotation3, levels = c("Hom-MG_1", "Hom-MG_2", "Hom-MG_3",
                                                                            "Act-MG_1", "Act-MG_2",  "Act-MG_3",  "Act-MG_4",
                                                                            "ncMo", "Mo", "intMoMphi",
                                                                            "Mphi_1", "Mphi_2", "Mphi_3"))
#NA values for Ccr2 in rep1
annotation[annotation$replicate=="r1", "Ccr2-prot-raw"]<-NA

MM.RNA <- annotation[,MMmarkers]
MM.prot <- annotation[annotation$Cell_annotation1 %in% c("MG", "Mo/Mphi"), grep("prot", colnames(panel_99))]
MM.prot<-MM.prot[,-c(3,10)]

col_annotation=c(rep("RNA", length(MMmarkers)), rep("protein", 10))

MM<-cbind(MM.RNA, MM.prot)

quantile(a, c(.01, .99))
#col_fun<-colorRamp2(c(-1.5,0,3), brewer_pal(palette = "BuPu")(3))
col_fun<-colorRamp2(c(-1.5,0,3), inferno(3))
#col_fun<-colorRamp2(c(-3,-1,0,1,3), diverge_hcl(5, palette="Vik"))
a<-scale(MM)

#perform clustering separately to apply dendrogram sorting
clust<-hclust(dist(a), method="ward.D2")
dend<-dendsort(clust, type="min")
color_code<-data.frame("Cell_annotation3" = names(table(panel_99$Cell_annotation3)),
                       "color"= hue_pal()(length(names(table(panel_99$Cell_annotation3)))))
color_code<-color_code[color_code$Cell_annotation3 %in% annotation$Cell_annotation3,]
color_code<-color_code[order(match(color_code$Cell_annotation3, levels(annotation$Cell_annotation3) )), ]
color<-color_code$color
names(color)<-color_code$Cell_annotation3

pdf(file="plots/Cell_heatmap_MMdendsort.pdf", width=12, height=6)
Heatmap(t(a), 
        col=col_fun,
        cluster_columns = dend,
        show_column_names = F,
        row_split=col_annotation,
        top_annotation = HeatmapAnnotation(#"Sex"=annotation$sex,
                                           #"Time point"=annotation$day,
                                           "Cell Type" = annotation$Cell_annotation3,
                                           col=list("Sex"=c("male"="blue", "female"="yellow"),
                                                    "Time point"=c("D0"="#BAA6B1", "D14"="#53424C", "D21"="#EB23B3"),
                                                    "Cell Type"=color)))
dev.off()


#generate plots
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(colorspace)
library(ggpubr)
library(dplyr)
library(ggridges)
library(latex2exp)
library(scales)
library(ggpubr)


#print sex UMAP

pdf(file="plots30/UMAP_sex.pdf", width=12, height=9)
ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2, color=sex))+
  geom_jitter(size=0.8, alpha=0.5)+
  #facet_grid(sex~day)+
  scale_color_manual(values=c( "#ffff00",  "#4FCAFF"))+
  theme_feature+
  guides(colour = guide_legend(override.aes = list(size=7)))

dev.off()


#Generate feature plots for transcriptomic profile groups

#Interferon response
genes<-colnames(panel_99)[grep("Ifi|Isg|Irf", colnames(panel_99))]

#Hypoxic signature
genes<-c("Bnip3", "Adam8", "Fam162a", "Mif", "Hilpda", "Arg1")

#canonical microglia
genes<-c( "Tmem119","Gpr34", "Fcrls", "P2ry12", "Cx3cr1", "Selplg", "Olfml3", "Siglech")

#Trascription factors/co-activators
genes<-c("Klf2", "Egr1", "Cited2", "Klf4", "Fos", "Atf3", "Fosb", "Ier2", 
      "Jun", "Junb", "Jund" )

#Tumor-activated
genes-c("Cd74", "H2_Aa", "H2_Ab1", "H2_Eb1",
                   "Stat1", "Ccl12", "C4b", "Lgals3bp",
                   "Bst2", "B2m", "Cd52", "Ifitm3")

#Cytokines
genes<-c("Tnf", "Ccl2", "Ccl3", "Ccl4", "Il1a", "Il1b")

#Phagocytosis
genes<-c("Apoe", "Lgals3", "Fabp5", "Gpnmb", "Spp1")

proliferation<-c("Stmn1", "Top2a", "Birc5", "Tubb5", "Tubb4b", "Tuba1b", "Ube2c","Hmgb2", 
                 "Pclaf", "Ccnb1", "Ccnb2", "Ccna2","Cks1b",  "Cdk1")


graphs<-list()

for (i in seq_along(genes)){
  max<-round(max(panel_99[,genes[i]]), digits=1)
  min<-0 #ceiling(min(panel_99[,genes[i]]))
  plot<-ggplot(panel_99, aes(x=UMAP_1, y=UMAP_2))+
    geom_jitter(size=1, alpha=0.5, aes_string(color=as.name(genes[i])))+
    labs(title=genes[i])+
    scale_color_gradientn(colors= (c(viridis(50),
                                     YlOrRd(50))))+
                          #breaks=c(min,max))+
                          #, labels=c("min", "max"))+
    xlab("UMAP_1")+
    ylab("UMAP_2")+
    #xlab("")+
    #ylab("")+
    theme_feature+
    theme(legend.position = "right", legend.text = element_text(size=20))
  #facet_grid(sex~day)
  graphs[[i]]<-plot
  print(paste(max,min)) 
}
a<-2
b<-4

png(file="Panel_Feature.png", width=500*a, height=500*b)
ggarrange(plotlist=graphs, ncol=a, nrow=b, common.legend = F)
dev.off()


  
