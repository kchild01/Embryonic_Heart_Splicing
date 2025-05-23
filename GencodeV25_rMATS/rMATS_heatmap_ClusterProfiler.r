library(tidyverse)
library("clusterProfiler")
library(enrichplot)
library("org.Hs.eg.db")
library("AnnotationHub")
library(ReactomePA)
library(DOSE)

#set path
directory <- "~/cotney_lab/heart/rMATS_heatmap"

setwd(directory)

##SE
cluster1 <- read.table("SE/GO/cluster1_genes.txt")
cluster2 <- read.table("SE/GO/cluster2_genes.txt")
cluster3 <- read.table("SE/GO/cluster3_genes.txt")

cluster1 <- gsub('\\..+$', '', cluster1$V1)
cluster2 <- gsub('\\..+$', '', cluster2$V1)
cluster3 <- gsub('\\..+$', '', cluster3$V1)

#DefaultLayer(cds1[["RNA"]]) <- 'data'
require(biomaRt)
load("~/Rscript/ensembl/ensembl_july2019.archive.RData")
mart <- ensembl

cluster1 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster1,
  uniqueRows=TRUE)

cluster2 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster2,
  uniqueRows=TRUE)

cluster3 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster3,
  uniqueRows=TRUE)

#gene ontology enrichment across clusters

#need to use Entrez id in clusterprofiler
genelist <- list("cluster1" = cluster1$entrezgene_id, 
                 "cluster2" = cluster2$entrezgene_id,
                 "cluster3" = cluster3$entrezgene_id
)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 10, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 10, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 10, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "SE/GO/ClusterProfiler/RMATSHM_cluster_enrichment_dotplot.pdf", w = 10, h = 12)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
dev.off()

p1 <- emapplot(BPclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "SE/GO/ClusterProfiler/rMATSHM_cluster_enrichment_emapplot.pdf", h = 8.5, w = 11)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
dev.off()

#per cluster enrichments
library('openxlsx')
go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "SE/GO/ClusterProfiler/gene_onology_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "SE/GO/ClusterProfiler/disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "SE/GO/ClusterProfiler/kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "SE/GO/ClusterProfiler/reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  
  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "SE/GO/ClusterProfiler/barplots_per_main_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist), function(x){
  
  name <- names(genelist)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

##A3SS
cluster1 <- read.table("A3SS/GO/cluster1_genes.txt")
cluster2 <- read.table("A3SS/GO/cluster2_genes.txt")
cluster3 <- read.table("A3SS/GO/cluster3_genes.txt")

cluster1 <- gsub('\\..+$', '', cluster1$V1)
cluster2 <- gsub('\\..+$', '', cluster2$V1)
cluster3 <- gsub('\\..+$', '', cluster3$V1)

#DefaultLayer(cds1[["RNA"]]) <- 'data'
require(biomaRt)
load("~/Rscript/ensembl/ensembl_july2019.archive.RData")
mart <- ensembl
#when ensembl is down
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE)

cluster1 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster1,
  uniqueRows=TRUE)

cluster2 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster2,
  uniqueRows=TRUE)

cluster3 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster3,
  uniqueRows=TRUE)

#need to use Entrez id in clusterprofiler
genelist <- list("cluster1" = cluster1$entrezgene_id, 
                 "cluster2" = cluster2$entrezgene_id,
                 "cluster3" = cluster3$entrezgene_id
)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 10, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 10, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 10, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "A3SS/GO/ClusterProfiler/A3SS_RMATSHM_cluster_enrichment_dotplot.pdf", w = 10, h = 12)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
dev.off()

p1 <- emapplot(BPclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "A3SS/GO/ClusterProfiler/A3SS_rMATSHM_cluster_enrichment_emapplot.pdf", h = 8.5, w = 11)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
dev.off()

#per cluster enrichments
library('openxlsx')
go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "A3SS/GO/ClusterProfiler/gene_onology_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "A3SS/GO/ClusterProfiler/disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "A3SS/GO/ClusterProfiler/kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "A3SS/GO/ClusterProfiler/reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  
  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "A3SS/GO/ClusterProfiler/barplots_per_main_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist), function(x){
  
  name <- names(genelist)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

##A5SS
cluster1 <- read.table("A5SS/GO/cluster1_genes.txt")
cluster2 <- read.table("A5SS/GO/cluster2_genes.txt")
cluster3 <- read.table("A5SS/GO/cluster3_genes.txt")

cluster1 <- gsub('\\..+$', '', cluster1$V1)
cluster2 <- gsub('\\..+$', '', cluster2$V1)
cluster3 <- gsub('\\..+$', '', cluster3$V1)

cluster1 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster1,
  uniqueRows=TRUE)

cluster2 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster2,
  uniqueRows=TRUE)

cluster3 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster3,
  uniqueRows=TRUE)

#need to use Entrez id in clusterprofiler
genelist <- list("cluster1" = cluster1$entrezgene_id, 
                 "cluster2" = cluster2$entrezgene_id,
                 "cluster3" = cluster3$entrezgene_id
)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 10, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 10, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 10, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "A5SS/GO/ClusterProfiler/A5SS_RMATSHM_cluster_enrichment_dotplot.pdf", w = 10, h = 12)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
dev.off()

p1 <- emapplot(BPclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "A5SS/GO/ClusterProfiler/A5SS_rMATSHM_cluster_enrichment_emapplot.pdf", h = 8.5, w = 11)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
dev.off()

#per cluster enrichments
library('openxlsx')
go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "A5SS/GO/ClusterProfiler/gene_onology_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "A5SS/GO/ClusterProfiler/disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "A5SS/GO/ClusterProfiler/kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "A5SS/GO/ClusterProfiler/reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  
  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "A5SS/GO/ClusterProfiler/barplots_per_main_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist), function(x){
  
  name <- names(genelist)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

##MXE
cluster1 <- read.table("MXE/GO/cluster1_genes.txt")
cluster2 <- read.table("MXE/GO/cluster2_genes.txt")
cluster3 <- read.table("MXE/GO/cluster3_genes.txt")

cluster1 <- gsub('\\..+$', '', cluster1$V1)
cluster2 <- gsub('\\..+$', '', cluster2$V1)
cluster3 <- gsub('\\..+$', '', cluster3$V1)

cluster1 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster1,
  uniqueRows=TRUE)

cluster2 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster2,
  uniqueRows=TRUE)

cluster3 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster3,
  uniqueRows=TRUE)

#need to use Entrez id in clusterprofiler
genelist <- list("cluster1" = cluster1$entrezgene_id, 
                 "cluster2" = cluster2$entrezgene_id,
                 "cluster3" = cluster3$entrezgene_id
)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 10, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 10, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 10, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "MXE/GO/ClusterProfiler/MXE_RMATSHM_cluster_enrichment_dotplot.pdf", w = 10, h = 12)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
dev.off()

p1 <- emapplot(BPclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "MXE/GO/ClusterProfiler/MXE_rMATSHM_cluster_enrichment_emapplot.pdf", h = 8.5, w = 11)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
dev.off()

#per cluster enrichments
library('openxlsx')
go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "MXE/GO/ClusterProfiler/gene_onology_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "MXE/GO/ClusterProfiler/disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "MXE/GO/ClusterProfiler/kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "MXE/GO/ClusterProfiler/reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  
  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "MXE/GO/ClusterProfiler/barplots_per_main_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist), function(x){
  
  name <- names(genelist)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

##RI
cluster1 <- read.table("RI/GO/cluster1_genes.txt")
cluster2 <- read.table("RI/GO/cluster2_genes.txt")
cluster3 <- read.table("RI/GO/cluster3_genes.txt")

cluster1 <- gsub('\\..+$', '', cluster1$V1)
cluster2 <- gsub('\\..+$', '', cluster2$V1)
cluster3 <- gsub('\\..+$', '', cluster3$V1)

cluster1 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster1,
  uniqueRows=TRUE)

cluster2 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster2,
  uniqueRows=TRUE)

cluster3 <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = cluster3,
  uniqueRows=TRUE)

#need to use Entrez id in clusterprofiler
genelist <- list("cluster1" = cluster1$entrezgene_id, 
                 "cluster2" = cluster2$entrezgene_id,
                 "cluster3" = cluster3$entrezgene_id
)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 10, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 10, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 10, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "RI/GO/ClusterProfiler/RI_RMATSHM_cluster_enrichment_dotplot.pdf", w = 10, h = 12)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:4])
dev.off()

p1 <- emapplot(BPclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))

cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
pdf(file = "RI/GO/ClusterProfiler/RI_rMATSHM_cluster_enrichment_emapplot.pdf", h = 8.5, w = 11)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:4])
dev.off()

#per cluster enrichments
library('openxlsx')
go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "RI/GO/ClusterProfiler/gene_onology_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "RI/GO/ClusterProfiler/disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "RI/GO/ClusterProfiler/kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "RI/GO/ClusterProfiler/reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  
  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "RI/GO/ClusterProfiler/barplots_per_main_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist), function(x){
  
  name <- names(genelist)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()
