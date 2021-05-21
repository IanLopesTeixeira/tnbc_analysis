setwd("C:/Users/TESTER/Desktop/Organizing/Projeto LBCM/Código")

# IMPORT WITH GENES ##########

library(tximport)


samples <- list.files(".", pattern = "kallisto")
files <- file.path(".", samples, "human_abundance.tsv")
names(files) <- substr(samples, 1, 8)
genes_h <- read.table("gene_human.tsv",header = T)
kallisto_h <- tximport(files, type = "kallisto", 
                       tx2gene = genes_h, txOut = FALSE)

# NORMALIZATION ############
library(edgeR)

cts <- kallisto_h$counts
group <- factor(c("Control","Control","Beva","Control","Beva"), 
                levels = c("Control", "Beva") )
y <- DGEList(counts=cts, group=group)

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
cpms <- edgeR::cpm(y, offset = y$offset, log = T)

#switch gene to symbol (cpms and y)
h_sep <- read.table("h_sep.tsv", header = T) #subtable of gene_total with human transcripts only
rownames(cpms) <- h_sep[match(rownames(cpms), h_sep[,2]),6]
cpms <- cpms[!duplicated(rownames(cpms)),]
rownames(y) <- h_sep[match(rownames(y), h_sep[,2]),6]
y <- y[!duplicated(rownames(y)),]


# SAMPLE LABELLING ################

pacman::p_load(pacman, tidyverse)
library("factoextra")
library(PCAtools)

# CLUSTERING #
values.cor <-as.dist(1 - cor(cpms))
dendgm_cpms <- hclust(values.cor)
plot(dendgm_cpms)

#Plot similarity between samples
fviz_dist(values.cor, gradient = list(low = "black",
                                      high = "white"))

#Plot dendogram with 2 groups
fviz_dend(dendgm_cpms, k = 2,
          cex = 1.2, # label size
          k_colors = c("darkblue","red"),
          color_labels_by_k = T,
          rect = T,
          main = "SampleType_Human")

# PCA #

#Prepare data
metadata <- data.frame(group)
rownames(metadata) <- colnames(cpms)
colnames(metadata) <- "SampleType"

#Run PCA
pca.res <- pca(cpms, metadata=metadata)

#Plot variance explained by each component
screeplot(pca.res,
          title = "SampleType_Human")

#Plot 2 selected components/eigenvectors
biplot(pca.res, colby="SampleType", hline = 0, vline = 0,legendPosition = 'top',
       title = "SampleType_Human") # Biplot with colors by sample type

#Plot several components
pairsplot(pca.res, colby="SampleType", pointSize = 2)

# Plot the component loadings and label genes most responsible for variation
plotloadings(pca.res, components = getComponents(pca.res, seq_len(2)), title = "SampleType_Human") #retaining 5% of the loadings per PC



# DIFFERENTIAL EXPRESSION ###################

design <- model.matrix(~group)
y <- estimateDisp(y, design)
plotBCV(y)

fit_q <- glmQLFit(y, design)
qlf <- glmQLFTest(fit_q,coef= 2)
plotQLDisp(qlf)
topTags(qlf)
is.de_q <- decideTestsDGE(qlf)
summary(is.de_q)
deq <- topTags(qlf, n = Inf, p = 0.05)$table
up_qlf <- row.names(deq[deq$logFC > 0,])
down_qlf <- row.names(deq[deq$logFC < 0,])


fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
plotMD.DGELRT(lrt)
topTags(lrt)
is.de_l <- decideTestsDGE(lrt)
summary(is.de_l)



#color data frame
dge<- data.frame(qlf)
dge$color <- "black"
dge$color[match(up_qlf,rownames(dge))] <- "red"
dge$color[match(down_qlf,rownames(dge))] <- "blue"


plot(x = dge[,1], y = -log10(dge[,4]), type = "p",col = dge$color, 
     xlab = "logFC", ylab = "-log10Pvalue", main = "DEG_Human",
     legend = legend("bottomright", legend=c("Up", "NotSig", "Down"),col=c("red","black", "blue"),pch = 20, cex=0.8))

#Heatmap

library(pheatmap)
sub_cpms <- data.matrix(cpms[match(rownames(deq)[1:20], rownames(cpms)),])
pheatmap(sub_cpms, main = "Top 20 DEG Human")



# PATHWAYS ##################

library(fgsea)

ranks <- qlf$table$logFC
names(ranks) <- rownames(qlf$table)
reactome <-  gmtPathways("c2.cp.reactome.v7.2.symbols.gmt")
fgseaRes <- fgsea(pathways = reactome,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)
topPathways <- fgseaRes[head(order(pval), n=10), pathway]

reactome_2all <-  gmtPathways("c2.all.v7.2.symbols.gmt")
fgseaRes_2all <- fgsea(pathways = reactome_2all,
                       stats    = ranks,
                       minSize  = 15,
                       maxSize  = 500)
topPathways_2all <- fgseaRes_2all[head(order(pval), n=10), pathway]

reactome_h <-  gmtPathways("h.all.v7.2.symbols.gmt")
fgseaRes_h <- fgsea(pathways = reactome_h,
                    stats    = ranks,
                    minSize  = 15,
                    maxSize  = 500)
topPathways_h <- fgseaRes_h[head(order(pval), n=10), pathway]

reactome_all <-  gmtPathways("msigdb.v7.2.symbols.gmt")
fgseaRes_all <- fgsea(pathways = reactome_all,
                      stats    = ranks,
                      minSize  = 15,
                      maxSize  = 500)
topPathways_all <- fgseaRes_all[head(order(pval), n=10), pathway]
comparative_df <- data.frame (topPathways_2all, topPathways_h,
                              topPathways_all)

plot.new()
plotGseaTable(reactome_2all[topPathways_2all], ranks, fgseaRes_2all, gseaParam=0.8)
collapsedPathways <- collapsePathways(fgseaRes_2all [order(padj)][padj < 0.05], reactome_2all, ranks)
mainPathways <- fgseaRes_2all[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
plot.new()
plotGseaTable(reactome_2all[mainPathways[1:10]], ranks, fgseaRes_2all, gseaParam = 0.8)
plot.new()
plotGseaTable(reactome_2all[tail(mainPathways, n=10)], ranks, fgseaRes_2all, gseaParam = 0.8)
down_pathways <- fgseaRes_2all[fgseaRes_2all$pathway %in% tail(mainPathways, n=7),]
down_pathways <- down_pathways[order(down_pathways$NES, decreasing = T),]
up_pathways <- fgseaRes_2all[fgseaRes_2all$pathway %in% head(mainPathways, n=10),]
up_pathways <- up_pathways[order(up_pathways$NES, decreasing = T),]
select_pathways <- rbind(up_pathways,down_pathways)
bp <- barplot(select_pathways$NES, cex.names= 1, cex.axis = 1, names.arg = "", xlim=range(pretty(c(min(select_pathways$NES), max(select_pathways$NES)))),
              horiz = T, col = c(rep("red",10),rep("blue",7)), xlab="Normalized Enrichment Score (NES)")
text(x=-0.1, y=bp[1:10],gsub("_"," ",sub("REACTOME_","",select_pathways$pathway[1:10])), cex=0.5, adj=1, font=2)
text(x=0.1, y=bp[11:17],gsub("_"," ",sub("REACTOME_","",select_pathways$pathway[11:17])), cex=0.5, adj=0, font=2)
