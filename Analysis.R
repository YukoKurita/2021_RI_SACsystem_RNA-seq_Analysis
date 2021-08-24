# This is the R script used for RNA-Seq analysis in the following paper.
# "Visualization of seasonal phosphorus re-translocation, and expression profile of phosphate transporters in a shortened annual cycle system of the deciduous poplar tree"
# Kurita et al. 

oldpar <- par(no.readonly = TRUE)  # Save the default values.

# read sample attribute, Gene descriptions, expected count matrix -------
tmp <- read.csv("input/SampleAttribute_lab.csv")
rownames(tmp) <- tmp$name
dim(tmp) # 44  8
at <- list(L = tmp[tmp$tissue == "Leaf",],
           S = tmp[tmp$tissue == "Stem",])


des <- read.csv("input/Pt3.1_GeneDescription.csv", row.names = 1)
dim(des) # 42950     8

tmp <- read.csv("input/excnt_lab.csv", row.names = 1)
dim(tmp) # 42950    44
excnts <- list(L = tmp[,rownames(at$L)],
               S = tmp[,rownames(at$S)])

n <- c("Leaf", "Stem")

# plot raw data counts ------------------------------------------------
dir.create("plot_sub")
pdf("plot_sub/transcript_couts_all.pdf", width = 8, height = 5)
par(mfcol = c(1,2), oma = c(0, 0, 1, 0))
for(i in 1:length(excnts)){
  tmp <- colSums(excnts[[i]])
  plot(tmp, ylim = c(0, 10^7), ylab = "counts", xlab = "sample",
       main = sprintf("%s, min: %s", n[i], min(tmp)))
  abline(a=10^6, b=0, lty="dashed") 
  cat(min(tmp),"\n")
}
mtext(text = "Transcript expected counts", side = 3, line = -0.5, 
      outer=T, cex=1.2, adj = 0)
dev.off()

# make DGElists -------------------------------------------------------
library(edgeR)
ds <- list(L = DGEList(counts = excnts$L, group = factor(at$L$stageID)),
           S = DGEList(counts = excnts$S, group = factor(at$S$stageID)))

# cut off non-expressed genes ------------------------------------------
ds2 <- NULL
for(i in 1:length(ds)){
  cpm <- cpm(ds[[i]])
  lcpm <- cpm(ds[[i]], log=TRUE)
  
  L <- mean(ds[[i]]$samples$lib.size)*(1e-6)
  M <- median(ds[[i]]$samples$lib.size)*(1e-6)
  cat(c(L, M),"\n")  # mean and median
  
  table(rowSums(ds[[i]]$counts==0)==ncol(ds[[i]]$counts))
  keep.exprs <- filterByExpr(ds[[i]], group=ds[[i]]$samples$group)
  d <- ds[[i]][keep.exprs,, keep.lib.sizes=FALSE]
  
  ## plot cutt off
  lcpm.cutoff <- log2(10/M + 2/L)
  library(RColorBrewer)
  nsamples <- ncol(d)
  col <- brewer.pal(nsamples, "Paired")
  
  png(filename = sprintf("plot_sub/cutoff_%s.png", n[i]), width = 400, height = 350)
    par(mfrow=c(1,2))    # before cut off
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(main="A. Raw data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (j in 2:nsamples){
      den <- density(lcpm[,j])
      lines(den$x, den$y, col=col[j], lwd=2)
    }
    legend("topright", colnames(d), text.col=col, bty="n")
    
    lcpm <- cpm(d, log=TRUE)   # after cut off
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(main="B. Filtered data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (j in 2:nsamples){
      den <- density(lcpm[,j])
      lines(den$x, den$y, col=col[j], lwd=2)
    }
    legend("topright", colnames(d), text.col=col, bty="n")
    mtext(text = sprintf("cut off: %s, exgenes: %s", names(ds)[i], nrow(d)),
          side = 3, line = -1.2, outer=T, cex=1.5, adj = 0)
  dev.off()
  
  ds2 <- c(ds2,list(d))
}

names(ds2) <- names(ds)
dir.create("data_Robj")
save(ds2, file = "data_Robj/ds2_cutoff_DGElist")
rm(ds, col)  # if need, run "# make DGElists"

# expressed gene list ---------------------------------------------------
egl <- list(L = rownames(ds2$L$counts),
            S = rownames(ds2$S$counts))

summary(egl) # Leaf 23646 genes, Stem 25151 genes

# Differential expression analysis --------------------------------------
efits <- NULL
vs <- NULL
for (i in 1:length(ds2)) {
  design <- model.matrix(~0+ds2[[i]]$samples$group)
  colnames(design) <- levels(ds2[[i]]$samples$group)
  design
  
  tmp <- colnames(design) # make contr.matrix label
  tmp2 <- NULL
  for (j in 1:(length(tmp)-1)) {
    tmp3 <- tmp[j]
    for(k in 1:(length(tmp)-j)){
      tmp2 <- c(tmp2, paste(tmp3,"-",tmp[j+k]))
    }
  }
  
  tmp2
  contr.matrix <- makeContrasts(contrasts=tmp2,levels = colnames(design))
  contr.matrix
  
  png(filename = sprintf("plot_sub/voom_%s_s4.png", n[i]), width = 400, height = 350)
  par(mfrow=c(1,2))
  v <- voom(ds2[[i]], design, plot=TRUE, save.plot=TRUE)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  plotSA(efit, main="Final model: Mean-variance trend")
  dev.off()
  
  vs <- c(vs, list(v))
  efits <- c(efits,list(efit))
}

par(oldpar)
names(vs) <- names(ds2)
save(vs, file = "data_Robj/vs")
names(efits) <- names(ds2)
save(efits, file = "data_Robj/efits")


# ANOVA, p < 0.01, log-fold changes > 1 -------------------------------------------
Flists<- NULL
for (i in 1:length(efits)) {
  tmp <- topTableF(efits[[i]], number = nrow(ds2[[i]]$counts), p.value = 0.01, lfc = 1)
  Flists <- c(Flists, list(tmp))
  cat(sprintf("%s: %s / %s genes \n",names(efits[i]),nrow(tmp), nrow(ds2[[i]]$counts)))  
  # L: 15743 / 23646 genes 
  # S: 17234 / 25151 genes 
}

names(Flists) <- names(ds2)

# clustering of exprettion trend ----------------------------------------------------
## take mean(log2(cpm))
mlcpms <- NULL
for (i in 1:length(ds2)) {
  lcpm <- cpm(ds2[[i]], log=TRUE) 
  tmp <- lcpm[rownames(Flists[[i]]),]
  cn <- levels(ds2[[i]]$samples$group)
  mlcpm <- matrix(NA, nrow = nrow(tmp), ncol = length(cn))
  rownames(mlcpm) <- rownames(tmp) 
  
  for(j in 1:nrow(tmp)){
    mlcpm[j,] <- tapply(tmp[j,], ds2[[i]]$samples$group, mean)
  }
  colnames(mlcpm) <- cn
  mlcpms <- c(mlcpms, list(mlcpm))
}

names(mlcpms) <- names(ds2)
str(mlcpms)  # L: num [1:15743, 1:5], S: num [1:17234, 1:8] 
rm(mlcpm)

# k-meand clustering -----------------------------------------------------------------
clus.No <- c(3,5) 
clusts <- NULL

for (i in 1:length(mlcpms)) {
  tmp <- scale(t(mlcpms[[i]]))
  tmp <- t(tmp)
  km <- kmeans(tmp, centers = clus.No[i], nstart = 1) # check warning
  clust <- km$cluster
  save(clust, file = sprintf("data_Robj/%s_clust%s", n[i], clus.No[i]))
  clusts <- c(clusts, list(clust))
}

names(clusts) <- names(ds2)
str(clusts)  # L: Named int [1:15743], S: Named int [1:17234]
rm(clust)

# NOTE: The numbering of the k-means clustering changes each time. -------------------
# We load and use the data we used in the paper to keep the names of the clusters the same as in the figure.
clusts <- NULL
load("input/Leaf_clust3") # name: clust
clusts <- c(clusts, list(clust))

load("input/Stem_clust5")
clusts <- c(clusts, list(clust))
names(clusts) <- names(ds2)
str(clusts) 
rm(clust)

# draw cluster plot Fig.6ab ----------------------------------------------------------
library(ggplot2)
library(tidyr)
dir.create("plot_Fig")
coltmp <- c("#00CC6620", "#FF773E20")

## set cluster plot function
fun.cluster_plot <- function(x,t,col){
  x <- as.data.frame(x)
  stage <- rownames(x)
  x <- cbind(stage, x)
  tmp2 <- tidyr::gather(x, key=gene, value = scaled_exp_value, -stage)
  
  g <- ggplot(tmp2, aes(x = stage, y = scaled_exp_value, group = gene)) + 
    geom_line(col = col) + xlab("") + ylab("scaled_value") +
    ggtitle(t) + 
    theme_bw(base_size=15) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15),
          axis.text.x = element_text(angle = 45, hjust = 1))
  plot(g)
}

## clustering plot for paper
lab <- list(L = c("Cluster C","Cluster A","Cluster B"), 
            S = c("Cluster C","Cluster E","Cluster D","Cluster A","Cluster B"))

for (i in 1:length(clusts)) {
  for(j in 1:clus.No[i]){
    cat(n[i], max(clusts[[i]]), j, "\n")
    tmp <- mlcpms[[i]][clusts[[i]] == j, ]
    tmp <- scale(t(tmp))
    title <- sprintf("%s, %s genes", lab[[i]][j], ncol(tmp))
    png(file=sprintf("plot_Fig/%s_clust%s_%s.png", n[i], clus.No[i], j),
        width = 320, height = 300)
    fun.cluster_plot(x=tmp, t=title, col = coltmp[i])
    dev.off()
  }
}

# make gene list, Supplementary Table 2 and 3 ----------------------------------------
exgenes <- list(L = rownames(des[des$NormalizationGroup == "data",]) %in% rownames(ds2$L[["counts"]]),
                S = rownames(des[des$NormalizationGroup == "data",]) %in% rownames(ds2$S[["counts"]]))

SIgenes <- list(L = rownames(des[des$NormalizationGroup == "data",]) %in% rownames(Flists$L),
                S = rownames(des[des$NormalizationGroup == "data",]) %in% rownames(Flists$S))

sum(exgenes$L) # 23646
sum(exgenes$S) # 25151
sum(SIgenes$L)  # 15743
sum(SIgenes$S)  # 17234

gls <- NULL
for(i in 1:length(exgenes)){　# It'll take a while.
  gl <- cbind(exgenes[[i]], SIgenes[[i]], matrix(NA, length(exgenes[[i]]),2))
  gl <- as.data.frame(gl)
  rownames(gl) <- rownames(des[des$NormalizationGroup == "data",])
  colnames(gl) <- c("Expressed gene", "SI gene", "adj.P.val", "Cluster")
  
  for (j in rownames(Flists[[i]])) {
    gl[j,3] <- Flists[[i]][j,]$adj.P.Val
  }
  
  for (j in names(clusts[[i]])){
    gl[j,4] <- clusts[[i]][j] 
  }
  
  gl$Cluster <- as.character(gl$Cluster)
  cat(mode(gl$Cluster), "\n")
  gls <- c(gls, list(gl))
}

names(gls) <- names(exgenes)
gls$L$Cluster <- chartr("123","CAB",gls$L$Cluster)
gls$S$Cluster <- chartr("12345","CEDAB",gls$S$Cluster)

dir.create("sup_tables")
write.csv(gls$L, file = "sup_tables/Leaf_genelist.csv")
write.csv(gls$S, file = "sup_tables/Stem_genelist.csv")

# GO analysis
source("input/GOanalysis_functions.R")
load("input/ulg.Potri_GObyTAIR_190903") # name: ulg
dim(ulg)  # 1897916       3

for (i in 1:length(ds2)) {
  # make subset of ulg
  ulg2 <- ulg[is.element(ulg[,"transcript_id"], rownames(ds2[[i]])),]
  cat(n[i], dim(ulg2), "\n")   
  colnames(ulg2) <- c("locus","GOid","Ath_locus")
  
  # fisher's exact test for each cluster
  for(j in 1:clus.No[i]){
    tmp <- mlcpms[[i]][clusts[[i]] == j, ]
    result <- ng.mft(cgt=ulg2, gn.test=rownames(tmp))
    result2 <- ng.prepGOtestOutTable(result, alpha = 0.05)
    result2 <- as.data.frame(result2, stringsAsFactors = F)
    result2 <- result2[!(result2$Description == "NA"),]
    
    result2$`Adjusted P value` <- as.numeric(result2$`Adjusted P value`)
    result2 <- result2[order(result2$`Adjusted P value`),]
    write.csv(result2, file = sprintf("sup_tables/%s_GO_cluster_%s.csv", n[i], j))
  }
}

# GO plot from REVIGO file, Fig. 6cd -------------------------------------------------------
# NOTE: The above GO analysis output files were summarized using REVIGO (http://revigo.irb.hr) (Supek et al., 2011) 
# allowed similarity = 0.5,  Others are default.
# The REVIGO output files(only biological process) need to be placed in a folder "REVIGOout"
# This script was written with reference to the previous protocol, Bonnot T, Gillard M, Nagel D. 2019. "A Simple Protocol for Informative Visualization of Enriched Gene Ontology Terms." Bio-Protocol e3429. 

## get GOlist  
ord <- list(L = c(2,3,1), S = c(4,5,1,3,2))
GOlists <- NULL

for (i in 1:length(ord)) {
  GOlist <- NULL
  for(j in ord[[i]]){　　　# check cluster order
    fn <- sprintf("sup_tables/%s_GO_cluster_%s.csv", n[i], j) # check cluster No
    tmp <- read.csv(fn, row.names = 1,stringsAsFactors = F) 
    fn <- sprintf("REVIGOout/%s_REVIGO_cluster%s.csv", n[i], j)
    tmp2 <- read.csv(fn, stringsAsFactors = F)
    
    tmp2 <- tmp2[tmp2$eliminated == 0,]
    tmp <- tmp[tmp2$term_ID,]
    tmp$clus <- rep(sprintf("Cluster%s", j), nrow(tmp))
    rownames(tmp) <- NULL
    tmp <- tmp[order(tmp$Adjusted.P.value),]
    
    GOlist <- rbind(GOlist, tmp)
  }
  
  colnames(GOlist) <- c("Adjusted.P.value", "ID", "Description", "Gene_number",
                        "A", "B", "U", "clus" )
  GOlist$mlog10P <- -log10(GOlist$Adjusted.P.value)
  GOlist$ylabel <- paste(GOlist$Description, ":", GOlist$ID)
  GOlist$ylabel <- factor(GOlist$ylabel, levels = unique(GOlist$ylabel))
  
  GOlists <- c(GOlists, list(GOlist))
}

names(GOlists) <- names(ord)

GOlists$L$clus <- factor(GOlists$L$clus, levels = c("Cluster2", "Cluster3", "Cluster1"))
GOlists$S$clus <- factor(GOlists$S$clus, levels = c("Cluster4", "Cluster5", "Cluster1", 
                                              "Cluster3", "Cluster2"))

## draw GO plot 
library(ggplot2)
xlabels <- list(L = c("Cluster A","Cluster B","Cluster C"),
                S = c("Cluster A","Cluster B","Cluster C","Cluster D", "Cluster E"))

ws <- c(5.5, 7.6)
hs <- c(8.3, 13)

for (i in 1:length(GOlists)) {
  tmp <- GOlists[[i]]
  
  fn <- sprintf("plot_Fig/%s_GOplot.pdf", n[i])
  pdf(fn, width = ws[i], height = hs[i]) 
  g <- ggplot(tmp, aes(x = ylabel, y = clus)) +
      geom_point(data = tmp,
                 aes(x = ylabel, y = clus, size = Gene_number, colour = mlog10P), 
                 alpha = .7)+
      scale_y_discrete(labels =xlabels[[i]])+
      scale_color_gradient(low = "cyan", high = "blue", limits=c(0, NA))+
      coord_flip()+
      theme_bw()+
      theme(axis.ticks.length=unit(-0.1, "cm"),
            axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1),
            axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
            axis.text = element_text(color = "black"),
            axis.title.x=element_blank(),
            legend.key.size = unit(0.4, "cm"),
            legend.text=element_text(size=8),
            legend.title=element_text(size=8))+
      xlab("GO biological processes")+
      labs(color=expression(atop("-log"[10],paste("(Adj. p)"))), size="Number\nof genes")
  plot(g)
  dev.off()
}

# expression heatmap, Fig. 7 ----------------------------------------------------------
library(gplots)

fun.heatmap <- function(mlcpm, tage, lab, tgl){
  cols <- colorRampPalette(c("navyblue","blue3","purple","orange","yellow"))
  mycol <- cols(1000)
  
  heatmap.2(mlcpm[tage,], scale="row",labCol=lab, 
            labRow = sprintf("%s, %s", tage, tgl[tage,]$name),
            col=mycol, trace="none", density.info="none", 
            margin=c(6,22), dendrogram="none", Colv=F, Rowv=F,
            cexCol = 1.2, keysize = 1)
}

# heatmap for Fig
tgl <- read.csv("input/TG_RI_paper.csv", stringsAsFactors = F)
rownames(tgl) <- tgl[,1]

## Leaf
tgl2 <- tgl[tgl$Leaf == T,]
tgl2 <- tgl2[rownames(tgl2) %in% rownames(Flists$L),]
tmp <- scale(t(mlcpms$L[rownames(tgl2),]))
tmp <- t(tmp)
tgl2 <- tgl2[order(tmp[,5],decreasing = F),]
tgl2 <-  tgl2[c(which(tgl2$group == "PHT1"),which(tgl2$group == "PHT5"),which(tgl2$group == "PHO"),
                which(tgl2$group == "Lu"), which(tgl2$group == "DPD1")),]
tage <- rownames(tgl2)

lab <- c("st1w4","st2w2","st2w4","st3w3","st3w5")

pdf("plot_Fig/Leaf_heatmap.pdf", width=5.7, height=5.8)
fun.heatmap(mlcpm = mlcpms$L, tage = tage, lab = lab, tgl = tgl)
dev.off()

## Stem
tgl2 <- tgl[tgl$Stem == T,]
tgl2 <- tgl2[rownames(tgl2) %in% rownames(Flists$S),]
tmp <- scale(t(mlcpms$S[rownames(tgl2),]))
tmp <- t(tmp)
tgl2 <- tgl2[order(tmp[,1],decreasing = T),]
tgl2 <- tgl2[c(which(tgl2$group == "PHT1"),which(tgl2$group == "PHT5"),which(tgl2$group == "PHO")),]
tage <- rownames(tgl2)

lab <- c("st1w4","st2w2","st2w4","st3w3","st3w5","st3w8","st1'w1","st1'w2")

pdf("plot_Fig/Stem_heatmap.pdf", width=6.2, height=5.7)
fun.heatmap(mlcpm = mlcpms$S, tage = tage, lab = lab, tgl = tgl)
dev.off()


# rawdata plot, Supplementary Figure S6-9 -----------------------------------------
library(ggplot2)
library(multcomp)

## Leaf
tgl2 <- tgl[tgl$Leaf == T,]
tage <- rownames(tgl2)
lcpm <- cpm(ds2$L, log=TRUE) 
tage <- tage[tage %in% rownames(lcpm)]
cl <- c("cluster C","cluster A","cluster B")

pdf("plot_Fig/Leaf_rawPlot.pdf", width=3.5, height=3.8)
for (i in 1:length(tage)){
  pn <- tage[i]
  log2cpm <- lcpm[pn,]
  stage <- as.character(ds2$L$samples$group)
  levels(stage) <- c("st1w4", "st2w2", "st2w4", "st3w3", "st3w5")
  x <- data.frame(stage, log2cpm)
  
  if(is.na(clusts$L[pn]) == F){
    res1 <- aov(log2cpm~stage,d=x)
    res2 <- glht(res1,linfct=mcp(stage="Tukey"))
    res3 <- cld(res2, level = 0.05)
    s <- res3$mcletters$Letters
  }else{
    s <- rep("", length(levels(stage)))
  }
  
  m <- tapply(x$log2cpm, x$stage, max)
  m <- m + 2
  m[m > 14] <- 14
  
  t <- sprintf("%s, %s", pn, tgl[pn,]$name)
  
  g <- ggplot(x, aes(x = stage, y = log2cpm)) + ylim(-2,14) +
    geom_point(shape = 1, size = 3) + xlab("") + ylab(expression(paste(log[2],"(cpm)"))) +
    ggtitle(t) + 
    theme_bw(base_size=10) + 
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=13),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.margin= unit(c(0.5, 2, 0, 0.5), "lines")) +
    annotate("text", x=4.8, y=13.5, label=cl[clusts$L[pn]]) + 
    annotate("text", x=names(m), y=m, label=s, col = "green3")
  plot(g)
}
dev.off()

## Stem
tgl2 <- tgl[tgl$Stem == T,]
tage <- rownames(tgl2)
lcpm <- cpm(ds2$S, log=TRUE) 
tage <- tage[tage %in% rownames(lcpm)]
cl <- c("cluster C","cluster E","cluster D","cluster A","cluster B")

pdf("plot_Fig/Stem_rawPlot.pdf", width=3.5, height=3.8)
for (i in 1:length(tage)){
  pn <- tage[i]
  log2cpm <- lcpm[pn,]
  stage <- ds2$S$samples$group
  levels(stage) <- c("st1w4", "st2w2", "st2w4", "st3w3", "st3w5", "st3w8", "st1'w1", "st1'w2")
  x <- data.frame(stage, log2cpm)
  
  if(is.na(clusts$S[pn]) == F){
    res1 <- aov(log2cpm~stage,d=x)
    res2 <- glht(res1,linfct=mcp(stage="Tukey"))
    res3 <- cld(res2, level = 0.05)
    s <- res3$mcletters$Letters
  }else{
    s <- rep("", length(levels(stage)))
  }
  
  m <- tapply(x$log2cpm, x$stage, max)
  m <- m + 2
  m[m > 14] <- 14
  
  t <- sprintf("%s, %s", pn, tgl[pn,]$name)
  
  g <- ggplot(x, aes(x = stage, y = log2cpm)) + ylim(-2,12) +
    geom_point(shape = 1, size = 3) + xlab("") + ylab(expression(paste(log[2],"(cpm)"))) +
    ggtitle(t) + 
    theme_bw(base_size=10) + 
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=13),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.margin= unit(c(0.5, 2, 0, 0.5), "lines")) +
    annotate("text", x=7.2, y=11.5, label=cl[clusts$S[pn]]) +
    annotate("text", x=names(m), y=m, label=s, col = "tan1")
  plot(g)
}
dev.off()


###### END ##########################################################################
# Package versions
sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] multcomp_1.4-12      TH.data_1.0-10       MASS_7.3-51.5        survival_3.1-11     
# [5] mvtnorm_1.1-0        gplots_3.0.3         GO.db_3.10.0         AnnotationDbi_1.48.0
# [9] IRanges_2.20.2       S4Vectors_0.24.4     Biobase_2.46.0       BiocGenerics_0.32.0 
# [13] tidyr_1.0.2          ggplot2_3.3.0        RColorBrewer_1.1-2   edgeR_3.28.1        
# [17] limma_3.42.2        
# 
# loaded via a namespace (and not attached):
#   [1] gtools_3.8.2       zoo_1.8-7          tidyselect_1.0.0   locfit_1.5-9.1    
# [5] purrr_0.3.3        splines_3.6.1      lattice_0.20-40    colorspace_1.4-1  
# [9] vctrs_0.2.4        blob_1.2.1         rlang_0.4.5        pillar_1.4.3      
# [13] glue_1.3.1         withr_2.1.2        DBI_1.1.0          bit64_0.9-7       
# [17] lifecycle_0.2.0    munsell_0.5.0      gtable_0.3.0       caTools_1.18.0    
# [21] codetools_0.2-16   memoise_1.1.0      labeling_0.3       Rcpp_1.0.3        
# [25] KernSmooth_2.23-16 scales_1.1.0       gdata_2.18.0       farver_2.0.3      
# [29] bit_1.1-15.2       digest_0.6.25      dplyr_0.8.5        grid_3.6.1        
# [33] tools_3.6.1        bitops_1.0-6       sandwich_2.5-1     magrittr_1.5      
# [37] tibble_2.1.3       RSQLite_2.2.0      crayon_1.3.4       pkgconfig_2.0.3   
# [41] ellipsis_0.3.0     Matrix_1.2-18      assertthat_0.2.1   rstudioapi_0.11   
# [45] R6_2.4.1           compiler_3.6.1 
