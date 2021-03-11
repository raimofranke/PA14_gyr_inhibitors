# Treatment of PA14 with gyrase inhibitors and cystobactamids, growth until OD 0.6, exp. phase
# Script by Raimo Franke

# load libraries
library(xcms)
library(readxl)
library(CAMERA)
library(RColorBrewer)
library(pheatmap)

my.dir <- "/Users/rmf/Dropbox/_Paper/Metabolomics_data/"

# Get mzxml files of MS1 data
mzxmls <- list.files(paste0(my.dir,"/mzXML"), full.names = T, recursive = T)
sample_name <- sub(basename(mzxmls), pattern = "_1ul.mzXML", replacement = "", fixed = TRUE)
sample_name <- sub(sample_name, pattern = "p_HM-188_", replacement = "", fixed = TRUE)

## Create a phenodata data.frame
pd <- data.frame(sample_name = sample_name,
                 sample_group = c(rep("AR351", 3), rep("Ciprofloxacin", 3),
                                  rep("Coumermycin", 3), rep("CN-DM-861", 3),
                                  rep("Levofloxacin", 3), rep("Lomefloxacin", 3), rep("Novobiocin", 3),
                                  rep("untreated", 3)), stringsAsFactors = FALSE) 

raw_data <- readMSData(files = mzxmls, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk") 

# Define colours for the 8 groups
group_colors <- brewer.pal(8, "Set3")[1:8]
names(group_colors) <- levels(as.factor(pd$sample_group))


## Get the total ion current by file
tc <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "intensity", main = "Total ion current", las = 2, cex.axis = 0.7,
        names = pd$sample_group) 


#log intensities
tc <- split(log2(tic(raw_data)), f = fromFile(raw_data))

postscript("boxplot_log2_TIC.eps", width = 8, height = 6,
           horizontal = T, onefile = FALSE, paper = "special")

# standard margin: c(5.1, 4.1, 4.1, 2.1)
par(mar = c(6.5, 4.1, 4.1, 2.1))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "log2 intensity", main = "Total ion current", las = 2, names = pd$sample_group) 
dev.off()

## Define the settings for the centwave peak detection
cwp <- CentWaveParam(peakwidth=c(5, 25), ppm=10, snthresh = 100, mzdiff = 0.01,
                     prefilter = c(2, 1000), noise = 100)
xdata <- findChromPeaks(raw_data, param = cwp)

## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[raw_data$sample_group],
        las = 2, names = pd$sample_group, cex.axis = 0.7,
        ylab = expression(log[2]~intensity), main = "Peak intensities")


## Obiwarp.
xdata_obi <- adjustRtime(xdata, param = ObiwarpParam())

#make xcmsSet object
xset3 <- as(xdata_obi, "xcmsSet")
class.vector <- c("AR351", "AR351","AR351",
                  "Ciprofloxacin", "Ciprofloxacin", "Ciprofloxacin",
                  "Coumermycin", "Coumermycin", "Coumermycin",
                  "CN-DM-861", "CN-DM-861", "CN-DM-861",
                  "Levofloxacin", "Levofloxacin", "Levofloxacin",
                  "Lomefloxacin", "Lomefloxacin", "Lomefloxacin",
                  "Novobiocin", "Novobiocin", "Novobiocin",
                  "untreated", "untreated", "untreated")

sampclass(xset3) <- class.vector
xset4 <- group(xset3, method = "density", bw = 5, mzwid = 0.015, minfrac = 0.5, minsamp = 1)
xset5 <- fillPeaks(xset4)

# Annotation with CAMERA
# Create an xsAnnotate object
an <- xsAnnotate(xset5)
# Group based on RT
anF <- groupFWHM(an, perfwhm = 0.6)
# Annotate isotopes
anI <- findIsotopes(anF, mzabs = 0.01)
# Verify grouping
anIC <- groupCorr(anI, cor_eic_th = 0.75)
#Annotate adducts
anFA <- findAdducts(anIC, polarity="positive")

my_peaklist <- getPeaklist(anFA, intval = "into")


# Prepare peaktable
# my_peaktable <- peakTable(xset5)
my_groupval <- groupval(xset5, value = "into")
my_peaktable <- cbind(compound = rownames(my_groupval), my_peaklist)

# Filter RT
my_peaktable.filtered <- my_peaktable[my_peaktable$rt > 49 & my_peaktable$rt < 1201, ]
temp_colnames <- colnames(my_peaktable.filtered)
new_colnames <- c(temp_colnames[1:16], pd$sample_name, temp_colnames[41:43])
colnames(my_peaktable.filtered) <- new_colnames

# plot internal standards to check alignment and intensities
compound_names <- c("Caffeine", "Nortriptyline", "Naproxen", "Phenylalanine", "Tryptophan", "Nicotinamide")
m_vect <- c(195.0879, 264.1749, 231.1018, 166.0866, 205.0976, 123.0555)
t_vect <- c(435, 672, 767, 248, 368, 105)
m_tol <- 1e-5 #1 mDa
t_tol <- 50

pdf(file="EICs_IS_t1000.pdf", paper="a4r")
for (i in 1:length(m_vect)) {
  
  mmat <- matrix(nrow=1, ncol=2)
  mmat[1,] <- c(m_vect[i]-m_tol, m_vect[i]+m_tol)
  tmat <- matrix(nrow=1, ncol=2)
  tmat[1,] <- c(t_vect[i]-t_tol, t_vect[i]+t_tol)
  
  myEIC.raw <- getEIC(xset5, mzrange=mmat, rtrange=tmat,  rt = "raw")
  myEIC.corrected <- getEIC(xset5, mzrange=mmat, rtrange=tmat,  rt = "corrected")
  par(mfrow=c(2,1))
  plot(myEIC.raw, mzdec = 4)
  mtext(compound_names[i], side=1, line=3, outer=F, adj=0, cex=1)
  plot(myEIC.corrected, mzdec = 4)
  mtext(compound_names[i], side=1, line=3, outer=F, adj=0, cex=1)
  
}

dev.off()

# Result: Alignment worked fine, Caffeine is not well aligned, Nortriptyline
# and Naproxen are well aligned, Nortriptyline shows the least variance

# Normalization with the internal standard Nortriptyline
nortriptyline <- my_peaktable.filtered[which(my_peaktable.filtered$compound == "264.2/684"),]
norm_nortrip <- as.numeric(nortriptyline[,17:40]) / max(as.numeric(nortriptyline[,17:40]))
peaktable_normalized_IS <- my_peaktable.filtered
for(i in 1:nrow(peaktable_normalized_IS)){
  peaktable_normalized_IS[i, 17:40] <- peaktable_normalized_IS[i, 17:40]/norm_nortrip
}


# Removal of internal standard features, that would lead to problems after normalization with OD
IS_features <- read_excel("IS_features.xlsx")

peaktable_normalized_noIS <- peaktable_normalized_IS
for (i in 1:length(IS_features$compound)){
  peaktable_normalized_noIS <- subset(peaktable_normalized_noIS, compound != IS_features$compound[i])
}

# Normalization with OD
samples_OD <- read_excel("HM188_OD.xlsx")
norm_OD <- samples_OD$OD / max(samples_OD$OD)

peaktable_normalized_OD <- peaktable_normalized_noIS
for(i in 1:nrow(peaktable_normalized_OD)){
  peaktable_normalized_OD[i, 17:40] <- peaktable_normalized_OD[i, 17:40]/norm_OD
}

# impute missing values using imputeRosMinRand
# replace 0 with NA 
peaktable_normalized_OD_nozero <- peaktable_normalized_OD
tempmat <- peaktable_normalized_OD_nozero[, 17:40]
tempmat[tempmat==0] <- NA
tempmat2 <- imputeRowMinRand(tempmat)
peaktable_normalized_OD_nozero[, 17:40] <- tempmat2

# peaktable annotation
# annotation is done externally with Bruker DataAnalysis
peaktable_iso_anno <- read_excel("peaktable_iso_anno.xlsx")

# removal of features coming from antibiotics treatment

novo_features <- read_excel("novo_features.xlsx")

for (i in 1:length(novo_features$compound)){
  peaktable_iso_anno <- subset(peaktable_iso_anno, compound != novo_features$compound[i])
}


# substitute NAs in annotation with mz-IDs
for (i in 1:nrow(peaktable_iso_anno)){
  ifelse (peaktable_iso_anno$anno[i]=="NA",
          peaktable_iso_anno$anno[i] <- as.character(peaktable_iso_anno$compound[i]),
          print(peaktable_iso_anno$anno[i])
  )
}


# de-isotoping

for (i in 1:length(peaktable_iso_anno$compound)){
  peaktable_deiso_anno <- subset(peaktable_iso_anno, iso_anno != 1)
}

write.csv(peaktable_deiso_anno, "peaktable_deiso_anno.csv")

# hierarchical clustering
# scaling by devision with max(x)
peaktable_scaled <- peaktable_deiso_anno
for (i in 1:nrow(peaktable_scaled)) {
  peaktable_scaled[i,18:41] <- peaktable_scaled[i,18:41] / max(peaktable_scaled[i,18:41])
}
features_OD_normalized <- dist(t(peaktable_scaled[,18:41]))
plot(hclust(features_OD_normalized, method = "ward.D"), hang = -1, main = "", sub = "")

# generate dendrogram figure
postscript("dendro_features_ODnorm.eps", width = 8, height = 6,
           horizontal = T, onefile = FALSE, paper = "special")

# standard margin: c(5.1, 4.1, 4.1, 2.1)
par(mar = c(6.5, 4.1, 4.1, 2.1))
plot(hclust(features_OD_normalized, method = "ward.D"), hang = -1, main = "", sub = "")
dev.off()


# generate heatmap
# max scaled
my_mat <- peaktable_scaled[,18:41]
my_mat <- as.data.frame(my_mat)
#colnames(my_mat) <- pd$sample_name
rownames(my_mat) <- peaktable_scaled$anno
pheatmap(my_mat, clustering_method = "ward.D", show_rownames = T, fontsize_row = 1, scale = "none")

# max scaled and centred
my_mat <- scale(t(peaktable_scaled[,18:41]), center = T, scale = F)
my_mat <- as.data.frame(t(my_mat))


my_color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100)
pheatmap(my_mat, clustering_method = "ward.D", show_rownames = T, fontsize_row = 2, fontsize_col = 8,
         scale = "none", color = my_color)

# z scaled and centred
my_mat <- peaktable_normalized_OD_nozero[,18:41]
my_mat <- as.data.frame(my_mat)

my_color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100)
pheatmap(my_mat, clustering_method = "ward.D", show_rownames = T, fontsize_row = 2, fontsize_col = 8,
         scale = "row", color = my_color)


# PCA
data <- t(peaktable_scaled[,18:41])
colnames(data) <- peaktable_scaled$anno

pca_pt_scaled <- prcomp(data, scale. = F, center = T)
pca_results <- summary(pca_pt_scaled)

scores <- data.frame(pca_pt_scaled$x[,c("PC1","PC2")])
scores$class <- pd$sample_group

lab <- rownames(scores)

## Plot scores
scores.plot <- ggplot(data=scores,aes(x=PC1, y=PC2, colour=pd$sample_group)) +
  labs(colour = "groups") +
  geom_point(alpha = I(0.7), size=3)+
  geom_text(data=scores, mapping=aes(x=PC1, y=PC2, label=lab),size=3, vjust=1, hjust=-0.15)+
  #geom_hline(yintercept = 0)+
  #geom_vline(xintercept = 0)+
  theme_minimal(base_size = 12)+
  #theme_cowplot(12)
  xlab(paste0("PC1 (", round(pca_results$importance[2,1] *100, 2), "%)")) +
  ylab(paste0("PC2 (", round(pca_results$importance[2,2] *100, 2), "%)")) +
  theme(legend.position = "bottom")

scores.plot

# generate PCA figure

postscript("PCA_scores.eps", width = 10, height = 6,
           horizontal = F, onefile = FALSE, paper = "special")
scores.plot

dev.off()


# Boxplots of all features
class.factor <- factor(pd$sample_group, levels = c("Lomefloxacin", "Levofloxacin", "Ciprofloxacin",
                                                   "Novobiocin", "Coumermycin",
                                                   "CN-DM-861", "AR351", "untreated"))

pdf(file="boxplots_all.pdf", paper="a4")
for (i in 1: nrow(peaktable_iso_anno)){
  my_name <- peaktable_iso_anno$anno[i]
  df0 <- peaktable_iso_anno[i, 18:41]
  
  df1 <- data.frame(compound = class.factor, peak_area = as.numeric(df0[1,]))
  
  
  p10 <- ggplot(df1, aes(x = compound, y = peak_area)) +
    geom_boxplot(colour = "black", fill = "lightgrey") +
    geom_jitter(width = 0.25) +
    ggtitle(my_name)
  p10 <- p10 + theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
  
  print(p10)
  
}
dev.off()


# heatmap only using annotated main features (monoisotopic ions), max-scaled and centred

pt_filtered <- peaktable_deiso_anno %>% filter(main_peak == 1)

my_mat <- pt_filtered[,18:41]

# scaling by devision with max(x)
for (i in 1:nrow(my_mat)) {
  my_mat[i,] <- my_mat[i,] / max(my_mat[i,])
}

my_mat <- scale(t(my_mat), scale = F, center = T)
my_mat <- as.data.frame(t(my_mat))
rownames(my_mat) <- pt_filtered$anno

my_color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100)
pheatmap(my_mat, clustering_method = "ward.D", show_rownames = T, fontsize_row = 2, fontsize_col = 8,
         scale = "none", color = my_color)


# heatmap only using annotated main features (monoisotopic ions), median, then max-scaled and centred
pt_filtered <- peaktable_annotated2 %>% filter(main_peak == 1)

my_mat <- pt_filtered[,17:40]

my_medians<-apply(my_mat, 1, function(x) tapply(x, classvector, mean))
my_clumat <- t(my_medians)

# scaling by devision with max(x)
for (i in 1:nrow(my_clumat)) {
  my_clumat[i,] <- my_clumat[i,] / max(my_clumat[i,])
}

my_clumat <- scale(t(my_clumat), scale = F, center = T)
my_clumat <- as.data.frame(t(my_clumat))
rownames(my_clumat) <- pt_filtered$anno



my_color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100)
pheatmap(my_clumat, clustering_method = "ward.D", show_rownames = T, fontsize_row = 5, fontsize_col = 8,
         scale = "none", color = my_color)


# Fold change and Welch´s t-test (unequal variances)

# Create the fold change function
D2 <- as.data.frame(peaktable_deiso_anno[,18:41])
classvector <- class.vector


means<-apply(D2, 1, function(x) tapply(x, classvector, mean))
means <- t(means)

AR351_mean <- means[,"AR351"];control_mean <- means[,"untreated"]
log2FC_AR351 <- log2(AR351_mean/control_mean)

Cipro_mean <- means[,"Ciprofloxacin"]
log2FC_Cipro <- log2(Cipro_mean/control_mean)

Coumer_mean <- means[,"Coumermycin"]
log2FC_Coumer <- log2(Coumer_mean/control_mean)

CN_DM_861_mean <- means[,"CN-DM-861"]
log2FC_CN_DM_861 <- log2(CN_DM_861_mean/control_mean)

Levo_mean <- means[,"Levofloxacin"]
log2FC_Levo <- log2(Levo_mean/control_mean)

Lome_mean <- means[,"Lomefloxacin"]
log2FC_Lome <- log2(Lome_mean/control_mean)

Novo_mean <- means[,"Novobiocin"]
log2FC_Novo <- log2(Novo_mean/control_mean)


# Welch´s t-test (unequal variances)

pval_AR351 <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,18:20]))
  pval_AR351 <- c(pval_AR351, t_result$p.value)
}
padj_AR351 <- p.adjust(pval_AR351, method = "fdr")

pval_Cipro <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,21:23]))
  pval_Cipro <- c(pval_Cipro, t_result$p.value)
}
padj_Cipro <- p.adjust(pval_Cipro, method = "fdr")

pval_Coumer <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,24:26]))
  pval_Coumer <- c(pval_Coumer, t_result$p.value)
}
padj_Coumer <- p.adjust(pval_Coumer, method = "fdr")

pval_CN_DM_861 <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,27:29]))
  pval_CN_DM_861 <- c(pval_CN_DM_861, t_result$p.value)
}
padj_CN_DM_861 <- p.adjust(pval_CN_DM_861, method = "fdr")

pval_Levo <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,30:32]))
  pval_Levo <- c(pval_Levo, t_result$p.value)
}
padj_Levo <- p.adjust(pval_Levo, method = "fdr")

pval_Lome <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,33:35]))
  pval_Lome <- c(pval_Lome, t_result$p.value)
}
padj_Lome <- p.adjust(pval_Lome, method = "fdr")

pval_Novo <- c()
for(i in 1:nrow(peaktable_deiso_anno)) {
  t_result <- t.test(as.numeric(peaktable_deiso_anno[i,39:41]), as.numeric(peaktable_deiso_anno[i,36:38]))
  pval_Novo <- c(pval_Novo, t_result$p.value)
}
padj_Novo <- p.adjust(pval_Novo, method = "fdr")

fc.res <- data.frame(cbind(control_mean, AR351_mean, log2FC_AR351, pval_AR351, padj_AR351,
                           Cipro_mean, log2FC_Cipro, pval_Cipro, padj_Cipro,
                           Coumer_mean, log2FC_Coumer, pval_Coumer, padj_Coumer,
                           CN_DM_861_mean, log2FC_CN_DM_861, pval_CN_DM_861, padj_CN_DM_861,
                           Levo_mean, log2FC_Levo, pval_Levo, padj_Levo,
                           Lome_mean, log2FC_Lome, pval_Lome, padj_Lome,
                           Novo_mean, log2FC_Novo, pval_Novo, padj_Novo))


fc.results.rounded <- apply(fc.res, 2, round, digits = 3)
peaktable_FC <- cbind(peaktable_deiso_anno, fc.res)

write.csv(peaktable_FC, "peaktable_FC.csv")

# determine regulated features
fc.res %>% filter(log2FC_Cipro < -log2(1.5)) %>% filter(pval_Cipro < 0.05) %>% dim()
fc.res %>% filter(log2FC_Cipro > log2(1.5)) %>% filter(pval_Cipro < 0.05) %>% dim()
fc.res %>% filter(log2FC_Levo < -log2(1.5)) %>% filter(pval_Levo < 0.05) %>% dim()
fc.res %>% filter(log2FC_Levo > log2(1.5)) %>% filter(pval_Levo < 0.05) %>% dim()
fc.res %>% filter(log2FC_Lome < -log2(1.5)) %>% filter(pval_Lome < 0.05) %>% dim()
fc.res %>% filter(log2FC_Lome > log2(1.5)) %>% filter(pval_Lome < 0.05) %>% dim()
fc.res %>% filter(log2FC_Novo < -log2(1.5)) %>% filter(pval_Novo < 0.05) %>% dim()
fc.res %>% filter(log2FC_Novo > log2(1.5)) %>% filter(pval_Novo < 0.05) %>% dim()
fc.res %>% filter(log2FC_Coumer < -log2(1.5)) %>% filter(pval_Coumer < 0.05) %>% dim()
fc.res %>% filter(log2FC_Coumer > log2(1.5)) %>% filter(pval_Coumer < 0.05) %>% dim()
fc.res %>% filter(log2FC_CN_DM_861 < -log2(1.5)) %>% filter(pval_CN_DM_861 < 0.05) %>% dim()
fc.res %>% filter(log2FC_CN_DM_861 > log2(1.5)) %>% filter(pval_CN_DM_861 < 0.05) %>% dim()
fc.res %>% filter(log2FC_AR351 < -log2(1.5)) %>% filter(pval_AR351 < 0.05) %>% dim()
fc.res %>% filter(log2FC_AR351 > log2(1.5)) %>% filter(pval_AR351 < 0.05) %>% dim()

# Volcano plots

library(EnhancedVolcano)

# AR351 vs. control

df_AR351 <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_AR351),
                       pval = as.numeric(peaktable_FC$pval_AR351),
                       anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_AR351 <- df_AR351[!is.na(df_AR351$log2FC),]
df_AR351 <- df_AR351[df_AR351$log2FC < 10 & df_AR351$log2FC > -10, ]

EnhancedVolcano(toptable = df_AR351, x = "log2FC", y = "pval", lab = df_AR351$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "AR351 vs. Control")

# cipro vs. control
df_Cipro <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_Cipro),
                       pval = as.numeric(peaktable_FC$pval_Cipro),
                       anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_Cipro <- df_Cipro[!is.na(df_Cipro$log2FC),]
df_Cipro <- df_Cipro[df_Cipro$log2FC < 10 & df_Cipro$log2FC > -10, ]

EnhancedVolcano(toptable = df_Cipro, x = "log2FC", y = "pval", lab = df_Cipro$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "Cipro vs. Control")

# Coumer vs. control
df_Coumer <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_Coumer),
                        pval = as.numeric(peaktable_FC$pval_Coumer),
                        anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_Coumer <- df_Coumer[!is.na(df_Coumer$log2FC),]
df_Coumer <- df_Coumer[df_Coumer$log2FC < 10 & df_Coumer$log2FC > -10, ]

EnhancedVolcano(toptable = df_Coumer, x = "log2FC", y = "pval", lab = df_Coumer$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "Coumer vs. Control")

#GT3 vs. control
df_GT3 <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_GT3),
                     pval = as.numeric(peaktable_FC$pval_GT3),
                     anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_GT3 <- df_GT3[!is.na(df_GT3$log2FC),]
df_GT3 <- df_GT3[df_GT3$log2FC < 10 & df_GT3$log2FC > -10, ]

EnhancedVolcano(toptable = df_GT3, x = "log2FC", y = "pval", lab = df_GT3$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "GT3 vs. Control")

# Levo vs. control
df_Levo <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_Levo),
                      pval = as.numeric(peaktable_FC$pval_Levo),
                      anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_Levo <- df_Levo[!is.na(df_Levo$log2FC),]
df_Levo <- df_Levo[df_Levo$log2FC < 10 & df_Levo$log2FC > -10, ]

EnhancedVolcano(toptable = df_Levo, x = "log2FC", y = "pval", lab = df_Levo$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "Levo vs. Control")


# Lome vs. control
df_Lome <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_Lome),
                      pval = as.numeric(peaktable_FC$pval_Lome),
                      anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_Lome <- df_Lome[!is.na(df_Lome$log2FC),]
df_Lome <- df_Lome[df_Lome$log2FC < 10 & df_Lome$log2FC > -10, ]

EnhancedVolcano(toptable = df_Lome, x = "log2FC", y = "pval", lab = df_Lome$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "Lome vs. Control")


# novo vs. control
df_novo <- data.frame(log2FC = as.numeric(peaktable_FC$log2FC_Novo),
                      pval = as.numeric(peaktable_FC$pval_Novo),
                      anno = peaktable_FC$anno)

# log2FC >10 < -10 NaN have to be removed 

df_novo <- df_novo[!is.na(df_novo$log2FC),]
df_novo <- df_novo[df_novo$log2FC < 10 & df_novo$log2FC > -10, ]

EnhancedVolcano(toptable = df_novo, x = "log2FC", y = "pval", lab = df_novo$anno, ylim = c(0, 5),
                transcriptLabSize = 2,
                title = "Novo vs. Control")


