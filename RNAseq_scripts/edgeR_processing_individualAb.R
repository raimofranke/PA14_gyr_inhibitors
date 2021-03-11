# load libraries
library(edgeR)
library(readxl)

C# load mappings
#setwd("/Volumes/CBIOCOM/RAIMO-FRANKE-RMF/GyrInh_Paper/DA_RNAseq")
setwd("/Users/raimofranke/Dropbox/GyrRNAseq/")

load("mapping_1A1.RData")
load("mapping_1B1.RData")
load("mapping_1C1.RData")
load("mapping_2A1.RData")
load("mapping_2B1.RData")
load("mapping_2C1.RData")
load("mapping_3A1.RData")
load("mapping_3B1.RData")
load("mapping_3C1.RData")
load("mapping_4A1.RData")
load("mapping_4B1.RData")
load("mapping_4C1.RData")
load("mapping_5A1.RData")
load("mapping_5B1.RData")
load("mapping_5C1.RData")
load("mapping_6A1.RData")
load("mapping_6B1.RData")
load("mapping_6C1.RData")
load("mapping_7A1.RData")
load("mapping_7B1.RData")
load("mapping_7C1.RData")
load("mapping_9A1.RData")
load("mapping_9B1.RData")
load("mapping_9C1.RData")


pa14_anno <- read_excel("Pseudomonas aeruginosa_PA14_full annotation.xlsx")


rawdata <- data.frame(ID = pa14_anno$PA14_ID, PAO1_ID = pa14_anno$PAO1_ID, Name = pa14_anno$Name,
                      Product_Name = pa14_anno$Product.Name,
                      GO_BP = pa14_anno$`GO_category:Biological Process`,
                      KEGG_Pathway = pa14_anno$`Pathway:KEGG pathway name`,
                      PseudoCAP_Pathway = pa14_anno$`Pathway:PseudoCAP`,
                      untreated.1 = as.numeric(mapping_9A1$counts),
                      untreated.2 = as.numeric(mapping_9B1$counts),
                      untreated.3 = as.numeric(mapping_9C1$counts),
                      ciprofloxacin.1 = as.numeric(mapping_1A1$counts),
                      ciprofloxacin.2 = as.numeric(mapping_1B1$counts),
                      ciprofloxacin.3 = as.numeric(mapping_1C1$counts),
                      levofloxacin.1 = as.numeric(mapping_2A1$counts),
                      levofloxacin.2 = as.numeric(mapping_2B1$counts),
                      levofloxacin.3 = as.numeric(mapping_2C1$counts),
                      lomefloxacin.1 = as.numeric(mapping_3A1$counts),
                      lomefloxacin.2 = as.numeric(mapping_3B1$counts),
                      lomefloxacin.3 = as.numeric(mapping_3C1$counts),
                      novobiocin.1 = as.numeric(mapping_4A1$counts),
                      novobiocin.2 = as.numeric(mapping_4B1$counts),
                      novobiocin.3 = as.numeric(mapping_4C1$counts),
                      coumermycin.1 = as.numeric(mapping_5A1$counts),
                      coumermycin.2 = as.numeric(mapping_5B1$counts),
                      coumermycin.3 = as.numeric(mapping_5C1$counts),
                      GT3_043.1 = as.numeric(mapping_6A1$counts),
                      GT3_043.2 = as.numeric(mapping_6B1$counts),
                      GT3_043.3 = as.numeric(mapping_6C1$counts),
                      AR351.1 = as.numeric(mapping_7A1$counts),
                      AR351.2 = as.numeric(mapping_7B1$counts),
                      AR351.3 = as.numeric(mapping_7C1$counts)
                      )

group <- c("untreated", "untreated", "untreated",
                    "cipro", "cipro", "cipro",
                    "levo", "levo", "levo",
                    "lome", "lome", "lome",
                    "novo", "novo", "novo",
                    "coumer", "coumer", "coumer",
                    "GT3", "GT3", "GT3",
                    "AR351", "AR351", "AR351")
group <- factor(group)

y <- DGEList(counts=rawdata[,8:31], genes=rawdata[,1:7], group = group)
keep <- filterByExpr(y)
summary(keep) #only 4 genes do not fulfill criterion, no need to filter

plotMD(cpm(y, log=TRUE), column=4)
abline(h=0, col="red", lty=2, lwd=2)
             

#Normalization
y_norm <- calcNormFactors(y)


plotMDS(y_norm)
plotMD(cpm(y_norm, log=TRUE), column=4)
abline(h=0, col="red", lty=2, lwd=2)



# design matrix
design <- model.matrix(~0 + group)
rownames(design) <- colnames(y_norm)
colnames(design) <- levels(group)

#estimate dispersion
y_disp <- estimateDisp(y_norm, design, robust=TRUE)
y_disp$common.dispersion
plotBCV(y_disp)

#differential expression
fit <- glmQLFit(y_disp, design, robust = TRUE)
plotQLDisp(fit)

# cipro vs. untreated
con <- makeContrasts(cipro - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_cipro_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# levo vs. untreated
con <- makeContrasts(levo - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_levo_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# lome vs. untreated
con <- makeContrasts(lome - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_lome_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# novo vs. untreated
con <- makeContrasts(novo - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_novo_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# coumer vs. untreated
con <- makeContrasts(coumer - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_coumer_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# GT3 vs. untreated
con <- makeContrasts(GT3 - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_GT3_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# AR351 vs. untreated
con <- makeContrasts(AR351 - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_AR351_untreated <- topTags(qlf, n=5979)
summary(decideTests(qlf, adjust.method ="fdr", p.value = 0.05, lfc = 1 ))
plotMD(qlf)

# generate big table
all_merged <- merge(results_cipro_untreated$table, results_lome_untreated$table, by = "ID")

results_all <- cbind(all_merged[,1:7],
                     log2FC_cipro = all_merged$logFC.x,
                     FDR_cipro = all_merged$FDR.x,
                     log2FC_lome = all_merged$logFC.y,
                     FDR_lome = all_merged$FDR.y)

all_merged <- merge(results_all, results_levo_untreated$table, by = "ID")
results_all <- cbind(results_all,
                     log2FC_levo = all_merged$logFC,
                     FDR_levo = all_merged$FDR)                     

all_merged <- merge(results_all, results_coumer_untreated$table, by = "ID")
results_all <- cbind(results_all,
                     log2FC_coumer = all_merged$logFC,
                     FDR_coumer = all_merged$FDR)   

all_merged <- merge(results_all, results_novo_untreated$table, by = "ID")
results_all <- cbind(results_all,
                     log2FC_novo = all_merged$logFC,
                     FDR_novo = all_merged$FDR)   

all_merged <- merge(results_all, results_GT3_untreated$table, by = "ID")
results_all <- cbind(results_all,
                     log2FC_GT3 = all_merged$logFC,
                     FDR_GT3 = all_merged$FDR)   

all_merged <- merge(results_all, results_AR351_untreated$table, by = "ID")
results_all <- cbind(results_all,
                     log2FC_AR351 = all_merged$logFC,
                     FDR_AR351 = all_merged$FDR)   



write.csv(results_all, "results_all.csv")


# processing for clustering etc.

logcpm <- cpm(y_disp, prior.count=2, log=TRUE) #for clustering, heatmap
my_cpm <- cpm(y_disp, prior.count=2, log=F)
my_cpm_export <- as.data.frame(t(my_cpm))
colnames(my_cpm_export) <- paste(pa14_anno$PA14_ID, "_", pa14_anno$Product.Name)
my_cpm_export <- cbind(colnames(rawdata[,7:30]), group, my_cpm_export)
write.csv(my_cpm_export, "my_cpm.csv")

distmat <- dist(t(logcpm))
plot(hclust(distmat, method = "ward.D"), hang = -1)

go <- order(lrt$table$PValue)
cpm(y_disp)[o[1:10],]

plotMD(lrt)
abline(h=c(-1, 1), col="blue")

# Annotation with EC numbers from Harvard
# http://pa14.mgh.harvard.edu/cgi-bin/pa14/annotation/downloads.cgi

EC_anno <- read.csv2("PA14_ECnumbers.csv")
rawdata_annotated <- cbind(rawdata[,1:2], EC = NA, rawdata[,3:28])

for(i in 1:length(EC_anno$LocusName)){
ma <- match(EC_anno$LocusName[i], rawdata_annotated$ID)
rawdata_annotated$EC[ma] <- as.character(EC_anno$ECNumber[i])
}

#GO annotation
PA14_GO <- read.csv("gene_ontology_csv.csv")

#pathview
write.csv2(lrt$table, file = "edgeR_results.csv")
edgeR_results <- read.csv2("edgeR_results.csv", dec = ",")

de_test <- data.frame(FC = edgeR_results$logFC)
de_test <- data.frame(FC = rnorm(5979, sd = 3))
rownames(de_test) <- as.character(df4$GeneID)                       
pv.out <- pathview(gene.data =  de_test, pathway.id = "pau02040", species = "pau",
                   gene.idtype = "KEGG", out.suffix = "test")


