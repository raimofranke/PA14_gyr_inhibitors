library(Rsubread)
buildindex(basename="PA14_index",reference="Pseudomonas_aeruginosa_UCBPP-PA14_109.fna")


my_filenames <- dir(pattern = "fastq")
my_filenames_edited <- sub(my_filenames, pattern = ".fastq.gz", replacement = "")

i <- -1
repeat{
  i = i+2
align(index="PA14_index",readfile1=paste0(my_filenames_edited[i],".fastq.gz"),
      readfile2=paste0(my_filenames_edited[i+1],".fastq.gz"),type="dna",
      output_format = "BAM", output_file = paste0(sub(my_filenames_edited[i],
                                pattern = "_R1", replacement = ""), ".bam"),
      minFragLength=20,maxFragLength=600)

if (i==(length(my_filenames_edited)-4)){
  break}
}


# calculate reads per gene

df <- read.delim("Pseudomonas_aeruginosa_UCBPP-PA14_109.gtf", header = F)
library(splitstackshape)
df2 <- cSplit(df, "V9", sep = ";")
df3 <- cSplit(df2, "V9_3", sep = " ")
df4 <- data.frame(GeneID = df3$V9_3_2, Chr = rep("gi|116048575|ref|NC_008463|pseudocap|138", 5979), 
                       Start = df3$V4, End = df3$V5, Strand = df3$V7)


mapping_8C1 <- featureCounts(files = "8C1.bam", annot.ext = df4, isPairedEnd = T, minFragLength = 20)
save(mapping_8C1, file = "mapping_8C1.RData")


library(edgeR)
rawdata <- data.frame(ID = rownames(mapping_1A1$counts), ID_names = df3$V9_4,
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
                      AR351.3 = as.numeric(mapping_7C1$counts),
                      CN861_2.1 = as.numeric(mapping_8B1$counts),
                      CN861_2.2 = as.numeric(mapping_8C1$counts)
                      )

group <- c("untreated", "untreated", "untreated",
               "cipro", "cipro", "cipro",
              "levo", "levo", "levo",
               "lome", "lome", "lome",
               "novo", "novo", "novo",
               "coumer", "coumer", "coumer",
               "GT3", "GT3", "GT3",
               "AR351", "AR351", "AR351",
               "CN861", "CN861")
group <- factor(group)

y <- DGEList(counts=rawdata[,3:28], genes=rawdata[,1:2], group = group)
keep <- filterByExpr(y)
summary(keep) #only 3 genes do not fulfill criterion

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

con <- makeContrasts(coumer - cipro, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
results_coumer_cipro_untreated <- topTags(qlf, n=5979)
write.csv2(results_coumer_cipro_untreated, file = "results_coumer_cipro.csv")
summary(decideTests(qlf))
 RsplotMD(qlf)




logcpm <- cpm(y_disp, prior.count=2, log=TRUE) #for clustering, heatmap
my_cpm <- cpm(y_disp, prior.count=2, log=F)
my_cpm_export <- as.data.frame(my_cpm)
my_cpm_export <- cbind(y_disp$genes$ID, my_cpm_export)
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


