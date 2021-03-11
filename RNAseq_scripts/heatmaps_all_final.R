# heatmaps
# read results table
results_all <- read_excel("results_all.xlsx")
results_all <- results_all[order(results_all$`#`),]


# RF-pyocin region
my_ind <- 631:666

# cipro
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_cipro < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_cipro < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_cipro < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), "")
  }
}



FC_pval_cipro <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_cipro <- c(FC_pval_cipro, x_temp)
}

# lome
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_lome < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_lome < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_lome < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), "")
  }
}


my_ind <- 631:666
FC_pval_lome <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_lome <- c(FC_pval_lome, x_temp)
}

# levo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_levo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_levo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_levo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), "")
  }
}


my_ind <- 631:666
FC_pval_levo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_levo <- c(FC_pval_levo, x_temp)
}

# coumer
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_coumer < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_coumer < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_coumer < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), "")
  }
}


my_ind <- 631:666
FC_pval_coumer <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_coumer <- c(FC_pval_coumer, x_temp)
}

# novo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_novo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_novo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_novo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), "")
  }
}


my_ind <- 631:666
FC_pval_novo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_novo <- c(FC_pval_novo, x_temp)
}

# GT3
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_GT3 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), "")
  }
}


my_ind <- 631:666
FC_pval_GT3 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_GT3 <- c(FC_pval_GT3, x_temp)
}

# AR351
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_AR351 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), "")
  }
}


FC_pval_AR351 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_AR351 <- c(FC_pval_AR351, x_temp)
}


my_df <- data.frame(ciprofloxacin = results_all[my_ind,]$log2FC_cipro,
                    levofloxacin = results_all[my_ind,]$log2FC_levo,
                    lomefloxacin = results_all[my_ind,]$log2FC_lome,
                    novobiocin = results_all[my_ind,]$log2FC_novo,
                    coumermycin = results_all[my_ind,]$log2FC_coumer,
                    GT3_043 = results_all[my_ind,]$log2FC_GT3,
                    AR351 = results_all[my_ind,]$log2FC_AR351)


row_anno <- c("PA14_07950_prtN", "PA14_07960_prtR", 
              results_all[my_ind,]$ID[3:20], "PA14_08160_lys", results_all[my_ind,]$ID[22:36])
rownames(my_df) <- row_anno
colnames(my_df) <- c("ciprofloxacin", "levofloxacin", "lomefloxacin",
                     "novobiocin", "coumermycin", "CN-DM-861", "AR351")

my_lables <- matrix(c(FC_pval_cipro, FC_pval_levo, FC_pval_lome,
                      FC_pval_novo, FC_pval_coumer,
                      FC_pval_GT3, FC_pval_AR351), 36, 7)


operon_df = data.frame("annotation" = factor(c("regulation", "regulation", "unknown", "unknown",
                                               "lysis_gene_cass",
                                               rep("R2_pyocin", 15),
                                               "lysis_gene_cass", "lysis_gene_cass", "lysis_gene_cass",
                                               "unknown",
                                               rep("F2_pyocin", 12))))
rownames(operon_df) <- rownames(my_df)
ann_colors <- list(annotation = c(regulation = "yellow", lysis_gene_cass = "green", unknown = "purple",
                                  R2_pyocin= "red", F2_pyocin = "blue"))


breaksList = seq(-7, 7, by = 1)
pheatmap(as.matrix(my_df), scale = "none",
         #color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         border_color = "white", cluster_rows = F, cluster_cols = F, main = "",
         display_numbers = my_lables, number_color = "black", fontsize_number = 10,
         annotation_row = operon_df,
         annotation_colors = ann_colors,
         #gaps_row = 18
)


# SOS regulon (only annotated genes, according to Cirz et al.), arp pathway and S pyocins

# SOS genes
my_ind <- grep("lexA", results_all$Name)
my_ind <- c(my_ind, grep("recA", results_all$Name))
my_ind <- c(my_ind, grep("recN", results_all$Name))
my_ind <- c(my_ind, grep("recX", results_all$Name))
my_ind <- c(my_ind, grep("dinG", results_all$Name))
my_ind <- c(my_ind, grep("dnaE2", results_all$Name))
my_ind <- c(my_ind, grep("PA0069", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA0670", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA0671", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA0922", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA1044", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA2288", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA3008", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA3413", results_all$PAO1_ID))
my_ind <- c(my_ind, grep("PA3414", results_all$PAO1_ID))


# S pyocin
my_ind <- c(my_ind, grep("pyoS3I", results_all$Name))
my_ind <- c(my_ind, grep("pyoS3A", results_all$Name))
my_ind <- c(my_ind, grep("PA14_59220", results_all$ID)) # pyoS5

# alp pathway
my_ind <- c(my_ind, grep("PA0906", results_all$PAO1_ID)) # alpR
my_ind <- c(my_ind, grep("PA0907", results_all$PAO1_ID)) # alpA
my_ind <- c(my_ind, grep("PA0908", results_all$PAO1_ID)) # alpB
my_ind <- c(my_ind, grep("PA0909", results_all$PAO1_ID)) # alpC
my_ind <- c(my_ind, grep("PA0910", results_all$PAO1_ID)) # alpD
my_ind <- c(my_ind, grep("PA0911", results_all$PAO1_ID)) # alpE

my_rownames <- c("PA3007, lexA", "PA3617, recA", "PA4763, recN", "PA3616, recX", "PA1045, dinG", "PA0669, dnaE2",
                 "PA0069, phl", "PA0670, imuB", "PA0671, sulA2", "PA0922", "PA1044", "PA2288",
                 "PA3008, sulA", "PA3413, yebG", "PA3414",
                 "pyoS3I", "pyoS3A", "pyoS5",
                 "alpR", "alpA", "alpB", "alpC", "alpD", "alpE")


# cipro
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_cipro < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_cipro < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_cipro < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), "")
  }
}



FC_pval_cipro <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_cipro <- c(FC_pval_cipro, x_temp)
}

# lome
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_lome < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_lome < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_lome < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), "")
  }
}

FC_pval_lome <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_lome <- c(FC_pval_lome, x_temp)
}

# levo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_levo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_levo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_levo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), "")
  }
}

FC_pval_levo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_levo <- c(FC_pval_levo, x_temp)
}

# coumer
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_coumer < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_coumer < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_coumer < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), "")
  }
}

FC_pval_coumer <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_coumer <- c(FC_pval_coumer, x_temp)
}

# novo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_novo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_novo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_novo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), "")
  }
}

FC_pval_novo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_novo <- c(FC_pval_novo, x_temp)
}

# GT3
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_GT3 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), "")
  }
}

FC_pval_GT3 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_GT3 <- c(FC_pval_GT3, x_temp)
}

# AR351
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_AR351 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), "")
  }
}

FC_pval_AR351 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_AR351 <- c(FC_pval_AR351, x_temp)
}


my_df <- data.frame(ciprofloxacin = results_all[my_ind,]$log2FC_cipro,
                    levofloxacin = results_all[my_ind,]$log2FC_levo,
                    lomefloxacin = results_all[my_ind,]$log2FC_lome,
                    novobiocin = results_all[my_ind,]$log2FC_novo,
                    coumermycin = results_all[my_ind,]$log2FC_coumer,
                    GT3_043 = results_all[my_ind,]$log2FC_GT3,
                    AR351 = results_all[my_ind,]$log2FC_AR351)


colnames(my_df) <- c("ciprofloxacin", "levofloxacin", "lomefloxacin",
                     "novobiocin", "coumermycin", "CN-DM-861", "AR351")
rownames(my_df) <- my_rownames

my_lables <- matrix(c(FC_pval_cipro, FC_pval_levo, FC_pval_lome,
                      FC_pval_novo, FC_pval_coumer,
                      FC_pval_GT3, FC_pval_AR351), 24, 7)


operon_df = data.frame("annotation" = factor(c(rep("SOS_genes", 15),
                                               rep("S_pyocins", 3),
                                               rep("alp_pathway", 6))))
rownames(operon_df) <- rownames(my_df)
ann_colors <- list(annotation = c(SOS_genes = "green", S_pyocins = "red", alp_pathway = "blue"))


breaksList = seq(-7, 7, by = 1)
pheatmap(as.matrix(my_df), scale = "none",
         #color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         border_color = "white", cluster_rows = F, cluster_cols = F, main = "",
         display_numbers = my_lables, number_color = "black", fontsize_number = 10,
         annotation_row = operon_df,
         annotation_colors = ann_colors,
         gaps_row = c(15, 18)
)


# phenazine genes

my_ind <- grep("phzA1", results_all$Name)
my_ind <- c(my_ind, grep("phzB1", results_all$Name))
my_ind <- c(my_ind, grep("phzC1", results_all$Name))
my_ind <- c(my_ind, grep("phzD1", results_all$Name))
my_ind <- c(my_ind, grep("phzE1", results_all$Name))
my_ind <- c(my_ind, grep("phzF1", results_all$Name))
my_ind <- c(my_ind, grep("phzG1", results_all$Name))
my_ind <- c(my_ind, grep("phzA2", results_all$Name))
my_ind <- c(my_ind, grep("phzB2", results_all$Name))
my_ind <- c(my_ind, grep("phzC2", results_all$Name))
my_ind <- c(my_ind, grep("phzD2", results_all$Name))
my_ind <- c(my_ind, grep("phzE2", results_all$Name))
my_ind <- c(my_ind, grep("phzF2", results_all$Name))
my_ind <- c(my_ind, grep("phzG2", results_all$Name))

# rhamnolipids
my_ind <- c(my_ind, grep("rhlA", results_all$Name))
my_ind <- c(my_ind, grep("rhlB", results_all$Name))
my_ind <- c(my_ind, grep("rhlC", results_all$Name))

# cipro
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_cipro < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_cipro < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_cipro < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), "")
  }
}



FC_pval_cipro <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_cipro <- c(FC_pval_cipro, x_temp)
}

# lome
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_lome < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_lome < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_lome < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), "")
  }
}

FC_pval_lome <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_lome <- c(FC_pval_lome, x_temp)
}

# levo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_levo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_levo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_levo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), "")
  }
}

FC_pval_levo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_levo <- c(FC_pval_levo, x_temp)
}

# coumer
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_coumer < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_coumer < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_coumer < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), "")
  }
}

FC_pval_coumer <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_coumer <- c(FC_pval_coumer, x_temp)
}

# novo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_novo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_novo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_novo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), "")
  }
}

FC_pval_novo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_novo <- c(FC_pval_novo, x_temp)
}

# GT3
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_GT3 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), "")
  }
}

FC_pval_GT3 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_GT3 <- c(FC_pval_GT3, x_temp)
}

# AR351
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_AR351 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), "")
  }
}

FC_pval_AR351 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_AR351 <- c(FC_pval_AR351, x_temp)
}


my_df <- data.frame(ciprofloxacin = results_all[my_ind,]$log2FC_cipro,
                    levofloxacin = results_all[my_ind,]$log2FC_levo,
                    lomefloxacin = results_all[my_ind,]$log2FC_lome,
                    novobiocin = results_all[my_ind,]$log2FC_novo,
                    coumermycin = results_all[my_ind,]$log2FC_coumer,
                    GT3_043 = results_all[my_ind,]$log2FC_GT3,
                    AR351 = results_all[my_ind,]$log2FC_AR351)



rownames(my_df) <- results_all[my_ind,]$Name
colnames(my_df) <- c("ciprofloxacin", "levofloxacin", "lomefloxacin",
                     "novobiocin", "coumermycin", "CN-DM-861", "AR351")

my_lables <- matrix(c(FC_pval_cipro, FC_pval_levo, FC_pval_lome,
                      FC_pval_novo, FC_pval_coumer,
                      FC_pval_GT3, FC_pval_AR351), 17, 7)


operon_df = data.frame("annotation" = factor(c(rep("phenazine_biosyn", 14),
                                               rep("rhamnolipid_biosyn", 3))))
rownames(operon_df) <- rownames(my_df)
ann_colors <- list(annotation = c(phenazine_biosyn = "blue", rhamnolipid_biosyn = "red"))


breaksList = seq(-7, 7, by = 1)
pheatmap(as.matrix(my_df), scale = "none",
         #color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         border_color = "white", cluster_rows = F, cluster_cols = F, main = "",
         display_numbers = my_lables, number_color = "black", fontsize_number = 10,
         annotation_row = operon_df,
         annotation_colors = ann_colors,
         gaps_row = 14
)

# pyoverdine and pyocheline synthesis
my_ind <- grep("pvd", results_all$Name)
my_ind2 <- grep("pch", results_all$Name)
my_ind <- c(my_ind, my_ind2)


# cipro
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_cipro < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_cipro < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_cipro < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), "")
  }
}



FC_pval_cipro <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_cipro <- c(FC_pval_cipro, x_temp)
}

# lome
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_lome < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_lome < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_lome < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), "")
  }
}

FC_pval_lome <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_lome <- c(FC_pval_lome, x_temp)
}

# levo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_levo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_levo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_levo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), "")
  }
}

FC_pval_levo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_levo <- c(FC_pval_levo, x_temp)
}

# coumer
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_coumer < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_coumer < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_coumer < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), "")
  }
}

FC_pval_coumer <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_coumer <- c(FC_pval_coumer, x_temp)
}

# novo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_novo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_novo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_novo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), "")
  }
}

FC_pval_novo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_novo <- c(FC_pval_novo, x_temp)
}

# GT3
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_GT3 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), "")
  }
}

FC_pval_GT3 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_GT3 <- c(FC_pval_GT3, x_temp)
}

# AR351
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_AR351 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), "")
  }
}

FC_pval_AR351 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_AR351 <- c(FC_pval_AR351, x_temp)
}


my_df <- data.frame(ciprofloxacin = results_all[my_ind,]$log2FC_cipro,
                    levofloxacin = results_all[my_ind,]$log2FC_levo,
                    lomefloxacin = results_all[my_ind,]$log2FC_lome,
                    novobiocin = results_all[my_ind,]$log2FC_novo,
                    coumermycin = results_all[my_ind,]$log2FC_coumer,
                    GT3_043 = results_all[my_ind,]$log2FC_GT3,
                    AR351 = results_all[my_ind,]$log2FC_AR351)



rownames(my_df) <- results_all[my_ind,]$Name

my_lables <- matrix(c(FC_pval_cipro, FC_pval_levo, FC_pval_lome,
                      FC_pval_novo, FC_pval_coumer,
                      FC_pval_GT3, FC_pval_AR351), 22, 7)


operon_df = data.frame("annotation" = factor(c(rep("pyoverdine_biosyn", 13),
                                               rep("pyocheline_biosyn", 9))))
rownames(operon_df) <- rownames(my_df)
ann_colors <- list(annotation = c(pyoverdine_biosyn = "blue", pyocheline_biosyn = "red"))


breaksList = seq(-7, 7, by = 1)
pheatmap(as.matrix(my_df), scale = "none",
         #color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         border_color = "white", cluster_rows = F, cluster_cols = F, main = "",
         display_numbers = my_lables, number_color = "black", fontsize_number = 10,
         annotation_row = operon_df,
         annotation_colors = ann_colors,
         gaps_row = 13
)


# heatshock genes and chaperones

my_ind <- grep("groEL", results_all$Name)
my_ind <- c(my_ind, grep("groES", results_all$Name))
my_ind <- c(my_ind, grep("grpE", results_all$Name))
my_ind <- c(my_ind, grep("dnaK", results_all$Name))
my_ind <- c(my_ind, grep("hslU", results_all$Name))
my_ind <- c(my_ind, grep("clpB", results_all$Name))
my_ind <- c(my_ind, grep("ibpA", results_all$Name))
my_ind <- c(my_ind, grep("hslV", results_all$Name))
my_ind <- c(my_ind, grep("dnaJ", results_all$Name))
my_ind <- c(my_ind, grep("htpG", results_all$Name))


# cipro
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_cipro < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_cipro < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_cipro < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), "")
  }
}



FC_pval_cipro <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_cipro <- c(FC_pval_cipro, x_temp)
}

# lome
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_lome < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_lome < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_lome < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), "")
  }
}

FC_pval_lome <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_lome <- c(FC_pval_lome, x_temp)
}

# levo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_levo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_levo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_levo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), "")
  }
}

FC_pval_levo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_levo <- c(FC_pval_levo, x_temp)
}

# coumer
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_coumer < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_coumer < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_coumer < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), "")
  }
}

FC_pval_coumer <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_coumer <- c(FC_pval_coumer, x_temp)
}

# novo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_novo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_novo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_novo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), "")
  }
}

FC_pval_novo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_novo <- c(FC_pval_novo, x_temp)
}

# GT3
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_GT3 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), "")
  }
}

FC_pval_GT3 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_GT3 <- c(FC_pval_GT3, x_temp)
}

# AR351
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_AR351 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), "")
  }
}

FC_pval_AR351 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_AR351 <- c(FC_pval_AR351, x_temp)
}


my_df <- data.frame(ciprofloxacin = results_all[my_ind,]$log2FC_cipro,
                    levofloxacin = results_all[my_ind,]$log2FC_levo,
                    lomefloxacin = results_all[my_ind,]$log2FC_lome,
                    novobiocin = results_all[my_ind,]$log2FC_novo,
                    coumermycin = results_all[my_ind,]$log2FC_coumer,
                    GT3_043 = results_all[my_ind,]$log2FC_GT3,
                    AR351 = results_all[my_ind,]$log2FC_AR351)



rownames(my_df) <- results_all[my_ind,]$Name

my_lables <- matrix(c(FC_pval_cipro, FC_pval_levo, FC_pval_lome,
                      FC_pval_novo, FC_pval_coumer,
                      FC_pval_GT3, FC_pval_AR351), 10, 7)


breaksList = seq(-7, 7, by = 1)
pheatmap(as.matrix(my_df), scale = "none",
         #color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         border_color = "white", cluster_rows = F, cluster_cols = F, main = "",
         display_numbers = my_lables, number_color = "black", fontsize_number = 10,
         
)

#glucose catabolism

my_ind <- grep("PA14_35250", results_all$ID)
my_ind <- c(my_ind, grep("PA14_35270", results_all$ID))
my_ind <- c(my_ind, grep("PA14_35290", results_all$ID))
my_ind <- c(my_ind, grep("PA14_35300", results_all$ID))
my_ind <- c(my_ind, grep("PA14_35320", results_all$ID))
my_ind <- c(my_ind, grep("PA14_35330", results_all$ID))
my_ind <- c(my_ind, grep("PA14_35340", results_all$ID))
my_ind <- c(my_ind, grep("PA14_35360", results_all$ID))



# cipro
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_cipro < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_cipro < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_cipro < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_cipro, 1), nsmall = 1), "")
  }
}



FC_pval_cipro <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_cipro <- c(FC_pval_cipro, x_temp)
}

# lome
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_lome < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_lome < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_lome < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_lome, 1), nsmall = 1), "")
  }
}

FC_pval_lome <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_lome <- c(FC_pval_lome, x_temp)
}

# levo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_levo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_levo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_levo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_levo, 1), nsmall = 1), "")
  }
}

FC_pval_levo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_levo <- c(FC_pval_levo, x_temp)
}

# coumer
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_coumer < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_coumer < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_coumer < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_coumer, 1), nsmall = 1), "")
  }
}

FC_pval_coumer <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_coumer <- c(FC_pval_coumer, x_temp)
}

# novo
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_novo < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_novo < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_novo < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_novo, 1), nsmall = 1), "")
  }
}

FC_pval_novo <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_novo <- c(FC_pval_novo, x_temp)
}

# GT3
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_GT3 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_GT3 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_GT3, 1), nsmall = 1), "")
  }
}

FC_pval_GT3 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_GT3 <- c(FC_pval_GT3, x_temp)
}

# AR351
annoFCpval <- function(my_ind) {
  if (results_all[my_ind,]$FDR_AR351 < 0.001 ) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ***")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.01) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " ** ")
  } else if(results_all[my_ind,]$FDR_AR351 < 0.05) {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), " *  ")
  } else  {
    paste0(format(round(results_all[my_ind,]$log2FC_AR351, 1), nsmall = 1), "")
  }
}

FC_pval_AR351 <- c()
for (i in 1:length(my_ind)){
  x_temp <-   annoFCpval(my_ind[i])
  FC_pval_AR351 <- c(FC_pval_AR351, x_temp)
}


my_df <- data.frame(ciprofloxacin = results_all[my_ind,]$log2FC_cipro,
                    levofloxacin = results_all[my_ind,]$log2FC_levo,
                    lomefloxacin = results_all[my_ind,]$log2FC_lome,
                    novobiocin = results_all[my_ind,]$log2FC_novo,
                    coumermycin = results_all[my_ind,]$log2FC_coumer,
                    CN_DM_861 = results_all[my_ind,]$log2FC_GT3,
                    AR351 = results_all[my_ind,]$log2FC_AR351)

colnames(my_df) <- c("ciprofloxacin", "levofloxacin", "lomefloxacin",
                     "novobiocin", "coumermycin", "CN-DM-861", "AR351")

rownames(my_df) <- c("PA14_35250", "PA14_35270", "PA14_35290_gad", "PA14_35300", "PA14_35320_kguD",
                     "PA14_35330_kguT", "PA14_35340_kguK", "PA14_35360_kguE")

my_lables <- matrix(c(FC_pval_cipro, FC_pval_levo, FC_pval_lome,
                      FC_pval_novo, FC_pval_coumer,
                      FC_pval_GT3, FC_pval_AR351), 8, 7)


operon_df = data.frame("annotation" = factor(c(rep("gad_operon", 4),
                                               rep("kgu_operon", 4))))
rownames(operon_df) <- rownames(my_df)
ann_colors <- list(annotation = c(gad_operon = "blue", kgu_operon = "red"))


breaksList = seq(-7, 7, by = 1)
pheatmap(as.matrix(my_df), scale = "none",
         #color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         border_color = "white", cluster_rows = F, cluster_cols = F, main = "",
         display_numbers = my_lables, number_color = "black", fontsize_number = 10,
         annotation_row = operon_df,
         annotation_colors = ann_colors,
         gaps_row = 4
)
