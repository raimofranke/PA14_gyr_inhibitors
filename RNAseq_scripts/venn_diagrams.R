library(VennDiagram)


fluoroquinolones_down <- results_all[results_all$log2FC_fluoroquinolone < -1 & results_all$FDR_fluoroquinolone < 0.05,]$ID
aminocoumarins_down <- results_all[results_all$log2FC_aminocoumarin < -1 & results_all$FDR_aminocoumarin < 0.05,]$ID
cystobactamid_down <- results_all[results_all$log2FC_cystobactamid < -1 & results_all$FDR_cystobactamid < 0.05,]$ID

fluoroquinolones_up <- results_all[results_all$log2FC_fluoroquinolone > 1 & results_all$FDR_fluoroquinolone < 0.05,]$ID
aminocoumarins_up <- results_all[results_all$log2FC_aminocoumarin > 1 & results_all$FDR_aminocoumarin < 0.05,]$ID
cystobactamid_up <- results_all[results_all$log2FC_cystobactamid > 1 & results_all$FDR_cystobactamid < 0.05,]$ID

fluoroquinolones_all <- c(as.character(fluoroquinolones_down), as.character(fluoroquinolones_up))
aminocoumarins_all <- c(as.character(aminocoumarins_down), as.character(aminocoumarins_up))
cystobactamid_all <- c(as.character(cystobactamid_down), as.character(cystobactamid_up))


venn.diagram(list(
  fluoroquinolones = fluoroquinolones_down,
  aminocoumarins = aminocoumarins_down,
  cystobactamids = cystobactamid_down),

file = "venn_down.tiff")

venn.diagram(list(
  fluoroquinolones = fluoroquinolones_up,
  aminocoumarins = aminocoumarins_up,
  cystobactamids = cystobactamid_up),
  file = "venn_up.tiff")

venn.diagram(list(
  fluoroquinolones = fluoroquinolones_all,
  aminocoumarins = aminocoumarins_all,
  cystobactamids = cystobactamid_all),
  file = "venn_all.tiff")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
#myCol <- brewer.pal(3, "RdYlBu")
myCol <- c("#4BB446", "#AF46B4", "#4682B4")

# Chart
venn.diagram(x = list(fluoroquinolones = fluoroquinolones_all,
    aminocoumarins = aminocoumarins_all,
    cystobactamids = cystobactamid_all),
  file = "venn_all.tiff",
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  height = 960 , 
  width = 1040 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


# calculate overlap
myoverlap_down <- calculate.overlap(list(fluoroquinolones = fluoroquinolones_down,
  aminocoumarins = aminocoumarins_down,
  cystobactamids = cystobactamid_down))

myoverlap_up <- calculate.overlap(list(fluoroquinolones = fluoroquinolones_up,
                                       aminocoumarins = aminocoumarins_up,
                                       cystobactamids = cystobactamid_up))


# filter results table
# downregulate only FQ
my_ind <- c()
for (i in 1:length(myoverlap_down$a1)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a1[i], results_all$ID))
}

mytable_down_FQonly <- results_all[my_ind,]
write.csv(mytable_down_FQonly, file = "mytable_down_FQ_only.csv")

# downregulate FQ and AC
my_ind <- c()
for (i in 1:length(myoverlap_down$a2)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a2[i], results_all$ID))
}

mytable_down_FQ_AC <- results_all[my_ind,]
write.csv(mytable_down_FQ_AC, file = "mytable_down_FQ_AC.csv")

# downregulate only AC
my_ind <- c()
for (i in 1:length(myoverlap_down$a3)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a3[i], results_all$ID))
}

mytable_down_AConly <- results_all[my_ind,]
write.csv(mytable_down_AConly, file = "mytable_down_AC_only.csv")

# downregulation only CY
my_ind <- c()
for (i in 1:length(myoverlap_down$a7)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a7[i], results_all$ID))
}

mytable_down_CYonly <- results_all[my_ind,]
write.csv(mytable_down_CYonly, file = "mytable_down_CY_only.csv")

## upregulation
# upregulation only FQ
my_ind <- c()
for (i in 1:length(myoverlap_up$a1)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a1[i], results_all$ID))
}

mytable_up_FQonly <- results_all[my_ind,]
write.csv(mytable_up_FQonly, file = "mytable_up_FQ_only.csv")

# upregulation only FQ and AC
my_ind <- c()
for (i in 1:length(myoverlap_up$a2)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a2[i], results_all$ID))
}

mytable_up_FQ_AC <- results_all[my_ind,]
write.csv(mytable_up_FQ_AC, file = "mytable_up_FQ_AC.csv")

# upregulation only AC
my_ind <- c()
for (i in 1:length(myoverlap_up$a3)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a3[i], results_all$ID))
}

mytable_up_AConly <- results_all[my_ind,]
write.csv(mytable_up_AConly, file = "mytable_up_AC_only.csv")


# upregulation FQ and CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a4)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a4[i], results_all$ID))
}

mytable_up_FQ_CY <- results_all[my_ind,]
write.csv(mytable_up_FQ_CY, file = "mytable_up_FQ_CY.csv")

# upregulation FQ and AC and CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a5)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a5[i], results_all$ID))
}

mytable_up_FQ_AC_CY <- results_all[my_ind,]
write.csv(mytable_up_FQ_AC_CY, file = "mytable_up_FQ_AC_CY.csv")

# upregulation AC and CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a6)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a6[i], results_all$ID))
}

mytable_up_AC_CY <- results_all[my_ind,]
write.csv(mytable_up_AC_CY, file = "mytable_up_AC_CY.csv")

# upregulation only CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a7)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a7[i], results_all$ID))
}

mytable_up_CYonly <- results_all[my_ind,]
write.csv(mytable_up_CYonly, file = "mytable_up_CY_only.csv")







# Filter for opposing regulation
FQ_up_AC_down <- results_all[results_all$log2FC_fluoroquinolone > 1 & results_all$FDR_fluoroquinolone < 0.05 &
                               results_all$log2FC_aminocoumarin < -1 & results_all$FDR_aminocoumarin < 0.05,]$ID


FQ_down_AC_up <- results_all[results_all$log2FC_fluoroquinolone < -1 & results_all$FDR_fluoroquinolone < 0.05 &
                               results_all$log2FC_aminocoumarin > 1 & results_all$FDR_aminocoumarin < 0.05,]$ID


FQ_up_CY_down <- results_all[results_all$log2FC_fluoroquinolone > 1 & results_all$FDR_fluoroquinolone < 0.05 &
                               results_all$log2FC_cystobactamid < -1 & results_all$FDR_cystobactamid < 0.05,]$ID

FQ_down_CY_up <- results_all[results_all$log2FC_fluoroquinolone < -1 & results_all$FDR_fluoroquinolone < 0.05 &
                               results_all$log2FC_cystobactamid > 1 & results_all$FDR_cystobactamid < 0.05,]$ID

AC_up_CY_down <- results_all[results_all$log2FC_aminocoumarin > 1 & results_all$FDR_aminocoumarin < 0.05 &
                               results_all$log2FC_cystobactamid < -1 & results_all$FDR_cystobactamid < 0.05,]$ID

AC_down_CY_up <- results_all[results_all$log2FC_aminocoumarin < -1 & results_all$FDR_aminocoumarin < 0.05 &
                               results_all$log2FC_cystobactamid > 1 & results_all$FDR_cystobactamid < 0.05,]$ID


  