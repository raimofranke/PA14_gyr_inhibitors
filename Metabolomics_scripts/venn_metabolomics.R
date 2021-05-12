library(VennDiagram)

fq_index_down <- peaktable_FC$log2FC_Cipro < -log2(1.5) & peaktable_FC$pval_Cipro < 0.05 &
  peaktable_FC$log2FC_Levo < -log2(1.5) & peaktable_FC$pval_Levo < 0.05 &
  peaktable_FC$log2FC_Lome < -log2(1.5) & peaktable_FC$pval_Lome < 0.05

ac_index_down <- peaktable_FC$log2FC_Novo < -log2(1.5) & peaktable_FC$pval_Novo < 0.05 &
  peaktable_FC$log2FC_Coumer < -log2(1.5) & peaktable_FC$pval_Coumer < 0.05

cy_index_down <- peaktable_FC$log2FC_CN_DM_861 < -log2(1.5) & peaktable_FC$pval_CN_DM_861 < 0.05 &
  peaktable_FC$log2FC_AR351 < -log2(1.5) & peaktable_FC$pval_AR351 < 0.05

fq_index_up <- peaktable_FC$log2FC_Cipro > log2(1.5) & peaktable_FC$pval_Cipro < 0.05 &
  peaktable_FC$log2FC_Levo > log2(1.5) & peaktable_FC$pval_Levo < 0.05 &
  peaktable_FC$log2FC_Lome > log2(1.5) & peaktable_FC$pval_Lome < 0.05

ac_index_up <- peaktable_FC$log2FC_Novo > log2(1.5) & peaktable_FC$pval_Novo < 0.05 &
  peaktable_FC$log2FC_Coumer > log2(1.5) & peaktable_FC$pval_Coumer < 0.05

cy_index_up <- peaktable_FC$log2FC_CN_DM_861 > log2(1.5) & peaktable_FC$pval_CN_DM_861 < 0.05 &
  peaktable_FC$log2FC_AR351 > log2(1.5) & peaktable_FC$pval_AR351 < 0.05


fluoroquinolones_down <- peaktable_FC[fq_index_down,]
aminocoumarins_down <- peaktable_FC[ac_index_down,]
cystobactamids_down <- peaktable_FC[cy_index_down,]

fluoroquinolones_up <- peaktable_FC[fq_index_up,]
aminocoumarins_up <- peaktable_FC[ac_index_up,]
cystobactamids_up <- peaktable_FC[cy_index_up,]
write.csv(fluoroquinolones_down, "FQdown_1_5.csv")
write.csv(fluoroquinolones_up, "FQup_1_5.csv")
write.csv(aminocoumarins_down, "ACdown_1_5.csv")
write.csv(aminocoumarins_up, "ACup_1_5.csv")
write.csv(cystobactamids_down, "CYdown_1_5.csv")
write.csv(cystobactamids_up, "CYup_1_5.csv")


fluoroquinolones_all <- c(as.character(fluoroquinolones_down$compound), as.character(fluoroquinolones_up$compound))
aminocoumarins_all <- c(as.character(aminocoumarins_down$compound), as.character(aminocoumarins_up$compound))
cystobactamids_all <- c(as.character(cystobactamids_down$compound), as.character(cystobactamids_up$compound))


venn.diagram(list(
  fluoroquinolones = fluoroquinolones_down$compound,
  aminocoumarins = aminocoumarins_down$compound,
  cystobactamids = cystobactamids_down$compound),
  
  file = "venn_down.tiff")

venn.diagram(list(
  fluoroquinolones = fluoroquinolones_up$compound,
  aminocoumarins = aminocoumarins_up$compound,
  cystobactamidss = cystobactamids_up$compound),
  file = "venn_up.tiff")

venn.diagram(list(
  fluoroquinolones = fluoroquinolones_all,
  aminocoumarins = aminocoumarins_all,
  cystobactamidss = cystobactamids_all),
  file = "venn_all.tiff")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
#myCol <- brewer.pal(3, "RdYlBu")
myCol <- c("#4BB446", "#AF46B4", "#4682B4")

# Chart
venn.diagram(x = list(fluoroquinolones = fluoroquinolones_all,
                      aminocoumarins = aminocoumarins_all,
                      cystobactamids = cystobactamids_all),
             file = "venn_all.tiff",
             output=TRUE,
             
             # Output features
             imagetype="tiff" ,
             height = 960 ,
             width = 1040 ,
             wcompoundth = 520 , 
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
             cat.pos = c(-25, 20, 145),
             cat.dist = c(0.06, 0.06, 0.06),
             cat.fontfamily = "sans",
             rotation = 1
)


# calculate overlap
myoverlap_down <- calculate.overlap(list(fluoroquinolones = fluoroquinolones_down$compound,
                                         aminocoumarins = aminocoumarins_down$compound,
                                         cystobactamidss = cystobactamids_down$compound))

myoverlap_up <- calculate.overlap(list(fluoroquinolones = fluoroquinolones_up$compound,
                                       aminocoumarins = aminocoumarins_up$compound,
                                       cystobactamidss = cystobactamids_up$compound))


# filter results table
# downregulate only FQ
my_ind <- c()
for (i in 1:length(myoverlap_down$a1)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a1[i], peaktable_FC$compound))
}

mytable_down_FQonly <- peaktable_FC[my_ind,]
write.csv(mytable_down_FQonly, file = "mytable_down_FQ_only.csv")

# downregulate FQ and AC
my_ind <- c()
for (i in 1:length(myoverlap_down$a2)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a2[i], peaktable_FC$compound))
}

mytable_down_FQ_AC <- peaktable_FC[my_ind,]
write.csv(mytable_down_FQ_AC, file = "mytable_down_FQ_AC.csv")

# downregulate only AC
my_ind <- c()
for (i in 1:length(myoverlap_down$a3)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a3[i], peaktable_FC$compound))
}

mytable_down_AConly <- peaktable_FC[my_ind,]
write.csv(mytable_down_AConly, file = "mytable_down_AC_only.csv")

# downregulation only CY
my_ind <- c()
for (i in 1:length(myoverlap_down$a7)) {
  my_ind <- c(my_ind, grep(myoverlap_down$a7[i], peaktable_FC$compound))
}

mytable_down_CYonly <- peaktable_FC[my_ind,]
write.csv(mytable_down_CYonly, file = "mytable_down_CY_only.csv")

## upregulation
# upregulation only FQ
my_ind <- c()
for (i in 1:length(myoverlap_up$a1)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a1[i], peaktable_FC$compound))
}

mytable_up_FQonly <- peaktable_FC[my_ind,]
write.csv(mytable_up_FQonly, file = "mytable_up_FQ_only.csv")

# upregulation only FQ and AC
my_ind <- c()
for (i in 1:length(myoverlap_up$a2)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a2[i], peaktable_FC$compound))
}

mytable_up_FQ_AC <- peaktable_FC[my_ind,]
write.csv(mytable_up_FQ_AC, file = "mytable_up_FQ_AC.csv")

# upregulation only AC
my_ind <- c()
for (i in 1:length(myoverlap_up$a5)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a5[i], peaktable_FC$compound))
}

mytable_up_AConly <- peaktable_FC[my_ind,]
write.csv(mytable_up_AConly, file = "mytable_up_AC_only.csv")


# upregulation FQ and CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a4)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a4[i], peaktable_FC$compound))
}

mytable_up_FQ_CY <- peaktable_FC[my_ind,]
write.csv(mytable_up_FQ_CY, file = "mytable_up_FQ_CY.csv")

# upregulation FQ and AC and CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a5)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a5[i], peaktable_FC$compound))
}

mytable_up_FQ_AC_CY <- peaktable_FC[my_ind,]
write.csv(mytable_up_FQ_AC_CY, file = "mytable_up_FQ_AC_CY.csv")

# upregulation AC and CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a6)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a6[i], peaktable_FC$compound))
}

mytable_up_AC_CY <- peaktable_FC[my_ind,]
write.csv(mytable_up_AC_CY, file = "mytable_up_AC_CY.csv")

# upregulation only CY
my_ind <- c()
for (i in 1:length(myoverlap_up$a7)) {
  my_ind <- c(my_ind, grep(myoverlap_up$a7[i], peaktable_FC$compound))
}

mytable_up_CYonly <- peaktable_FC[my_ind,]
write.csv(mytable_up_CYonly, file = "mytable_up_CY_only.csv")







# Filter for opposing regulation
FQ_up_AC_down <- peaktable_FC[fq_index_up & ac_index_down,]$compound


FQ_down_AC_up <- peaktable_FC[peaktable_FC$log2FC_fluoroquinolone < -1 & peaktable_FC$pval_fluoroquinolone < 0.05 &
                               peaktable_FC$log2FC_aminocoumarin > 1 & peaktable_FC$pval_aminocoumarin < 0.05,]$compound


FQ_up_CY_down <- peaktable_FC[peaktable_FC$log2FC_fluoroquinolone > 1 & peaktable_FC$pval_fluoroquinolone < 0.05 &
                               peaktable_FC$log2FC_cystobactamids < -1 & peaktable_FC$pval_cystobactamids < 0.05,]$compound

FQ_down_CY_up <- peaktable_FC[peaktable_FC$log2FC_fluoroquinolone < -1 & peaktable_FC$pval_fluoroquinolone < 0.05 &
                               peaktable_FC$log2FC_cystobactamids > 1 & peaktable_FC$pval_cystobactamids < 0.05,]$compound

AC_up_CY_down <- peaktable_FC[peaktable_FC$log2FC_aminocoumarin > 1 & peaktable_FC$pval_aminocoumarin < 0.05 &
                               peaktable_FC$log2FC_cystobactamids < -1 & peaktable_FC$pval_cystobactamids < 0.05,]$compound

AC_down_CY_up <- peaktable_FC[peaktable_FC$log2FC_aminocoumarin < -1 & peaktable_FC$pval_aminocoumarin < 0.05 &
                               peaktable_FC$log2FC_cystobactamids > 1 & peaktable_FC$pval_cystobactamids < 0.05,]$compound


