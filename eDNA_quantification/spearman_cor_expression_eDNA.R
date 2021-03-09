# calculate spearman correlation between expression values in cpm
# and amount of released eDNA

lys_cpm <- c(8.22, 1549.82, 113.93, 212.67)

# untreated, cipro, novo, CN-DM-861
eDNA <- c(18.86, 261.28, 21.84, 51.48)
cor(lys_cpm, eDNA, method = "spearman")
cor.test(lys_cpm, eDNA, method = "spearman")

my_data <- data.frame(eDNA = eDNA, lys = lys_cpm)

library("ggpubr")
ggscatter(my_data, x = "lys", y = "eDNA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "lys expression [cpm]", ylab = "eDNA [µg/µl]") +
  theme(text = element_text(size=14))
        
