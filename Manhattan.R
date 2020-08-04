library(data.table)
library(qqman

data <- data.frame(fread("adjusted_Re_AD.assoc.logistic",header=T,sep=" "))

manhattan(data, main = "Manhattan Plot", ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"),chrlabs = c(1:20, "X", "Y"))

qq(data$P, main = "Q-Q plot of GWAS p-values",xlim = c(0, 7), ylim = c(0,12), pch = 18, col = "blue4", cex = 1.5, las = 1)

