install.packages("SKAT")
library(SKAT)

TEST = read.table("merged.SetID", header=F)
test = TEST[-which(duplicated(TEST$V2)),]
write.table(test, "merged.SetID", quote = FALSE, sep = '\t', row.names=F, col.names=F)

File.Bed <- ("merged.bed")
File.Bim <- ("merged.bim")
File.Fam <- ("merged.fam")
File.Cov <- ("merged.cov")    
File.SetID <- ("merged.SetID")
File.SSD <- ("merged.SSD")
File.Info <- ("merged.info")
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
SSD.INFO <- Open_SSD(File.SSD, File.Info)

phenotypes = read.table(File.Cov, header = T)
y = phenotypes$Phenotype
Gender = phenotypes$Sex
out.skat <- SKAT.SSD.All(SSD.INFO, SKAT_Null_Model(y~1, out_type="C"))
out.skato <- SKAT.SSD.All(SSD.INFO, SKAT_Null_Model(y~1+Gender, out_type="C"), method="optimal")
out.burden <- SKAT.SSD.All(SSD.INFO, SKAT_Null_Model(y~1, out_type="C"), r.corr=1)

skat <- out.skat$results
pval_skat <- skat$P.value
pval_skat.bon <- p.adjust(pval_skat, method = "bonferroni")
skato <-out.skato$results
pval_skato <- skato$P.value
pval_skato.bon <- p.adjust(pval_skato, method = "bonferroni")
burden <-out.burden$results
pval_burden <- burden$P.value
pval_burden.bon <- p.adjust(pval_burden, method = "bonferroni")

data <- cbind(skat,pval_skat.bon,skato,pval_skato.bon,burden,pval_burden.bon)
write.table(data,"SKAT_results.txt",col.names=T,row.names=F,quote=F,sep="\t")

