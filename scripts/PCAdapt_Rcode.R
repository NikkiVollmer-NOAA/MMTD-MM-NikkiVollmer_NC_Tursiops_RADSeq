#https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
#install the package and load it
install.packages("pcadapt")
library(pcadapt)

#to convert structure file format into pcadapt-accepted format can use the program "LEA". To install LEA need to install/
#load devtools
library(devtools)

#use devtools to install LEA and then load it
devtools::install_github("bcm-uga/LEA")
library(LEA)

#set wd to where your structure file is to convert
setwd("C:/Users/nicole.vollmer/Desktop")

#use struct2geno function to convert structure format to 'lfmm' format (latent factor mixed models). My structure (.txt) file contains all of 
#my 5767snp individs with pop info for each and sorted with all individs in each pop clumped together so I can do the poplist 
#function below. Also any non-polymorphic loci (across all pops) will be flagged and you will have to remove them for conversion
#to work. Missing data must be -9
struct2geno("5760pops.txt", ploidy = 2, FORMAT = 2, extra.row = 1, extra.column = 2)

#now that you have data in lfmm format you can use it in pcadapt (also use same input file for OutFLANK)
filename <- read.pcadapt("6475pops.txt.lfmm", type = "lfmm")

#want to start off testing a large number of PC's (aka K), manual says 20 is large, i went with 100
x <- pcadapt(input = filename, K = 100)
#plots scree plot, and recommend to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell's rule)
plot(x, option = "screeplot")

#alternatively, can use score plot to choose # PCs, here I am assigning the 1st 15 individs as the darkgreen pop, 
#the next 14 as lightgreen and so on.
poplist.int <- c(rep("orange", 84), rep("red", 21), rep("lightgreen", 14), rep("darkgreen", 15))
#double check it looks correct
poplist.int
#you can plot different PCs and when you stop seeing the structure you expect that is where your PC cutoff is
plot(x, option = "scores", pop = poplist.int)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.int)

#once you have your PC # can run the pcadapt function again with it
x <- pcadapt(filename, K = 4)
#the following shows numerical quantities obtained after performing a PCA on the genotype matrix
summary(x)
#to look at any category in summary do x$name, e.g., x$pvalues; and to export a table do write.table(x$pvalues, "pvalues.txt", sep="\t")
x$pvalues
write.table(x$pvalues, "pvalues.txt", sep="\t")

#A Manhattan plot displays ???log10 of the p-values
plot(x , option = "manhattan")
#Can also check the expected uniform distribution of the p-values using a Q-Q plot
plot(x, option = "qqplot")
#An histogram of p-values confirms that most of the p-values follow an uniform distribution
#The excess of small p-values indicates the presence of outliers.
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#The presence of outliers is also visible when plotting a histogram of the test statistic Dj
plot(x, option = "stat.distribution")

#To provide a list of outliers and choose a cutoff for outlier detection, 
#there are several methods that are listed below from the less conservative one to the more conservative one

#A: The R package qvalue, transforms p-values into q-values. see here about install for R 4.1 http://www.bioconductor.org/packages/release/bioc/html/qvalue.html 
#load the package
library(qvalue)
#SNPs with q-values less than ?? will be considered as outliers with an expected false discovery rate bounded by ??. 
#The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs.
#here set a FDR cutoff of 5%
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
outliers
#to export list of qvalues for all loci
write.table(qval, "qval_test.txt", sep="\t")


#B: Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
outliers

#C: Bonferroni Correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
outliers

