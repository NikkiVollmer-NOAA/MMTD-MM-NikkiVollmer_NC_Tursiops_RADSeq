#https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
#https://github.com/whitlock/OutFLANK/blob/master/OutFLANK%20readme.pdf

#load/install all the packages

library(devtools)
library(OutFLANK)
library(qvalue)

if (!("vcfR" %in% installed.packages())){install.packages("vcfR")}
library(vcfR)

#set your working directory
setwd("C:/Users/nicole.vollmer/Desktop/NC_Tursiops/Analyses/OutFLANK")

#can use the same input file created for PCAdapt, need to load it into a dataframe
SNPdata<-read.table("6475pops.txt.lfmm")

#check dataframe to make sure looks ok
head(SNPdata[,1:20])

#need table of locus names in same order as loci in SNPdata, and load those into dataframe
locusNames<-read.table("locusNames.txt")

#need table of population names in same order as samples in SNPdata, and load those into dataframe
popNames<-read.table("popNames.txt")

#this function used the dataframes just loaded to create the appropriate input dataframe to run OutFLANK
FstDataFrame<-MakeDiploidFSTMat(SNPdata, locusNames, popNames)
head(FstDataFrame)

#plots to corrected vs uncorrected Fst, Loci with unusual sample sizes will be outliers in this plot. Outliers should be removed
plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, 
     xlim=c(-0.05,0.5), ylim=c(-0.05,0.5),
     pch=20)
abline(0,1)

#plots heterozygosity vs uncorrected FST
plot(FstDataFrame$He, FstDataFrame$FSTNoCorr, pch=20, col="grey")

#this is the OutFLANK procedure; defaults are 0.05 for right and left trim, 0.1 for threshold for expected heterozygosity,
#0.05 for qthreshold (=FDR 0.05); and for #ofSamples that is your number of pops
out1 <- OutFLANK(FstDataFrame, NumberOfSamples=4)

#this uses the out1 from previous to plot the distribution of FSTNoCorr. Want to make sure the fit is good here.
OutFLANKResultsPlotter(out1, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)

#to zoom in on the right tail to double check fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom = TRUE, RightZoomFraction = 0.05, titletext = NULL)

#plots a histogram of the right tail, it should be flat to indicate a good fit
hist(out1$results$pvaluesRightTail)

#looks at the overview of results
str(out1)

#tells you the # of outliers
sum(out1$results$qvalues<0.10, na.rm=TRUE)


#if you have outliers this plots them highlighted in blue
plot(out1$results$He, out1$results$FST, pch=20, col="grey")
points(out1$results$He[out1$results$qvalues<0.10], out1$results$FST[out1$results$qvalues<0.10], pch=21, col="blue")

#list top candidates/outliers
top_candidates <- out1$results$qvalues<0.10 & out1$results$He>0.1
topcan <- out1$results[top_candidates,]
topcan[order(topcan$LocusName),]



##there is more you can do if you have outliers...see http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html

