library(dartR)
library(dplyr)

#may also need to run the following, see "Welcome to dartR" notes once you load dartR 
gl.install.vanilla.dartR()

#Import genepop file to genind object
Fourpop <- read.genepop(file = "6251snps_top50CvO.gen", ncode=3L)

#genind to genlight
Fourpop.gl <- gi2gl(Fourpop, verbose=5)
Fourpop.gl_2 <- gl.recalc.metrics(Fourpop.gl)
Fourpop.gl_2$other$loc.metrics <- as.data.frame(Fourpop.gl_2$other$loc.metrics)
Fourpop.gl_3 <- gl.compliance.check(Fourpop.gl_2) #make sure you get no errors with this line

#####HWE#####
#Reports departure of Hardy-Weinberg-Equilibrium for every loci per population or overall, if input file
#has multiple pops need to use subset = "each" so that HWE is calc for each pop separately
#below am doing Exact test (which I think is most similar to Genepop) and including sequential bonferroni (=holm)
HWE<-gl.report.hwe(Fourpop.gl_3, subset ="each", method_sig = "Exact", multi_comp=TRUE, multi_comp_method = "BH", 
                   alpha_val = 0.05, pvalue_type = "midp", verbose=5)
write.csv(HWE,"HWE_BH.csv", row.names = TRUE)


#####LD#####
#Import genepop file to genind object - will need to import each pop separately and run this code for each
Redpop <- read.genepop(file = "5767_red.gen", ncode=3L)

#genind to genlight
Redpop.gl <- gi2gl(Redpop, verbose=5)
Redpop.gl_2 <- gl.recalc.metrics(Redpop.gl)
Redpop.gl_2$other$loc.metrics <- as.data.frame(Redpop.gl_2$other$loc.metrics)
Redpop.gl_3 <- gl.compliance.check(Redpop.gl_2) #make sure you get no errors with this line

#Calculates pairwise population based Linkage Disequilibirum across all loci using the specified number of cores
LD<-gl.report.ld(Redpop.gl_3, save = TRUE, outpath = "./", nchunks = 1, ncores = 4, probar=TRUE, verbose=5)

#the above LD function will produce a dataframe with potentially millions of rows (too many for excel to read)
#so need to sort dataframe by p column (may need to load LDallp.rdata file first)
Red_LDsorted<-LDallp[order(LDallp$p),]

#can check first 6 rows
head(Red_LDsorted)

#create a subset of only the first 500,000 rows of the sorted table
subset<-head(Red_LDsorted, 500000)

#export subset as csv to open in excel, note this may take a few seconds
write.csv(subset, "Red_LDsorted500K.csv", row.names = TRUE)
                 
#####MonomorphicLoci#####
#to get number of monomorphic loci, should import pops as separate files to do this
monomorphLoci<-gl.report.monomorphs(Fourpop.gl_3)

#to get a list of monomorphic loci, should import pops as separate files to do this
x <- Fourpop.gl
na.counter <- 0
loc.list <- array(NA, nLoc(x))
mat <- as.matrix(x)
lN <- locNames(x)
for (i in 1:nLoc(x)) {
  row <- mat[, i]
  if (all(row == 0, na.rm = TRUE) | all(row == 2, na.rm = TRUE) |
      all(is.na(row))) {
    loc.list[i] <- lN[i]
    if (all(is.na(row))) {
      na.counter = na.counter + 1
    }
  }
}
loc.list <- loc.list[!is.na(loc.list)]

# display list of monomorphic loci
print(loc.list)


#####FixedDifferences#####
#This script takes SNP data or sequence tag P/A data grouped into populations in a genlight object (DArTSeq) and generates 
#a matrix of fixed differences between populations taken pairwise
Fixdiff<-gl.fixed.diff(Fourpop.gl_3, tloc=0.05)
Fixdiff


#####Principal Coordinates Analysis#####
#Import structure file to genind object
data <- read.structure("run9_noMAF_K2.stru", n.ind=237, n.loc=51741,onerowperin=FALSE, col.lab=1, col.pop=2,row.marknames=1)

#convert the genind from above to a genlight
data2 <- gi2gl(data, parallel = FALSE)

#need to calculate a distance between individs - below is for 
#scaled Euclidean (if scale=FALSE then is regular Euclidean). Note unscaled Euclidean is severely affected by missing values adn should only be used for complete data 
D <- gl.dist.ind(data2, method="Euclidean", scale=TRUE, plot.out=TRUE)
#Simple Mismatch Distance
#D <- gl.dist.ind(data2, method="Simple", scale=FALSE, plot.out=TRUE)
#Absolute Mismatch Distance
#D <- gl.dist.ind(data2, method="Absolute", scale=FALSE, plot.out=TRUE)
#Czekanowski (Manhattan) Distance
#D <- gl.dist.ind(data2, method="Manhattan", scale=FALSE, plot.out=TRUE)

#to do the PCoA
pc <- gl.pcoa(D)

#plot the PCoA. for a 2D plot remove the zaxis part.
PCoA <- gl.pcoa.plot(pc, data2, scale=FALSE, ellipse = FALSE, pop.labels = 'pop',interactive=TRUE, zaxis = 3)
