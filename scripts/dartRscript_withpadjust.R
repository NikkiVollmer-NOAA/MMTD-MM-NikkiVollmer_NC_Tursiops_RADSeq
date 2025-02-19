library(dartR)
library(dplyr)
library(stats)

#may also need to run the following, see "Welcome to dartR" notes once you load dartR 
gl.install.vanilla.dartR()

#Import genepop file to genind object
Fourpop <- read.genepop(file = "6479snps_orangepop.gen", ncode=3L)

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


#####LD with test for multile corrections#####
#Import genepop file to genind object - will need to import each pop separately and run this code for each
orangepop <- read.genepop(file = "6479snps_orangepop.gen", ncode=3L)

#genind to genlight
orangepop.gl <- gi2gl(orangepop, verbose=5)
orangepop.gl_2 <- gl.recalc.metrics(orangepop.gl)
orangepop.gl_2$other$loc.metrics <- as.data.frame(orangepop.gl_2$other$loc.metrics)
orangepop.gl_3 <- gl.compliance.check(orangepop.gl_2) #make sure you get no errors with this line

#Calculates pairwise population based Linkage Disequilibirum across all loci using the specified number of cores
LD<-gl.report.ld(orangepop.gl_3, save = TRUE, outpath = "./", nchunks = 1, ncores = 4, probar=TRUE, verbose=5)

#above does not carry the locus names from the input file to the output LD. So below is how you attach those names
LD$loc1name <- locNames(orangepop.gl_3)[LD$loc1]
LD$loc2name <- locNames(orangepop.gl_3)[LD$loc2]


#in order to do the BY correction on the LD results, need to use another program. Here am using p.adjust in 'stats'.
#this function takes a vector of p-values as the input and will output all the adjusted p-values for each.
#Then you can look at the adjusted values, determine what FDR cutoff you want to use, and see which comparisons 
#are or are not significant. NOTE p.adjust also can do BH, holm, hochberg, hommel and bonferroni.

#DONT HAVE TO DO THIS: first extract the locus columns (1 and 2) and the p-value column (9) from the dataframe using column names
#datasubset<-LD[,c("loc1name", "loc2name", "p")]

#run p.adjust on the p-values, this creates a column of adjusted p-values
BYvalues<-p.adjust(LD$p, method="BY")

#add the BY values back onto the LD dataframe (e.g. adding a new column with the BY adjusted p-values)
LD_new<-cbind(LD, BYvalues)

#want to know how many are less than my cutoff, and what those pairs are, so first apply the cutoff
is_lessthan_0.05<-LD_new$BYvalues < 0.05

#then isolate only those belwo the cutoff
result<-LD_new[is_lessthan_0.05,]

#then save results, only those less than cut off, into a csv file
write.csv(result, "resultsnew.csv", row.names = TRUE)




#####below is code I used before I was applying any correction, and is probably not needed anymore#####

#the above LD function will produce a dataframe with potentially millions of rows (too many for excel to read)
#so need to sort dataframe by p column (may need to load LDallp.rdata file first)
orange_LDsorted<-LDallp[order(LDallp$p),]


#can check first 6 rows
head(orange_LDsorted)

#create a subset of only the first 500,000 rows of the sorted table
subset<-head(orange_LDsorted, 500000)

#export subset as csv to open in excel, note this may take a few seconds
write.csv(subset, "Orange_new_LDsorted500K.csv", row.names = TRUE)
                 



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

