##How to run the different estimators on your data (see below for doing simulations)

#loads the package and required dependencies
library(related)
Loading required package: tools
Loading required package: ggplot2
Use suppressPackageStartupMessages() to eliminate package startup messages.

#loads the data (formatted according to Tutorial), make sure you have navigated to the correct place and file
input<-readgenotypedata("~/Desktop/related/XXXPop/XXX.txt")

#runs all 7 relatedness indices, but not doing 95% CI for any (would make “=2” if you wanted to do that). 
#NOTE for a pop with 543 samples, this took ~1hr 15min to run on my MMMGL desktop with 32GB RAM. A pop with 42 samples took 19 seconds.
output<-coancestry(input$gdata, dyadml=1, lynchli=1, lynchrd=1, quellergt=1, ritland=1, trioml=1, wang=1)

#exports each dataframe (if “=1”) into a .txt file. Can cut/paste the whole blue chunk into R. 
#Make sure to change the destination name for each line. 
#There is code here for 19 loci for exporting freq data, would need to add/subtract lines of code if using a different number of loci. 
#If doing 95% CI (“=2”) need to also do: 
#write.table(output$relatedness.ci95, "~/Desktop/related/BluePop/delta7.txt", sep="\t")
#write.table(output$delta7.ci95, "~/Desktop/related/BluePop/delta7.txt", sep="\t")
#write.table(output$delta8.ci95, "~/Desktop/related/BluePop/delta7.txt", sep="\t")
#write.table(output$inbreeding.ci95, "~/Desktop/related/BluePop/inbreeding.txt", sep="\t")

write.table(output$relatedness, "~/Desktop/related/BluePop/relatedness.txt", sep="\t")
write.table(output$delta7, "~/Desktop/related/BluePop/delta7.txt", sep="\t")
write.table(output$delta8, "~/Desktop/related/BluePop/delta8.txt", sep="\t")
write.table(output$inbreeding, "~/Desktop/related/BluePop/inbreeding.txt", sep="\t")
write.table(output$freqs$locus1, "~/Desktop/related/BluePop/locus1freq.txt", sep="\t")
write.table(output$freqs$locus2, "~/Desktop/related/BluePop/locus2freq.txt", sep="\t")
write.table(output$freqs$locus3, "~/Desktop/related/BluePop/locus3freq.txt", sep="\t")
write.table(output$freqs$locus4, "~/Desktop/related/BluePop/locus4freq.txt", sep="\t")
write.table(output$freqs$locus5, "~/Desktop/related/BluePop/locus5freq.txt", sep="\t")
write.table(output$freqs$locus6, "~/Desktop/related/BluePop/locus6freq.txt", sep="\t")
write.table(output$freqs$locus7, "~/Desktop/related/BluePop/locus7freq.txt", sep="\t")
write.table(output$freqs$locus8, "~/Desktop/related/BluePop/locus8freq.txt", sep="\t")
write.table(output$freqs$locus9, "~/Desktop/related/BluePop/locus9freq.txt", sep="\t")
write.table(output$freqs$locus10, "~/Desktop/related/BluePop/locus10freq.txt", sep="\t")
write.table(output$freqs$locus11, "~/Desktop/related/BluePop/locus11freq.txt", sep="\t")
write.table(output$freqs$locus12, "~/Desktop/related/BluePop/locus12freq.txt", sep="\t")
write.table(output$freqs$locus13, "~/Desktop/related/BluePop/locus13freq.txt", sep="\t")
write.table(output$freqs$locus14, "~/Desktop/related/BluePop/locus14freq.txt", sep="\t")
write.table(output$freqs$locus15, "~/Desktop/related/BluePop/locus15freq.txt", sep="\t")
write.table(output$freqs$locus16, "~/Desktop/related/BluePop/locus16freq.txt", sep="\t")
write.table(output$freqs$locus17, "~/Desktop/related/BluePop/locus17freq.txt", sep="\t")
write.table(output$freqs$locus18, "~/Desktop/related/BluePop/locus18freq.txt", sep="\t")
write.table(output$freqs$locus19, "~/Desktop/related/BluePop/locus18freq.txt", sep="\t")
write.table(output$freqs$locus18, "~/Desktop/related/BluePop/locus18freq.txt", sep="\t")
write.table(output$freqs$locus19, "~/Desktop/related/BluePop/locus19freq.txt", sep="\t")

#How to run simulations based on your data and determine which estimator is the best to use

#If running this on the SERVER the usual COPY/PASTE shortcuts DO NOT work: instead to COPY use CONTROL+INSERT and to PASTE use SHIFT+INSERT

#first input your data
input<-readgenotypedata(“C:/Users/nicole.vollmer/Desktop/related/RedPop/Red.txt”)

#the “compareestimators” command uses your genotype data (“input”) and simulates pairs of individs 
#(below set as 100 pairs) for each type of relationship. The function generates simulated genotypes of known relatedness
#(PO, FS, HS, unrelated), calculated relatedness using 4 moment-based estimators
#(lynchli, lynchrd, quellergt, wang), and then plots the data to view.  
#Also it estimates and prints out Pearson correlation coefficients between observed and expected relatedness values so you can compare
#the estimators and see which correlates best with expected values (you want to use the estimator that produces the highest number).
compareestimators(input,100)

#to do the above with some of the other estimators first load the data
input<-readgenotypedata(“C:/Users/nicole.vollmer/Desktop/related/RedPop/Red.txt”)

#generate simulated individs based on your data (here I simulated 100 of each type: 100 PO, 100 FS, 100 HS, 100 related)
simdata<-familysim(input$freqs,100)

#now estimate relatedness on the simulated data using the desired estimators, 
#here I am testing the 2 likelihood methods (dyadml, trioml) and 4 moment methods (queller&GN, wang, lynchli, lynchrd). 
#NOTE this takes a long time (130 individs at 19 loci on PC desktop with 8GB RAM took ~2.5 hrs, 265 individs at 19 loci on PC with 
#8GB RAM took ~5 hrs., 465 individs took 21hrs.)
output<-coancestry(simdata, dyadml=1, trioml = 1, quellergt=1, wang=1, lynchli=1, lynchrd=1)

# the resulting file with contain ALL pairwise estimates of relatedness, not just the ones we are interested in. 
#To reduce to only the desired values use cleanuprvals command. NOTE if you simulated 100 pairs with simdata command, 
#then you need to use 100 here too. If you simulated a different number with simdata, then need to make it the same here.
simrel<-cleanuprvals (output$relatedness, 100)

# to get ALL relatedness values calculated from simulation. You will have to move over the top row one cell to the right (so cell A1 is blank, 
#and cell L1 says dyadml). Note that this file may be too large (i.e., have too many rows) to open in excel. 
#NLV has python code to work with such a file (VARCalcRelated.ipynb saved on NLV personal laptop home directory)
write.table(output$relatedness, "All_simR.txt", sep="\t")

#OR for files with many rows that you will analyze in python do csv

write.csv(output$relatedness, "All_simR.csv")

# to get only the relatedness values calculated after the clean up. You will have to move over the top row one cell to the right
#(so cell A1 is blank, and cell L1 says dyadml).
write.table(simrel, "CleanUp_simR.txt", sep="\t")

# next the data needs to be parsed based on relatedness and by estimator used. 
#Since 100 individs were simulated then 100 of each type were created: 100 PO, 100 FS, 100 HS, 100 related. 
#And each estimator will be in a different, but estimator-specific, column (see related manual).
#Can cut and past the following code to parse results. NOTE it is specific for the 6 estimators estimated above and would need to be edited, 
#both name and column info, if different estimators are estimated. 
triomlpo <- simrel [1:100 , 5] 
triomlfs <- simrel [(100 + 1) : (2 * 100) , 5] 
triomlhs <- simrel [((2 * 100) + 1) : (3 * 100) , 5] 
triomlur <- simrel [((3 * 100) + 1) : (4 * 100) , 5] 
wangpo <- simrel [1:100 , 6]
wangfs <- simrel [(100 + 1) : (2 * 100) , 6] 
wanghs <- simrel [((2 * 100) + 1) : (3 * 100) , 6]
wangur <- simrel [((3 * 100) + 1) : (4 * 100) , 6]
lynchlipo <- simrel [1:100 , 7]
lynchlifs <- simrel [(100 + 1) : (2 * 100) , 7] 
lynchlihs <- simrel [((2 * 100) + 1) : (3 * 100) , 7]
lynchliur <- simrel [((3 * 100) + 1) : (4 * 100) , 7] 
lynchrdpo <- simrel [1:100 , 8]
lynchrdfs <- simrel [(100 + 1) : (2 * 100) , 8] 
lynchrdhs <- simrel [((2 * 100) + 1) : (3 * 100) , 8]
lynchrdur <- simrel [((3 * 100) + 1) : (4 * 100) , 8]  
quellergtpo <- simrel [1:100 , 10] 
quellergtfs <- simrel [(100 + 1) : (2 * 100) , 10] 
quellergths <- simrel [((2 * 100) + 1) : (3 * 100) , 10] 
quellergtur <- simrel [((3 * 100) + 1) : (4 * 100) , 10] 
dyadmlpo <- simrel [1:100 , 11] 
dyadmlfs <- simrel [(100 + 1) : (2 * 100) , 11] 
dyadmlhs <- simrel [((2 * 100) + 1) : (3 * 100) , 11]
dyadmlur <- simrel [((3 * 100) + 1) : (4 * 100) , 11]

# the following creates a list of labels for the different estimators, each repeated the appropriate number of times (i.e., 100)
trioml <- rep ("tri", 100)
wang <- rep ("W", 100)
lynchli <- rep("lyli", 100)
lynchrd <- rep("lyrd", 100)
quellergt <- rep ("QG", 100)
dyadml <- rep ("di", 100)
estimator2 <- c( trioml , wang , lynchli, lynchrd, quellergt , dyadml )
Estimator <- rep (estimator2 , 4)

#creates a list of labels for the different relationships
po <- rep ("Parent - Offspring ", (6 * 100) )
fs <- rep ("Full - Sibs ", (6 * 100) )
hs <- rep ("Half - Sibs ", (6 * 100) )
ur <- rep (" Unrelated ", (6 * 100) )
relationship <- c(po , fs , hs , ur )

# combines the different values for each estimator based on relatedness type, as lists
relatednesspo <- c( triomlpo , lynchlipo, lynchrdpo, wangpo , quellergtpo , dyadmlpo ) 
relatednessfs <- c( triomlfs , lynchlifs, lynchrdfs, wangfs , quellergtfs , dyadmlfs ) 
relatednesshs <- c( triomlhs , lynchlihs, lynchrdhs, wanghs , quellergths , dyadmlhs )
relatednessur <- c( triomlur , lynchliur, lynchrdur, wangur , quellergtur , dyadmlur )
Relatedness_Value <- c( relatednesspo , relatednessfs , relatednesshs , relatednessur )

#combines the data
combineddata <- as.data.frame(cbind (Estimator,relationship,Relatedness_Value))
combineddata$Relatedness_Value <-  
as.numeric(as.character(combineddata$Relatedness_Value))

#plots the data. NOTE the ylim might need to be changed depending on the data
ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 , 1.0) ) +  geom_boxplot () + facet_wrap (~ relationship )

#If running R on the SERVER type this in after running the ggplot (or ANY plot) 
#function (this isn’t necessary if just running on your desktop RStudio or R)
dev.off()

#need to manually calculate Pearson correlation coefficient this time, 
#need to first create vectors containing the appropriate relatedness values the appropriate number of times (here 100)
urval <- rep (0 , 100)
hsval <- rep (0.25 , 100)
fsval <- rep (0.5 , 100)
poval <- rep (0.5 , 100)
relvals <- c( poval , fsval , hsval , urval)

#now need to code for comparing observed and expected values. 
#To do this need to call the appropriate column from the simrel data dataframe for each estimator (see related manual table on pg 8)
cor ( relvals , simrel [ , 5])
cor ( relvals , simrel [ , 6])
cor ( relvals , simrel [ , 7])
cor ( relvals , simrel [ , 8])
cor ( relvals , simrel [ , 10])
cor ( relvals , simrel [ , 11])



#to get data used to make boxplots (ymin, lower, middle, upper, ymax, outliers, etc.) use the code below. 
#Output will be labeled 1-24 corresponding to the first box on the left in the top plot on the left as #1, and the last plot on the right
#in the bottom plot on the right as #24. Type “?geom_boxplot” to get more info on what each variable is. 
#The “boxplot_data.csv” will be saved in your working directory.

gg_bp <- ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 , 1.0) ) +  geom_boxplot () + facet_wrap (~ relationship )

data <- apply((layer_data(gg_bp)), 2, as.character)

write.csv(data, "boxplot_data.csv")
