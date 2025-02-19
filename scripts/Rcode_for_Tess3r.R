#to load tess3r need to use devtools
#If you don't have devtools do 
#install.packages("devtools")
#library(devtools)

#then do
#devtools::install_github("bcm-uga/TESS3_encho_sen")

library(tess3r)

#set your working directory to where you want it

#Need a tab-delimited file with 1st column being longitude, 2nd latitude, then 3rd column starts the genotypes.
#Each sample is 2 rows (like Structure format) but there should be no header columns or rows. 
#The following converts the .txt file to a dataframe, which is required for Tess3r.
snp6479<-read.delim("snps6479.txt", header = FALSE)

#This converts the info in the dataframe so it can be used by Tess3r. 
#May take a long few seconds to complete
data<-tess2tess3(dataframe=snp6479, TESS=TRUE, diploid = TRUE, FORMAT=2, extra.row = 0, extra.column = 0)

#This labels the genotype and coord data in your file so you can call it in subsequent code. 
genotype <-data$X
coordinates <-data$coord

##################### What to do if using Geographic Distance Matrix for W parameter #########################
#used same input file as here, opened in TESS 2.3 GUI, followed directions in Bk2 pg 15 to calculate Great Circle Distance matrix (did 
#not check the box for data stored in one line per individ though)

#opened matrix file in notepad++ and had to delete trailing white space (appeared at the end of each line) by doing:
#Edit -> Blank Operations -> Trim Trailing Space

#saved file as txt

#in R need to turn this text into a matrix

distance_matrix<-as.matrix(read.table("C:/Users/nicole.vollmer/Documents/NC_Tursiops/Analyses/Tess3r/6479snps/TestwGCdistancesAsW/GC_distances", sep = " "))

#and need to make sure to remove all column and row names - do this just to be sure otherwise will get error in Tess3 code
colnames(distance_matrix)<-NULL
rownames(distance_matrix)<-NULL

#now run the tess3 code
tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:10, method = "projected.ls", ploidy = 2, rep=100, W=distance_matrix)


#####################
#This tests K1-10 and calculates cross-validation score for each K. 
#There are 2 methods to choose from, projected.ls and qp: 
#the former is fast (6479 snps took seconds for 1 rep, and 1hr 40min for 100 reps) 
#and the latter is slower (6479 snps took 15min for 1 rep). 
#qp = alternating quadratic programming, 
#projected = alternating projected least squares. 
#Caye et al. 2017 shows that the projected algorithm is a good approximation to the solutions of 
#the qp algorithm.
tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:10, method = "projected.ls", ploidy = 2, rep=100)

#This plots the cross-validation score for each K. 
#There are 2 options for this CV-score that you can plot: 1) cross-entropy error  
#and 2) root mean square error (RMSE). The latter is the default. 
#The interpretation of this plot is similar to the cross-entropy plot of LEA or the cross-validation plot of
#ADMIXTURE. The cross-validation criterion is based on the prediction of a fraction of masked genotypes via 
#matrix completion, and comparison with masked values considered as the truth. Smaller values of the 
#cross-validation criterion mean better runs. The best choice for the K value is when the cross-validation 
#curve exhibits a plateau or starts increasing. But NOTE you may not see a plateau or a minimum value that 
#makes sense, but you should still keep going and look at the results of different K's to see if there are 
#any biologically meaningful populations. In the previous Tess3 step there are also reps you can do per K 
#with 'rep=' (default is 1). 
plot(tess3.obj, pch = 19, col = "blue", xlab = "Number of ancestral populations", ylab = "Cross-validation score")

#Makes a q matrix for what ever K you want.
q.matrix <- qmatrix(tess3.obj, K = 4)

#to export your qmatrix as a csv
write.csv(q.matrix, "K4_distances.csv", row.names = TRUE)

#Makes a Structure-like barplot for whatever K you chose for the q matrix. 
#There is a 'sort.by.Q' option (default is TRUE).  To see what order things are plotted/sorted can just type
#'bp'.  See the GitHub page about changing the color palette for the barplot.
barplot(q.matrix, border = NA, space = 0, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp

#To add labels to x-axis (matches order of 'bp')
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

#To plot the ancestry coefficients on a map you first need to download an appropriate ascii or raster or tiff
#of the geographic area you need. One option is https://maps.ngdc.noaa.gov/viewers/wcs-client/ 
#(Grid Extract Tool). Draw a box around the area you want, can keep layer as ETOPO1 and hit Download Data 
#button. Note this website used to let you export as an ascii file but as of March 2021 this was no longer 
#an option and you can only download as GeoTiff. After you download the GeoTiff you need to convert it to 
#ascii. Here is some R code to do that. Need to install/load Raster package. 
#[Another option is to use the 'tess_for_seas.r' code that NLV has in her NCTt Tess3r folder. 
#This option downloads a map from NOAA, but the map quality is much poor that the method below.]
library(raster)
f <- "path/to/downloaded/file.tif"
r <- raster(f)
ra <- aggregate(r, fact=2)  ## By default aggregates using mean (fun=mean)
writeRaster(ra, "path/to/outfile.asc", format="ascii")

#Now you have your ascii, but in order for the ancestry coefficients to plot on the water, and not the land,
#you need to manually change the signs off all the numbers from line 7 and beyond. So positive numbers
#become negative and vice versa. NOTE there may be a much better way to do this that I haven't found yet! 
#Open outfile.asc in NotePad++ and do the following:
  
#1)	Copy and paste the first 6 lines into a new file so you don't mess those up
#2)	Look at all the 1st number in each line and see if there are any that are negative and note that you will have to delete those negative signs individually at the very end
#3)	Find [space] and replace with [space][space][space]
#4)	Find [space][space][space] -  and replace with [space]+
#5)	Find [space][space][space] and replace with [space] -
#6)	Find + and replace with nothing
#7)	Manually insert a - before the 1st number in each line (except the ones that should be positive, see step 2)
#8)	Paste back in the 1st 6 lines and save

#This plots the ancestry values from the qmatrix on the water. For a crisper image (specifically the 
#coloring of ancestry coefficient colors) can increase resolution to 5000 for both numbers but this will 
#increase the processing time. To change the size of the sample dots change the cex #
plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude",  resolution = c(300,300), cex = 1, raster.filename="C:/Users/nicole.vollmer/Documents/NC_Tursiops/Analyses/Tess3r/wNAcoastalineforNCTT/outfile.asc")

#Can color the clusters how you want with the following code (need to figure out the correct order to match 
#colors to pops you want) 
my.colors <- c("blue", "red", "darkorange", "darkgreen", "darkolivegreen3", "purple", "yellow", "pink","grey")
#my.colors <-c("red","darkorange", "green")
my.palette <- CreatePalette(my.colors, 9)

#Then add col.palette=my.palette to the plot function
plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude",  resolution = c(300,300), cex = 1, raster.filename="C:/Users/nicole.vollmer/Documents/NC_Tursiops/Analyses/Tess3r/wNAcoastalineforNCTT/outfile.asc", col.palette = my.palette)

