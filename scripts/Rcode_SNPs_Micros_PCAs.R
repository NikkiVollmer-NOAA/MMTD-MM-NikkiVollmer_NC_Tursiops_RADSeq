library(adegenet)
library(factoextra)
library(plotly)
library(ggplot2)
#library(adegraphics)This is another package that can be used to make nice plots for PCAs, but some functions, like s.class, are already
#in adegenet, and if you load both packages the s.class in adegraphics gets masked and I couldn't figure out how to get it to work even when
#trying adegraphics::s.class

############Loading Data into Genind Object###########
#this can work for smaller micro datasets and/or SNP datasets that are a few thousand loci.
#consider using a Genlight object for very large SNP datasets

#First need to import data file
#For input file (.stru file - this needs to be .stru and not .txt), make sure to put in correct number of individs and loci, 
#and make sure there are ONLY column labels for the loci, do not have headers for samples or pop column(s)
#double check if data is in one row per individ or not and change code accordingly
#col.lab is the column where your sample labels are; col.pop is the the column where you have your 
#identifier (sex, extraction type, structure pop, etc)
#row.marknames is an integer giving the index of the row containing the names of the markers
data <- read.structure("6479snps_run4_K4structpops.stru", n.ind=141, n.loc=6479,onerowperin=FALSE, col.lab=1, col.pop=2,row.marknames=1)

#need to scale data and can replace missing (NA) values, if set to 'mean' it uses the mean allele frequency. Can also choose 'asis' to 
#keep NAs, or 'zero' to replace all NA with 0
#Note if you get this message when running 
#"Warning message: In .local(x, ...) : Some scaling values are null. Corresponding alleles are removed." 
#it is because there are non-polymorphic loci that are being excluded/removed.
data.scaled <- scaleGen(data,NA.method="mean")


###########Running a PCA with dudi.pca#############
#run the PCA on the scaled data, if you want to see the scree plot set scannf = TRUE, 
#if scannf set to FALSE set nf to number of axes to keep, below I am keeping 3
pca1 <- dudi.pca(data.scaled,center=FALSE,scale=FALSE,scannf=FALSE,nf=4)


###########Scree Plots##########
#to look at scree plot, this is nice cause there is a lot of detail
fviz_eig(pca1)

#Below is another scree plot option from adegenet
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))


##########Visualizing PCA with fviz##########
#to look at basic PCA with each point labeled according to input data column 1
#You also can choose which axes to plot
#this also prints the estimated contribution of the PC's, which is something I haven't figured out how to do with s.class (see below)
#to show points and labels make geom = c("point", "text") [this is actually the default]
fviz_pca_ind(pca1, axes = c(1,2), geom = c("point", "text"), repel=TRUE)

#to look at color coded PCA with each point labeled according to input data column 1. You can color code the data based on 
#"cos2" which is the quality of representation of rows/columns from the results of the PCA
#"contrib" which is the contribution of rows/columns from the results of the PCA
# and see ?fviz_pca_ind for other options
#repel = TRUE makes it so labels don't overlap but this doesn't look good if you have a tight cluster of points
fviz_pca_ind(pca1, col.ind="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = FALSE)

##########Visualizing PCA with s.class##########
#another way for plotting PCA results that provides options not available (I think) with fviz
#s.class has lots of options to play around with, one negative (compared to fviz) is that I can't figure out how to label x/y axes
#with %PCA explained by each
#xlim=c(-70,100)

s.class(pca1$li, xax =1, yax = 2, fac=pop(data), pch = c(15, 16, 17, 18, 19) [pop(data)], col=transp(c("darkorange","red","darkgreen","darkolivegreen3","black"), .6), axesell = FALSE, cstar=0, cpoint=3, label = NULL, addaxes = TRUE, xlim=c(-30,70))

#adds the pcas/eigenvalues plot as a subplot in the pca
#I am plotting the first 20 pcs, coloring in the first 4 with the first two as black (represented axes) and the
#the third gray because I retained a total of 3 axes
add.scatter.eig(pos = "topright",pca1$eig[1:20],4,3,4,ratio=.2)

#legend("bottomleft", pch=19,col=c("yellow", "red", "blue2"), c("Unclustered", "Coastals", "Offshores"), bty="o", box.col="black", cex=1.2)
#Tried adding legend function after add.scatter.eig but it doesn't plot correctly

#######Color-coding Basic Scatter Plot Based on PCA Results########
#there aren't many options with the 'plot' function (like I don't think you can draw ellipses), 
#but this works for a simple color-coded scatter plot of PCA. You would want something a bit more fancy like the above stuff
#for any publication
#to use this function load data and run PCA (dudi.pca) as above
#then assign colors to whatever ever categories you are interested in, for example:
pop.color<-data@pop
pop.color = as.character(pop.color)
pop.color[pop.color == '1'] = 'red'
pop.color[pop.color == '2'] = 'lightgreen'
pop.color[pop.color == '3'] = 'darkgreen'
pop.color[pop.color == '4'] = 'orange'
pop.color[pop.color == '0'] = 'yellow'
col=transp(pop.color)

#below makes your plot
plot(pca1$li[,1],pca1$li[,2],xlab='gPC1',ylab='gPC2',bg=col,col='black',
     pch=21,cex=2,xlim=c(-50,90),ylim=c(-35,40), main = "PCA colored based on K4")

#this is to add a legend for the plot
legend("bottomright", pch=19,col=c("red", "lightgreen", "darkgreen", "orange", "yellow"), 
       c("Red_21", "LightGreen_14", "DarkGreen_17", "Orange_84", "Unassigned_5"), bty="o", box.col="black", cex=1.2)


#########for interactive and 3D plots#############

#colors need to be in hex number format starting with #
# #008000 - green; #FF0000 - red; #0000FF - blue; #808080 - grey; #FFFF00 - yellow, #800080 - purple

p<-plot_ly(pca1$li, x= ~Axis1, y= ~Axis2, z = ~Axis3, color= I("black"), mode = 'markers', marker=list(size=10))%>%
  add_markers(marker=list(line=list(color='black',width=2,size = 2)), showlegend=T,
              symbol = ~pop(data), symbols = c('circle','circle', 'circle', 'circle', 'x'), color = ~pop(data),
              colors = c("#FF8000", "#FF0000","#006600", "#CCFFCC", "#FFFF00")) %>%
  layout(scene=list(xaxis=list(title="PC1"), yaxis=list(title="PC2"), 
                    zaxis=list(title="PC3")), legend=list(font = list(size = 10))) 

p


