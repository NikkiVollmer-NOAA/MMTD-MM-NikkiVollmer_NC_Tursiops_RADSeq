library(adegenet)

#set wd

#For Structure file
#data <- import2genind("unrelneut_6367loci.stru", n.ind=128, n.loc=6367,onerowperin=FALSE, col.lab=1, col.pop=2,row.marknames=1)

#For Genepop file
data <- import2genind("6479snps.gen",ncode=3L)

data


#NA.method can be asis, mean (replace with mean allele freq) or zero, note xval won't work if you leave NAs in
#data_scaledZ <- tab(data, NA.method="zero")
data_scaled <- tab(data,NA.method="mean")

grp<-find.clusters(data_scaled, max.n.clust=40)


#Make sure to type "groupMean" as you see it (don't capitalize the g)
xval<-xvalDapc(data_scaled, grp$grp, n.pca.max=250,training.set=0.9, result = "groupMean", center=TRUE,
               scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
xval

dapc1<-dapc(data_scaled, grp$grp)

scatter(dapc1)

#to change orientation of graph:
dapc1$ind.coord[,1] <- -dapc1$ind.coord[,1]
dapc1$loadings[,1] <- -dapc1$loadings[,1]
scatter(dapc1)


#to get more info type: ?scatter.dapc
#cex= size of marker
#solid = making marker transparent
#clab = 0 makes cluster label no appear 
#cstar=0 removes lines connecting center of cluster to markers
#pch changes the type of marker (20 = circle)
#cleg is size of legend
#ratio is size of insets
#to turn on individ labels use label.inds = list(air=0.5), the air = is how far away from dot to print label,
#if the distance is too large it will prevent close by labels from not printing

#col<-c("red", "darkorange", "green")
col=c("darkgreen", "darkorange", "green","red")
#scatter(dapc1, cex = 3, solid = .4, clab=0, cstar=0, pch = 20, posi.da="topright", ratio.da = 0.2, 
        #scree.pca = TRUE, ratio.pca= .20, posi.pca ="bottomleft", col=col,leg = TRUE, cleg = 1.5,
        #posi.leg ="bottomright", txt.leg=paste(c("Red", "Orange", "Green")))

scatter(dapc1, cex = 3, solid = .6, clab=0, cstar=0, pch = c(17, 15, 18, 16), posi.da="topright", ratio.da = 0.2, 
scree.pca = TRUE, ratio.pca= .20, posi.pca ="bottomright", col=col, grid = TRUE)

#to see membership probabilities
#dapc1$posterior

#to export membership probabilities into a csv file in home directory
write.csv(dapc1$posterior, "MembershipProbabilitiestest2_thia.csv", row.names=TRUE)

#to get summary info like grp size
summary(dapc1)


