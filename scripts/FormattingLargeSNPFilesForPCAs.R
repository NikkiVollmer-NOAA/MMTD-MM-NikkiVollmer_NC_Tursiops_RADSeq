library(dplyr)
library(tidyverse)
library(reshape2)
#read in the txt file after initially editing in NOTEPAD++ (see NLV Bk 11 pg 1). Note if the file is too large to convert from vcf to stru in the PGDSpider
#GUI then just use PGDSpider using the command line and all 'should' work fine
data<-read.table("as1pop.txt", header = TRUE, sep = "\t")
data

#remove the popdata column that pgd spider inserts = column X.1
data_clean<-subset(data, select = -c(X.1))
data_clean

#read in a tab delimited txt file that includes a column of ID and Location. Make sure this file is in the
#same order as the IDs in the data file so they match up correctly
locations<-read.delim("K3basin_51741snps.txt", header = TRUE, sep = "\t")
locations

#add the "Location" column from the locations dataframe to the end of the data_clean dataframe
data_all<-cbind(data_clean, locations[c("K3_basin")])

#move the Location column to be right after the column with the ID names ("X") [relocate function is part of dplyr]
data_all_new<-relocate(data_all, "K3_basin", .after= "X")

#export this new table as an .stru file to use in PCA code. Will need to open again in NOTEPAD++ to remove
#quote marks, replace the headers for columns 1 and 2 with tabs [so should be tabtabfirstlocusname], and remove line breaks
write.table(data_all_new,"test.stru",sep="\t", row.names = FALSE)

#the code below takes the data_clean file created above (structure file without popdata/info column) and counts up the 
#number of each allele and -9's(missing data) for each individual which can be used to calculate % missing data per individual
count<-data_clean %>%
  gather(var, val, -X) %>% #Transforming the data from wide to long format
  group_by(val, X) %>% #Grouping 
  summarise(count = n()) %>% #Performing the count
  dcast(X~val, value.var = "count") #Reshaping the data

#this creates a .txt file with the above info
write.table(count,"count.txt",sep="\t", row.names = FALSE)
