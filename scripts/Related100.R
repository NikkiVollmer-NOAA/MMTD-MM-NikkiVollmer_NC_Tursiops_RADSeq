library(related)

input<-readgenotypedata("/home/nicole.vollmer@nmfs.local/bluepop.txt")

output<-coancestry(input$gdata, dyadml=1, lynchli=1, lynchrd=1, quellergt=1, ritland=1, trioml=1, wang=1)

write.table(output$relatedness, "~/Desktop/related/BluePop/relatedness.txt", sep="\t")

simdata<-familysim(input$freqs,100)

output2<-coancestry(simdata, dyadml=1, trioml = 1, quellergt=1, wang=1, lynchli=1, lynchrd=1)

simrel<-cleanuprvals (output2$relatedness, 100)

write.table(simrel, "CleanUp_simR.txt", sep="\t")

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


trioml <- rep ("tri", 100)
wang <- rep ("W", 100)
lynchli <- rep("lyli", 100)
lynchrd <- rep("lyrd", 100)
quellergt <- rep ("QG", 100)
dyadml <- rep ("di", 100)
estimator2 <- c( trioml , wang , lynchli, lynchrd, quellergt , dyadml )
Estimator <- rep (estimator2 , 4)


po <- rep ("Parent - Offspring ", (6 * 100) )
fs <- rep ("Full - Sibs ", (6 * 100) )
hs <- rep ("Half - Sibs ", (6 * 100) )
ur <- rep (" Unrelated ", (6 * 100) )
relationship <- c(po , fs , hs , ur )


relatednesspo <- c( triomlpo , lynchlipo, lynchrdpo, wangpo , quellergtpo , dyadmlpo ) 
relatednessfs <- c( triomlfs , lynchlifs, lynchrdfs, wangfs , quellergtfs , dyadmlfs ) 
relatednesshs <- c( triomlhs , lynchlihs, lynchrdhs, wanghs , quellergths , dyadmlhs )
relatednessur <- c( triomlur , lynchliur, lynchrdur, wangur , quellergtur , dyadmlur )
Relatedness_Value <- c( relatednesspo , relatednessfs , relatednesshs , relatednessur )

combineddata <- as.data.frame(cbind (Estimator,relationship,Relatedness_Value))
combineddata$Relatedness_Value <- as.numeric(as.character(combineddata$Relatedness_Value))


ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 , 1.0) ) +  geom_boxplot () + facet_wrap (~ relationship )
dev.off()

urval <- rep (0 , 100)
hsval <- rep (0.25 , 100)
fsval <- rep (0.5 , 100)
poval <- rep (0.5 , 100)
relvals <- c( poval , fsval , hsval , urval)


cor ( relvals , simrel [ , 5])
cor ( relvals , simrel [ , 6])
cor ( relvals , simrel [ , 7])
cor ( relvals , simrel [ , 8])
cor ( relvals , simrel [ , 10])
cor ( relvals , simrel [ , 11])


gg_bp <- ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 , 1.0) ) +  geom_boxplot () + facet_wrap (~ relationship )

data <- apply((layer_data(gg_bp)), 2, as.character)

write.csv(data, "boxplot_data.csv")




