library(related)

input<-readgenotypedata("/home/nicole.vollmer@nmfs.local/related/orange.txt")

output<-coancestry(input$gdata, dyadml=1, lynchli=1, lynchrd=1, quellergt=1, ritland=1, trioml=1, wang=1)

write.table(output$relatedness, "./OrangePoprelatedness.txt", sep="\t")

simdata<-familysim(input$freqs,500)

output2<-coancestry(simdata, dyadml=1, trioml = 1, quellergt=1, wang=1, lynchli=1, lynchrd=1)

simrel<-cleanuprvals (output2$relatedness, 500)

write.table(simrel, "CleanUp_OrangePop_simR.txt", sep="\t")

triomlpo <- simrel [1:500 , 5] 
triomlfs <- simrel [(500 + 1) : (2 * 500) , 5] 
triomlhs <- simrel [((2 * 500) + 1) : (3 * 500) , 5] 
triomlur <- simrel [((3 * 500) + 1) : (4 * 500) , 5] 
wangpo <- simrel [1:500 , 6]
wangfs <- simrel [(500 + 1) : (2 * 500) , 6] 
wanghs <- simrel [((2 * 500) + 1) : (3 * 500) , 6]
wangur <- simrel [((3 * 500) + 1) : (4 * 500) , 6]
lynchlipo <- simrel [1:500 , 7]
lynchlifs <- simrel [(500 + 1) : (2 * 500) , 7] 
lynchlihs <- simrel [((2 * 500) + 1) : (3 * 500) , 7]
lynchliur <- simrel [((3 * 500) + 1) : (4 * 500) , 7] 
lynchrdpo <- simrel [1:500 , 8]
lynchrdfs <- simrel [(500 + 1) : (2 * 500) , 8] 
lynchrdhs <- simrel [((2 * 500) + 1) : (3 * 500) , 8]
lynchrdur <- simrel [((3 * 500) + 1) : (4 * 500) , 8]  
quellergtpo <- simrel [1:500 , 10] 
quellergtfs <- simrel [(500 + 1) : (2 * 500) , 10] 
quellergths <- simrel [((2 * 500) + 1) : (3 * 500) , 10] 
quellergtur <- simrel [((3 * 500) + 1) : (4 * 500) , 10] 
dyadmlpo <- simrel [1:500 , 11] 
dyadmlfs <- simrel [(500 + 1) : (2 * 500) , 11] 
dyadmlhs <- simrel [((2 * 500) + 1) : (3 * 500) , 11]
dyadmlur <- simrel [((3 * 500) + 1) : (4 * 500) , 11]


trioml <- rep ("tri", 500)
wang <- rep ("W", 500)
lynchli <- rep("lyli", 500)
lynchrd <- rep("lyrd", 500)
quellergt <- rep ("QG", 500)
dyadml <- rep ("di", 500)
estimator2 <- c( trioml , wang , lynchli, lynchrd, quellergt , dyadml )
Estimator <- rep (estimator2 , 4)


po <- rep ("Parent - Offspring ", (6 * 500) )
fs <- rep ("Full - Sibs ", (6 * 500) )
hs <- rep ("Half - Sibs ", (6 * 500) )
ur <- rep (" Unrelated ", (6 * 500) )
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

urval <- rep (0 , 500)
hsval <- rep (0.25 , 500)
fsval <- rep (0.5 , 500)
poval <- rep (0.5 , 500)
relvals <- c( poval , fsval , hsval , urval)


cor ( relvals , simrel [ , 5])
cor ( relvals , simrel [ , 6])
cor ( relvals , simrel [ , 7])
cor ( relvals , simrel [ , 8])
cor ( relvals , simrel [ , 10])
cor ( relvals , simrel [ , 11])


gg_bp <- ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 , 1.0) ) +  geom_boxplot () + facet_wrap (~ relationship )

data <- apply((layer_data(gg_bp)), 2, as.character)

write.csv(data, "OrangePop_boxplot_data.csv")




