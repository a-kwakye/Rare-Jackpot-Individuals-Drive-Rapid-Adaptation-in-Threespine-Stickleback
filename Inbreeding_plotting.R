library(ggplot2)
library(tidyr)
library(readr)
library(vioplot)

Inbreed_SC_13 <- read_table("SC13_SC13AF.res1", col_names = FALSE)
Inbreed_SC_14 <- read_table("SC14_SC13AF.res1", col_names = FALSE)
Inbreed_SC_15 <- read_table("SC15_SC13AF.res1", col_names = FALSE)
Inbreed_SC_17 <- read_table("SC17_SC13AF.res1", col_names = FALSE)





SC_13_vec=c(Inbreed_SC_13$X2)
SC_14_vec=c(Inbreed_SC_14$X2)
SC_15_vec=c(Inbreed_SC_15$X2)
SC_17_vec=c(Inbreed_SC_17$X2)

max_length <- max( length(SC_13_vec), length(SC_14_vec),
                  length(SC_15_vec),length(SC_17_vec))


SC_13_vec_2=c(SC_13_vec , rep(NA, max_length - length(SC_13_vec)))
SC_14_vec_2=c(SC_14_vec , rep(NA, max_length - length(SC_14_vec)))
SC_15_vec_2=c(SC_15_vec , rep(NA, max_length - length(SC_15_vec)))
SC_17_vec_2=c(SC_17_vec , rep(NA, max_length - length(SC_17_vec)))




Inbreed_scout_df=data.frame( SC_13_vec_2, SC_14_vec_2,SC_15_vec_2, 
                               SC_17_vec_2)
Inbreed_scout_df <-Inbreed_scout_df[-1, ]

labels =  c( "SC2013", "SC2014", "SC2015", "SC2017")
colnames(Inbreed_scout_df)=labels 

All_data=list( SC_13_vec, SC_14_vec, SC_15_vec , SC_17_vec )


cols=c("#9F2305", "#edea18", "#819d4e", "#155084")




pdf("Scout_inbreeding.pdf", width = 12.4, height = 8)
par(mar = c(5, 4, 4, 2))
vioplot(All_data, names = c( "SC2013", "SC2014", "SC2015", "SC2017"), 
        col = cols, 
        xlab ="Year", 
        ylab="Inbreeding coefficient (F)", 
        main = "Estimated inbreeding coefficient for each Scout timepoint", 
        cex.lab=34)

dev.off()


############# ngsrelate results ########

Inbreed_SC_13_ngs <- read_table("scout_13_inbreeding", col_names = FALSE)
Inbreed_SC_14_ngs <- read_table("scout_14_inbreeding", col_names = FALSE)
Inbreed_SC_15_ngs <- read_table("scout_15_inbreeding", col_names = FALSE)
Inbreed_SC_17_ngs <- read_table("scout_17_inbreeding", col_names = FALSE)
Inbreed_SC_20_ngs <- read_table("/gpfs/home/akwakye/inbreeding/Scout_2020_using_af_calculated_from_gatk", col_names = FALSE)
Inbreed_Lob_21_ngs <- read_table("/gpfs/home/akwakye/inbreeding/Loberg_21_using_af_calculated_from_gatk", col_names = FALSE)


SC_13_vec_ngs=c(Inbreed_SC_13_ngs$X3)
SC_14_vec_ngs=c(Inbreed_SC_14_ngs$X3)
SC_15_vec_ngs=c(Inbreed_SC_15_ngs$X3)
SC_17_vec_ngs=c(Inbreed_SC_17_ngs$X3)
SC_20_vec_ngs=c(Inbreed_SC_20_ngs$X3)
Lob_21_vec_ngs=c(Inbreed_Lob_21_ngs$X3)

All_data_ngs=list( SC_13_vec_ngs, SC_14_vec_ngs, SC_15_vec_ngs , SC_17_vec_ngs, SC_20_vec_ngs, Lob_21_vec_ngs )


cols=c("#9F2305", "#edea18", "#819d4e", "#155084", '#A4568B', '#8896AB')




pdf("Scout_inbreeding_ngs_relate_gatk_af.pdf", width = 12.4, height = 8)
par(mar = c(5, 4, 4, 2))
vioplot(All_data_ngs, names = c( "SC2013", "SC2014", "SC2015", "SC2017", 'SC2020','LOB21'), 
        col = cols, 
        xlab ="Year", 
        ylab="Inbreeding coefficient (F)", 
        main = "Estimated inbreeding coefficient for each Scout timepoint", 
        cex.lab=34)

dev.off()


############# ngsrelate outbred results ########

Inbreed_SC_13_ngs <- read_table("scout_13_inbreeding_outbred", col_names = FALSE)
Inbreed_SC_14_ngs <- read_table("scout_14_inbreeding_outbred", col_names = FALSE)
Inbreed_SC_15_ngs <- read_table("scout_15_inbreeding_outbred", col_names = FALSE)
Inbreed_SC_17_ngs <- read_table("scout_17_inbreeding_outbred", col_names = FALSE)


SC_13_vec_ngs=c(Inbreed_SC_13_ngs$X3)
SC_14_vec_ngs=c(Inbreed_SC_14_ngs$X3)
SC_15_vec_ngs=c(Inbreed_SC_15_ngs$X3)
SC_17_vec_ngs=c(Inbreed_SC_17_ngs$X3)

All_data_ngs=list( SC_13_vec_ngs, SC_14_vec_ngs, SC_15_vec_ngs , SC_17_vec_ngs )


cols=c("#9F2305", "#edea18", "#819d4e", "#155084")




pdf("Scout_inbreeding_ngs_relate_outbred.pdf", width = 12.4, height = 8)
par(mar = c(5, 4, 4, 2))
vioplot(All_data_ngs, names = c( "SC2013", "SC2014", "SC2015", "SC2017"), 
        col = cols, 
        xlab ="Year", 
        ylab="Inbreeding coefficient (F)", 
        main = "Estimated inbreeding coefficient for each Scout timepoint", 
        cex.lab=34)

dev.off()


####

############# ngsrelate using sc2013 AFs ########
library(ggplot2)
library(tidyr)
library(readr)
library(vioplot)

Inbreed_SC_13_ngs <- read_table("/gpfs/home/akwakye/inbreeding/SC_2013_af_filtered_inbreeding", col_names = FALSE)
Inbreed_SC_14_ngs <- read_table("/gpfs/home/akwakye/inbreeding/SC_2014_af_filtered_inbreeding", col_names = FALSE)
Inbreed_SC_15_ngs <- read_table("/gpfs/home/akwakye/inbreeding/SC_2015_af_filtered_inbreeding", col_names = FALSE)
Inbreed_SC_17_ngs <- read_table("/gpfs/home/akwakye/inbreeding/SC_2017_af_filtered_inbreeding", col_names = FALSE)
Inbreed_SC_20_ngs <- read_table("/gpfs/home/akwakye/inbreeding/Scout_2020_inbreeding_oubred", col_names = FALSE)
#Inbreed_Lob_21_ngs <- read_table("/gpfs/home/akwakye/inbreeding/Loberg_21_inbreeding_oubred", col_names = FALSE)
#Inbreed_Lob_21_ngs <- read_table("/gpfs/home/akwakye/inbreeding/Loberg_21_using_af_calculated_from_gatk", col_names = FALSE)


SC_13_vec_ngs=c(Inbreed_SC_13_ngs$X3)
SC_14_vec_ngs=c(Inbreed_SC_14_ngs$X3)
SC_15_vec_ngs=c(Inbreed_SC_15_ngs$X3)
SC_17_vec_ngs=c(Inbreed_SC_17_ngs$X3)
SC_20_vec_ngs=c(Inbreed_SC_20_ngs$X3)
Lob_21_vec_ngs=c(Inbreed_Lob_21_ngs$X3)

All_data_ngs=list( SC_13_vec_ngs, SC_14_vec_ngs, SC_15_vec_ngs , SC_17_vec_ngs, SC_20_vec_ngs )


cols=c("#9F2305", "#edea18", "#819d4e", "#155084", '#A4568B')


pdf("Scout_inbreeding_ngs_relate_sc13_AF_all_2.pdf", width = 8.4, height = 8.4)
par(mar = c(5, 4, 4, 2),  family = "Helvetica")
vioplot(All_data_ngs, names = c( "SC2013", "SC2014", "SC2015", "SC2017",  'SC2020'), 
        col = cols, 
        xlab ="Year", 
        ylab="Inbreeding coefficient (F)", 
        main = "Estimated inbreeding coefficient for each Scout timepoint", 
        cex.lab = 2,
        cex.axis = 1.2 
        )

dev.off()



