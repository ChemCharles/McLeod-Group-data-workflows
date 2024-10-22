# This data workflow was made by and is the intellectual property of C.C.J.Fitzgerald and the McLeod group (Australian National University, Canberra, Australia, 2602)

# For citation: if you use this code in any form please ask the appropriate authors for premission of use. 

# Disclaimer of liability.
# The use of this code is at the sole risk of the participant.
#======================================================================================================================
# Export files for Metaboanalyst 5.0 ####
# This code is designed to take the output from the Kmeans data and re-format it so it can be used in MetaboAnalyst 5.0. 
#======================================================================================================================
# Load Packages ####
library(tidyverse)
library(readxl)
#======================================================================================================================
# Load the data ####
dir()
data.all<-read_excel("data_R_H5.xlsx", sheet=1)
#======================================================================================================================
# Format and Structure data ####

# Subset peak intensity data
# This is the data that will go into MetaboAnalyst 5.0. 
DataMA <- subset(data.all[c(7:62)])

# Specify the names of the data groups for MetaboAnalyst 5.0 below:
Names <- c("Ta/-48h",	"Tb/-24h",	"Tc/0h",	"Td/4h",	"Te/8h","Tf/36h",	"Tj/72h","Tk/168h",	"Tl/192h",	"Tm/240h","Tn/288h","To/336h","Tp/432","Tq/Con1","Tr/Con2","Ts/QC")

# Specify the number of replicates per grouping that appears above:  
Replicates <- c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,11)

# specify the structure of the data set e.g. replicates. 
Groupings<-as.vector(factor(
  rep(Names, times=c(Replicates),levels = c(Names)), levels = Names))

#======================================================================================================================
#Put data in MetaboAnalyst 5.0 Format  ####
# condense data set to peak area info and concatenate the compound row_id, Retention time and m/z, this is used as a unique identifier (foreign key) for the data.  

DataMA <- rbind(Groupings,DataMA)
Samples <-  c("Group", as.vector(paste(data.all$row_id,"/",data.all$`RT [min]`,"/",data.all$`m/z`)))
DataMA$Samples <- Samples
# Makes 'Samples' the first column
DataMA_out <- DataMA[ , c("Samples",
                      names(DataMA)[names(DataMA) != "Samples"])]
#======================================================================================================================
# Write '.csv' file that you can put straight into MetaboAnalyst 5.0 web-browser.This can be found @ https://www.metaboanalyst.ca/ 

write.csv(Data_out_MA,"H5_data_all_trial.csv", row.names=FALSE)
#======================================================================================================================
# End of Script ####
#======================================================================================================================
# Post Script ####

#Instructions for MetaboAnalyst 
# The output data file is formatted to 'fit' into MetaboAnalyst5.0's Statistical Analysis (One factor) GUI on their website. 
# In this GUI mutliple analysis' can be done such as, PCA, PLS, OPLS-DA, Heatmaps, dendrograms, volcano-plots, and T-test.    
#======================================================================================================================
