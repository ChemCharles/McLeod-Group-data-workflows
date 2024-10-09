#last accessed 9 March 2022
library(readxl)
library(ggplot2)
library(limma)
library(Glimma)
library(lattice)
library(wrapr)
library(dplyr)
library(factoextra)
library(ggfortify)
library("factoextra")
library("scatterplot3d")
library(plotly)
library(plyr)
library(tibble)
library(openxlsx)

#Set up directory, input data in excel sheet format####
#getwd()
setwd("C:/Users/cjjfi/OneDrive - Australian National University/2022_Jan_ARFL/Jan2022_H1/EUS_H1_UTWF/EUS_H1_IR_MA")
#dir()

#Py_1 Add file name of the output of intensity alignment####

data1<-read_excel("20220306_H1_Py1_output.xlsx", sheet=1)
#data1[is.na(data1)] <- 0

#Py_2 Add file name of the output of normalised intensity alignment #### 
data2<-read_excel("20220306_H1_Py2_output.xlsx", sheet=1)
#data2[is.na(data2)] <- 0

#Data wrangling####

#make sure they are data frames
data1 <-as.data.frame(data1)
data2 <-as.data.frame(data2)

#name first column row_id
colnames(data1)[1] <- "row_id"

colnames(data2)[1] <- "row_id"

#order data by row_id_ this way everything will align later on, so very important step!
data1 <- data1[order(data1$row_id),]

data2 <- data2[order(data2$row_id),]

#I like to check in excel whether everything is sorted

#write.csv(data1,"checkrow_id_data1.csv")

#write.csv(data2,"checkrow_id_data2.csv")

#________Py1__________####
#Calculation of intensity ratio####
colnames(data1)

#subset aligned MS2 data 
data1a<- subset(data1, select = 97:1430) #aligned MS2 data here
data1b<-subset(data1, select = 1:97) #MS-Dial output here
names(data1a) #check


#input names of all the column headers for one frame of ms2 data
Coltitles <- c("row_id",
              "AnalyteFullAnnotation",
              "DataFile",
              "RT [min]",
              "Molecular Weight",
              "m/z",
              "z",
              "Precursor Intensity",
              ".SO3-_intensity",
              "HSO3-_intensity",
              ".SO4-_intensity",
              "HSO4-_intensity",
              "Pre Int",
              "IonLoss_.SO3-_intensity",
              "IonLoss_HSO3-_intensity",
              "IonLoss_.SO4-_intensity",
              "IonLoss_HSO4-_intensity",
              "NeutralLoss_SO3_intensity",
              "NeutralLoss_H2SO4_intensity",
              "Sulfate Int",
              "Total Reporter Intensity (with INL)",
              "Total Int",
              "Internal Neutral Loss Transitions")

#give all columns the same name for each aligned MS2 run
Group<-factor(rep(c(Coltitles), times =58), #times indicates how many mgf files you have aligned
              levels = c(Coltitles))

colnames(data1a)<- Group #name all the columns the same thing for each set in the aligned files

#split all MS2 data into a list containing n MS2 spectra
data1c<-lapply(seq(1, ncol(data1a), by=23), function(i) # by = n indicates the number of columns per files
  data1a[i: pmin((i+22), ncol(data1a))]) #i+(n-1) so it knows how many to iterate through before stopping

#view(data1c[[3]]) #check data

#IR calculation#### 
data1_IR=NULL #make a matrix that we can pour the data into later

#this function will search through the list data1c and then run the below function based on the column names, 
#it then will out put the calculations into a new matrix

for (i in data1c) {
  i$ratio=((i$"Sulfate Int"/(i$"Total Int"-i$"Pre Int"))*100)
  data1_IR<-cbind(data1_IR,i$ratio)
  }

#Average ratio####
#by making all 0's NA the means ignore the 0's getting a 'truer' idea of the ratio

data1_IR[is.na(data1_IR)] <- 0
data1_IR<-as.data.frame(data1_IR)
data1_IR[data1_IR == 0] <- NA
#below will take the average of all the ratios of real values i.e. it will ignore 0's or in this case Na's
data1_IR$IRaverage <- rowMeans(data1_IR, na.rm = TRUE)
data1_IR[is.na(data1_IR)] <- 0

data1_IR$row_id<-factor(data1b$row_id)

data1_ex<-data1_IR[,c("row_id","IRaverage")]

#write.csv(data1_IR, "check_ave_ratio_matrix.csv")

#data1b$IR_Average<-data1_IR$IRaverage


#_________Py2__________####
#Calculate max abundance from normalised intensities####

#subset to only get MS2 data
data2a<- subset(data2, select = 97:1430)

data2b <- data2a

#names all col names by the group factor established earlier
colnames(data2b)<- Group

#makes all characters into numeric, will return NA on all character columns. this is fine as we are not interested in them

data2b <- sapply(data2b, as.numeric)
data2b <- as.data.frame(data2b)
data2b[is.na(data2b)] <- 0
data2b[data2b == 0] <- NA

#split all MS2 data into a list containing n MS2 spectra
data2c<-lapply(seq(1, ncol(data2b), by=23), function(i) 
  data2b[i: pmin((i+22), nrow(data2b))])

#this will get rid of all those annoying character columns.

dataT <- data2b[ , -which(names(data2b) %in% c("AnalyteFullAnnotation",
                                               "DataFile",
                                               "RT [min]",
                                               "Molecular Weight",
                                               "m/z",
                                               "z",
                                               "Internal Neutral Loss Transitions"))]

#split all MS2 data into a list containing n MS2 spectra
data2c<-lapply(seq(1, ncol(dataT), by=16), function(i) 
  dataT[i: pmin((i+15), nrow(dataT))])

check5 <- data2c[[58]]

#https://stackoverflow.com/questions/31465415/combine-multiple-data-frames-and-calculate-average 

#this function takes the mean over all of the list data2c. 
data2_ave <- aaply(laply(data2c, as.matrix), c(2, 3), function (x) mean(x, na.rm = TRUE))

data2d <- as.data.frame(data2_ave)

data2d<-data2d[ , -which(names(data2d) %in% c("row_id"))] #get rid of the incorrect row_id's

#put in the row id's from MSdial

data2d$row_id<-data2$row_id
data2d[is.na(data2d)] <- 0
data2d<-data2d[,c(ncol(data2d),1:(ncol(data2d)-1))]

#calculate the max count 
data2d$Max_Count <- pmax(data2d$`HSO3-_intensity`,data2d$`HSO4-_intensity`,
                          data2d$`.SO3-_intensity`,data2d$`.SO4-_intensity`,
                          data2d$NeutralLoss_H2SO4_intensity,data2d$NeutralLoss_SO3_intensity)

#also get rid of the ionloss data this doesn't really tell us anything
data2d<-data2d[ , -which(names(data2d) %in% c("IonLoss_.SO3-_intensity",
                                              "IonLoss_HSO3-_intensity",
                                              "IonLoss_.SO4-_intensity",
                                              "IonLoss_HSO4-_intensity"))]

#________Output_________####
#put it all back together####

#merge everything based the the MS-dial row_id

output<-merge.data.frame(data1b,data2d, by = "row_id")

output<-merge.data.frame(output,data1_ex, by = "row_id")

#write the output as the below csv file

#write.csv(output,"data_R_H1.csv")

write.xlsx(output, 'data_R_H1.xlsx')

#