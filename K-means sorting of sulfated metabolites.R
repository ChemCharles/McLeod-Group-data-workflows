# This data workflow was made by and is the intellectual property of C.C.J.Fitzgerald and the McLeod group (Australian National University, Canberra, Australia, 2602)

#This workflow describes the use of k-means clustering to differentiate putative sulfated and non-sulfated features in an un-targeted metabolomics mass spectrometry assay. For further details please see https://doi.org/10.3389/fmolb.2022.829511

# For citation: if you use this code in any form please ask the appropriate authors for premission of use. 

# Disclaimer of liability.
# The use of this code is at the sole risk of the participant.

#=============================================================================
#=============================================================================

#Set up data###

# Load in packages ####
library(readxl)
library(ggplot2)
library(factoextra)
#=============================================================================
#Set up directory, input data in excel sheet format

getwd()

dir() 

#enter the names of the excel file below from the output

data.all<-read_excel(".xlsx", sheet=1)

data.all[is.na(data.all)] <- 0

#=============================================================================
#=============================================================================

# k-means clustering of sulfated and non-sulfated features ####

#this line first subsets the 8 variables  into a data frame called Kmdata that will be used to cluster sulfated and non-sulfated features. Specifically these are, the intensities for the two sulfate product ions and the two neutral losses, and the calculated variables for IR average and max count

Kmdata <- data.all[,c(".SO3-_intensity","HSO3-_intensity",
                      ".SO4-_intensity","HSO4-_intensity",
                      "NeutralLoss_SO3_intensity","NeutralLoss_H2SO4_intensity",
                      "IRaverage","Max_Count")]

#this line executes the k-means clustering using the k-means function in base r, and puts the results into a data frame called Kmd. There are a couple of adjustable parameters, but in general you want to set clustering to 2 and maximum iterations to 10000. 

Kmd <- kmeans(Kmdata,2,iter.max = 10000)

# this then adds the clustering values (grouped as either 1 or 2) to your original data set. 

data.all$cluster=as.character(Kmd$cluster)

# this will generate a scatter plot of all your results. You can use this to quickly check the clustering across all of your variables. 

plot(Kmdata, col = data.all$cluster)
points(data.all$centers, col = 1:2, pch = 8, cex = 2)


#scatter plot with K-means clustering####

Kmeansplot<-
  
  ggplot(data=data.all, aes(data.all$Max_Count,data.all$IRaverage,color = factor(cluster)))+
  geom_point()+
    scale_color_manual(name="Compound type", values=c("Black","orange"), labels = c('Non-sulfates','Sulfates'))+
  labs(x="Maximum Abundance (%)", y ="Intensity Ratio (%)" )+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5, size = 14)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"))+
  expand_limits(x=c(0,100), y=c(0, 100))


print (Kmeansplot)


# Export k-means clustering data as a .csv file.  

write.csv(data.all,"H5_data_kmeans_clustered.csv")

#=============================================================================
#=============================================================================

# k-means checking, descriptive statistics and violin plots of k-means clustering ####

#Average silhouette method####

# This is a good check to see if you have the optimal number of clusters in k-means for the given data set. 

#The average silhouette approach measures the quality of a clustering, that is it determines how well each object lies within its cluster.

#It computes the average silhouette of observations for different values of k. The optimal number of clusters k is the one that maximizes the avearage silhouette.
#over a range of possible values of k.

#Kaufman, L. and Rousseeuw, P.J., 2009. Finding groups in data: an introduction to cluster analysis (Vol. 344). John Wiley & Sons.

fviz_nbclust(Kmdata, kmeans, method = "silhouette")

# source: https://uc-r.github.io/kmeans_clustering 

#=============================================================================

# Summary of descriptive stats for sulfates and non-sulfates#### 

Des_stats<- data.all %>% group_by(cluster) %>% 
  select(Max_Count,IRaverage) %>% 
  summarise_all(list(mean = mean,
                     median = median,
                     std = sd, 
                     min = min, 
                     max = max))


# Interquartile range for sulfates and non-sulfates 

# For maximum abundance (MA)

tbl_MA <- data.all %>%  group_by(cluster) %>% 
  summarise(SD = sd(Max_Count),
            Mean = mean(Max_Count),
            Median = as.numeric(median(Max_Count)),
            "Trimmed Mean" = mean(Max_Count, trim = 0.2),
            "Geometric Mean" = psych::geometric.mean(Max_Count),
            "Harmonic Mean" = psych::harmonic.mean(Max_Count),
            IQR = IQR(Max_Count),
            "%25 Q" = quantile(Max_Count, .25),
            "%50 Q" = quantile(Max_Count, .5),
            "%75 Q" = quantile(Max_Count, .75))

# For intensity ratio (IR)

tbl_IR <- data.all %>%  group_by(cluster) %>% 
  summarise(SD = sd(IRaverage),
            Mean = mean(IRaverage),
            Median = as.numeric(median(IRaverage)),
            "Trimmed Mean" = mean(IRaverage, trim = 0.2),
            "Geometric Mean" = psych::geometric.mean(IRaverage),
            "Harmonic Mean" = psych::harmonic.mean(IRaverage),
            IQR = IQR(IRaverage),
            "%25 Q" = quantile(IRaverage, .25),
            "%50 Q" = quantile(IRaverage, .5),
            "%75 Q" = quantile(IRaverage, .75))

#Export descriptive statistics ####

write.xlsx(tbl_IR,".xlsx")
write.xlsx(tbl_MA,".xlsx")


# Violin plot of data####

#Max_Count Product Ion abundance
P_MA <- ggplot(data.all, aes(x= data.all$cluster, y=data.all$`Max_Count`, fill=data.all$cluster)) +
  geom_violin(width=0.4, trim = TRUE) +
  geom_boxplot(width=0.1, alpha = .1)  +
  theme(legend.position="none",
        plot.title = element_text(size=11)) +
  ggtitle("Product Ion Abundance") +
  ylab ("Maximum Abundance (%)")+
  xlab("") +
  scale_x_discrete(breaks=c("1","2"),
                   labels=c("Non-Sulfate\nn=", "Sulfate\nn=")) +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=c("grey", "orange"),
                    name="Compound type",
                    breaks=c("1", "2"),
                    labels=c("Non-Sulfate", "Sulfate"))+
  coord_flip()

P_MA

#theme(legend.position = "none")

#Intensity Ratio
P_IR<- ggplot(data.all, aes(x= data.all$cluster, y=data.all$IRaverage, fill=data.all$cluster)) +
  geom_violin(width = .5, trim = TRUE) +
  geom_boxplot(width=.1, alpha = .1)  +
  theme(legend.position="none",
        plot.title = element_text(size=11)) +
  ggtitle("Proportion of Sulfate Ions") +
  ylab ("Intensity Ratio (%)")+
  xlab("") +
  scale_x_discrete(breaks=c("1","2"),
                   labels=c("Non-Sulfate\nn=", "Sulfate\nn=")) +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=c("grey", "orange"),
                    breaks=c("1", "2"),
                    labels=c("Non-Sulfate", "Sulfate"))
P_IR



