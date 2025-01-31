# This data workflow was made by and is the intellectual property of C.C.J.Fitzgerald and the McLeod group (Australian National University, Canberra, Australia, 2602)
# For citation: if you use this code in any form please ask the appropriate authors for premission of use. 

# Disclaimer of liability.
# The use of this code is at the sole risk of the participant.
#========================================================================================
#Compound Searcher#### 

# this function searches through every m/z and looks to see if it is with in 5 ppm error of the given list of compounds

# ppm error; put in the error amount you want

ppm <- 5

# Lower-bound PPM function PPM = 5

NegPPM <- function(z){ 
  ((z*-ppm)/(10^6))+(z)
}

# Upper-bound PPM function PPM = 5
PosPPM <- function(z){ 
  ((z*ppm)/(10^6))+(z)
}

#========================================================================================
# List of possible Steroid sulfate####
Steroid_Theoretical_mz = c(#Put in steroid masses with a minimum of 4.d.p)

# Set list as a single column data frame 
# Note you can also just import a '.csv' list of compounds

datamz <- data.frame(Steroid_Theoretical_mz)

# Find lower bound, then add it as a column
datamz$lowerbound <- apply(datamz,c(1,2),NegPPM) #applys NegPPM function to a given dataframe

# Find upper bound, then add it as a column
datamz$upperbound <- apply(datamz,c(1,2),PosPPM) #applys PosPPM function to a given dataframe

# Make Steroid_datamz into a readable (conformable) data set

m <- as.matrix.data.frame(datamz)
colnames(m)[c(1,2,3)] <- c("Theoretical_mz","lowerbound","upperbound") #naming the three columns
m<-subset(m,select = -c(upperbound.lowerbound)) # an extra row is made that needs to be deleted 
m<-as.data.frame(m)
#========================================================================================
# Rename your data that you want to match to the above list 
n <- data.all

# Apply's function to all listed m/z in the data set and returns "Match" if their theoretical m/z is found
n$Theoretical_mz <- ifelse(sapply(n$`m/z`, function(p) 
  any(m$lowerbound <= p & m$upperbound >= p)),"Match", NA)

# Subsets list only including values where a match was found
Steroids <- subset(n,(!is.na(n$Theoretical_mz))) 

#========================================================================================
# Output of list of matched steroids 
write.csv(Steroids,"Matched_molecules.csv")

#========================================================================================
#End of script 
#========================================================================================
