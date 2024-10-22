# This data workflow was made by and is the intellectual property of C.C.J.Fitzgerald and the McLeod group (Australian National University, Canberra, Australia, 2602)
# For citation: if you use this code in any form please ask the appropriate authors for premission of use. 

# Disclaimer of liability.
# The use of this code is at the sole risk of the participant.
#========================================================================================
# Limma analysis####

x <- c("Tm48h","Tm24h","T0h",	"T4h",	"T8h",	"T36h",	"T72h",
       "T168h",	"T192h","T240h",	"T288h",	"T336h", "T432h","Con_1","Con_2","QC")

Group<-factor(
  rep(x,
      times=c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,11),levels = c(x)), levels = x)

levels(Group)

design1<-model.matrix(~0+Group)
colnames(design1)<-x

fit <- lmFit(datalmScaled, design=design1)

#========================================================================================
# Look at contrasts between all groups relative to the zero hour time point. 
cont.matrix <- makeContrasts(T0hvsTm48h  =Tm48h-T0h,
                             T0hvsTm24h  =Tm24h-T0h,
                             T0hvsT4h  =T4h-T0h,
                             T0hvsT8h  =T8h-T0h,
                             T0hvsT36h  =T36h-T0h,
                             T0hvsT72h =T72h-T0h,
                             T0hvsT168h =T168h-T0h,
                             T0hvsT192h =T192h-T0h,
                             T0hvsT240h =T240h-T0h,
                             T0hvsT288h =T288h-T0h,
                             T0hvsT336h=T336h-T0h,
                             T0hvsT432h=T432h-T0h,
                             T0hvsCon_1=Con_1-T0h,
                             T0hvsCon_2=Con_2-T0h,
                             T0hvsQC =QC-T0h,
                             levels=design1)
#========================================================================================
# Isolate results as fit1
fit1 <- contrasts.fit(fit, cont.matrix)
fit1 <- eBayes(fit1)
summa.fit1<-decideTests(fit1)
summary(summa.fit1)
#========================================================================================
# Add p-value and mz's as columns to fit1
topTable(fit1, coef=1) 
fit1$p.value_log<- -1*log10(fit1$p.value)
fit1$mz<-data.all$`m/z`
str(fit1)
names(fit1)
#========================================================================================

