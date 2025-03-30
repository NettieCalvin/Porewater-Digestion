setwd()

#load libraries
library("PMCMRplus")
library("rstatix")
citation("PMCMRplus")
citation("rstatix")

#read in original data for EF OCT
EF.friedman.OCT=read.csv("EF OCT Trial Table.csv", header = F, nrow = 4)
names(EF.friedman.OCT)=c("Trial6", "Trial7", "Trial8")
rownames(EF.friedman.OCT)=c("4", "6", "8","10")
EF.friedman.matrix.OCT=as.matrix(EF.friedman.OCT)
print(EF.friedman.OCT)

friedmanTest(EF.friedman.matrix.OCT) #Friedman chi-squared = 6.5, df = 2, p-value = 0.03877
EF.OCT.Nemenyi = frdAllPairsNemenyiTest(EF.friedman.matrix.OCT)
       #Trial6 Trial7
#Trial7 0.759  -     
#Trial8 0.181  0.036 
print(EF.OCT.Nemenyi)
print(EF.OCT.Nemenyi$statistic)


EF.effsize =read.csv("EF.effsize.OCT.csv", header =F)
names(EF.effsize)=c("pM", "trial", "depth")
print(EF.effsize)

EF.effsize %>% friedman_effsize(pM ~ trial | depth, ci=T, ci.type = "perc", conf.level = 0.95, nboot = 1000)
#.y.       n effsize conf.low conf.high method    magnitude
  #pM        4   0.812     0.75         1 Kendall W large   

#histograms of EF OCT data by trial
png("Distributions Trials 6-8 at EF.png", units='in',width=28,height=5,res=200)
par(mfrow=c(1,3))
hist(EF.friedman.OCT$Trial6,col="light pink",freq=FALSE, breaks =  5, main="Distribution EF OCT Trial 6")
hist(EF.friedman.OCT$Trial7,col="light blue",freq=FALSE, breaks =5, main="Distribution EF OCT Trial 7")
hist(EF.friedman.OCT$Trial8,col="light green",freq=FALSE, breaks =5, main="Distribution EF OCT Trial 8")

dev.off()

#source: http://www.sthda.com/english/wiki/paired-samples-wilcoxon-test-in-r


####ZERO substituted for BLOD data
EF.friedman.OCT.0=read.csv("EF OCT Trial Table ZEROS.csv", header = F, nrow = 4)
names(EF.friedman.OCT.0)=c("Trial6", "Trial7", "Trial8")
rownames(EF.friedman.OCT.0)=c("4", "6", "8","10")
EF.friedman.matrix.OCT.0=as.matrix(EF.friedman.OCT.0)
EF.friedman.matrix.OCT.0

median(EF.friedman.OCT.0$Trial6)
#1.025
median(EF.friedman.OCT.0$Trial7)
#0
median(EF.friedman.OCT.0$Trial8)
#15.91
friedmanTest(EF.friedman.matrix.OCT.0)
#Friedman chi-squared = 6.5333, df = 2, p-value = 0.03813
EF.Nemenyi.ZEROS = frdAllPairsNemenyiTest(EF.friedman.matrix.OCT.0)
#        Trial6 Trial7
#Trial7 0.933  -     
#Trial8 0.126  0.056 
print(EF.Nemenyi.ZEROS$statistic)

#       Trial6 Trial7
#Trial7   0.50     NA
#Trial8   2.75   3.25

######effect size for EF OCT ZERO substituted for BLOD data ######################
EF.effsize.ZEROS =read.csv("EF.effsize.ZEROS.csv", header = F)
names(EF.effsize.ZEROS)=c("pM", "trial", "depth")
print(EF.effsize.ZEROS)

EF.effsize.ZEROS %>% friedman_effsize(pM ~ trial | depth, ci=T, ci.type = "perc", conf.level = 0.95, nboot = 1000)

#.y.       n effsize conf.low conf.high method    magnitude
#pM        4   0.817     0.75         1 Kendall W large  

# start section for HALF LOD substituted for BLOD data
EF.friedman.OCT.H=read.csv("EF OCT Trial Table HALFS.csv", header = F, nrow = 4)
names(EF.friedman.OCT.H)=c("Trial6", "Trial7", "Trial8")
rownames(EF.friedman.OCT.H)=c("4", "6", "8","10")
EF.friedman.matrix.OCT.H=as.matrix(EF.friedman.OCT.H)
EF.friedman.matrix.OCT.H

median(EF.friedman.OCT.H$Trial6)
#1.025
median(EF.friedman.OCT.H$Trial7)
#1.85
median(EF.friedman.OCT.H$Trial8)
#15.91

boxplot(EF.friedman.OCT.H$Trial6,EF.friedman.OCT.H$Trial7,EF.friedman.OCT.H$Trial8,names=c("Trial 6","Trial 7","Trial 8"))

friedmanTest(EF.friedman.matrix.OCT.H) #Friedman chi-squared = 8, df = 2, p-value = 0.01832
EF.Nemenyi.HALFS = frdAllPairsNemenyiTest(EF.friedman.matrix.OCT.H)
print(EF.Nemenyi.HALFS$p.value)
#      Trial6 Trial7
#Trial7 0.333  -     
#Trial8 0.013  0.333 

print(EF.Nemenyi.HALFS$statistic)
#         Trial6 Trial7
#Trial7      2     NA
#Trial8      4      2

######effect size for EF OCT HALF LOD substituted for BLOD data
EF.effsize.HALFS =read.csv("EF.effsize.HALFS.csv", header = F)
names(EF.effsize.HALFS)=c("pM", "trial", "depth")

EF.effsize.HALFS %>% friedman_effsize(pM ~ trial | depth, ci=T, ci.type = "perc", conf.level = 0.95, nboot = 10)
#"All values of t are equal to  1 \n Cannot calculate confidence intervals"

######## FULL LODS  ######
######## # start section for FULL LOD substituted for BLOD data
EF.friedman.OCT.L=read.csv("EF OCT Trial Table LODS.csv", header = F, nrow = 4)
print(EF.friedman.OCT.L)
names(EF.friedman.OCT.L)=c("Trial6", "Trial7", "Trial8")
rownames(EF.friedman.OCT.L)=c("4", "6", "8","10")
EF.friedman.matrix.OCT.L=as.matrix(EF.friedman.OCT.L)
print(EF.friedman.matrix.OCT.L)

median(EF.friedman.OCT.L$Trial6)
#1.025
median(EF.friedman.OCT.L$Trial7)
#1.9
median(EF.friedman.OCT.L$Trial8)
#15.91

boxplot(EF.friedman.OCT.L$Trial6,EF.friedman.OCT.L$Trial7,EF.friedman.OCT.L$Trial8,names=c("Trial 6","Trial 7","Trial 8"))

friedmanTest(EF.friedman.matrix.OCT.L) #Friedman chi-squared = 8, df = 2, p-value = 0.01832
EF.Nemenyi.FULL.LOD = frdAllPairsNemenyiTest(EF.friedman.matrix.OCT.L)
print(EF.Nemenyi.FULL.LOD$p.value)
#      Trial6 Trial7
#Trial7 0.333  -     
#Trial8 0.013  0.333 
print(EF.Nemenyi.FULL.LOD$statistic)


EF.effsize.FULL.LOD =read.csv("EF.effsize.FULL.LOD.csv", header = F)
names(EF.effsize.FULL.LOD)=c("pM", "trial", "depth")
print(EF.effsize.FULL.LOD)

EF.effsize.FULL.LOD %>% friedman_effsize(pM ~ trial | depth, ci=T, ci.type = "perc", conf.level = 0.95, nboot = 1000)
# "All values of t are equal to  1 \n Cannot calculate confidence intervals"

####################################################################################


