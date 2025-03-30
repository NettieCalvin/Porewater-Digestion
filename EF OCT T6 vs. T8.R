setwd()

library("ggplot2")
library("ggpubr")
library("lawstat")
library("rstatix")
library("stats")
library("tidyverse")
#

#use stat test:
EF.wilcox.OCT=read.csv("T6vT8.csv", header = F)
names(EF.wilcox.OCT)=c("Trial6", "Trial8") #name columns

#use stat test to check symmetry as a proxy for normal distibution:
#a large p value allows us to accept NULL which is "symmetric enough"
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"

EF.oct.diff <- EF.wilcox.OCT %>% mutate(differences = Trial6 - Trial8)
gghistogram(EF.oct.diff, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

#use stat test to check symmetry which is a proxy for normal distribution:
##Distribution of Differences: While the test does not assume a specific form for the overall distribution of measurements, it does concern itself with the distribution of differences between pairs. 
#The central tendency of these differences is what the test seeks to evaluate, under the hypothesis that the median difference is zero.
#source: https://www.statisticssolutions.com/free-resources/directory-of-statistical-analyses/assumptions-of-the-wilcox-sign-test/
#a large p value allows us to accept NULL which is "symmetric enough"
#Symmetry test by Miao, Gel, and Gastwirth (2006)
#do the MGG symmetry test on the differences
EF.oct.matrix=as.matrix(EF.oct.diff)
SYM1=symmetry.test(EF.oct.matrix, boot = F)
SYM1
#Test statistic = 0.32487, p-value = 0.7453

##for EF OCT T6 vs.T8 HgT concentrations "as calculated"
wilcox.test(EF.wilcox.OCT$Trial6[1:5],EF.wilcox.OCT$Trial8[1:5], 
            paired = T, 
            alternative = "less",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)
#V = 0, p-value = 0.03125

median(EF.wilcox.OCT$Trial6[1:5])
#1.05

median(EF.wilcox.OCT$Trial8[1:5])
#13.1

###effect size of as calculated values############################
# Transform into long data: 
# gather the before and after values in the same column
EF.6v8.long <- EF.wilcox.OCT %>%
  gather(key = "trial", value = "pM", Trial6, Trial8)
print(EF.6v8.long )

EF.T6vT8.effsize =  EF.6v8.long  %>%
  wilcox_effsize(pM ~ trial, paired = TRUE,  alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)
print(alt.effsize )

#.y.   group1 group2 effsize    n1    n2 conf.low conf.high magnitude
# pM    Trial6 Trial8   0.905     5     5      0.9      0.95 large  

###effect size of ZEROS values############################
EF.wilcox.OCT.0=read.csv("T6vT8.ZEROS.csv", header = F)
names(EF.wilcox.OCT.0)=c("Trial6", "Trial8") #name columns

#make a plot of the symmetry!
EF.oct0.diff <- EF.wilcox.OCT.0 %>% mutate(differences = Trial6 - Trial8)
gghistogram(EF.oct0.diff, x = "differences", y = "..density..", 
            fill = "pink",bins = 5, add_density = TRUE)

#MGG symmetry test:
EF.oct0.matrix=as.matrix(EF.oct0.diff)
SYM2=symmetry.test(EF.oct0.matrix, boot = F)

SYM2
#Test statistic = 0.30548, p-value = 0.76

##for EF OCT with ZERO substituted for BLOD data
wilcox.test(EF.wilcox.OCT.0$Trial6[1:5],EF.wilcox.OCT.0$Trial8[1:5], 
            paired = T, 
            alternative = "less",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)
#V = 1, p-value = 0.0625

median(EF.wilcox.OCT.0$Trial6[1:5])
#1.05

median(EF.wilcox.OCT.0$Trial8[1:5])
#13.1

#convert to long data:
## gather the before and after values in the same column
EF.6v8.ZEROS.long <- EF.wilcox.OCT.0 %>%
  gather(key = "trial", value = "pM", Trial6, Trial8)
print(EF.6v8.ZEROS.long)

EF.ZEROS.effsize=EF.6v8.ZEROS.long %>%
  wilcox_effsize(pM ~ trial, paired = TRUE,  alternative ="less", ci=TRUE, ci.type = "perc",  conf.level =0.95)

print(EF.ZEROS.effsize)
#.y.   group1 group2 effsize    n1    n2 conf.low conf.high magnitude
#pM    Trial6 Trial8   0.784     5     5     0.18      0.93 large  

###effect size of HALFS values############################
EF.wilcox.OCT.H=read.csv("T6vT8.HALFS.csv", header = F)
names(EF.wilcox.OCT.H)=c("Trial6", "Trial8") #name columns

#make a plot of the symmetry!
EF.octH.diff <- EF.wilcox.OCT.H %>% mutate(differences = Trial6 - Trial8)
gghistogram(EF.octH.diff, x = "differences", y = "..density..", 
            fill = "lightblue",bins = 5, add_density = TRUE)

##do the MGG symmetry test on the differences
EF.octH.matrix=as.matrix(EF.octH.diff)
SYM3=symmetry.test(EF.octH.matrix, boot = F)
SYM3
#Test statistic = 0.22447, p-value = 0.8224

##for EF OCT HALF LOD substituted for BLOD data
wilcox.test(EF.wilcox.OCT.H$Trial6[1:5],EF.wilcox.OCT.H$Trial8[1:5], 
            paired = T, 
            alternative = "less",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)
#V = 1, p-value = 0.0625

median(EF.wilcox.OCT.H$Trial6[1:5])
#1.05

median(EF.wilcox.OCT.H$Trial8[1:5])
#13.1

#convert to long data:
## gather the before and after values in the same column
EF.6v8.HALFS.long <- EF.wilcox.OCT.H %>%
  gather(key = "trial", value = "pM", Trial6, Trial8)
print(EF.6v8.HALFS.long)

EF.HALFS.effsize = EF.6v8.HALFS.long  %>%
  wilcox_effsize(pM ~ trial, paired = TRUE,  alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)

print(EF.HALFS.effsize)
#.y.   group1 group2 effsize    n1    n2 conf.low conf.high magnitude
  #pM    6      8        0.784     5     5     0.18      0.93 large   


###effect size of FULL LODS values############################
EF.wilcox.OCT.LOD=read.csv("T6vT8.LOD.csv", header = F)
names(EF.wilcox.OCT.LOD)=c("Trial6", "Trial8") #name columns

#make a  plot of the symmetry of the differences!
EF.octLOD.diff <- EF.wilcox.OCT.LOD %>% mutate(differences = Trial6 - Trial8)
gghistogram(EF.octLOD.diff, x = "differences", y = "..density..", 
            fill = "yellow",bins = 5, add_density = TRUE)

#do the MGG symmetry test on the differences
EF.octLOD.matrix=as.matrix(EF.octLOD.diff)
SYM4=symmetry.test(EF.octLOD.matrix, boot = F)
SYM4
#Test statistic = 0.32779, p-value = 0.7431
#
##for EF OCT LOD itself substituted for BLOD data
wilcox.test(EF.wilcox.OCT.LOD$Trial6[1:5],EF.wilcox.OCT.LOD$Trial8[1:5], 
            paired = T, 
            alternative = "less",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)
#V = 0, p-value = 0.03125

median(EF.wilcox.OCT.LOD$Trial6[1:5])
#1.05

median(EF.wilcox.OCT.LOD$Trial8[1:5])
#13.1

#convert to long data:
## gather the before and after values in the same column
EF.6v8.LOD.long <- EF.wilcox.OCT.LOD %>%
  gather(key = "trial", value = "pM", Trial6, Trial8)
print(EF.6v8.LOD.long)

EF.LOD.effsize=EF.6v8.LOD.long %>%
  wilcox_effsize(pM ~ trial, paired = TRUE,  alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)

print(EF.LOD.effsize)
#.y.   group1 group2 effsize    n1    n2 conf.low conf.high magnitude
#pM    Trial6 Trial8   0.905     5     5      0.9      0.95 large    
