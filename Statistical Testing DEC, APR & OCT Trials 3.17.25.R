setwd()
#install pkgs and load libraries
install.packages("ggplot2") #use binaries not source code
install.packages("ggpubr")
install.packages(c("tidyverse", "broom", "AICcmodavg"))
install.packages("PMCMRplus")
#install rstatix for Wilcoxon testing options
install.packages("rstatix")
install.packages("stats")
install.packages("lawstat")
install.packages('coin')
library('coin')
library("lawstat")
library("rstatix")
library("stats")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("PMCMRplus")
library("tidyverse")
library("broom") 
library("AICcmodavg")

#BZ DEC friedman test, compare three conditions 
BZ.friedman=read.csv("BZ DEC Trial Table.csv", header = F)
names(BZ.friedman)=c("Trial1", "Trial2", "Trial3")
#histograms for each trial to evaluate normality
hist(BZ.friedman$Trial1)

hist(BZ.friedman$Trial2)

hist(BZ.friedman$Trial3)


#will throw error w/o rownames specified
rownames(BZ.friedman)=c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
BZ.friedman.matrix=as.matrix(BZ.friedman)
BZ.friedman.matrix
friedman.test(BZ.friedman.matrix)

median(BZ.friedman$Trial1)
#13.9
median(BZ.friedman$Trial2)
#12.6
median(BZ.friedman$Trial3)
#22.4

boxplot(BZ.friedman$Trial1,BZ.friedman$Trial2,BZ.friedman$Trial3,names=c("Trial 1","Trial 2","Trial 3"))

#source: stcp-marquier-Friedman R test .pdf page 3 (in stat files on desktop)

#function no longer available
#posthoc.friedman.nemenyi.test(BZ.friedman.matrix)

#use this instead also from PMCMRplus:
nemenyi.BZ.DEC=frdAllPairsNemenyiTest(BZ.friedman.matrix) 
#now can do post hoc testing to find out which pairs are sig. different!
#print the test statistic
nemenyi.BZ.DEC
print(nemenyi.BZ.DEC$statistic)
#          Trial1   Trial2
#Trial2 1.507557       NA
#Trial3 2.412091 3.919647

BZ.friedman.matrix

# Trial 2 is diff. than 3, p = 0.02, df = 2

#move on to EF DEC, HgT conc. from overlying water (0cm) to 20cm, n = 11
EF.friedman=read.csv("EF DEC Trial Table.csv", header = F)
names(EF.friedman)=c("Trial1", "Trial2", "Trial3")

hist(EF.friedman$Trial1)

hist(EF.friedman$Trial2)

hist(EF.friedman$Trial3)


#histograms for both sites all in one figure 
#
png("Distributions Trials 1-3 EF and BZ.png", units='in',width=28,height=10,res=200)
par(mfrow=c(2,3))
hist(EF.friedman$Trial1,col="light pink",freq=FALSE, breaks =  5, main=" Distribution EF DEC Trial 1")
hist(EF.friedman$Trial2,col="light blue",freq=FALSE, breaks =5, main=" Distribution EF DEC Trial 2")
hist(EF.friedman$Trial3,col="light green",freq=FALSE, breaks =5, main="Distribution EF DEC Trial3")
hist(BZ.friedman$Trial1,col="light pink",freq=FALSE, breaks =  5, main=" Distribution BZ DEC Trial 1")
hist(BZ.friedman$Trial2,col="light blue",freq=FALSE, breaks =5, main="Distribution  BZ DEC Trial 2")
hist(BZ.friedman$Trial3,col="light green",freq=FALSE, breaks =5, main="Distribution  BZ DEC Trial 3")
dev.off()

rownames(EF.friedman)=c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
EF.friedman.matrix=as.matrix(EF.friedman)
friedman.test(EF.friedman.matrix)
EF.friedman.matrix

median(EF.friedman$Trial1)
#28.75
median(EF.friedman$Trial2)
#34.85
median(EF.friedman$Trial3)
#38.49

boxplot(EF.friedman$Trial1,EF.friedman$Trial2,EF.friedman$Trial3,names=c("Trial 1","Trial 2","Trial 3"))

#post hoc nemenyi test to tell which group means are stat. sig. diff.
nemenyi.EF.DEC = frdAllPairsNemenyiTest(EF.friedman.matrix)

nemenyi.EF.DEC
#P value adjustment method: single-step
#print the test statistic
print(nemenyi.EF.DEC$statistic)
        #Trial1    Trial2
#Trial2 1.507557       NA
#Trial3 3.919647 2.412091
print(nemenyi.EF$dist)
#"q"
# p-values are computed from the studentized range distribution

#move on to EF OCT, n=4, HgT conc. from depth 4-10cm
EF.friedman.OCT=read.csv("EF OCT Trial Table.csv", header = F)
names(EF.friedman.OCT)=c("Trial6", "Trial7", "Trial8")
hist(EF.friedman.OCT$Trial6)

hist(EF.friedman.OCT$Trial7)

hist(EF.friedman.OCT$Trial8)

rownames(EF.friedman.OCT)=c("4", "6", "8", "10")
EF.friedman.matrix.OCT=as.matrix(EF.friedman.OCT)
friedman.test(EF.friedman.matrix.OCT)


median(EF.friedman.OCT$Trial6)
#1.025
median(EF.friedman.OCT$Trial7)
#-3
median(EF.friedman.OCT$Trial8)
#15.91

#posthoc test for EF OCT
nemenyi.EF = frdAllPairsNemenyiTest(EF.friedman.matrix.OCT)
nemenyi.EF
#P value adjustment method: single-step
#print the test statistic
print(nemenyi.EF$statistic)
#         Trial6  Trial7
#Trial7    1.0     NA
#Trial8    2.5    3.5
#print the critical value
print(nemenyi.EF$dist)
#"q"
# p-values are computed fromthe studentized range distribution
#source: http://www.sthda.com/english/wiki/paired-samples-wilcoxon-test-in-r
##OK move on to non-parametric Wilcoxon Signed-Rank test 

# BZ April read in dataframe
BZ.wilcox=read.csv("BZ APR Trial Table.csv", header = F)
names(BZ.wilcox)=c("Trial4", "Trial5")

#compare two conditions, two sided alternative 
##for BZ April, non-parametric Wilcoxon Signed-Rank test on BZ April data, n=11
wilcox.test(BZ.wilcox$Trial4,BZ.wilcox$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)

median(BZ.wilcox$Trial4)
median(BZ.wilcox$Trial5)
mean(BZ.wilcox$Trial4)
mean(BZ.wilcox$Trial5)


##for EF April dataframe read in
EF.wilcox=read.csv("EF Trial Table APR ORIG.csv", header = F)
names(EF.wilcox)=c("Trial4", "Trial5")
EF.wilcox.matrix=as.matrix(EF.wilcox)

##for EF April, Wilcoxon signed rank test 
wilcox.test(EF.wilcox$Trial4,EF.wilcox$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)


# Transform into long data: 
# gather the before and after values in the same column
EF.APR.long <- EF.wilcox %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(EF.APR.long, 10)
#wow it worked!!!!
#summary stats?
EF.APR.long %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")
#not working
bxp <- ggpaired(EF.APR.long, x = "trial", y = "pM", 
                order = c("Trial4", "Trial5"),
                ylab = "[pM]", xlab = "April East Fork Trials")
bxp


#
EF.Apr.diff <- EF.wilcox %>% mutate(differences = Trial5 - Trial4)
gghistogram(EF.Apr.diff, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)
#make a matrix just to make R happy
EF.APR.matrix=as.matrix(EF.Apr.diff)

#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
symmetry.test(EF.APR.matrix, boot = F) #Symmetry test by Miao, Gel, and Gastwirth (2006)

#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
#fine to use Wilcoxon signed rank in this case then: Test statistic = 1.0891, p-value = 0.2761

#do a Shapiro-Wilk test also on the differences to check the distribution
# p > 0.05 means it is normally distributed
shapiro.test(EF.APR.matrix)
# not normally distributed , use non-parametric test
# W = 0.95987, p-value = 0.2563
# 
stat.test=
  wilcox.test(x=EF.wilcox$Trial4, y=EF.wilcox$Trial5, paired = TRUE) %>%
  add_significance()
stat.test

EF.APR.long %>%
  wilcox_effsize(pM ~ trial, paired = TRUE)

#####BZ section########################
#####
BZ.wilcox=read.csv("BZ APR Trial Table.csv", header = F)
names(BZ.wilcox)=c("Trial4", "Trial5")


# Transform into long data: 
# gather the before and after values in the same column
BZ.APR.long <- BZ.wilcox %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(BZ.APR.long, 10)
#wow it worked!!!!
#summary stats?
BZ.APR.long %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")

bxp <- ggpaired(BZ.APR.long, x = "trial", y = "pM", 
                order = c("Trial4", "Trial5"),
                ylab = "[pM]", xlab = "April Beach Zone Trials")
bxp

#
BZ.Apr.diff <- BZ.wilcox %>% mutate(differences = Trial5 - Trial4)
gghistogram(BZ.Apr.diff, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

BZ.APR.matrix=as.matrix(BZ.Apr.diff)
#make a matrix 
#
#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
symmetry.test(BZ.APR.matrix, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
#fine to use Wilcoxon signed rank in this case

#check distribution with Shapiro-Wilk also
shapiro.test(BZ.APR.matrix)
# not normally distributed , use non-parametric test
# W = 0.97364,  p-value = 0.5868

stat.test=
  wilcox.test(x=BZ.wilcox$Trial4, y=BZ.wilcox$Trial5, paired = TRUE) %>%
  add_significance()
stat.test

BZ.effsize=BZ.APR.long %>%
  wilcox_effsize(pM ~ trial, paired = TRUE)


