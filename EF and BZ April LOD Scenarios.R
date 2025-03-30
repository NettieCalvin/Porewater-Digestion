setwd()

library("tidyverse")
library("ggplot2")
library("ggpubr")
library("rstatix")
library("stats")
library("lawstat")
library('coin')

#source: http://www.sthda.com/english/wiki/paired-samples-wilcoxon-test-in-r
##OK move on to non-parametric Wilcoxon Signed-Rank test 


##ORIGINAL calculated concentrations
####with all 11 depths incl.################
##for EF April df read in

EF.wilcox=read.csv("EF Trial Table APR ORIG.csv", header = F)
names(EF.wilcox)=c("Trial4", "Trial5")
EF.wilcox.matrix=as.matrix(EF.wilcox)
##for EF April
wilcox.test(EF.wilcox$Trial4,EF.wilcox$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)

#Wilcoxon signed rank exact test: V = 1, p-value = 0.001953
#
# Transform into long data: 
# gather the before and after values in the same column
EF.APR.long <- EF.wilcox %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(EF.APR.long, 11)
#wow it worked!!!!
#summary stats?
EF.APR.long %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")
#not working
bxp.EF.APR.ORIG <- ggpaired(EF.APR.long, x = "trial", y = "pM",
                            order = c("Trial4", "Trial5"),
                            ylab = "[pM]", xlab = "April East Fork Trials")
bxp.EF.APR.ORIG


#check the differences of the paired trials for normal distribution 
EF.Apr.diff <- EF.wilcox %>% mutate(differences = Trial5 - Trial4)
gghistogram(EF.Apr.diff, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

EF.APR.matrix=as.matrix(EF.Apr.diff)
EF.diff=EF.Apr.diff$differences

#use stat test to check symmetry which is a proxy for normal distribution:
##Distribution of Differences: While the test does not assume a specific form for the overall distribution of measurements, it does concern itself with the distribution of differences between pairs. 
#The central tendency of these differences is what the test seeks to evaluate, under the hypothesis that the median difference is zero.
#source: https://www.statisticssolutions.com/free-resources/directory-of-statistical-analyses/assumptions-of-the-wilcox-sign-test/
#a large p value allows us to accept NULL which is "symmetric enough"
#Symmetry test by Miao, Gel, and Gastwirth (2006)
SYM1=symmetry.test(EF.APR.matrix, boot = F)
SYM1$statistic #1.089069 
SYM1$p.value # 0.2761233 

#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"

symmetry.test(EF.diff, boot = F)
#fine to use Wilcoxon signed rank in this case

stat.test=
  wilcox.test(x=EF.wilcox$Trial4, y=EF.wilcox$Trial5, paired = TRUE) %>%
  add_significance()
stat.test

EF.APR.eff.table= EF.APR.long %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)

print(EF.APR.eff.table)

##ZEROS SUBSTITUTED for BLOD concentrations
##for EF April df read in
EF.wilcox.0=read.csv("EF Apr Trial Table ZEROS.csv", header = F)
names(EF.wilcox.0)=c("Trial4", "Trial5")

##for EF April
wilcox.test(EF.wilcox.0$Trial4,EF.wilcox.0$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)

# Transform into long data: 
# gather the before and after values in the same column
EF.APR.long.0 <- EF.wilcox.0 %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(EF.APR.long.0, 11)
#wow it worked!!!!
#summary stats?
EF.APR.long.0 %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")

bxp.EF.APR.ZEROS <- ggpaired(EF.APR.long.0, x = "trial", y = "pM",
                             order = c("Trial4", "Trial5"),
                             ylab = "[pM]", xlab = "April East Fork Trials ZEROS")
bxp.EF.APR.ZEROS

#
EF.Apr.diff.0 <- EF.wilcox.0 %>% mutate(differences = Trial5 - Trial4)
gghistogram(EF.Apr.diff.0, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

EF.APR.matrix.0=as.matrix(EF.Apr.diff.0)

#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
SYM2 = symmetry.test(EF.APR.matrix.0, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM2$statistic #1.083249 
SYM2$p.value #0.278698 

#fine to use Wilcoxon signed rank in this case
stat.test=
  wilcox.test(x=EF.wilcox.0$Trial4, y=EF.wilcox.0$Trial5, paired = TRUE) %>%
  add_significance()
stat.test

#V = 1, p-value = 0.001953

EF.APR.long.0 %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)


##HALF the LOD SUBSTITUTED for BLOD concentrations
##for EF April df read in
EF.wilcox.H=read.csv("EF Apr Trial Table HALFS.csv", header = F)
names(EF.wilcox.H)=c("Trial4", "Trial5")

##for EF April
wilcox.test(EF.wilcox.H$Trial4,EF.wilcox.H$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)

# Transform into long data: 
# gather the before and after values in the same column
EF.APR.long.H <- EF.wilcox.H %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(EF.APR.long.H, 11)
#wow it worked!!!!
#summary stats?
EF.APR.long.H %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")


bxp.EF.APR.HALFS <- ggpaired(EF.APR.long.H, x = "trial", y = "pM",
                             order = c("Trial4", "Trial5"),
                             ylab = "[pM]", xlab = "April East Fork Trials 1/2*LOD")
bxp.EF.APR.HALFS


#
EF.Apr.diff.H <- EF.wilcox.H %>% mutate(differences = Trial5 - Trial4)
gghistogram(EF.Apr.diff.H, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

EF.APR.matrix.H=as.matrix(EF.Apr.diff.H)

#a large p value allows us to accept NULL which is "symmetric enough"
SYM3 = symmetry.test(EF.APR.matrix.H, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM3$statistic #1.309556 
SYM3$p.value #0.1903461 

stat.test=
  wilcox.test(x=EF.wilcox.H$Trial4, y=EF.wilcox.H$Trial5, paired = TRUE) %>%
  add_significance()
stat.test
#V = 1, p-value = 0.001953

EF.APR.long.H %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)


##LOD itself SUBSTITUTED for BLOD concentrations
##for EF April df read in
EF.wilcox.L=read.csv("EF Apr Trial Table LODS.csv", header = F)
names(EF.wilcox.L)=c("Trial4", "Trial5")

##for EF April
wilcox.test(EF.wilcox.L$Trial4,EF.wilcox.L$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)


# Transform into long data: 
# gather the before and after values in the same column
EF.APR.long.L <- EF.wilcox.L %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(EF.APR.long.L, 11)
#wow it worked!!!!
#summary stats?
EF.APR.long.L %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")

bxp.EF.APR.LODS <- ggpaired(EF.APR.long.L, x = "trial", y = "pM",
                            order = c("Trial4", "Trial5"),
                            ylab = "[pM]", xlab = "April East Fork Trials LOD")
bxp.EF.APR.LODS



EF.Apr.diff.L <- EF.wilcox.L %>% mutate(differences = Trial5 - Trial4)
gghistogram(EF.Apr.diff.L, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

EF.APR.matrix.L=as.matrix(EF.Apr.diff.L)

#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
SYM4 = symmetry.test(EF.APR.matrix.L, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM4$statistic #1.3263
SYM4$p.value #0.1847

stat.test=
  wilcox.test(x=EF.wilcox.L$Trial4, y=EF.wilcox.L$Trial5, paired = TRUE) %>%
  add_significance()
stat.test
#V = 1, p-value = 0.001953

EF.APR.long.L %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)

#############SWITCH TO BEACH ZONE #####################
#############
##ORIGINAL CALCULATED CONCENTRATIONS
##for BZ April df read in
BZ.wilcox=read.csv("BZ APR Trial Table.csv", header = F)
names(BZ.wilcox)=c("Trial4", "Trial5")
BZ.wilcox.matrix=as.matrix(BZ.wilcox)
##for BZ April
wilcox.test(BZ.wilcox$Trial4,BZ.wilcox$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)

#
# Transform into long data: 
# gather the before and after values in the same column
BZ.APR.long <- BZ.wilcox %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(BZ.APR.long, 11)
#wow it worked!!!!
#summary stats?
BZ.APR.long %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")

bxp.BZ.APR.ORIG <- ggpaired(BZ.APR.long, x = "trial", y = "pM",
                            order = c("Trial4", "Trial5"),
                            ylab = "[pM]", xlab = "April Beach Zone Trials")
bxp.BZ.APR.ORIG


#
BZ.Apr.diff <- BZ.wilcox %>% mutate(differences = Trial5 - Trial4)
gghistogram(BZ.Apr.diff, x = "differences", y = "..density..", 
            fill = "lightblue",bins = 5, add_density = TRUE)

BZ.APR.matrix=as.matrix(BZ.Apr.diff)

#a large p value allows us to accept NULL which is "symmetric enough"
SYM5 = symmetry.test(BZ.APR.matrix, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM5$statistic #0.02454369 
SYM5$p.value #0.9804189 

stat.test=
  wilcox.test(x=BZ.wilcox$Trial4, y=BZ.wilcox$Trial5, paired = TRUE) %>%
  add_significance()
stat.test
#V = 4, p-value = 0.006836

BZ.effsize=BZ.APR.long %>%
  wilcox_effsize(pM ~ trial, paired = TRUE,  alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)

print(BZ.effsize)

##ZEROS SUBSTITUTED for BLOD concentrations
##for BZ April df read in
BZ.wilcox.0=read.csv("BZ APR Trial Table ZEROS.csv", header = F)
names(BZ.wilcox.0)=c("Trial4", "Trial5")

##for BZ April
wilcox.test(BZ.wilcox.0$Trial4,BZ.wilcox.0$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)


# Transform into long data: 
# gather the before and after values in the same column
BZ.APR.long.0 <- BZ.wilcox.0 %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(BZ.APR.long.0, 11)  #make sure to change this or wont incl. all values (# means how many rows to incl.)

#wow it worked!!!!
#summary stats?
BZ.APR.long.0 %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")


bxp.BZ.APR.ZEROS <- ggpaired(BZ.APR.long.0, x = "trial", y = "pM",
                             order = c("Trial4", "Trial5"),
                             ylab = "[pM]", xlab = "April Beach Zone Trials ZEROS")
bxp.BZ.APR.ZEROS



BZ.Apr.diff.0 <- BZ.wilcox.0 %>% mutate(differences = Trial5 - Trial4)
gghistogram(BZ.Apr.diff.0, x = "differences", y = "..density..", 
            fill = "lightgreen",bins = 5, add_density = TRUE)

BZ.APR.matrix.0=as.matrix(BZ.Apr.diff.0)

#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
SYM6 = symmetry.test(BZ.APR.matrix.0, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM6$statistic #-0.1709171 
SYM6$p.value #0.8642889 

#fine to use Wilcoxon signed rank in this case
stat.test=
  wilcox.test(x=BZ.wilcox.0$Trial4, y=BZ.wilcox.0$Trial5, paired = TRUE) %>%
  add_significance()
stat.test
#V = 4, p-value = 0.006836

BZ.APR.long.0 %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)


##HALF the LOD SUBSTITUTED for BLOD concentrations
##for BZ April df read in
BZ.wilcox.H=read.csv("BZ APR Trial Table HALFS.csv", header = F)
names(BZ.wilcox.H)=c("Trial4", "Trial5")

##for BZ April
wilcox.test(BZ.wilcox.H$Trial4,BZ.wilcox.H$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)

# Transform into long data: 
# gather the before and after values in the same column
BZ.APR.long.H <- BZ.wilcox.H %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(BZ.APR.long.H, 11)
#wow it worked!!!!
#summary stats?
BZ.APR.long.H %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")
#not working
bxp.BZ.APR.HALFS <- ggpaired(BZ.APR.long.H, x = "trial", y = "pM",
                             order = c("Trial4", "Trial5"),
                             ylab = "[pM]", xlab = "April Beach Zone Trials 1/2*LOD")
bxp.BZ.APR.HALFS


#
BZ.Apr.diff.H <- BZ.wilcox.H %>% mutate(differences = Trial5 - Trial4)
gghistogram(BZ.Apr.diff.H, x = "differences", y = "..density..", 
            fill = "lightblue",bins = 5, add_density = TRUE)

BZ.APR.matrix.H=as.matrix(BZ.Apr.diff.H)


#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
SYM7 = symmetry.test(BZ.APR.matrix.H, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM7$statistic
SYM7$p.value

stat.test=
  wilcox.test(x=BZ.wilcox.H$Trial4, y=BZ.wilcox.H$Trial5, paired = TRUE) %>%
  add_significance()
stat.test

#V = 4, p-value = 0.006836
#
BZ.APR.long.H %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)


##LOD itself SUBSTITUTED for BLOD concentrations
##for BZ April df read in
BZ.wilcox.L=read.csv("BZ APR Trial Table LODS.csv", header = F)
names(BZ.wilcox.L)=c("Trial4", "Trial5")

##for BZ April
wilcox.test(BZ.wilcox.L$Trial4,BZ.wilcox.L$Trial5, 
            paired = T, 
            alternative = "two.sided",
            mu = 0, 
            exact = T, 
            correct = TRUE, 
            conf.int = T,
            conf.level = 0.95)


# Transform into long data: 
# gather the before and after values in the same column
BZ.APR.long.L <- BZ.wilcox.L %>%
  gather(key = "trial", value = "pM", Trial4, Trial5)
head(BZ.APR.long.L, 11)
#wow it worked!!!!
#summary stats?
BZ.APR.long.L %>%
  group_by(trial) %>%
  get_summary_stats(pM, type = "median_iqr")

bxp.BZ.APR.LODS <- ggpaired(BZ.APR.long.L, x = "trial", y = "pM",
                            order = c("Trial4", "Trial5"),
                            ylab = "[pM]", xlab = "April East Fork Trials LOD")
bxp.BZ.APR.LODS

BZ.Apr.diff.L <- BZ.wilcox.L %>% mutate(differences = Trial5 - Trial4)
gghistogram(BZ.Apr.diff.L, x = "differences", y = "..density..", 
            fill = "lightblue",bins = 5, add_density = TRUE)

BZ.APR.matrix.L=as.matrix(BZ.Apr.diff.L)

#use stat test:
#a large p value allows us to accept NULL which is "symmetric enough"
SYM8 = symmetry.test(BZ.APR.matrix.L, boot = F)
#asymptotic distribution of the chosen statistic if you choose FALSE for "boot"
SYM8$statistic #0.3145264 
SYM8$p.value #0.7531212 


stat.test=
  wilcox.test(x=BZ.wilcox.L$Trial4, y=BZ.wilcox.L$Trial5, paired = TRUE) %>%
  add_significance()
stat.test
##V = 7, p-value = 0.01855

BZ.APR.long.L %>%
  wilcox_effsize(pM ~ trial, paired = TRUE, alternative ="two.sided", ci=TRUE, ci.type = "perc",  conf.level =0.95)


