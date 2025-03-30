setwd()

install.packages("matrixStats")
install.packages("car")
install.packages("datarium")
library(matrixStats)
library(rstatix)
library(datarium)
library(tidyverse)
library(ggpubr)
library(car)
library(dplyr)
library(dunn.test)
#install forcats since using a fancy new way of doing the level specification
library(forcats)

#start with Beach Zone section
##read in 2 df's
pore_data = read.csv("Collected Porewater Results 12.20.23.csv", header = T)
pore_data  #inspect orientation
#specify "levels" for the treatments so it doesn't force alphabetize them which is uneccessarily confusing
#To order the levels in the order they appear in the dataset, use fct_inorder.
##no need to specify levels in order with this command, 
#should use the order in which they are given in dataset!
pore_data$Treatment = fct_inorder(pore_data$Treatment)
pore_data$Treatment

# Convert depth to a factor variable 
pore_data$Depth_cm = as.factor(pore_data$Depth_cm)
pore_data$Depth_cm

pore_data$Sample.ID = as.factor(pore_data$Sample.ID)
pore_data$Sample.ID
#now that depth is a factor, we can specify the order of the levels by hand
##I want them in the reverse of the dataframe, so forcats won't help here
#list by hand with 
pore_data$Depth_cm <- factor(pore_data$Depth_cm, levels = c("20", "10", "2"))
pore_data$Depth_cm 

#example for filtering dataframe:
#georoc[georoc$rock.type == "Basalt" & georoc$tectonic.setting == "Convergent margin", ] 
#need to give full data_frame$column format each time; this creates a new data frame because the part before the [ ] is georoc

#subset the BZ data only 
BZ_pore = pore_data[pore_data$Site =="Beach Zone",]
BZ_pore 
class(BZ_pore$Depth_cm)
BZ_pore$Depth_cm
#create vectors of each BZ treatment group
bz_epa = BZ_pore$pM[(BZ_pore$Depth =="2" & BZ_pore$Treatment =="EPA 1631")]
bz_wb = BZ_pore$pM[(BZ_pore$Depth =="2" & BZ_pore$Treatment =="Heated Water Bath")]
bz_14 = BZ_pore$pM[(BZ_pore$Depth =="2" & BZ_pore$Treatment =="14 Day Digestion")]
bz_24 = BZ_pore$pM[(BZ_pore$Depth =="2" & BZ_pore$Treatment =="24 Day Digestion")]
#create a vector of the BZ treatment groups
bz.mean.object = rbind(bz_epa, bz_wb, bz_14, bz_24)
#calc means of the vectors
rowMeans(bz.mean.object, na.rm =T)
rowSds(bz.mean.object, na.rm =T)
#bz_epa 6.997330 + 0.06094896
#bz_wb 10.209597 + 3.15481313
#bz_14 5.540743 + 0.23591493
#bz_24 9.999591 + 2.88183193
#
#calculate medians for BZ 2 cm groups by row
apply(bz.mean.object, 1, median)
#bz_epa    bz_wb    bz_14    bz_24 
#6.982525 9.747350 5.549600 8.995873 

#a water bath worked no different than 24 day digestion for this matrix
#pM conc. was higher for wb and 24, compared to epa method and 14 day digestion.
#
#move on to BZ 2% 20 cm depth
##create vectors of each BZ treatment group
bz_epa20 = BZ_pore$pM[(BZ_pore$Depth =="20" & BZ_pore$Treatment =="EPA 1631")]
bz_wb20 = BZ_pore$pM[(BZ_pore$Depth =="20" & BZ_pore$Treatment =="Heated Water Bath")]
bz_1420 = BZ_pore$pM[(BZ_pore$Depth =="20" & BZ_pore$Treatment =="14 Day Digestion")]
bz_2420 = BZ_pore$pM[(BZ_pore$Depth =="20" & BZ_pore$Treatment =="24 Day Digestion")]

#create a vector of the BZ treatment groups, 20 cm depth
bz.mean.object20 = rbind(bz_epa20, bz_wb20, bz_1420, bz_2420)
#calc means of the vectors
rowMeans(bz.mean.object20, na.rm =T)
rowSds(bz.mean.object20, na.rm =T)

#bz_epa20 9.326254 +0.3960132 
#bz_wb20 13.574645 +1.9304913
#bz_1420 7.382826 + 0.2613113
#bz_2420 19.457374 +7.0675327

#HISTOGRAMS OF BZ 2 and 20
bz_2 = BZ_pore[(BZ_pore$Depth =="2" ),]
bz_20 = BZ_pore[(BZ_pore$Depth =="20" ),]
bz_2 

hist(bz_2$pM)
hist(bz_20$pM)
#shapiro tests not really normal?
shapiro.test(bz_2$pM)
#W = 0.8325, p-value = 0.007609
shapiro.test(bz_20$pM)
#W = 0.81744, p-value = 0.00621

#not normal, natural log transform
lnbz2=log(bz_2$pM)
lnbz20=log(bz_20$pM)
#create histograms to check distribution after transform
hist(lnbz2)
hist(lnbz20)
#assign the ln-transformed data back to the vectors
bz_2$pM = lnbz2
bz_20$pM = lnbz20

#check normality assumption by treatment group

bz_2 %>%
  group_by(Treatment) %>%
  shapiro_test(pM)
bz_20 %>%
  group_by(Treatment) %>%
  shapiro_test(pM)
#ok, normal enough to use ANOVA


#check variances to see if they are equal before ANOVA
# #pick a stat test for comparing if 2 or more variances are different
# use Levene's test as it is more robust for non normal data than an F test and can assess variances across more than 2 groups
#leveneTest(response variable ~ group variable, data = data, center = median)
#if the p value is greater than 0.05 then accept null hypothesis:
# The variance among the groups is equal.
# https://www.statology.org/levenes-test-r/
levene_test(pM ~ Treatment, data = bz_2, center = median)
#p = 0.04, variances are not equal enough to perform ANOVA
##df1   df2 statistic     p
#3    12      3.83 0.0390
levene_test(pM ~ Treatment, data = bz_20, center = median)
#df1   df2 statistic     p     # ok to use ANOVA for BZ 20 cm samples
#<int> <int>     <dbl> <dbl>
#3    11      1.35 0.308

#convert ID and Treatment to factor variables (did it above at the beginning before subset BZ
class(bz_2$Sample.ID)
class(bz_2$Treatment)
#
#assign not transformed data back to BZ 2 cm, will proceed with non-parametric test
bz_2= BZ_pore[(BZ_pore$Depth =="2" ),]


##Kruskal-Wallis test for BZ 2 cm, , does not meet assumption of equal variances
kruskal.test(bz_2$pM ~ bz_2$Treatment, data = bz_2)
#Kruskal-Wallis chi-squared = 12.728, df = 3, p-value = 0.005263

# Dunn's post-hoc test of multiple comparisons:
BZ.dunn.2 = dunn.test(bz_2$pM, g=bz_2$Treatment, method="bonferroni", kw=TRUE, label=TRUE, 
                       wrap=F, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)


#Comparison of x by group                            
#(Bonferroni)                                  
#Col Mean-|
 # Row Mean |   14 Day D   24 Day D   EPA 1631
#---------+---------------------------------
 # 24 Day D |  
 #            -3.044703
#              0.0070*
  
 # EPA 1631 |  -1.188177   1.856526
#               0.7043     0.1901

 # Heated W |  -2.896181   0.148522  -1.708004
#              0.0113*     1.0000     0.2629

##ANOVA test for BZ 20 cm 
res.aovBZ_20 =  aov(pM ~ Treatment, data = bz_20)
summary(res.aovBZ_20) #0.000196
#tukey hsd
TukeyHSD(res.aovBZ_20, conf.level=.95)
#24 Day Digestion-EPA 1631   0.0013838
## 14 Day Digestion-Heated Water Bath 0.0068966
### 24 Day Digestion-14 Day Digestion  0.0002327

#move on to EF 2 and 10 cm samples
#subset the EF 2% BrCl data only 
EF_pore_2 = pore_data[pore_data$Site =="East Fork" & pore_data$BrCl. == "2",]
EF_pore_2

#further subset EF 2 %  to be 2cm and 10cm
EF_pore_2_2cm = EF_pore_2[EF_pore_2$Depth_cm == "2",]
EF_pore_2_10cm = EF_pore_2[EF_pore_2$Depth_cm == "10",]

EF_pore_2_10cm 

group_by(EF_pore_2_2cm, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(pM, na.rm =T),
    sd =sd(pM, na.rm =T)
  )
#   treatment         count  mean  sd
#1 EPA 1631              4  21.9  6.34
#2 Heated Water Bath     4  64.7 16.6 
#3 14 Day Digestion      4  14.4 13.3 
#4 24 Day Digestion      4  20.5 16.5 
#wow, thanks dplyr made that much faster!

#calculate medians by treatment 
group_by(EF_pore_2_2cm, Treatment) %>%
  summarise(
    count = n(),
    median = median(pM, na.rm =T),
  )
#Treatment             count median

# 1 EPA 1631              4   24.3
# 2 Heated Water Bath     4   58.4
# 3 14 Day Digestion      4   12.1
# 4 24 Day Digestion      4   16.5

group_by(EF_pore_2_10cm, Treatment) %>%
  summarise(
    count = n(),
    median = median(pM, na.rm =T),
    
  )
#Treatment             count median

# 1 EPA 1631              4  0.471
# 2 Heated Water Bath     4 29.4  
# 3 14 Day Digestion      4  0    
# 4 24 Day Digestion      4  0  

group_by(EF_pore_2_10cm, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(pM, na.rm =T),
    
  )
#Treatment             count mean  sd
#1 EPA 1631              4  0.471 0.666
#2 Heated Water Bath     4 29.4   0.654
#3 14 Day Digestion      4  0     0    
#4 24 Day Digestion      4  0     0

#check normality:
hist(EF_pore_2_2cm$pM) #noooo
shapiro.test(EF_pore_2_2cm$pM) #p-value = 0.1128 shapiro.wilk is ok with it but the histogram looks lognormal
hist(EF_pore_2_10cm$pM)
shapiro.test(EF_pore_2_10cm$pM) #definitely not normal W = 0.59113, p-value = 3.385e-05


#ln transforms
ln_ef_2_2 = log(EF_pore_2_2cm$pM)
#assign new ln transformed values back into new column of df
EF_pore_2_2cm$pM_ln = ln_ef_2_2 

ln_ef_2_10 = log(EF_pore_2_10cm$pM)
#assign new ln transformed values back into column of df
EF_pore_2_10cm$pM_ln = ln_ef_2_10 
#check normality again
hist(ln_ef_2_2 ) #better
hist(ln_ef_2_10 ) #no
shapiro.test(ln_ef_2_2 ) #W = 0.6951, p-value = 0.0001491, not normally dist.

shapiro.test(ln_ef_2_10 ) #too many zeros, p value is NA
#assign it back to not transformed, use non-parametric tests for EF 2 and 10 cm with 2% BrCl
EF_pore_2_10cm = EF_pore_2[EF_pore_2$Depth_cm == "10",]
EF_pore_2_10cm$pM
EF_pore_2_2cm = EF_pore_2[EF_pore_2$Depth_cm == "2",]
EF_pore_2_2cm$pM

#use Kruskal-Wallis for EF 2%_2cm as data not become normal with ln-transformation
kruskal.test(pM ~ Treatment, data = EF_pore_2_2cm)
#Kruskal-Wallis chi-squared = 9.1544, df = 3, p-value = 0.02731

#followed by Dunn's post-hoc test
EF.dunn.2cm = dunn.test(EF_pore_2_2cm$pM, g=EF_pore_2_2cm$Treatment, method="bonferroni", kw=TRUE, label=TRUE, 
                                  wrap=F, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)



#Comparison of x by group                            
#(Bonferroni)                                  
#Col Mean-|
 # Row Mean |   14 Day D   24 Day D   EPA 1631
#---------+---------------------------------
# 24 Day D |  -0.297044
#             1.0000

#EPA 1631    -0.816871  -0.519827
            # 1.0000     1.0000

#  Heated W |  -2.747659  -2.450615  -1.930787
#              0.0180*     0.0428     0.1605

#alpha = 0.05
#Reject Ho if p <= alpha/2


#use Kruskal_Wallis for EF 2%_10cm
kruskal.test(pM ~ Treatment, data = EF_pore_2_10cm)
#Kruskal-Wallis chi-squared = 11.642, df = 3, p-value = 0.008717

# Dunn's test of multiple comparisons:
EF.dunn.10 = dunn.test(EF_pore_2_10cm$pM, g=EF_pore_2_10cm$Treatment, method="bonferroni", kw=TRUE, label=TRUE, 
                    wrap=F, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)


#Comparison of x by group                            
#(Bonferroni)                                  
#Col Mean-|
 # Row Mean |   14 Day D   24 Day D   EPA 1631
#---------+---------------------------------
 # 24 Day D |   0.000000
#               1.0000

 # EPA 1631 |  -0.804217  -0.804217
#              1.0000     1.0000

  #Heated W |  -2.954884  -2.954884  -1.608435
#              0.0094*    0.0094*     0.3232

#alpha = 0.05
#Reject Ho if p <= alpha/2

#with Dunn's, HWB is stat diff. than 14 Day digestion (p = 0.0094). HWB is diff. than 24 day digestion (p = 0.0094*).
#HWB is not diff than 1631 (p = 0.32) 24D is not diff from 14 Day, EPA not diff from 14 or 24 day
#HWB was contaminated, all others had insufficient BrCl for complete oxidation
#
#subset the EF 5% data only 
EF_pore_5 = pore_data[ pore_data$Site =="East Fork" & pore_data$BrCl. == "5",]
EF_pore_5

#further subset EF 5 %  to be 2cm and 10cm
EF_pore_5_2cm = EF_pore_5[EF_pore_5$Depth_cm == "2",]
EF_pore_5_10cm = EF_pore_5[EF_pore_5$Depth_cm == "10",]

#summary stats, 2cm, 5% BrCl, calc.means + sd
group_by(EF_pore_5_2cm, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(pM, na.rm =T),
    sd =sd(pM, na.rm =T)
  )
 #  Treatment            Count mean sd
#1 EPA 1631              4  32.1 5.56 
#2 Heated Water Bath     4  34.3 0.327
#3 14 Day Digestion      4  46.6 0.499
#4 24 Day Digestion      4  49.3 1.79 

#calculate medians, 2cm, 5% BrCl
group_by(EF_pore_5_2cm, Treatment) %>%
  summarise(
    count = n(),
    median = median(pM, na.rm =T),
  )
#Treatment            count median

#1 EPA 1631              4   31.8
#2 Heated Water Bath     4   34.2
#3 14 Day Digestion      4   46.4
#4 24 Day Digestion      4   48.6


#summary stats, 10cm, 5% BrCl, calc.means + sd
group_by(EF_pore_5_10cm, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(pM, na.rm =T),
    sd =sd(pM, na.rm =T)
  )
#  Treatment            Count mean sd
#1 EPA 1631              4  32.8 4.33 
#2 Heated Water Bath     4  16.6 2.40 
#3 14 Day Digestion      4  27.5 0.780
#4 24 Day Digestion      4  29.2 2.03 


#calculate medians, 10cm, 5% BrCl
group_by(EF_pore_5_10cm, Treatment) %>%
  summarise(
    count = n(),
    median = median(pM, na.rm =T),
  )
#Treatment             count median

#1 EPA 1631              4   32.8
#2 Heated Water Bath     4   17.2
#3 14 Day Digestion      4   27.0
#4 24 Day Digestion      4   28.7


#check normality

hist(EF_pore_5_2cm$pM) #bimodal?
hist(EF_pore_5_10cm$pM) #okayish
shapiro.test(EF_pore_5_2cm$pM) # good enough W = 0.90404, p-value = 0.1519
shapiro.test(EF_pore_5_10cm$pM)# good enough W = 0.91564, p-value = 0.1652

#now check variances with levene's test
#EF 5% 2 cm
levene_test(pM ~ Treatment, data = EF_pore_5_2cm, center = median)
#       df1   df2 statistic      p
#      3     9      3.97     0.0469   #variances are not equal enough to perform ANOVA

#EF 5% 10 cm
levene_test(pM ~ Treatment, data = EF_pore_5_10cm, center = median)
#df1   df2 statistic       p
#<int> <int>     <dbl>   <dbl>
#  3    11      7.12 0.00631  #variances are not equal enough to perform ANOVA


# use non-paramentric Kruskal-Wallis for both as it does not assume equal variances
kruskal.test(pM ~ Treatment, data = EF_pore_5_2cm)

#Dunn's test
dunn.test(EF_pore_5_2cm$pM, g=EF_pore_5_2cm$Treatment, method="bonferroni", kw=TRUE, label=TRUE, 
          wrap=F, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

#Comparison of x by group                            
#(Bonferroni)                                  
#Col Mean-|
 # Row Mean |   14 Day D   24 Day D   EPA 1631
#---------+---------------------------------
#  24 Day D |  -0.943456
#              1.0000

 # EPA 1631 |   1.933144   2.941742
#               0.1597    0.0098*
 
 # Heated W |   1.257941   2.201398  -0.588348
#              0.6252     0.0831     1.0000


#14 Day Digestion-EPA 1631    p = 0.1597      
#24 Day Digestion-EPA 1631    p =  0.0098*    
#14 Day Digestion-Heated Water Bath p = 0.6252 
#24 Day Digestion-Heated Water Bath p = 0.0831 

#no sig. diff. between HWB and 1631, not sig. diff. btwn. 14 and 24 day digestions
#for this site and depth, 14 and 24 day digestions gave agreement and had higher medians than HWB or 1631.
#longer digestion with 5% was beneficial in supporting oxidation.
#heat and 1631 (RT) performed similarly (p = 1.0), but revealed lower HgT conc. than extended digestions 


#and EF 5% 10 cm
# use non-parametric Kruskal-Wallis  as it does not assume equal variances
kruskal.test(pM ~ Treatment, data = EF_pore_5_10cm)
#Kruskal-Wallis chi-squared = 11.167, df = 3, p-value = 0.01086

#Dunn's test
dunn.test(EF_pore_5_10cm$pM, g=EF_pore_5_10cm$Treatment, method="bonferroni", kw=TRUE, label=TRUE, 
          wrap=F, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

#Comparison of x by group                            
#(Bonferroni)                                  
#Col Mean-|
#  Row Mean |   14 Day D   24 Day D   EPA 1631
#---------+---------------------------------
#  24 Day D |  -0.975900
#               0.9873

#  EPA 1631 |  -1.707825  -0.790569
#                 0.2630     1.0000

#  Heated W |   1.219875   2.371708   3.162277
#                0.6675     0.0531    0.0047*
  

#Heated Water Bath-EPA 1631     p =  0.0047* 
#14 Day Digestion-Heated Water Bath  not significant p = 0.6675 
#24 Day Digestion-Heated Water Bath  p = 0.0531 


#1631, 14 day, and 24 day all performed similarly.  HWB revealed the lowest mean pM value.  
#for this matrix, heating resulted in less Hg revealed even with glass jars.  Not clear if Hg0 evaded
#with 5% BrCl, it seems that incubation from 12 hours to 24 days gave similar results.
