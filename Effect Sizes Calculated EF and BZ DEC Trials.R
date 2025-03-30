setwd()

#the purpose of this script is to produce the effect sizes from the Friedman test on Dec HgT concentrations from Trials 1-3 at BZ and EF

#BEACH ZONE SECTION
BZ.effsize = read.csv("BZ.effsize.csv", header = F)
names(BZ.effsize)=c("pM", "trial", "depth")
print(BZ.effsize)

BZ.effsize %>% friedman_effsize(pM ~ trial | depth, ci=T, ci.type = "perc", conf.level = 0.95, nboot = 1000)
#.y.       n effsize conf.low conf.high method    magnitude
 
#pM       11   0.355     0.07      0.77 Kendall W moderate 

############ EAST FORK SECTION  DEC   ##########################

#move on to EF DEC

EF.effsize.DEC = read.csv("EF.effsize.DEC.csv", header = F)
names(EF.effsize.DEC)=c("pM", "trial", "depth")
print(EF.effsize.DEC)

#.y.       n effsize conf.low conf.high method    magnitude
  
#pM       11   0.355      0.1       0.8 Kendall W moderate 

EF.effsize.DEC%>% friedman_effsize(pM ~ trial | depth, ci=T, ci.type = "perc", conf.level = 0.95, nboot = 1000)

