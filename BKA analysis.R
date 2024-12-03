############################################################
##BKA analysis 
##written by Andrea Ayala last change 11/26/24
##Reference sites:
##https://www.statmethods.net/management/merging.html
##https://www.statmethods.net/management/subset.html
##https://statsandr.com/blog/descriptive-statistics-in-r/
##https://www.marsja.se/how-to-rename-column-or-columns-in-r-with-dplyr/
##https://stackoverflow.com/questions/7531868/how-to-rename-a-single-column-in-a-data-frame
##http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
##http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
##https://stackoverflow.com/questions/41384075/r-calculate-and-interpret-odds-ratio-in-logistic-regression
##https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
##https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html
##https://stackoverflow.com/questions/1296646/sort-order-data-frame-rows-by-multiple-columns
##http://homepage.stat.uiowa.edu/~luke/classes/STAT4580-2022/proportions.html
##https://statsandr.com/blog/chi-square-test-of-independence-in-r/
##https://bookdown.org/danieljcarter/r4steph/two-way-frequency-tables.html
##http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
##https://stats.stackexchange.com/questions/41270/nonparametric-equivalent-of-ancova-for-continuous-dependent-variables
##https://stackoverflow.com/questions/36459967/ggplot-glm-fitted-curve-without-interaction
##https://danstich.github.io/worst-r/13.1-intro13.html
##https://www.r-bloggers.com/2021/06/oddsplotty-has-landed-on-cran/
##https://stats.stackexchange.com/questions/168928/categorical-variables-in-lasso-regression
#####################################################################################

rm(list=ls()) #this clears the workspace to make sure no leftover variables are left. Not strictly needed but often a good idea
#graphics.off(); #close all graphics windows, in case there are still some open from previous work performed

#Loading necessary packages
library(readr)
library(ggplot2)
library(googlesheets4)
library(httpuv)
library(googledrive)
library(gsheet)
library(dplyr)
library(pastecs)
library(plyr)
library(data.table)
library(tidyr)
library(ggpubr)
library(TTAinterfaceTrendAnalysis)
library(stargazer)
library(RColorBrewer)
library(viridis)
library(glue)
library(hrbrthemes)
library(vcd)
library(ggstatsplot)
library(epiDisplay)
library(tidyverse)
library(magrittr)
library(sm)
library(aod)
library(MASS)
library(leaps)
library(caret)
library(sandwich)
library(broom)
library(boot)
library(pscl)
library(lattice)
library(lmtest)
library(lme4)
library(dplyr)
library(epitools)
library(Epi)
library(pwr)



#Downloading the google sheet from the web - we have only one sheet, so we are going to end up with one dataframe

bka0<-gsheet2tbl('https://docs.google.com/spreadsheets/d/1hxsSShuI08kNSe7LOHiCOdw3k-VE7hEUlzNOR6jOlS4/edit#gid=0')

#Removing columns by name that are not useful to the analysis

bka1 = subset(bka0, select = -c(Tarsus_mm, Serum, Zymo, Culture))

#Imputing lead values below the LOD based on https://github.com/nx10/lnormimp-r/blob/master/docs/example-usage.R
library(fitdistrplus)
library(lnormimp)

#Removing the 999999's which are birds that did not get bled
bka_I <- bka1[-which(bka1$Lead_Imputation=='999999'),] 

#Starting with the data manipulation

#Let's visualize our distribution with a probability density curve - the lognormal can be used based on the distribution

#Now doing this for concentration - first I input zeroes where the NA's were

bka_I$Lead_Imputation[is.na(bka_I$Lead_Imputation)] <- 0

plot(density(bka_I$Lead_Imputation), col = "blue", lty = "dashed", main = "Probability-density")

#Now to censor the data we define lower and upper cutoffs

lower_cutoff <- 2.5
upper_cutoff <- 30

#Which will be visualized by drawing vertical lines on our density curve

abline(v = lower_cutoff, col = "darkgray")
abline(v = upper_cutoff, col = "darkgray")

below_lower <- bka_I$Lead_Imputation < lower_cutoff
above_upper <- bka_I$Lead_Imputation > upper_cutoff

n_below <- sum(below_lower)
n_above <- sum(above_upper)

data_censored <- bka_I$Lead_Imputation[!(below_lower | above_upper)]

#Let's have a quick look at how much data got censored:

data.frame(n = length(bka_I$Lead_Imputation),
           n_below = n_below,
           n_above = n_above,
           n_censored = length(data_censored))

#Now we can visualize the censored distribution:

lines(density(data_censored), col = "red", lty = "dotted")

data_imputed <- lnormimp(
  data_censored,
  censn = c(n_below, n_above),
  cutoff = c(lower_cutoff, upper_cutoff)
)

#Note that if a plausible measurement range is known for real data, 
#it can be specified with the optional range parameter to obtain even more realistic results. The measurement range is set by default to 0-Infinity.

#Let's visualize the results one last time:

lines(density(data_imputed), col = "darkgreen")

#Visually inspecting the data, the imputed distribution (straight green) resembles the original data 
#(dashed blue) much closer than the censored distribution (dotted red).

ks.test(bka_I$Lead_Imputation, data_censored)$statistic
#D 
#0.5

ks.test(bka_I$Lead_Imputation, data_imputed)$statistic
#D 
#0.5

print(data_censored)
print(data_imputed)

#Creating a dataframe of the signal values

# Finding maximum length 
max_ln <- max(c(length(data_censored), length(data_imputed))) 
blood_lead<- data.frame(col1 = c(data_censored,rep(NA, max_ln - length(data_censored))), 
                               col2 = c(data_imputed,rep(NA, max_ln - length(data_imputed)))) 
blood_lead
is.data.frame((blood_lead)) 

#Now taking the mean value of inputted data - and getting rid of all actual lead values

#dropping col1

blood_lead1 = subset(blood_lead, select = -c(col1))

blood_lead2 <- subset(blood_lead1, col2<2.5)

mean(blood_lead2$col2)
# 1.197163 ==== 1.2 is the value I will impute

#Replacing all the zeroes in the imputed dataframe with the mean of 1.2

bka_I$Lead_Imputation[bka_I$Lead_Imputation == 0] <- 1.2

###############################################Analyzing the actual values now, starting with BKA

###Now removing all the values of 999999 for BKA with missing values

bka2 <- subset(bka_I,!(Mean_BKA == 999999))
bka3 <- subset(bka2,!(Mass_g == 999999))


#Now doing a Shapiro-Wilk's tests on the continuous variables

#Starting with mass

shapiro.test(bka3$Mass_g)

#data:  bka3$Mass_kg
#W = 0.96722, p-value = 0.6227 #normal distribution 

library(fitdistrplus) #Identifying the distribution
descdist(bka3$Mass_g)

#Closest to normal distribution

#Next the Bka values

shapiro.test(bka3$Mean_BKA)

#data:  bka3$Mean_BKA
#W = 0.90455, p-value = 0.03141

descdist(bka3$Mean_BKA)
#Between the beta and uniform distribution

shapiro.test(bka3$Lead_Imputation)

#data:  bka3$Lead_Imputation
#W = 0.72952, p-value = 3.466e-05

descdist(bka3$Lead_Imputation)
#Beta distribution

#Now working with the categorical variables
# Counts for each factor
table(bka3$Age)

#AHY  HY   
#18   5    

table(bka3$Sex)

#F  M  
#8 15  


####################################################Bivariate analyses between independent variables
#Lead and weight in grams

corr1 <- cor.test(x=bka3$Lead_Imputation, y=bka3$Mass_g, method = 'spearman')
corr1

#Lead by Age

#Boxplot of Lead by Age

boxplot(Lead_Imputation ~ Age,data=bka3, main="Lead by Age",
        xlab="Age", ylab="Lead")


# PERFORM MANN WHITNEY U TEST
age_result <- wilcox.test(Lead_Imputation ~ Age, data = bka3, exact = FALSE)
age_result

#data:  Lead by Age

#W = 64, p-value = 0.1679
#alternative hypothesis: true location shift is not equal to 0

#Lead by Sex

#Boxplot of Lead by Sex

boxplot(Lead_Imputation ~ Sex,data=bka3, main="Lead by Sex",
        xlab="Sex", ylab="Lead")


# PERFORM MANN WHITNEY U TEST
sex_result <- wilcox.test(Lead_Imputation ~ Sex, data = bka3, exact = FALSE)
sex_result

###data:  Lead_Imputation by Sex
##W = 56, p-value = 0.8109
#alternative hypothesis: true location shift is not equal to 0

# PERFORM Independent samples t-test Mass vs. Age
age_result1 <- t.test (Mass_g ~ Age, var.equal=TRUE, data = bka3)
age_result1


# PERFORM Independent samples t-test Mass vs. Sex
sex_result1 <- t.test (Mass_g ~ Sex, var.equal=TRUE, data = bka3)
sex_result1

#Chi-Square of Age and Sex

table(bka3$Age, bka3$Sex)

chisq.test(bka3$Age, bka3$Sex, correct=FALSE)


#######################################Generalized Linear Model: BKA versus Pb and Age####################################################

bka3$Age <- factor(bka3$Age)

bka3$Sex <- factor(bka3$Sex)

#https://www.r-bloggers.com/2021/10/analysis-of-covariance-ancova-using-r/

#Get summary statistics based on dependent variable and covariate,
library(rstatix)

#summary statistics for dependent variable yield 
bka3 %>% group_by(Age) %>%  get_summary_stats(Mean_BKA, type="common")

#Age   variable     n   min   max median   iqr  mean    sd    se    ci
#<fct> <fct>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 AHY   Mean_BKA    18 0.142  1.14  0.758 0.511 0.768 0.321 0.076 0.16 
#2 HY    Mean_BKA     5 0.268  1.05  0.389 0.372 0.543 0.327 0.146 0.406
 
# summary statistics for covariate height
bka3 %>% group_by(Age) %>%  get_summary_stats(Lead_Imputation, type="common")

#Age   variable            n   min   max median   iqr  mean    sd    se    ci
#<fct> <fct>           <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 AHY   Lead_Imputation    18   1.2 29.5    3.83  6.98  6.50  7.68  1.81  3.82
#2 HY    Lead_Imputation     5   1.2  9.45   1.2   7.35  4.32  4.28  1.92  5.32

#Since it's not parametric, running the glm

glm.bka <- glm(Mean_BKA ~ Lead_Imputation + Age, data = bka3, family = gaussian())
anova(glm.bka)
summary(glm.bka)


library(effectsize)

eta_squared(glm.bka, partial = FALSE)

#Plot of BKA values by Individual

bka3$Bird_ID <- as.character(bka3$Bird_ID) 

ggplot(data = bka3, aes(x = Mean_BKA)) +
  geom_histogram() +
  labs(x ='BKA Capacity', y='Number of Mottled Ducks', title = 'Mottled Duck BKA Distributions')


#Doing a power analysis  #https://www.r-bloggers.com/2021/05/power-analysis-in-statistics-with-r/



#Choosing the power test for the GLM
#For a one-way GLM comparing 2 groups, calculate the sample size needed 
#in each group to obtain a power of 0.80, when the effect size is moderate (0.25) 
#and a significance level of 0.05 is employed.
#https://cran.r-project.org/web/packages/pwr/pwr.pdf
#https://stats.stackexchange.com/questions/523092/numerator-degrees-of-freedom-in-power-analysis-for-regression-r-vs-gpower

u = length(coef(glm.bka))-1

(v <- pwr.f2.test(u = u, v = , f2 = 0.25, sig.level = .05, power = .80)$v)

# Determine required sample size per Age group
ceiling(v + u + 1)
