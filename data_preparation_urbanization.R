##############################################################
# Authors: 
# Alfredo Sanchez-Tojar and Endika Blanco 

# Script first created in Nov 2018

##############################################################
# Description of script and Instructions
##############################################################

# This script is to prepare the dataset needed for a 
# meta-analysis on the effect of urbanization on biodiversity.

# Endika Blanco-Urdillo, Alfredo Sánchez-Tójar, Juan D. 
# Ibáñez-Álamo. Methodological approaches to study the effect 
# of urbanization on biodiversity: a review and meta-analysis.


##############################################################
# Packages needed
##############################################################

library(metafor)
library(plyr)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# function to convert chi-squared to Pearson's r (Nakagawa and Cuthill 2007, Table 2)
chi.to.r<-function(chi,n){
  r<-sqrt(chi/n)
}

# function to convert t values from a continuous predictor variable to Pearson's r 
# (Nakagawa and Cuthill 2007,pg 598)
t.to.r<-function(t,N){
  df<-N-2 #this assumes that there are no other fixed effects included in the models from which t was obtained
  r<-t/sqrt(((t^2)+df))
}

# function to convert F-test to Pearson's r (Lajeunesse 2013)
F.to.r<-function(f,N){
  r<-sqrt(f/(f+N-2))
}

# function to convert R2 to Pearson's r (adjusted) (Nakagawa and Cuthill 2007)
R2.to.r<-function(R2,N,k){
  r<-sqrt(1-(((N-1)*(1-R2))/(N-k-1)))
}

# function to convert Spearman's rho to Pearson's r (Lajeunesse 2013)
spearman_to_pearson <- function(rho){
  r <- 2*sin((pi*rho)/6)
}

# function to convert Kendall's tau to Pearson's r (Lajeunesse 2013)
kendall_to_pearson <- function(tau){
  r <- sin((pi*tau)/2)
}

# function to obtain variance of r (Viechtbauer 2018, Meta-analysis course, Maastrich, Dec 2018)
Vr <- function(r,N){
  Vr<-((1-r^2)^2)/(N-1)
}


#####################
# Preparing the data
#####################

meta <- read.table("meta_AST_v2.csv",header=T,sep=",")

#a bit of information on database history:

#meta_AST differs from meta in that meta_AST has a variable with the prefer estimate
#for the conversion to Zr for each of the studies. Additionally, the signs of some
#correlations R were fixed.

#meta_AST_v2 differs from meta_AST in how es_type is coded, and also some typos have
#been fixed (e.g. Wang et al. 2000 N = 43 instead of 63, Slawski et al 2008 R2 instead 
#of chi-square and F test added too, Kurta & Teramino 1992 N = 77 instead of 42...)

#Note that the database we provide with the manuscript is the final database after
#calculating effect sizes, assigning the correct signs, recoding the necessary variables,
#etc. However, the database used in this script (i.e. meta_AST_v2.csv) can be found in 
# the github repository: https://github.com/ASanchez-Tojar/meta-analysis_urbanization


#####################
# re-coding variables
#####################

#creating a unique identifier for each study
meta$studyID <- as.factor(paste(meta$Author,meta$Year,sep=" "))

#since sample sizes were to small for the original idea, we have recoded the factor
#taxa so that the sample sizes are decent and still biologically meaningful
meta$taxa.2 <- as.factor(ifelse(meta$TAXA=="AMPHIBIANS" |
                                  meta$TAXA=="BIRDS" | 
                                  meta$TAXA=="FISH" | 
                                  meta$TAXA=="MAMMALS",
                                "VERTEBRATES",as.character(meta$TAXA)))

#table(meta$taxa.2)
#additionally, we reorder the levels to set the reference level (intercept) to vertebrates
meta$taxa.2 <- factor(meta$taxa.2,levels(meta$taxa.2)[c(3,1,2)])


#since sample sizes were to small for the original idea, we have recoded the factor
#so that the sample sizes are decent and biologically meaningful
meta$continent.2 <- as.factor(ifelse(meta$CONTINENT=="ASIA" |
                                       meta$CONTINENT=="AUSTRALASIA" | 
                                       meta$CONTINENT=="SOUTH AMERICA",
                                     "OTHERS",as.character(meta$CONTINENT)))

#table(meta$continent.2)
#additionally, we reorder the levels to set the reference level (intercept) to North America
meta$continent.2 <- factor(meta$continent.2,levels(meta$continent.2)[c(2,1,3)])


#since sample sizes were to small for the original idea, we have recoded the factor
#so that the sample sizes are decent and biologically meaningful
# but first, we need to rename SCALE using the right terminology
meta$SCALE <- revalue(meta$SCALE, c("INTERMEDIATE"="LOCAL", "LOCAL"="CITY"))

meta$scale.2 <- as.factor(ifelse(meta$SCALE=="LOCAL" |
                                   meta$SCALE=="REGIONAL",
                                 "LOCAL/REGIONAL",as.character(meta$SCALE)))

#table(meta$scale.2)
#additionally, we reorder the levels to set the reference level (intercept) to Local
meta$scale.2 <- factor(meta$scale.2,levels(meta$scale.2)[c(3,2,1)])


#table(meta$METHOD)
#for method we reorder the levels to set the reference level (intercept) to Urban Gradient
meta$METHOD <- factor(meta$METHOD,levels(meta$METHOD)[c(4,1,3,2)])


#table(meta$index)
#for index we reorder the levels to set the reference level (intercept) to Richness
meta$index <- factor(meta$index,levels(meta$index)[c(2,1)])


#some papers need to be excluded, at least for the time being until we find out if 
#we can actually obtain some data from them.
exclude.for.time.being <- c("Roy et al 1999","?opucki et al 2013",
                            "Trentanovi et al 2013","Faggi & Dadon 2011")

meta <- meta[!(meta$studyID%in%exclude.for.time.being),]


#############################
# re-estimating effect sizes 
#############################

# strategy = we transform all group comparisons (means, sd) to biserial correlations
# and all the rest to Pearson's correlation coefficient. This allows us to compare data
# obtained from comparing urban vs. rural, to studies where urbanization was measured 
# as a continuous variable. Since we use biserial correlations, we do not transform
# Pearson's r to Zr. Note that some biserial correlations might be larger than 1 or
# smaller than -1. To aid interpretation, we will set this values to -1, in which case
# our results should be considered a bit conservative.

# we apply a different procedure for each original estimate
# we subset the database according to the original estimate and then we put
# everything together once the effect sizes are estimated

# For those studies where more than one type of effect size was reported, we chose
# only one based on the following order of priority:
# 1. Pearson's correlation coefficient (r) + N's
# 2. Group means + SD's + N's
# 3. R2 + N's + k (number of predictors included in the model)
# 4. other estimates + N's. The estimates where t, F and chi-squared, and no order of preference was needed to chose among them

# We follow the following rule:
# If the number of species is lower in urban than in rural, the sign of the effect size
# is assigned as (-) negative. In other words, urbanization has a negative effect on biodiversity

# More information in the materials and methods of the study.

# For studies reporting a correlation coefficient, we assumed it to be Pearson's correlation
# coefficient r unless stated otherwise.

# The following are the correlations that according to the methods of their respective papers are not
# Pearson's, and thus, need renaming.

meta$es_type <- as.character(meta$es_type)

meta[meta$studyID=="Orde?ana et al 2010","es_type"] <- "spearmans"
meta[meta$studyID=="Potapova et al 2005","es_type"] <- "kendalls"
meta[meta$studyID=="Kinzig et al 2005" | 
       meta$studyID=="Slawski et al 2008" | 
       meta$studyID=="Zanette et al 2005","es_type"] <- "R2"


# For Riem et al. 2012, we estimated a Spearman's rho using the data of Figure 3, becuase we
# were not entirely sure on what results were presented in the paper 
meta[meta$studyID=="Riem et al 2012","es_type"] <- "spearmans"

meta[meta$studyID=="Riem et al 2012","R"] <- cor(c(4,4,7,7,4,2),c(1,2,3,4,5,6),
                                                 method = c("spearman"))

meta$es_type <- factor(meta$es_type)


#############
# two groups
#############
twogroups <- meta[meta$es_type=="twogroups",]

#we use the escalc() function from metafor to estimate biserial correlations and their sampling variance
twogroups <- escalc(measure="RBIS", m1i=twogroups$UMEAN, m2i=twogroups$RMEAN,
                    sd1i=twogroups$USD,sd2i=twogroups$RSD,
                    n1i=twogroups$UN,n2i=twogroups$RN,data=twogroups,
                    var.names=c("r","Vr"))

#obtaining final N for the effect size as N1+N2
twogroups$N.final <- c(twogroups$RN+twogroups$UN)

#comments: Tzortzakaki et al 2017 R value refers to regression coefficient, so we used the twogroups to
#estimate biserial correlations


##############
# Pearson's r
##############
r <- meta[meta$es_type=="r",]

r$r <- r$R
r$Vr <- Vr(r$r,r$N)
r$N.final <- r$N


#################
# Spearman's rho
#################
spearmans <- meta[meta$es_type=="spearmans",]

spearmans$r <- spearman_to_pearson(spearmans$R)
spearmans$Vr <- Vr(spearmans$r,spearmans$N)
spearmans$N.final <- spearmans$N


#################
# Kendall's tau
#################
kendalls <- meta[meta$es_type=="kendalls",]

kendalls$r <- kendall_to_pearson(kendalls$R)
kendalls$Vr <- Vr(kendalls$r,kendalls$N)
kendalls$N.final <- kendalls$N


#################
# R2
#################

# k is the number of predictors included in the models from which R2 was obtained
# k will be used when estimating r from R2
kvalues <- data.frame(studyID=as.factor(c("Kinzig et al 2005","Slawski et al 2008","Zanette et al 2005")),
                      k=as.integer(c(0,1,0)),
                      stringsAsFactors=FALSE)

R2 <- meta[meta$es_type=="R2",]
R2.2 <- merge(R2,kvalues,by="studyID",all.x=T) # to add each k value to the correct 

R2$r <- R2.to.r(R2.2$R,R2.2$N,R2.2$k)
R2$Vr <- Vr(R2$r,R2$N)
R2$N.final <- R2$N


#############
# t-test
#############
t_test <- meta[meta$es_type=="t-test",]

t_test$r <- t.to.r(t_test$t,t_test$N)
t_test$Vr <- Vr(t_test$r,t_test$N)
t_test$N.final <- t_test$N


#############
# F-test
#############
f_test <- meta[meta$es_type=="F-test",]

f_test$r <- F.to.r(f_test$F,f_test$N)
f_test$Vr <- Vr(f_test$r,f_test$N)
f_test$N.final <- f_test$N


#############
# chi-squared
#############
chi_squared <- meta[meta$es_type=="chi-squared",]

chi_squared$r <- chi.to.r(chi_squared$chi_squared,chi_squared$N)
chi_squared$Vr <- Vr(chi_squared$r,chi_squared$N)
chi_squared$N.final <- chi_squared$N


#################
# FINAL DATABASE
#################

meta.final <- rbind(as.data.frame(twogroups),r,spearmans,kendalls,R2,t_test,f_test,chi_squared)

#################
# revising signs
#################

# We followed the following logic for all effect sizes:
# If the number of species is lower in urban than in rural, the sign of the effect size
# is assigned as (-) negative


# To do so, the following effect sizes need to be negative 
should.be.negative <- c("Kinzig et al 2005","Hill & Wood 2014","Johnson et al 2012",
                        "Kurta & Teramino 1992","Penone et al 2012","Su et al 2011",
                        "Wang et al 2000","Zanette et al 2005","Slawski et al 2008",
                        "Fortel et al 2014","Jokimaki & Suhonen 1993")

# Which is done with the following
meta.final$r <- ifelse(meta.final$studyID%in%should.be.negative,
                       -meta.final$r,meta.final$r)


# are there any r values < -1 or > 1? If so, we will set them to -1 and 1, respectively
# to aid interpretation (more details in the materials and methods of the study). 
meta.final[meta.final$r<(-1),]
meta.final[meta.final$r>1,]


# In our case, there are 5 values < -1, ranging from -1.22 to -1.02, we are going to set
# them to -1 to aid interpretation. Keep in mind that our results should therefore be 
# considered slightly conservative
meta.final[meta.final$r<(-1),"r"] <- -1


# The following additional variable is needed so that the models account for both between-study
# and within-study heterogeneity. 
meta.final$obsID <- 1:nrow(meta.final)

# Saving a reduced version of the final database to make double-checking easier
write.csv(meta.final[,c("studyID","Year","r","Vr","N.final","METHOD","taxa.2",
                        "continent.2","scale.2","obsID")],
          "output_final_db/meta_final_reduced.csv",row.names=FALSE)