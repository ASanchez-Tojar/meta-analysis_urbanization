##############################################################
# Authors: 
# Alfredo Sanchez-Tojar and Endika Blanco 

# Script first created in Nov 2018

##############################################################
# Description of script and Instructions
##############################################################

# This script is for a meta-analysis on the effect of urbanization
# on biodiversity.


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
#meta_AST differs from meta in that meta_AST has a variable with the prefer estimate
#for the conversion to Zr for each of the studies. Additionally, the signs of some
#correlations R have now been fixed.

#meta_AST_v2 differs from meta_AST in how es_type is coded, and also some typos have
#been fixed (e.g. Wang et al. 2000 N = 43 instead of 63, Slawski et al 2008 R2 instead 
#of chi-square and F test added too, Kurta & Teramino 1992 N = 77 instead of 42...)

#Note that the database we provide with the manuscript is the final database after
#calculating effect sizes, assigning the correct signs, recoding the necessary variables,
#etc. Please go to line 342 of this script to be able to import and analyze that final
#database.


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
#additionally, we reorder the levels to set the reference level (intercept) to Urban Gradient
meta$METHOD <- factor(meta$METHOD,levels(meta$METHOD)[c(4,1,3,2)])


#table(meta$index)
#additionally, we reorder the levels to set the reference level (intercept) to Urban Gradient
meta$index <- factor(meta$index,levels(meta$index)[c(2,1)])

#some papers are excluded, at least for the time being until we find out if we can actually obtain
#some data from them.
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


# If the number of species is lower in urban than in rural, the sign of the effect size
# is assigned as (-) negative. In other words, urbanization has a negative effect on biodiversity

# For studies reporting a correlation coefficient, we assumed it to be Pearson's correlation
# coefficient r unless stated otherwise.

# These are the correlations that according to the methods of their respective papers are not
# Pearson's, and thus, need renaming.

meta$es_type <- as.character(meta$es_type)

meta[meta$studyID=="Orde?ana et al 2010","es_type"] <- "spearmans"
meta[meta$studyID=="Potapova et al 2005","es_type"] <- "kendalls"
meta[meta$studyID=="Kinzig et al 2005" | 
       meta$studyID=="Slawski et al 2008" | 
       meta$studyID=="Zanette et al 2005","es_type"] <- "R2"

# For Riem et al. 2012, we estimated a Spearman's rho using the data of Figure 3, as we
# were not entirely sure on what results were presented in the paper 
meta[meta$studyID=="Riem et al 2012","es_type"] <- "spearmans"

meta[meta$studyID=="Riem et al 2012","R"] <- cor(c(4,4,7,7,4,2),c(1,2,3,4,5,6),
                                                method = c("spearman"))

meta$es_type <- factor(meta$es_type)

#############
# two groups
#############
twogroups <- meta[meta$es_type=="twogroups",]

#we use the escalc() function from metafor to estimate biserial correlations and their variance
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
# to aid interpretation. 
meta.final[meta.final$r<(-1),]
meta.final[meta.final$r>1,]

# In our case, there are 5 values < -1, ranging from -1.22 to -1.02, we are going to set
# them to -1 to aid interpretation. Keep in mind that our results are therefore slightly 
# conservative
meta.final[meta.final$r<(-1),"r"] <- -1


# The following additional variable is needed so that the models account for both between-study
# and within-study heterogeneity. To do so, we include a nested random effect in all models
meta.final$obsID <- 1:nrow(meta.final)

# Saving a reduced version of the final database to make double-checking easier
write.csv(meta.final[,c("studyID","Year","r","Vr","N.final","METHOD","taxa.2",
                        "continent.2","scale.2","obsID")],
          "output_final_db/meta_final_reduced.csv",row.names=FALSE)

# meta.final <- read.table("finaldatabasename.csv",
#                          header=T, sep=",")

################
# META-ANALYSIS
################

# MAIN MODEL = multilevel meta-analysis
# running the model
model <- rma.mv(r, Vr, random = ~ 1 | studyID/obsID, data=meta.final)

# printing the summary results of the model
print(model, digits=3)

# printing the results again, but adding the credibility/prediction interval, which
# uses the heterogeneity to generate a 95% interval that should contain the effect
# sizes of any future or unknown studies 
predict(model, digits=3)

# looking at the foresplot
forest(model)

# estimating the variances for both random effects, and their 95% CI, which we
# use to estimate I2 (see below)
confint(model)

# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(model)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that Potapova et al. 2005 and Wang et al. 2009 are potential outliers. 
# We should consider running a sensitive analysis without these data points (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17 | meta.final$obsID==24 | meta.final$obsID==25,] 
  
# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(model),type="b")

# This plot tell us that Wang et al. 2009 is a very influencial data point. We should
# consider running a sensitive analysis without this data point (perhaps overall, without
# any plant data points)
# meta.final[meta.final$obsID==17,] 


# Estimating heterogeneity as I2 (Nakagawa and Santos 2012). We used the approach 
# offered by Wolfgang Viechtbauer in his webpage:
#http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
W <- diag(1/meta.final$Vr)
X <- model.matrix(model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

# mean I2
I2 <- 100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
round(I2,0) 

# To estimate the 95% CI for the overall I2
bounds <- confint(model) 

# extracting the upper and lower bounds for each variance component
between.lb <- bounds[[1]]$random["sigma^2.1","ci.lb"]
between.ub <- bounds[[1]]$random["sigma^2.1","ci.ub"]
within.lb <- bounds[[2]]$random["sigma^2.2","ci.lb"]
within.ub <- bounds[[2]]$random["sigma^2.2","ci.ub"]

total.lb <- between.lb+within.lb
total.ub <- between.ub+within.ub

I2.lower <- 100 * total.lb / (total.lb + (model$k-model$p)/sum(diag(P)))
round(I2.lower,0) 
I2.upper <- 100 * total.ub / (total.ub + (model$k-model$p)/sum(diag(P)))
round(I2.upper,0) 

# overall I2 and 95% CI = 94.3% (80.3,97.5). That is, high levels of relative heterogeneity

###################
# PUBLICATION BIAS
###################

# since Vr includes r itself (see equation above), Vr and r are correlated, which means 
# that we should not use Vr for the Egger's regression, and instead we are gonna use N
# meta.final$Precision<-sqrt(1/meta.final$Vr)
# meta.final$zresid<-resid*meta.final$Precision

eggers.model <- rma.mv(r, Vr, mods = ~ N.final, random = ~ 1 | studyID/obsID, data=meta.final)
print(eggers.model, digits=3)
confint(model)

# For each meta-regression, one should estimate pseudo R^2 to know how much heterogeneity is explained.
# For that we follow the code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
# This way of estimating R2 is not very precise, therefore sometimes R2 is negative. In those case, we set R2
# to 0%, and intepret it like that
(sum(model$sigma2) - sum(eggers.model$sigma2)) / sum(model$sigma2)

# this indicates that there is some evidence for publication bias but it is not too strong

#########
# PLOT
#########

eggers.plotting <- predict(eggers.model)

newdat <- data.frame(N.final=meta.final$N.final,
                     fit=eggers.plotting$pred,
                     upper=eggers.plotting$ci.ub,
                     lower=eggers.plotting$ci.lb,
                     stringsAsFactors=FALSE)

newdat <- newdat[order(newdat$N.final),]


tiff("plots/Figure_S5_Eggers_regression.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=800)

xaxis <- meta.final$N.final
yaxis <- meta.final$r


op <- par(mar = c(4.5,4.5,1,1)) #bottom, left, top, and right. 

plot(xaxis,yaxis,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,1),
     xlim=c(0,220))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,220,30),
     cex.axis=0.8,tck=-0.02)

axis(2,
     at=round(seq(-1,1,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(xlab = "study sample size", 
      ylab = "effect size (r)",
      line = 2.75, cex.lab=1.4)


lines(newdat$N.final, newdat$fit, lwd=2.75,col="black") 

polygon(c(newdat$N.final,rev(newdat$N.final)),
        c(newdat$lower,rev(newdat$upper)),
        border=NA,col=rgb(0,0,0, 0.25))

points(jitter(xaxis,2),yaxis,
       bg=rgb(0,0,0, 0.35),
       pch=21,
       cex=1)


dev.off()



#################
# time-lag bias
#################

# First one to test for time-lag bias:
# I will first z-transform year, to aid interpretation (this does not change the conclusions)
# more in: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00012.x
meta.final$Year.z <- scale(meta.final$Year)

timelag <- rma.mv(r, Vr, mods = ~ Year.z, random = ~ 1 | studyID/obsID, data=meta.final)
print(timelag, digits=3)
confint(timelag)

# For each meta-regression, one should estimate pseudo R^2 to know how much heterogeneity is explained.
# For that we follow the code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
# This way of estimating R2 is not very precise, therefore sometimes R2 is negative. In those case, we set R2
# to 0%, and intepret it like that
(sum(model$sigma2) - sum(timelag$sigma2)) / sum(model$sigma2)

# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(timelag)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that Potapova et al. 2005 and Wang et al. 2009 are potential outliers. 
# We should consider running a sensitive analysis without these data points (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17 | meta.final$obsID==24 | meta.final$obsID==25,] 

# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(timelag),type="b")

# This plot tell us that Wang et al. 2009 is a very influencial data point. We should
# consider running a sensitive analysis without this data point (perhaps overall, without
# any plant data points)
# meta.final[meta.final$obsID==17,] 


#########
# PLOT
#########

timelag.plotting <- predict(timelag)

newdat <- data.frame(Year=meta.final$Year,
                     fit=timelag.plotting$pred,
                     upper=timelag.plotting$ci.ub,
                     lower=timelag.plotting$ci.lb,
                     stringsAsFactors=FALSE)

newdat <- newdat[order(newdat$Year),]


tiff("plots/Figure_S6_Time_lag_bias.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=800)

xaxis <- meta.final$Year
yaxis <- meta.final$r
cex.study <- meta.final$N.final/40


op <- par(mar = c(4.5,4.5,1,1)) #bottom, left, top, and right. 

plot(xaxis,yaxis,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,1),
     xlim=c(1982,2018))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(1982,2018,6),
     cex.axis=0.8,tck=-0.02)

axis(2,
     at=round(seq(-1,1,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(xlab = "year of publication", 
      ylab = "effect size (r)",
      line = 2.75, cex.lab=1.4)


lines(newdat$Year, newdat$fit, lwd=2.75,col="black") 

polygon(c(newdat$Year,rev(newdat$Year)),
        c(newdat$lower,rev(newdat$upper)),
        border=NA,col=rgb(0,0,0, 0.25))

points(jitter(xaxis,2),yaxis,
       bg=rgb(0,0,0, 0.35),
       pch=21,
       cex=cex.study)


cex.legend <- c(4/40,30/40,100/40)

legend(1982,1,
       c("N =   4",
         "N =  30",
         "N = 100"),
       pt.bg=rgb(0,0,0, 0.75),
       pt.cex=cex.legend,
       pch=21,
       inset=c(0,0),
       y.intersp=1,x.intersp=1)


dev.off()


##################
#META-REGRESSIONs
##################

#########
# method
#########

#Second one, to test if the method influence the effect size
metareg.method <- rma.mv(r, Vr, mods = ~ METHOD, random = ~ 1 | studyID/obsID, data=meta.final)
print(metareg.method, digits=3)
table(meta.final$METHOD)
confint(metareg.method)

#for each meta-regression, one should estimate the R2marginal. It seems that one can do it
#in metafor with the following code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
(sum(model$sigma2) - sum(metareg.method$sigma2)) / sum(model$sigma2)


# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(metareg.method)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that mostly Wang et al. 2009 is potential outlier. 
# We should consider running a sensitive analysis without this data point (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17,] 

# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(metareg.method),type="b")

# This plot tell as that Wang et al. 2009 is a very influencial data point. We should
# consider running a sensitive analysis without this data point (perhaps overall, without
# any plant data points)
# meta.final[meta.final$obsID==17,] 


#########
# PLOT
#########

method.plotting <- predict(metareg.method)

newdat <- data.frame(method=meta.final$METHOD,
                     fit=method.plotting$pred,
                     upper=method.plotting$ci.ub,
                     lower=method.plotting$ci.lb,
                     stringsAsFactors=FALSE)

# as one can see fit upper and lower are the same for all the values
# correspoding to a specific level of the factor, that, we reduce
# the database to have only four rows.
newdat <- newdat[c(2,3,1,21),]


#############
# Actual plot

tiff("plots/Figure_3_metaregression_method.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)


op <- par(mar = c(3,4.5,1,1)) #bottom, left, top, and right. 

plot(c(0,1,2,3),newdat$fit,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,0.5),
     xlim=c(-0.4,3.4))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,3,1),
     labels=c("Urban\ngradient","City\ncomparison","Urban-Rural\ncomparison","Temporal\ngradient"),
     cex.axis=0.7,tck=-0.02)

axis(2,
     at=round(seq(-1,0.5,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(ylab = "effect size (r)",
      line = 2.75, cex.lab=1.75)


arrows(0,newdat$lower[1],
       0,newdat$upper[1],
       angle=90,code=3,
       col=rgb(0,0,139/255,0.85),
       length = 0,lwd=3.5)

arrows(1,newdat$lower[2],
       1,newdat$upper[2],
       angle=90,code=3,
       col=rgb(0,1,1, 0.85),
       length = 0,lwd=3.5)

arrows(2,newdat$lower[3],
       2,newdat$upper[3],
       angle=90,code=3,
       col=rgb(0,104/255,139/255,0.85),
       length = 0,lwd=3.5)

arrows(3,newdat$lower[4],
       3,newdat$upper[4],
       angle=90,code=3,
       col=rgb(1,0,0, 0.85),
       length = 0,lwd=3.5)

points(0,newdat$fit[1],
       col=rgb(0,0,139/255,0.95),
       pch=19,
       cex=2)

points(1,newdat$fit[2],
       col=rgb(0,1,1, 0.95),
       pch=19,
       cex=2)

points(2,newdat$fit[3],
       col=rgb(0,104/255,139/255,0.95),
       pch=19,
       cex=2)

points(3,newdat$fit[4],
       col=rgb(1,0,0, 0.95),
       pch=19,
       cex=2)


text(0.04,0.3,"k = 23",adj=0.5,cex=0.9)
text(1.04,0.3,"k = 13",adj=0.5,cex=0.9)
text(2.04,0.3,"k =  6",adj=0.5,cex=0.9)
text(3.04,0.3,"k =  2",adj=0.5,cex=0.9)


dev.off()



#################
# Taxa
#################

#third one, to test if the taxa influence the effect size
metareg.taxa <- rma.mv(r, Vr, mods = ~ taxa.2, random = ~ 1 | studyID/obsID, data=meta.final)
print(metareg.taxa, digits=3)
confint(metareg.taxa)
table(meta.final$taxa.2)


#for each meta-regression, one should estimate the R2marginal. It seems that one can do it
#in metafor with the following code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
round(((sum(model$sigma2) - sum(metareg.taxa$sigma2)) / sum(model$sigma2))*100,0)


# this funnel plot seems to indicate that there might be negative effects
# missing, which would be the opposite of publication bias :) But this might
# be the consequence of including plants, and having heterogeneity
funnel(metareg.taxa)

# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(metareg.taxa)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that mostly Wang et al. 2009 but also Kinzig et al 2005 are potential outliers. 
# We should consider running a sensitive analysis without these data points (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17 | meta.final$obsID==29,] 

# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(metareg.taxa),type="b")

# This plot tell as that Wang et al. 2009 and mostly Potapova et al 2005 are very influencial 
# data points. We should consider running a sensitive analysis without this data point 
# (perhaps overall, without any plant data points)
# meta.final[meta.final$obsID==17 | meta.final$obsID==25 ,] 


#########
# PLOT
#########

taxa.plotting <- predict(metareg.taxa)

newdat <- data.frame(method=meta.final$taxa.2,
                     fit=taxa.plotting$pred,
                     upper=taxa.plotting$ci.ub,
                     lower=taxa.plotting$ci.lb,
                     stringsAsFactors=FALSE)

# as one can see fit upper and lower are the same for all the values
# correspoding to a specific level of the factor, that, we reduce
# the database to have only four rows.
newdat <- newdat[c(1,5,10),]


#############
# Actual plot

tiff("plots/Figure_6_metaregression_taxa.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)


op <- par(mar = c(3,4.5,1,1)) #bottom, left, top, and right. 

plot(c(0,1,2),newdat$fit,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,0.5),
     xlim=c(-0.4,2.4))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,2,1),
     labels=c("Vertebrates\n","Invertebrates\n","Plants\n"),
     cex.axis=0.7,tck=-0.02)

axis(2,
     at=round(seq(-1,0.5,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(ylab = "effect size (r)",
      line = 2.75, cex.lab=1.75)


arrows(0,newdat$lower[1],
       0,newdat$upper[1],
       angle=90,code=3,
       col=rgb(238/255,197/255,145/255,0.85),
       length = 0,lwd=3.5)

arrows(1,newdat$lower[2],
       1,newdat$upper[2],
       angle=90,code=3,
       col=rgb(139/255,115/255,85/255, 0.85),
       length = 0,lwd=3.5)

arrows(2,newdat$lower[3],
       2,newdat$upper[3],
       angle=90,code=3,
       col=rgb(0,100/255,0,0.85),
       length = 0,lwd=3.5)


points(0,newdat$fit[1],
       col=rgb(238/255,197/255,145/255,0.95),
       pch=19,
       cex=2)

points(1,newdat$fit[2],
       col=rgb(139/255,115/255,85/255, 0.95),
       pch=19,
       cex=2)

points(2,newdat$fit[3],
       col=rgb(0,100/255,0,0.95),
       pch=19,
       cex=2)



text(0.04,0.3,"k = 27",adj=0.5,cex=0.9)
text(1.04,0.3,"k = 10",adj=0.5,cex=0.9)
text(2.04,0.3,"k =  7",adj=0.5,cex=0.9)


dev.off()


#################
# continent
#################

#fourth one, to test if the continent influence the effect size
metareg.continent <- rma.mv(r, Vr, mods = ~ continent.2, random = ~ 1 | studyID/obsID, data=meta.final)
print(metareg.continent, digits=3)
confint(metareg.continent)
table(meta.final$continent.2)


#for each meta-regression, one should estimate the R2marginal. It seems that one can do it
#in metafor with the following code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
round(((sum(model$sigma2) - sum(metareg.continent$sigma2)) / sum(model$sigma2))*100,0)


# this funnel plot seems to indicate that there might be negative effects
# missing, which would be the opposite of publication bias :) But this might
# be the consequence of including plants, and having heterogeneity
funnel(metareg.continent)

# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(metareg.continent)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that again Potapova et al. 2005 and  Wang et al. 2009 are potential outliers. 
# We should consider running a sensitive analysis without these data points (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17 | meta.final$obsID==24 | meta.final$obsID==25,] 

# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(metareg.continent),type="b")

# This plot tell as that gain Wang et al. 2009 and  Potapova et al 2005 are very influencial 
# data points, but also Leveau 2017. 
# We should consider running a sensitive analysis without this data point 
# (perhaps overall, without any plant data points)
# meta.final[meta.final$obsID==13 | meta.final$obsID==17 | meta.final$obsID==26,] 


#########
# PLOT
#########

continent.plotting <- predict(metareg.continent)

newdat <- data.frame(method=meta.final$continent.2,
                     fit=continent.plotting$pred,
                     upper=continent.plotting$ci.ub,
                     lower=continent.plotting$ci.lb,
                     stringsAsFactors=FALSE)

# as one can see fit upper and lower are the same for all the values
# correspoding to a specific level of the factor, that, we reduce
# the database to have only four rows.
newdat <- newdat[c(1,7,11),]


#############
# Actual plot

tiff("plots/Figure_4_metaregression_continent.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)


op <- par(mar = c(3,4.5,1,1)) #bottom, left, top, and right. 

plot(c(0,1,2),newdat$fit,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,0.5),
     xlim=c(-0.4,2.4))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,2,1),
     labels=c("North\nAmerica","Europe\n","Others\n"),
     cex.axis=0.7,tck=-0.02)

axis(2,
     at=round(seq(-1,0.5,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(ylab = "effect size (r)",
      line = 2.75, cex.lab=1.75)


arrows(0,newdat$lower[1],
       0,newdat$upper[1],
       angle=90,code=3,
       col=rgb(139/255,0,0,0.85),
       length = 0,lwd=3.5)

arrows(1,newdat$lower[2],
       1,newdat$upper[2],
       angle=90,code=3,
       col=rgb(0,0,1,0.85),
       length = 0,lwd=3.5)

arrows(2,newdat$lower[3],
       2,newdat$upper[3],
       angle=90,code=3,
       col=rgb(153/255,153/255,153/255,0.85),
       length = 0,lwd=3.5)


points(0,newdat$fit[1],
       col=rgb(139/255,0,0,0.95),
       pch=19,
       cex=2)

points(1,newdat$fit[2],
       col=rgb(0,0,1,0.95),
       pch=19,
       cex=2)

points(2,newdat$fit[3],
       col=rgb(153/255,153/255,153/255,0.95),
       pch=19,
       cex=2)



text(0.04,0.3,"k = 25",adj=0.5,cex=0.9)
text(1.04,0.3,"k = 13",adj=0.5,cex=0.9)
text(2.04,0.3,"k =  6",adj=0.5,cex=0.9)


dev.off()



#################
# scale
#################

#fifth one, to test if the scale influence the effect size
metareg.scale <- rma.mv(r, Vr, mods = ~ scale.2, random = ~ 1 | studyID/obsID, data=meta.final)
print(metareg.scale, digits=3)
confint(metareg.scale)
table(meta.final$scale.2)

#for each meta-regression, one should estimate the R2marginal. It seems that one can do it
#in metafor with the following code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
round(((sum(model$sigma2) - sum(metareg.scale$sigma2)) / sum(model$sigma2))*100,0)


# this funnel plot seems to indicate that there might be negative effects
# missing, which would be the opposite of publication bias :) But this might
# be the consequence of including plants, and having heterogeneity
funnel(metareg.scale)

# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(metareg.scale)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that again Potapova et al. 2005 and  Wang et al. 2009 are potential outliers. 
# We should consider running a sensitive analysis without these data points (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17 | meta.final$obsID==24 | meta.final$obsID==25,] 

# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(metareg.scale),type="b")

# This plot tell as that gain Wang et al. 2009 is a very influencial data point
# We should consider running a sensitive analysis without this data point 
# (perhaps overall, without any plant data points)
# meta.final[meta.final$obsID==17,] 


#########
# PLOT
#########

scale.plotting <- predict(metareg.scale)

newdat <- data.frame(method=meta.final$scale.2,
                     fit=scale.plotting$pred,
                     upper=scale.plotting$ci.ub,
                     lower=scale.plotting$ci.lb,
                     stringsAsFactors=FALSE)

# as one can see fit upper and lower are the same for all the values
# correspoding to a specific level of the factor, that, we reduce
# the database to have only four rows.
newdat <- newdat[c(1,3,6),]


#############
# Actual plot

tiff("plots/Figure_7_metaregression_scale.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)


op <- par(mar = c(3,4.5,1,1)) #bottom, left, top, and right. 

plot(c(0,1,2),newdat$fit,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,0.5),
     xlim=c(-0.4,2.4))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,2,1),
     labels=c("City\n","Local/Regional\n","Global\n"),
     cex.axis=0.7,tck=-0.02)

axis(2,
     at=round(seq(-1,0.5,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(ylab = "effect size (r)",
      line = 2.75, cex.lab=1.75)


arrows(0,newdat$lower[1],
       0,newdat$upper[1],
       angle=90,code=3,
       col=rgb(1,230/255,140/255,0.85),
       length = 0,lwd=3.5)

arrows(1,newdat$lower[2],
       1,newdat$upper[2],
       angle=90,code=3,
       col=rgb(1,165/255,0, 0.85),
       length = 0,lwd=3.5)

arrows(2,newdat$lower[3],
       2,newdat$upper[3],
       angle=90,code=3,
       col=rgb(205/255,55/255,0,0.85),
       length = 0,lwd=3.5)


points(0,newdat$fit[1],
       col=rgb(1,230/255,140/255,0.95),
       pch=19,
       cex=2)

points(1,newdat$fit[2],
       col=rgb(1,165/255,0, 0.95),
       pch=19,
       cex=2)

points(2,newdat$fit[3],
       col=rgb(205/255,55/255,0,0.95),
       pch=19,
       cex=2)



text(0.04,0.3,"k = 25",adj=0.5,cex=0.9)
text(1.04,0.3,"k = 13",adj=0.5,cex=0.9)
text(2.04,0.3,"k =  6",adj=0.5,cex=0.9)


dev.off()



#################
# index
#################

#fifth one, to test if the scale influence the effect size
metareg.index <- rma.mv(r, Vr, mods = ~ index, random = ~ 1 | studyID/obsID, data=meta.final)
print(metareg.index, digits=3)
confint(metareg.index)
table(meta.final$index)

#for each meta-regression, one should estimate the R2marginal. It seems that one can do it
#in metafor with the following code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
round(((sum(model$sigma2) - sum(metareg.index$sigma2)) / sum(model$sigma2))*100,0)


# this funnel plot seems to indicate that there might be negative effects
# missing, which would be the opposite of publication bias :) But this might
# be the consequence of including plants, and having heterogeneity
funnel(metareg.index)

# we can look for potential outliers based on the residuals of the model
# any points below -2 or above 2 could be considered outliers
resid<-rstandard(metareg.index)
plot(meta.final$obsID,resid$z,type="b")
abline(h = 2)
abline(h = -2)

# This plot tell us that again Potapova et al. 2005 and  Wang et al. 2009 are potential outliers. 
# We should consider running a sensitive analysis without these data points (perhaps overall, 
# without any plant data points?)
# meta.final[meta.final$obsID==17 | meta.final$obsID==24 | meta.final$obsID==25,] 

# we can look for potential influential effect sizes by plotting Cook's distances
# and looking for spikes
plot(meta.final$obsID,cooks.distance.rma.mv(metareg.index),type="b")

# This plot tell as that gain Wang et al. 2009 is a very influencial data point
# We should consider running a sensitive analysis without this data point 
# (perhaps overall, without any plant data points)
# meta.final[meta.final$obsID==17,] 


#########
# PLOT
#########

index.plotting <- predict(metareg.index)

newdat <- data.frame(method=meta.final$index,
                     fit=index.plotting$pred,
                     upper=index.plotting$ci.ub,
                     lower=index.plotting$ci.lb,
                     stringsAsFactors=FALSE)

# as one can see fit upper and lower are the same for all the values
# correspoding to a specific level of the factor, that, we reduce
# the database to have only four rows.
newdat <- newdat[c(1,10),]


#############
# Actual plot

tiff("plots/Figure_5_metaregression_index.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)


op <- par(mar = c(3,4.5,1,1)) #bottom, left, top, and right. 

plot(c(0,1),newdat$fit,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,0.5),
     xlim=c(-0.4,1.4))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,1,1),
     labels=c("Richness\n","Diversity\nindex"),
     cex.axis=0.7,tck=-0.02)

axis(2,
     at=round(seq(-1,0.5,0.5),1),
     cex.axis=0.8,las=2,tck=-0.02)

title(ylab = "effect size (r)",
      line = 2.75, cex.lab=1.75)


arrows(0,newdat$lower[1],
       0,newdat$upper[1],
       angle=90,code=3,
       col=rgb(0,0,0,0.85),
       length = 0,lwd=3.5)

arrows(1,newdat$lower[2],
       1,newdat$upper[2],
       angle=90,code=3,
       col=rgb(127/255,127/255,127/255, 0.85),
       length = 0,lwd=3.5)


points(0,newdat$fit[1],
       col=rgb(0,0,0,0.95),
       pch=19,
       cex=2)

points(1,newdat$fit[2],
       col=rgb(127/255,127/255,127/255, 0.95),
       pch=19,
       cex=2)



text(0.04,0.3,"k = 39",adj=0.5,cex=0.9)
text(1.04,0.3,"k =  5",adj=0.5,cex=0.9)


dev.off()

