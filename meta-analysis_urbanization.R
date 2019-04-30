##############################################################
# Authors: 
# Alfredo Sanchez-Tojar and Endika Blanco 

# Script first created in Nov 2018

##############################################################
# Description of script and Instructions
##############################################################

# This script is analyze the data for a meta-analysis on the 
# effect of urbanization on biodiversity.

# Endika Blanco-Urdillo, Alfredo Sánchez-Tójar, Juan D. 
# Ibáñez-Álamo. Under review. Methodological approaches 
# to study the effect of urbanization on biodiversity: a 
# review and meta-analysis.

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

# none

#####################
# Importing the data
#####################

meta.final <- read.table("output_final_db/meta_final_reduced.csv",
                         header=T, sep=",")

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

eggers.model <- rma.mv(r, Vr, mods = ~ N.final, random = ~ 1 | studyID/obsID, data=meta.final)
print(eggers.model, digits=3)
confint(model)

# For each meta-regression, one should estimate pseudo R^2 to know how much heterogeneity is explained.
# For that we follow the code from: https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-September/000232.html
# This way of estimating R2 is not very precise, therefore sometimes R2 is negative. In those case, we set R2
# to 0%, and interpret it like that
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
     labels=c("Richness\n","Others\n"),
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

