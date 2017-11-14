# Source: Based on example from https://www.stat.ubc.ca/~lang/Stat527a/ex4.r

library(nlme)
library(ggplot2)


#### Example:  NLME models for the AIDS data ####
# We return to the HIV dataset described in Chapter 1. We only consider 
# the first 90 day data. The viral load data "lgcopy" (response) is 
# log10-transformed RNA (to make RNA data more normally distributed).  

dat0 <- read.table("https://www.stat.ubc.ca/~lang/Stat527a/aids.dat2",head=T)

dat1 <- dat0[dat0$day<=90, ]   # use only first 90-day data

dat2 <- dat1[!apply(is.na(dat1),1,any),]  # remove missing data 

dat2 <- merge(dat2, data.frame(patid = unique(dat2$patid), trt = c("A","B"))) # Randomly assign a trt just for demonstration purposes
dat2$lgcopy <- ifelse(dat2$trt == "A", dat2$lgcopy , dat2$lgcopy*0.8) # let's just add some fake effect of treatment (B treatment is 20% lower)

# A nonlinear function: a bi-exponential viral dynamic model ################################

logexp2 <- function(p1,b1,p2,b2,t) log10(exp(p1-b1*t)+exp(p2-b2*t)) 

# Nonlinear model fit, assuming i.i.d. data #################################################

start0 <- c(p1=10,b1=0.5,p2=6,b2=0.005)  # starting value

nls.fit <- nls(lgcopy ~ logexp2(p1,b1,p2,b2,day),
               data =dat2,
               start=start0)

summary(nls.fit) # returns the estimates and significance for the non-linear coefficients 

# plot response curve
plot(dat2$day,dat2$lgcopy)
lines(0:90,predict(nls.fit, newdata = data.frame(day = 0:90)),
      col = "blue")

# Plot residuals 
plot(nls.fit) 


# Next, let's treat the data as longitudinal (or grouped) data  #############################
aids.dat <- groupedData(lgcopy ~ day | patid, data=dat2)

# patid is the patient id
# lgcopy is the log of the numbers is the response variable
# day is the independent variable

plot(aids.dat) # the response for each subject

# A NLME model fit, with random effects on all 4 parameters #################################

start <-  c(10,0.5,6,0.005)# starting value 

nlme.fit <- nlme(lgcopy ~ logexp2(p1,b1,p2,b2,day),
                 fixed = p1 + b1 + p2 + b2 ~ 1, # 'global' parameters estimates
                 random = pdDiag(p1 + b1 + p2 + b2 ~ 1), # we ingnore correlation among parameter estimates 
                 data =aids.dat,
                 start=c(start))

summary(nlme.fit) # the estimates are similar to nls.fit but we see that standard errors of the estimates are much reduced. 

plot(random.effects(nlme.fit))

plot(augPred(nlme.fit,level = c(0,1))) # fixed is the 'average' curve 

# Now we add the effect of trt ###############################################################

start2 <-  rep(c(10,0.5,6,0.005),each = 2)  # starting value 

nlme.fit2 <- nlme(lgcopy ~ logexp2(p1,b1,p2,b2,day),
                 fixed = p1 + b1 + p2 + b2 ~ trt, # parameters are calculated for each treatment
                 random = pdDiag(p1 + b1 + p2 + b2 ~ 1), # we ingnore correlation among parameter estimates 
                 data =aids.dat,
                 start=c(start2))

anova(nlme.fit2) # there is an effect of trt pm p1, b1 and p2, but not b2

# And pooling the effect of trt on b2 #########################################################

start3 <- c(rep(c(10,0.5,6),each = 2),0.005)

nlme.fit3 <- nlme(lgcopy ~ logexp2(p1,b1,p2,b2,day),
                  fixed = list(p1 + b1 + p2 ~ trt,
                               b2 ~ 1), # parameters are calculated for each treatment
                  random = pdDiag(p1 + b1 + p2 + b2 ~ 1), # we ingnore correlation among parameter estimates 
                  data =aids.dat,
                  start=c(start3))

anova(nlme.fit3) # there is an effect of trt pm p1, b1 and p2, but not b2


summary(nlme.fit3)
plot(augPred(nlme.fit3,level = c(0,1)))

# Final plot ###############################################################################

expand.grid(day = 0:90, trt = c("A","B")) -> plotdat
plotdat$pred <- predict(nlme.fit3, newdata = plotdat, level=0)  

ggplot() + 
  geom_point(aes(day,lgcopy, fill = trt),shape = 21, data = dat2) +
  geom_line(aes(day,pred, colour = trt),size = 1.2, data = plotdat) + 
  theme_classic() + theme(legend.position = c(0.8,0.8))