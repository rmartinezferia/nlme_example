# Rgronomists! ################################################################
# Fitting non-linear mixed effects models using `nlme` package ################
# By R. Martinez-Feria, Nov, 2017 #############################################

# Load package
library(nlme)
library(ggplot2)

# We will be working with the Oats data set included in the nlme package -----
oats <- as.data.frame(nlme::Oats)

# I don't like the units so I convert to kg/ha
oats$nitro <- round(oats$nitro*100*1.12,0)
oats$yield <- oats$yield*32*1.12

# We need to save as groupedData() to help us during the fitting process
oats <- groupedData(yield ~ nitro | Block ,
                    data = oats, 
                    units = list(x = "N rate (kg/ha)",
                                 y = "Yield (kg/ha)"))

plot(oats) # vizualize data 

# A word about mixed-effect models --------------------------------------------

# Fit simple linear model without random effects
fit.lm <- lm(yield~Variety*nitro, data = oats)
anova(fit.lm)
summary(fit.lm)

# Fit a linear mixed-effects model with subplot as a random effect
fit.lme <- lme(yield~Variety*nitro, random = ~1|Block/Variety, data = oats)
anova(fit.lme)
summary(fit.lme)
coefficients(fit.lme)

expand.grid(Block = unique(oats$Block),
            Variety = unique(oats$Variety),
            nitro = 1:70) -> newOats

newOats$predict <- predict(fit.lme, newdata = newOats)
newOats <- merge(newOats, data.frame(oats),all = T)

ggplot(newOats,
       aes(y= predict, x = nitro, colour = Variety)) +
  geom_line() + 
  geom_point(aes(y=yield)) + 
  facet_grid(Block~.)

ggplot(newOats,
       aes(y= predict, x = nitro, colour = Block)) +
  geom_line() + 
  geom_point(aes(y=yield)) + 
  facet_grid(Variety~.)


# Why a non-linear model? ------------------------------------------------------

# Many biological processes are not linear. For intance, can we expect for yield
# to increase linearly with N fertilizer indefinitively?

# Non-linear hypothesis: There is a level where N that no longer increases yield
# Here we use a linear plateau model: 

linearPlateau <- function(x,a,b,c) a + b * (x - c) * (x <= c)

# The nlsList function uses the group data estrcuture to fit a curve for each
# group (i.e. block)

fit.nls <- nlsList(yield ~ linearPlateau(nitro,a,b,c),
                   start = list(a = 30, b = 50, c = 60),
                   data = oats)

summary(fit.nls)
plot(fit.nls)

# No we define our random effects on the model parameters. Here nlme knows the 
# data structure we passed as groupedData() so it knows the random effect is
# the block

fit.nlme <- nlme(fit.nls, random = pdDiag(a + b + c ~ 1))
plot(ranef(fit.nlme))
random.effects(fit.nlme) # look at the scale for b and c
anova(fit.nlme)
coef(fit.nlme)
(fe <- fixed.effects(fit.nlme))

# Now we can ask: Is there an effect of variety? ------------------------------
update(fit.nlme, 
       fixed = list(a + b + c ~ Variety),
       start = c(fe[1],5,5, #Here is where you'll spend most of the time
                 fe[2],0,0,
                 fe[3],0,0)) -> fit.nlme2

anova(fit.nlme2)
coef(fit.nlme2)
plot(augPred(fit.nlme2,level = c(0,1)))

# Final plot ------------------------------------------------------------------
expand.grid(Block = unique(oats$Block),
            Variety = unique(oats$Variety),
            nitro = 1:70) -> newOats

newOats$predict <- predict(fit.nlme, newdata = newOats)
newOats <- merge(newOats, data.frame(oats),all = T)

ggplot(newOats,
       aes(y= predict, x = nitro, colour = Variety)) +
  geom_line() + 
  geom_point(aes(y=yield)) + 
  facet_grid(Block~.)

ggplot(newOats,
       aes(y= predict, x = nitro, colour = Block)) +
  geom_line() + 
  geom_point(aes(y=yield)) + 
  facet_grid(Variety~.)
