
# Contaminant Exposure Shapes Mitochondrial Bioenergetics in a Wild Seabird----
#Date: 24/03/2025
#Author: Guadalupe Lopez Nava (Max Planck Institute of Biological Intelligence)

## Reading data ----
# Set the working directory
#knitr::opts_knit$set(root.dir = getwd())

data <- read.csv("data/dataset_MitoCont_Shearwaters.csv")

# Display the first few rows of the dataset
head(data)

# Display summary statistics
summary(data)

## Loading the libraries----
library(ggplot2)# Visualizing results
library(lme4) # to run the mixed model 
library(arm) # Bayesian analysis
library(plyr) #calculate means and SD 
library(tidyverse)
library(gridExtra)
library(dplyr)
library(ggtext)
library(ggpubr)
library(cowplot)
library(rstanarm)
library(lubridate)
library(ggtext) 
library(corrplot)
library(ggExtra)
library(patchwork)
library(car) 
library(performance)
library(grid)
library(patchwork)
library(gt)


#Libraries for map Linosa 
library(sf)
library(rnaturalearth)
library(ggspatial)
library(osmdata)

# load Libraries to print tables
library(sjPlot)
library(sjmisc)
library(sjlabelled)

#To import image of shearwater and convert it to raster 
library(ggplotify)
library(magick)



## Setting a theme for all plots----

# Define a custom theme
custom_theme <- theme_bw() +
  theme(text = element_text(family = "Arial"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        title = element_text(size=10),
        axis.title.y = element_markdown())

### Set seed for reproducibility (needed for the simulations)-----
set.seed(123)

## Data curation----

# Convert specific columns to factors
data <- data %>%
  mutate(across(c(year, nest, ring, age, sex, gps, DateMito, TimeMito, Notes), as.factor))


# Ensure the age variable is numeric 
data$age <- as.numeric(as.character(data$age))


# Retaining only two decimals for d13C and d15N
data$d13C <- round(data$d13C, 2)
data$d15N <- round(data$d15N, 2)


# Display the modified data
head(data)

#Scaling the predictors 
data$d15N.Z<-scale(data$d15N)
data$d13C.Z<-scale(data$d13C)
data$age.Z<- scale(data$age)
data$bodymass.Z<-scale(data$bodymass)
data$hg.Z<-scale(data$hg)
data$SUMPFAS.Z<-scale(data$SUMPFAS)


#Statistical Models -----

## Table 1 -----
### Table 1a.  Mercury differences by sex, year and age----

mod <- lmer(hg ~ year + sex +age.Z+ bodymass.Z + (1 | nest), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
check_overdispersion(mod)
MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2

#Extracting the ICC
summary(mod)$varcor  # Extract variance components
icc <- 0.00 / (0.0 + 1.53)  #the residual changes with every model
print(icc)  

#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab1a.html")#Extracting p values from frequentist results

# Save the output:


## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)

### Table 1b. SumPFAS differences by sex, year and age -----

mod <- glmer(SUMPFAS ~ sex + age.Z + bodymass.Z + (1 | nest),
             data = data,       
             family = Gamma(link = "log"),  # Specify Gamma distribution with log link
             na.action = na.exclude)

#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab1b.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)


## Table 2 -----
### Table 2a.  Mercury differences by d15n AND d13C----

mod <- lmer(hg ~ d15N.Z + d13C.Z + (1 | nest), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
check_overdispersion(mod)
MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2

#Extracting the ICC
summary(mod)$varcor  # Extract variance components
icc <- 0.00 / (0.0 + 1.16)  #the residual changes with every model
print(icc)  

#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab2a.html") #Extracting p values from frequentist results

## Assessing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)


### Table 2b.  SUMPFAS differences by d15n AND d13C----

mod <- glm(SUMPFAS ~ d15N.Z + d13C.Z,
           data = data,         
           family = Gamma(link = "log"),  # Specify Gamma distribution with log link
           na.action = na.exclude)


MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2

#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab2b.html")

## Assessing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))


# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

## Table 3 -----
### Table 3a CMR and Mercury (+ predictors)-----

mod <- lmer(CMR ~ hg.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab3a.html")

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 3b OXPHOS and Mercury (+ predictors)-----

mod <- lmer(OXPHOS ~ hg.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab3b.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusions

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)


### Table 3c ETS and Mercury (+ predictors)-----

mod <- lmer(ETS ~ hg.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab3c.html")

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)


### Table 3d LEAK and Mercury (+ predictors)-----

mod <- lmer(LEAK ~ hg.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab3d.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)


### Table 3e FCR1 and Mercury (+ predictors)-----

mod <- lmer(FCR1 ~ hg.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab3e.html")#Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

## Table 4 -----
### Table 4a CMR and SUMPFAS (+ predictors)-----

mod <- lmer(CMR ~ SUMPFAS.Z + sex + age.Z + bodymass.Z + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab4a.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 4b OXPHOS and SUMPFAS (+ predictors)-----

mod <- lmer(OXPHOS ~ SUMPFAS.Z + sex + age.Z + bodymass.Z + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab4b.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)


### Table 4c ETS and SUMPFAS (+ predictors)-----

mod <- lmer(ETS ~ SUMPFAS.Z + sex + age.Z + bodymass.Z + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab4c.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 4d LEAK and SUMPFAS (+ predictors)-----

mod <- lmer(LEAK ~ SUMPFAS.Z + sex + age.Z + bodymass.Z + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab4d.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 4e FCR1 and SUMPFAS (+ predictors)-----

mod <- lmer(FCR1 ~ SUMPFAS.Z + sex + age.Z + bodymass.Z + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab4e.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

## Table 5 -----
### PFOS and mito bioenergetic traits-----
mod <-glm(LEAK~ PFOS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFOS_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

### PFHXS and mito bioenergetic traits-----
mod <-glm(LEAK~ PFHXS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFHXS_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

### PFNA and mito bioenergetic traits-----
mod <-glm(LEAK~ PFNA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFNA_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

### PFDODA and mito bioenergetic traits-----
mod <-glm(LEAK~ PFDODA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFDODA_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

### PFOA and mito bioenergetic traits-----
mod <-glm(LEAK~ PFOA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFOA_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


mod <-glm(FCR1~ PFOA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFOA_FCR1.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


### PFHPS and mito bioenergetic traits-----
mod <-glm(LEAK~ PFHPS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFHPS_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


mod <-glm(FCR1~PFHPS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFHPS_FCR1.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


### PFDA and mito bioenergetic traits-----
mod <-glm(CMR~ PFDA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFDA_CMR.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


mod <-glm(LEAK~PFDA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFDA_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

### PFDS and mito bioenergetic traits-----
mod <-glm(CMR~ PFDS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFDS_CMR.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


mod <-glm(ETS~PFDS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFDS_ETS.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

mod <-glm(LEAK~PFDS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 

r2(mod)
tab_model(mod, file = "tables/Tab5_PFDS_LEAK.html")

#Drawing conclusions 
# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))


## Table 6 -----
### Table 6a CMR and d15N (+ predictors)-----

mod <- lmer(CMR ~ d15N.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab6a.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 6b OXPHOS and d15N (+ predictors)-----

mod <- lmer(OXPHOS ~ d15N.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab6b.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 6c ETS and d15N (+ predictors)-----

mod <- lmer(ETS ~ d15N.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab6c.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 6d LEAK and d15N (+ predictors)-----

mod <- lmer(LEAK ~ d15N.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab6d.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")


qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)

### Table 6e FCR1 and d15N (+ predictors)-----

mod <- lmer(FCR1 ~ d15N.Z + sex + age.Z + bodymass.Z + year + 
              (1 | nest) +
              (1| TimeMito), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

MuMIn::r.squaredGLMM(mod) #Extracting conditional and marginal R2


#Summary of the model 
summary(mod) 
anova(mod)
tab_model(mod, file = "tables/Tab6e.html") #Extracting p values from frequentist results

## Assesing assumptions 

# Tukey-Ascombe plot: 
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")  

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$nest[, 1], main = "Normal Q-Q plot, random effect: nest")
qqline(ranef(mod)$nest[, 1])  # Reference line for QQ plot of random effect

# Normal Q-Q plot of random effect (random_factor) from the model
qqnorm(ranef(mod)$TimeMito[, 1], main = "Normal Q-Q plot, random effect: Timemito")
qqline(ranef(mod)$TimeMito[, 1])  # Reference line for QQ plot of random effect

# Residuals and fitted values
residuals <- resid(mod)
fitted_values <- fitted(mod)

# Plot residuals
plot(fitted_values, residuals)
abline(h = 0, col = "red")

vif(mod) # Extracting multicollinearity values below If vif is higher than 5, then multicolinearity is though to be highly problematic.

## Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Obtain the estimates and their 95 CrI of fixed effects 
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Obtain the estimates and their 95 CrI of random intercepts 
round(quantile(as.vector(bsim@ranef$nest), prob = c(0.025, 0.5, 0.975)),2)
round(quantile(as.vector(bsim@ranef$TimeMito), prob = c(0.025, 0.5, 0.975)),2)



## Figure 1: Map Linosa  -----

# Get the world map with medium resolution as sf object
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define coordinates of Linosa Island
linosa_coords <- data.frame(
  lon = 12.867, 
  lat = 35.867
)

# Convert Linosa coords to sf point
linosa_sf <- st_as_sf(linosa_coords, coords = c("lon", "lat"), crs = 4326)

# Crop map to Mediterranean bounding box (roughly)
med_bbox <- st_bbox(c(xmin = -1, xmax = 28, ymin = 30, ymax = 44), crs = st_crs(4326))

# Plot map zoomed to Mediterranean Sea and add Linosa point
linosa_medit <-ggplot() +
  geom_sf(data = world, fill = "antiquewhite") +
  geom_sf(data = linosa_sf, color = "purple", size = 4, shape = 18) +
  coord_sf(xlim = c(med_bbox["xmin"], med_bbox["xmax"]),
           ylim = c(med_bbox["ymin"], med_bbox["ymax"]),
           expand = FALSE) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering) +
  labs(title = "Linosa Island in the Mediterranean Sea") +
  theme_minimal()

linosa_medit


# Query OpenStreetMap for Linosa Island polygon
linosa_sf <- opq(bbox = "Linosa, Italy") %>%
  add_osm_feature(key = "place", value = "island") %>%
  osmdata_sf()

# Extract polygon geometry for the island
linosa_polygon <- linosa_sf$osm_multipolygons

# And manarazza_coords is a data.frame with longitude and latitude for Manarazza

manarazza_coords <- data.frame(
  lon = 12.8618, # example longitude
  lat = 35.8747  # example latitude
)

# Convert Manarazza coords to sf object
manarazza_sf <- st_as_sf(manarazza_coords, coords = c("lon", "lat"), crs = st_crs(linosa_polygon))

# Plot full island with Manarazza point
linosa <- ggplot() +
  geom_sf(data = linosa_polygon, fill = "antiquewhite") +
  geom_sf(data = manarazza_sf, color = "purple", size = 4, shape=18) +
  coord_sf(xlim = c(12.84, 12.90), ylim = c(35.85, 35.88)) +
  scale_y_continuous(breaks = seq(35.85, 35.88, by = 0.01)) +  # adjust y-axis breaks
  scale_x_continuous(breaks = seq(12.84, 12.90, by = 0.02)) +
  labs(title = "Sampling area on Linosa Island") +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering) +
  theme_minimal()

linosa

#Import a phot of the shearwaters 

# Import image from a file path 
shearimage <- image_read("images/shearwateradult.jpeg")

# Print or plot the image
print(shearimage)

# Convert shearimage (e.g., a raster image) to ggplot
shear_gg <- as.ggplot(~ plot(shearimage))


combined_plot<-ggarrange(
  linosa_medit, linosa, shear_gg,
  ncol = 3,        # All three plots in one row,
  labels = c("A", "B", "C"),
  label.x = 0.05,            # Moves label closer to left edge
  label.y = 0.95,  
  widths = c(1, 1, 0.7),
  font.label = list(size = 20, face = "bold", color = "black"),
  align = "hv" 
)


# Save the combined plot
ggsave(
  filename = "figures/Fig1.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 16,                       # Width in inches
  height = 4
  ,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)

## Figure 2 -----

mod <- lmer(hg ~ year + sex + age + bodymass + (1|nest), 
            data = data, 
            na.action = na.omit)


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~ year + sex + age + bodymass, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))



color_palette <- c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9")

# Calculate the count of observations per sex and year
count_data <- data %>%
  group_by(sex, year) %>%
  summarise(N = n(), .groups = "drop")


# Create the plot
hg_sex_year_plot <- ggplot(data = data, aes(x = sex, y = hg, shape=sex)) +
  geom_jitter(position = position_jitter(width = 0.08), 
              aes(color = factor(sex)), size = 2.7, alpha = 0.3, show.legend = "right") + 
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  facet_wrap(~year) +
  labs(
    x = "Sex",        # Rename x-axis
    y = "<b>Mercury</b> µg/mg"  
  ) +
  scale_color_manual(values = color_palette) + 
  geom_point(data = newdat, aes(x = sex, y = fit, color = sex),
             position = position_dodge(width = 0.5), stroke = 1.5,  size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat, aes(x = sex, y = fit, ymin = lwr, ymax = upr, color = sex),
                width = 0.2, color = "black", position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_text(data = count_data, aes(x = sex, y = min(data$hg, na.rm = TRUE) - 0.5, label = paste("n =", N)), 
            position = position_dodge(width = 1), size = 3, alpha = 0.5, color = "black", fontface = "bold", show.legend = FALSE) +  # Add count labels
  scale_y_continuous(limits = c(0, max(data$hg +1, na.rm = TRUE))) +
  custom_theme 
  


# View the plot
hg_sex_year_plot

## Plotting from the model: age

mod <- lmer(hg ~ year + sex + age + bodymass + (1|nest), 
            data = data, 
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

## Fitting the model

# Calculate the range of values for inddependent variable
range_age <- range(data$age, na.rm = TRUE)
range_age


newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = seq(6, 18, length = 100),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~ year + sex + age + bodymass, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))


newdat<-newdat %>% filter(year=="2020") %>%  droplevels() #leave only the references of the other
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels() #leave only the references of the other
head(newdat)

#Plot the model age 

# Age vs Mercury with credible intervals and no legend

total_count <- nrow(data[!is.na(data$age) & !is.na(data$LEAK), ])

color_palette <- c("#E69F00", "#56B4E9")

# Update plot
p_age_hg <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = age, y = hg, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Age (years)",         # Rename x-axis
    y = "<b>Mercury</b> µg/mg", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "right")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = age, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = age, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = age, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = age, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$hg +2, na.rm = TRUE)))



p_age_hg 


## Plotting the model: d15N


mod <- lmer(hg ~ d15N + d13C + (1 | nest), 
            data = data,
            na.action = na.exclude,  # Exclude missing data only for this model
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

## Fitting the model to plot d15N

# Calculate the range of values for inddependent variable
range_age <- range(data$d15N, na.rm = TRUE)
range_age

newdat <- expand.grid(d15N = seq(7.42, 9.98,length=100),
                      d13C = mean(data$d13C))

Xmat <- model.matrix(~ d15N + d13C, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)


#Plot the model d15N

# Age vs Mercury with credible intervals and no legend

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

# Calculate the total count of observations
total_count <- nrow(data)
# Update plot
p_d15N_hg <- ggplot() + # Start without default data or aes
  geom_point(data = data, aes(x = d15N, y = hg, color = sex, fill = sex, shape= sex), alpha = 0.7, size = 2.7) +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "δ¹⁵N",         # Rename x-axis
    y = "<b>Mercury</b> µg/mg", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape = "Sex"
  )  +
  scale_x_continuous(breaks = seq(min(data$d15N, na.rm = TRUE), max(data$d15N, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) +
  custom_theme +
  theme(legend.position = "none")+
  geom_line(data = newdat, aes(x = d15N, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = d15N, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = d15N, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = d15N, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1)+
  scale_y_continuous(limits = c(0, max(data$hg + 2, na.rm = TRUE)))  

p_d15N_hg 

## Plotting the model: d13C

## Fitting the modelto d13C

# Calculate the range of values for independent variable
range_d13c <- range(data$d13C, na.rm = TRUE)
range_d13c

newdat <- expand.grid(d13C = seq(-19.39, -18.18,length=100),
                      d15N = mean(data$d15N))

Xmat <- model.matrix(~ d15N + d13C, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

#Plot the model d13C

# Age vs Mercury with credible intervals and no legend

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

# Calculate the total count of observations
total_count <- nrow(data)
# Update plot
p_d13C_hg <- ggplot() + # Start without default data or aes
  geom_point(data = data, aes(x = d13C, y = hg, color = sex, fill = sex, shape = sex), alpha = 0.7, size = 2.7) +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "δ¹³C",         # Rename x-axis
    y = "<b>Mercury</b> µg/mg", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape= "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5
  ) +
  geom_line(data = newdat, aes(x = d13C, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = d13C, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = d13C, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = d13C, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_y_continuous(limits = c(0, max(data$hg + 2, na.rm = TRUE)))



p_d13C_hg  


#Plotting both figures. 
# Combine plots with labels
combined_plot <- p_age_hg + hg_sex_year_plot + p_d15N_hg +  p_d13C_hg + 
  plot_annotation(
    tag_levels = "a"  # Automatically labels A, B
  ) +
  theme(
    plot.tag = element_text(face = "bold")  # Make the tags bold
  )

# Save the combined plot
ggsave(
  filename = "figures/Fig2.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 7,                       # Width in inches
  height = 7,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)

## Figure 3 Mercury -----
mod <- lmer(LEAK ~ hg + sex +age + bodymass + year + (1|nest) + (1|TimeMito ), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$hg, na.rm = TRUE)

# Display the range
range_hg

newdat <- expand.grid(hg= seq(1.91, 11.77,length=100),
                      sex =levels(data$sex),
                      age= mean(data$age),
                      bodymass = mean(data$bodymass),
                      year=levels(data$year))

Xmat <- model.matrix(~ hg + sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


## Plotting the model

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

total_count <- nrow(data[!is.na(data$hg) & !is.na(data$LEAK), ])
# Update plot
p_leak_hg <-ggplot() + 
  geom_point(data = data, aes(x = hg, y = LEAK, color = sex, shape= sex), alpha = 0.7, size = 2.7) +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +
  labs(
    x = "Mercury (µg/mg)",        
    y = "<b>LEAK</b> pmol O2sec-1ml-1", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape= "Sex"
  ) +
  custom_theme +
  geom_line(data = newdat, aes(x = hg, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = hg, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = hg, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = hg, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  scale_x_continuous(breaks = seq(min(data$hg, na.rm = TRUE), max(data$hg, na.rm = TRUE), by = 3)) + # Set age ticks by 3
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5
  )

p_leak_hg


## Figure 3 PFAS-----

#PFOS and LEAK

mod <-glm(LEAK~ PFOS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)



### Drawing conclusion

# Number of simulations
nsim <- 1000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim #0.989

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] > 0) / nsim #0.592


# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFOS, na.rm = TRUE)
print(range) # Display the calculated range: 21.89 125.55

# Create a new data frame 'newdat' with 'predictor' ranging from 21.89 to 125.55with 100 points
newdat <- expand.grid(PFOS = seq(21.89 ,125.55, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFOS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

# Create a new data frame 'newdatf' (possibly as a backup of 'newdat')
newdatf <- newdat

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

total_count <- nrow(data[!is.na(data$PFOS) & !is.na(data$LEAK), ])
# Update plot
p1<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFOS, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFOS (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFOS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFOS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = PFOS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFOS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFOS, na.rm = TRUE), max(data$PFOS, na.rm = TRUE), by = 30)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p1


#PFHXS and LEAK


mod <-glm(LEAK~ PFHXS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)



### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFHXS, na.rm = TRUE)
print(range) # Display the calculated range: 21.89 125.55

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFHXS = seq(0.10, 7.95, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFHXS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

# Create a new data frame 'newdatf' (possibly as a backup of 'newdat')
newdatf <- newdat

#Plotting the model 


# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFHXS) & !is.na(data$LEAK), ])
# Update plot
p2<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFHXS, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFHXS (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFHXS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFHXS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFHXS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFHXS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFHXS, na.rm = TRUE), max(data$PFHXS, na.rm = TRUE), by = 2.5)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p2

#PFNA and LEAK

mod <-glm(LEAK~ PFNA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)



### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFNA, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFNA = seq( 0.10, 9.58, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFNA, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 


# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFNA) & !is.na(data$LEAK), ])
# Update plot
p3<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFNA, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFNA (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFNA, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFNA, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFNA, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFNA, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFNA, na.rm = TRUE), max(data$PFNA, na.rm = TRUE), by = 3)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p3


#PFDODA and LEAK


mod <-glm(LEAK~ PFDODA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim #

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim #

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFDODA, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFDODA = seq(1.92, 11.38, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFDODA, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFDODA) & !is.na(data$LEAK), ])
# Update plot
p4<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFDODA, y = LEAK, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFDODA (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFDODA, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFDODA, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFDODA, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFDODA, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFDODA, na.rm = TRUE), max(data$PFDODA, na.rm = TRUE), by = 3)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p4

#PFOA and LEAK

mod <-glm(LEAK~ PFOA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion


# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFOA, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFOA = seq(0.10, 7.43, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFOA, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFOA) & !is.na(data$LEAK), ])
# Update plot
p5<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFOA, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFOA (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFOA, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFOA, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFOA, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFOA, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFOA, na.rm = TRUE), max(data$PFOA, na.rm = TRUE), by = 2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p5


#PFHPS and LEAK


mod <-glm(LEAK~ PFHPS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion


# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFHPS, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFHPS = seq(0.04, 0.48, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFHPS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFHPS) & !is.na(data$LEAK), ])
# Update plot
p6<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFHPS, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFHPS (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFHPS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFHPS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFHPS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFHPS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFHPS, na.rm = TRUE), max(data$PFHPS, na.rm = TRUE), by = 0.20)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p6

#PFDA and LEAK

mod <-glm(LEAK~ PFDA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 
# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 
# Calculate the range of predictor
range <- range(data$PFDA, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFDA = seq(0.02, 6.23, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFDA, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

total_count <- nrow(data[!is.na(data$PFDA) & !is.na(data$LEAK), ])
# Update plot
p7<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFDA, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFDA (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFDA, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFDA, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFDA, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFDA, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFDA, na.rm = TRUE), max(data$PFDA, na.rm = TRUE), by = 2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p7

#PFDS and LEAK

mod <-glm(LEAK~ PFDS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFDS, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFDS = seq(0.02, 0.67, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFDS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'


#Plotting the model 


# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFDS) & !is.na(data$LEAK), ])
# Update plot
p8<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFDS, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFDS (ng/g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFDS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFDS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFDS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFDS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFDS, na.rm = TRUE), max(data$PFDS, na.rm = TRUE), by = 0.2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))

p8


#Plotting Fig 3. 


combined_plot <- wrap_plots(
  list(
    p_leak_hg, plot_spacer(), plot_spacer(), plot_spacer(),
    p1, p2, p3, p4,
    p5, p6, p7, p8
  ),
  design = "
abcd
efgh
ijkl
") +
  plot_annotation(tag_levels = "a", 
                  theme = theme(plot.tag = element_text(face = "bold")))

print(combined_plot)



# Save the combined plot
ggsave(
  filename = "figures/Fig3.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 10,                       # Width in inches
  height = 7,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)

## Figure 4 -----

#PFOA and FCR1

mod <-glm(FCR1~ PFOA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFOA, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFOA = seq(0.10, 7.43, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFOA, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'


#Plotting the model 


# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFOA) & !is.na(data$FCR1), ])
# Update plot
p1<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFOA, y = FCR1, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFOA (ng/g)",         # Rename x-axis
    y = "<b>FCR1</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFOA, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFOA, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFOA, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFOA, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFOA, na.rm = TRUE), max(data$PFOA, na.rm = TRUE), by = 2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$FCR1, na.rm = TRUE)))

p1

#PFHPS and FCR1

mod <-glm(FCR1~ PFHPS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFHPS, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFHPS = seq(0.04, 0.48, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFHPS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFHPS) & !is.na(data$FCR1), ])
# Update plot
p2<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFHPS, y = FCR1, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFHPS (ng/g)",         # Rename x-axis
    y = "<b>FCR1</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape= "Sex"
  ) +
  custom_theme +
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFHPS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFHPS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFHPS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFHPS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFHPS, na.rm = TRUE), max(data$PFHPS, na.rm = TRUE), by = 0.20)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$FCR1, na.rm = TRUE)))

p2

#Plotting Fig 4.

# Combine plots with labels
combined_plot <- p1 + p2 + 
  plot_annotation(
    tag_levels = "a"  # Automatically labels A, B
  ) +
  theme(
    plot.tag = element_text(face = "bold")  # Make the tags bold
  )
# Ensure three columns
combined_plot <- combined_plot + plot_layout(ncol = 2)

# Save the combined plot
ggsave(
  filename = "figures/Fig4.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 6,                       # Width in inches
  height = 3,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)


## Figure 5 -----
#PFDA and CMR

mod <-glm(CMR~ PFDA,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFDA, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFDA = seq(0.02, 6.23, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFDA, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFDA) & !is.na(data$CMR), ])
# Update plot
p1<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFDA, y = CMR, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFDA (ng/g)",         # Rename x-axis
    y = "<b>CMR</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFDA, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFDA, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFDA, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFDA, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFDA, na.rm = TRUE), max(data$PFDA, na.rm = TRUE), by = 2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$CMR, na.rm = TRUE)))

p1

#PFDS and LEAK


mod <-glm(CMR~ PFDS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 

# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFDS, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFDS = seq(0.02, 0.67, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFDS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFDS) & !is.na(data$CRM), ])
# Update plot
p2<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFDS, y = CMR, color = sex, fill = sex, shape =sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFDS (ng/g)",         # Rename x-axis
    y = "<b>CMR</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFDS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFDS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFDS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFDS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFDS, na.rm = TRUE), max(data$PFDS, na.rm = TRUE), by = 0.2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$CMR, na.rm = TRUE)+100))

p2

#PFDS and ETS

mod <-glm(ETS~ PFDS,
          data = data, 
          family = Gamma(link = "log"), 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod)


### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values greater than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


# Fitting the model 

# Calculate the range of predictor
range <- range(data$PFDS, na.rm = TRUE)
print(range) # Display the calculated range: 

# Create a new data frame 'newdat' with 'predictor' ranging from with 100 points
newdat <- expand.grid(PFDS = seq(0.02, 0.67, length = 100))
str(newdat) # Display the structure of the 'newdat' data frame

# Generate the model matrix 'Xmat' based on 'scaled_mass_index'
Xmat <- model.matrix(~PFDS, data = newdat)

# Create an empty matrix 'fitmat' to store the model fit for each simulation
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim@coef[i,])

# Calculate lower and upper quantiles for the model fits
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

# Calculate the fitted values using the fixed effects from 'modfb'
newdat$fit <- exp(Xmat %*% coef(mod))

# Display the content of 'newdat'
newdat
head(newdat) # Display the first few rows of 'newdat'

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$PFDS) & !is.na(data$ETS), ])
# Update plot
p3<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = PFDS, y = ETS, color = sex, fill = sex, shape=sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "PFDS (ng/g)",         # Rename x-axis
    y = "<b>ETS</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape= "Sex"
  ) +
  custom_theme +
  theme(legend.position = "right")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = PFDS, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = PFDS, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x =PFDS, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = PFDS, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$PFDS, na.rm = TRUE), max(data$PFDS, na.rm = TRUE), by = 0.2)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$ETS, na.rm = TRUE)+100))

p3

#Plotting Fig 5.

# Combine plots with labels
combined_plot <- p1 + p2 + p3+
  plot_annotation(
    tag_levels = "a"  # Automatically labels A, B
  ) +
  theme(
    plot.tag = element_text(face = "bold")  # Make the tags bold
  )
# Ensure three columns
combined_plot <- combined_plot + plot_layout(ncol = 3)

# Save the combined plot
ggsave(
  filename = "figures/Fig5.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 10,                       # Width in inches
  height = 3,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)

## Figure 6 -----


## Fitting the model

mod <- lmer(LEAK ~ sex+ d15N+ bodymass + age + year +  (1|nest) + (1|TimeMito ), 
            data = data)

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Calculate the range of values for inddependent variable
range_d15N <- range(data$d15N, na.rm = TRUE)

# Display the range
range_d15N


newdat <- expand.grid(sex=levels(data$sex),
                      d15N= seq(7.42, 9.98,length=100),
                      bodymass = mean(data$bodymass),
                      age= mean(data$age),
                      year=levels(data$year))


Xmat <- model.matrix(~sex+ d15N  + bodymass + age + year, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels() #leave only the references of the other
head(newdat)

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

# Calculate the total count of observations
total_count <- nrow(data[!is.na(data$d15N) & !is.na(data$LEAK), ])

# Update plot
p_d15N_LEAK <- ggplot() + # Start without default data or aes
  geom_point(data = data, aes(x = d15N, y = LEAK, color = sex, fill = sex, shape = sex), alpha = 0.7, size = 2.7) +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "δ¹⁵N",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape = "Sex"
  )  +
  scale_x_continuous(breaks = seq(min(data$d15N, na.rm = TRUE), max(data$d15N, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) +
  custom_theme +
  theme(legend.position = "none")+
  geom_line(data = newdat, aes(x = d15N, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = d15N, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = d15N, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = d15N, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1)+
  scale_y_continuous(limits = c(0, max(data$LEAK + 2, na.rm = TRUE)))  

p_d15N_LEAK



mod <- lmer(FCR1 ~ sex+ d15N+ bodymass + age + year +  (1|nest) + (1|TimeMito ), 
            data = data)

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Calculate the range of values for inddependent variable
range_d15N <- range(data$d15N, na.rm = TRUE)

# Display the range
range_d15N


newdat <- expand.grid(sex=levels(data$sex),
                      d15N= seq(7.42, 9.98,length=100),
                      bodymass = mean(data$bodymass),
                      age= mean(data$age),
                      year=levels(data$year))


Xmat <- model.matrix(~sex+ d15N  + bodymass + age + year, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels() #leave only the references of the other
head(newdat)

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

# Calculate the total count of observations
total_count <- nrow(data[!is.na(data$d15N) & !is.na(data$FCR1), ])

# Update plot
p_d15N_FCR1 <- ggplot() + # Start without default data or aes
  geom_point(data = data, aes(x = d15N, y = FCR1, color = sex, fill = sex, shape = sex), alpha = 0.7, size = 2.7) +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "δ¹⁵N",         # Rename x-axis
    y = "<b>FCR1</b> pmol O2sec-1ml-1", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape = "Sex"
  )  +
  scale_x_continuous(breaks = seq(min(data$d15N, na.rm = TRUE), max(data$d15N, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) +
  custom_theme +
  theme(legend.position = "right")+
  geom_line(data = newdat, aes(x = d15N, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = d15N, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = d15N, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = d15N, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1)+
  scale_y_continuous(limits = c(0, max(data$FCR1 +0.1, na.rm = TRUE)))  

p_d15N_FCR1

# Combine plots with labels
combined_plot <- p_d15N_LEAK + p_d15N_FCR1 + 
  plot_annotation(
    tag_levels = "a"  # Automatically labels A, B
  ) +
  theme(
    plot.tag = element_text(face = "bold")  # Make the tags bold
  )
# Ensure two columns
combined_plot <- combined_plot + plot_layout(ncol = 2)

# Save the combined plot
ggsave(
  filename = "figures/Fig6.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 6,                       # Width in inches
  height = 3,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)

# Supporting information-----
 
## Table S1----

# Function to calculate summary statistics
summarize_data <- function(data, variable) {
  data %>%
    group_by(sex, year) %>%
    summarise(
      n = n(),
      median = round(median(!!sym(variable), na.rm = TRUE), 2),
      mean_sd = paste0(round(mean(!!sym(variable), na.rm = TRUE), 2), 
                       " ± ", 
                       round(sd(!!sym(variable), na.rm = TRUE), 2)),
      min = round(min(!!sym(variable), na.rm = TRUE),2),
      max = round(max(!!sym(variable), na.rm = TRUE),2),
      .groups = "drop"
    ) %>%
    arrange(sex, year)
}

# List of variables to summarize
variables <- c("hg", "d15N", "d13C", "SUMPFAS")

# Apply function to each variable and store results
results_list <- lapply(variables, function(var) {
  summary_table <- summarize_data(data, var)
  summary_table$Variable <- var  # Add variable name
  return(summary_table)
})

# Combine all results into one table
final_summary <- do.call(rbind, results_list)

# Print the summary table
print(final_summary)

# Save table for publication 

gt_table <- final_summary %>%
  gt() %>%  
  tab_header(title = "TableS1: Summary Statistics of Variables by Sex and Year") %>% 
  opt_table_font(font = list(google_font("Arial"))) %>%
  opt_row_striping() %>%
  tab_options(table.width = pct(80))  %>% # narrower table
  tab_style(
    style = cell_text(weight = "bold", align = "left"),
    locations = cells_column_labels(everything())
  ) %>% 
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  )

# Save as HTML
gtsave(gt_table, filename = "suppinfo/TabS1.html")

## Table S2 (additional analyses)----------
###PFOS and bodymass------

mod <-glm(bodymass~ PFOS *sex,
          data = data, 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFOS.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


###PFHXS------

mod <-glm(bodymass~ PFHXS *sex,
          data = data, 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFHXS.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


###PFNA------

mod <-glm(bodymass~ PFNA *sex,
          data = data, 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFNA.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 



###PFDODA------

mod <-glm(bodymass~ PFDODA *sex,
          data = data, 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFDODA.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


###PFOA------

mod <-glm(bodymass~ PFOA *sex,
          data = data, 
          na.action = na.exclude)

summary(mod)
round(coef(mod),3) 
anova(mod)
1 - mod$deviance/mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFOA.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot") # Title for the plot indicating strong shrinkage

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)


### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim) # Display structure of the simulated data

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 




###PFHPS ------
mod <- glm(bodymass ~ PFHPS *sex,
           data = data, 
           na.action = na.exclude)

summary(mod)
round(coef(mod), 3) 
anova(mod)
1 - mod$deviance / mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFHPS.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim)

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 




#### PFDA ------

mod <- glm(bodymass ~ PFDA *sex,
           data = data, 
           na.action = na.exclude)

summary(mod)
round(coef(mod), 3) 
anova(mod)
1 - mod$deviance / mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFDA.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim)

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 


#### PFDS ------

mod <- glm(bodymass ~ PFDS *sex,
           data = data, 
           na.action = na.exclude)

summary(mod)
round(coef(mod), 3) 
anova(mod)
1 - mod$deviance / mod$null.deviance

r2(mod)
tab_model(mod, file = "suppinfo/TabS2_PFDS.html")

### Assessing model assumptions 

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2))

# Tukey-Ascombe plot: Residuals vs Fitted values with a reference line at 0
scatter.smooth(fitted(mod), resid(mod))
abline(h = 0, lty = 2)
title("Tukey-Ascombe plot")

# Normal Q-Q plot of residuals
qqnorm(resid(mod), main = "Normal Q-Q plot, residuals")
qqline(resid(mod))

# Residuals' square root vs Fitted values
scatter.smooth(fitted(mod), sqrt(abs(resid(mod)))) 

plot(mod)

### Drawing conclusion

# Number of simulations
nsim <- 10000

# Simulate from the model 'modfb'
bsim <- sim(mod, n.sim = nsim)
str(bsim)

# Use exp to transform effects from log-link function to natural values
apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Proportion of values greater than 0 in the first column of bsim@fixef
sum(bsim@coef[, 1] > 0) / nsim 

# Proportion of values lower than 0 in the second column of bsim@fixef
sum(bsim@coef[, 2] < 0) / nsim 




## Figure S1-----

##CMR

mod <- lmer(CMR ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

## Plotting the model

 color_palette <- c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9")

# Calculate the count of observations per sex and year
count_data <- data %>%
  group_by(sex, year) %>%
  summarise(N = n(), .groups = "drop")


# Create the plot
p1 <- ggplot(data = data, aes(x = sex, y = CMR, shape= sex, shape = sex)) +
  geom_jitter(position = position_jitter(width = 0.08), aes(color = factor(sex)), size = 2.7, alpha = 0.3, show.legend = "right") + 
  facet_wrap(~year) +
  labs(
    x = "Sex",        # Rename x-axis
    y = "<b>CMR</b> pmol O2sec-1ml-1"  
  ) +
  scale_color_manual(values = color_palette) + 
  scale_x_discrete(labels = c("F" = "Females", "M" = "Males")) +
  geom_point(data = newdat, aes(x = sex, y = fit, color = sex),
             position = position_dodge(width = 0.5), stroke = 1.5,  size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat, aes(x = sex, y = fit, ymin = lwr, ymax = upr, color = sex),
                width = 0.2, color = "black", position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_text(data = count_data, aes(x = sex, y = min(data$CMR, na.rm = TRUE) - 80, label = paste("n =", N)), 
            position = position_dodge(width = 1), size = 3, alpha = 0.5, color = "black", fontface = "bold", show.legend = FALSE) +  # Add count labels
  scale_y_continuous(limits = c(0, max(data$CMR, na.rm = TRUE))) +
  custom_theme 

# View the plot
p1


## Fitting the model: AGE

mod <- lmer(CMR ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = seq(6, 18, length = 100),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$age) & !is.na(data$CMR), ])
# Update plot
p2 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = age, y = CMR, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Age (years)",         # Rename x-axis
    y = "<b>CMR</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = age, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = age, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = age, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = age, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$CMR, na.rm = TRUE)))

p2

## Fitting the model: bodymass

mod <- lmer(CMR ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$bodymass, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = seq(445, 772, length = 100) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

total_count <- nrow(data[!is.na(data$bodymass) & !is.na(data$CMR), ])
# Update plot
p3 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = bodymass, y = CMR, color = sex, fill = sex, shape = sex ), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Body mass (g)",         # Rename x-axis
    y = "<b>CMR</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = bodymass, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = bodymass, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = bodymass, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = bodymass, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$bodymass, na.rm = TRUE), max(data$bodymass, na.rm = TRUE), by = 50)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$CMR, na.rm = TRUE)))



p3

# 1.1 - OXPHOS

mod <- lmer(OXPHOS ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
#newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
#newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

## Plotting the model


 color_palette <- c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9")

# Calculate the count of observations per sex and year
count_data <- data %>%
  group_by(sex, year) %>%
  summarise(N = n(), .groups = "drop")


# Create the plot
p4 <- ggplot(data = data, aes(x = sex, y = OXPHOS, shape= sex)) +
  geom_jitter(position = position_jitter(width = 0.08), aes(color = factor(sex)), size = 2.7, alpha = 0.3, show.legend = "right") + 
  facet_wrap(~year) +
  labs(
    x = "Sex",        # Rename x-axis
    y = "<b>OXPHOS</b> pmol O2sec-1ml-1"  
  ) +
  scale_color_manual(values = color_palette) + 
  scale_x_discrete(labels = c("F" = "Females", "M" = "Males")) +
  geom_point(data = newdat, aes(x = sex, y = fit, color = sex),
             position = position_dodge(width = 0.5), stroke = 1.5,  size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat, aes(x = sex, y = fit, ymin = lwr, ymax = upr, color = sex),
                width = 0.2, color = "black", position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_text(data = count_data, aes(x = sex, y = min(data$CMR, na.rm = TRUE) - 80, label = paste("n =", N)), 
            position = position_dodge(width = 1), size = 3, alpha = 0.5, color = "black", fontface = "bold", show.legend = FALSE) +  # Add count labels
  scale_y_continuous(limits = c(0, max(data$CMR, na.rm = TRUE))) +
  custom_theme 

# View the plot
p4


## Fitting the model: AGE

mod <- lmer(OXPHOS ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = seq(6, 18, length = 100),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$age) & !is.na(data$OXPHOS), ])
# Update plot
p5 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = age, y = OXPHOS, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Age (years)",         # Rename x-axis
    y = "<b>OXPHOS</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = age, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = age, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = age, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = age, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$OXPHOS, na.rm = TRUE)))



p5

## Fitting the model: bodymass

mod <- lmer(OXPHOS ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$bodymass, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = seq(445, 772, length = 100) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$bodymass) & !is.na(data$OXPHOS), ])
# Update plot
p6 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = bodymass, y = OXPHOS, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Body mass (g)",         # Rename x-axis
    y = "<b>OXPHOS</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = bodymass, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = bodymass, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = bodymass, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = bodymass, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$bodymass, na.rm = TRUE), max(data$bodymass, na.rm = TRUE), by = 50)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$OXPHOS, na.rm = TRUE)))



p6

# 1.1 - ETS

mod <- lmer(ETS~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
#newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
#newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


## Plotting the model

 color_palette <- c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9")

# Calculate the count of observations per sex and year
count_data <- data %>%
  group_by(sex, year) %>%
  summarise(N = n(), .groups = "drop")


# Create the plot
p7 <- ggplot(data = data, aes(x = sex, y = ETS, shape = sex)) +
  geom_jitter(position = position_jitter(width = 0.08), aes(color = factor(sex)), size = 2.7, alpha = 0.3, show.legend = "right") + 
  facet_wrap(~year) +
  labs(
    x = "Sex",        # Rename x-axis
    y = "<b>ETS</b> pmol O2sec-1ml-1"  
  ) +
  scale_color_manual(values = color_palette) + 
  scale_x_discrete(labels = c("F" = "Females", "M" = "Males")) +
  geom_point(data = newdat, aes(x = sex, y = fit, color = sex),
             position = position_dodge(width = 0.5), stroke = 1.5,  size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat, aes(x = sex, y = fit, ymin = lwr, ymax = upr, color = sex),
                width = 0.2, color = "black", position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_text(data = count_data, aes(x = sex, y = min(data$ETS, na.rm = TRUE) - 5, label = paste("n =", N)), 
            position = position_dodge(width = 1), size = 3, alpha = 0.5, color = "black", fontface = "bold", show.legend = FALSE) +  # Add count labels
  scale_y_continuous(limits = c(0, max(data$ETS, na.rm = TRUE))) +
  custom_theme 

# View the plot
p7


## Fitting the model: AGE
mod <- lmer(ETS ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = seq(6, 18, length = 100),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$age) & !is.na(data$ETS), ])
# Update plot
p8 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = age, y = ETS, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Age (years)",         # Rename x-axis
    y = "<b>ETS</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = age, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = age, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = age, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = age, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$ETS, na.rm = TRUE)))

p8

## Fitting the model: bodymass

mod <- lmer(ETS ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$bodymass, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = seq(445, 772, length = 100) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$bodymass) & !is.na(data$ETS), ])
# Update plot
p9 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = bodymass, y = ETS, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Body mass (g)",         # Rename x-axis
    y = "<b>ETS</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = bodymass, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = bodymass, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = bodymass, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = bodymass, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$bodymass, na.rm = TRUE), max(data$bodymass, na.rm = TRUE), by = 50)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$ETS, na.rm = TRUE)))



p9

# 1.1 - LEAK

## Fitting the model: year and sex

mod <- lmer(LEAK~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
#newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
#newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

## Plotting the model

 color_palette <- c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9")

# Calculate the count of observations per sex and year
count_data <- data %>%
  group_by(sex, year) %>%
  summarise(N = n(), .groups = "drop")


# Create the plot
p10 <- ggplot(data = data, aes(x = sex, y = LEAK, shape = sex)) +
  geom_jitter(position = position_jitter(width = 0.08), aes(color = factor(sex)), size = 2.7, alpha = 0.3, show.legend = "right") + 
  facet_wrap(~year) +
  labs(
    x = "Sex",        # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1"  
  ) +
  scale_color_manual(values = color_palette) + 
  scale_x_discrete(labels = c("F" = "Females", "M" = "Males")) +
  geom_point(data = newdat, aes(x = sex, y = fit, color = sex),
             position = position_dodge(width = 0.5), stroke = 1.5,  size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat, aes(x = sex, y = fit, ymin = lwr, ymax = upr, color = sex),
                width = 0.2, color = "black", position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_text(data = count_data, aes(x = sex, y = min(data$LEAK, na.rm = TRUE) - 5, label = paste("n =", N)), 
            position = position_dodge(width = 1), size = 3, alpha = 0.5, color = "black", fontface = "bold", show.legend = FALSE) +  # Add count labels
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE))) +
  custom_theme 

# View the plot
p10


## Fitting the model: AGE

mod <- lmer(LEAK ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = seq(6, 18, length = 100),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$age) & !is.na(data$LEAK), ])
# Update plot
p11 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = age, y = LEAK, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Age (years)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = age, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = age, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = age, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = age, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))



p11
# Save the plot

## Fitting the model: bodymass

mod <- lmer(LEAK ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$bodymass, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = seq(445, 772, length = 100) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$bodymass) & !is.na(data$LEAK), ])
# Update plot
p12 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = bodymass, y = LEAK, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Body mass (g)",         # Rename x-axis
    y = "<b>LEAK</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = bodymass, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = bodymass, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = bodymass, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = bodymass, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$bodymass, na.rm = TRUE), max(data$bodymass, na.rm = TRUE), by = 50)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$LEAK, na.rm = TRUE)))



p12

# 1.1 - FCR1

## Fitting the model: year and sex

mod <- lmer(FCR1~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)


## Plotting the model

color_palette <- c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9")


# Calculate the count of observations per sex and year
count_data <- data %>%
  group_by(sex, year) %>%
  summarise(N = n(), .groups = "drop")


# Create the plot
p13 <- ggplot(data = data, aes(x = sex, y = FCR1, shape= sex )) +
  geom_jitter(position = position_jitter(width = 0.08), aes(color = factor(sex)), size = 2.7, alpha = 0.3, show.legend = "right") + 
  facet_wrap(~year) +
  labs(
    x = "Sex",        # Rename x-axis
    y = "<b>FCR1</b> pmol O2sec-1ml-1"  
  ) +
  scale_color_manual(values = color_palette) + 
  scale_x_discrete(labels = c("F" = "Females", "M" = "Males")) +
  geom_point(data = newdat, aes(x = sex, y = fit, color = sex),
             position = position_dodge(width = 0.5), stroke = 1.5,  size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat, aes(x = sex, y = fit, ymin = lwr, ymax = upr, color = sex),
                width = 0.2, color = "black", position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_text(data = count_data, aes(x = sex, y = min(data$FCR1, na.rm = TRUE) - 0.01, label = paste("n =", N)), 
            position = position_dodge(width = 1), size = 3, alpha = 0.5, color = "black", fontface = "bold", show.legend = FALSE) +  # Add count labels
  scale_y_continuous(limits = c(0, max(data$FCR1, na.rm = TRUE))) +
  custom_theme 

# View the plot
p13


## Fitting the model: AGE
mod <- lmer(FCR1 ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$age, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = seq(6, 18, length = 100),
                      bodymass = mean(data$bodymass) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$age) & !is.na(data$FCR1), ])
# Update plot
p14 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = age, y = FCR1, color = sex, fill = sex, shape= sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Age (years)",         # Rename x-axis
    y = "<b>FCR1</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = age, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = age, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = age, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = age, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$FCR1, na.rm = TRUE)))



p14
# Save the plot

## Fitting the model: bodymass

mod <- lmer(FCR1 ~sex +age + bodymass + year + (1|nest), 
            data = data,
            na.action = na.omit)


# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

# Transform effects from log-link function to natural values using 'exp()'
round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$bodymass, na.rm = TRUE)

# Display the range
range_hg
newdat <- expand.grid(year = levels(data$year),
                      sex= levels(data$sex),
                      age = mean(data$age),
                      bodymass = seq(445, 772, length = 100) )

Xmat <- model.matrix(~  sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other


#Plotting the model 

# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")


total_count <- nrow(data[!is.na(data$bodymass) & !is.na(data$FCR1), ])
# Update plot
p15 <- ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = bodymass, y = FCR1, color = sex, fill = sex, shape = sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "Body mass (g)",         # Rename x-axis
    y = "<b>FCR1</b> pmol O2sec-1ml-1" , # Rename y-axis
    color = "Sex",
    fill = "Sex"
  ) +
  custom_theme +
  theme(legend.position = "none")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = bodymass, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = bodymass, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = bodymass, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = bodymass, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$bodymass, na.rm = TRUE), max(data$bodymass, na.rm = TRUE), by = 50)) + # Set age ticks by 3
  scale_y_continuous(limits = c(0, max(data$FCR1, na.rm = TRUE)))



p15


#Plotting aLL figures. 


# Combine plots with labels
combined_plot <- p1 + p2 + p3 + p4 +p5 + p6 + p10+p11+p12+p13+p14+p15+
  theme(
    plot.tag = element_text(face = "bold")  # Make the tags bold
  ) +
  plot_annotation(
    tag_levels = "a"  # Automatically labels A, B
  )
# Ensure three columns
combined_plot <- combined_plot + plot_layout(ncol = 3)

# Save the combined plot
ggsave(
  filename = "suppinfo/FigS1.png",  # Save as a file
  plot = combined_plot,            # The combined plot object
  width = 8,                       # Width in inches
  height = 10,                      # Height in inches
  units = "in",# Units for dimensions
  dpi = 300   
)


## Figure S2 ----

# Correlation between trophic levels and foraging areas 
mod<-lm(d15N~d13C, data=data)


summary(mod)$sigma

nsim<-10000
bsim<-sim(mod, n.sim=nsim)

round(apply(bsim@coef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_age <- range(data$d13C, na.rm = TRUE)
range_age

newdat <- expand.grid(d13C = seq(-19.39, -18.18,length=100))

Xmat <- model.matrix(~d13C, data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@coef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%coef(mod))

head(newdat)


# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

total_count <- nrow(data[!is.na(data$d13C) & !is.na(data$d15N), ])
# Update plot
ss2<-ggplot() + # Start without default data or aes
  geom_jitter(data = data, aes(x = d13C, y =d15N, color = sex, fill = sex, shape=sex), 
              width = 0, height = 0,  # Adjust width to spread points slightly
              alpha = 0.7, size = 2.7)  +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +  # Use custom fill palette
  labs(
    x = "δ¹³C",         # Rename x-axis
    y = "δ¹⁵N", # Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape= "Sex"
  ) +
  custom_theme +
  theme(legend.position = "right")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5) + 
  geom_line(data = newdat, aes(x = d13C, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = d13C, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  geom_line(data = newdat, aes(x = d13C, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = d13C, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = seq(min(data$d13C, na.rm = TRUE), max(data$age, na.rm = TRUE), by = 1)) + # Set age ticks by 3
  scale_y_continuous(limits = c(7, max(data$d15N, na.rm = TRUE)))



ss2 
# Save the plot
ggsave("suppinfo/FigS2.png", ss2 ,width = 4, height = 3, units = "in",   dpi = 300   )


## Figure S3----

## Fitting the model

mod <- lmer(CMR ~ hg + sex +age + bodymass + year + (1|nest) + (1|TimeMito ), 
            data = data,
            na.action = na.omit)

# Set seed for reproducibility
set.seed(123)

# Number of simulations
nsim <- 10000

# Simulate data from the model  using 'sim()' function
bsim <- sim(mod, n.sim = nsim)
str(bsim)  # Show structure of the simulated data

round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 2)

# Calculate the range of values for inddependent variable
range_hg <-range(data$hg, na.rm = TRUE)

# Display the range
range_hg

newdat <- expand.grid(hg= seq(1.91, 11.77,length=100),
                      sex =levels(data$sex),
                      age= mean(data$age),
                      bodymass = mean(data$bodymass),
                      year=levels(data$year))

Xmat <- model.matrix(~ hg + sex + age + bodymass + year 
                     , data=newdat)

fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- Xmat%*%bsim@fixef[i,]
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
newdat$fit <- (Xmat%*%fixef(mod))

head(newdat)

# LEAK vs Mercury with credible intervals and no legend
newdat<-newdat %>% filter(sex=="Male") %>%  droplevels()
newdat<-newdat %>% filter(year=="2021") %>%  droplevels() #leave only the references of the other

## Plotting the model


# Define custom color palette
color_palette <- c("#E69F00", "#56B4E9")

total_count <- nrow(data[!is.na(data$hg) & !is.na(data$CMR), ])
# Update plot
p_CMR_hg <-ggplot() + 
  geom_point(data = data, aes(x = hg, y = CMR, color = sex, shape= sex), alpha = 0.7, size = 2.7) +
  scale_color_manual(values = color_palette) + # Use custom color palette for points
  scale_fill_manual(values = color_palette) +
  labs(
    x = "Mercury (µg/mg)",        
    y = "<b>CMR</b> pmol O2sec-1ml-1" ,# Rename y-axis
    color = "Sex",
    fill = "Sex",
    shape= "Sex"
  ) +
  custom_theme +
  geom_line(data = newdat, aes(x = hg, y = fit), linetype = "solid", color = "black", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = hg, y = lwr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_line(data = newdat, aes(x = hg, y = upr), linetype = "dashed", color = "grey", alpha = 0.5, size = 1) +
  geom_ribbon(data = newdat, aes(x = hg, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.2) + # Shade between lwr and upr
  scale_x_continuous(breaks = seq(min(data$hg, na.rm = TRUE), max(data$hg, na.rm = TRUE), by = 3)) + # Set age ticks by 3
  theme(legend.position = "right")+
  geom_text(
    aes(x = Inf, y = Inf, label = paste("n =", total_count)),
    size = 3, color = "black", fontface = "bold", hjust = 1.5, vjust = 1.5
  )

p_CMR_hg
# Save the plot
ggsave("suppinfo/FigS3.png", p_CMR_hg, width = 4, height = 3, units = "in",   dpi = 300  )

