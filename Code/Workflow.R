
########################################################################################################################
##### More information #####
########################################################################################################################

## For more details on the SEM approach see Hindle et al (2018) MEE (https://doi.org/10.1111/2041-210X.13085) and the 
## example code at https://doi.org/10.5281/zenodo.1312855

## For more details on the FLM approach including example code see Teller et al (2016) MEE (https://doi.org/10.1111/2041-210X.12486)

## For general details on fitting smoothers in JAGS see Wood et al (2016) Journal of Statistical Software (https://doi.org/10.18637/jss.v075.i07)

########################################################################################################################
##### Load packages and data #####
########################################################################################################################

set.seed(110189)

### Load required packages ###
requiredPackageNames <- c(
  "dplyr", "lme4", "ggplot2", "rjags", 
  "tidyr", "plotrix", "mgcv",
  "runjags", "snow", "foreach", 
  "doParallel", "data.table", "gridExtra", "cowplot")
# install (if we don't have it) and load up the package 
sapply(requiredPackageNames, function(name) {
  if(!require(name, quietly=TRUE, character.only=TRUE)) {
    install.packages(name, repo = 'http://cran.uk.r-project.org')
    library(name, quietly=TRUE, character.only=TRUE)
  }
})

### Load Soay sheep data ###
load(file = "./Data/SoayData.rda")
# NB: data is in a list for importing into JAGs
# Density is log population size and Time is year of study
# The remaining elements are the demographic data e.g. surv.L is the number of surviving ewe lambs each year and 
# n.surv.L is the total number of ewe lambs each year

### Load climate data ###
load("./Data/PrecipData.rda")
# pcovar is a matrix, where pcovar[t,w] of precipitation in year t and week w.
# weeks is a matrix giving the week indices for every year.

########################################################################################################################
##### Formatting for JAGS  #####
########################################################################################################################

# Thanks to Wood (2016) the jagam function can be used to turn a GAM specified in mgcv format to the JAGS
# code and data.

precipdata$latent <- rnorm(nrow(precipdata$pcovar), 0, 1)
precipdata$Year <- factor(1:30)

examp <- jagam(latent ~ s(weeks, k=8, by=pcovar, bs="cr") + s(Year, bs="re") -1, 
              data=precipdata, file="FLM.jags", sp.prior = "log.uniform")

dataUseList$X <- examp$jags.data$X
dataUseList$S1 <- examp$jags.data$S1
dataUseList$zero <- rep(0, 8)

# The example code from the file created here ("FLM.jags") was used to write the JAGS model ("SEM-FLM.jags").  

########################################################################################################################
##### Running model in JAGS #####
########################################################################################################################

## Specify which parameters to save

par <- c("asymp.", "beta.env.", "thresh.", "beta.time.asymp.")
stages <- c("L", "Y", "A", "RL", "RY", "RA")
rate <- c("surv.", "rep.", "twn.")

parKeep <- c(do.call(paste0, expand.grid(par, rate, stages)), 
             "b", "rho", "lambda", "env.effect", "clim.ef",
             "beta.time.thresh", "sig.latent", "latef", "cor.lat", 
             "beta.repro.rep.L", "beta.repro.rep.A", "beta.repro.twn.A")

## Give some initial values
inits <- list(b=c(examp$jags.ini$b[1:8], rep(NA, 30)), rho = examp$jags.ini$rho[1:2], 
              beta.repro.rep.L = 1, beta.repro.rep.A = 1, beta.repro.twn.A=1)

## Warning: this takes hours rather than minutes to run
## "Output.rda" file is provided if you don't want to run this part yourself
postSample <- run.jags("./Code/SEM-FLM.jags",
                       adapt=5000, burnin= 1e5, sample= 1e3, thin=2e3, n.chains=2,
                       monitor=parKeep, data=dataUseList, inits = inits,
                       modules="glm",
                       method="rjags")


save(postSample, file = "Output.rda")

########################################################################################################################
##### Initial checks #####
########################################################################################################################

## Carry out usual checks of convergence etc 

par(mfrow=c(3,2))
plot(postSample)

gelman.diag(postSample)


########################################################################################################################
##### Plotting output #####
########################################################################################################################

## This rest of this file includes code to make some of the figures from the manuscript including 
# Figure 3a) - Plot of vital rates against population density
# Figure 4) - Plot of observed vs. predicted vital rates
# Figure 5) - Climatic effects using FLMs

source("./Code/PlotFunctions.R")
load("./Data/Observed.rda")

base <- do.call(rbind, postSample$mcmc); bmeds <- apply(base, 2, median)

# NB: there may be minor differences between the plots here and those in the manuscript as we've used 
# the base model (with no climate predictors) for Figs 3 and 4 in the manuscript and here we use one of 
# the FLM models

########################################################################################################################
##### Posterior checks (e.g. fig 4 in manuscript) #####
########################################################################################################################

## We recommend carrying out a range of posterior checks to check various aspects of the model e.g. here we plot the 
## observed vs predicted vital rates.

prds <- data.frame(do.call("rbind", lapply(0:29, predsur, bmeds, enviro = pdata, clim = examp$jags.data$X[,1:8]))) %>% 
  mutate(obsY=1985:2014) %>%
  gather(param, pred, c(SuL, SuY, SuA, RepL, RepY, RepA, TwnY, TwnA, RSuL, RSuY, RSuA)) %>% 
  mutate(Age= ifelse(grepl("L", param), "Lambs",ifelse(grepl("Y", param), "Yearlings", "Adults")),
         Vital=ifelse(grepl("Rep", param), "Reproduction", ifelse(grepl("Su", param), "Survival", "Twinning")),
         Sex = ifelse(grepl("RS", param), "Rams", "Ewes")) %>% dplyr::select(Age, Vital, Sex, obsY, pred)

toplot <- pdata %>% inner_join(prds, by=c("Age", "Sex", "obsY", "Vital")) %>% ungroup() %>% 
  mutate(Age = factor(Age, levels = c("Lambs", "Yearlings", "Adults")))

## Example plot
ggplot(filter(toplot, Vital=="Survival"), aes(x=Prop, y=pred)) + scale_color_manual(values=c("#56B4E9", "#999999", "#E69F00")) + 
  geom_point(aes(colour=Age, shape=Sex), alpha = 0.5) + geom_abline(linetype=2) + labs(x="Observed survival", y="Predicted survival") + 
  theme(aspect.ratio = 1, text = element_text(size = 9), axis.text= element_text(size = 8)) + xlim(c(0, 1)) + ylim(c(0, 1))

## We also recommend the use of cross validation to compare the predictive performance of different models  (e.g. comparing 
## each model with a climatic driver to a base model with no drivers) to reduce the risk of overfitting

########################################################################################################################
##### Plot density vs. vital rates (fig. 3a in manuscript) #####
########################################################################################################################

## Predict vital rates (for mid year of study)
pred <- data.frame(do.call(rbind, mapply(createpvec, dens = seq(2.25, 2.85, length = 100), 
                                         MoreArgs = list(jagsout = bmeds, obstime = 14, norand = TRUE, 
                                                         climate = as.data.frame(examp$jags.data$X[,1:8])), ## Using mid-year from the study
                                         SIMPLIFY = FALSE))) %>% select(-contains("Thr"), -c(latent, envef)) %>% 
  mutate(density = seq(2.25, 2.85, length = 100)) %>% gather(Vital, Prop, SuL:TwnA) %>% 
  mutate(Age = ifelse(grepl("L", Vital), "Lambs", ifelse(grepl("Y", Vital), "Yearlings", "Adults")), 
         Sex = ifelse(grepl("RS", Vital), "Rams", "Ewes"), 
         Vital = ifelse(grepl("Su", Vital), "Survival", ifelse(grepl("Rep", Vital), "Reproduction", "Twinning"))) 

## Example plot
ggplot(filter(pdata, Vital =="Survival"), aes(x=density, y = Prop, colour = Age, group = interaction(Age, Sex))) +
  geom_point(aes(shape = Sex, size = obsY), alpha = 0.65) + 
  theme(aspect.ratio=1, text = element_text(size = 9), axis.text = element_text(size = 8),
        legend.title = element_text(size=8), legend.key.size = unit(0.5, 'lines'), 
        strip.text = element_text(hjust = 0, vjust = -1)) + 
  scale_size_continuous(breaks=seq(1985, 2010, by = 10), range = c(0.5, 1.75)) + 
  geom_line(data = filter(pred, Vital=="Survival"), aes(linetype = Sex)) + scale_color_manual(values=c("#56B4E9", "#999999", "#E69F00")) + 
  labs(x='Density (log scale)', y='Proportion surviving', size = "Year")

########################################################################################################################
##### Plotting FLMs (fig. 5 in manuscript) #####
########################################################################################################################

par(mfrow=c(1,1))
plotflm(postSample, precipdata)

## This plots an example FLM output (i.e. Figure 5 in the manuscript)
## The median is in black and samples from the posterior in grey (use nsamp argument to change number of samples)

