
## Function to predict vital rates from JAGS output ##

createpvec <- function(jagsout, dens, obstime, norand = FALSE, climate = FALSE){
  p.vec <- numeric() 
  get_par <- function(name) {
    nameTry <- paste(name, suffix, sep=".")
    if (name %in% names(jagsout)) jagsout[name] else jagsout[nameTry]
  }
  ## Latent year effects ##
  if(norand){
    latent <- c(0, 0)
  } else {
    latent <- c(jagsout[paste0("latef[", obstime+1, ",1]")], jagsout[paste0("latef[", obstime+1, ",2]")])
    }
  p.vec["latent"] <- latent[1]
  ## Environmental effect ##
  suffix <- "surv.L"
  envef <- dens - latent[1] - get_par("beta.time.thresh") * obstime 
  if(is.data.frame(climate)){
    b.clim <- jagsout[paste0("b[", 1:8, "]")]  ## extract climate coefficients
    mu.clim <- rowSums(b.clim * climate[(obstime + 1), ])
    envef <- envef - mu.clim
  }
  p.vec["envef"] <- envef
  ## Ewe Lamb survival ##
  p.vec[c("ThrL", "SuL")] <- survival("surv.L", jagsout, envef, obstime)
  ## Ewe Yearling survival ## 
  p.vec[c("ThrY", "SuY")] <- survival("surv.Y", jagsout, envef, obstime)
  ## Ewe Adult survival ##
  p.vec[c("ThrA", "SuA")] <- survival("surv.A", jagsout, envef, obstime)
  ## Ram Lamb survival ##
  p.vec["RSuL"] <- fecund("surv.RL", jagsout, envef, obstime) 
  ## Ram Yearling survival ##
  p.vec[c("RThrY", "RSuY")] <- survival("surv.RY", jagsout, envef, obstime)
  ## Ram Adult survival ##
  p.vec[c("RThrA", "RSuA")] <- survival("surv.RA", jagsout, envef, obstime)
  ## Lamb reproduction
  p.vec["RepL"] <-fecund("rep.L", jagsout, envef, obstime, latent[2])
  ## Yearling reproduction ##
  p.vec["RepY"] <- fecund("rep.Y", jagsout, envef, obstime)
  ## Adult reproduction ##
  p.vec["RepA"] <- fecund("rep.A", jagsout, envef, obstime, latent[2])
  ## Yearling twinning ##
  p.vec["TwnY"] <- fecund("twn.Y", jagsout, envef, obstime)
  ## Adult twinning ##
  p.vec["TwnA"] <- fecund("twn.A", jagsout, envef, obstime, latent[2])
  return(p.vec)
}

# Function for survival (threshold models) #
survival <- function(suffix, pvec, env, yr, both=T){
  get_par <- function(name) {
    nameTry <- paste(name, suffix, sep=".")
    if (name %in% names(pvec)) pvec[name] else pvec[nameTry]
  }
  thr <- env >= get_par("thresh")
  sur <- get_par("asymp") + get_par("beta.time.asymp") * yr - thr * get_par("beta.env") * (env - get_par("thresh"))
  if(both==T) {return(c(thr, logistic(sur)))
  } else { return(logistic(sur))}
}

# Function for repro/twinning (mixed models) models #
fecund  <- function(suffix, pvec, env, yr, replat = F){
  get_par <- function(name) {
    nameTry <- paste(name, suffix, sep=".")
    if (name %in% names(pvec)) pvec[name] else pvec[nameTry]
  }
  fec <- get_par("asymp") + get_par("beta.time.asymp") * yr - get_par("beta.env") * env
  if(is.numeric(replat)) fec <- fec - get_par("beta.repro") * replat
  return(logistic(fec))
}

logistic <- function(p.vec)
{                
  p <- 1/(1+exp(-p.vec))                                          
  return(p)
}


predsur <- function(yr, pvec, enviro, clim){
  inp <- enviro %>% filter(obsY==yr + 1985) 
  p.vec <- createpvec(pvec, unique(inp$density), yr, climate = as.data.frame(clim))
  return(p.vec)
}

########################################################################################################################
##### Function for plotting smooth splines #####
########################################################################################################################

# jagspost is the posterior samples
# climdata is a list containing 1) a matrix of climate data and 2) a matrix of weeks
# nsamp is the number of samples from the posterior to plot
# nweek is the number of weeks (/fortnights/months etc) of climate data

plotflm <- function(jagspost, climdata, nsamp = 100, knots = 8, nweek = 41){
  ## Calculate posterior medians
  post <- do.call(rbind, jagspost$mcmc)
  betas <- post[,paste0("b[", 1:38, "]")] 
  betamodes <- apply(betas, 2, median)
  ## Sample from posterior
  betasample <- betas[sample(1:nrow(betas), nsamp),]
  ## Coefficients
  jagstemp <- jags.samples(as.jags(jagspost), variable.names=c("b", "rho"), n.iter= 2, thin=1)
  jagstemp$b[,,1] <- betamodes
  ## Set up jagam
  dataout <- data.frame(Year = factor(1985:2014), latpar = rnorm(30, 0, 1))
  dataout$pcovar <- climdata$pcovar
  dataout$weeks <- climdata$weeks
  dataout$Year <- factor(dataout$Year)
  mod <- jagam(latpar ~ s(weeks, k=knots, by=pcovar, bs="cr") + s(Year, bs="re")-1, 
               data=dataout, file="template.jags")
  ## Modes
  jamtemp <- sim2jam(jagstemp, mod$pregam)
  newdat <- data.frame(weeks = 1:nweek, pcovar =1, Year = 1985)
  Xp <- predict(jamtemp, type = "lpmatrix", newdata = newdat)[,1:knots]
  modes <- Xp %*% jagstemp$b[1:knots,1,1]
  ## Posterior samples
  jagstemp <- jags.samples(as.jags(jagspost), variable.names=c("b", "rho"), n.iter= nsamp, thin=1)
  jagstemp$b[,,1] <- t(betasample)
  Xp <- predict(jamtemp, type="lpmatrix", newdata=newdat)[,1:knots]
  p <- matrix(NA, ncol = nsamp, nrow = nweek)
  for (i in 1:nsamp){
    p[,i] <- Xp %*% jagstemp$b[1:knots,i,1]
  }
  ## Plotting
  plot(modes, type="l", xlab="Fortnight", ylab="Precipitation coefficient", lwd=2, 
       ylim=c(min(p), max(p)), cex.axis = 0.7, cex.lab=0.8)
  for(i in 1:nsamp){
    lines(p[,i], col=adjustcolor("grey", alpha=0.3))
  }
  abline(h=0, lty=2, col="red")
}