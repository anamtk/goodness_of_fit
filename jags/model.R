model {
  
  #-------------------------------------## 
  # WHWO Nest survival model
  # Based on model from Kozma et al. 2017
  # adapted by Ana Miller-ter Kuile
  # November 3, 2021
  
  # this model tests the importance
  # of forest treatment + other environmental covariates on nest survival
  # of white-headed woodpeckers
  
  # General model attributes:
  # - hierarchical structure of nest within transect
  # - crossed random effect of year
  # - Treatment covariates at the nest
  # - visit interval covariates related to nest stage
  # - additional covariates for the multiple scales of environmental factors (
  # nest, local, and landscape)
  
  # - this model includes bayesian p-value estimation for important covariates
  #-------------------------------------##
  
  #-------------------------------------## 
  #Model of nest survival ###
  #-------------------------------------##
  
  for(i in 1:n.nests) { #for each nest
    for(j in 1:n.t[i]){ #and each interval in which the nest was surveyed
      
      #observed values of y are of a 1/0 bernoulli distribution based on s,
      # the survival probability for that nest in that survey period of time t
      y[i, j] ~ dbern(mu.s[i,j]) 
      
      # set of covariates, z, as a function of random effects and covariates
      # based on nest location and interval 
      
      z[i, j] <- #hierarchical structure of intercept with hierarchical
        # centering to fix identifiability issues
        b0.nest[Nest.num[i]] + #this encapsulates multiple spatial hierarchies
        #coded into the priors for this - see below
        #Crossed random effect of year:
        b0.year[Year.num[i]] + #this is summed to zero for identifiabilty
        # see in priors below
        #categorical covariates
        # Interval categorical covariate:
        b1StageID[StageID[i, j]] + 
        # Treatment categorical covariates
        b2TreatmentID[TreatmentID[i]] +
        b3TrtTime[TrtTime[i]] +
        #nest categorical covariate
        b4SpeciesID[SpeciesID[i]] +
        #continuous covariates
        #Treatment continuous
        b[5]*NTrt[i] +
        #Nest continuouse covariates
        b[6]*NestHt[i]+
        b[7]*cosOrientation[i] +
        b[8]*InitDay[i]+
        #local continuouse covariates
        b[9]*Trees50[i] +
        b[10]*Trees2550[i] +
        b[11]*PercPonderosa[i] +
        #landscape continuous covariates
        b[12]*Tmax[i] +
        b[13]*Tmax[i]^2 +
        b[14]*ForestCV[i] +
        b[15]*Contag[i] +
        b[16]*OpenNm[i] +
        b[17]*LandHa[i] +
        b[18]*LandBu[i] +
        #interacting covariates
        b[19]*Trees50[i]*PercPonderosa[i] +
        b[20]*Trees2550[i]*PercPonderosa[i] +
        b[21]*Trees50[i]*Tmax[i] +
        b[22]*Trees2550[i]*Tmax[i] +
        b[23]*LandHa[i]*LandBu[i]
      
      #SURVIVAL PROBABILITY
      #survival probability, s, for each nest and time period
      #(break into three lines of code to ensure models plays nice)
      
      #1. daily survival:
      s.daily[i, j] <- (exp(z[i,j]))/(1+exp(z[i,j]))
      
      #2. log transform that and multiply by t (rather than doing ^t,
      # this is what can break things)
      #daily survival taken to the ^t power for the number of 
      #days in that survey period:
      log.s[i,j] <- t[i,j]*log(s.daily[i,j])
      
      #3.  retransform s back to non-log scale
      s[i, j] <- exp(log.s[i,j])
      
      #And then keep s away from 1 and 0
      mu.s[i,j] <- min(0.999, max(0.001, s[i, j]))
      
      # the resulting mu.s[i,j] value gives an estimate
      # of survival probability for that nest for that survey period of time t
      
      #-------------------------------------## 
      # Model Goodness-of-fit objects ###
      #-------------------------------------##
      
      #Create replicated data for gof
      yrep[i, j] ~ dbern(mu.s[i,j])
      
      #Residuals
      resid[i,j] <- y[i,j] - mu.s[i,j]
      
    }
    
    #-------------------------------------## 
    # Imputing missing data ###
    #-------------------------------------##
    
    #Some covariate data are msising, so use the following to model those 
    # missing data
    #Basing these distributions off of the distributions of the 
    # data for each variable
    Trees2550[i] ~ dnorm(mu.t25, tau.t25)
    Trees50[i] ~ dnorm(mu.t50, tau.t50)
    PercPonderosa[i] ~ dnorm(mu.pp, tau.pp)
    InitDay[i] ~ dnorm(mu.init, tau.init)
    cosOrientation[i] ~ dnorm(mu.orient, tau.orient)
    
    #temp is dependent on forest location
    Tmax[i] ~ dnorm(mu.tmax[Forest.num[i]], tau.tmax[Forest.num[i]])
    
  }
  
  #-------------------------------------## 
  #Prior distributions on parameters ###
  #-------------------------------------##
  
  ## HIERARCHICAL EFFECTS PRIORS
  # Random covariates:
  # Hierarchical spatial random effects
  #each level depends on the level higher than it
  #Nested spatial random structure with hierarchical centering: 
  for(n in 1:n.nests){ #nests
    #nests in points effect
    b0.nest[n] ~ dnorm(b0.transect[Transect.num[n]], tau.nest)
  } 
  
  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0, tau.transect)
  }
  
  #Crossed effect for year
  #this effect is summed to zero for identifiability issues
  
  #for every year but the last one:
  for(y in 1:(n.years-1)){
    b0.year[y] ~ dnorm( 0, tau.year)
  }
  #set the last year to be the -sum of all other years so the 
  # overall fo all year levels == 0
  b0.year[n.years] <- -sum(b0.year[1:(n.years-1)])
  
  #Random and intercept priors
  b0 ~ dnorm(0, 1E-2)
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.transect ~ dunif(0, 10)
  sig.year ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
  tau.year <- 1/pow(sig.year, 2)
  
  tau.nest ~ dgamma(0.1, 0.1)
  sig.nest <- 1/sqrt(tau.nest)
  
  #FIXED COVARIATE PRIORS
  #Categorical variables
  #this is all in relation to first treatment
  #Ensure treatment == 1 has the most observations!!
  for(s in 2:n.stages){
    b1StageID[s] ~ dnorm(0, 1E-2)
  }
  b1StageID[1] <- 0

  for(tt in 2:n.trt){
    b2TreatmentID[tt] ~ dnorm(0, 1E-2)
  }
  b2TreatmentID[1] <- 0
  
  for(t in 2:n.times){
    b3TrtTime[t] ~ dnorm(0, 1E-2)
  }
  b3TrtTime[1] <- 0
  
  for(s in 2:n.species){
    b4SpeciesID[s] ~ dnorm(0, 1E-2)
  }
  b4SpeciesID[1] <- 0
  
  for(i in 5:23){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #MISSING DATA IMPUTING PRIORS
  #Priors for mean and tau of missing covariates in the model
  mu.t25 ~ dunif(-10, 10)
  sig.t25 ~ dunif(0, 20)
  tau.t25 <- pow(sig.t25, -2)
  mu.t50 ~ dunif(-10, 10)
  sig.t50 ~ dunif(0, 20)
  tau.t50 <- pow(sig.t50, -2)
  mu.pp ~ dunif(-10, 10)
  sig.pp ~ dunif(0, 20)
  tau.pp <- pow(sig.pp, -2)
  mu.init ~ dunif(-10, 10)
  sig.init ~ dunif(0, 20)
  tau.init <- pow(sig.init, -2)
  mu.orient ~ dunif(-10, 10)
  sig.orient ~ dunif(0, 20)
  tau.orient <- pow(sig.orient, -2)

  #these need to be indexed by forest ID
  for(f in 1:n.forests){
    mu.tmax[f] ~ dunif(-10, 10)
    sig.tmax[f] ~ dunif(0, 20)
    tau.tmax[f] <- pow(sig.tmax[f], -2)
  }
  
  #-------------------------------------## 
  # Covariate P-values ###
  #-------------------------------------##
  
  #generate a 1-0 vector for each covariate
  #such that 1 = + in that iteration, 0 = - in that iteration
  # the mean of this value will tell us whether something is mostly positive
  # (high mean posterior value), mostly negative (low mean posterior value)
  # or somewhree in the middle (often 0, so 0.5 mean posterior)
  
  #generates per level of categorical variables
  zi.b1 <- step(b1StageID)
  zi.b2 <- step(b2TreatmentID)
  zi.b3 <- step(b3TrtTime)
  zi.b4 <- step(b4SpeciesID)
  
  #generate p-values for all continuous covariates
  for(i in 5:23){
    zi[i] <- step(b[i])
  }
  
}

    
      