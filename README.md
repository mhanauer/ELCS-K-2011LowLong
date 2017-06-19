# Jags-Ymet-XmetSsubj-MrobustHier.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , x10Name = "X10", x11Name = "X11", x12Name = "X12", x13Name="X13" , x14Name = "X14", x15Name = "X15", x16Name = "X16", x17Name = "X17", yName="y" , sName="s" ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault) { 
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x10 = data[,x10Name]
  x11 = data[,x11Name]
  x12 = data[,x12Name]
  x13 = data[,x13Name]
  x14 = data[,x14Name]
  x15 = data[,x15Name]
  x16 = data[,x16Name]
  x17 = data[,x17Name]
  # Convert sName to consecutive integers:
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x10 = x10,
    x11 = x11,
    x12 = x12,
    x13 = x13,
    x14 = x14,
    x15 = x15,
    x16 = x16,
    x17 = x17,
    y = y ,
    s = s ,
    Nsubj = max(s)  # should equal length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
  Ntotal <- length(y)
  x10m <- mean(x10)
  x11m <- mean(x11)
  x12m <- mean(x12)
  x13m <- mean(x13)
  x14m <- mean(x14)
  x15m <- mean(x15)
  x16m <- mean(x16)
  x17m <- mean(x17)
  ym <- mean(y)
  x10sd <- sd(x10)
  x11sd <- sd(x11)
  x12sd <- sd(x12)
  x13sd <- sd(x13)
  x14sd <- sd(x14)
  x15sd <- sd(x15)
  x16sd <- sd(x16)
  x17sd <- sd(x17)
  
  ysd <- sd(y)
  for ( i in 1:length(y) ) {
  zx10[i] <- ( x10[i] - x10m ) / x10sd
  zx11[i] <- ( x11[i] - x11m ) / x11sd
  zx12[i] <- ( x12[i] - x12m ) / x12sd
  zx13[i] <- ( x13[i] - x13m ) / x13sd
  zx14[i] <- ( x14[i] - x14m ) / x14sd
  zx15[i] <- ( x15[i] - x15m ) / x15sd
  zx16[i] <- ( x16[i] - x16m ) / x16sd
  zx17[i] <- ( x17[i] - x17m ) / x17sd
  zy[i] <- ( y[i] - ym ) / ysd
  }
  }
  # Specify the model for standardized data:
  model {
  for ( i in 1:Ntotal ) {
  zy[i] ~ dt( zbeta0[s[i]] + zbeta10[s[i]] * zx10[i] + zbeta11[s[i]] * zx11[i] + zbeta12[s[i]] * zx12[i] + zbeta13[s[i]] * zx13[i] + zbeta14[s[i]] * zx14[i] + zbeta15[s[i]] * zx15[i] + zbeta16[s[i]] * zx16[i] + zbeta17[s[i]] * zx17[i], 1/zsigma^2 , nu )
  }
  for ( j in 1:Nsubj ) {
  zbeta10[j] ~ dnorm( zbeta10mu , 1/(zbeta10sigma)^2 )
  zbeta11[j] ~ dnorm( zbeta11mu , 1/(zbeta11sigma)^2 )
  zbeta12[j] ~ dnorm( zbeta12mu , 1/(zbeta12sigma)^2 )
  zbeta13[j] ~ dnorm( zbeta13mu , 1/(zbeta13sigma)^2 )
  zbeta14[j] ~ dnorm( zbeta14mu , 1/(zbeta14sigma)^2 )
  zbeta15[j] ~ dnorm( zbeta15mu , 1/(zbeta15sigma)^2 )
  zbeta16[j] ~ dnorm( zbeta16mu , 1/(zbeta16sigma)^2 )
  zbeta17[j] ~ dnorm( zbeta17mu , 1/(zbeta17sigma)^2 )
  }
  # Priors vague on standardized scale:
  zbeta10mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta11mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta12mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta13mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta14mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta15mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta16mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta17mu ~ dnorm( 0 , 1/(10)^2 )
  
  zsigma ~ dnorm( 1.0E-3 , 1.0E+3 )
  zbeta10sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta11sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta12sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta13sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta14sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta15sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta16sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta17sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  
  nu ~ dexp(1/30.0)
  # Transform to original scale:
  for ( j in 1:Nsubj ) {
  beta10[j] <- zbeta10[j] * ysd / x10sd
  beta11[j] <- zbeta11[j] * ysd / x11sd
  beta12[j] <- zbeta12[j] * ysd / x12sd
  beta13[j] <- zbeta13[j] * ysd / x13sd 
  beta14[j] <- zbeta14[j] * ysd / x14sd
  beta15[j] <- zbeta15[j] * ysd / x15sd
  beta16[j] <- zbeta16[j] * ysd / x16sd
  beta17[j] <- zbeta17[j] * ysd / x17sd
  beta0[j] <- zbeta0[j] * ysd  + ym - zbeta10[j] * x10m * ysd / x10sd + zbeta11[j] * x11m * ysd / x11sd + zbeta12[j] * x12m * ysd / x12sd + ym - zbeta13[j] * x13m * ysd / x13sd + zbeta14[j] * x14m * ysd / x14sd + zbeta15[j] * x15m * ysd / x15sd + zbeta16[j] * x16m * ysd / x16sd + zbeta17[j] * x17m * ysd / x17sd 
  }

  beta10mu <- zbeta10mu * ysd / x10sd
  beta11mu <- zbeta11mu * ysd / x11sd
  beta12mu <- zbeta12mu * ysd / x12sd
  beta13mu <- zbeta13mu * ysd / x13sd
  beta14mu <- zbeta14mu * ysd / x14sd
  beta15mu <- zbeta15mu * ysd / x15sd
  beta16mu <- zbeta16mu * ysd / x16sd
  beta17mu <- zbeta17mu * ysd / x17sd
  
  
  beta0mu <- zbeta0mu * ysd  + ym -  zbeta10mu * x10m * ysd / x10sd + zbeta11mu * x11m * ysd / x11sd + zbeta12mu * x12m * ysd / x12sd + zbeta12mu * x12m * ysd / x12sd + zbeta13mu * x13m * ysd / x13sd + zbeta14mu * x14m * ysd / x14sd + zbeta15mu * x15m * ysd / x15sd + zbeta16mu * x16m * ysd / x16sd + zbeta17mu * x17m * ysd / x17sd
  sigma <- zsigma * ysd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c("beta0"  "beta10",  "beta11",  "beta12", "beta13" ,  "beta14","beta15" , "beta16","beta17", "beta0mu" , "beta10mu", "beta11mu", "beta12mu","beta13mu" , "beta14mu" , "beta15mu", "beta16mu","beta17mu", "zbeta0", "zbeta10", "zbeta11", "zbeta12", "zbeta13" , "zbeta14" , "zbeta15", "zbeta16", "zbeta17",  "zbeta0mu" ,  "zbeta10mu", "zbeta11mu", "zbeta12mu", "zbeta13mu" , "zbeta14mu" ,"zbeta15mu", "zbeta16mu", "zbeta17mu", "zsigma", "sigma", "nu" , "zbeta0sigma" , "zbeta10sigma" , "zbeta11sigma" , "zbeta12sigma", "zbeta13sigma" , "zbeta14sigma","zbeta15sigma", "zbeta16sigma","zbeta17sigma")
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 2000
  runJagsOut <- run.jags( method="parallel" ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}
