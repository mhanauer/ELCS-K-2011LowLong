# Jags-Ymet-XmetSsubj-MrobustHier.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="X" , x2Name = "X2", x3Name = "X3", x4Name = "X4", x5Name = "X5", x6Name = "X6", x7Name = "X7", x8Name = "X8", yName="y" , sName="s" ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault) { 
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  x2 = data[,x2Name]
  x3 = data[,x3Name]
  x4 = data[,x4Name]
  x5 = data[,x5Name]
  x6 = data[,x6Name]
  x7 = data[,x7Name]
  x8 = data[,x8Name]
  # Convert sName to consecutive integers:
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5 = x5,
    x6 = x6,
    x7 = x7,
    x8 = x8,
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
  xm <- mean(x)
  x2m <- mean(x2)
  x3m <- mean(x3)
  x4m <- mean(x4)
  x5m <- mean(x5)
  x6m <- mean(x6)
  x7m <- mean(x7)
  x8m <- mean(x8)
  ym <- mean(y)
  xsd <- sd(x)
  x2sd <- sd(x2)
  x3sd <- sd(x3)
  x4sd <- sd(x4)
  x5sd <- sd(x5)
  x6sd <- sd(x6)
  x7sd <- sd(x7)
  x8sd <- sd(x8)
  
  ysd <- sd(y)
  for ( i in 1:length(y) ) {
  zx[i] <- ( x[i] - xm ) / xsd
  zx2[i] <- ( x2[i] - x2m ) / x2sd
  zx3[i] <- ( x3[i] - x3m ) / x3sd
  zx4[i] <- ( x4[i] - x4m ) / x4sd
  zx5[i] <- ( x5[i] - x5m ) / x5sd
  zx6[i] <- ( x6[i] - x6m ) / x6sd
  zx7[i] <- ( x7[i] - x7m ) / x7sd
  zx8[i] <- ( x8[i] - x8m ) / x8sd
  zy[i] <- ( y[i] - ym ) / ysd
  }
  }
  # Specify the model for standardized data:
  model {
  for ( i in 1:Ntotal ) {
  zy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] + zbeta2[s[i]] * zx2[i] + zbeta3[s[i]] * zx3[i] + zbeta4[s[i]] * zx4[i]  + zbeta5[s[i]] * zx5[i] + zbeta6[s[i]] * zx6[i] + zbeta7[s[i]] * zx7[i] + zbeta8[s[i]] * zx8[i], 1/zsigma^2 , nu )
  }
  for ( j in 1:Nsubj ) {
  zbeta0[j] ~ dnorm( zbeta0mu , 1/(zbeta0sigma)^2 )  
  zbeta1[j] ~ dnorm( zbeta1mu , 1/(zbeta1sigma)^2 )
  zbeta2[j] ~ dnorm( zbeta2mu , 1/(zbeta2sigma)^2 )
  zbeta3[j] ~ dnorm( zbeta3mu , 1/(zbeta3sigma)^2 )
  zbeta4[j] ~ dnorm( zbeta4mu , 1/(zbeta4sigma)^2 )
  zbeta5[j] ~ dnorm( zbeta5mu , 1/(zbeta5sigma)^2 )
  zbeta6[j] ~ dnorm( zbeta6mu , 1/(zbeta6sigma)^2 )
  zbeta7[j] ~ dnorm( zbeta7mu , 1/(zbeta7sigma)^2 )
  zbeta8[j] ~ dnorm( zbeta8mu , 1/(zbeta8sigma)^2 )

  }
  # Priors vague on standardized scale:
  zbeta0mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta1mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta2mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta3mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta4mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta5mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta6mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta7mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta8mu ~ dnorm( 0 , 1/(10)^2 )
  
  zsigma ~ dnorm( 1.0E-3 , 1.0E+3 )
  zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta2sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta3sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta4sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta5sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta6sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta7sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta8sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  
  
  nu ~ dexp(1/30.0)
  # Transform to original scale:
  for ( j in 1:Nsubj ) {
  beta1[j] <- zbeta1[j] * ysd / xsd 
  beta2[j] <- zbeta2[j] * ysd / x2sd
  beta3[j] <- zbeta3[j] * ysd / x3sd
  beta4[j] <- zbeta4[j] * ysd / x4sd
  beta5[j] <- zbeta5[j] * ysd / x5sd
  beta6[j] <- zbeta6[j] * ysd / x6sd
  beta7[j] <- zbeta7[j] * ysd / x7sd
  beta8[j] <- zbeta8[j] * ysd / x8sd
  beta0[j] <- zbeta0[j] * ysd  + ym - zbeta1[j] * xm * ysd / xsd + zbeta2[j] * x2m * ysd / x2sd + zbeta3[j] * x3m * ysd / x3sd + zbeta4[j] * x4m * ysd / x4sd + zbeta5[j] * x5m * ysd / x5sd + zbeta6[j] * x6m * ysd / x6sd + zbeta7[j] * x7m * ysd / x7sd + zbeta8[j] * x8m * ysd / x8sd
  }
  beta1mu <- zbeta1mu * ysd / xsd
  beta2mu <- zbeta2mu * ysd / x2sd
  beta3mu <- zbeta3mu * ysd / x3sd
  beta4mu <- zbeta4mu * ysd / x4sd
  beta5mu <- zbeta5mu * ysd / x5sd
  beta6mu <- zbeta6mu * ysd / x6sd
  beta7mu <- zbeta7mu * ysd / x7sd
  beta8mu <- zbeta8mu * ysd / x8sd
  
  
  beta0mu <- zbeta0mu * ysd  + ym - zbeta1mu * xm * ysd / xsd + zbeta2mu * x2m * ysd / x2sd + zbeta3mu * x3m * ysd / x3sd + zbeta4mu * x4m * ysd / x4sd + zbeta5mu * x5m * ysd / x5sd + zbeta6mu * x6m * ysd / x6sd + zbeta7mu * x7m * ysd / x7sd + zbeta8mu * x8m * ysd / x8sd 
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
  parameters = c( "beta0" ,  "beta1","beta2" , "beta3",  "beta4", "beta5", "beta6",  "beta7",  "beta8",  "beta0mu" , "beta1mu" , "beta2mu", "beta3mu", "beta4mu", "beta5mu", "beta6mu", "beta7mu", "beta8mu", 
                  "zbeta0" , "zbeta1" , "zbeta2", "zbeta3", "zbeta4", "zbeta5", "zbeta6", "zbeta7", "zbeta8",  "zbeta0mu" , "zbeta1mu" ,"zbeta2mu", "zbeta3mu", "zbeta4mu", "zbeta5mu", "zbeta6mu", "zbeta7mu", "zbeta8mu", 
                  "zsigma", "sigma", "nu" , "zbeta0sigma" , "zbeta1sigma","zbeta2sigma", "zbeta3sigma", "zbeta4sigma" , "zbeta5sigma" , "zbeta6sigma" , "zbeta7sigma" , "zbeta8sigma")
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
