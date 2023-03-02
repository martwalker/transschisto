#########################################
## utility functions
########################################

funcs <- list(
# mating probability (monogamous, May 1977)
matingf =  function(W,k){
  
  W <- max(0, W)
  
  alpha= W/(k+W)
  
  integrand <- function(theta)
  {(1-cos(theta))/(((1+alpha*cos(theta))^(1+k)))}
  intsol <- integrate(integrand, lower=0, upper=2*pi)
  
  int <- intsol$value
  I= (((1-alpha)^(1+k))/(2*pi))
  I1= 1-I*int
  
  return(I1)
},

ddred = function(n, b){
    n^(b-1) 
},

prob = function(n, W, k, rho) {
  dnbinom(x= n, size = k, mu= W*rho)
}, 

# density-dependent reduction in worm fecundity
densdep = function(W, k, b, rho){
  
  max <- max(qnbinom(0.999, size=k, mu=I(W*rho)),1)
  
  probs <- dnbinom(seq(1,max), size = k, mu= W*rho)
    
  sumsol <- sum( funcs$ddred(1:max, b=b)*probs ) / (1-(1+(W*rho)/k)^(-k))

  return(as.numeric(sumsol))
},

netred = function(W, k, b) {
  funcs$densdep2(W=W,k=k,b=b)*funcs$matingf(W=W,k=k)
},

densdep2 = function(W, k, b) {
  A=W/(W+k)
  B=exp(-b)
  p0 = (1+W/k)^(-k)
  (1-A)^k/(B*(1-p0))*((1-A*B)^(-k)-1)
},

# treament event function
eventfun = function(t, y, par) {
  with(as.list(y), {
    
    ymat <- matrix(y, ncol=1, nrow=parameters["nw"])
    ymat[,1] <- ymat[,1] - ymat[,1]*parameters["coverage"]*(parameters["efficacy"])
    
    return(c(ymat))
  })
  
},

# main function for running the simple model
runmod = function(parameters, W0=2) {
  
  
  if (parameters["dotx"]==1) {
  if(parameters["stop.t"] <= parameters["start.tx"] + parameters["n.tx"]*parameters["freq.tx"]) {
    stop("treatments run longer than simulation time")
  }
  }
  
  inits <- rep(as.numeric(W0/parameters["nw"]), parameters["nw"])
  
  t <- seq(0,parameters["stop.t"],by=parameters["dt"])
  
  if(parameters["dotx"]==1) {
  events = list(func= funcs$eventfun, time= seq(parameters["start.tx"],
                                          parameters["start.tx"]+parameters["n.tx"]*parameters["freq.tx"],
                                          parameters["freq.tx"]), par= parameters)
  
  
  out<-lsoda(y=inits,times=t, func=mod, par=parameters,  events = events)
  } else {
    
    out <- lsoda(y=inits,times=t, func=mod, par=parameters)
  
    }
  
  out <- as.data.frame(out)
  out
}, 

derivs = function(W, par) {
  mod(t=0, y=W, par=par)[[1]]
},

findbreak = function(par)
{
  W<-cntrlpar$accbreak
  while(W<cntrlpar$maxbreak)
  {
    W0 <- W
    W1 <- W0+cntrlpar$accbreak
    dW0 <- funcs$derivs(W0, par=parameters)
    dW1 <- funcs$derivs(W1, par=parameters)
    i <- diff(sign(c(dW0, dW1)))
    if(i!=0) {
      return(W0+(W1-W0)/2)
      break
    }
    W <- W + cntrlpar$accbreak
  }
  return(NA)
},

findstable = function(Wbreak, par) {
  if(is.na(Wbreak)) {
    return(NA)
  } else {
    W <- Wbreak+cntrlpar$accbreak/2
    while(W<cntrlpar$maxendem) {
      W0 <- W
      W1 <- W0+cntrlpar$accendem
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0+(W1-W0)/2)
        break
      }
      W <- W + cntrlpar$accendem
    }
    return(NA)
  }
},

findstable2 = function(par) {
  maxW <- max(as.numeric((parameters["R0"]/parameters["rho"])*parameters["mu2"]-parameters["mu2"]/
    (parameters["R0hs"]*parameters["mu1"]*parameters["N1N2"])),1)
  {
    W <- maxW-cntrlpar$accendem/2
    while(I(W - cntrlpar$accendem)>0) {
      W0 <- W
      W1 <- W0-cntrlpar$accendem
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0-(W1-W0)/2)
        break
      }
      W <- W - cntrlpar$accendem
      #print(W)
      
    }
    return(NA)
  }
},

findstable3 =function(par) {
  tmp <- funcs$runmod(parameters = par)
  W <- tmp[nrow(tmp), "W"]
  if (W<cntrlpar$accendem) {
    return(NA)
  } else {
    return(W)
  }
},

findbreak2 = function(par, Wstar)
{
  if(is.na(Wstar)) {
    return(NA)
  } else {
    W<-1.0E-12
    while(W<Wstar)
    {
      W0 <- W
      W1 <- W0+cntrlpar$accbreak
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0+(W1-W0)/2)
        break
      }
    W <- W + cntrlpar$accbreak
    }
  }
  return(NA)
},

PRCC=function(outcome,covariates){
  
  rank_outcome<-rank(outcome)
  rank_covariates<-as.data.frame(apply(covariates, 2, rank))
  
  PRCC_out<-c()
  for (par in 1:ncol(covariates)){
    xx<-rank_covariates[,par]		  
    xy<-rank_covariates[,-par]	
    
    xj<-xx-predict(lm(xx~.,data=xy))
    yy<-rank_outcome-predict(lm(rank_outcome~.,data=xy))
    
    PRCC_out[par]<-cov(xj,yy)/(sqrt(var(xj)*var(yy)))	
  }
  names(PRCC_out)<-colnames(covariates)
  PRCC_out
}



)


