#################################################
##simple (no age structure) transmission model
################################################

mod <- function(t,y,par){
  with(as.list(c(y,par)),{
    
    # initialise matrices & vectors
    ymat<- matrix(y, ncol=1, nrow=nw)
    dymat <- matrix(0, ncol=1, nrow=nw)
    
    ## mean worm burden
    W = sum(ymat[,1]) 
    
    ## density-dependent functions
    mmat <- funcs$matingf(W,k)
    dd <- funcs$densdep(W, k, b, rho)
    
    ## eggoutput
    epgout= rho*W*dd*a
    
    # prevalence of (detectable) female infection
    prev= (1-(1+(W*rho)/k)^-k)
      
    # prevalence of (detactable) heavy infection
    prev_h = (1-pnbinom((z/a)^(1/b), size = k, mu = W*rho)) 

    # model ODEs    
    for (i in 1:nw) { 
        if (i==1) { 
            dymat[i,1] = ((R0/rho)*mu2*W*mmat*dd) / 
              (R0hs*N1N2*W*mmat*dd + mu2/(nw*mu1)) - (nw*mu1)*ymat[i,1]
          } else if (i>1) {
            dymat[i,1] = nw*mu1*ymat[i-1,1] - (nw*mu1)*ymat[i,1]
          }
        }
      
    
    return(list(rbind(dymat),
                W=W, E=epgout, prev=prev,prevh=prev_h))
  })
}


