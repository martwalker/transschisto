#########################################
## default parameters
########################################
parameters= c(R0=5,                   # basic reproduction number
              R0hs= 10,               # R0 host -> snails
              N1N2 = 1/10,            # host-to-snail population density
              mu1= 1/5,               # adult worm mortality rate (1/mu1=average lifespan)
              mu2= 1/(1/12),          # snail mortality rate (1/mu2=average snail lifespan)
              k= 0.5,                 # overdispersion adult worms among hosts
              coverage= 0.50,         # population treatment coverage
              efficacy= 0.94,         # efficacy of PZQ
              start.tx=1,             # treatment start time
              n.tx=25,                # number of treatments
              freq.tx=1,              # frequency of treatments
              nw=1,                   # adult worm structure (not used)
              dt=1/12,                # integration time step (fixed)
              stop.t=51,              # simulation stop time
              a= exp(2.92),       # maximum egg output from 1 female (Neves et al 2021)
              b=0.4,              # density-dependent constraint on fecundity (Neves et al 2021)
              z=50,                   # threshold egg density for heavy infection
              dotx=1, 
              rho=0.5)                 # toggle treatments on/off

cntrlpar = list(accbreak=0.0125, accendem=0.5, block=10, nsamp=50000)
