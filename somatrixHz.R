somatrixHz <- function(pK,sigma){
  #SOMATRIX function to generate 4x4 projection matrix for sea otters
  #  (2 age classes for 2 sexes = 4stage) using proportional hazards
  #  with density dependenze, possible additional additive hazards and 
  #  environmental stochasticity
  #
  # EXPLANATION:
  # Uses proportional hazards formulation to estimate survival
  #   gamma = exp(sum(zeta0 + zeta1 + zeta2... + AiHz + eps))
  #   Sx = exp(-gamma)
  #     where AiHz = user supplied additional additive hazards 
  #     (which can be = 0 if no additional hazards)
  #     and eps = stochasticity, eps~rnorm(1,0,sigma)
  #     (if user-supplied sigma=0 then no stochasticity)
  # NOTE: user-supplied sigma represents std dev in log haz,
  #  and is ~5x larger than associated sd of log(lambda) 
  #
  # Produces a 2-sex matrix for 2 age classes,
  #  Sub-adults: 0.5-2.5 yr olds (2 year stage duration)
  #  Adults: 2.5-19.5 year olds (17 year stage duration)
  #  (assumes first pup born to females on 3rd birthday)
  # INPUT: vector of 8 elements containing vital rates:
  #  br = birth rates for adult females
  #  wr = weaning success rate for adult females
  #  fs = female survival rates (stages J and A)
  #  ms = male survival rates (stages J and A)
  # OUTPUT:
  # lambda = algebraic value of lambda
  # ssd = associated stable stage distribution  
  # M = 4 by 4 projection matrix
  # vrts = vector of vital rates
  #
  if(sigma==0){
    eps = 0
  }else{
    eps = max(-.5,min(.5,rnorm(1,0,sigma))) # env. stochasticity
  }
  # Rho parameter, pup survival proportional hazards function (log hazards)
  rho = 1.1482881         # log haz ratio, density-dependent hazards for pups (multiplied by prpn K) 
  zeta = numeric()
  # Zeta parameters for survival proportional hazards function (log hazards)
  zeta[1] = 0.0312688  # Base hazards (adult female), as instantaneous death rate 
  zeta[2] = 1.2023682	# log haz ratio, density-dependent hazards (multiplied by prpn K)
  zeta[3] = 0.4569099  # log haz ratio, for juveniles/sub-adults relative to adults 
  zeta[4] = 1.0864169  # log haz ratio, for aged adults relative to adults
  zeta[5] = 0.3289521    # Log haz ratio for sub-adult males relative to sub-adult females
  zeta[6] = 0.4188282    # Log haz ratio for adult males relative to adult females
  pupDbase = .15  # Base hazards for pups (as inst./ death rate), low pop densities
  BR = c(.98,.85) # Birth rates for adults and aged-adults
  wr = exp(-pupDbase*exp(rho*pK + eps))
  Rad = 0.5*BR[1]*wr
  Raa = 0.5*BR[2]*wr
  fsSA = exp(-zeta[1]*exp(zeta[2]*pK + zeta[3] + eps))
  fsA = exp(-zeta[1]*exp(zeta[2]*pK + eps))
  fsAA = exp(-zeta[1]*exp(zeta[2]*pK + zeta[4] + eps))
  msSA = exp(-zeta[1]*exp(zeta[2]*pK + zeta[3] + zeta[5] + eps))
  msA = exp(-zeta[1]*exp(zeta[2]*pK + zeta[6] + eps))
  msAA = exp(-zeta[1]*exp(zeta[2]*pK + zeta[4] + zeta[6] + eps))
  Dur = c(3,7) # stage durations for sub-adult and adult
  lm = max(.5,min(1.25,-1.06 + 2.3*fsA)); # Initial lambda estimate 
  # Parameterize matrix and solve (use iterations to stabalize lambda)
  for (i in 1:3){
    gf1 = (min(.999,fsSA/lm)^Dur[1] - min(.999,fsSA/lm)^(Dur[1]-1))/(min(.999,fsSA/lm)^Dur[1] -1)
    gm1 = (min(.999,msSA/lm)^Dur[1] - min(.999,msSA/lm)^(Dur[1]-1))/(min(.999,msSA/lm)^Dur[1] -1)
    gf2 = (min(.999,fsA/lm)^Dur[2] - min(.999,fsA/lm)^(Dur[2]-1))/(min(.999,fsA/lm)^Dur[2] -1)
    gm2 = (min(.999,msA/lm)^Dur[2] - min(.999,msA/lm)^(Dur[2]-1))/(min(.999,msA/lm)^Dur[2] -1)    
    M = matrix(c(fsSA*(1-gf1),  fsA*Rad,      fsAA*Raa,   0,            0,            0, 
                 fsSA*gf1,      fsA*(1-gf2),  0,          0,            0,            0, 
                 0,             fsA*gf2,      fsAA,       0,            0,            0,    
                 0,             fsA*Rad,      fsAA*Raa,   msSA*(1-gm1), 0,            0,
                 0,             0,            0,          msSA*gm1,     msA*(1-gm2),  0,
                 0,             0,            0,          0,            msA*gm2,      msAA),
              byrow=T,ncol=6) 
    lm=eigen(M)$values[1]    # lambdas=vector of eigenvalues
  }
  lambda = lm 
  W=eigen(M)$vectors          # W=matrix of right eigenvectors 
  w=abs(W[,1])					      # w=stable distribution, unscaled
  ssd = w/sum(w)                # w=stable distribution, scaled
  result <- list(lam=lambda,SSD = ssd,M=M,
                 vrts = c(wr,fsSA,fsA,fsAA,msSA,msA,msAA))
  return(result)
}
