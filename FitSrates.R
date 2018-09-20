 # Fit Survival rates
require(runjags)
require(parallel)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(reshape2)
require(readxl)
# Load data ------------------------------------------------------------------
# 
rm(list = c())
dfSx = read_xlsx("../Data/Sx_logit_studies.xlsx")
dfWr = read_xlsx("../Data/Wr_logit_studies.xlsx")
dfPK = read_xlsx("../Data/Est_pK_Areas.xlsx")
dfLam = read_xlsx("../Data/Est_Lam_Areas.xlsx")
#
# Process data for running JAGS ----------------------------------------------
#
Nareas = dim(dfLam)[1]
Nyrs =  dim(dfLam)[2]-2
Years = variable.names(dfLam); Years = Years[3:dim(dfLam)[2]]
Years = as.numeric(gsub("Y", "", Years))
MinYr = min(Years)
YearsN = Years - MinYr + 1
# Proportion k, pK
pK = as.matrix(dfPK[,3:dim(dfPK)[2]])
# Log lambda data from surveys
LogLam = numeric()
Asrv = numeric()
Ysrv = numeric()
ii = 0
for (i in 1:Nareas){
  for (j in 1:Nyrs){ 
    ii = ii+1
    Asrv[ii] = i
    Ysrv[ii] = j
    LogLam[ii] = as.numeric(log(dfLam[i,j+2]))
  }
}
# Extreme log-lambda values at 1% and 100% K
Nsurv = length(LogLam)
LogLamAv = numeric()
LogLamAv[1] = 0.182
LogLamAv[2] = 0
# 
# Survival estimates
NSxEst = dim(dfSx)[1]
SxEst = dfSx$lgt_mn
tauESx = 1/(0.9*dfSx$lgt_sd)^2
Sstage = dfSx$Stage
Sarea = dfSx$Area
Syear = dfSx$Year - MinYr + 1
#
# Wean rate estimates
NWrEst = dim(dfWr)[1]
WrEst = dfWr$lgt_mn
tauEWr = 1/(0.9*dfWr$lgt_sd)^2
Warea = dfWr$Area
Wyear = dfWr$Year - MinYr + 1
#
# Create dummy pop vectors corresponding to ssd for low and high densities
n1 = matrix(c(235, 226,  66, 230, 200,  43,
              228, 246,  76, 219, 194,  37),
            byrow=T,ncol=6)
n1 = t(n1)
# Set up Jags inputs -----------------------------------------------------
#
fitmodel = c("FitSrates.jags")
#  
jags.data <- list(Nareas=Nareas,Nyrs=Nyrs,LogLam=LogLam,Asrv=Asrv,
                  Ysrv=Ysrv,Nsurv=Nsurv,LogLamAv=LogLamAv,
                  NSxEst=NSxEst,SxEst=SxEst,tauESx=tauESx,
                  Sstage=Sstage,Sarea=Sarea,Syear=Syear,
                  NWrEst=NWrEst,WrEst=WrEst,tauEWr=tauEWr,
                  Warea=Warea,Wyear=Wyear,n1=n1,pK=pK)
#                  Kguess=Kguess,NKguess=NKguess,tauKg=tauKg 
#
inits <- function() list(sigS=runif(1, .1, .5))
# sigK=runif(1, .1, .5),
params <- c("sigS","rho","zeta","eps") # 
#
nsamples <- 1000
nthin <- 10
nadapt <- 500
nburnin <- 2000
cores = detectCores()
ncore = min(20,cores-1)
#cl <- makeCluster(ncore)
nc <- ncore
#
# Run JAGS to fit model---------------------------------------------
#
out <- run.jags(data = jags.data, 
                monitor = params, 
                model = fitmodel,
                adapt=nadapt,
                inits = inits,
                n.chains = nc, 
                thin = nthin, 
                sample = nsamples, 
                burnin = nburnin,
                method="parallel") #inits = inits, 
#
# Diagnostic plots -------------------------------------------------
#
mcmc <- as.mcmc.list(out)
vn = varnames(out[["mcmc"]])
Nsims = length(as.matrix(mcmc)[,1])
sumstats = summary(out)

plot(out, vars = "sig", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "rho", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "zeta", plot.type = c("trace", "histogram"),
          layout = c(1,2))
BR = numeric(); Dur = numeric()
Dur[1] <- 3      # Stage duration for sub-adult 
Dur[2] <- 7      # Stage duration for adult
BR[1] <- .98     # Birth rate for adults
BR[2] <- .85     # Birth rate for aged-adults
pupDbase <- .15  # Base hazards for pups (as inst./ death rate), low pop densities
#
S1lg_m = matrix(nrow=5,ncol = Nyrs)
S2lg_m = matrix(nrow=5,ncol = Nyrs)
S3lg_m = matrix(nrow=5,ncol = Nyrs)
S4lg_m = matrix(nrow=5,ncol = Nyrs)
S5lg_m = matrix(nrow=5,ncol = Nyrs)
S6lg_m = matrix(nrow=5,ncol = Nyrs)
Wrlg_m = matrix(nrow=5,ncol = Nyrs)
S1lg_s = matrix(nrow=5,ncol = Nyrs)
S2lg_s = matrix(nrow=5,ncol = Nyrs)
S3lg_s = matrix(nrow=5,ncol = Nyrs)
S4lg_s = matrix(nrow=5,ncol = Nyrs)
S5lg_s = matrix(nrow=5,ncol = Nyrs)
S6lg_s = matrix(nrow=5,ncol = Nyrs)
Wrlg_s = matrix(nrow=5,ncol = Nyrs)
for (i in 1:(Nareas-1)){
  for (j in 1:Nyrs){
    eps = as.matrix(mcmc)[,which(vn==paste0("eps[",i,",",j,"]"))]
    zeta = as.matrix(mcmc)[,which(startsWith(vn,"zeta"))]
    rho =  as.matrix(mcmc)[,which(startsWith(vn,"rho"))]
    Wr = exp(-pupDbase*exp(rho*pK[i,j] + eps))
    S1 =  exp(-zeta[,1]*exp(zeta[,2]*pK[i,j] + zeta[,3] + eps))
    S2 =  exp(-zeta[,1]*exp(zeta[,2]*pK[i,j] + eps))
    S3 =  exp(-zeta[,1]*exp(zeta[,2]*pK[i,j] + zeta[,4] + eps))
    S4 =  exp(-zeta[,1]*exp(zeta[,2]*pK[i,j] + zeta[,3] + zeta[,5] + eps))
    S5 =  exp(-zeta[,1]*exp(zeta[,2]*pK[i,j] + zeta[,6] + eps))
    S6 =  exp(-zeta[,1]*exp(zeta[,2]*pK[i,j] + zeta[,4] + zeta[,6] + eps))	
    Wrlg_m[i,j] = mean(logit(Wr)); Wrlg_s[i,j] = sd(logit(Wr))
    S1lg_m[i,j] = mean(logit(S1)); S1lg_s[i,j] = sd(logit(S1))
    S2lg_m[i,j] = mean(logit(S2)); S2lg_s[i,j] = sd(logit(S2))
    S3lg_m[i,j] = mean(logit(S3)); S3lg_s[i,j] = sd(logit(S3))
    S4lg_m[i,j] = mean(logit(S4)); S4lg_s[i,j] = sd(logit(S4))
    S5lg_m[i,j] = mean(logit(S5)); S5lg_s[i,j] = sd(logit(S5))
    S6lg_m[i,j] = mean(logit(S6)); S6lg_s[i,j] = sd(logit(S6))
  }
}

VRmats = list(Wrlg_m=Wrlg_m,Wrlg_s=Wrlg_s,S1lg_m=S1lg_m,S1lg_s=S1lg_s,
              S2lg_m=S2lg_m,S2lg_s=S2lg_s,S3lg_m=S3lg_m,S3lg_s=S3lg_s,
              S4lg_m=S4lg_m,S4lg_s=S4lg_s,S5lg_m=S5lg_m,S5lg_s=S5lg_s,
              S6lg_m=S6lg_m,S6lg_s=S6lg_s)  
save(VRmats,file="../Data/Vrates_est.rdata")
save.image(file="../Results/FitSrates_Results_1.rdata")
