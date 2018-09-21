# Fit COD impacts on stage-speficic survival
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
dfCOD = read_xlsx("../Data/COD_data.xlsx")
CODdefs = read_xlsx("../Data/COD_Groups.xlsx")
dfPK = read_xlsx("../Data/Est_pK_Areas.xlsx")
load("../Data/Vrates_est.rdata")
#
# Process data for running JAGS ----------------------------------------------
#
Nareas = max(dfCOD$Area)
Nyrs =  max(dfCOD$YearN)
Nstage = 6
Nqtr = 4
Ncod = max(CODdefs$COD_N)
Years = variable.names(dfPK); Years = Years[3:dim(dfPK)[2]]
Years = as.numeric(gsub("Y", "", Years))
MinYr = min(Years)
YearsN = Years - MinYr + 1
pK = as.matrix(dfPK[,3:dim(dfPK)[2]])
Ncasedist = dim(dfCOD)[1]
NSrates = Nareas*Nyrs*Nstage
#
Cstage = as.numeric(dfCOD$Class)
Carea = as.numeric(dfCOD$Area)
Cyear = as.numeric(dfCOD$YearN)
Cqtr = as.numeric(dfCOD$Quarter)
ncase = as.numeric(dfCOD$Ncarcs)
Casedist = as.matrix(dfCOD[,7:15])
Casedist[is.na(Casedist)] = 0
#
Sstage = numeric()
Sarea = numeric()
Syear = numeric()
Slg = numeric()
tauS = numeric()
ii = 0
for (i in 1:Nareas){
  for (j in 1:Nyrs){
    for (s in 1:Nstage){
      ii = ii+1
      Sstage[ii] = s
      Sarea[ii] = i
      Syear[ii] = j
      txt1 = paste0("Slg[ii] = VRmats$S",s,"lg_m[",i,",",j,"]")   
      eval(parse(text = txt1))  
      txt2 = paste0("sd = VRmats$S",s,"lg_s[",i,",",j,"]")  
      eval(parse(text = txt2))  
      tauS[ii] = 1/sd^2
    }
  }
}
#
# Set up Jags inputs -----------------------------------------------------
#
fitmodel = c("FitCOD.jags")
#  
jags.data <- list(Nareas=Nareas,Nyrs=Nyrs,Nstage=Nstage,Nqtr=Nqtr,
                  Ncod=Ncod,Ncasedist=Ncasedist,NSrates=NSrates,
                  Casedist=Casedist,Cstage=Cstage,Carea=Carea,
                  Cyear=Cyear,Cqtr=Cqtr,ncase=ncase,
                  Slg=Slg,Sstage=Sstage,Sarea=Sarea,
                  Syear=Syear,tauS=tauS,pK=pK)
#
inits <- function() list(sigQ=runif(1, .1, .5))
params <- c("sigQ","sigZA","sigZY","phi","B","Q","Zmn","Z") # 
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
  