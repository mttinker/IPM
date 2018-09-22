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
Stagefilter = matrix(1,nrow=Ncod,ncol = Nstage)
Stagefilter[2,1] = 0; Stagefilter[2,4:6] = 0 

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
                  Syear=Syear,tauS=tauS,pK=pK,
                  Stagefilter=Stagefilter)
#
inits <- function() list(sigQ=runif(1, .1, .5),
                         sigA=runif(1, .1, .4),sigT=runif(1, .02, .08))
params <- c("sigQ","sigA","sigT","phi","B","Q","Zmn","Z") # 
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
plot(out, vars = "phi", plot.type = c("trace", "histogram"),
     layout = c(1,2))
# plot(out, vars = "B", plot.type = c("trace", "histogram"),
#      layout = c(1,2))
#   
Nrandsamp = 5000
randsamp = sample(Nsims,Nrandsamp)
Femmort = matrix(nrow = Nrandsamp,ncol = Ncod)
for (i in 1:Ncod){
  tmp1=as.matrix(mcmc)[randsamp,which(vn==paste0("B[",i,",",2,"]"))]
  tmp2=as.matrix(mcmc)[randsamp,which(vn==paste0("phi[",i,"]"))]
  Femmort[,i] = exp(tmp1+tmp2*.90)
}
boxplot(Femmort,ylab="Log Hazard ratio",
        main = "Relative Impacts of COD, Adult Females, 90% K",
        xlab="",xaxt="n",col = "lightgray",outline=FALSE)
axis(1, labels = FALSE)
text(x =  seq_along(CODdefs$COD_Group), y = par("usr")[3] - .02, srt = 45, adj = 1,
     labels = CODdefs$COD_Group, xpd = TRUE)

rm(Femmort)
save.image(file="../Results/FitCOD_Results_2.rdata")
