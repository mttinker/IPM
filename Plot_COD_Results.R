# Plot COD results
# Fit Survival rates
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(reshape2)
require(readxl)
require(gridExtra)
# Load data ------------------------------------------------------------------
# 
load("../Results/FitCOD_Results_2.rdata")
# 
post = as.matrix(mcmc)
Nrandsamp = 1000
randsamp = sample(Nsims,Nrandsamp)
post =post[randsamp,] 
Stages = factor(c("Fem_Subadult","Fem_Adult","Fem_AgedAd",
           "Male_Subadult","Male_Adult","Male_AgedAd"),
           levels = c("Fem_Subadult","Fem_Adult","Fem_AgedAd",
                      "Male_Subadult","Male_Adult","Male_AgedAd"))
CODnames = factor(CODdefs$COD_Group, levels = CODdefs$COD_Group)
# 
# Survival vs density-----------------------------------------------

df_Slg = data.frame(Proportion_K = as.numeric(pK[1:5,]), Logit_AdFemSx = as.numeric(VRmats$S2lg_m),
                Logit_SubAdFemSx = as.numeric(VRmats$S1lg_m),Logit_WeanRate = as.numeric(VRmats$Wrlg_m))
plt1 = ggplot(df_Slg, aes(x = Proportion_K, y = Logit_AdFemSx)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  ggtitle("Female Survival vs Density Relative to K")
print(plt1)
plt2 = ggplot(df_Slg, aes(x = Proportion_K, y = Logit_WeanRate)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  ggtitle("Wean Success vs Density Relative to K")
print(plt2)

# Violin plot of phi by COD, showing degree of density-dependence------------- 
tempdf = data.frame(post[,which(vn=="phi[1]")])
for (i in 2:Ncod){
  tempdf = cbind(tempdf,data.frame(post[,which(vn==paste0("phi[",i,"]"))]))
}
colnames(tempdf) = CODnames
dfphi = melt(tempdf)
rm(tempdf)
plt1 = ggplot(dfphi, aes(x = variable, y = value)) +
  geom_violin(fill = "light blue", colour = "black",
              alpha = 0.7,scale = "width",trim = TRUE) +
  scale_y_continuous(name = "Parameter Value") +
  scale_x_discrete(name = "Density-dependent variation (phi) by Cause of Death") +
  ggtitle("Posterior distributions, density-dependent variation (phi) for cause-specific hazards") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt1)

# Impacts of COD by stage ---------------------------------------------
#
tmp1 = data.frame(Stage = rep(Stages,Ncod),
                          Sex = rep(c("F","F","F","M","M","M"),Ncod),
                          COD = rep(CODnames,1,each=Nstage))
tmp2 = matrix(0, nrow = Nstage*Ncod,ncol = 3); tmp2 = as.data.frame(tmp2)
colnames(tmp2) = c("Mean","CIlo","CIhi")
CODHazRepsHD =  cbind(tmp1,tmp2)
pK = 0.95
cc = 0
for (x in 1:Ncod){
  for (i in 1:Nstage){
    cc = cc+1
    if (Stagefilter[x,i]==0){
      CODHazRepsHD$Mean[cc] = 0
      CODHazRepsHD$CIlo[cc] = 0
      CODHazRepsHD$CIhi[cc] = 0
    }else{
      ii = which(vn==paste0("B[",x,",",i,"]"))
      tmp1 = post[,ii]
      # tmp1 = rnorm(10000,sumstats[ii,4],sumstats[ii,5])
      ii = which(vn==paste0("phi[",x,"]"))
      #tmp2 = rnorm(10000,sumstats[ii,4],sumstats[ii,5])
      tmp2 = post[,ii]
      tmp3 = exp(tmp1 + tmp2*pK)
      CODHazRepsHD$Mean[cc] = mean(tmp3)
      CODHazRepsHD$CIlo[cc] = quantile(tmp3,.05)
      CODHazRepsHD$CIhi[cc] = quantile(tmp3,.95)      
    }
  }
}
CODHazRepsLD = CODHazRepsHD
pK = 0.25
cc = 0
for (x in 1:Ncod){
  for (i in 1:Nstage){
    cc = cc+1
    if (Stagefilter[x,i]==0){
      CODHazRepsLD$Mean[cc] = 0
      CODHazRepsLD$CIlo[cc] = 0
      CODHazRepsLD$CIhi[cc] = 0
    }else{
      ii = which(vn==paste0("B[",x,",",i,"]"))
      tmp1 = post[,ii]
      # tmp1 = rnorm(10000,sumstats[ii,4],sumstats[ii,5])
      ii = which(vn==paste0("phi[",x,"]"))
      #tmp2 = rnorm(10000,sumstats[ii,4],sumstats[ii,5])
      tmp2 = post[,ii]
      tmp3 = exp(tmp1 + tmp2*pK)
      CODHazRepsLD$Mean[cc] = mean(tmp3)
      CODHazRepsLD$CIlo[cc] = quantile(tmp3,.05)
      CODHazRepsLD$CIhi[cc] = quantile(tmp3,.95)      
    }
  }
}
# Females
ii = which(CODHazRepsHD$Sex=="F")
pd <- position_dodge(0.5) # move them .05 to the left and right
plt1Fhd = ggplot(CODHazRepsHD[ii,], 
               aes(x=COD, y=Mean, group=Stage, color=Stage)) + 
  geom_pointrange(aes(ymin=CIlo, ymax=CIhi),position=pd) +
  labs(x = "Cause of Death", y = "Log Hazard Ratio") + 
  ggtitle("Relative Impacts on Survival by Stage, Females, 95% K") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt1Fhd)

plt1Fld = ggplot(CODHazRepsLD[ii,], 
                 aes(x=COD, y=Mean, group=Stage, color=Stage)) + 
  geom_pointrange(aes(ymin=CIlo, ymax=CIhi),position=pd) +
  labs(x = "Cause of Death", y = "Log Hazard Ratio") + 
  ggtitle("Relative Impacts on Survival by Stage, Females, 25% K") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt1Fld)

# Males
ii = which(CODHazRepsHD$Sex=="M")
pd <- position_dodge(0.5) # move them .05 to the left and right
plt1Mhd = ggplot(CODHazRepsHD[ii,], 
                 aes(x=COD, y=Mean, group=Stage, color=Stage)) + 
  geom_pointrange(aes(ymin=CIlo, ymax=CIhi),position=pd) +
  labs(x = "Cause of Death", y = "Log Hazard Ratio") + 
  ggtitle("Relative Impacts on Survival by Stage, Males, 95% K") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt1Mhd)

plt1Mld = ggplot(CODHazRepsLD[ii,], 
                 aes(x=COD, y=Mean, group=Stage, color=Stage)) + 
  geom_pointrange(aes(ymin=CIlo, ymax=CIhi),position=pd) +
  labs(x = "Cause of Death", y = "Log Hazard Ratio") + 
  ggtitle("Relative Impacts on Survival by Stage, Males, 25% K") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt1Mld)

# Variance-covariance in relative frequency ------------------------------

CODfrq = matrix(0,nrow = Nareas*Nyrs,ncol = Ncod)
for (a in 1:Nareas){
  jj = seq((a-1)*Nyrs+1,(a-1)*Nyrs+Nyrs)
  for (x in 1:Ncod){
    ii = which(startsWith(vn,paste0("Z[",x,",",a)))
    CODfrq[jj,x] = sumstats[ii,4]
  }
}

CODfrq = as.data.frame(CODfrq)
colnames(CODfrq) = CODnames

splom(CODfrq)
covmat = cov(CODfrq)
cormat = cor(CODfrq)
