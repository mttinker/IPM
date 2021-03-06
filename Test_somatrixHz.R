# Explore SOmatrixHz at various proportion K and env. stochasticity 
require(ggplot2)
source("somatrixHz.R")
#
dns = c(.01,seq(.05,1.05,by = .05))
Ndns = length(dns)
AFmS = numeric()
JFmS = numeric()
PupS = numeric()
lams = numeric()
SSD = matrix(nrow=6,ncol=Ndns)
sig = .35
for (d in 1:Ndns){
  PrpK = dns[d]
  rslt = somatrixHz(PrpK,0) # no stochasticity 
  lams[d] = rslt$lam
  vrates = rslt$vrts
  SSD[,d] = rslt$SSD
  PupS[d] = vrates[1]
  JFmS[d] = vrates[2]
  AFmS[d] = vrates[3]
}
df = data.frame(Density = dns, Lambda = lams, Survive_AF = AFmS,
                Survive_JF = JFmS,WeanRate = PupS)

plt1 = (ggplot(df, aes(x= Density))+
          geom_line(aes(y=Survive_AF, colour = "Survive_AF"),size=1) +
          geom_line(aes(y=Survive_JF, colour = "Survive_JF"),size=1) +
          geom_line(aes(y=WeanRate, colour = "WeanRate"),size=1) +
          xlab("Population Density (Proportion of K)") +
          ylab("Vital Rate") +
          ggtitle("Density-dependent variation in vital rates, proportional hazards model")+
          scale_colour_discrete(name="Vital rates: ",
                                breaks=c("Survive_AF", "Survive_JF", "WeanRate"),
                                labels=c("Adult female survival", "Subadult female survival", 
                                         "Wean success rate")) +
          theme(legend.position="bottom"))
print(plt1)

plt2 = (ggplot(df, aes(x= Density,y=Lambda))+
          geom_line(size=1) +
          xlab("Population Density (Proportion of K)") +
          ylab("Annual Growth Rate (Lambda)") +
          ggtitle("Density-dependent variation in population growth, proportional hazards model"))  
print(plt2)

print("Linear function of Lambda ~ Adult female survival:")
print(summary(lm(lams ~ AFmS)))

if (sig>0){
lamvar = numeric()
vrates = matrix(nrow=1000,ncol = 7)
for (i in 1:10000){
  rslt = somatrixHz(1,sig) # with stochasticity 
  lamvar[i] = rslt$lam
  vrates[i,] = as.numeric(rslt$vrts)
}
 
boxplot(lamvar,ylab="Lambda",main = "Lambda at K with Environmental Stochasticity")
labels = c("Pups","Subadult Fem","Adult Fem","Aged-adl Fem",
                                  "Subadult Male","Adult Male","Aged-adl Male")
boxplot(vrates,ylab="Annual Survival Rates",main = "Vital rates at K with Environmental Stochasticity",
        xlab="",xaxt="n",col = "lightgray")
axis(1, labels = FALSE)
text(x =  seq_along(labels), y = par("usr")[3] - .02, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
print(paste0("Geometric mean lambda = ",exp(mean(log(lamvar)))))
}