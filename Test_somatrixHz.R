# Explore SOmatrixHz at various proportion K and env. stochasticity 
require(ggplot2)
source("somatrixHz.R")
#
dns = seq(.05,1.05,by = .05)
Ndns = length(dns)
AFmS = numeric()
JFmS = numeric()
PupS = numeric()
lams = numeric()
SSD = matrix(nrow=6,ncol=Ndns)
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
                                labels=c("Adult female survival", "Juvenile female survival", 
                                         "Wean success rate")) +
          theme(legend.position="bottom"))
print(plt1)

plt2 = (ggplot(df, aes(x= Density,y=Lambda))+
          geom_line(size=1) +
          xlab("Population Density (Proportion of K)") +
          ylab("Annual Growth Rate (Lambda)") +
          ggtitle("Density-dependent variation in population growth, proportional hazards model"))  
print(plt2)

summary(lm(lams ~ AFmS))
