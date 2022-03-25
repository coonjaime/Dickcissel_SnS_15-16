####..............................####
#METADATA----
####..............................####

#Attributions
#>> Code accompanies the following in review manuscript:
#Coon, J.J., S.B. Maresh Nelson, R.C. Daughtridge, and J.R. Miller
#Title....
#>> Code was written by J. Coon and can be found at: https://github.com/coonjaime/Dickcissel_SnS_16

#Packages----
#install.packages("easypackages")
#remotes::install_github("pcdjohnson/GLMMmisc")

library(easypackages)

packages("remotes","tidyverse", "ggeffects", "unmarked", "patchwork","lme4","MASS")
#"AICcmodavg",

#Importing & Setup----
#setwd('/cloud/project/Data')
setwd("~/Dropbox/_Manuscripts/Dissertation/Ch. 2/Dickcissel_SnS_16/Data")

ArthAbund  = read_csv("ArthAbund.csv")
Clutch     = read_csv("ClutchSize_SnS.csv")
DICK_data  = read_csv("DICK_RData2.csv")
Mass       = read_csv("NestlingMass_SnS.csv")
Parasitism = read_csv("Parasitism_SprayPasturesOnly_EarlyDeathsRemoved_6.8.18.csv")
PatchData  = read_csv("PatchNestData.csv") #note to self: this was dataset #10
Prov       = read_csv("Provisioning_byClip.csv")
OrthSize   = read_csv("DICKMeasuring_SnS_Orth.csv")
AranSize   = read_csv("DICKMeasuring_SnS_Aran.csv")
LepiSize   = read_csv("DICKMeasuring_SnS_Lepi.csv")
Survival   = read_csv("NestMonitoring3.csv")

dodge <- position_dodge(width=0.9)

#ggplot themes----
theme_bar_SnS_noleg <- function () { 
  theme(text=element_text(size=10),
        axis.title=element_text(face="bold", size=10),
        axis.text=element_text(size=10,color="black"),
        axis.line=element_line(color="black",size=1),
        panel.background=element_rect("snow2"),
        panel.grid=element_blank(),
        legend.position="none",
        legend.text=element_text(size=10, color="black"))}

theme_bar_SnS_leg <- function () { 
  theme(text=element_text(size=10),
        axis.title=element_text(face="bold", size=10),
        axis.text=element_text(size=10,color="black"),
        axis.line=element_line(color="black",size=1),
        panel.background=element_rect("snow2"),
        panel.grid=element_blank(),
        legend.position="right",
        legend.text=element_text(size=10, color="black"))}

#functions----
Coefficients =function(model) {
  coefs=c(coef(model, type = "state"))
  coefs=data.frame(coefs)
} #extracting coefficients from final models
ConfidenceIntervals = function(model,coefs) {
  CIs = as.data.frame(confint(model, type = "state", level = 0.85))
  colnames(CIs)=c("LCL", "UCL") #renames columns
  CIs
  coefs$LCL=CIs$LCL
  coefs$UCL=CIs$UCL
  coefs
}
####..............................####
#1. PASTURE SCALE----
####..............................####
####1a. Male abundance####

abundance_by_visit=DICK_data[,4:8]
siteCovs=data.frame(list(Past_Pat_Year=DICK_data[,1], Year=DICK_data[,2], Pasture=DICK_data[,3], Area=DICK_data[,34], FireTreat=DICK_data[,35], TSF=DICK_data[,36],TSH=DICK_data[,37],HerbTreat=DICK_data[,38],GrazingYesNo=DICK_data[,39],HerbYesNo=DICK_data[,40],GrazingTreat=DICK_data[,41]))
obsCovs=list(wind=DICK_data[,24:28],DOY=DICK_data[,14:18],clouds=DICK_data[,9:13],starttime=DICK_data[,19:23],obs=DICK_data[,29:33])
DICK_PCount=unmarkedFramePCount(abundance_by_visit,
                                siteCovs,
                                obsCovs)

DICK_Abund_Global=pcount(~obs+wind+clouds+DOY+starttime ~HerbTreat_14.18+GrazingYesNo+offset(log(Area_ha)), mixture="P",data=DICK_PCount, K=100)
summary(DICK_Abund_Global)

exp(coef(DICK_Abund_Global, type = "state"))
exp(confint(DICK_Abund_Global, type = "state", level = 0.90))

DICK_Abund_Final=pcount(~obs+wind+clouds ~HerbTreat_14.18+GrazingYesNo+offset(log(Area_ha)), mixture="P",data=DICK_PCount, K=100)
summary(DICK_Abund_Final)

#predict abundance 
PredictedValues=function(model,newdata) {
  abundance_estimates = as.data.frame(predict(model, type = "state", newdata = newdata, level=0.85,appendData = T))
  abundance_estimates} 

DICK_Abund_newdata = data.frame(GrazingYesNo = c("Yes", "Yes", "Yes","No", "No","No"),
                                HerbTreat_14.18 = c("Con","SnS","Spr"),
                                Area_ha=1)
#print(DICK_Abund_newdata)
DICK_Abund_Predicted=PredictedValues(DICK_Abund_Final,DICK_Abund_newdata)


#reordering factors
DICK_Abund_Predicted$HerbTreat_14.18=factor(DICK_Abund_Predicted$HerbTreat_14.18,levels=c("Con","Spr","SnS"))

#plot
DICK_Abund_Plot=ggplot(data=DICK_Abund_Predicted, y=Predicted, x=GrazingYesNo)+  
  geom_bar(aes(x=GrazingYesNo, y=Predicted,fill=HerbTreat_14.18), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_x_discrete(labels=c("Ungrazed","Grazed"))+
  scale_y_continuous(limits=c(0,11.2),expand = c(0, 0)) +
  geom_errorbar(aes(x = GrazingYesNo, ymin = lower, ymax = upper, group=HerbTreat_14.18),position = dodge, width = 0.2)+
  labs(y = "Dickcissel Males per Ha", x=" ")+
  theme_bar_SnS_noleg()
DICK_Abund_Plot


####1b. Total Nests####
TotNestsGlobal  = glmmTMB (TotalNests ~  HerbTreat + GrazingYesNo + offset(log(Patchsize_ha))+ (1|Pasture), family=poisson, data=PatchData, REML=FALSE) 
summary(TotNestsGlobal)

#Final model
TotNestsFinal  = glmmTMB (TotalNests ~ HerbTreat + GrazingYesNo + as.factor(offset(log(Patchsize_ha)))+ (1|Pasture), REML="FALSE", family=nbinom2, data=PatchData) 
summary(TotNestsFinal)

#producing predicted values
TotNests_Pred = as.data.frame(ggpredict(TotNestsFinal,c("HerbTreat", "GrazingYesNo"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) 
colnames(TotNests_Pred)=c("HerbTreat", "Predicted","SE","Lower","Upper", "Grazing") #renames columns
print(TotNests_Pred) 
confint(TotNestsFinal, level = 0.85)

c_hat(TotNestsFinal, method = "pearson") #overdispersion >2 (2.66)

#reordering factors
TotNests_Pred$HerbTreat=factor(TotNests_Pred$HerbTreat,levels=c("Con","Spr","SnS"),labels=c("Control","Spray","Spray & Seed"))

#plot
TotNests_Plot=ggplot(data=TotNests_Pred, y=Predicted, x=Grazing)+  
  geom_bar(aes(x=Grazing, y=Predicted,fill=HerbTreat), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_x_discrete(labels=c("Ungrazed","Grazed"))+
  scale_y_continuous(limits=c(0,40),expand = c(0, 0)) +
  geom_errorbar(aes(x = Grazing, ymin = Lower, ymax = Upper, group=HerbTreat),position = dodge, width = 0.2)+
  labs(y = "Total Nests per Patch", x=" ")+
  theme_bar_SnS_leg()+
  theme(
    legend.title=element_blank())
TotNests_Plot

####1c. Fledgling Production per patch####
TotFledgedGlobal  = glmmTMB (TotalDICKFledged ~  HerbTreat + GrazingYesNo + offset(log(Patchsize_ha))+ (1|Pasture), REML="FALSE", family=poisson, data=PatchData) 
summary(TotFledgedGlobal)

c_hat(TotFledgedGlobal, method = "pearson") #overdispersion <2 (2)

TotFledgedFinal  = glmmTMB (TotalDICKFledged ~ HerbTreat +(TotalNests)+ (1|Pasture), REML="FALSE", family=poisson, data=PatchData) 
summary(TotFledgedFinal)


#producing predicted values
TotFledged_Pred = as.data.frame(ggpredict(TotFledgedFinal,c("HerbTreat"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) 
colnames(TotFledged_Pred)=c("HerbTreat", "Predicted","SE","Lower","Upper") #renames columns
print(TotFledged_Pred) 
confint(TotFledgedFinal, level = 0.85)

#reordering factors
TotFledged_Pred$HerbTreat=factor(TotFledged_Pred$HerbTreat,levels=c("Con","Spr","SnS"))

#plot
TotFledged_Plot=ggplot(data=TotFledged_Pred, y=Predicted, x=HerbTreat)+  
  geom_bar(aes(x=HerbTreat, y=Predicted,fill=HerbTreat), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_x_discrete(labels=c("Control","Spray","Spray & Seed"))+
  scale_y_continuous(breaks=c(2,4,6,8),limits=c(0,8),expand = c(0, 0)) +
  geom_errorbar(aes(x = HerbTreat, ymin = Lower, ymax = Upper),position = dodge, width = 0.2)+
  labs(y = "Fledged per Patch", x="Herbicide Treatment")+
  theme_bar_SnS_noleg()

TotFledged_Plot
####..............................####
####2. NEST SCALE####
####..............................####

#####2a. Nest survival####
#analysis conducted in SAS using logistic exposure method

#some data manipulaton to add herbicide treatment to this data, merged from the Clutch dataset

#logistic exposure function
logexp <- function(exposure = 1) {
  get_exposure <- function() {
    if (exists("..exposure", env=.GlobalEnv))
      return(get("..exposure", envir=.GlobalEnv))
    exposure
  }
  linkfun <- function(mu) qlogis(mu^(1/get_exposure()))
  linkinv <- function(eta) plogis(eta)^get_exposure()
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta)>30,.Machine$double.eps,
           exp(eta)/(1+exp(eta))^2)
  }
  mu.eta <- function(eta) {       
    get_exposure() * plogis(eta)^(get_exposure()-1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

Survive_Global = glm (Survive ~ HerbTreat + GrazingYesNo + JulianDate + Stage + Pasture, data=Survival,family=binomial(link=logexp(exposure=Survival$ExpDays))) # binomial
summary(Survive_Global)

Survive_Final = glm (Survive ~ HerbTreat  + JulianDate + Pasture, data=Survival,family=binomial(link=logexp(exposure=Survival$ExpDays))) # binomial
summary(Survive_Final)
#https://rpubs.com/bbolker/logregexp
#also https://www.perrywilliams.us/wp-content/uploads/2018/03/Crimmins2016factors.pdf

summary(Survival$HerbTreat)

#producing predicted values
..exposure <- 1
Survive_Pred = as.data.frame(ggpredict(Survive_Final,terms=c("HerbTreat","Pasture"),
                                       ci.lvl=0.85, 
                                       back.transform=TRUE, 
                                       append=TRUE)) #turns predictions into a dataframe that we can more easily manipulate
colnames(Survive_Pred)=c("HerbTreat", "Predicted","SE","Lower","Upper","Pasture") #renames columns
print(Survive_Pred) 
rm(..exposure)
confint(Survive_Final, level = 0.85)

Survive_Pred_Sum =Survive_Pred%>%
  group_by(HerbTreat)%>%
  summarise_at(vars(Predicted,SE,Lower,Upper),mean)


#plotg
Survive_Plot=ggplot(data=Survive_Pred_Sum, y=Predicted, x=HerbTreat)+  
  geom_bar(aes(x=HerbTreat, y=Predicted,fill=HerbTreat), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  #scale_x_discrete(labels=c("Control","Spray","Spray & Seed"))+
  scale_y_continuous(limits=c(0,1),expand = c(0, 0)) +
  geom_errorbar(aes(x = HerbTreat, ymin = Lower, ymax = Upper,group=HerbTreat),position = dodge, width = 0.2)+
  labs(y = "Daily Nest Survival",x="Herbicide Treatment")+
  scale_x_discrete(labels=c("Control","Spray","Spray & Seed"))+
  
  theme_bar_SnS_noleg()
Survive_Plot

#####2b. Parasitism Rates####
ParaGlobal  = glmmTMB (Parasitized ~ InitiationJulian + SprayTreat  + Grazed + (1|Pasture), REML="FALSE", family="binomial",   data=Parasitism) # binomial
summary(ParaGlobal)

ParaFinal= glmmTMB (Parasitized ~  InitiationJulian +(1|Pasture), REML="FALSE", family="binomial",   data=Parasitism) # binomial
summary(ParaFinal)

#no plot

#####2c. Parasitism Intensity ####
Parasitism$SprayTreat=recode(Parasitism$SprayTreat,"None"="Con")
Parasitism$SprayTreat=recode(Parasitism$SprayTreat,"Spray"="Spr")
Parasitism$SprayTreat=recode(Parasitism$SprayTreat,"SpraySeed"="SnS")

summary(as.factor(Parasitism$SprayTreat))

#global model
ParaIntGlobal  = glmmTMB (ParasitismIntensity ~ InitiationJulian + SprayTreat  + Grazed + (1|Pasture), REML="FALSE", family=poisson,   data=Parasitism) # 
summary(ParaIntGlobal)
c_hat(ParaIntGlobal) #1.33
#final model
ParaIntFinal= glmmTMB(ParasitismIntensity ~ InitiationJulian + SprayTreat   + (1|Pasture), REML="FALSE", family=poisson,   data=Parasitism) # binomial # 
summary(ParaIntFinal)

#producing predicted values
ParaInt_Pred = as.data.frame(ggpredict(ParaIntFinal,c("SprayTreat"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) #turns predictions into a dataframe that we can more easily manipulate
colnames(ParaInt_Pred)=c("HerbTreat", "Predicted","SE","Lower","Upper") #renames columns
View(ParaInt_Pred) 
confint(ParaIntFinal, level = 0.85)

#reordering factors
ParaInt_Pred$HerbTreat=factor(ParaInt_Pred$HerbTreat,levels=c("Con","Spr","SnS"),labels=c("Control","Spray","Spray & Seed"))

#plot
ParaInt_Plot=ggplot(data=ParaInt_Pred, y=Predicted, x=HerbTreat)+  
  geom_bar(aes(x=HerbTreat, y=Predicted,fill=HerbTreat), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_x_discrete(labels=c("Control","Spray","Spray & Seed"))+
  scale_y_continuous(limits=c(0,2.4),expand = c(0, 0)) +
  geom_errorbar(aes(x = HerbTreat, ymin = Lower, ymax = Upper),position = dodge, width = 0.2)+
  labs(y = "Parasitism Intensity (#BHCO)", x="Herbicide Treatment")+
  theme_bar_SnS_noleg()
ParaInt_Plot

#####2d. Clutch Size####
ClutchGlobal  = glmmTMB (ClutchSizeAdjust ~ InitiationJulian + HerbTreat  + GrazingYesNo, REML="FALSE", family=poisson, data=Clutch) # 
summary(ClutchGlobal)
c_hat(ClutchGlobal) #.22 - underdispersed

ClutchFinal  = glmmTMB (ClutchSizeAdjust ~ 1, REML="FALSE", family=poisson,   data=Clutch) # 
summary(ClutchFinal)

#####2e. Nestling Mass####

MassGlobal    = glmmTMB (NestlingMass ~ OrdinalMeasured  + TimeOfDay  + NestlingAge + HerbTreat  + GrazingYesNo + (1|NestID), REML="FALSE", family=gaussian, data=Mass) 
summary(MassGlobal)

MassFinal = glmmTMB (NestlingMass ~ NestlingAge + (1|NestID), REML="FALSE", family=gaussian, data=Mass) 
summary(MassFinal)


#####2f. Provisioning Rates per hr####

ProvGlobal    = glmmTMB (InstancesPerHr ~ NestlingAge_Days + TotalNestlingNum + HerbTreat +GrazingYesNo +(1|NestID), REML="FALSE", family="tweedie", data=Prov)  
summary(ProvGlobal)

ProvFinal    = glmmTMB (InstancesPerHr ~ NestlingAge_Days + TotalNestlingNum + (1|NestID), REML="FALSE", family="tweedie", data=Prov)  
summary(ProvFinal)


#####2g. Orthoptera Abundance####
OrthGlobal  = glmmTMB (Orth ~ OrdinalDate + Method + HerbTreat_1418 + GrazingYesNo + (1|NestID), REML="FALSE", family=nbinom2,   data=ArthAbund) #negative binomial
summary(OrthGlobal)
#c_hat(OrthGlobal, method = "pearson") overdispersion >4

OrthFinal= glmmTMB (Orth ~ OrdinalDate + Method + (1|NestID), REML="FALSE", family=nbinom2,   data=ArthAbund) #negative binomial
summary(OrthFinal)

#####2h. Araneae Abundance####
AranGlobal  = glmmTMB (Aran~OrdinalDate+Method+HerbTreat_1418 + GrazingYesNo + (1|NestID),REML="FALSE", family=nbinom2, data=ArthAbund) #negative binomial
summary(AranGlobal)

AranFinal  = glmmTMB (Aran~  GrazingYesNo + (1|NestID),REML="FALSE", family=nbinom2, data=ArthAbund) #negative binomial
summary(AranFinal)

#producing predicted values
Aran_Pred = as.data.frame(ggpredict(AranFinal,c("GrazingYesNo"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) #turns predictions into a dataframe that we can more easily manipulate
colnames(Aran_Pred)=c("Grazing", "Predicted","SE","Lower","Upper") #renames columns
print(Aran_Pred) 
confint(AranFinal, level = 0.85)

#plot
Aran_Plot=ggplot(data=Aran_Pred, y=Predicted, x=Grazing)+  
  geom_bar(aes(x=Grazing, y=Predicted), position=dodge, stat="identity")+
  scale_x_discrete(labels=c("Ungrazed","Grazed"))+
 scale_y_continuous(limits=c(0,20.2),expand = c(0, 0)) +
  geom_errorbar(aes(x = Grazing, ymin = Lower, ymax = Upper),position = dodge, width = 0.2)+
  labs(y = "Spider per Sample", x="Grazing Status")+
  theme_bar_SnS_noleg()
Aran_Plot

#####2i. Lepidopteran Larvae Abundance####
Lepid_CatGlobal  = glmmTMB (Lepid_Cat~OrdinalDate+Method+HerbTreat_1418 + GrazingYesNo,REML="FALSE", family=poisson,data=ArthAbund) 
summary(Lepid_CatGlobal)
c_hat(Lepid_CatGlobal)

Lepid_CatFinal  = glmmTMB (Lepid_Cat~ OrdinalDate + Method + GrazingYesNo,REML="FALSE", family=poisson,data=ArthAbund) 
summary(Lepid_CatFinal)

#producing predicted values
Lepid_Cat_Pred = as.data.frame(ggpredict(Lepid_CatFinal,c("GrazingYesNo"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) #turns predictions into a dataframe that we can more easily manipulate
colnames(Lepid_Cat_Pred)=c("Grazing", "Predicted","SE","Lower","Upper") #renames columns
print(Lepid_Cat_Pred) 
confint(Lepid_CatFinal, level = 0.85)

#plot
Lepid_Cat_Plot=ggplot(data=Lepid_Cat_Pred, y=Predicted, x=Grazing)+  
  geom_bar(aes(x=Grazing, y=Predicted), position=dodge, stat="identity")+
  scale_x_discrete(labels=c("Ungrazed","Grazed"))+
  scale_y_continuous(limits=c(0,5.2),expand = c(0, 0)) +
  geom_errorbar(aes(x = Grazing, ymin = Lower, ymax = Upper),position = dodge, width = 0.2)+
  labs(y = "Caterpillars per Sample", x=" ")+
  theme_bar_SnS_noleg()+
  theme(
    axis.text.x=element_blank()
  )
Lepid_Cat_Plot

#####2j. Orth Biomass####
OrthSizeGlobal  = glmmTMB (log(Biomass_mg)~OrdinalDate+HerbTreat_1418 + GrazingYesNo +(1|NestID),REML="FALSE", family=gaussian,data=OrthSize) 
summary(OrthSizeGlobal)

OrthSizeFinal  = glmmTMB (log(Biomass_mg)~OrdinalDate +(1|NestID),REML="FALSE", family=gaussian,data=OrthSize) 
summary(OrthSizeFinal)

#####2k. Aran Biomass####
AranSizeGlobal  = glmmTMB (log(Biomass_mg)~OrdinalDate+HerbTreat_1418 + GrazingYesNo +(1|NestID),REML="FALSE", family=gaussian,data=AranSize) 
summary(AranSizeGlobal)

AranSizeFinal  = glmmTMB (log(Biomass_mg)~HerbTreat_1418 +(1|NestID),REML="FALSE", family=gaussian,data=AranSize) 
summary(AranSizeFinal)

#producing predicted values
AranSize_Pred = as.data.frame(ggpredict(AranSizeFinal,c("HerbTreat_1418"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) #turns predictions into a dataframe that we can more easily manipulate
colnames(AranSize_Pred)=c("HerbTreat", "Predicted","SE","Lower","Upper") #renames columns
print(AranSize_Pred) 
confint(ParaIntFinal, level = 0.85)

#reordering factors
AranSize_Pred$HerbTreat=factor(AranSize_Pred$HerbTreat,levels=c("Con","Spr","SnS"))

#plot
AranSize_Plot=ggplot(data=AranSize_Pred, y=Predicted, x=HerbTreat)+  
  geom_bar(aes(x=HerbTreat, y=Predicted,fill=HerbTreat), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_x_discrete(breaks=c("Con","Spr","SnS"),labels=c("Control","Spray","Spray & Seed"))+
  scale_y_continuous(limits=c(0,11),expand = c(0, 0)) +
  geom_errorbar(aes(x = HerbTreat, ymin = Lower, ymax = Upper),position = dodge, width = 0.2)+
 labs(y = "Spider Size (mg)", x="Herbicide Treatment")+
  theme_bar_SnS_noleg()+
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  )
AranSize_Plot

#####2l. Lepi Biomass####
LepidSizeGlobal  = glmmTMB (log(Biomass_mg)~OrdinalDate+HerbTreat_1418 + GrazingYesNo +(1|NestID),REML="FALSE", family=gaussian,data=LepiSize) 
summary(LepidSizeGlobal)

LepidSizeFinal  = glmmTMB (log(Biomass_mg)~HerbTreat_1418  +(1|NestID),REML="FALSE", family=gaussian,data=LepiSize) 
summary(LepidSizeFinal)

#producing predicted values
LepidSize_Pred = as.data.frame(ggpredict(LepidSizeFinal,c("HerbTreat_1418"),ci.lvl=0.85, back.transform=TRUE, append=TRUE)) #turns predictions into a dataframe that we can more easily manipulate
colnames(LepidSize_Pred)=c("HerbTreat", "Predicted","SE","Lower","Upper") #renames columns
print(LepidSize_Pred) 
confint(ParaIntFinal, level = 0.85)

#reordering factors
LepidSize_Pred$HerbTreat=factor(LepidSize_Pred$HerbTreat,levels=c("Con","Spr","SnS"))

#plot
LepidSize_Plot=ggplot(data=LepidSize_Pred, y=Predicted, x=HerbTreat)+  
  geom_bar(aes(x=HerbTreat, y=Predicted,fill=HerbTreat), position=dodge, stat="identity")+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_x_discrete(breaks=c("Con","Spr","SnS"),labels=c("Control","Spray","Spray & Seed"))+
  scale_y_continuous(limits=c(0,13),expand = c(0, 0)) +
  geom_errorbar(aes(x = HerbTreat, ymin = Lower, ymax = Upper),position = dodge, width = 0.2)+
  labs(y = "Caterpillar Size (mg)", x=" ")+
  theme_bar_SnS_noleg()+
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  )
LepidSize_Plot

####..............................####
####3. Graphing/layout####
####..............................####
print(PatchData)
PatchData$HerbTreat=factor(PatchData$HerbTreat,levels=c("Con","Spr","SnS"))

#plot sample size
SampSize   =PatchData%>%
  #filter(GrazingYesNo == "Yes")%>%
  select(Pasture_Patch_Year,TotalNests,SuccessfulNests,HerbTreat,GrazingYesNo)%>%
  group_by(HerbTreat,GrazingYesNo) %>%
  summarise_at(vars(TotalNests, SuccessfulNests), sum)%>%
  gather(Nests,Number,TotalNests:SuccessfulNests)

SampSize$Nests <- factor(SampSize$Nests, 
                               levels = c("SuccessfulNests", "TotalNests") )

SampSize$GrazingYesNo <- factor(SampSize$GrazingYesNo, levels = c("No", "Yes"),
                  labels = c("Ungrazed", "Grazed"))

SampSize$HerbTreat <- factor(SampSize$HerbTreat, levels = c("Con", "Spr","SnS"),
                                labels = c("Control", "Spray", "Spray & Seed"))

Fig3 = ggplot(data=SampSize,aes(x=HerbTreat, y=Number, group=Nests))+  
  geom_bar(aes(fill=HerbTreat,alpha=factor(Nests)),  position = "identity", stat="identity")+
  geom_text(aes(label=Number), position="identity", vjust=-0.25,size=3)+
  facet_wrap(~GrazingYesNo, strip.position="top")+
  scale_alpha_manual(values = c("SuccessfulNests"=1, "TotalNests"=0.6), guide='none')+
  scale_fill_manual(values=c("goldenrod3","darkseagreen4","darkslategray"))+
  scale_y_continuous(limits=c(0,30),expand = c(0, 0)) +
  theme_bar_SnS_noleg()+
  labs(y = "# Nests (Total & Successful)", x="Grazing")+
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text=element_text(color="black",size=10, face="bold"),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank())
       #legend.position=c(0.85, 0.8))
Fig3

setwd('/cloud/project/Figures')

ggsave("Fig3.jpeg", Fig3,unit="in", width=6,height=4, dpi=600)

Fig4 =  
  TotNests_Plot + 
  guide_area() +
  DICK_Abund_Plot + 
  TotFledged_Plot + 
  plot_layout(guides = "collect",
              nrow=2)+ 
  plot_annotation(tag_levels = 'A')
Fig4

ggsave("Fig4.jpeg", Fig4,unit="in", width=8,height=6, dpi=600)


Fig5 =  
  AranSize_Plot +
  Aran_Plot + 
  LepidSize_Plot + 
  Lepid_Cat_Plot + 
  ParaInt_Plot+
  Survive_Plot+
  plot_layout(guides = "collect",
              nrow=3)+ 
  plot_annotation(tag_levels = 'A')
Fig5

ggsave("Fig5.jpeg", Fig5, unit="in", width=8,height=9, dpi=600)


