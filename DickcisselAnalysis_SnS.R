####Packages####
#install.packages(c("tidyverse", "glmmTMB","unmarked","ggeffects","aiccmodavg","check_overdispersion"))
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(AICcmodavg)
library(unmarked)

####Importing####
setwd('/cloud/project/Data')
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


####..............................####
####1. PASTURE SCALE####
####..............................####
####1a. Male abundance####

abundance_by_visit=DICK_data[,4:8]
siteCovs=data.frame(list(Past_Pat_Year=DICK_data[,1], Year=DICK_data[,2], Pasture=DICK_data[,3], Area=DICK_data[,34], FireTreat=DICK_data[,35], TSF=DICK_data[,36],TSH=DICK_data[,37],HerbTreat=DICK_data[,38],GrazingYesNo=DICK_data[,39],HerbYesNo=DICK_data[,40],GrazingTreat=DICK_data[,41]))
obsCovs=list(wind=DICK_data[,24:28],DOY=DICK_data[,14:18],clouds=DICK_data[,9:13],starttime=DICK_data[,19:23],obs=DICK_data[,29:33])
DICK_PCount=unmarkedFramePCount(abundance_by_visit,
                                siteCovs,
                                obsCovs)

DICK_Abund_Global=pcount(~obs+wind+clouds+DOY+starttime ~HerbTreat_14.18+GrazingYesNo+offset(log(Area_ha)), data=DICK_PCount, K=100)
summary(DICK_Abund_Global)

exp(coef(DICK_Abund_Global, type = "state"))
exp(confint(DICK_Abund_Global, type = "state", level = 0.90))

DICK_Abund_Final=pcount(~obs+wind+clouds ~HerbTreat_14.18+GrazingYesNo+offset(log(Area_ha)), data=DICK_PCount, K=100)
summary(DICK_Abund_Final)

#predict abundance over newdata values
DICK_Abund_newdata = data.frame(GrazingYesNo = c("Yes", "Yes", "Yes","No", "No","No"),
                     HerbTreat_14.18 = c("Con","SnS","Spr"),
                     Area_ha=1)
print(DICK_Abund_newdata)

predict(DICK_Abund_Final, type = "state", newdata = DICK_Abund_newdata, appendData = T, level=0.9)

abundance_estimates = as.data.frame(predict(DICK_Abund_Final, type = "state", newdata = DICK_Abund_newdata, appendData = T, level=0.9))


####1b. Total Nests####
TotNestsGlobal  = glmmTMB (TotalNests ~  HerbTreat + GrazingYesNo + offset(log(Patchsize_ha))+ (1|Pasture), REML="FALSE", family=poisson, data=PatchData) 
summary(TotNestsGlobal)

c_hat(TotNestsGlobal, method = "pearson") #overdispersion <4 (2.66)

TotNestsFinal  = glmmTMB (TotalNests ~ HerbTreat + GrazingYesNo + +offset(log(Patchsize_ha))+ (1|Pasture), REML="FALSE", family=poisson, data=PatchData) 
summary(TotNestsFinal)


####1c. Fledgling Production per patch####
TotFledgedGlobal  = glmmTMB (TotalDICKFledged ~  HerbTreat + GrazingYesNo + offset(log(Patchsize_ha))+ (1|Pasture), REML="FALSE", family=poisson, data=PatchData) 
summary(TotFledgedGlobal)

c_hat(TotFledgedGlobal, method = "pearson") #overdispersion <4 (2)

TotNestsFinal  = glmmTMB (TotalNests ~ HerbTreat + GrazingYesNo +offset(log(Patchsize_ha))+ (1|Pasture), REML="FALSE", family=poisson, data=PatchData) 
summary(TotNestsFinal)

####..............................####
####2. NEST SCALE####
####..............................####

#####2a. Nest survival####
#analysis conducted in SAS using logistic exposure method


#####2b. Parasitism Rates####
ParaGlobal  = glmmTMB (Parasitized ~ InitiationJulian + SprayTreat  + Grazed + (1|Pasture), REML="FALSE", family="binomial",   data=Parasitism) # binomial
summary(ParaGlobal)

ParaFinal= glmmTMB (Parasitized ~  (1|Pasture), REML="FALSE", family="binomial",   data=Parasitism) # binomial
summary(ParaFinal)

#####2c. Parasitism Intensity ####
ParaIntGlobal  = glmmTMB (ParasitismIntensity ~ InitiationJulian + SprayTreat  + Grazed + (1|Pasture), REML="FALSE", family=poisson,   data=Parasitism) # 
summary(ParaIntGlobal)
c_hat(ParaIntGlobal)

ParaIntFinal= glmmTMB(ParasitismIntensity ~ InitiationJulian + SprayTreat   + (1|Pasture), REML="FALSE", family=poisson,   data=Parasitism) # binomial # 
summary(ParaIntFinal)


#####2d. Clutch Size####
ClutchGlobal  = glmmTMB (ClutchSizeAdjust ~ InitiationJulian + HerbTreat  + GrazingYesNo, REML="FALSE", family=poisson,   data=Clutch) # 

summary(ClutchGlobal)
c_hat(ClutchGlobal)

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

AranFinall  = glmmTMB (Aran~  GrazingYesNo + (1|NestID),REML="FALSE", family=nbinom2, data=ArthAbund) #negative binomial
summary(AranFinal)


#####2i. Lepidopteran Larvae Abundance####
Lepid_CatGlobal  = glmmTMB (Lepid_Cat~OrdinalDate+Method+HerbTreat_1418 + GrazingYesNo,REML="FALSE", family=poisson,data=ArthAbund) 
summary(Lepid_CatGlobal)
c_hat(Lepid_CatGlobal)

Lepid_CatFinal  = glmmTMB (Lepid_Cat~ OrdinalDate + Method + GrazingYesNo,REML="FALSE", family=poisson,data=ArthAbund) 
summary(Lepid_CatFinal)

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

#####2l. Lepi Biomass####
LepidSizeGlobal  = glmmTMB (log(Biomass_mg)~OrdinalDate+HerbTreat_1418 + GrazingYesNo +(1|NestID),REML="FALSE", family=gaussian,data=LepiSize) 
summary(LepidSizeGlobal)

LepidSizeFinal  = glmmTMB (log(Biomass_mg)~HerbTreat_1418  +(1|NestID),REML="FALSE", family=gaussian,data=LepiSize) 
summary(LepidSizeFinal)


####..............................####
####3. Graphing/layout####
####..............................####
