####Packages####
#install.packages(c("tidyverse", "glmmTMB","unmarked","ggeffects","aiccmodavg","check_overdispersion"))
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(AICcmodavg)
library(unmarked)

####Importing####
PatchData=read_csv("PatchNestData.csv") #note to self: this was dataset #10
ArthAbund=read_csv("ArthAbund.csv")
Clutch=read_csv("ClutchSize_SnS.csv")
DICKData=read_csv("DICK_RData2.csv")

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

# Selecting the best model for detection
null=pcount(~1~1, data=DICK_PCount, K=100)
model_1=pcount(~wind ~1, data=DICK_PCount, K=100)
model_2=pcount(~DOY ~1, data=DICK_PCount, K=100)
model_3=pcount(~clouds ~1, data=DICK_PCount, K=100)
model_4=pcount(~obs ~1, data=DICK_PCount, K=100)
model_5=pcount(~clouds+wind ~1, data=DICK_PCount, K=100)
model_6=pcount(~wind+obs+clouds ~1, data=DICK_PCount, K=100)


null

model_selection_det <- fitList("Null"=null,
                               "det(wind)abund(.)"=model_1,
                               "det(DOY)abund(.)"=model_2,
                               "det(clouds)abund(.)"=model_3,
                               "det(obs)abund(.)"=model_4,
                               "det(clouds+wind)abund(.)"=model_5,
                               "det(obs+wind)abund(.)"=model_6)
modSel(model_selection_det)

# Selecting the best model for abundance
null2=pcount(~obs+wind ~offset(log(Area)), data=DICK_PCount, K=100)
model_7=pcount(~obs+wind ~HerbTreat+offset(log(Area)), data=DICK_PCount, K=100)
model_8=pcount(~obs+wind ~GrazingYesNo+offset(log(Area)), data=DICK_PCount, K=100)
model_9=pcount(~obs+wind ~GrazingYesNo+HerbTreat+offset(log(Area)), data=DICK_PCount, K=100)
model_10=pcount(~obs+wind ~HerbYesNo+offset(log(Area)), data=DICK_PCount, K=100)
model_11=pcount(~obs+wind ~HerbYesNo+GrazingYesNo+offset(log(Area)), data=DICK_PCount, K=100)

model_selection_abund <- fitList("det(obs+wind)abund(.)"=null2,
                                 "det(obs+wind)abund(HerbTreat)"=model_7,
                                 "det(obs+wind)abund(GrazingYesNo)"=model_8,
                                 "det(obs+wind)abund(GrazingYesNo+HerbTreat)"=model_9,
                                 "det(obs+wind)abund(HerbYesNo)"=model_10,
                                 "det(obs+wind)abund(HerbYesNo+GrazingYesNo)"=model_11)
modSel(model_selection_abund)

summary (model_9)
LRT(model_9, model_7)
LRT(model_9,model_8)
LRT(model_7,null2)

#coefficients and confidence intervals in log scale

coef(model_9, type = "state")
backTransform(model_9,type='det')
confint(model_9, type = "state", level = 0.90)

#odds ratios (exponentiated beta coefficients)

exp(coef(model_9, type = "state"))
exp(confint(model_9, type = "state", level = 0.90))

deviance (model_9)

#predict abundance over newdata values

newdata = data.frame(GrazingYesNo = c("Yes", "Yes", "Yes","No", "No","No"),
                     HerbTreat = c("Con","SnS","Spr"),
                     Area=1)
print(newdata)

mean((predict(model_11, type = "det"))[,1])
#mean detection rate at all surveys for all observers is 0.29

predict(model_9, type = "state", newdata = newdata, appendData = T)

abundance_estimates = as.data.frame(predict(model_9, type = "state", newdata = newdata, appendData = T))

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


#####2a. Orthoptera Abundance####
OrthGlobal  = glmmTMB (Orth ~ OrdinalDate + Method + HerbTreat_1418 + GrazingYesNo + (1|NestID), REML="FALSE", family=nbinom2,   data=ArthAbund) #negative binomial
summary(OrthGlobal)
#c_hat(OrthGlobal, method = "pearson") overdispersion >4

OrthFinal= glmmTMB (Orth ~ OrdinalDate + Method + (1|NestID), REML="FALSE", family=nbinom2,   data=ArthAbund) #negative binomial
summary(OrthFinal)
#Global model - spiders (Aran)

AranGlobal  = glmmTMB (Aran~OrdinalDate+Method+HerbTreat_1418 + GrazingYesNo + (1|NestID),REML="FALSE", family=nbinom2, data=ArthAbund) #negative binomial
summary(AranGlobal)
