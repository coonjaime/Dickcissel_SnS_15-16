####Packages####
install.packages(c("tidyverse"))

library(c(tidyverse))

####Importing####
ArthAbund=read_csv("ArthAbund.csv")
Clutch=read_csv("ClutchSize_SnS.csv")

#Global model - Orthoptera
#model Orth = ordinaldate  method  HerbTreat_1418 GrazingYesNo 

