/*######################################################################################
ARTHROPOD ABUNDANCE
######################################################################################*/

/*importing arthropod abundance dataset*/
  PROC IMPORT OUT= WORK.ArthAbund 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\Analysis_ArthAbund\ArthAbund.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  /*global model ORTHOPTERAN ABUNDANCE*/
    
  Proc glimmix data=ArthAbund method=laplace plots=residualpanel /*orth*/;
  class NestID Past HerbYesNo_alltime HerbTreat_1418 GrazingYesNo method;
  model Orth = ordinaldate  method  HerbTreat_1418 GrazingYesNo / dist=nb solution;
  *lsmeans / ilink cl alpha=0.15 ;
  random int /subject=NestID type=vc G solution;
  
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  run;
  
  /*final model ORTHOPTERAN ABUNDANCE*/
    Proc glimmix data=ArthAbund method=laplace plots=residualpanel /*orth*/;
  class NestID Past HerbYesNo_alltime HerbTreat_1418 GrazingYesNo method;
  model Orth = ordinaldate  method / dist=nb solution;
  *lsmeans / ilink cl alpha=0.15 ;
  random int /subject=NestID type=vc G solution;
  
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  /*global model ARANEAE ABUNDANCE*/
    Proc glimmix data=ArthAbund method=laplace plots=residualpanel /*aran*/;
  class NestID Past HerbYesNo_alltime HerbTreat_1418 GrazingYesNo method;
  model Aran =  GrazingYesNo HerbTreat_1418 ordinaldate method / dist=nb solution;
  random int /subject=NestID  type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  /*final model ARANEAE ABUNDANCE*/
    Proc glimmix data=ArthAbund method=laplace plots=residualpanel /*aran*/;
  class NestID Past HerbYesNo_alltime HerbTreat_1418 GrazingYesNo method;
  model Aran =  GrazingYesNo / dist=nb solution;
  random int /subject=NestID  type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  estimate 'No grazing'	intercept 1 	GrazingYesNo 1 0 	 	/ cl 	alpha=0.1 	ilink;
  estimate 'Yes grazing'	intercept 1 	GrazingYesNo 0 1	 	/ cl 	alpha=0.1 	ilink;
  
  run;
  
  /*global model & final model LEPIDOPTERAN LARVAE ABUNDANCE*/
    Proc glimmix data=ArthAbund method=laplace plots=residualpanel /*lepi*/;
  class NestID Past HerbYesNo_alltime HerbTreat_1418 GrazingYesNo method;
  model Lepid_Cat = ordinaldate HerbTreat_1418 GrazingYesNo method / dist=p solution;
  *lsmeans / ilink cl alpha=0.15 ;
  random int/subject=NestID type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  estimate ' Con ' 	intercept 1 	HerbTreat_1418 1 0 0 	ordinaldate 195 	GrazingYesNo 0.5 0.5 	method 0.5 0.5 / 	cl 	alpha=0.1 	ilink;
  estimate ' SnS ' 	intercept 1 	HerbTreat_1418 0 1 0	ordinaldate 195 	GrazingYesNo 0.5 0.5 	method 0.5 0.5 /	cl 	alpha=0.1 	ilink;
  estimate ' Spr ' 	intercept 1 	HerbTreat_1418 0 0 1 	ordinaldate 195 	GrazingYesNo 0.5 0.5 	method 0.5 0.5 /	cl 	alpha=0.1 	ilink;
  
  run;
  
  /*######################################################################################
  CLUTCH SIZE
  ######################################################################################*/
  
  /*importing clutch size dataset*/
    PROC IMPORT OUT= WORK.Clutch 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\Analysis_ClutchSize\ClutchSize_SnS.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  
  /*global model CLUTCH SIZE*/
    Proc glimmix data=clutch method=laplace plots=residualpanel /*clutch*/;
  class NestID Var2 Pasture_Patch Pasture_patch_year HerbYesNo HerbTreat GrazingYesNo;
  model ClutchSizeAdjust = HerbTreat GrazingYesNo InitiationJulian  / dist=p solution;
  random int/subject=var2 type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  run;
  
  /*final model CLUTCH SIZE*/
    Proc glimmix data=clutch method=laplace plots=residualpanel /*clutch*/;
  class NestID Var2 Pasture_Patch Pasture_patch_year HerbYesNo HerbTreat GrazingYesNo;
  model ClutchSizeAdjust =  / dist=p solution;
  random int/subject=var2 type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  run;
  
  /*######################################################################################
  ORTHOPTERAN BIOMASS
  ######################################################################################*/
  
  /*importing ORTHOPTERAN BIOMASS dataset*/
    PROC IMPORT OUT= WORK.ArthSize 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\Analysis_NestArthropodSize_Round1\DICKMeasuring_SnS_Orth.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  
  /*global & final model ORTHOPTERAN BIOMASS*/
    proc means data=arthsize;
  var ordinaldate;
  run;
  Proc glimmix data=ArthSize method=laplace plots=residualpanel /*orthbio*/;
  class Category NestId var10 Pasture_Patch Pasture_Patch_year HerbYesNo_alltime HerbTreat_1418 GrazingYesNo;
  model Biomass_mg = HerbTreat_1418 GrazingYesNo OrdinalDate / dist=gamma solution;
  random int / subject=var10;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  estimate 'Con Grazing Yes'	intercept 1 	HerbTreat_1418 1 0 0 	GrazingYesNo 0 1 	ordinaldate 201	 	/ cl 	alpha=0.1 	ilink ;
  estimate 'Spr Grazing Yes'	intercept 1 	HerbTreat_1418 0 0 1 	GrazingYesNo 0 1   	ordinaldate 201		/ cl 	alpha=0.1	ilink;
  estimate 'SnS Grazing Yes' 	intercept 1 	HerbTreat_1418 0 1 0 	GrazingYesNo 0 1  	ordinaldate 201		/ cl 	alpha=0.1 	ilink;
  
  estimate 'Con Grazing No'	intercept 1 	HerbTreat_1418 1 0 0 	GrazingYesNo 1 0 	ordinaldate 201 	/ cl 	alpha=0.1 	ilink ;
  estimate 'Spr Grazing No'	intercept 1 	HerbTreat_1418 0 0 1 	GrazingYesNo 1 0   	ordinaldate 201		/ cl 	alpha=0.1	ilink;
  estimate 'SnS Grazing No' 	intercept 1 	HerbTreat_1418 0 1 0 	GrazingYesNo 1 0  	ordinaldate 201		/ cl 	alpha=0.1 	ilink;
  
  run;
  
  /*######################################################################################
  ARANEAE BIOMASS
  ######################################################################################*/
  /*importing ARANEAE BIOMASS dataset*/
    PROC IMPORT OUT= WORK.AranSize 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\Analysis_NestArthropodSize_Round1\DICKMeasuring_SnS_Aran.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  
  /*global model ARANEAE BIOMASS*/
    Proc glimmix data=AranSize method=laplace plots=residualpanel /*aranbio*/;
  class Category NestId var10 Pasture_Patch Pasture_Patch_year HerbYesNo_alltime HerbTreat_1418 GrazingYesNo;
  model Biomass_mg = HerbTreat_1418 GrazingYesNo ordinaldate method / dist=logn solution;
  random int / subject=var10;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  run;
  
  /*final model ARANEAE BIOMASS*/
    Proc glimmix data=AranSize method=laplace plots=residualpanel/*aranbio*/;
  class Category NestId var10 Pasture_Patch Pasture_Patch_year HerbYesNo_alltime HerbTreat_1418 GrazingYesNo;
  model Biomass_mg = HerbTreat_1418 GrazingYesNo / dist=logn solution;
  random int / subject=var10;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  estimate 'Con Grazing Yes'	intercept 1 	HerbTreat_1418 1 0 0 	GrazingYesNo 0 1 	 		/ cl 	alpha=0.1 	ilink ;
  estimate 'Spr Grazing Yes'	intercept 1 	HerbTreat_1418 0 0 1 	GrazingYesNo 0 1   	 		/ cl 	alpha=0.1	ilink;
  estimate 'SnS Grazing Yes' 	intercept 1 	HerbTreat_1418 0 1 0 	GrazingYesNo 0 1  	 		/ cl 	alpha=0.1 	ilink;
  
  estimate 'Con Grazing No'	intercept 1 	HerbTreat_1418 1 0 0 	GrazingYesNo 1 0 	 		/ cl 	alpha=0.1 	ilink ;
  estimate 'Spr Grazing No'	intercept 1 	HerbTreat_1418 0 0 1 	GrazingYesNo 1 0   	 		/ cl 	alpha=0.1	ilink;
  estimate 'SnS Grazing No' 	intercept 1 	HerbTreat_1418 0 1 0 	GrazingYesNo 1 0  	 		/ cl 	alpha=0.1 	ilink;
  
  run;
  
  /*######################################################################################
  LEPIDOPTERAN LARVAE BIOMASS
  ######################################################################################*/
  /*importing LEPIDOPTERAN LARVAE BIOMASS dataset*/
    /*LEPIDOPTERA LARVAE*/
    PROC IMPORT OUT= WORK.LepiSize 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\Analysis_NestArthropodSize_Round1\DICKMeasuring_SnS_Lepi.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  
  /*global model LEPIDOPTERAN LARVAE BIOMASS*/
    Proc glimmix data=LepiSize method=laplace plots=residualpanel /*lepibio*/;
  class Category NestId Var10 Pasture_Patch Pasture_Patch_year HerbYesNo_alltime HerbTreat_1418 GrazingYesNo;
  model Biomass_mg = HerbTreat_1418 GrazingYesNo ordinaldate   / dist=expo solution;
  random int / subject=var10 type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  
  /*final model LEPIDOPTERAN LARVAE BIOMASS*/
    Proc glimmix data=LepiSize method=laplace plots=residualpanel /*lepibio*/;
  class Category NestId Var10 Pasture_Patch Pasture_Patch_year HerbYesNo_alltime HerbTreat_1418 GrazingYesNo;
  model Biomass_mg =    / dist=expo solution;
  random int / subject=var10 type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  
  /*######################################################################################
  NESTLING MASS
  ######################################################################################*/
  /*importing NESTLING MASS BIOMASS dataset*/
    PROC IMPORT OUT= WORK.Mass 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\Analysis_NestlingMass\NestlingMass_SnS.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  
  /*global model NESTLING MASS BIOMASS*/
    Proc glimmix data=mass method=laplace plots=residualpanel /*nestlingmass*/;
  class NestID _Pasture Pasture_Patch Pasture_patch_year HerbYesNo HerbTreat GrazingYesNo;
  model Nestlingmass = HerbTreat GrazingYesNo NestlingAge OrdinalMeasured ChicksInNest TimeofDay / dist=expo solution;
  random int / subject=_Pasture type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  /*final model NESTLING MASS BIOMASS*/
    Proc glimmix data=mass method=laplace plots=residualpanel /*nestlingmass*/;
  class NestID _Pasture Pasture_Patch Pasture_patch_year HerbYesNo HerbTreat GrazingYesNo;
  model Nestlingmass =  / dist=expo solution;
  random int / subject=_Pasture type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  /*######################################################################################
  PROVISIONING (BY CLIP)
  ######################################################################################*/
  
  /*importing PROVISIONING (BY CLIP) dataset*/
    PROC IMPORT OUT= WORK.DICKProvByClip 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest metrics\Analysis_Provisioning_Round1\ProvDataSetRound2\Provisioning_byClip.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  /*global model PROVISIONING (BY CLIP)*/
    Proc glimmix data=DICKProvByClip method=laplace plots=residualpanel /*prov*/;
  class Nest_ID Nest_ID_Session Nest_ID Behavior Pasture HerbTreat HerbYesNo GrazingYesNo;
  model InstancesPerHr = Nestling_Age_Days Total_Nestling__ GrazingYesNo HerbTreat ordinaldate / dist=expo solution;
  random int/subject=Nest_ID_session type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  /*final model PROVISIONING (BY CLIP)*/
    Proc glimmix data=DICKProvByClip method=laplace plots=residualpanel /*prov*/;
  class Nest_ID Nest_ID_Session Nest_ID Behavior Pasture HerbTreat HerbYesNo GrazingYesNo;
  model InstancesPerHr = Nestling_Age_Days Total_Nestling__ / dist=expo solution;
  random int/subject=Nest_ID_session type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  run;
  
  /*######################################################################################
  PATCH SCALE DATA: TOTAL NESTS + FLEDGLINGS
  ######################################################################################*/
  
  /*importing PATCH LEVEL data*/
    PROC IMPORT OUT= WORK.DICKPATCH 
  DATAFILE= "C:\Users\jjcoon2\Dropbox\Manuscripts\Dissertation
\Ch. 2\Dickcissel Collaborative Manuscripts\#3 Management and all nest m
etrics\PatchNestData_ForAnalysis10.csv" 
  DBMS=CSV REPLACE;
  GETNAMES=YES;
  DATAROW=2; 
  RUN;
  /*global model TOTALNESTS*/
    Proc glimmix data=dickpatch method=quad plots=residualpanel /*totalnests*/;
  class pasture_year Pasture_Patch Pasture Year FireTreat HerbYesNo HerbTreat GrazingYesNo GrazingTreat ;
  model TotalNests = HerbTreat GrazingYesNo / dist=p solution offset=Ln_Patchsize_ha ;
  random int/subject=pasture type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  run;
  
  /*final model TOTALNESTS*/
    Proc glimmix data=dickpatch method=laplace plots=residualpanel /*totalnests*/;
  class pasture_year Pasture_Patch Pasture Year FireTreat HerbYesNo HerbTreat GrazingYesNo GrazingTreat ;
  model TotalNests = HerbTreat GrazingYesNo/ dist=p link=log solution offset=ln_Patchsize_ha;
  random int/subject=pasture type=vc G solution;
  
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  estimate 'Con herbicide grazing no'		intercept 1 	HerbTreat 1 0 0		grazingyesno 1 0  	/ cl 	alpha=0.1 	ilink;
  estimate 'Spr herbicide grazing no' 	intercept 1 	HerbTreat 0 0 1		grazingyesno 1 0	/ cl 	alpha=0.1 	ilink;
  estimate 'SnS herbicide grazing no' 	intercept 1 	HerbTreat 0 1 0		grazingyesno 1 0	/ cl 	alpha=0.1 	ilink;
  
  estimate 'Con herbicide grazing yes'	intercept 1 	HerbTreat 1 0 0		grazingyesno 0 1  	/ cl 	alpha=0.1 	ilink;
  estimate 'Spr herbicide grazing yes' 	intercept 1 	HerbTreat 0 0 1		grazingyesno 0 1	/ cl 	alpha=0.1 	ilink;
  estimate 'SnS herbicide grazing yes' 	intercept 1 	HerbTreat 0 1 0		grazingyesno 0 1	/ cl 	alpha=0.1 	ilink;
  
  run;
  
  
  /*global model DICK FLEDGED*/
    Proc glimmix data=dickpatch method=quad plots=residualpanel /*fledglings*/;
  class pasture_year Pasture_Patch Pasture Year FireTreat HerbYesNo HerbTreat GrazingYesNo GrazingTreat ;
  model TotalDICKFledged = TotalNests HerbTreat GrazingYesNo/ dist=p solution;
  random int/subject=pasture type=vc G solution;
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  run;
  
  /*final model DICK FLEDGED*/
    Proc glimmix data=dickpatch method=quad plots=residualpanel /*fledglings*/;
  class pasture_year Pasture_Patch Pasture Year FireTreat HerbYesNo HerbTreat GrazingYesNo GrazingTreat ;
  model TotalDICKFledged = TotalNests HerbTreat / dist=p solution;
  random int/subject=pasture type=vc G solution;
  
  covtest 'glm|indep' indep;
  covtest 'zerog' zerog;
  covtest 'homogeneity' homogeneity;
  
  estimate 'No herbicide '		intercept 1 	HerbTreat 1 0 0	TotalNests 5. 	/ cl 	alpha=0.1 	ilink;
  estimate 'Spr herbicide '		intercept 1 	HerbTreat 0 0 1	TotalNests 5 	/ cl 	alpha=0.1	ilink;
  estimate 'SnS herbicide '		intercept 1 	HerbTreat 0 1 0	TotalNests 5 	/ cl 	alpha=0.1	ilink;
  
  run;
  
  
  
  
  
  
  