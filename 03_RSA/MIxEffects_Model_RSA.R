#Mass univariate mixed effects model analysis for each experiment time point predicting RT from RSA t-values 

#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
fsep<-.Platform$file.sep;
Dir_R<-path.expand("~/DataAnalysis/R")
Dir_RSA<-path.expand("~/DataAnalysis/RSA/Conjunction")
Dir_EDATA<-path.expand("~/DataAnalysis/Matlab/EEG")
Dir_BDATA<-path.expand("~/DataAnalysis/Matlab/BEH/w_ALLGRAND")
Dir_GRAND<-paste0(Dir_EDATA,"/w_ALLGRAND")

Dir_flash <- "/Volumes/RESEARCH/Flanker"
# Load libraries
library(data.table)
library(dplyr)
library(broom)
library(rhdf5)
library(caret)
library(foreach)
library(doMC)
library(tidyr)
library(binhf)
library(GGally)
library(RColorBrewer)
library(RSQLite)
library(lazyeval)
library(stringr)
library(psych)
library(lme4)
library(RcppEigen)
library(rio)
library(gridExtra)

# non analysis related libraries
library(parallel)
# Source Files
setwd(Dir_R)
source('basic_lib.R')

#Define basic info
srate<-4
timeL<-seq(-200,1200,4);#Stim_LK:seq(-600,600,4),Resp_LK:seq(-800,0,4),NA=SEG,JMB,CRS DATA
timeS<-seq(-200,1200,12)#c(0,600)#seq(-600,600,4)#c(begining,end,window)
refRLV <- NA;# Realignment: re-epoch data referencing refRLV timing

# Local functions, kudos to Melissa
if (!is.na(refRLV)){adj<-which(timeL==0)}else{adj<-0}
cutF <- function(t){as.numeric(cut(t,breaks=timeS,include.lowest=T))}
idxF<-function(i) which.min(abs(timeL-i))-adj # normal

# Database= PP: REG / RSAFits: FIT_RSA_TSCONJ,FIT_RSA_TSCONJ_RL
setwd(Dir_GRAND)
expN<-"FlankerConj";
#freqband<-"theta"
#tableN<-paste0("FIT_RSA_MCONDconj_",freqband);
tableN<-"FIT_RSA_MCONDconj";
depV<-c("TARG_RSA","FLANK_RSA","CONJ_RSA", "CONG_RSA", "CTR_RSA"); # "CONJ_RSA","TARG_RSA","FLANK_RSA","TFCONF_RSA","CONG_RSA","TARG_UNBOUND_RSA","FLANK_UNBOUND_RSA","CTR_RSA"
grpV<-c("SUBID","time","RT", "BLOCK", "TRIAL", "prevConfTrial", "ConfTrial")
#grpV<-c("SUBID","time","RT", "BLOCK", "TRIAL")


# Open database connection
ds_cj<- fread("probe_by_flank_sub_RTs.txt")# Load Individuals' Conjunction RT

#Fit_RSA <- "FIT_RSA_MCONDconj";#"FlankerConj_FIT_RSA_MCONDconj"
dbCon=dbConnect(dbDriver("SQLite"),paste0("FlankerConj_",tableN));#paste0(Dir_GRAND,fsep,Fit_RSA));#open db connection and pass path to db
#tablename <- dbListTables(dbCon);
dscheck=as.data.table(dbGetQuery(dbCon,paste0('SELECT * FROM ',tableN,' LIMIT 10000')));
dbCL=dbListFields(dbCon,dbListTables(dbCon)); #list column/variable names of db

subs <- as.data.table(dbGetQuery(dbCon,paste0('SELECT DISTINCT SUBID FROM ',tableN)))#check subids

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AGGREGATION: grouping SUBID x gvar in database with flexible query sentences
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregation with SQ laguage
tIDX<-sapply(range(timeS),idxF);whereQ="";print(tIDX)
selectQ<-paste0('SELECT ',paste(unique(c(grpV)),collapse = ", "))
compQ<-paste(sprintf(", AVG(%s) AS %s",depV,depV),collapse = "")
if (!any(is.na(timeL))){whereQ <- sprintf(" WHERE %s >= %d AND %s <= %d",grpV[2],tIDX[1],grpV[2],tIDX[2])}
groupQ<-paste0(' GROUP BY ',paste(unique(c(grpV)),collapse = " , "))
querySTW<-paste0(selectQ,compQ,sprintf(" FROM %s",tableN),whereQ,groupQ);
print(querySTW);
dsA=as.data.table(dbGetQuery(dbCon,querySTW));



# Keep constant time labels
timeC<-colnames(dsA)[grepl("time",colnames(dsA))]
if (!any(is.na(timeL)) & diff(unique(dsA[,get(timeC)]))[1]!=srate){
  dsA[,(timeC):=lapply(.SD,function(t)timeL[t]),.SDcols=timeC]
}


dsA$BLOCK<-as.numeric(scale(dsA$BLOCK,center=T,scale=F));
dsA$TRIAL<-as.numeric(scale(dsA$TRIAL,center=T,scale=F));

#calculate residuals for RT 
ds_resid <- dsA %>% mutate(logRT = log(RT)) %>% 
  group_by(SUBID, BLOCK, TRIAL) %>% 
  dplyr::summarize(RTlog = mean(logRT, na.rm = TRUE)) %>% 
  ungroup()
rt_resid <- summary(lm(RTlog~BLOCK+TRIAL+I(BLOCK^2)+I(TRIAL^2), ds_resid))
ds_resid$RTlog_resid <- rt_resid$residuals

dsA <- merge(dsA, ds_resid, by = c("SUBID","BLOCK","TRIAL"))

# # <Models with class probability>
depV<-"RTlog_resid"

#formula_pw<-paste0(depV,'~BLOCK+(1|SUBID)')

#dsA <- dsA %>% dplyr::mutate(FLANK_CONJ_RSA_DIF = FLANK_RSA - CONJ_RSA)
#formula<-paste0(depV,'~1+TARG_RSA+FLANK_CONJ_RSA_DIF+CONG_RSA+CTR_RSA+(1|SUBID)');
formula<-paste0(depV,'~1+FLANK_RSA+TARG_RSA+CONJ_RSA+CONG_RSA+CTR_RSA+(1|SUBID)');

vars<-colnames(dsA)
dsA$RT<-log(dsA$RT)

dsAll <- dsA

#dsA <- dsAll[dsAll$prevConfTrial==0 & dsAll$ConfTrial==1,]

dsConf <- dsAll[dsAll$ConfTrial==1,]#Incongruent trials
dsnoConf <- dsAll[dsAll$ConfTrial==0,]#Congruent trials
dspConf <- dsAll[dsAll$prevConfTrial==1,]#n-1 incongruent trials
dspnoConf <- dsAll[dsAll$prevConfTrial==0,]#n-1 congruent trials

confCon <- c("dsConf","dsnoConf","dspConf","dspnoConf");
conflist <- setNames(vector("list",length(confCon)), confCon)
# Run time-series regression
registerDoMC(4)#3max for laptop
ptm<- proc.time()
#============================================================
#*****Rerun this section for conflict and no conflict trials
#============================================================
#loop through each conflict condition grouped data set 
for (dsC in 1:length(confCon)){
  cl <- names(conflist[dsC]);
  r<-
  foreach(t=unique(get(cl)[[grpV[2]]]))%dopar% {
    #Load all data ahead (its okay for data < 3GB)
    ds<-get(cl)[get(grpV[2])==t,]
    #ds<-dsA[time>=t-aveW & time<=t+aveW,lapply(.SD,mean),keyby=c("SUBID","BLOCK","TRIAL")]#moving ave
    
    #Prewitten DV!
    #pw<-lmer(formula_pw,data=ds,control=lmerControl(calc.derivs=FALSE));#
    #ds[[depV]]<-residuals(pw);
    
    #Fit regression!!===========================================================
    controlSL<-lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6))#bobyqa,Nelder_Mead,L-BFGS-B
    m<-lme4::lmer(formula,data=ds,control=controlSL);#Fitting
    #r<-tidy(m,effects=c("fixed"),conf.int=T,conf.method="Wald");
    
    # Summarize results
    #print(r)
    m2 <- summary(m)
    #std_beta_MLM(m)@fixef
    
    #Summarize results
    #r<-tidy(m,effects=c("fixed"));
    #r<-std_beta_MLM(m)@fixef %>% data.table()
    r <- data.table(m2$coefficients)
    r$predictor <- rownames(m2$coefficients)
    r$time<-t
    return(r)
  }
  rc <- paste0('r',cl)
  assign(rc,r)
  #get(rc)
  print(dsC)
 
}
proc.time() - ptm



# results <- rbindlist(r); results_backup <- results
# results <- results[results$predictor != "(Intercept)" & results$predictor != "CTR_RSA",]

#ds_list <- mget(confCon,envir=.GlobalEnv, ifnotfound=list(NULL))

resultsConf <- rbindlist(rdsConf);
resultsnoConf <- rbindlist(rdsnoConf);#check why theres no CONG_RSA t-value
resultspConf <- rbindlist(rdspConf);
resultspnoConf <- rbindlist(rdspnoConf);

#combined lists as dataframes for each conflict condition dataset
results1 <- resultsConf[resultsConf$predictor != "(Intercept)" & resultsConf$predictor != "CTR_RSA",]
results2 <- resultsnoConf[resultsnoConf$predictor != "(Intercept)" & resultsnoConf$predictor != "CTR_RSA",]
results3 <- resultspConf[resultspConf$predictor != "(Intercept)" & resultspConf$predictor != "CTR_RSA",]
results4 <- resultspnoConf[resultspnoConf$predictor != "(Intercept)" & resultspnoConf$predictor != "CTR_RSA",]

#dir_fig <- ("/Volumes/Data2/EEG projects/UM_flanker_eeg/DataAnalysis/Figures");
setwd(Dir_flash)
save("results1",file="rIncongruent.rds")
save("results2",file="rCongruent.rds")
save("results3",file="rpIncongruent.rds")
save("results4",file="rpCongruent.rds")



  


