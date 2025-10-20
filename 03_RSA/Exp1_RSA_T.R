# Exp1_L1Model_RSA
# load pred file (contains single-trial accuracy & confidence)
# fit RSA models to decoding result of the conjunction variable
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
fsep<-.Platform$file.sep;
Dir_R<-path.expand("~/DataAnalysis/R")
Dir_RSA<-path.expand("~//DataAnalysis/RSA/Conjunction")
Dir_EDATA<-path.expand("~//DataAnalysis/Matlab/EEG")
Dir_BDATA<-path.expand("~/DataAnalysis/Matlab/BEH/w_ALLGRAND")
Dir_GRAND<-paste0(Dir_EDATA,"/w_ALLGRAND")

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

# Source Files
setwd(Dir_R)
source('basic_lib.R')

# Setting for Database
setwd(Dir_EDATA);
subs<-list.files( pattern="^A[[:digit:]]{3}");
subs<-subs[subs!="A010" & subs!="A018" & subs!="A025" & subs!="A035"] #exclude subs 10, 18, 25, and 35 for >33% AR
subs<-subs[subs!="A029"] #exclude sub 29 because marker file was not filled
subs<-subs[subs!="A012" & subs!="A026" & subs!="A036"] #exclude sub 12, 26, and 36 for chance-level decoding
fileN<-"FlankerConj"
modelV<-"MCOND"#TSCONJ/L1TSCONJ/TCTSCONJ/L1TSCONJ_TSCONJ
modelVs<-unlist(str_split(modelV,"_"))
f2Load<-paste0(modelV,'_nb_pred.rds');
dbName<-paste0("FIT_RSA_",modelV,"conj");
saveOn<-T;
dataL<-list();

# # Load beh data
setwd(Dir_BDATA)
ds_beh<-fread("FlankerConj_BehP.txt");
ds_beh[is.nanM(ds_beh)]<-NA

ds_beh <- ds_beh %>%
  dplyr::mutate(MCOND_cat = as.factor(paste0(Target, "_", Flanker))) %>%
  dplyr::mutate(MCOND = as.numeric(MCOND_cat)) %>%
  dplyr::rename(BLOCK = Block,
                TRIAL = Trial,
                SUBID = ID) %>% 
  dplyr::group_by(SUBID, BLOCK) %>% 
  dplyr::mutate(l1Error = lag(Error), prevConfTrial = lag(ConfTrial)) %>% 
  ungroup %>% 
  dplyr::filter(Error == 0 & l1Error == 0 & TRIAL > 1 & Practice == 0) %>% 
  dplyr::filter(RT < 2000) #remove extreme RTs

ds_beh <-data.table(ds_beh)
ds_beh[, (c("RSI_JIT","PairedKeyResp","L1Error")):=NULL]


# # Load Individuals' Conjunction RT
setwd(Dir_GRAND)
ds_cj<- fread("probe_by_flank_sub_RTs.txt") #mean RT for all class and participants
#======================================================================================================
#Load RSA Models (adjusted for different decoding variables)
#======================================================================================================
# # Adjust RSA model sources!
setwd(Dir_RSA)
mlist <- list.files(pattern = ".csv")
models<-lapply(as.list(mlist),function(f){fread(f)})
for (m in 1:length(mlist)){assign(gsub(".csv","",mlist)[m],models[[m]])}

# read in all .csv models (8x8) to be used in RSA
setwd(Dir_RSA)
flist <- list.files(pattern = ".csv")

for (matrx in 1:length(flist)) {
  name <- strsplit(flist[matrx], ".csv", fixed = TRUE)[[1]]
  df <- read.csv(flist[matrx], header = FALSE)
  #assign(name, apply(df,2,rev)) #only need for original model set (because I originally flipped the y axis)
  assign(name, df)
}

#======================================================================================================
#Fit RSA
#======================================================================================================

# Open Database  (erase old database of the same datatype, if it exists!!)
setwd(Dir_GRAND);
registerDoMC(4)
resultsA<-list();
if (file.exists(paste0(fileN,'_',dbName))){file.remove(list.files(pattern=dbName))}
dbCon=dbConnect(dbDriver("SQLite"),paste0(fileN,'_',dbName))# Open database connection

for (s in 1:length(subs)){#
  # Load data & set properties
  Dir_Data_i<-paste0(Dir_EDATA,fsep,subs[s],fsep,"DATASETS");
  setwd(Dir_Data_i);f2LoadL<-paste0(subs[s],"_",f2Load);
  pred<-data.table(readRDS(f2LoadL));#pred<-rio::import(f2LoadL)
  pred$SUBID <- pred$ID; pred$ID<-NULL
  subN <- gsub('[[:alpha:]]+', '', subs[s])
  mesVs<-colnames(pred)[colnames(pred) %in% LETTERS_ex(50)]
  idVs<-colnames(pred)[!colnames(pred) %in% c(mesVs,"acc")]
  nmbVs<-idVs[!idVs %in% c("FBAND","ELEC","CV")]# avoid character columns in idVs
  timeVs<-nmbVs[grepl('time',nmbVs)]# time coding columns
  castF1<-paste0("CLASS~",paste0(idVs[!idVs %in% c("SUBID")],collapse="+"))
  castF2<-paste0(paste0(idVs,collapse="+"),"~vars")
  
  # Step1: Melt data from wide (A:L) to long format (BLOCK,TRIAL,time)
  pred<-data.table::melt(pred,#this is way faster than gather/spread
                         id.vars = idVs,
                         measure.vars=mesVs,
                         variable.name="CLASS",value.name="PP")
  
  # Step1':Logit transform
  pred[PP==0,PP:=1e-10];pred[PP==1,PP:=1-1e-10]
  pred$PP<-psych::logit(pred$PP);
  
  # Check if all classes have complete cases 
  if (any(c(is.infinite(pred$PP) | is.nan(pred$PP)))){
    pred<-pred[!c(is.infinite(pred$PP) | is.nan(pred$PP)),]
    pred[,cmp:=uniqueN(CLASS),by=c("BLOCK","TRIAL",timeVs)]
    pred<-pred[cmp==length(mesVs),];pred[,cmp:=NULL]
  }
  
  sub_in_beh <- as.numeric(subN) %in% as.numeric(unique(ds_beh$SUBID))
  
  # Step2: Add condition of Modeled variable from behavior (if necessary!)
  if (!any(grepl("obs",names(pred)))){
    ds_behi<-ds_beh[as.numeric(SUBID)==as.numeric(subN),c(c("SUBID","BLOCK","TRIAL"),modelV),with=F]
    pred<-merge(pred,ds_behi,by=c("SUBID","BLOCK","TRIAL"));setnames(pred,modelV,"obs")
    }
  
  # Step3-1:Long format (BLOCK,TRIAL,time) to another wide format (A:L to all combination of BLOCK,TRIAL,time)
  # Then, add model terms (e.g.,TASK,STIM,RESP,CONJ)
  if(is.factor(pred$obs)){pred[,obs:=as.numeric(obs)]}
  clist<-which(mesVs %in% LETTERS_ex(50))
  
  # Step (4): Prepare control vector of RT
  CJCTR<-scale(ds_cj[ds_cj$SUBID==as.numeric(subN),]$subRT)[,1];
  
  if (sub_in_beh){
    results<-
      foreach(cc=clist) %dopar% {
        d<-pred[obs==cc,];
        predF<-data.table::dcast(d,castF1, value.var="PP")
        
        # This part could be more abstract with "eval" 
        # TSCONJ (Rule x Resp/Stim) with CJ RT control
        predF$TARG_RSA<-TargetID_M[,cc];#TARG MODEL
        predF$FLANK_RSA<-FlankerID_M[,cc];#FLANK MODEL
        predF$CONJ_RSA<-Conjunction_M[,cc];#CONJ MODEL
        predF$CONG_RSA<-Congruency_M[,cc];#CONGRUENCY MODEL
        #predF$TFCONF_RSA<-TF_Confusion_M[,cc];#TF/FT CONFUSION MODEL
        #predF$TARG_UNBOUND_RSA<-TargetUnbound_M[,cc];#Target Unbound MODEL
        #predF$FLANK_UNBOUND_RSA<-FlankerUnbound_M[,cc];#Flanker Unbound MODEL
        predF$CTR_RSA<-CJCTR;#CONTROL MODEL
        dvIdx<-grepl("CLASS|FLANK_RSA|TARG_RSA|CONJ_RSA|CONG_RSA|CTR_RSA",colnames(predF))
        
        # Step3-2:Regression all at once!
        Y<-as.matrix(predF[,which(!dvIdx),with=FALSE]);#DVs(this indexing takes time...)
        X<-as.matrix(predF[,which(dvIdx)[-1],with=FALSE]);#IVs(first var is CLASS!)
        # r<-coeff(.lm.fit(cbind(1,X),Y));#fastest, but only gives unstandardized coefficients
        r<-ls.print(lsfit(X,Y),print.it=F);#includes intercept!
        
        # Step3-3:Summarize results!
        names(r$coef.table)<-colnames(predF)[!dvIdx]
        estList<-rownames(r$coef.table[[1]]);#list of vars for estimates
        r<-rbindlist(lapply(r$coef.table,as.data.frame),idcol=TRUE)
        setnames(r,old=c("t-value","Pr(>|t|)","Std.Err"),new=c("tvalue","pvalue","SE"))
        
        #Rename variables and add basic identifier variables
        r[,c("SUBID",modelV):=list(pred$SUBID[1],cc)]
        r[,vars:=base::rep(estList,(dim(r)[1]/length(estList)))]
        
        #Restore labels from list names, then convert to wide format
        bltrtime<-str_split_fixed(r$.id,"_",n=Inf);#no assumption for .id tokens
        r[,c(idVs[-4]):=narray::split(bltrtime,along=2)]# assume SUBID is the fourth column
        r[,(nmbVs):=lapply(.SD,as.numeric),.SDcols=nmbVs]
        r<-r[vars!="Intercept",c(idVs,"vars","tvalue"),with=F]
        
        # Need wide format(dcast) // long-wide(dcast) & wide-long(melt)
        r_wide<-dcast(r,castF2,value.var="tvalue")
        return(r_wide)
      }
    
    
    # Summarize all results!!!!
    results <- rbindlist(results); results_backup <- results
    setorderv(results,nmbVs)
    depV<-colnames(results)[!colnames(results) %in% c(idVs,"CTR_RSA")]
    
    print(sprintf(subs[s]))
    
    # Detect outliers within each time point and replace with NA
    for (cc in depV){
      results[,std:=lapply(.SD,sd),by=c(timeVs),.SDcols=cc]
      results[,bad:=abs(results[[cc]]) > std*5] # outlier based on std
      results[,(cc):=ifelse(bad,NA,results[[cc]])]
      print(sprintf("For %s, %f2 percent of rows were removed",cc,(length(which(results$bad))/dim(results)[1])*100))
      results<-results[bad!=T,] # remove bad trials altogether??
      #results[,bad:=outlier(results[[cc]],logical=T)] # outlier based on package(only one value??)
    }
    
    # Merge to behavioral template and put into DB
    results[,c("std","bad"):=NULL]# reduce data!
    predB<-merge(results,ds_beh,by=c("SUBID","BLOCK","TRIAL"));#use ds_beh(containing all vars)
    dbWriteTable(dbCon,name=dbName,value=predB,row.names=FALSE,append=TRUE);#Try to append!
    
  } else {
    print(paste0("Didn't run ", subN))
  }
}
