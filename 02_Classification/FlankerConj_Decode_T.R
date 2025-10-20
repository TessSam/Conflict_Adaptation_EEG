# FlankerConj_Decode
# peforms basic time-series decoding
# 1. Keep this script (FlankerConj_Decode) most updated version!
# 2. DO NOT USE TOO MANY ABSTRACT SELF-DEFINED FUNCTIONS!
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
fsep <- .Platform$file.sep;
#AK comp:
#Dir_R<-path.expand("~/Desktop/FlankerConj")
#Dir_EDATA<-Dir_R
#Dir_BDATA<-Dir_R
#Dir_GRAND<-Dir_R
Dir_R <- path.expand("~/Dropbox (University of Oregon)/P.FlankerConj/DataAnalysis/R")
Dir_EDATA <- path.expand("~/Dropbox (University of Oregon)/P.FlankerConj/DataAnalysis/Matlab/EEG")
Dir_BDATA <- path.expand("~/Dropbox (University of Oregon)/P.FlankerConj/DataAnalysis/Matlab/BEH/w_ALLGRAND")
Dir_GRAND <- paste0(Dir_EDATA,"/w_ALLGRAND")


# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(rhdf5)
library(caret)
library(foreach)
library(doMC)
library(binhf)
library(pROC)
library(RColorBrewer)
library(scales)
library(forcats)
library(stringr)



# Source Files
setwd(Dir_R)
#source('basic_theme.R')
source('basic_lib.R')

#Load behavioral data
setwd(Dir_BDATA)
ds_b <- fread("FlankerConj_BehP.txt") #was originally ConflictFlanker_BehP.txt
ds_b[is.nanM(ds_b)] <- NA;

#ds_b <- ds_b %>% dplyr::filter(Target != Flanker) #use only incongruent trials

ds_b <- ds_b %>%
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

sub_rt <- ds_b %>% 
  dplyr::group_by(SUBID, MCOND) %>% 
  dplyr::summarize(subRT = mean(RT, na.rm = T))
setwd(Dir_GRAND)
write.table(sub_rt, "probe_by_flank_sub_RTs.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)
ds_b <-data.table(ds_b)

#Analysis Setting
setwd(Dir_EDATA)
subs <- list.files(pattern = "^A[[:digit:]]{3}");
subs<-subs[subs!="A010" & subs!="A018" & subs!="A025" & subs!="A035"] #exclude subs 10, 18, 25, and 35 for >30% AR
subs<-subs[subs!="A029"] #exclude sub 29 because marker file was not filled

modelV<-"MCOND" #was originally just "MCOND"
sessN<-"nb";#nb = regular decoding, nppr = shuffled label
balanceV<-"MCOND" #MCOND
cutV=NA;# 
saveACC<-T; #(T)
savePRED<-T; #(T)
saveCM<-F; #confusion matrix (T)
saveIMP<-F; #importance of different features for classifier (F)
saveROC<-F; #more standard way to get decoder accuracy with only 2 classes (area under curve) (F)
useFreq<-F; #collapse across all frequency bands, unless T

# Classification Setting
#"lda" = linear discriminant analysis
#"knn" = k nearest means
#"svmRadial" = support vector machine
#"nb" = naive bayes
#"rf" = random forest
#"svmLinear3" = L2 Regularized Support Vector Machine (dual) with Linear Kernel
#"pda" = penalized linear discriminant analysis
method <- 'pda';
formula <- as.formula(paste(modelV,' ~ .'))
metric <- "Accuracy";
control <- caret::trainControl(method = "repeatedcv",
                      number = 5, repeats = 10, p = 0.75,
                      selectionFunction = "oneSE",
                      classProbs = TRUE, allowParallel = TRUE, savePredictions = TRUE)

#Feature labels
freqL <- c("Delta","Theta","Alpha","Beta","Gamma");
#elecL <- fread(paste0(Dir_R,fsep,"AS-32_NO_REF.bvef"))
elecL <- 2:31
varL <- expand.grid(elec = as.vector(elecL),freq = freqL)
varL <- str_c(varL$freq,"_",varL$elec)

#subs <- subs[12:length(subs)]
#s <- subs[1]

#do for each frequency band separately
freqs <- c("_delta", "_theta", "_alpha", "_beta")
freqIDX <- 1 #delta (will need to run for each frequency band: "Delta","Theta","Alpha","Beta")

# IN THE LOOP FOR SUBJECT!!
for (s in subs) {
  
  sprintf("Subject %s:",s)
  
  # STEP 1: Merging data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load EEG data via HDF5 file
  Dir_Data_i <- paste0(Dir_EDATA,fsep,s,fsep,"DATASETS")
  setwd(Dir_Data_i); H5close();
  f2load <- list.files(path = Dir_Data_i, pattern = '*_SINGLE_allFB.h5', full.names = TRUE)
  #f2load<-list.files(path=Dir_Data_i,pattern='*RESPLOCK.h5',full.names=TRUE)
  ds_eeg <- h5read(f2load, "eegpower") #ds_info<-h5ls(f2load); #NOTE TO SELF: need to figure out how to match these data to BLOCK*TRIAL from partial data
  btIDX <- data.table(h5read(f2load ,"IDX")) %>% dplyr::rename(BLOCK = V1, TRIAL = V2);
  btIDX$count <- 1:length(btIDX$BLOCK)
  
  # Get Individual data to match to EEG data (this will reflect EEG artifact rejection)
  # ds_bl is typically smaller than original ds_b because btIDX reflects AR...
  subN <- as.numeric(gsub("A","",s));# Filter data and keep indexes
  #ds_bl<-merge(btIDX,ds_b[SUBID==subN],by=c('BLOCK','TRIAL')) #original line
  ds_bl <- merge(ds_b[SUBID == subN], btIDX,by = c('BLOCK','TRIAL')) #line I made for BR only decoding
  
  ds_eeg <- ds_eeg[ds_bl$count[],,,] #only save rows (BLOCK*TRIAL) that correspond to rows in the behavioral data
  
  modelV <- "MCOND"
  
  # Applies 1) filtering, 2) cutting, 3) balancing, and 4) factorization
  prepSet <- list(
    filtIDX = !is.na(ds_bl[[modelV]]) & ds_bl$Error == 0 & ds_bl$L1Error == 0 & ds_bl$TRIAL > 1 & ds_bl$BLOCK > 1,
    modelV = modelV, balanceV = balanceV, cutV = cutV)
  
  pp <- preprop_decode(ds_bl, ds_eeg, prepSet);
  ds_bl <- pp$ds_beh; ds_eeg <- pp$ds_eeg
  
  sprintf("Subject %d!", subN)
  
  # STEP (2):Aggregation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # STEP 2: Classification~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  numCores <- 3; # was originally 4
  registerDoMC(numCores);
  ptm<- proc.time();
  dims<-dim(ds_eeg);
  timeIDX<-1:dims[2];elecIDX<-1:dims[4];
  sprintf("Classification starts!! Retained data are %d trials- %d samples- %d fband- %d electrodes",dims[1],dims[2],dims[3],dims[4])
  
  if (!useFreq) {
    freqIDX<-1:dims[3];
  }
  
  results<-
    foreach(t=timeIDX[1:350]) %dopar% {
    #for (t in timeIDX) {
      
      #Get data of one point!
      d<-ds_eeg[,t,freqIDX,elecIDX];#pick frequency band to include!
      dim(d)<-c(dim(d)[1],prod(dim(d)[2:3],na.rm=T))#adjust dimensions !
      #dim(d)<-c(prod(dim(d)[1:2],na.rm=T),dim(d)[3])#adjust dimensions !
      
      #Assigne label
      C <- ds_bl[[modelV]];if(grepl("nppr",sessN)){C<-sample(C)} 
      d<-data.table(C=C,data.table(d));#make sure to data.table(d) before appending!
      setnames(d,"C",modelV)
      
      #Run classification
      m<-caret::train(formula, data=d, method=method, metric=metric, trControl=control);
      
      # #Get all results
      r <- summary_decode(m,varL,list(CM=saveCM,IMP=saveIMP,ROC=saveROC))
      r <- r[sapply(r, function(x) dim(x)[1]) > 0];# remove empty slots
      r_f <- lapply(r, cbind,ID=subN,time=timeIDX[t])# don't lapply r_f repeatedly
      
      r_f$r_prob <- r_f$r_prob %>% mutate(BLOCK=ds_bl$BLOCK,TRIAL=ds_bl$TRIAL) # should be safe!
      # r_f$r_prob <- r_f$r_prob %>% dplyr::mutate(FBAND = rep(freqL[1:4], each = 1680),
      #                                            BLOCK=rep(ds_bl$BLOCK,4),TRIAL=rep(ds_bl$TRIAL,4)) # should be safe!
      
      return(r_f)
    }
  proc.time() - ptm
  
  # Summarize all results
  rALL<-as.list(as.data.frame(do.call(rbind, results)))
  
  if (useFreq) {
    modelV <- paste0(modelV,freqs[freqIDX])
  }
    
  # Result 1:Overall accuracy
  setwd('/Users/tesufuai/Documents/Masters/EEG Project')
  r_accG<-rbindlist(rALL$r_acc);
  if (saveACC){saveRDS(r_accG,paste0(s,"_",modelV,"_",sessN,"_accG.rds"))}
  
  # Result 2:Single-trial accuracy and confidence
  # Both acc and prob would be 0~1 since its cross-validated!
  r_pred<-rbindlist(rALL$r_prob) %>% dplyr::select(-rowIndex);
  if (savePRED){saveRDS(r_pred,paste0(s,"_",modelV,"_",sessN,"_pred.rds"))}
  
  # Result 3:Confusion matrix
  r_cm<-rbindlist(rALL$r_cm)
  if (saveCM){saveRDS(r_cm,paste0(s,"_",modelV,"_",sessN,"_cm.rds"))}
  
  # Result 4:Importance map
  r_imp<-rbindlist(rALL$r_imp)
  if (saveIMP){saveRDS(r_imp,paste0(s,"_",modelV,"_",sessN,"_imp.rds"))}
  
  # STEP 3: Quickly plotting?~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot related setting (for quick check!)
  #theme_set(theme_bw(base_size = 20))#32/28
  CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
  chance<-1/length(unique(ds_bl[[modelV]]))
  
  # Overall Accuracy
  quartz(width=7.5,height=4.5); #theme_set(theme_bw(base_size = 20))#32/28
  print(ggplot(data=r_accG,aes(x=time,y=Accuracy,ymin=Accuracy-AccuracySD,ymax=Accuracy+AccuracySD)) +
          geom_vline(xintercept=51,linetype=1,size=1)+annotate("text",x=85, y=chance-0.05,label="Flanker")+
          geom_vline(xintercept=160,linetype=1,size=1)+annotate("text",x=165, y=chance-0.05,label="Stimulus")+
          geom_hline(yintercept=chance,linetype=1,size=1)+
          geom_ribbon(alpha=0.2,linetype=0,fill="red")+#
          ggtitle(s)+
          geom_line(size=1.5,color="red"))
}


# Send email to me 
#notifyM(paste0("ConflictFlanker_",date()))

