# Conflict Adaptation project pipeline workflow

## MATLAB Preprocessing pipeline** 
<u>Necessary Toolbox: EEGLAB, ERPLAB</u>
  
**01.PrepAll_EEGLAB_Conj_Analysis_MM.m
-Parameter settings for artifact rejection and TF decomposition.
-Handles the main preprocessing loop. 

**02.TimeFreq_EEGLAB_Preprocess_Conj_MM.m**
-Merge EEG, header, and marker files.
-Epoch data.
-Filter data for above threshold muscle noise and eye movements. 

**03.TimeFreq_EEGLAB_RunTF_Conj_MM.m**
-Time-Frequency analysis: wavelet convolution. 
-Frequency decomposed data available through this function but also handled below for transfer to R. 

**04.FlankConj_BehAnalysis**

**05.FlankerConj_Analysis**
-Line up and merge EEG and behavioral datasets with tf_meargeEEGBEH_MM.m.
-Convert merged dataset to HDF5 format to transfer to R.
-HDF5 file generated in tf_trans2env_MM.m.

## R, Classification and RSA Scripts  
**01.FlankerConj_Decode_T.R**
-Script that handles the decoding analysis.
-Output:_nb_pred.rds & MCOND_nb_pred.rds:decoder accuracy and classification probabilities for each class across all subids, electrodes, 
timepoints, and trials.
 
**02.Exp1_RSA_T.R**
-OLS regression predicting classification probability from RSA models for each subid and class.
 
**03.MixEffects_Model_RSA.R**
-Mass univariate mixed effects model analysis for each conflict condition and trial time point predicting RT from OLS regression 
t-values.
-Predicts how classification probability for each RSA model class at each timepoint predicts RT. 
  
**04.MixEffects_Model_Plots.R**
-Plots results from OLS regression, mixed effects model, RT as a function of conflict, and  RSA model t-values as a function conflict.
