# Conflict Adaptation - Project Pipeline Workflow

## MATLAB Preprocessing pipeline 
<ins>Necessary Toolbox: EEGLAB, ERPLAB</ins>
  
### 1. PrepAll_EEGLAB_Conj_Analysis_MM.m
- Parameter settings for artifact rejection and TF decomposition. 
- Handles the main preprocessing loop. 

### 2. TimeFreq_EEGLAB_Preprocess_Conj_MM.m
- Merge EEG, header, and marker files. 
- Epoch data. 
- Filter data for above threshold muscle noise and eye movements. 

### 3. TimeFreq_EEGLAB_RunTF_Conj_MM.m
- Time-Frequency analysis: wavelet convolution. 
- Frequency decomposed data available through this function but also handled below for transfer to R. 

### 4. FlankConj_BehAnalysis

### 5. FlankerConj_Analysis
- Line up and merge EEG and behavioral datasets with tf_meargeEEGBEH_MM.m. 
- Convert merged dataset to HDF5 format to transfer to R. 
- HDF5 file generated in tf_trans2env_MM.m. 

## R, Classification and RSA Scripts  
### 1. FlankerConj_Decode_T.R
- Script that handles the decoding analysis. 
- Output:_nb_pred.rds & MCOND_nb_pred.rds:decoder accuracy and classification probabilities for each class across all subids, electrodes, 
timepoints, and trials. 
 
### 2. Exp1_RSA_T.R
- OLS regression predicting classification probability from RSA models for each subid and class. 
 
### 3. MixEffects_Model_RSA.R
- Mass univariate mixed effects model analysis for each conflict condition and trial time point predicting RT from OLS regression 
t-values. 
- Predicts how classification probability for each RSA model class at each timepoint predicts RT. 
  
### 4. MixEffects_Model_Plots.R
- Plots results from OLS regression, mixed effects model, RT as a function of conflict, and  RSA model t-values as a function conflict. 
