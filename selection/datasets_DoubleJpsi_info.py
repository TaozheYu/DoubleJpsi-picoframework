#!/usr/bin/env python
#Insert dataset you want to analyze
#The format is that of a list of array in python
#Eeach array has 8 entries
#0. "year": year of the dataset
#1. "datasetFlag": may stands for "SigMC", "BkgMC", the name of the dataset (e.g. "JetHT", "MET")
#2. "processName": a shortcut of the datasetName. If a dataset is an "ext", add this label and its number (e.g. DYJetsToLL_M-50_HT-100to200ext1) 
#3. "datasetName": the datasetName
#4. xSec in pb
#5. The total number of events of the dataset. To find it, use the command: dasgoclient --query="file dataset=dataset name | sum(file.nevents)". For data you should run 1 time and see the num as in data we filter the evt with the json file (to make sure you run over all evt see e.g. the tot num of files processed) 
#6. "sumgenWeight": the  sum of the genWeight of a MC dataset. Set 0 for data dataset.
#7. Label containing the analyses that need to analyze the dataset 
#Note
#1. The order of insertion of each year is first Data and then MC and must respect to avoid problems in the analysis chain
#2. The first time your run you do not know what the value of sumgenWeight is. So you need to be careful that your run over the entire dataset to get this value (see the num of files or the tot num of abs evt). Up to then, set sumgenWeight = 1 

datasets_info = [
#####
##   Type:MC
#####  
#####
##   Background
#####
#####
##  Year: 2016
#####
#SPS 
#[2016,"BkgMC","SPS","/SPS_ToJPsiJPsi_TuneCP5_Inclusive_13TeV-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",1,7.2867e+05,1,"DoubleJpsi"],
#DPS
#[2016,"BkgMC","DPS","/DPS_ToJPsiJPsi_TuneCP5_Inclusive_13TeV-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",1,1.973008e+06,1,"DoubleJpsi"],
#####
##  Year: 2017
#####
#SPS 
[2017,"BkgMC","SPS","/SPS_ToJPsiJPsi_TuneCP5_Inclusive_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM",1,1.198173e+06,1,"DoubleJpsi"],
#DPS
[2017,"BkgMC","DPS","/DPS_ToJPsiJPsi_TuneCP5_Inclusive_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM",1,2.568719e+06,1,"DoubleJpsi"],
#####
##  Year: 2018
#####
#SPS 
[2018,"BkgMC","SPS","/SPS_ToJPsiJPsi_TuneCP5_DP_13TeV-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM",1,2.079682e+06,1,"DoubleJpsi"],
#DPS
[2018,"BkgMC","DPS","/DPS_ToJPsiJPsi_TuneCP5_DP_13TeV-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM",1,4.544333e+06,1,"DoubleJpsi"],
#####
##   Type:data
#####   
###### 2016 data ######
[2016,"Charmonium","Run2016B-02Apr2020_ver2-v1","/Charmonium/Run2016B-02Apr2020_ver2-v1/NANOAOD",1,5.1713106e+07,0,"DoubleJpsi"],
[2016,"Charmonium","Run2016C-02Apr2020-v1","/Charmonium/Run2016C-02Apr2020-v1/NANOAOD",1,1.9535374e+07,0,"DoubleJpsi"],
[2016,"Charmonium","Run2016D-02Apr2020-v1","/Charmonium/Run2016D-02Apr2020-v1/NANOAOD",1,3.2176529e+07,0,"DoubleJpsi"],
[2016,"Charmonium","Run2016E-02Apr2020-v1","/Charmonium/Run2016E-02Apr2020-v1/NANOAOD",1,2.6549568e+07,0,"DoubleJpsi"],
[2016,"Charmonium","Run2016F-02Apr2020-v1","/Charmonium/Run2016F-02Apr2020-v1/NANOAOD",1,2.1339582e+07,0,"DoubleJpsi"],
[2016,"Charmonium","Run2016G-02Apr2020-v1","/Charmonium/Run2016G-02Apr2020-v1/NANOAOD",1,6.2080523e+07,0,"DoubleJpsi"],
[2016,"Charmonium","Run2016H-02Apr2020-v1","/Charmonium/Run2016H-02Apr2020-v1/NANOAOD",1,7.1189433e+07,0,"DoubleJpsi"],
###### 2017 data ######
[2017,"Charmonium","Run2017B-02Apr2020-v1","/Charmonium/Run2017B-02Apr2020-v1/NANOAOD",1,1.9863860e+07,0,"DoubleJpsi"],
[2017,"Charmonium","Run2017C-02Apr2020-v1","/Charmonium/Run2017C-02Apr2020-v1/NANOAOD",1,7.8721035e+07,0,"DoubleJpsi"],
[2017,"Charmonium","Run2017D-02Apr2020-v1","/Charmonium/Run2017D-02Apr2020-v1/NANOAOD",1,3.2714414e+07,0,"DoubleJpsi"],
[2017,"Charmonium","Run2017E-02Apr2020-v1","/Charmonium/Run2017D-02Apr2020-v1/NANOAOD",1,6.0891440e+07,0,"DoubleJpsi"],
[2017,"Charmonium","Run2017F-02Apr2020-v1","/Charmonium/Run2017D-02Apr2020-v1/NANOAOD",1,7.4190347e+07,0,"DoubleJpsi"],
###### 2018 data ######
[2018,"Charmonium","Run2018A-02Apr2020-v1","/Charmonium/Run2018A-02Apr2020-v1/NANOAOD",1,9.3082605e+07,0,"DoubleJpsi"],
[2018,"Charmonium","Run2018B-02Apr2020-v1","/Charmonium/Run2018B-02Apr2020-v1/NANOAOD",1,4.6847388e+07,0,"DoubleJpsi"],
[2018,"Charmonium","Run2018C-02Apr2020-v1","/Charmonium/Run2018C-02Apr2020-v1/NANOAOD",1,4.5597510e+07,0,"DoubleJpsi"],
[2018,"Charmonium","Run2018D-02Apr2020-v1","/Charmonium/Run2018D-02Apr2020-v1/NANOAOD",1,2.17135785e+08,0,"DoubleJpsi"],
]
