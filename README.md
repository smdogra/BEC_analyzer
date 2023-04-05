# BEC_analyzer
Analyzer for BEC

Instructions how to run 

cmsrel CMSSW_12_5_0

cd CMSSW_12_5_0/src

cmsenv

git clone https://github.com/smdogra/BEC_analyzer.git

cd BEC_analyzer

root -l  'beccorrelation_analyzer.C("inputfiles/pPb_8160/Data_MB/Pbgoing/MB_PD1_Pbgoing_part0.txt","HiForest_pPb_Data_MB_PD1_Pbgoing_job_0",0,1)' # This will run mixing  

root -l  'beccorrelation_analyzer.C("inputfiles/pPb_8160/Data_MB/Pbgoing/MB_PD1_Pbgoing_part0.txt","HiForest_pPb_Data_MB_PD1_Pbgoing_job_0",0,0)' # This will not run mixing 
