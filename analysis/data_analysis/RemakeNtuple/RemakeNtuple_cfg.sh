#!/bin/bash

echo "Starting job on " `date` #Date/time of start of job
echo "Runing on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
#export X509_USER_PROXY=/afs/cern.ch/user/t/tayu/x509up_u93266
cd /afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/analysis/data_analysis/RemakeNtuple  
#cmsenv
eval `scram runtime -sh`

root -b -q RemakeNtuple.cc

 

