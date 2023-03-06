#!/bin/bash
/bin/hostname
gcc -v
pwd
cd /workfs2/cms/yutz/CMSSW_10_2_9_Trial/src/DoubleJpsi-picoframework/analysis/data_analysis/RemakeNtuple

root -l -b -q RemakeNtuple.cc

 
