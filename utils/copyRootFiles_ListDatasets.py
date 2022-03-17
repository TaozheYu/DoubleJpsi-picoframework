#!/usr/bin/env python
from commands import getoutput
import os, sys
copy_path = "/eos/user/t/tayu/DoubleJpsi/"
path_fw = os.environ['CMSSW_BASE']+"/src/PicoFramework/selection/"
path_utils = path_fw+"utils/"
sys.path.append(path_utils)
passwd = os.popen('cat /afs/cern.ch/user/t/tayu/private/vomsproxy').read().strip()
os.system('echo %s | voms-proxy-init --valid 192:00 -voms cms -rfc' % (passwd)) #192 means keep the credentials valid for 192h
from datasets_DoubleJpsi_info import *
#User input
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-ys',  '--years',       dest='years',       action='store', type=str, default='2016_2017_2018') #The years you want to process (separated by _!)
parser.add_argument('-df',  '--datasetFlag', dest='datasetFlag', action='store', type=str, default='All') #See datasets_info[d][1] in datasets_info.py. "All" stands for all he cases.
parser.add_argument('-i',   '--iterations',  dest='iterations',  action='store', type=int, default='10')
args = parser.parse_args()
years       = args.years
datasetFlag = args.datasetFlag
iterations  = args.iterations

def copyRootFiles(fileName,copy_path,year,copy_files):
 with open(fileName+".txt",'r') as file_to_read:
   while True:
      line = file_to_read.readline().split("\n")[0]
      if not line:
         break
         pass 
      print >> copy_files, "gfal-copy root://cms-xrd-global.cern.ch/"+line+" "+copy_path+year+'/'+fileName

#Main
if not os.path.exists("copy_rootFilesList"): os.makedirs("copy_rootFilesList")
if not os.path.exists("rootFilesList"): os.makedirs("rootFilesList")
os.chdir("copy_rootFilesList")
work_path = (os.popen("pwd").read()).strip() #.strip() removes empty characters
os.chdir("../rootFilesList")
path = (os.popen("pwd").read()).strip() #.strip() removes empty characters
#print work_path
num_datasets_of_interest = 0
finishedDatasets = []
keep_checking_count = 0
while(keep_checking_count<iterations):
 keep_checking_count = keep_checking_count+1
 print 'Iteration %s' % (keep_checking_count)
 for d in range(0,len(datasets_info)):
  os.chdir(path)
  year = str(datasets_info[d][0])
  if not year in years: continue
  if not (datasetFlag==datasets_info[d][1] or datasetFlag=="All"): continue
  if not os.path.exists(year): os.makedirs(year) 
  os.chdir(work_path)
  if not os.path.exists(year): os.makedirs(year)
  os.chdir(copy_path)
  if not os.path.exists(year): os.makedirs(year) 
  os.chdir(path+"/"+year)
  #os.chdir(year)
  fileName = datasets_info[d][2]
  os.chdir(copy_path+year)
  if not os.path.exists(fileName): os.makedirs(fileName)
  os.chdir(path+"/"+year)
  if os.path.isfile(datasets_info[d][2]+".txt"): continue #Helpful if you run this script more than one time and it had previously finished with missing processed datasets
  if keep_checking_count==1: num_datasets_of_interest = num_datasets_of_interest+1
  if fileName in finishedDatasets: continue
  files = getoutput('dasgoclient --query="file dataset=%s"' % (datasets_info[d][3]))
  #os.system("sleep 1s")
  rootFilesListFile = file(fileName+".txt","w")
  print >> rootFilesListFile, "%s" % (files)
  rootFilesListFile.close()
  rootFilesListFile = file(fileName+".txt","r")
  numRootFiles = int(getoutput('cat %s.txt | grep root | wc -l' % (fileName)))
  rootFilesListFile.close()
  copy_files = file(fileName+".sh","w")
  copyRootFiles(fileName,copy_path,year,copy_files)
  os.popen('mv %s/%s/%s.sh %s/%s' % (path,year,fileName,work_path,year))
  if numRootFiles>=1: #It assumes that if dasgoclient got 1 root file it got all of them (and there are datasets produced with just 1 file)
   finishedDatasets.append(fileName)
  else:
   os.system("rm %s.txt" % (fileName))
   if keep_checking_count==iterations: print '%s %s is empty!' % (year,fileName)
 if keep_checking_count>1 and len(finishedDatasets)==num_datasets_of_interest: break
print "Ended at iteration %s. If Iteration is 10, you must run it again" % (keep_checking_count)

  




