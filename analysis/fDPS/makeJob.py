import sys
import os
import glob
import string
import subprocess
import ROOT
from commands import getoutput


workpath = os.environ['CMSSW_BASE']+"/src/PicoFramework/analysis/fDPS/"
Job_path = workpath+"Job/"
cross_section_txt = 'the Jpsi pair production cross section is'
cross_section_error_txt = 'the Jpsi pair production cross section error is'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-s', '--sampletype', dest='sampletyple', action='store', type=str, default='DPS')
parser.add_argument('-w', '--whichJpsi', dest='whichJpsi', action='store', type=str, default='Jpsi1')
parser.add_argument('-c', '--cutcategory', dest='cutcategory', action='store', type=str, default='JpsiPair_mass')

args = parser.parse_args()
whichJpsi = args.whichJpsi
cutcategory = args.cutcategory
sampletype = args.sampletyple

if not (sampletype=='DPS' or sampletype=='SPS' or sampletype=='data' or sampletype=='All'):
   raise ValueError('"sampletype" must be "DPS" or "SPS" or "data" or "All"')

if cutcategory=='JpsiPair_mass':
   bins = [5.,10.,15.,20.,25.,30.]
elif cutcategory=='JpsiPair_DeltaEta':
   bins = [0.,0.5,1.0,1.5,2.0,2.5]
else:
   raise ValueError('"cutcategory" must be"JpsiPair_mass" or "JpsiPair_DetaEta" ')

h_CS_DPS = ROOT.TH1F("cross section of DPS", "h_CS_DPS",len(bins)-1, bins[0], bins[len(bins)-1])
h_CS_SPS = ROOT.TH1F("cross section of SPS", "h_CS_SPS",len(bins)-1, bins[0], bins[len(bins)-1])
h_CS_data = ROOT.TH1F("cross section of data", "h_CS_data",len(bins)-1, bins[0], bins[len(bins)-1])

if (sampletype=='DPS' or sampletype=='SPS'):
 codeFile_pre = sampletype+'_'+whichJpsi+'_'+cutcategory+'.sh' 
elif (sampletype=='data'):
 codeFile_pre = sampletype+'_'+cutcategory+'.sh' 

rootfilename = 'cross_section_'+cutcategory+'.root'
SF =1000


def main():
 if (sampletype == 'DPS' or  sampletype == 'SPS' or  sampletype == 'data'):
  os.chdir(Job_path)
  codeFile = file(codeFile_pre,"w")
 for i in range(len(bins)-1):
  minbin=bins[i]
  maxbin=bins[i+1]
  CS=[]
  CS_err=[]
  if (sampletype == 'DPS' or  sampletype == 'SPS' or  sampletype == 'data'):
   makecommand(workpath, codeFile, sampletype, whichJpsi, cutcategory, minbin, maxbin)
  if sampletype == 'All':  
   # read cross sestion in each bin
   readCS(workpath, 'DPS', cutcategory, minbin, maxbin, cross_section_txt, cross_section_error_txt, CS, CS_err)
   readCS(workpath, 'SPS', cutcategory, minbin, maxbin, cross_section_txt, cross_section_error_txt, CS, CS_err)
   readCS(workpath, 'data', cutcategory, minbin, maxbin, cross_section_txt, cross_section_error_txt, CS, CS_err)
   print "JpsiPair cross section of DPS in %s bin %f-%f is %f +- %f" % (cutcategory,minbin,maxbin,CS[0],CS_err[0])
   print "JpsiPair cross section of SPS in %s bin %f-%f is %f +- %f" % (cutcategory,minbin,maxbin,CS[1],CS_err[1])
   print "JpsiPair cross section of data in %s bin %f-%f is %f +- %f" % (cutcategory,minbin,maxbin,CS[2],CS_err[2])
   print ""
   h_CS_DPS.SetBinContent(i+1,SF*CS[0]) 
   h_CS_DPS.SetBinError(i+1,SF*CS_err[0]) 

   h_CS_SPS.SetBinContent(i+1,SF*CS[1]) 
   h_CS_SPS.SetBinError(i+1,SF*CS_err[1]) 

   h_CS_data.SetBinContent(i+1,SF*CS[2]) 
   h_CS_data.SetBinError(i+1,SF*CS_err[2]) 
 if (sampletype == 'DPS' or  sampletype == 'SPS' or  sampletype == 'data'):
  os.popen('chmod +x %s' % (codeFile_pre))
 os.chdir(workpath)
 if sampletype == 'All':  
  Norm_DPS = 0.9*h_CS_data.Integral()/h_CS_DPS.Integral()
  Norm_SPS = 0.1*h_CS_data.Integral()/h_CS_SPS.Integral()
  h_CS_DPS.Sumw2()
  h_CS_DPS.Scale(Norm_DPS)
  h_CS_SPS.Sumw2()
  h_CS_SPS.Scale(Norm_SPS)
  Workspace = ROOT.TFile.Open(rootfilename, "RECREATE")
  h_CS_DPS.SetName("DPS")
  h_CS_SPS.SetName("SPS")
  h_CS_data.SetName("data_obs")
  h_CS_DPS.Write() 
  h_CS_SPS.Write() 
  h_CS_data.Write() 

 
 #path_at = (os.popen("pwd").read()).strip()
 #print path_at
 #os.system('sh %s/%s' % (Job_path, codeFile_pre)) 
 #subprocess.Popen(['sh %s/%s' % (Job_path, codeFile_pre)], shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)  
 #print output

def makecommand(workpath, codeFile, sampletype, whichJpsi, cutcategory, minbin, maxbin):
 if (sampletype == 'DPS' or sampletype == 'SPS'): 
  print >> codeFile, "root -b -q %s/RooFit_%s_Jpsi.cc'(\"%s\",\"%s\",%f,%f)'" % (workpath, sampletype, whichJpsi, cutcategory, minbin, maxbin)
 elif sampletype == 'data':
  print >> codeFile, "root -b -q %s/RooFit2D_fDPS_%s.cc'(\"%s\",%f,%f)'" % (workpath, sampletype, cutcategory, minbin, maxbin)
  #os.system("root -b -q %s/RooFit_%s_Jpsi.cc'(\"%s\",\"%s\",%f,%f)'" % (workpath, sampletype, whichJpsi, cutcategory, minbin, maxbin))
 
def readCS(workpath, sampletype, cutcategory, minbin, maxbin, cross_section_txt, cross_section_error_txt, CS, CS_err):
 if (sampletype == 'data'): 
  os.system("root -b -q %s/RooFit2D_fDPS_%s.cc'(\"%s\",%f,%f)'>%s%s_%s.txt"  % (workpath, sampletype, cutcategory, minbin, maxbin, workpath, sampletype, cutcategory))
 if (sampletype == 'DPS' or sampletype == 'SPS'): 
  os.system("root -b -q %s/RooFit_%s_Jpsi.cc'(\"Jpsi1\",\"%s\",%f,%f)'>%s%s_%s.txt"  % (workpath, sampletype, cutcategory, minbin, maxbin, workpath, sampletype, cutcategory))
 cross_section = float(getoutput("cat %s%s_%s.txt | grep \"%s\" | cut -d ' ' -f 8" % (workpath, sampletype, cutcategory, cross_section_txt)))  
 cross_section_err = float(getoutput("cat %s%s_%s.txt | grep \"%s\" | cut -d ' ' -f 9" % (workpath, sampletype, cutcategory, cross_section_error_txt)))  
  #print cross_section
  #print cross_section_err
 CS.append(cross_section)
 CS_err.append(cross_section_err)


def Saveworkspace():
  pass    

 
if __name__=='__main__':
 main()


