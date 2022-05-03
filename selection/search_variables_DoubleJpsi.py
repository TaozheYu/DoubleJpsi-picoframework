#!/usr/bin/env python3
import os, sys
import ROOT
import math 
import numpy as np 

class search_variables(object):
 def __init__(self, name):

  #Create file 
  inputFile = name
  outputFileName = os.path.basename(str(inputFile)).split(".root", 1)[0]+"_Skim.root"
  compression = "LZMA:9"
  ROOT.gInterpreter.ProcessLine("#include <Compression.h>")
  (algo, level) = compression.split(":")
  compressionLevel = int(level)
  if   algo == "LZMA": compressionAlgo  = ROOT.ROOT.kLZMA
  elif algo == "ZLIB": compressionAlgo  = ROOT.ROOT.kZLIB
  else: raise RuntimeError("Unsupported compression %s" % algo)
  self.outputfile = ROOT.TFile(outputFileName, 'RECREATE',"",compressionLevel)
  self.outputfile.SetCompressionAlgorithm(compressionAlgo)

  #All entries 
  self.evtree = ROOT.TTree('evtree','evtree')
  self.add_float(self.evtree,"sumNumEvt")
  self.add_float(self.evtree,"sumgenWeight")

  #Efficiency/Num events
  self.effevt = ROOT.TTree('effevt','effevt')
  self.add_float(self.effevt,"passCut")
  self.add_float(self.effevt,"eegenWeight")
  self.add_float(self.effevt,"eelumiWeight")
        
  #Common variables 
  self.Events = ROOT.TTree('Events','Events')
  self.add_float(self.Events,"run")
  self.add_float(self.Events,"luminosityBlock")
  self.add_float(self.Events,"event")
  self.add_float(self.Events,"genWeight")
  self.add_float(self.Events,"lumiWeight")

  self.add_float(self.Events,"Jpsi1_pt")
  self.add_float(self.Events,"Jpsi1_eta")
  self.add_float(self.Events,"Jpsi1_phi")
  self.add_float(self.Events,"Jpsi1_mass")

  self.add_float(self.Events,"Jpsi2_pt")
  self.add_float(self.Events,"Jpsi2_eta")
  self.add_float(self.Events,"Jpsi2_phi")
  self.add_float(self.Events,"Jpsi2_mass")

  self.add_float(self.Events,"JpsiPair_mass")
  self.add_float(self.Events,"JpsiPair_DeltaEta")

  self.add_float(self.Events,"muon1_pt")
  self.add_float(self.Events,"muon1_eta")
  self.add_float(self.Events,"muon1_phi")

  self.add_float(self.Events,"muon2_pt")
  self.add_float(self.Events,"muon2_eta")
  self.add_float(self.Events,"muon2_phi")

  self.add_float(self.Events,"muon3_pt")
  self.add_float(self.Events,"muon3_eta")
  self.add_float(self.Events,"muon3_phi")

  self.add_float(self.Events,"muon4_pt")
  self.add_float(self.Events,"muon4_eta")
  self.add_float(self.Events,"muon4_phi")

  self.add_float(self.Events,"w_acc_Jpsi1")
  self.add_float(self.Events,"w_reco_Jpsi1")
  self.add_float(self.Events,"w_eff_Jpsi1")

  self.add_float(self.Events,"w_acc_Jpsi1_up")
  self.add_float(self.Events,"w_reco_Jpsi1_up")
  self.add_float(self.Events,"w_eff_Jpsi1_up")

  self.add_float(self.Events,"w_acc_Jpsi1_do")
  self.add_float(self.Events,"w_reco_Jpsi1_do")
  self.add_float(self.Events,"w_eff_Jpsi1_do")

  self.add_float(self.Events,"w_acc_Jpsi2")
  self.add_float(self.Events,"w_reco_Jpsi2")
  self.add_float(self.Events,"w_eff_Jpsi2")

  self.add_float(self.Events,"w_acc_Jpsi2_up")
  self.add_float(self.Events,"w_reco_Jpsi2_up")
  self.add_float(self.Events,"w_eff_Jpsi2_up")

  self.add_float(self.Events,"w_acc_Jpsi2_do")
  self.add_float(self.Events,"w_reco_Jpsi2_do")
  self.add_float(self.Events,"w_eff_Jpsi2_do")
  
  self.add_float(self.Events,"w_trig_Jpsi12")
  self.add_float(self.Events,"w_trig_Jpsi12_up")
  self.add_float(self.Events,"w_trig_Jpsi12_do")

 def add_float(self,tree,name,dtype=np.dtype(float)):
  if hasattr(self,name):
   print('ERROR! SetBranchAddress of name "%s" already exists!' % (name))
   exit(1)
  setattr(self,name,np.full((1),-99999999999999999999999999999999999999999999999999,dtype=dtype)) #1 elem w/ inizialization '-99999999999999999999999999999999999999999999999999'
  tree.Branch(name,getattr(self,name),'{0}/D'.format(name)) #The types in root (/D in this example) are defined here https://root.cern/root/html528/TTree.html

 def add_string(self,tree,name,dtype=np.dtype('S100')): #It assumes a string of max 100 characters
  if hasattr(self,name):
   print('ERROR! SetBranchAddress of name "%s" already exists!' % (name))
   exit(1)
  setattr(self,name,np.full((1),'noVal',dtype=dtype)) #1 elem w/ inizialization 'noval'
  tree.Branch(name,getattr(self,name),'{0}/C'.format(name))

 def add_vectorFloat(self,tree,name):
  if hasattr(self,name):
   print('ERROR! SetBranchAddress of name "%s" already exists!' % (name))
   exit(1)
  setattr(self,name,ROOT.std.vector('float')()) #No inizialization
  tree.Branch(name,getattr(self,name))

 def add_vectorString(self,tree,name):
  if hasattr(self,name):
   print('ERROR! SetBranchAddress of name "%s" already exists!' % (name))
   exit(1)
  setattr(self,name,ROOT.std.vector('string')()) #No inizialization
  tree.Branch(name,getattr(self,name))
