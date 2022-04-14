#!/usr/bin/env python3
import os, sys
import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from search_variables_DoubleJpsi import *
from search_functions_DoubleJpsi import *

class declareVariables(search_variables):
 def __init__(self, name):
  super(declareVariables, self).__init__(name)

class Producer(Module):
 def __init__(self, **kwargs):
  #User inputs
  self.channel     = kwargs.get('channel') 
  self.isData      = kwargs.get('dataType')=='data'
  self.year        = kwargs.get('year') 
  self.maxNumEvt   = kwargs.get('maxNumEvt')
  self.prescaleEvt = kwargs.get('prescaleEvt')
  self.lumiWeight  = kwargs.get('lumiWeight')

  #Analysis quantities
  if self.year==2018:
   #Trigger 
   print "not yet"
  elif self.year==2017:
   #Trigger
   if self.channel=="mumu":
    self.trigger = lambda e: e.HLT_Dimuon0_Jpsi_Muon

  elif self.year==2016:
   #Trigger
   if self.channel=="mumu":
    self.trigger = lambda e: e.HLT_Dimuon0_Jpsi_Muon

  else:
   raise ValueError('Year must be above 2016 (included).')

  #ID
  
  #Corrections

  #Cut flow table
        
 def beginJob(self):
  print "Here is beginJob"
  #pass
        
 def endJob(self):
  print "Here is endJob"
  #pass
        
 def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
  print "Here is beginFile"
  self.sumNumEvt = 0
  self.sumgenWeight = 0
  self.out = declareVariables(inputFile)
  self.correct_file = ROOT.TFile.Open("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/Scale_factor/corrector.root")
  self.h_w_acc = self.correct_file.Get("h_w_acc")
  self.h_w_reco= self.correct_file.Get("h_w_reco")
  self.h_w_eff = self.correct_file.Get("h_w_eff")
  self.h_w_trig= self.correct_file.Get("h_w_trig")
  #pass
        
 def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):        
  print "Here is endFile"
  self.out.sumNumEvt[0] = self.sumNumEvt
  self.out.sumgenWeight[0] = self.sumgenWeight
  self.out.evtree.Fill()
  self.out.outputfile.Write()
  self.out.outputfile.Close()
  #pass
        
 def analyze(self, event):
  """process event, return True (go to next module) or False (fail, go to next event)"""
  #For all events
  if(self.sumNumEvt>self.maxNumEvt and self.maxNumEvt!=-1): return False
  self.sumNumEvt = self.sumNumEvt+1
  if not self.isData: self.sumgenWeight = self.sumgenWeight+(event.genWeight/abs(event.genWeight))
  if not self.sumNumEvt%self.prescaleEvt==0: return False
  passCut = 0 

  #Primary vertex (loose PV selection)
  '''
  if not event.PV_npvs>0:
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 1
  '''
  #Trigger        
  if not self.trigger(event):
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 2  
  #Muons
  #print "run:lumi:evt %s:%s:%s" % (event.run,event.luminosityBlock,event.event)
  muons = Collection(event, 'Muon')
  selectedMusIdx = []
  selectMus(event,selectedMusIdx) 
  #print len(selectedMusIdx)
  if len(selectedMusIdx)<4:
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 3
  #print selectedMusIdx
  
  #Reconstruct Jpsi by 2 muons
  Jpsi_muindex1=[]
  Jpsi_muindex2=[]
  selectedJpsi=[]
  selectJpsi(event, selectedMusIdx, Jpsi_muindex1, Jpsi_muindex2, selectedJpsi)
  if len(selectedJpsi)<2:
    return False
  #print Jpsi_muindex1
  #print Jpsi_muindex2
  #print selectedJpsi[0].Pt()

  #print "dimuon1_mu1:%s dimuon1_mu2:%s dimuon2_mu1:%s dimuon2_mu2:%s" % (muons[Jpsi_muindex1[0]].pt, muons[Jpsi_muindex1[1]].pt, muons[Jpsi_muindex2[0]].pt, muons[Jpsi_muindex2[1]].pt)

  #select Jpsi pair
  fourmuon = []
  JpsiPair = []
  selectJpsiPair(event, Jpsi_muindex1, Jpsi_muindex2, selectedJpsi, fourmuon, JpsiPair) 
  if (len(JpsiPair)<2 or len(fourmuon)<4):
    return False
  Jpsi1 = JpsiPair[0]
  Jpsi2 = JpsiPair[1]
  mu1 = fourmuon[0] 
  mu2 = fourmuon[1] 
  mu3 = fourmuon[2] 
  mu4 = fourmuon[3]


  # add corrector in Jpsi (pt,eta) bin
  Jpsi1_corr_list=[] #0=acc, 1=reco, 2=eff, 
  Jpsi1_corr_list_up=[] #0=acc, 1=reco, 2=eff, 
  Jpsi1_corr_list_do=[] #0=acc, 1=reco, 2=eff, 

  Jpsi2_corr_list=[] #0=acc, 1=reco, 2=eff, 3=trig
  Jpsi2_corr_list_up=[] #0=acc, 1=reco, 2=eff, 3=trig
  Jpsi2_corr_list_do=[] #0=acc, 1=reco, 2=eff, 3=trig

  AddCorrector(JpsiPair, Jpsi1_corr_list, Jpsi1_corr_list_up, Jpsi1_corr_list_do, Jpsi2_corr_list, Jpsi2_corr_list_up, Jpsi2_corr_list_do, self.h_w_acc, self.h_w_reco, self.h_w_eff, self.h_w_trig)
   
  #print "mu1:%s mu2:%s mu3:%s mu4:%s" % (mu1.Pt(),mu2.Pt(),mu3.Pt(),mu4.Pt())
  '''
  mu0 = muons[selectedMusIdx[0]].p4()
  mu1 = muons[selectedMusIdx[1]].p4() 
  if not mu0.Pt()>=30:
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 4
  if not mu0.Pt()-mu1.Pt()>=30: 
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 5
  if not mu0.DeltaR(mu1)>=0.3:
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 6
  if not (mu0+mu1).M()>=200:
   effevt(passCut,self,event)
   return False
  passCut = passCut+1 #passCut = 7
  effevt(passCut,self,event) 
  '''
  #print "mu0 idx pT |eta| ID Iso %s %s %s %s %s" % (selectedMusIdx[0],event.Muon_pt[selectedMusIdx[0]],abs(event.Muon_eta[selectedMusIdx[0]]),event.Muon_tightId[selectedMusIdx[0]],event.Muon_pfRelIso04_all[selectedMusIdx[0]])
  #print "mu1 idx pT |eta| ID Iso %s %s %s %s %s" % (selectedMusIdx[1],event.Muon_pt[selectedMusIdx[1]],abs(event.Muon_eta[selectedMusIdx[1]]),event.Muon_tightId[selectedMusIdx[1]],event.Muon_pfRelIso04_all[selectedMusIdx[1]])

  print "run:lumi:evt %s:%s:%s" % (event.run,event.luminosityBlock,event.event)
  if not self.isData: print "genWeights event.Pileup_nTrueInt %s %s" % (event.genWeight/abs(event.genWeight), event.Pileup_nTrueInt) 
  print ""

  #Event
  self.out.run[0] = event.run
  self.out.luminosityBlock[0] = event.luminosityBlock
  self.out.event[0] = event.event #& 0xffffffffffffffff
  self.out.Jpsi1_pt[0]   = Jpsi1.Pt()
  self.out.Jpsi1_eta[0]  = Jpsi1.Eta()
  self.out.Jpsi1_phi[0]  = Jpsi1.Phi()
  self.out.Jpsi1_mass[0] = Jpsi1.M()
   
  self.out.Jpsi2_pt[0]   = Jpsi2.Pt()  
  self.out.Jpsi2_eta[0]  = Jpsi2.Eta()
  self.out.Jpsi2_phi[0]  = Jpsi2.Phi()
  self.out.Jpsi2_mass[0] = Jpsi2.M()   

  self.out.muon1_pt[0]  = mu1.Pt() 
  self.out.muon1_eta[0] = mu1.Eta()
  self.out.muon1_phi[0] = mu1.Phi()

  self.out.muon2_pt[0]  = mu2.Pt() 
  self.out.muon2_eta[0] = mu2.Eta()
  self.out.muon2_phi[0] = mu2.Phi()

  self.out.muon3_pt[0]  = mu3.Pt() 
  self.out.muon3_eta[0] = mu3.Eta()
  self.out.muon3_phi[0] = mu3.Phi()

  self.out.muon4_pt[0]  = mu4.Pt() 
  self.out.muon4_eta[0] = mu4.Eta()
  self.out.muon4_phi[0] = mu4.Phi()
 
  self.out.w_acc_Jpsi1[0] = Jpsi1_corr_list[0]
  self.out.w_reco_Jpsi1[0] = Jpsi1_corr_list[1]
  self.out.w_eff_Jpsi1[0] = Jpsi1_corr_list[2]
 
  self.out.w_acc_Jpsi1_up[0] = Jpsi1_corr_list_up[0]
  self.out.w_reco_Jpsi1_up[0] = Jpsi1_corr_list_up[1]
  self.out.w_eff_Jpsi1_up[0] = Jpsi1_corr_list_up[2]
 
  self.out.w_acc_Jpsi1_do[0] = Jpsi1_corr_list_do[0]
  self.out.w_reco_Jpsi1_do[0] = Jpsi1_corr_list_do[1]
  self.out.w_eff_Jpsi1_do[0] = Jpsi1_corr_list_do[2]
 
  self.out.w_acc_Jpsi2[0] = Jpsi2_corr_list[0]
  self.out.w_reco_Jpsi2[0] = Jpsi2_corr_list[1]
  self.out.w_eff_Jpsi2[0] = Jpsi2_corr_list[2]
 
  self.out.w_acc_Jpsi2_up[0] = Jpsi2_corr_list_up[0]
  self.out.w_reco_Jpsi2_up[0] = Jpsi2_corr_list_up[1]
  self.out.w_eff_Jpsi2_up[0] = Jpsi2_corr_list_up[2]
 
  self.out.w_acc_Jpsi2_do[0] = Jpsi2_corr_list_do[0]
  self.out.w_reco_Jpsi2_do[0] = Jpsi2_corr_list_do[1]
  self.out.w_eff_Jpsi2_do[0] = Jpsi2_corr_list_do[2]

  self.out.w_trig_Jpsi12[0] = Jpsi2_corr_list[3]
  self.out.w_trig_Jpsi12_up[0] = Jpsi2_corr_list_up[3]
  self.out.w_trig_Jpsi12_do[0] = Jpsi2_corr_list_do[3]
  #Save tree
  self.out.Events.Fill() 
  return True
