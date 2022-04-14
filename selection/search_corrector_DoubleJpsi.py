#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True


class ExampleAnalysis(Module):
  def __init__(self):
    self.writeHistFile = True

  def beginJob(self, histFile=None, histDirName=None):
    Module.beginJob(self, histFile, histDirName)
    self.sumNumEvt = 0
    xbins= [5.,6.,7.,8.,9.,10.,12.,15.,20.,25.,40.]
    ybins = [-2.2,-2.,-1.75,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,1.75,2.,2.2]
    pt_bins = [5.,7.,9.,12.,15.,25.,40.] 

    N = len(xbins)
    M = len(ybins)
    L = len(pt_bins)

    self.h_total_Jpsi = ROOT.TH2F("h_total_Jpsi","h_total_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_acc_Jpsi = ROOT.TH2F("h_acc_Jpsi","h_acc_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_reco_Jpsi = ROOT.TH2F("h_reco_Jpsi","h_reco_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_eff_Jpsi = ROOT.TH2F("h_eff_Jpsi","h_eff_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 

    self.h_betri_Jpsi = ROOT.TH2F("h_betri_Jpsi","h_betri_Jpsi",L-1 , array('d',pt_bins), L-1, array('d',pt_bins)) 
    self.h_patri_Jpsi = ROOT.TH2F("h_patri_Jpsi","h_patri_Jpsi",L-1 , array('d',pt_bins), L-1, array('d',pt_bins)) 
    self.addObject(self.h_total_Jpsi)
    self.addObject(self.h_acc_Jpsi)
    self.addObject(self.h_reco_Jpsi)
    self.addObject(self.h_eff_Jpsi)
    self.addObject(self.h_betri_Jpsi)
    self.addObject(self.h_patri_Jpsi)

  def analyze(self, event):
    """process event, return True (go to next module) or False (fail, go to next event)"""
    if self.sumNumEvt>10000000: return False
    self.sumNumEvt = self.sumNumEvt+1
    gJpsi_pt = -999
    gJpsi_eta= -999
    genparts = Collection(event, 'GenPart')
    muons = Collection(event, 'Muon')
    ''' 
    for igen in range(event.nGenPart):
      if abs(event.GenPart_pdgId[igen])==13:
        gmotherID = event.GenPart_genPartIdxMother[igen]
        if gmotherID>0:
	  if event.GenPart_pdgId[gmotherID] == 443 and event.GenPart_eta[gmotherID]<=2:
	    gJpsi_pt = event.GenPart_pt[gmotherID]
	    gJpsi_eta= event.GenPart_eta[gmotherID]
	    #if (0.<gJpsi_pt<=1.0) and (-2<gJpsi_eta<=-1.75):
	    self.h_total_Jpsi.Fill(gJpsi_pt,gJpsi_eta)
	  else:
	    continue  
    return True
    '''
    selectedJpsi_pt = []
    JpsiPair = []
    for igen in range(event.nGenPart):
      if not (event.GenPart_pdgId[igen]==443 and event.GenPart_eta[igen]<=2.2): continue
      for jgen in range(igen+1, event.nGenPart):
        if not (abs(event.GenPart_pdgId[jgen])==13 and event.GenPart_genPartIdxMother[jgen]==igen): continue
        for kgen in range(jgen+1, event.nGenPart):
          if not (abs(event.GenPart_pdgId[kgen])==13 and event.GenPart_genPartIdxMother[kgen]==igen): continue
          self.h_total_Jpsi.Fill(event.GenPart_pt[igen],event.GenPart_eta[igen])
          if not (abs(event.GenPart_eta[kgen])<=2.4 and abs(event.GenPart_eta[jgen])<=2.4): continue
          self.h_acc_Jpsi.Fill(event.GenPart_pt[igen],event.GenPart_eta[igen])
          for imu in range(event.nMuon):
            for jmu in range(imu+1, event.nMuon):
              imuIdx = event.Muon_genPartIdx[imu]
              jmuIdx = event.Muon_genPartIdx[jmu]
              if not ((imuIdx==jgen and jmuIdx==kgen) or (imuIdx==kgen and jmuIdx==jgen)): continue
              self.h_reco_Jpsi.Fill(event.GenPart_pt[igen],event.GenPart_eta[igen])
              if not (event.Muon_mediumId[imu] and event.Muon_mediumId[jmu]): continue 
              if not muons[imu].charge+muons[jmu].charge==0: continue
              mui = muons[imu].p4()
              muj = muons[jmu].p4() 
              Jpsi= mui+muj
              if not Jpsi.Pt()>=5: continue
              if not abs(Jpsi.Eta())<=2.2: continue
              if abs(Jpsi.Eta())<=1:
                if not Jpsi.Pt()>=6: continue
              if not abs(Jpsi.M()-3.092)<=0.3: continue              
              self.h_eff_Jpsi.Fill(Jpsi.Pt(),Jpsi.Eta())
              selectedJpsi_pt.append(Jpsi.Pt())

    
    if len(selectedJpsi_pt)<2: return False
    #print len(selectedJpsi_pt)
    if not len(selectedJpsi_pt)==2: return False
    #for iJpsi in range(len(selectedJpsi)):
      #for jJpsi in range(iJpsi+1,len(selectedJpsi)):
    #This aim to assign the Jpsi1 and Jpsi2 randomly
    gRandom = ROOT.TRandom2(0)
    random = gRandom.Uniform(0,1)
    if random >=0.5:
    #if(selectedJpsi_pt[0]>=selectedJpsi_pt[1]):
      Jpsi1_pt= selectedJpsi_pt[0]
      Jpsi2_pt= selectedJpsi_pt[1]
    else:
      Jpsi1_pt= selectedJpsi_pt[1]
      Jpsi2_pt= selectedJpsi_pt[0]
    
    self.h_betri_Jpsi.Fill(Jpsi1_pt,Jpsi2_pt)
    if not event.HLT_Dimuon0_Jpsi_Muon==1: return False
    self.h_patri_Jpsi.Fill(Jpsi1_pt,Jpsi2_pt) 

    return True
   
 

files = ["/eos/user/t/tayu/DoubleJpsi/data_and_MC/2016/DPS/DPS.root"]
p = PostProcessor(".", files, branchsel=None, modules=[
                  ExampleAnalysis()], noOut=True, histFileName="histOut.root", histDirName="plots") 
p.run()
