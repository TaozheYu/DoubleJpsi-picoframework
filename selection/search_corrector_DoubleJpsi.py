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
    N = 10 
    xbins= [5.,6.,7.,8.,10.,12.,15.,20.,30.,40.]
    M = 19
    ybins = [-2.2,-2.,-1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.,2.2]
    self.h_total_Jpsi = ROOT.TH2F("h_total_Jpsi","h_total_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_acc_Jpsi = ROOT.TH2F("h_acc_Jpsi","h_acc_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_reco_Jpsi = ROOT.TH2F("h_reco_Jpsi","h_reco_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_eff_Jpsi = ROOT.TH2F("h_eff_Jpsi","h_eff_Jpsi",N-1 , array('d',xbins), M-1, array('d',ybins)) 
    self.h_tri_Jpsi = ROOT.TH2F("h_tri_Jpsi","h_tri_Jpsi",N-1 , array('d',xbins), N-1, array('d',xbins)) 
    self.addObject(self.h_total_Jpsi)
    self.addObject(self.h_acc_Jpsi)
    self.addObject(self.h_reco_Jpsi)
    self.addObject(self.h_eff_Jpsi)

  def analyze(self, event):
    """process event, return True (go to next module) or False (fail, go to next event)"""
    if self.sumNumEvt>100000: return False
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
              #if not Jpsi.Pt()>=3.5: continue
              if not Jpsi.Eta()<=2.2: continue
              if not abs(Jpsi.M()-3.092)<=0.3: continue              
              self.h_eff_Jpsi.Fill(event.GenPart_pt[igen],event.GenPart_eta[igen])
    return True
   
    '''
    for igen in range(event.nGenPart):
      if event.GenPart_pdgId[igen]==443:
        for jgen in range(event.nGenPart):
          if event.GenPart_pdgId[jgen]==13 and event.GenPart_genPartIdxMother[jgen]==igen:
            for kgen in range(event.nGenPart):
              if event.GenPart_pdgId[kgen]==-13 and event.GenPart_genPartIdxMother[kgen]==igen:
                if abs(event.GenPart_eta[kgen])<=2.4 and abs(event.GenPart_eta[jgen])<=2.4:
                  for imu in range(event.Muon):
                    for jmu in range(event.Muon):
                      if (Muon_genPartIdx[imu]==jgen and Muon_genPartIdx[jmu]==kgen) or (Muon_genPartIdx[jmu]==jgen)
    ''' 
 

files = ["/eos/user/t/tayu/DoubleJpsi/data_and_MC/2016/DPS/DPS.root"]
p = PostProcessor(".", files, branchsel=None, modules=[
                  ExampleAnalysis()], noOut=True, histFileName="histOut.root", histDirName="plots") 
p.run()
