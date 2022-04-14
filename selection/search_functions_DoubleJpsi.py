#!/usr/bin/env python3
import ROOT
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection


def effevt(passCut,self,event):
 passCut = passCut+1
 self.out.passCut[0] = passCut
 self.out.eelumiWeight[0] = self.lumiWeight
 if self.isData:
  self.out.eegenWeight[0] = 1
 else:
  self.out.eegenWeight[0] = event.genWeight/abs(event.genWeight)
 self.out.effevt.Fill()

def selectMus(event,selectedMusIdx):
 muons = Collection(event, 'Muon') 
 for imu in range(event.nMuon):
  #print "mu idx pT |eta| ID Iso %s %s %s %s %s" % (imu,event.Muon_pt[imu],abs(event.Muon_eta[imu]),event.Muon_tightId[imu],event.Muon_pfRelIso04_all[imu])
  if not event.Muon_pt[imu]>=3: continue
  if not abs(event.Muon_eta[imu])<=2.4: continue
  if not event.Muon_mediumId[imu]: continue
  #if not event.Muon_pfRelIso04_all[imu]<=0.15: continue
  selectedMusIdx.append(imu)

def selectJpsi(event, selectedMusIdx, Jpsi_muindex1, Jpsi_muindex2, selectedJpsi):
   muons = Collection(event, 'Muon')
   Jpsi = ROOT.TLorentzVector(-99,-99,-99,-99)
   mui = ROOT.TLorentzVector(-99,-99,-99,-99)
   muj = ROOT.TLorentzVector(-99,-99,-99,-99)
   for imu in range(len(selectedMusIdx)):
     for jmu in range(imu+1,len(selectedMusIdx)):
       if not muons[selectedMusIdx[imu]].charge + muons[selectedMusIdx[jmu]].charge==0: continue
       mui.SetPtEtaPhiM(muons[selectedMusIdx[imu]].pt,muons[selectedMusIdx[imu]].eta,muons[selectedMusIdx[imu]].phi,muons[selectedMusIdx[imu]].mass)
       muj.SetPtEtaPhiM(muons[selectedMusIdx[jmu]].pt,muons[selectedMusIdx[jmu]].eta,muons[selectedMusIdx[jmu]].phi,muons[selectedMusIdx[jmu]].mass)
       Jpsi = mui+muj
       if not Jpsi.Pt()>=5: continue
       if not abs(Jpsi.Eta())<=2.2: continue
       if abs(Jpsi.Eta())<=1:
         if not Jpsi.Pt()>=6: continue
       if not abs(Jpsi.M()-3.092)<=0.3: continue
       Jpsi_muindex1.append(imu)
       Jpsi_muindex2.append(jmu)
       selectedJpsi.append(Jpsi)        

def selectJpsiPair(event, Jpsi_muindex1, Jpsi_muindex2, selectedJpsi, fourmuon, JpsiPair): 
  muons = Collection(event, 'Muon')
  mu1 = ROOT.TLorentzVector(-99,-99,-99,-99)
  mu2 = ROOT.TLorentzVector(-99,-99,-99,-99)
  mu3 = ROOT.TLorentzVector(-99,-99,-99,-99)
  mu4 = ROOT.TLorentzVector(-99,-99,-99,-99)
  Jpsi1 = ROOT.TLorentzVector(-99,-99,-99,-99)
  Jpsi2 = ROOT.TLorentzVector(-99,-99,-99,-99)
  minkaisquare=9999
  foundJpsiPair = False
  for iJpsi in range(len(selectedJpsi)):
    for jJpsi in range(iJpsi+1,len(selectedJpsi)):
      imu1 = Jpsi_muindex1[iJpsi]
      imu2 = Jpsi_muindex2[iJpsi]
      jmu1 = Jpsi_muindex1[jJpsi] 
      jmu2 = Jpsi_muindex2[jJpsi]
      if (imu1==jmu1 or imu1==jmu2 or imu2==jmu1 or imu2==jmu2): continue
      if not (selectedJpsi[iJpsi].M()-3.092)*(selectedJpsi[iJpsi].M()-3.092)+(selectedJpsi[jJpsi].M()-3.092)*(selectedJpsi[jJpsi].M()-3.092)<minkaisquare: continue 
      #This aim to assign the Jpsi1 and Jpsi2 randomly
      gRandom = ROOT.TRandom2(0)
      random = gRandom.Uniform(0,1)
      if random >=0.5:
      #if(selectedJpsi[iJpsi].Pt()>=selectedJpsi[jJpsi].Pt()):
        Jpsi1= selectedJpsi[iJpsi]
        Jpsi2= selectedJpsi[jJpsi]
        mu1.SetPtEtaPhiM(muons[imu1].pt, muons[imu1].eta, muons[imu1].phi,muons[imu1].mass)
        mu2.SetPtEtaPhiM(muons[imu2].pt, muons[imu2].eta, muons[imu2].phi,muons[imu2].mass)
        mu3.SetPtEtaPhiM(muons[jmu1].pt, muons[jmu1].eta, muons[jmu1].phi,muons[jmu1].mass)
        mu4.SetPtEtaPhiM(muons[jmu2].pt, muons[jmu2].eta, muons[jmu2].phi,muons[jmu2].mass)
        minkaisquare = (selectedJpsi[iJpsi].M()-3.092)*(selectedJpsi[iJpsi].M()-3.092)+(selectedJpsi[jJpsi].M()-3.092)*(selectedJpsi[jJpsi].M()-3.092)      
      else:
        Jpsi1= selectedJpsi[jJpsi]
        Jpsi2= selectedJpsi[iJpsi]
        mu1.SetPtEtaPhiM(muons[jmu1].pt, muons[jmu1].eta, muons[jmu1].phi, muons[jmu1].mass)
        mu2.SetPtEtaPhiM(muons[jmu2].pt, muons[jmu2].eta, muons[jmu2].phi, muons[jmu2].mass)
        mu3.SetPtEtaPhiM(muons[imu1].pt, muons[imu1].eta, muons[imu1].phi, muons[imu1].mass)
        mu4.SetPtEtaPhiM(muons[imu2].pt, muons[imu2].eta, muons[imu2].phi, muons[imu2].mass)
      foundJpsiPair = True
  if foundJpsiPair:
    JpsiPair.append(Jpsi1)
    JpsiPair.append(Jpsi2)
    fourmuon.append(mu1)
    fourmuon.append(mu2)
    fourmuon.append(mu3)
    fourmuon.append(mu4)

def AddCorrector(JpsiPair, Jpsi1_corr_list, Jpsi1_corr_list_up, Jpsi1_corr_list_do, Jpsi2_corr_list, Jpsi2_corr_list_up, Jpsi2_corr_list_do, h_w_acc, h_w_reco, h_w_eff, h_w_trig):
  Jpsi1 = ROOT.TLorentzVector(-99,-99,-99,-99)
  Jpsi2 = ROOT.TLorentzVector(-99,-99,-99,-99)
  Jpsi1 = JpsiPair[0]
  Jpsi2 = JpsiPair[1]
  if 5<Jpsi1.Pt()<40 and abs(Jpsi1.Eta())<2.2: 
    pt_bin_Jpsi1 = h_w_acc.GetXaxis().FindFixBin(Jpsi1.Pt())
    eta_bin_Jpsi1 = h_w_acc.GetYaxis().FindFixBin(Jpsi1.Eta())

    w_acc_Jpsi1 = h_w_acc.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1)
    w_reco_Jpsi1 = h_w_reco.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1)
    w_eff_Jpsi1 = h_w_eff.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1)
 
    w_acc_Jpsi1_up = h_w_acc.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1) + h_w_acc.GetBinError(pt_bin_Jpsi1, eta_bin_Jpsi1)
    w_reco_Jpsi1_up = h_w_reco.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1) + h_w_reco.GetBinError(pt_bin_Jpsi1, eta_bin_Jpsi1)
    w_eff_Jpsi1_up = h_w_eff.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1) + h_w_eff.GetBinError(pt_bin_Jpsi1, eta_bin_Jpsi1)

    w_acc_Jpsi1_do = h_w_acc.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1) - h_w_acc.GetBinError(pt_bin_Jpsi1, eta_bin_Jpsi1)
    w_reco_Jpsi1_do = h_w_reco.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1) - h_w_reco.GetBinError(pt_bin_Jpsi1, eta_bin_Jpsi1)
    w_eff_Jpsi1_do = h_w_eff.GetBinContent(pt_bin_Jpsi1, eta_bin_Jpsi1) - h_w_eff.GetBinError(pt_bin_Jpsi1, eta_bin_Jpsi1)

    if w_acc_Jpsi1==0 :   
       w_acc_Jpsi1=1
       w_acc_Jpsi1_up=1
       w_acc_Jpsi1_do=1

    if w_reco_Jpsi1==0 :
       w_reco_Jpsi1=1
       w_reco_Jpsi1_up=1
       w_reco_Jpsi1_do=1

    if w_eff_Jpsi1==0 : 
       w_eff_Jpsi1=1
       w_eff_Jpsi1_up=1
       w_eff_Jpsi1_do=1
 
    Jpsi1_corr_list.append(w_acc_Jpsi1)
    Jpsi1_corr_list.append(w_reco_Jpsi1)
    Jpsi1_corr_list.append(w_eff_Jpsi1)

    Jpsi1_corr_list_up.append(w_acc_Jpsi1_up)
    Jpsi1_corr_list_up.append(w_reco_Jpsi1_up)
    Jpsi1_corr_list_up.append(w_eff_Jpsi1_up)

    Jpsi1_corr_list_do.append(w_acc_Jpsi1_do)
    Jpsi1_corr_list_do.append(w_reco_Jpsi1_do)
    Jpsi1_corr_list_do.append(w_eff_Jpsi1_do)
  else:
    Jpsi1_corr_list.append(1)
    Jpsi1_corr_list.append(1)
    Jpsi1_corr_list.append(1)

    Jpsi1_corr_list_up.append(1)
    Jpsi1_corr_list_up.append(1)
    Jpsi1_corr_list_up.append(1)

    Jpsi1_corr_list_do.append(1)
    Jpsi1_corr_list_do.append(1)
    Jpsi1_corr_list_do.append(1)
          
  if 5<Jpsi1.Pt()<40 and  5<Jpsi2.Pt()<40 and abs(Jpsi1.Eta())<2.2: 
    pt_bin_Jpsi2 = h_w_acc.GetXaxis().FindFixBin(Jpsi2.Pt())
    eta_bin_Jpsi2 = h_w_acc.GetYaxis().FindFixBin(Jpsi2.Eta())

    w_acc_Jpsi2 = h_w_acc.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_reco_Jpsi2 = h_w_reco.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_eff_Jpsi2 = h_w_eff.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) 
    w_trig_Jpsi12 = h_w_trig.GetBinContent(pt_bin_Jpsi1,pt_bin_Jpsi2)
 
    w_acc_Jpsi2_up = h_w_acc.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) + h_w_acc.GetBinError(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_reco_Jpsi2_up = h_w_reco.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) + h_w_reco.GetBinError(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_eff_Jpsi2_up = h_w_eff.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) + h_w_eff.GetBinError(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_trig_Jpsi12_up = h_w_trig.GetBinContent(pt_bin_Jpsi1,pt_bin_Jpsi2) + h_w_trig.GetBinError(pt_bin_Jpsi1,pt_bin_Jpsi2)

    w_acc_Jpsi2_do = h_w_acc.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) - h_w_acc.GetBinError(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_reco_Jpsi2_do = h_w_reco.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) - h_w_reco.GetBinError(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_eff_Jpsi2_do = h_w_eff.GetBinContent(pt_bin_Jpsi2, eta_bin_Jpsi2) - h_w_eff.GetBinError(pt_bin_Jpsi2, eta_bin_Jpsi2)
    w_trig_Jpsi12_do = h_w_trig.GetBinContent(pt_bin_Jpsi1,pt_bin_Jpsi2) - h_w_trig.GetBinError(pt_bin_Jpsi1,pt_bin_Jpsi2)

    if w_acc_Jpsi2==0 : 
       w_acc_Jpsi2=1
       w_acc_Jpsi2_up=1
       w_acc_Jpsi2_do=1

    if w_reco_Jpsi2==0 : 
       w_reco_Jpsi2=1
       w_reco_Jpsi2_up=1
       w_reco_Jpsi2_do=1

    if w_eff_Jpsi2==0 : 
       w_eff_Jpsi2=1
       w_reco_Jpsi2_up=1
       w_reco_Jpsi2_do=1

    if w_trig_Jpsi12==0 : 
       w_trig_Jpsi12=1
       w_trig_Jpsi12_up=1
       w_trig_Jpsi12_do=1

    Jpsi2_corr_list.append(w_acc_Jpsi2)
    Jpsi2_corr_list.append(w_reco_Jpsi2)
    Jpsi2_corr_list.append(w_eff_Jpsi2)
    Jpsi2_corr_list.append(w_trig_Jpsi12)

    Jpsi2_corr_list_up.append(w_acc_Jpsi2_up)
    Jpsi2_corr_list_up.append(w_reco_Jpsi2_up)
    Jpsi2_corr_list_up.append(w_eff_Jpsi2_up)
    Jpsi2_corr_list_up.append(w_trig_Jpsi12_up)

    Jpsi2_corr_list_do.append(w_acc_Jpsi2_do)
    Jpsi2_corr_list_do.append(w_reco_Jpsi2_do)
    Jpsi2_corr_list_do.append(w_eff_Jpsi2_do)
    Jpsi2_corr_list_do.append(w_trig_Jpsi12_do)
  else:
    Jpsi2_corr_list.append(1)
    Jpsi2_corr_list.append(1)
    Jpsi2_corr_list.append(1)
    Jpsi2_corr_list.append(1)

    Jpsi2_corr_list_up.append(1)
    Jpsi2_corr_list_up.append(1)
    Jpsi2_corr_list_up.append(1)
    Jpsi2_corr_list_up.append(1)

    Jpsi2_corr_list_do.append(1)
    Jpsi2_corr_list_do.append(1)
    Jpsi2_corr_list_do.append(1)
    Jpsi2_corr_list_do.append(1)
    #print Jpsi1.Pt()
    #print pt_bin_Jpsi1 
    #print Jpsi1.Eta()
    #print eta_bin_Jpsi1
'''
def selectJets(jetType,event,selectedWJetsIdx,selectedMusIdx,selectedJetsIdx):
 muons = Collection(event, 'Muon')
 fatjets = Collection(event, 'FatJet')
 jets = Collection(event, 'Jet')
 for ijet in range(event.nJet):
  #print "jets idx pT eta phi ID %s %s %s %s %s" % (ijet,event.Jet_pt[ijet],abs(event.Jet_eta[ijet]),event.Jet_phi[ijet],event.Jet_jetId[ijet])
  #for imu in range(len(selectedMusIdx)):
   #print "dR(jet,mu) %s %s %s" % (ijet,imu,muons[selectedMusIdx[imu]].p4().DeltaR(jets[ijet].p4()))
  #Kinematic
  if not event.Jet_pt[ijet]>=20: continue
  if not abs(event.Jet_eta[ijet])<=5.0: continue
  if not (jetType==1 or jetType==2 or jetType==3):
   if not event.Jet_pt[ijet]>=30: continue
  if not jetType==4:
   if not abs(event.Jet_eta[ijet])<=2.4: continue
  #ID
  if not event.Jet_jetId[ijet]>=2: continue
  #b-jet
  if jetType==2:
   if not event.Jet_btagCSVV2[ijet]>0.8838: continue #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X 
  #Cleaning
  isNotFatJetLep = True
  for ifatjet in range(len(selectedWJetsIdx)):
   if fatjets[selectedWJetsIdx[ifatjet]].p4().DeltaR(jets[ijet].p4())<0.8:
    isNotFatJetLep = False
  if not isNotFatJetLep: continue
  for imu in range(len(selectedMusIdx)):
   if muons[selectedMusIdx[imu]].p4().DeltaR(jets[ijet].p4())<0.4:
    isNotFatJetLep = False
  if not isNotFatJetLep: continue
  selectedJetsIdx.append(ijet)

def selectFatJets(jetType,event,selectedMusIdx,selectedFatJetsIdx):
 muons = Collection(event, 'Muon')
 fatjets = Collection(event, 'FatJet')
 for ifatjet in range(event.nFatJet):
  #print "fat jets idx pT eta phi ID msoftdrop tau2/tau1 %s %s %s %s %s %s %s" % (ifatjet,event.FatJet_pt[ifatjet],abs(event.FatJet_eta[ifatjet]),event.FatJet_phi[ifatjet],event.FatJet_jetId[ifatjet],event.FatJet_msoftdrop[ifatjet],event.FatJet_tau2[ifatjet]/event.FatJet_tau1[ifatjet])
  #Kinematic
  if not event.FatJet_pt[ifatjet]>=180: continue
  if not abs(event.FatJet_eta[ifatjet])<2.4: continue #This is the recommendation for all the fat jets (there are not reconstructed forward fat jets)
  #ID
  if not event.FatJet_jetId[ifatjet]>=2: continue
  if(jetType==0):
   if not (65<event.FatJet_msoftdrop[ifatjet] and event.FatJet_msoftdrop[ifatjet]<105): continue
   if not event.FatJet_tau2[ifatjet]/event.FatJet_tau1[ifatjet]<0.55: continue
  elif(jetType==1):
   if not event.FatJet_pt[ifatjet]>=400: continue
   if not (105<event.FatJet_msoftdrop[ifatjet] and event.FatJet_msoftdrop[ifatjet]<220): continue
   if not event.FatJet_tau3[ifatjet]/event.FatJet_tau2[ifatjet]<0.81: continue
  #Cleaning
  isNotLep = True
  for imu in range(len(selectedMusIdx)):
   #print "dR(fat jet,mu) %s %s" % (selectedMusIdx[imu],muons[selectedMusIdx[imu]].p4().DeltaR(fatjets[ifatjet].p4()))
   if muons[selectedMusIdx[imu]].p4().DeltaR(fatjets[ifatjet].p4())<0.8:
    #print "mu pT,eta, phi %s %s %s" % (muons[selectedMusIdx[imu]].p4().Pt(),muons[selectedMusIdx[imu]].p4().Eta(),muons[selectedMusIdx[imu]].p4().Phi())
    isNotLep = False
  if not isNotLep: continue
  selectedFatJetsIdx.append(ifatjet)
'''
