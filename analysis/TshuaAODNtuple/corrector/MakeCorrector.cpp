#define MakeCorrector_cxx
#include "MakeCorrector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MakeCorrector(){
  char openTree[500];   sprintf(openTree, "rootuple/oniaTree"); 
  TFile *file = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/TsHAODNtuple/DPSToJJ_13TeV_pythia8_official.root");
  TTree *Tree = (TTree*)file->Get(openTree);
  Init(Tree);
  if (fChain == 0) return;
  const int N = 11;  Double_t xbins[N] = {5.,6.,7.,8.,9.,10.,12.,15.,20.,25.,40.};
  const int M = 11;   Double_t ybins[M] = {-2.,-1.75,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,1.75,2.};
  TH2F *h_total_Jpsi = new TH2F("h_total_Jpsi","h_total_Jpsi",N-1 , xbins, M-1, ybins);
  TH2F *h_acc_Jpsi = new TH2F("h_acc_Jpsi","h_acc_Jpsi",N-1 , xbins, M-1, ybins);
  TH2F *h_reco_Jpsi = new TH2F("h_reco_Jpsi","h_reco_Jpsi",N-1 , xbins, M-1, ybins);
  TH2F *h_eff_Jpsi = new TH2F("h_eff_Jpsi","h_eff_Jpsi",N-1 , xbins, M-1, ybins);
  TH2F *h_vtx_Jpsi = new TH2F("h_vtx_Jpsi","h_vtx_Jpsi",N-1 , xbins, M-1, ybins);

  const int L = 7; Double_t pt_bins[L] = {5.,7.,9.,12.,15.,25.,40.};
  TH2F *h_JpsiPair = new TH2F("h_JpsiPair","h_JpsiPair",L-1 , pt_bins, L-1, pt_bins);
  TH2F *h_HLT_JpsiPair = new TH2F("h_HLT_JpsiPair","h_HLT_JpsiPair",L-1 , pt_bins, L-1, pt_bins);
  TH2F *h_evt_JpsiPair = new TH2F("h_evt_JpsiPair","h_evt_JpsiPair",L-1 , pt_bins, L-1, pt_bins);


  Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //nb = fChain->GetEntry(jentry);   nbytes += nb;
    fChain->GetEntry(jentry);
    for (int i=0; i<GENups_pt->size(); i++){
      if(!(fabs(GENups_y->at(i))<2)) continue;
      //if(!(GENups_DaughtersId->at(i)==13)) continue; 
      if(!((GENups_DaughtersId->at(2*i)==13)&&(GENups_DaughtersId->at(2*i+1)==13))) continue; 
      h_total_Jpsi->Fill(GENups_pt->at(i),GENups_y->at(i));
      if(!((fabs(GENups_Daughter_mueta->at(2*i))<2.4)&&(fabs(GENups_Daughter_mueta->at(2*i+1))<2.4))) continue; 
      //if(!((GENups_Daughter_mupt->at(2*i)>3.5)&&(GENups_Daughter_mupt->at(2*i+1)>3.5))) continue; 
      h_acc_Jpsi->Fill(GENups_pt->at(i),GENups_y->at(i));    
      }
    if (ups1_pt_GenMatched>0) h_reco_Jpsi->Fill(ups1_pt_GenMatched,ups1_y_GenMatched);
    if (ups2_pt_GenMatched>0) h_reco_Jpsi->Fill(ups2_pt_GenMatched,ups2_y_GenMatched);
    
    if (ups1_pt_GenMatched_ID>0) h_eff_Jpsi->Fill(ups1_pt_GenMatched_ID,ups1_y_GenMatched_ID);
    if (ups2_pt_GenMatched_ID>0) h_eff_Jpsi->Fill(ups2_pt_GenMatched_ID,ups2_y_GenMatched_ID);

    if (ups1_pt_GenMatched_ID_OS_VTX>0) h_vtx_Jpsi->Fill(ups1_pt_GenMatched_ID_OS_VTX,ups1_y_GenMatched_ID_OS_VTX);
    if (ups2_pt_GenMatched_ID_OS_VTX>0) h_vtx_Jpsi->Fill(ups2_pt_GenMatched_ID_OS_VTX,ups2_y_GenMatched_ID_OS_VTX);

    if (fourMuFit_ups1_pt->size()>0&&fourMuFit_ups2_pt->size()>0) h_JpsiPair->Fill(fourMuFit_ups1_pt->at(0),fourMuFit_ups2_pt->at(0)) ;
      
    if (!((trigger&16)==0)){
       if (fourMuFit_ups1_pt->size()>0&&fourMuFit_ups2_pt->size()>0) h_HLT_JpsiPair->Fill(fourMuFit_ups1_pt->at(0),fourMuFit_ups2_pt->at(0)) ;
       if (fourMuFit_ups1_pt->size()>0&&fourMuFit_ups2_pt->size()>0&&fourMuFit_VtxProb->at(0)>0.05) h_evt_JpsiPair->Fill(fourMuFit_ups1_pt->at(0),fourMuFit_ups2_pt->at(0));
      
    }
    
      // if (Cut(ientry) < 0) continue;
  }
  char NewFileName[500]; sprintf(NewFileName, "Gen_Event.root");
  TFile f(NewFileName,"new");
  h_total_Jpsi->Write();
  h_acc_Jpsi->Write();
  h_reco_Jpsi->Write();
  h_eff_Jpsi->Write();
  h_vtx_Jpsi->Write();
  h_JpsiPair->Write();
  h_HLT_JpsiPair->Write();
  h_evt_JpsiPair->Write();
  f.Close();
  cout<<"Job Finish"<<endl;
  
}
