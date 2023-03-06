#include "RemakeNtuple.h"
#include <TRandom2.h>

void RemakeNtuple(){
  bool data = true;
  //char openFile[500];   sprintf(openFile, "/publicfs/cms/user/yutz/Double_JPci/BBarToJJ_Ntuple_lxy/bDecay18.root");
  //char openFile[500];   sprintf(openFile, "/publicfs/cms/user/yutz/Double_JPci/SPSToJJ_Ntuple_lxy/SPS18.root");
  //char openFile[500];   sprintf(openFile, "/publicfs/cms/user/yutz/Double_JPci/DPSToJJ_Ntuple_lxy/DPS18.root");
  //char openFile[500];   sprintf(openFile, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018.root");
  char openFile[500];   sprintf(openFile, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018B_UL_triggerMatch.root");
  //char openFile[500];   sprintf(openFile, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018C_UL_triggerMatch.root");
  char openTree[500];   sprintf(openTree, "rootuple/oniaTree");
  TFile *file = TFile::Open(openFile);
  TTree *Tree = (TTree*)file->Get(openTree);
  Init(Tree, data);
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  //Creat output root file
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_BBbar18.root");
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_No4mucut_BBbar18.root");
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_No4mucut_SPS18.root");
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_No4mucut_DPS18.root");
  //char NewFileName[500]; sprintf(NewFileName, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/RemakeNtuple_Run2018.root");
  //char NewFileName[500]; sprintf(NewFileName, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/RemakeNtuple_Run2018B_triggerMatch.root");
  char NewFileName[500]; sprintf(NewFileName, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/RemakeNtuple_No4mucut_Run2018B_triggerMatch.root");
  //char NewFileName[500]; sprintf(NewFileName, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/RemakeNtuple_No4mucut_Run2018C_triggerMatch.root");
  TFile f(NewFileName,"new");
  TTree *NewTree = new TTree("Events","Events");
  SetNewTree(NewTree);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    fChain->GetEntry(jentry);
    initialVal();
    if (fourMuFit_ups1_mass->size()<0.5) continue;
    if (fourMuFit_ups2_mass->size()<0.5) continue;
    if (fourMuFit_ups1_LxyPVSig_MC->size()<0.5) continue;
    if (fourMuFit_ups2_LxyPVSig_MC->size()<0.5) continue;
    //if (fourMuFit_VtxProb->size()<0.5) continue;
    //This aim to assign the Jpsi1 and Jpsi2 randomly
    TRandom2 *gRandom = new TRandom2(0);
    float random = gRandom->Uniform(0,1);
    if (random > 0.5){
    Jpsi1_mass = fourMuFit_ups1_mass->at(0);
    Jpsi2_mass = fourMuFit_ups2_mass->at(0);
    Jpsi1_pt = fourMuFit_ups1_pt->at(0);
    Jpsi2_pt = fourMuFit_ups2_pt->at(0);
    Jpsi1_rapidity = fourMuFit_ups1_rapidity->at(0);
    Jpsi2_rapidity = fourMuFit_ups2_rapidity->at(0);
    if (fourMuFit_ups1_VtxProb->size()>0.5) Jpsi1_vtxProb  = fourMuFit_ups1_VtxProb->at(0);
    if (fourMuFit_ups2_VtxProb->size()>0.5) Jpsi2_vtxProb  = fourMuFit_ups2_VtxProb->at(0);

    Jpsi1_Lxy = fourMuFit_ups1_Lxy_MC->at(0);
    Jpsi2_Lxy = fourMuFit_ups2_Lxy_MC->at(0);

    Jpsi1_LxyPV = fourMuFit_ups1_LxyPV_MC->at(0);
    Jpsi2_LxyPV = fourMuFit_ups2_LxyPV_MC->at(0);

    Jpsi1_cTau = fourMuFit_ups1_cTau_MC->at(0);
    Jpsi2_cTau = fourMuFit_ups2_cTau_MC->at(0);

    Jpsi1_LxyPVSig = fourMuFit_ups1_LxyPVSig_MC->at(0);
    Jpsi2_LxyPVSig = fourMuFit_ups2_LxyPVSig_MC->at(0);
    } else {
    Jpsi1_mass = fourMuFit_ups2_mass->at(0);
    Jpsi2_mass = fourMuFit_ups1_mass->at(0);
    Jpsi1_pt = fourMuFit_ups2_pt->at(0);
    Jpsi2_pt = fourMuFit_ups1_pt->at(0);
    Jpsi1_rapidity = fourMuFit_ups2_rapidity->at(0);
    Jpsi2_rapidity = fourMuFit_ups1_rapidity->at(0);
    if ( fourMuFit_ups2_VtxProb->size() > 0.5 ) Jpsi1_vtxProb  = fourMuFit_ups2_VtxProb->at(0);
    if ( fourMuFit_ups1_VtxProb->size() > 0.5 ) Jpsi2_vtxProb  = fourMuFit_ups1_VtxProb->at(0);

    Jpsi1_Lxy = fourMuFit_ups2_Lxy_MC->at(0);
    Jpsi2_Lxy = fourMuFit_ups1_Lxy_MC->at(0);

    Jpsi1_LxyPV = fourMuFit_ups2_LxyPV_MC->at(0);
    Jpsi2_LxyPV = fourMuFit_ups1_LxyPV_MC->at(0);

    Jpsi1_cTau = fourMuFit_ups2_cTau_MC->at(0);
    Jpsi2_cTau = fourMuFit_ups1_cTau_MC->at(0);

    Jpsi1_LxyPVSig = fourMuFit_ups2_LxyPVSig_MC->at(0);
    Jpsi2_LxyPVSig = fourMuFit_ups1_LxyPVSig_MC->at(0);    

    }
    Jpsi_DistanceSig = fourMuFit_DistanceSig_MC->at(0);
    
    MaxLxyPVSig = (fabs(fourMuFit_ups1_LxyPVSig_MC->at(0)+fourMuFit_ups2_LxyPVSig_MC->at(0))+fabs(fourMuFit_ups1_LxyPVSig_MC->at(0)-fourMuFit_ups2_LxyPVSig_MC->at(0)))/2;
    if (fourMuFit_VtxProb->size() > 0.5) {
       fourMu_VtxProb = fourMuFit_VtxProb->at(0);
    }
    trigger_matched = Trigger_Matched;
    NewTree->Fill();

    }
    NewTree->Write();
    f.Close();
    cout<<"Job Finish"<<endl;

}

void SetNewTree(TTree* NewTree){
  NewTree->Branch("Jpsi1_mass",         &Jpsi1_mass,      "Jpsi1_mass/F"      );
  NewTree->Branch("Jpsi2_mass",         &Jpsi2_mass,      "Jpsi2_mass/F"      );
  NewTree->Branch("Jpsi1_pt",         &Jpsi1_pt,      "Jpsi1_pt/F"      );
  NewTree->Branch("Jpsi2_pt",         &Jpsi2_pt,      "Jpsi2_pt/F"      );
  NewTree->Branch("Jpsi1_rapidity",         &Jpsi1_rapidity,      "Jpsi1_rapidity/F"      );
  NewTree->Branch("Jpsi2_rapidity",         &Jpsi2_rapidity,      "Jpsi2_rapidity/F"      );
  NewTree->Branch("Jpsi1_vtxProb",         &Jpsi1_vtxProb,      "Jpsi1_vtxProb/F"      );
  NewTree->Branch("Jpsi2_vtxProb",         &Jpsi2_vtxProb,      "Jpsi2_vtxProb/F"      );

  NewTree->Branch("Jpsi1_Lxy",         &Jpsi1_Lxy,      "Jpsi1_Lxy/F"      );
  NewTree->Branch("Jpsi2_Lxy",         &Jpsi2_Lxy,      "Jpsi2_Lxy/F"      );
  NewTree->Branch("Jpsi1_LxyPV",         &Jpsi1_LxyPV,      "Jpsi1_LxyPV/F"      );
  NewTree->Branch("Jpsi2_LxyPV",         &Jpsi2_LxyPV,      "Jpsi2_LxyPV/F"      );
  NewTree->Branch("Jpsi1_cTau",         &Jpsi1_cTau,      "Jpsi1_cTau/F"      );
  NewTree->Branch("Jpsi2_cTau",         &Jpsi2_cTau,      "Jpsi2_cTau/F"      );
  NewTree->Branch("Jpsi1_LxyPVSig",         &Jpsi1_LxyPVSig,      "Jpsi1_LxyPVSig/F"      );
  NewTree->Branch("Jpsi2_LxyPVSig",         &Jpsi2_LxyPVSig,      "Jpsi2_LxyPVSig/F"      );

  NewTree->Branch("MaxLxyPVSig",         &MaxLxyPVSig,      "MaxLxyPVSig/F"      );
  NewTree->Branch("Jpsi_DistanceSig",         &Jpsi_DistanceSig,      "Jpsi_DistanceSig/F"      );
  NewTree->Branch("fourMu_VtxProb",         &fourMu_VtxProb,      "fourMu_VtxProb/F"      );
  NewTree->Branch("trigger_matched",         &trigger_matched,      "trigger_matched/O"      );
}

void initialVal(){
  
  Jpsi1_mass = -99;
  Jpsi2_mass = -99;
  Jpsi1_pt = -99;
  Jpsi2_pt = -99;
  Jpsi1_rapidity = -99;
  Jpsi2_rapidity = -99;
  Jpsi1_vtxProb = -99;
  Jpsi2_vtxProb = -99;

  Jpsi1_Lxy = -99;
  Jpsi2_Lxy = -99;
  Jpsi1_LxyPV = -99;
  Jpsi2_LxyPV = -99;
  Jpsi1_cTau = -99;
  Jpsi2_cTau = -99;
  Jpsi1_LxyPVSig = -99;
  Jpsi2_LxyPVSig = -99;

  MaxLxyPVSig = -99;
  Jpsi_DistanceSig = -99;
  fourMu_VtxProb = -99;
  trigger_matched = false;
}
