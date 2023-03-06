#include "RemakeNtuple.h"
#include <TRandom2.h>

void RemakeNtuple(){
  char openFile[500];   sprintf(openFile, "/eos/user/t/tayu/DoubleJpsi/Run2016C_UL.root");
  //char openFile[500];   sprintf(openFile, "/eos/user/t/tayu/DoubleJpsi/Run2016D_UL.root");
  //char openFile[500];   sprintf(openFile, "/eos/user/t/tayu/DoubleJpsi/Run2016E_UL.root");
  //char openFile[500];   sprintf(openFile, "/eos/user/t/tayu/DoubleJpsi/Run2016FH_UL.root");
  char openTree[500];   sprintf(openTree, "rootuple/oniaTree");
  TFile *file = TFile::Open(openFile);
  TTree *Tree = (TTree*)file->Get(openTree);
  Init(Tree);
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  //Creat output root file
  char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_Run2016C_UL.root");
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_Run2016D_UL.root");
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_Run2016E_UL.root");
  //char NewFileName[500]; sprintf(NewFileName, "RemakeNtuple_Run2016FH_UL.root");
  //char NewFileName[500]; sprintf(NewFileName, "/publicfs/cms/user/yutz/Double_JPci/Run2018UL/RemakeNtuple_Run2018.root");
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
    if (fourMuFit_VtxProb->size()<0.5) continue;
    //This aim to assign the Jpsi1 and Jpsi2 randomly
    TRandom2 *gRandom = new TRandom2(0);
    float random = gRandom->Uniform(0,1);
    //cout<< random <<endl;
    if (random > 0.5){
    Jpsi1_mass = fourMuFit_ups1_mass->at(0);
    Jpsi2_mass = fourMuFit_ups2_mass->at(0);
    Jpsi1_pt = fourMuFit_ups1_pt->at(0);
    Jpsi2_pt = fourMuFit_ups2_pt->at(0);
    Jpsi1_rapidity = fourMuFit_ups1_rapidity->at(0);
    Jpsi2_rapidity = fourMuFit_ups2_rapidity->at(0);

    Jpsi1_Lxy = fourMuFit_ups1_Lxy_MC->at(0);
    Jpsi2_Lxy = fourMuFit_ups2_Lxy_MC->at(0);

    Jpsi1_LxyPV = fourMuFit_ups1_LxyPV_MC->at(0);
    Jpsi2_LxyPV = fourMuFit_ups2_LxyPV_MC->at(0);

    Jpsi1_cTau = fourMuFit_ups1_cTau_MC->at(0);
    Jpsi2_cTau = fourMuFit_ups2_cTau_MC->at(0);

    Jpsi1_LxyPVSig = fourMuFit_ups1_LxyPVSig_MC->at(0);
    Jpsi2_LxyPVSig = fourMuFit_ups2_LxyPVSig_MC->at(0);
    } else{
    Jpsi1_mass = fourMuFit_ups2_mass->at(0);
    Jpsi2_mass = fourMuFit_ups1_mass->at(0);
    Jpsi1_pt = fourMuFit_ups2_pt->at(0);
    Jpsi2_pt = fourMuFit_ups1_pt->at(0);
    Jpsi1_rapidity = fourMuFit_ups2_rapidity->at(0);
    Jpsi2_rapidity = fourMuFit_ups1_rapidity->at(0);

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
    fourMu_VtxProb = fourMuFit_VtxProb->at(0);
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
}

void initialVal(){
  
  Jpsi1_mass = -99;
  Jpsi2_mass = -99;
  Jpsi1_pt = -99;
  Jpsi2_pt = -99;
  Jpsi1_rapidity = -99;
  Jpsi2_rapidity = -99;

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
}
