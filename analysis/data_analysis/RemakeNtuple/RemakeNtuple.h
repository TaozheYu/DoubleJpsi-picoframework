#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"

void SetNewTree(TTree* NewTree);
void initialVal();

TTree  *fChain; //pointer to analyzed TTree
Int_t  fCurrent; //!current Tree number

// Declaration of leaf types
vector<float> *fourMuFit_ups1_mass;
vector<float> *fourMuFit_ups2_mass;
vector<float> *fourMuFit_ups1_pt;
vector<float> *fourMuFit_ups2_pt;
vector<float> *fourMuFit_ups1_rapidity;
vector<float> *fourMuFit_ups2_rapidity;

vector<float> *fourMuFit_ups1_Lxy_MC;
vector<float> *fourMuFit_ups2_Lxy_MC;
vector<float> *fourMuFit_ups1_LxyPV_MC;
vector<float> *fourMuFit_ups2_LxyPV_MC;
vector<float> *fourMuFit_ups1_cTau_MC;
vector<float> *fourMuFit_ups2_cTau_MC;
vector<float> *fourMuFit_ups1_LxyPVSig_MC;
vector<float> *fourMuFit_ups2_LxyPVSig_MC;
vector<float> *fourMuFit_DistanceSig_MC;
vector<float> *fourMuFit_VtxProb;

// List of branches
TBranch *b_fourMuFit_ups1_mass;
TBranch *b_fourMuFit_ups2_mass;
TBranch *b_fourMuFit_ups1_pt;
TBranch *b_fourMuFit_ups2_pt;
TBranch *b_fourMuFit_ups1_rapidity;
TBranch *b_fourMuFit_ups2_rapidity;

TBranch *b_fourMuFit_ups1_Lxy_MC;
TBranch *b_fourMuFit_ups2_Lxy_MC;
TBranch *b_fourMuFit_ups1_LxyPV_MC;
TBranch *b_fourMuFit_ups2_LxyPV_MC;
TBranch *b_fourMuFit_ups1_cTau_MC;
TBranch *b_fourMuFit_ups2_cTau_MC;
TBranch *b_fourMuFit_ups1_LxyPVSig_MC;
TBranch *b_fourMuFit_ups2_LxyPVSig_MC;
TBranch *b_fourMuFit_DistanceSig_MC;
TBranch *b_fourMuFit_VtxProb;

float Jpsi1_mass;
float Jpsi2_mass;
float Jpsi1_pt;
float Jpsi2_pt;
float Jpsi1_rapidity;
float Jpsi2_rapidity;

float Jpsi1_Lxy;
float Jpsi2_Lxy;
float Jpsi1_LxyPV;
float Jpsi2_LxyPV;
float Jpsi1_cTau;
float Jpsi2_cTau;
float Jpsi1_LxyPVSig;
float Jpsi2_LxyPVSig;
float MaxLxyPVSig;
float Jpsi_DistanceSig;
float fourMu_VtxProb;

void Init(TTree *tree)
{
  fourMuFit_ups1_mass = 0;
  fourMuFit_ups2_mass = 0;
  fourMuFit_ups1_pt = 0;
  fourMuFit_ups2_pt = 0;
  fourMuFit_ups1_rapidity = 0;
  fourMuFit_ups2_rapidity = 0;

  fourMuFit_ups1_Lxy_MC = 0;
  fourMuFit_ups2_Lxy_MC = 0;
  fourMuFit_ups1_LxyPV_MC = 0;
  fourMuFit_ups2_LxyPV_MC = 0;
  fourMuFit_ups1_cTau_MC = 0;
  fourMuFit_ups2_cTau_MC = 0;
  fourMuFit_ups1_LxyPVSig_MC = 0;
  fourMuFit_ups2_LxyPVSig_MC = 0;
  fourMuFit_DistanceSig_MC = 0;
  fourMuFit_VtxProb = 0;
  
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1); 
  
  fChain->SetBranchAddress("fourMuFit_ups1_mass", &fourMuFit_ups1_mass, &b_fourMuFit_ups1_mass);
  fChain->SetBranchAddress("fourMuFit_ups2_mass", &fourMuFit_ups2_mass, &b_fourMuFit_ups2_mass);
  fChain->SetBranchAddress("fourMuFit_ups1_pt", &fourMuFit_ups1_pt, &b_fourMuFit_ups1_pt);
  fChain->SetBranchAddress("fourMuFit_ups2_pt", &fourMuFit_ups2_pt, &b_fourMuFit_ups2_pt);
  fChain->SetBranchAddress("fourMuFit_ups1_rapidity", &fourMuFit_ups1_rapidity, &b_fourMuFit_ups1_rapidity);
  fChain->SetBranchAddress("fourMuFit_ups2_rapidity", &fourMuFit_ups2_rapidity, &b_fourMuFit_ups2_rapidity);

  fChain->SetBranchAddress("fourMuFit_ups1_Lxy_MC", &fourMuFit_ups1_Lxy_MC, &b_fourMuFit_ups1_Lxy_MC);
  fChain->SetBranchAddress("fourMuFit_ups2_Lxy_MC", &fourMuFit_ups2_Lxy_MC, &b_fourMuFit_ups2_Lxy_MC);
  fChain->SetBranchAddress("fourMuFit_ups1_LxyPV_MC", &fourMuFit_ups1_LxyPV_MC, &b_fourMuFit_ups1_LxyPV_MC);
  fChain->SetBranchAddress("fourMuFit_ups2_LxyPV_MC", &fourMuFit_ups2_LxyPV_MC, &b_fourMuFit_ups2_LxyPV_MC);
  fChain->SetBranchAddress("fourMuFit_ups1_cTau_MC", &fourMuFit_ups1_cTau_MC, &b_fourMuFit_ups1_cTau_MC);
  fChain->SetBranchAddress("fourMuFit_ups2_cTau_MC", &fourMuFit_ups2_cTau_MC, &b_fourMuFit_ups2_cTau_MC);
  fChain->SetBranchAddress("fourMuFit_ups1_LxyPVSig_MC", &fourMuFit_ups1_LxyPVSig_MC, &b_fourMuFit_ups1_LxyPVSig_MC);
  fChain->SetBranchAddress("fourMuFit_ups2_LxyPVSig_MC", &fourMuFit_ups2_LxyPVSig_MC, &b_fourMuFit_ups2_LxyPVSig_MC);
  fChain->SetBranchAddress("fourMuFit_DistanceSig_MC", &fourMuFit_DistanceSig_MC, &b_fourMuFit_DistanceSig_MC);
  fChain->SetBranchAddress("fourMuFit_VtxProb", &fourMuFit_VtxProb, &b_fourMuFit_VtxProb);
}



