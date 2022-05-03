TFile *file = TFile::Open("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/selection/SPS_Skim.root");
TTree *Tree = (TTree*)file->Get("Events");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);

const int N = 6;  Double_t mass_OFbins[N] = {5,10,15,20,25,999999}; Double_t mass_bins[N] = {5,10,15,20,25,30};

const int M = 6;  Double_t DeltaEta_OFbins[M] = {0,0.5,1.0,1.5,2.0,99.0}; Double_t DeltaEta_bins[N] = {0,0.5,1.0,1.5,2.0,2.5};

