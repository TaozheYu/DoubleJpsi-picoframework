//TFile *file_DPS2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/DPSToJJ_Ntuple_lxy/DPSToJJ_Ntuple_MuMuLxy/DPS_2018_MuMuLxy_v2.root");
//TFile *file_BBbarToJpsi2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/BBarToJJ_lxy/BBarToJJ_2018_Ntuple_MuMulxy_v2.root");

TFile *file_DPS2018 = TFile::Open("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/analysis/RemakeNtuple/RemakeNtuple_DPS18.root");
TFile *file_SPS2018 = TFile::Open("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/analysis/RemakeNtuple/RemakeNtuple_SPS18.root");
TFile *file_BBbarToJpsi2018 = TFile::Open("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/analysis/RemakeNtuple/RemakeNtuple_bDecay18.root");

TTree *Tree_DPS2018 = (TTree*)file_DPS2018->Get("Events");
TTree *Tree_SPS2018 = (TTree*)file_SPS2018->Get("Events");
TTree *Tree_BBbarToJpsi2018 = (TTree*)file_BBbarToJpsi2018->Get("Events");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
