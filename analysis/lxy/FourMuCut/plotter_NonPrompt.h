//TFile *file_DPS2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/DPSToJJ_Ntuple_lxy/DPSToJJ_Ntuple_MuMuLxy/DPS_2018_MuMuLxy_v2.root");
//TFile *file_BBbarToJpsi2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/BBarToJJ_lxy/BBarToJJ_2018_Ntuple_MuMulxy_v2.root");

TFile *file_DPS2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/DPSToJJ_Ntuple_lxy/DPSToJJ_NewNtuple/DPS18.root");
TFile *file_SPS2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/SPSToJJ_Ntuple_lxy/SPS18.root");
TFile *file_BBbarToJpsi2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/BBarToJJ_lxy/bDecay18.root");

TTree *Tree_DPS2018 = (TTree*)file_DPS2018->Get("rootuple/oniaTree");
TTree *Tree_SPS2018 = (TTree*)file_SPS2018->Get("rootuple/oniaTree");
TTree *Tree_BBbarToJpsi2018 = (TTree*)file_BBbarToJpsi2018->Get("rootuple/oniaTree");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
