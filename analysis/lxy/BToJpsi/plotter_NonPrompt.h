TFile *file_DPS2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/DPSToJJ_Ntuple_lxy/DPSToJJ_Ntuple_MuMuLxy/DPS_2018_MuMuLxy_v2.root");
TFile *file_BToJpsi2018 = TFile::Open("/eos/user/t/tayu/DoubleJpsi/BToJpsi_lxy/BsToJPsiPhi_MuMulxy_v2.root");

TTree *Tree_DPS2018 = (TTree*)file_DPS2018->Get("rootuple/oniaTree");
TTree *Tree_BToJpsi2018 = (TTree*)file_BToJpsi2018->Get("rootuple/oniaTree");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
