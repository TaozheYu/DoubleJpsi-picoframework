TFile *file = TFile::Open("/eos/user/t/tayu/DoubleJpsi/BToJpsi_lxy/BsToJPsiPhi_MuMulxy_v2.root");
TTree *Tree = (TTree*)file->Get("rootuple/oniaTree");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
