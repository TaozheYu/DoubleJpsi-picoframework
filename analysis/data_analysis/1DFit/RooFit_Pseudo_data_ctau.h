TFile *file_Run2018A_UL = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018A_UL_v4.root");

TTree *Tree_Run2018A_UL = (TTree*)file_Run2018A_UL->Get("rootuple/oniaTree");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
