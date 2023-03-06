TFile *file = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018A_UL_new.root");
TTree *Tree = (TTree*)file->Get("rootuple/oniaTree");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
void GetHisto2D(char CUT[1000], TTree *Tree, TH2F* & histo, const char *plotX,const char *plotY, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX);
void GetHisto3D(char CUT[1000], TTree *Tree, TH3F* & histo, const char *plotX, const char *plotY, const char *plotZ, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX, int ZBIN, float ZMIN, float ZMAX);
void MakeHistoErrors(int m, TH1F* &histo, TH1F* histo_up, TH1F* histo_do);
