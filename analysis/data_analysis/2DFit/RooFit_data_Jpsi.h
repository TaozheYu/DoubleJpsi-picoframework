//TFile *file = TFile::Open("/eos/user/t/tayu/DoubleJpsi/Skim/2016/Run2016.root");
//TFile *file = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018A_UL_v4.root");
//TTree *Tree = (TTree*)file->Get("rootuple/oniaTree");

TFile *file = TFile::Open("input/RemakeNtuple_Run2018.root");
TTree *Tree = (TTree*)file->Get("Events");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
void GetHisto2D(char CUT[1000], TTree *Tree, TH2F* & histo, const char *plotX,const char *plotY, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX);
void MakeHistoErrors(int m, TH1F* &histo, TH1F* histo_up, TH1F* histo_do);
