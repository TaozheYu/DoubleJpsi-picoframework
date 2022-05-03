//TFile *file = TFile::Open("/eos/user/t/tayu/DoubleJpsi/Skim/2016/Run2016.root");
TFile *file = TFile::Open("/eos/user/t/tayu/DoubleJpsi/Skim/2016/Run2016.root");
TTree *Tree = (TTree*)file->Get("Events");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
void GetHisto2D(char CUT[1000], TTree *Tree, TH2F* & histo, const char *plotX,const char *plotY, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX);
void MakeHistoErrors(int m, TH1F* &histo, TH1F* histo_up, TH1F* histo_do);


const int N = 6;  Double_t mass_OFbins[N] = {5,10,15,20,25,999999}; Double_t mass_bins[N] = {5,10,15,20,25,30};

const int M = 6;  Double_t DeltaEta_OFbins[M] = {0,0.5,1.0,1.5,2.0,99.0}; Double_t DeltaEta_bins[N] = {0,0.5,1.0,1.5,2.0,2.5};

