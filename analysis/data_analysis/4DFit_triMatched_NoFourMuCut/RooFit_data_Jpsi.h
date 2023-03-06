//TFile *file = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/Run2018UL/Run2018.root");
//TTree *Tree = (TTree*)file->Get("rootuple/oniaTree");

//TFile *file_DPS2018 = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/DPSToJJ_Ntuple_lxy/DPS18.root");
//TFile *file_SPS2018 = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/SPSToJJ_Ntuple_lxy/SPS18.root");
//TFile *file_BBbarToJpsi2018 = TFile::Open("/publicfs/cms/user/yutz/Double_JPci/BBarToJJ_Ntuple_lxy/bDecay18.root");

//TTree *Tree_DPS2018 = (TTree*)file_DPS2018->Get("rootuple/oniaTree");
//TTree *Tree_SPS2018 = (TTree*)file_SPS2018->Get("rootuple/oniaTree");
//TTree *Tree_BBbarToJpsi2018 = (TTree*)file_BBbarToJpsi2018->Get("rootuple/oniaTree");

TFile *file = TFile::Open("./input/RemakeNtuple_No4mucut_Run2018BC_triggerMatch.root");
TTree *Tree = (TTree*)file->Get("Events");

//TFile *file_DPS2018 = TFile::Open("./input/RemakeNtuple_DPS18.root");
//TFile *file_SPS2018 = TFile::Open("./input/RemakeNtuple_SPS18.root");
//TFile *file_BBbarToJpsi2018 = TFile::Open("./input/RemakeNtuple_BBbar18.root");

//TTree *Tree_DPS2018 = (TTree*)file_DPS2018->Get("Events");
//TTree *Tree_SPS2018 = (TTree*)file_SPS2018->Get("Events");
//TTree *Tree_BBbarToJpsi2018 = (TTree*)file_BBbarToJpsi2018->Get("Events");


void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
void GetHisto2D(char CUT[1000], TTree *Tree, TH2F* & histo, const char *plotX,const char *plotY, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX);
void GetHisto3D(char CUT[1000], TTree *Tree, TH3F* & histo, const char *plotX, const char *plotY, const char *plotZ, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX, int ZBIN, float ZMIN, float ZMAX);
void MakeHistoErrors(int m, TH1F* &histo, TH1F* histo_up, TH1F* histo_do);
