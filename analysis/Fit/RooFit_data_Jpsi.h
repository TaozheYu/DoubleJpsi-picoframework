TFile *file = TFile::Open("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/selection/DoubleJpsi/test/2016/Run2016.root");
TTree *Tree = (TTree*)file->Get("Events");

void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX);
