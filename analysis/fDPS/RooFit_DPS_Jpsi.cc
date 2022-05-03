#include "RooFit_DPS_Jpsi.h"
#include "CalculateCS.h"
//#ifndef __CINT__
//#include "RooGlobalFunc.h"
//#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFFTConvPdf.h"
//#include "RooCBShape.h"
//#include "RooGaussian.h"
//#include "TCanvas.h"
//#include "RooPlot.h"
//#include "TTree.h"
//#include "TH1D.h"
//#include "TRandom.h"
using namespace RooFit;

void RooFit_DPS_Jpsi(const char *WhichJpsi, const char *CutCategory, float min_cut, float max_cut){
  /*float min_cut, float max_cut*/
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location
  gStyle->SetTitleY(0.96); //title Y location
  gStyle->SetPaintTextFormat(".2f");
  gErrorIgnoreLevel = kError;
  using namespace std;
  TCanvas* canvas = new TCanvas("canvas","canvas",0,0,900,600);
  TH1F *DPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  RooRealVar x("x", "x", 2.80, 3.40);
  x.setBins(1000,"cache");
  RooRealVar PDGmZ("PDGmZ","PDGmZ",3.092);
  RooRealVar PDGwZ("PDGwZ","PDGwZ",0.05,0.01,0.08);
  RooBreitWigner PDGBW("PDGBW","PDGBW",x,PDGmZ,PDGwZ);
  PDGmZ.setConstant(kTRUE);
  PDGwZ.setConstant(kTRUE);
  RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi mean", 3.092,2,4); 
  RooRealVar Jpsi_width("Jpsi_width", "Jpsi width",0.05,0,1);
  RooRealVar alpha_1("alpha_1", "",1,0,10);
  RooRealVar alpha_2("alpha_2", "",1,0,10);
  RooRealVar n_1("n_1","",10,0,100);
  RooRealVar n_2("n_2","",10,0,500);
  RooDoubleCB DCB("DCB","DCB",x,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);

  string  WhichJpsi_= WhichJpsi;  string CutCategory_=CutCategory;
  
  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  if(WhichJpsi_.find("Jpsi1")!=string::npos) {name.push_back("Jpsi1_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{12}[GeV]");}
  if(WhichJpsi_.find("Jpsi2")!=string::npos) {name.push_back("Jpsi2_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{34}[GeV]");}


  //for(int nbin=0; nbin<N-1; nbin++){
    //float min = mass_bins[nbin];
    //float max = mass_bins[nbin+1];
    float min = min_cut;
    float max = max_cut;
     
    for(int i=0; i<name.size(); i++){
      const char *plot = name[i];
      char CUTpre[1000];
      char CUT[1000];
      char MassCUT[1000]; sprintf(MassCUT,"1");
      char DeltaEtaCUT[1000]; sprintf(DeltaEtaCUT,"1"); 
      //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
      sprintf(CUTpre,"(((Jpsi1_pt>6 && abs(Jpsi1_eta)>1)||(Jpsi1_pt>7 && abs(Jpsi1_eta)<1)) && ((Jpsi2_pt>6 && abs(Jpsi2_eta)>1)||(Jpsi2_pt>7 && abs(Jpsi2_eta)<1)))");
      if(CutCategory_.find("JpsiPair_mass")!=string::npos) {sprintf(MassCUT,"(JpsiPair_mass>%f && JpsiPair_mass<%f)",min,max);}
      if(CutCategory_.find("JpsiPair_DeltaEta")!=string::npos) {sprintf(DeltaEtaCUT,"(JpsiPair_DeltaEta>%f && JpsiPair_DeltaEta<%f)",min,max);}
      //sprintf(CUT,"1* %s", CUTpre);
      //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12)) *%s*%s ", CUTpre, MassCUT);

      sprintf(CUT,"1. *%s*%s*%s", CUTpre, MassCUT, DeltaEtaCUT);
      GetHisto(CUT, Tree, DPS ,plot,bin[i],Min[i],Max[i]);

      cout<<DPS->Integral()<<endl;; 

      for(int j=1; j<DPS->GetNbinsX()+1; j++){
	      //data->SetBinContent(j,data_pre->GetBinContent(j));

	      DPS->SetBinContent(j,DPS->GetBinContent(j));
	      DPS->SetBinError(j,sqrt(DPS->GetBinContent(j)));
      }

      TH1F *h_DCB; h_DCB = (TH1F*)DPS->Clone(); h_DCB->SetLineColor(kBlue);

      RooDataHist roohist_DPS("roohist_DPS","roohist_DPS",x,DPS);
      RooPlot* frame_DPS = x.frame();
      DCB.fitTo(roohist_DPS,SumW2Error(kTRUE));
      roohist_DPS.plotOn(frame_DPS);
      DCB.plotOn(frame_DPS,LineColor(kBlue)); 
      frame_DPS->SetTitle(""); 
      frame_DPS->GetXaxis()->SetTitle(axis[i]); 
      frame_DPS->GetYaxis()->SetTitle("Events/(0.01 GeV)"); 
      frame_DPS->Draw();

      Double_t data_num = DPS->Integral();
      CalculateCS("DPS", 1., 0., data_num);

      TPad *pad = new TPad("pad","pad",0.01,0.01,0.99,0.99);
      gPad->RedrawAxis();
      TString channelText = "";
      float channelTextFont   = 42;
      float channelTextSize   = 0.060;
      TString cmsText     = "CMS";
      float cmsTextFont   = 61;  // default is helvetic-bold
      bool writeExtraText = true;
      TString extraText   = "Simulation";
      //TString extraText   = "";
      float extraTextFont = 52;  // default is helvetica-italics
      // text sizes and text offsets with respect to the top frame in unit of the top margin size
      float lumiTextSize     = 0.7;
      float lumiTextOffset   = 0.2;
      float cmsTextSize      = 0.75;
      float cmsTextOffset    = 0.1;  // only used in outOfFrame version
      float relPosX    = 0.045;
      float relPosY    = 0.035;
      float relExtraDY = 1.2;
      // ratio of "CMS" and extra text size
      float extraOverCmsTextSize  = 0.65;
      TString lumi_13TeV;
      lumi_13TeV = "36.3 fb^{-1}";
      TString lumiText;
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
      float t = pad->GetTopMargin();
      float b = pad->GetBottomMargin();
      float r = pad->GetRightMargin();
      float l = pad->GetLeftMargin();
      TLatex latex;
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(kBlack);
      float extraTextSize = extraOverCmsTextSize*cmsTextSize;
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(lumiTextSize*t);
      latex.DrawLatex(1-r,0.92,lumiText);
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11);
      latex.SetTextSize(cmsTextSize*t);
      latex.DrawLatex(l+0.08, 1-t+lumiTextOffset*t-0.12,cmsText);//0.075
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.DrawLatex(l+0.08, 1-t+lumiTextOffset*t-0.19, extraText);
      //latex.SetTextFont(channelTextFont);
      //latex.SetTextSize(channelTextSize);
      //latex.DrawLatex(l+0.165, 1-t+lumiTextOffset*t-0.19, channelText1);
      //latex.DrawLatex(l+0.165, 1-t+lumiTextOffset*t-0.28, channelText2);

      TLegend *pl2 = new TLegend(0.55,0.70,0.88,0.87);
      pl2->SetTextSize(0.030);
      pl2->SetTextFont(62);
      pl2->SetFillColor(0);
      TLegendEntry *ple2 = pl2->AddEntry(&roohist_DPS, "JpsiJpsi simulated events",  "PE");
      pl2->AddEntry(h_DCB, "fit", "L");
      pl2->SetBorderSize(0);
      pl2->Draw();

      TString NAME = name[i];
      char mass_min[10];  sprintf(mass_min,"%f",   min); TString MIN = mass_min;
      char mass_max[10];  sprintf(mass_max,"%f",   max); TString MAX = mass_max;
      cout<<"Finished "<<NAME<<endl;

      canvas->SaveAs("/afs/cern.ch/work/t/tayu/CMSSW_10_2_9_DoubleJpci/src/PicoFramework/analysis/fDPS/figures/DPS_fit_"+NAME+"_"+CutCategory+"_"+MIN+"_"+MAX+".pdf");
    }
  //}
}


void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX){
  char input[50]; sprintf(input,"%s>>h(%i,%f,%f)",plot,BIN,MIN,MAX);
  TH1F *temp = new TH1F("", "", BIN,MIN,MAX);
  if(Tree->Draw(input,CUT)){ Tree->Draw(input,CUT); TH1F* h=(TH1F*)gDirectory->Get("h"); histo = (TH1F*)h->Clone(); delete h; } else { histo = (TH1F*)temp->Clone(); }
  histo->SetDirectory(0);
  histo->SetName("histo");
  delete temp;
} 
