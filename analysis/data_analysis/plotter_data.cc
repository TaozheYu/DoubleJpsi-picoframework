#include "plotter_data.h"
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

void plotter_data(){
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
  canvas->SetLogy();
  TH1F *Run2018A_UL; 

  float twomuon_cut = 0.05;
  float fourmuon_cut = 0.01;
  float Jpsi_mass_do = 2.90;
  float Jpsi_mass_up = 3.30;


  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  //name.push_back("fourMuFit_lxyPV_noMC"); bin.push_back(150); Min.push_back(-0.03); Max.push_back(0.03); axis.push_back("fourMuFit_lxyPV_noMC");
  //name.push_back("fourMuFit_ctauPV_noMC"); bin.push_back(100); Min.push_back(-0.05); Max.push_back(0.05); axis.push_back("fourMuFit_ctauPV_noMC");
  //name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(150); Min.push_back(-0.1); Max.push_back(0.15); axis.push_back("L_{xy}(J/#psi) [cm]");
  //name.push_back("fourMuFit_ups1_mass"); bin.push_back(60); Min.push_back(Jpsi_mass_do); Max.push_back(Jpsi_mass_up); axis.push_back("m_{12}[GeV]");
  //name.push_back("fourMuFit_ups2_mass"); bin.push_back(60); Min.push_back(Jpsi_mass_do); Max.push_back(Jpsi_mass_up); axis.push_back("m_{34}[GeV]");
  name.push_back("abs(fourMuFit_ups1_cTau_MC)"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.2); axis.push_back("J/#psi c#tau [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    char CUT_sideband[1000];
    //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
    //sprintf(CUT,"v_mumufit_VtxCL_noMC>%f && fourMuFit_VtxProb>%f && fourMuFit_ups1_mass>%f && fourMuFit_ups1_mass<%f && fourMuFit_ups2_mass>%f && fourMuFit_ups2_mass<%f",twomuon_cut, fourmuon_cut, Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
    sprintf(CUT,"fourMuFit_ups1_mass>%f && fourMuFit_ups1_mass<%f && fourMuFit_ups2_mass>%f && fourMuFit_ups2_mass<%f", Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
    GetHisto(CUT, Tree_Run2018A_UL, Run2018A_UL ,plot,bin[i],Min[i],Max[i]);
       
    cout<<Run2018A_UL->Integral()<<endl;; 
   
    for(int j=1; j<Run2018A_UL->GetNbinsX()+1; j++){

      Run2018A_UL->SetBinContent(j,Run2018A_UL->GetBinContent(j));
      Run2018A_UL->SetBinError(j,sqrt(Run2018A_UL->GetBinContent(j)));

    }

    //TH1F *h_DCB; h_DCB = (TH1F*)Run2018A_UL->Clone(); h_DCB->SetLineColor(kBlue);
    Run2018A_UL->SetLineColor(1);   

    Run2018A_UL->SetTitle(""); 
    Run2018A_UL->GetXaxis()->SetTitle(axis[i]); 
    Run2018A_UL->GetYaxis()->SetTitle("Events"); 
    Run2018A_UL->Draw("samehisto");
    //Run2018A_UL->SetMaximum(1);
    //hs->SetMinimum(3);

 

    TPad *pad = new TPad("pad","pad",0.01,0.01,0.99,0.99);
    gPad->RedrawAxis();
    TString channelText = "";
    float channelTextFont   = 42;
    float channelTextSize   = 0.060;
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    bool writeExtraText = true;
    TString extraText   = "preliminary";
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
    lumi_13TeV = "";
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



    TLegend *pl2 = new TLegend(0.60,0.70,0.88,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    //TLegendEntry *ple2 = pl2->AddEntry(&roohist_Run2018A_UL, "JpsiJpsi simulated events",  "PE");
    TLegendEntry *ple2 = pl2->AddEntry(Run2018A_UL, "data 2018",  "L");
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[i];
    cout<<"Finished "<<NAME<<endl;

    canvas->SaveAs("plots/"+NAME+".pdf");
  }
}


void GetHisto(char CUT[1000], TTree *Tree, TH1F* & histo, const char *plot, int BIN, float MIN, float MAX){
  char input[50]; sprintf(input,"%s>>h(%i,%f,%f)",plot,BIN,MIN,MAX);
  TH1F *temp = new TH1F("", "", BIN,MIN,MAX);
  if(Tree->Draw(input,CUT)){ Tree->Draw(input,CUT); TH1F* h=(TH1F*)gDirectory->Get("h"); histo = (TH1F*)h->Clone(); delete h; } else { histo = (TH1F*)temp->Clone(); }
  histo->SetDirectory(0);
  histo->SetName("histo");
  delete temp;
} 
