#include "plotter_NonPrompt.h"
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

void plotter_NonPrompt(){
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
  TH1F *DPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *BBToJpsi; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);


  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  //name.push_back("fourMuFit_lxyPV_noMC"); bin.push_back(150); Min.push_back(-0.03); Max.push_back(0.03); axis.push_back("fourMuFit_lxyPV_noMC");
  //name.push_back("fourMuFit_ctauPV_noMC"); bin.push_back(100); Min.push_back(-0.05); Max.push_back(0.05); axis.push_back("fourMuFit_ctauPV_noMC");
  //name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(300); Min.push_back(-0.2); Max.push_back(0.5); axis.push_back("L_{xy}(J/#psi) [cm]");
  name.push_back("MuMuFit_ctauPV_MC"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.1); axis.push_back("J/#psi c#tau [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
    //sprintf(CUT,"(fourMuFit_VtxProb>0.01)");
    sprintf(CUT,"");
    //sprintf(CUT,"fourMuFit_VtxProb>0.01");
    GetHisto(CUT, Tree_DPS2018, DPS ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT, Tree_BBToJpsi2018, BBToJpsi ,plot,bin[i],Min[i],Max[i]);
       
    cout<<DPS->Integral()<<endl;; 
   
    for(int j=1; j<DPS->GetNbinsX()+1; j++){
      DPS->SetBinContent(j,DPS->GetBinContent(j));
      DPS->SetBinError(j,sqrt(DPS->GetBinContent(j)));

      BBToJpsi->SetBinContent(j,BBToJpsi->GetBinContent(j));
      BBToJpsi->SetBinError(j,sqrt(BBToJpsi->GetBinContent(j)));

    }

    //TH1F *h_DCB; h_DCB = (TH1F*)DPS->Clone(); h_DCB->SetLineColor(kBlue);
    DPS->SetLineColor(kBlue);   BBToJpsi->SetLineColor(kRed);

    DPS->SetTitle(""); 
    DPS->GetXaxis()->SetTitle(axis[i]); 
    DPS->GetYaxis()->SetTitle("Events"); 
    DPS->Draw("samehisto");
    DPS->SetMaximum(100000);
    //hs->SetMinimum(3);

    BBToJpsi->Draw("samehisto");
 

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



    TLegend *pl2 = new TLegend(0.55,0.70,0.88,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    //TLegendEntry *ple2 = pl2->AddEntry(&roohist_DPS, "JpsiJpsi simulated events",  "PE");
    TLegendEntry *ple2 = pl2->AddEntry(DPS, "DPS",  "L");
    pl2->AddEntry(BBToJpsi, "BBar to Jpsi Jpsi", "L");
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[i];
    cout<<"Finished "<<NAME<<endl;

    canvas->SaveAs(NAME+".pdf");
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
