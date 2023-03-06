#include "RooFit_prompt_and_nonprompt.h"
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

void RooFit_prompt_and_nonprompt_ctau(){
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
  TH1F *BBbar; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  RooRealVar x("x", "x", -0.2 , 0.2);
  x.setBins(1000,"cache");
  //RooRealVar mean("mean", "mean",0.00065,-0.002,0.002); 
  //RooRealVar sigma("sigma", "sigma",0.0037,0,0.01);
  //RooRealVar rlife("rlife", "rlife",0.016,0,0.5);
  /* 
  //DPS fit function parameter
  RooRealVar DPS_mean("DPS_mean", "DPS mean", 0,-0.02,0.02); 
  RooRealVar DPS_width("DPS_width", "DPS width",0.01,0,0.08);
  RooRealVar DPS_alpha_1("DPS_alpha_1", "",1,0,10);
  RooRealVar DPS_alpha_2("DPS_alpha_2", "",1,0,10);
  RooRealVar DPS_n_1("DPS_n_1","",10,0,500);
  RooRealVar DPS_n_2("DPS_n_2","",10,0,500);
  RooDoubleCB DPS_DCB("DCB","DCB",x,DPS_mean,DPS_width,DPS_alpha_1,DPS_n_1,DPS_alpha_2,DPS_n_2);
  */
  //DPS double gauss fit
  RooRealVar DPS_mean1("DPS_mean1", "DPS mean1", 0.000055728, -0.1, 0.1); 
  RooRealVar DPS_width1("DPS_width1", "DPS width1",0.00174684);

  //RooRealVar DPS_mean2("DPS_mean2", "DPS mean2", 0,-0.02,0.04); 
  RooRealVar DPS_width2("DPS_width2", "DPS width2",0.00369649);
  
  RooGaussian DPS_gauss1("gauss1","gaussian1 PDF",x,DPS_mean1,DPS_width1);
  RooGaussian DPS_gauss2("gauss2","gaussian2 PDF",x,DPS_mean1,DPS_width2);
  RooRealVar gfrac_gauss1("gfrac_gauss1","",0.801);
  RooAddPdf DPS_model("DPS_model","DPS_model",RooArgList(DPS_gauss1,DPS_gauss2),RooArgList(gfrac_gauss1));


  //BBbarToJpsiJpsi fit function parameter
  RooRealVar BBbar_mean("BBbar_mean", "BBbar mean", 0.005,0,0.02); 
  RooRealVar BBbar_width("BBbar_width", "BBbar width",0.01,0,0.06);
  RooRealVar BBbar_alpha_1("BBbar_alpha_1", "",1,0,10);
  RooRealVar BBbar_alpha_2("BBbar_alpha_2", "",1,0,10);
  RooRealVar BBbar_n_1("BBbar_n_1","",10,0,500);
  RooRealVar BBbar_n_2("BBbar_n_2","",10,0,500);
  RooDoubleCB BBbar_DCB("BBbar_DCB","BBbar_DCB",x,BBbar_mean,BBbar_width,BBbar_alpha_1,BBbar_n_1,BBbar_alpha_2,BBbar_n_2);


  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,DCB);
  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,CB_1);
  //RooGExpModel CtauModel("ctauModel","ctauModel",x,mean,sigma,rlife);

  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  name.push_back("MuMuFit_ctauPV_MC"); bin.push_back(50); Min.push_back(-0.02); Max.push_back(0.04); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(100); Min.push_back(-0.1); Max.push_back(0.12); axis.push_back("L_{xy}(J/#psi) [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
    //sprintf(CUT,"");
    sprintf(CUT,"fourMuFit_VtxProb>0.01");
    GetHisto(CUT, Tree_DPS2018, DPS ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT, Tree_BBbarToJpsi2018, BBbar ,plot,bin[i],Min[i],Max[i]);
       
    cout<<DPS->Integral()<<endl;; 
   
    for(int j=1; j<DPS->GetNbinsX()+1; j++){
      //data->SetBinContent(j,data_pre->GetBinContent(j));

      DPS->SetBinContent(j,DPS->GetBinContent(j));
      DPS->SetBinError(j,sqrt(DPS->GetBinContent(j)));

      BBbar->SetBinContent(j,BBbar->GetBinContent(j));
      BBbar->SetBinError(j,sqrt(BBbar->GetBinContent(j)));
    }

    DPS->SetLineColor(kBlue); DPS->SetLineWidth(2); DPS->SetMarkerColor(kBlue); DPS->SetMarkerSize(1.3); DPS->SetMarkerStyle(20);
    BBbar->SetLineColor(kRed); BBbar->SetLineWidth(2); BBbar->SetMarkerColor(kRed); BBbar->SetMarkerSize(1.3); BBbar->SetMarkerStyle(20);

    RooDataHist roohist_DPS("roohist_DPS","roohist_DPS",x,DPS);
    RooDataHist roohist_BBbar("roohist_BBbar","roohist_BBbar",x,BBbar);

    RooPlot* frame_DPS = x.frame();

    DPS_model.fitTo(roohist_DPS,SumW2Error(kTRUE));
    roohist_DPS.plotOn(frame_DPS,LineColor(kBlue),MarkerColor(kBlue));
    DPS_model.plotOn(frame_DPS,LineColor(kBlue)); 

    BBbar_DCB.fitTo(roohist_BBbar,SumW2Error(kTRUE));
    roohist_BBbar.plotOn(frame_DPS,LineColor(kRed),MarkerColor(kRed));
    BBbar_DCB.plotOn(frame_DPS,LineColor(kRed)); 

    frame_DPS->SetTitle(""); 
    frame_DPS->GetXaxis()->SetTitle(axis[i]); 
    frame_DPS->GetYaxis()->SetTitle("Events"); 
    frame_DPS->SetMaximum(100000);;
    frame_DPS->SetMinimum(0.1);;
    frame_DPS->Draw();


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
    //lumi_13TeV = "36.3 fb^{-1}";
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




    TLegend *pl2 = new TLegend(0.53,0.70,0.80,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(DPS, "DPS simulation events",  "PE");
    pl2->AddEntry(BBbar, "b#bar{b} #rightarrow J/#psi J/#psi simulation events",  "PE");
    pl2->AddEntry(DPS, "DPS fit", "L");
    pl2->AddEntry(BBbar, "b#bar{b} #rightarrow J/#psi J/#psi fit", "L");
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[i];
    cout<<"Finished "<<NAME<<endl;

    canvas->SaveAs("Fit_plots/fit_"+NAME+".pdf");
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
