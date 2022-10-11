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
  TH1F *BToJpsi; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *BBbar; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  RooRealVar x("x", "x", -0.3 , 1);
  x.setBins(1000,"cache");
  //RooRealVar mean("mean", "mean",0.00065,-0.002,0.002); 
  //RooRealVar sigma("sigma", "sigma",0.0037,0,0.01);
  //RooRealVar rlife("rlife", "rlife",0.016,0,0.5);
  
  //DPS fit function parameter
  RooRealVar DPS_mean("DPS_mean", "DPS mean", 0); 
  RooRealVar DPS_width("DPS_width", "DPS width",0.01,0,0.08);
  RooRealVar DPS_alpha_1("DPS_alpha_1", "",1,0,10);
  RooRealVar DPS_alpha_2("DPS_alpha_2", "",1,0,10);
  RooRealVar DPS_n_1("DPS_n_1","",10,0,100);
  RooRealVar DPS_n_2("DPS_n_2","",10,0,500);
  RooDoubleCB DPS_DCB("DCB","DCB",x,DPS_mean,DPS_width,DPS_alpha_1,DPS_n_1,DPS_alpha_2,DPS_n_2);


  //BBbarToJpsiJpsi fit function parameter
   
  RooRealVar BBbar_mean("BBbar_mean", "BBbar mean", 0.0012, 0, 0.03); 
  //RooRealVar BBbar_mean("BBbar_mean", "BBbar mean", 0.0023); 
  RooRealVar BBbar_width("BBbar_width", "BBbar width",0.01,0, 0.02);
  RooRealVar BBbar_alpha_1("BBbar_alpha_1", "",0.1,0,1);
  RooRealVar BBbar_alpha_2("BBbar_alpha_2", "",0.1,0,1);
  RooRealVar BBbar_n_1("BBbar_n_1","",10,0,500);
  RooRealVar BBbar_n_2("BBbar_n_2","",10,0,500);
  RooDoubleCB BBbar_DCB("BBbar_DCB","BBbar_DCB",x,BBbar_mean,BBbar_width,BBbar_alpha_1,BBbar_n_1,BBbar_alpha_2,BBbar_n_2);


  //BToJpsi fit function parameter

  RooRealVar BToJpsi_mean("BToJpsi_mean", "BToJpsi mean", 0.0012,0,0.03); 
  //RooRealVar BToJpsi_mean("BToJpsi_mean", "BToJpsi mean", 0.0023); 
  RooRealVar BToJpsi_width("BToJpsi_width", "BToJpsi width",0.001,0,0.01);
  RooRealVar BToJpsi_alpha_1("BToJpsi_alpha_1", "",1,0,10);
  RooRealVar BToJpsi_alpha_2("BToJpsi_alpha_2", "",1,0,10);
  RooRealVar BToJpsi_n_1("BToJpsi_n_1","",10,0,500);
  RooRealVar BToJpsi_n_2("BToJpsi_n_2","",10,0,500);
  RooDoubleCB BToJpsi_DCB("BToJpsi_DCB","BToJpsi_DCB",x,BToJpsi_mean,BToJpsi_width,BToJpsi_alpha_1,BToJpsi_n_1,BToJpsi_alpha_2,BToJpsi_n_2);

  //Non-prompt function for time
  //Exponential
  RooRealVar x_NonPrompt_DPS_Tau("x_NonPrompt_DPS_Tau", "x_NonPrompt_DPS_Tau", -3, -10, -1);      //1
  RooRealVar x_NonPrompt_BBbar_Tau("x_NonPrompt_BBbar_Tau", "x_NonPrompt_BBbar_Tau", -25.6, -30, -20);      //1
  RooRealVar x_NonPrompt_BToJpsi_Tau("x_NonPrompt_BToJpsi_Tau", "x_NonPrompt_BToJpsi_Tau", -27.9, -35, -25);      //1
  RooRealVar x_NonPrompt_Tau("x_NonPrompt_Tau", "x_NonPrompt_Tau", -27.8, -35, -25);      //1

  //Step Function
  RooRealVar x_StepValue("x_StepValue", "x_StepValue", 1);
  RooRealVar x_LowerLimit("x_LowerLimit", "x_LowerLimit", 0);
  RooRealVar x_HigherLimit("x_HigherLimit", "x_HigherLimit", 1);

  //Gaus
  RooRealVar x_NonPrompt_Mean("x_NonPrompt_Mean", "x_NonPrompt_Mean",0.0001,0.,0.1);
  RooRealVar x_NonPrompt_Sigma("x_NonPrompt_Sigma", "x_NonPrompt_Sigma", 0.0001, 0., 0.03);


  //Convolution
  //StepFunction
  RooAbsReal* Stepx = new RooStepFunction("Stepx", "Stepx", x, RooArgList(x_StepValue), RooArgList(x_LowerLimit, x_HigherLimit));
  //Generic Exponential 
  RooAbsPdf* Expx_NonPrompt = new RooGenericPdf("Expx_NonPrompt", "Expx_NonPrompt", "Stepx * exp(x * x_NonPrompt_Tau)", RooArgSet(*Stepx, x, x_NonPrompt_Tau));
  RooAbsPdf* Expx_NonPrompt_DPS = new RooGenericPdf("Expx_NonPrompt_DPS", "Expx_NonPrompt_DPS", "Stepx * exp(x * x_NonPrompt_DPS_Tau)", RooArgSet(*Stepx, x, x_NonPrompt_DPS_Tau));
  RooAbsPdf* Expx_NonPrompt_BBbar = new RooGenericPdf("Expx_NonPrompt_BBbar", "Expx_NonPrompt_BBbar", "Stepx * exp(x * x_NonPrompt_BBbar_Tau)", RooArgSet(*Stepx, x, x_NonPrompt_BBbar_Tau));
  RooAbsPdf* Expx_NonPrompt_BToJpsi = new RooGenericPdf("Expx_NonPrompt_BToJpsi", "Expx_NonPrompt_BToJpsi", "Stepx * exp(x * x_NonPrompt_BToJpsi_Tau)", RooArgSet(*Stepx, x, x_NonPrompt_BToJpsi_Tau));

  //Gaus for resolution
  RooAbsPdf* Gausx_NonPrompt = new RooGaussian("Gausx_NonPrompt", "Gausx_NonPrompt", x, x_NonPrompt_Mean, x_NonPrompt_Sigma);

  //Covolution
  RooFFTConvPdf NonPromptx("NonPromptx", "NonPromptx", x, *Expx_NonPrompt, *Gausx_NonPrompt);
  RooFFTConvPdf NonPromptx_DPS("NonPromptx_DPS", "NonPromptx_DPS", x, *Expx_NonPrompt_DPS, DPS_DCB);
  RooFFTConvPdf NonPromptx_BBbar("NonPromptx_BBbar", "NonPromptx_BBbar", x, *Expx_NonPrompt_BBbar, BBbar_DCB);
  RooFFTConvPdf NonPromptx_BToJpsi("NonPromptx_BToJpsi", "NonPromptx_BToJpsi", x, *Expx_NonPrompt_BToJpsi, BToJpsi_DCB);

  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  name.push_back("MuMuFit_ctauPV_MC"); bin.push_back(200); Min.push_back(-0.2); Max.push_back(0.6); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(100); Min.push_back(-0.2); Max.push_back(0.5); axis.push_back("L_{xy}(J/#psi) [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
    //sprintf(CUT,"");
    //sprintf(CUT,"fourMuFit_VtxProb>0.01");
    sprintf(CUT,"");
    GetHisto(CUT, Tree_DPS2018, DPS ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT, Tree_BToJpsi2018, BToJpsi ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT, Tree_BBbarToJpsi2018, BBbar ,plot,bin[i],Min[i],Max[i]);
       
    cout<<DPS->Integral()<<endl;; 
   
    for(int j=1; j<DPS->GetNbinsX()+1; j++){
      //data->SetBinContent(j,data_pre->GetBinContent(j));

      DPS->SetBinContent(j,DPS->GetBinContent(j));
      DPS->SetBinError(j,sqrt(DPS->GetBinContent(j)));

      BToJpsi->SetBinContent(j,BToJpsi->GetBinContent(j));
      BToJpsi->SetBinError(j,sqrt(BToJpsi->GetBinContent(j)));

      BBbar->SetBinContent(j,BBbar->GetBinContent(j));
      BBbar->SetBinError(j,sqrt(BBbar->GetBinContent(j)));
    }


    DPS->SetLineColor(kBlue); DPS->SetLineWidth(2); DPS->SetMarkerColor(kBlue); DPS->SetMarkerSize(1.3); DPS->SetMarkerStyle(20);
    BToJpsi->SetLineColor(kGreen); BToJpsi->SetLineWidth(2); BToJpsi->SetMarkerColor(kGreen); BToJpsi->SetMarkerSize(1.3); BToJpsi->SetMarkerStyle(20);
    BBbar->SetLineColor(kRed); BBbar->SetLineWidth(2); BBbar->SetMarkerColor(kRed); BBbar->SetMarkerSize(1.3); BBbar->SetMarkerStyle(20);

    RooDataHist roohist_DPS("roohist_DPS","roohist_DPS",x,DPS);
    RooDataHist roohist_BToJpsi("roohist_BToJpsi","roohist_BToJpsi",x,BToJpsi);
    RooDataHist roohist_BBbar("roohist_BBbar","roohist_BBbar",x,BBbar);

    RooPlot* frame_DPS = x.frame();

    DPS_DCB.fitTo(roohist_DPS,SumW2Error(kTRUE));
    roohist_DPS.plotOn(frame_DPS,LineColor(kBlue),MarkerColor(kBlue));
    DPS_DCB.plotOn(frame_DPS,LineColor(kBlue)); 
    //NonPromptx_DPS.fitTo(roohist_DPS,SumW2Error(kTRUE));
    //roohist_DPS.plotOn(frame_DPS,LineColor(kRed),MarkerColor(kRed));
    //NonPromptx_DPS.plotOn(frame_DPS,LineColor(kRed)); 


    NonPromptx_BBbar.fitTo(roohist_BBbar,SumW2Error(kTRUE));
    roohist_BBbar.plotOn(frame_DPS,LineColor(kRed),MarkerColor(kRed));
    NonPromptx_BBbar.plotOn(frame_DPS,LineColor(kRed)); 

    NonPromptx_BToJpsi.fitTo(roohist_BToJpsi,SumW2Error(kTRUE));
    roohist_BToJpsi.plotOn(frame_DPS,LineColor(kGreen),MarkerColor(kGreen));
    NonPromptx_BToJpsi.plotOn(frame_DPS,LineColor(kGreen)); 

    frame_DPS->SetTitle(""); 
    frame_DPS->GetXaxis()->SetTitle(axis[i]); 
    frame_DPS->GetYaxis()->SetTitle("Events"); 
    frame_DPS->SetMaximum(1000000);;
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
    //latex.DrawLatex(l+0.08, 1-t+lumiTextOffset*t-0.12,cmsText);//0.075
    latex.SetTextFont(extraTextFont);
    latex.SetTextSize(extraTextSize*t);
    //latex.DrawLatex(l+0.08, 1-t+lumiTextOffset*t-0.19, extraText);
    //latex.SetTextFont(channelTextFont);
    //latex.SetTextSize(channelTextSize);
    //latex.DrawLatex(l+0.165, 1-t+lumiTextOffset*t-0.19, channelText1);
    //latex.DrawLatex(l+0.165, 1-t+lumiTextOffset*t-0.28, channelText2);



    
    TLegend *pl2 = new TLegend(0.55,0.70,0.83,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(DPS, "DPS",  "PE");
    pl2->AddEntry(BToJpsi, "b #rightarrow J/#psi simulation events",  "PE");
    pl2->AddEntry(BBbar, "b#bar{b} #rightarrow J/#psi J/#psi simulation events",  "PE");
    pl2->AddEntry(DPS, "DPS", "L");
    pl2->AddEntry(BToJpsi, "b #rightarrow J/#psi fit", "L");
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
