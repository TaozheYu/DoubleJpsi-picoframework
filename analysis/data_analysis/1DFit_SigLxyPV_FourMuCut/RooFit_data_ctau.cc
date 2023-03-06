#include "RooFit_data_ctau.h"
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

void RooFit_data_ctau(){
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
  TH1F *Total; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *Run2018A_UL; 

  float twomuon_cut = 0.05;
  float fourmuon_cut = 0.01;
  float Jpsi_mass_do = 2.8;
  float Jpsi_mass_up = 3.4;
  RooRealVar x("x", "x", 0 , 10);
  x.setBins(1000,"cache");
  //Non-prompt and prompt pdf
    //Step Function
  RooRealVar x_StepValue("x_StepValue", "x_StepValue", 1);
  RooRealVar x_LowerLimit("x_LowerLimit", "x_LowerLimit", 0);
  RooRealVar x_HigherLimit("x_HigherLimit", "x_HigherLimit", 1000000);
  RooAbsReal* Stepx = new RooStepFunction("Stepx", "Stepx", x, RooArgList(x_StepValue), RooArgList(x_LowerLimit, x_HigherLimit));
  /*
  RooRealVar SPS_mean("SPS_mean", "SPS mean", 0);
  RooRealVar SPS_width("SPS_width", "SPS width",0.0017,0.0015,0.0025);
  RooRealVar SPS_alpha_1("SPS_alpha_1", "",1.2,0.7,2);
  RooRealVar SPS_alpha_2("SPS_alpha_2", "",1.2,0.7,2);
  RooRealVar SPS_n_1("SPS_n_1","",5.8,2,10);
  RooRealVar SPS_n_2("SPS_n_2","",4.2,3,5);
  RooDoubleCB prompt_model("DCB","DCB",x,SPS_mean,SPS_width,SPS_alpha_1,SPS_n_1,SPS_alpha_2,SPS_n_2);
  */
   //prompt
     //Generic Exponential 
   RooRealVar SPS_Tau("SPS_Tau", "SPS_Tau", -1.317, -2.0, -0.5);      //1
   RooAbsPdf* Expx_SPS = new RooGenericPdf("SPS", "SPS", "Stepx * exp(x * SPS_Tau)", RooArgSet(*Stepx, x, SPS_Tau));
   //Gauss
   RooRealVar SPS_Mean("SPS_Mean", "SPS_Mean",0.1693, 0.08, 2);
   RooRealVar SPS_Sigma("SPS_Sigma", "SPS_Sigma", 0.22525, 0.05, 2);
   RooAbsPdf* Gauss_SPS = new RooGaussian("Gauss_SPS", "Gauss_SPS", x, SPS_Mean, SPS_Sigma);
   
   //Convolution
   RooFFTConvPdf prompt_model("ConvPdf_SPS", "ConvPdf_SPS", x, *Expx_SPS, *Gauss_SPS);

  //Non-prompt
    //Generic Exponential
  RooRealVar BBbar_Tau("BBbar_Tau", "BBbar_Tau", -0.8418);      //1

  RooAbsPdf* Expx_BBbar = new RooGenericPdf("BBbar", "BBbar", "Stepx * exp(x * BBbar_Tau)", RooArgSet(*Stepx, x, BBbar_Tau));

    //Gauss
  RooRealVar BBbar_Mean("BBbar_Mean", "BBbar_Mean",0.5374);
  RooRealVar BBbar_Sigma("BBbar_Sigma", "BBbar_Sigma", 0.40001);
  RooAbsPdf* Gauss_BBbar = new RooGaussian("Gauss_BBbar", "Gauss_BBbar", x, BBbar_Mean, BBbar_Sigma);
    //Convolution
  RooFFTConvPdf non_prompt_model("ConvPdf_BBbar", "ConvPdf_BBbar", x, *Expx_BBbar, *Gauss_BBbar);

  RooRealVar gfrac1("gfrac1","",0.5,0.2,0.85);
  RooAddPdf model_ctau("model_ctau","prompt_model+non_prompt_model",RooArgList(prompt_model,non_prompt_model),RooArgList(gfrac1));

  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  //name.push_back("MuMuFit_ctauPV_noMC"); bin.push_back(50); Min.push_back(-0.02); Max.push_back(0.04); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(100); Min.push_back(-0.1); Max.push_back(0.12); axis.push_back("L_{xy}(J/#psi) [cm]");
  //name.push_back("fourMuFit_ups1_cTau_MC"); bin.push_back(200); Min.push_back(-0.02); Max.push_back(0.2); axis.push_back("J/#psi c#tau [cm]");
  name.push_back("fourMuFit_ups1_LxyPVSig_MC"); bin.push_back(100); Min.push_back(0); Max.push_back(5); axis.push_back("Significance of L_{xy}(J/#psi) [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    sprintf(CUT,"fourMuFit_VtxProb>%f && fourMuFit_ups1_mass>%f && fourMuFit_ups1_mass<%f && fourMuFit_ups2_mass>%f && fourMuFit_ups2_mass<%f", fourmuon_cut, Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
    //sprintf(CUT,"fourMuFit_ups1_mass>%f && fourMuFit_ups1_mass<%f && fourMuFit_ups2_mass>%f && fourMuFit_ups2_mass<%f", Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
    GetHisto(CUT, Tree_Run2018A_UL, Run2018A_UL ,plot,bin[i],Min[i],Max[i]);
       
    cout<<Run2018A_UL->Integral()<<endl;; 
   
    for(int j=1; j<Run2018A_UL->GetNbinsX()+1; j++){

      Run2018A_UL->SetBinContent(j,Run2018A_UL->GetBinContent(j));
      Run2018A_UL->SetBinError(j,sqrt(Run2018A_UL->GetBinContent(j)));

    }

    DPS = (TH1F*)Run2018A_UL->Clone(); 
    BBbar = (TH1F*)Run2018A_UL->Clone(); 
    Total = (TH1F*)Run2018A_UL->Clone(); 

    DPS->SetLineColor(kBlue); DPS->SetLineWidth(2); DPS->SetMarkerColor(kBlue); DPS->SetMarkerSize(1.3); DPS->SetMarkerStyle(20);
    BBbar->SetLineColor(kRed); BBbar->SetLineWidth(2); BBbar->SetMarkerColor(kRed); BBbar->SetMarkerSize(1.3); BBbar->SetMarkerStyle(20);
    Total->SetLineColor(kGreen); Total->SetLineWidth(2); Total->SetMarkerColor(kGreen); Total->SetMarkerSize(1.3); Total->SetMarkerStyle(20);
    
    RooDataHist roohist_Run2018A("roohist_Run2018A","roohist_Run2018A",x,Run2018A_UL);
    
    double total_num = Run2018A_UL->Integral();

    RooPlot* frame_Run2018A = x.frame();

    model_ctau.fitTo(roohist_Run2018A,SumW2Error(kTRUE));
    roohist_Run2018A.plotOn(frame_Run2018A);
    model_ctau.plotOn(frame_Run2018A,Components(prompt_model),LineColor(kBlue)); 
    model_ctau.plotOn(frame_Run2018A,Components(non_prompt_model),LineColor(kRed)); 
    model_ctau.plotOn(frame_Run2018A,LineColor(kGreen)); 

    Double_t frac1 = gfrac1.getVal();
    Double_t frac1_err = gfrac1.getError();
    cout<<"number of event is "<<total_num<<endl;
    cout<<"DPS frac is "<<frac1<<" +- "<<frac1_err<<endl;
    cout<<"DPS number is "<<total_num*frac1<<" +- "<<total_num*frac1_err<<endl;

    frame_Run2018A->SetTitle(""); 
    frame_Run2018A->GetXaxis()->SetTitle(axis[i]); 
    frame_Run2018A->GetYaxis()->SetTitle("Events"); 
    frame_Run2018A->SetMaximum(100000);;
    frame_Run2018A->SetMinimum(0.5);;
    frame_Run2018A->Draw();


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

    TLegend *pl2 = new TLegend(0.50,0.70,0.80,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(&roohist_Run2018A, "data",  "PE");
    pl2->AddEntry(DPS, "Prompt J/#psi fit", "L");
    pl2->AddEntry(BBbar, "Non-prompt J/#psi fit", "L");
    pl2->AddEntry(Total, "total fit", "L");
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[i];
    cout<<"Finished "<<NAME<<endl;
    canvas->SaveAs("Fit_plots/"+NAME+".pdf");
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

