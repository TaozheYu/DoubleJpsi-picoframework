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
  TH1F *SPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  //TH1F *BToJpsi; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
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
  

  //DPS double gauss fit
  
  RooRealVar DPS_mean1("DPS_mean1", "DPS mean1", 0.000055728); 
  RooRealVar DPS_width1("DPS_width1", "DPS width1",0.00174684, 0, 0.04);

  //RooRealVar DPS_mean2("DPS_mean2", "DPS mean2", 0,-0.02,0.04); 
  RooRealVar DPS_width2("DPS_width2", "DPS width2",0.00369649, 0, 0.05);
  
  RooGaussian DPS_gauss1("gauss1","gaussian1 PDF",x,DPS_mean1,DPS_width1);
  RooGaussian DPS_gauss2("gauss2","gaussian2 PDF",x,DPS_mean1,DPS_width2);
  RooRealVar gfrac_DPS_gauss("gfrac_DPS_gauss","",0.801,0.5,0.9);
  RooAddPdf DPS_model("DPS_model","DPS_model",RooArgList(DPS_gauss1,DPS_gauss2),RooArgList(gfrac_DPS_gauss));
  
  //SPS double gauss fit
  
  RooRealVar SPS_mean1("SPS_mean1", "SPS mean1", 0.000042); 
  RooRealVar SPS_width1("SPS_width1", "SPS width1",0.00116, 0, 0.004);

  //RooRealVar SPS_mean2("SPS_mean2", "SPS mean2", 0,-0.02,0.04); 
  RooRealVar SPS_width2("SPS_width2", "SPS width2",0.003367, 0, 0.005);
  
  RooGaussian SPS_gauss1("gauss1","gaussian1 PDF",x,SPS_mean1,SPS_width1);
  RooGaussian SPS_gauss2("gauss2","gaussian2 PDF",x,SPS_mean1,SPS_width2);
  RooRealVar gfrac_SPS_gauss("gfrac_SPS_gauss","",0.703,0.5,0.9);
  RooAddPdf SPS_model("SPS_model","SPS_model",RooArgList(SPS_gauss1,SPS_gauss2),RooArgList(gfrac_SPS_gauss));
  
  //SPS fit function parameter
  
  RooRealVar SPS_mean("SPS_mean", "SPS mean", 0); 
  RooRealVar SPS_width("SPS_width", "SPS width",0.01,0,0.08);
  RooRealVar SPS_alpha_1("SPS_alpha_1", "",1,0,10);
  RooRealVar SPS_alpha_2("SPS_alpha_2", "",1,0,10);
  RooRealVar SPS_n_1("SPS_n_1","",10,0,100);
  RooRealVar SPS_n_2("SPS_n_2","",10,0,500);
  RooDoubleCB SPS_DCB("DCB","DCB",x,SPS_mean,SPS_width,SPS_alpha_1,SPS_n_1,SPS_alpha_2,SPS_n_2);
  

  //Step Function
  RooRealVar x_StepValue("x_StepValue", "x_StepValue", 1);
  RooRealVar x_LowerLimit("x_LowerLimit", "x_LowerLimit", 0);
  RooRealVar x_HigherLimit("x_HigherLimit", "x_HigherLimit", 1000000);
  RooAbsReal* Stepx = new RooStepFunction("Stepx", "Stepx", x, RooArgList(x_StepValue), RooArgList(x_LowerLimit, x_HigherLimit));

  //BBbarToJpsiJpsi fit function parameter
    //Generic Exponential 
  RooRealVar BBbar_Tau("BBbar_Tau", "BBbar_Tau", -27.5);      //1

  RooAbsPdf* Expx_BBbar = new RooGenericPdf("BBbar", "BBbar", "Stepx * exp(x * BBbar_Tau)", RooArgSet(*Stepx, x, BBbar_Tau));
  
    //Gauss
  RooRealVar BBbar_Mean("BBbar_Mean", "BBbar_Mean",0.00171);
  RooRealVar BBbar_Sigma("BBbar_Sigma", "BBbar_Sigma", 0.00235);
  RooAbsPdf* Gauss_BBbar = new RooGaussian("Gauss_BBbar", "Gauss_BBbar", x, BBbar_Mean, BBbar_Sigma);
   
    //Convolution
  RooFFTConvPdf ConvPdf_BBbar("ConvPdf_BBbar", "ConvPdf_BBbar", x, *Expx_BBbar, *Gauss_BBbar);

  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  //name.push_back("MuMuFit_ctauPV_MC"); bin.push_back(200); Min.push_back(-0.2); Max.push_back(0.6); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(100); Min.push_back(-0.2); Max.push_back(0.5); axis.push_back("L_{xy}(J/#psi) [cm]");
  name.push_back("fourMuFit_ups1_cTau_MC"); bin.push_back(200); Min.push_back(-0.02); Max.push_back(0.2); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("fourMuFit_ups1_cTau_MC"); bin.push_back(80); Min.push_back(-0.02); Max.push_back(0.02); axis.push_back("J/#psi c#tau [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
    //sprintf(CUT,"");
    //sprintf(CUT,"fourMuFit_VtxProb>0.01");
    sprintf(CUT,"");
    GetHisto(CUT, Tree_DPS2018, DPS ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT, Tree_SPS2018, SPS ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT, Tree_BBbarToJpsi2018, BBbar ,plot,bin[i],Min[i],Max[i]);
       
    cout<<DPS->Integral()<<endl;; 
   
    for(int j=1; j<DPS->GetNbinsX()+1; j++){
      //data->SetBinContent(j,data_pre->GetBinContent(j));

      DPS->SetBinContent(j,DPS->GetBinContent(j));
      DPS->SetBinError(j,sqrt(DPS->GetBinContent(j)));

      SPS->SetBinContent(j,SPS->GetBinContent(j));
      SPS->SetBinError(j,sqrt(SPS->GetBinContent(j)));

      BBbar->SetBinContent(j,BBbar->GetBinContent(j));
      BBbar->SetBinError(j,sqrt(BBbar->GetBinContent(j)));
    }


    DPS->SetLineColor(kBlue); DPS->SetLineWidth(2); DPS->SetMarkerColor(kBlue); DPS->SetMarkerSize(1.3); DPS->SetMarkerStyle(20);
    SPS->SetLineColor(kGreen); SPS->SetLineWidth(2); SPS->SetMarkerColor(kGreen); SPS->SetMarkerSize(1.3); SPS->SetMarkerStyle(20);
    BBbar->SetLineColor(kRed); BBbar->SetLineWidth(2); BBbar->SetMarkerColor(kRed); BBbar->SetMarkerSize(1.3); BBbar->SetMarkerStyle(20);

    RooDataHist roohist_DPS("roohist_DPS","roohist_DPS",x,DPS);
    RooDataHist roohist_SPS("roohist_SPS","roohist_SPS",x,SPS);
    RooDataHist roohist_BBbar("roohist_BBbar","roohist_BBbar",x,BBbar);

    RooPlot* frame_DPS = x.frame();

    //DPS_model.fitTo(roohist_DPS,SumW2Error(kTRUE));
    //roohist_DPS.plotOn(frame_DPS,LineColor(kBlue),MarkerColor(kBlue));
    //DPS_model.plotOn(frame_DPS,LineColor(kBlue)); 

    DPS_DCB.fitTo(roohist_DPS,SumW2Error(kTRUE));
    roohist_DPS.plotOn(frame_DPS,LineColor(kBlue),MarkerColor(kBlue));
    DPS_DCB.plotOn(frame_DPS,LineColor(kBlue)); 
    
    //NonPromptx_DPS.fitTo(roohist_DPS,SumW2Error(kTRUE));
    //roohist_DPS.plotOn(frame_DPS,LineColor(kRed),MarkerColor(kRed));
    //NonPromptx_DPS.plotOn(frame_DPS,LineColor(kRed)); 

    //SPS_model.fitTo(roohist_SPS,SumW2Error(kTRUE));
    //roohist_SPS.plotOn(frame_DPS,LineColor(kGreen),MarkerColor(kGreen));
    //SPS_model.plotOn(frame_DPS,LineColor(kGreen)); 

    SPS_DCB.fitTo(roohist_SPS,SumW2Error(kTRUE));
    roohist_SPS.plotOn(frame_DPS,LineColor(kGreen),MarkerColor(kGreen));
    SPS_DCB.plotOn(frame_DPS,LineColor(kGreen)); 

    //NonPromptx_BToJpsi.fitTo(roohist_BToJpsi,SumW2Error(kTRUE));
    //roohist_BToJpsi.plotOn(frame_DPS,LineColor(kGreen),MarkerColor(kGreen));
    //NonPromptx_BToJpsi.plotOn(frame_DPS,LineColor(kGreen)); 

    ConvPdf_BBbar.fitTo(roohist_BBbar,SumW2Error(kTRUE));
    roohist_BBbar.plotOn(frame_DPS,LineColor(kRed),MarkerColor(kRed));
    ConvPdf_BBbar.plotOn(frame_DPS,LineColor(kRed)); 


    frame_DPS->SetTitle(""); 
    frame_DPS->GetXaxis()->SetTitle(axis[i]); 
    frame_DPS->GetYaxis()->SetTitle("Events"); 
    frame_DPS->SetMaximum(1000000);;
    frame_DPS->SetMinimum(1);;
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
    pl2->AddEntry(SPS, "SPS",  "PE");
    //pl2->AddEntry(BToJpsi, "b #rightarrow J/#psi simulation events",  "PE");
    pl2->AddEntry(BBbar, "b#bar{b} #rightarrow J/#psi J/#psi simulation events",  "PE");
    pl2->AddEntry(DPS, "DPS", "L");
    pl2->AddEntry(SPS, "SPS", "L");
    //pl2->AddEntry(BToJpsi, "b #rightarrow J/#psi fit", "L");
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
