#include "RooFit_data.h"
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

void RooFit_data_mass(){
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
  TCanvas* canvas = new TCanvas("canvas","canvas",0,0,1100,900);
  //canvas->SetLogy();
  TH1F *DPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *SPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *BBbar; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *Total; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *Run2018; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  RooRealVar x("x", "x", 2.8 , 3.4);
  x.setBins(1000,"cache");
  //RooRealVar mean("mean", "mean",0.00065,-0.002,0.002); 
  //RooRealVar sigma("sigma", "sigma",0.0037,0,0.01);
  //RooRealVar rlife("rlife", "rlife",0.016,0,0.5);
  float twomuon_cut = 0.05;
  float fourmuon_cut = 0.01;
  float Jpsi_mass_CR_do = 2.7;
  float Jpsi_mass_SR_do = 2.8;
  float Jpsi_mass_SR_up = 3.4;
  float Jpsi_mass_CR_up = 3.5;
  //Jpsi mass 
  RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi mean", 3.092); 
  RooRealVar Jpsi_width("Jpsi_width", "Jpsi width",0.032786, 0.01, 0.05);
  RooRealVar alpha_1("alpha_1", "",1.28, 0.3, 3);
  RooRealVar alpha_2("alpha_2", "",1.2424, 0.3, 3);
  RooRealVar n_1("n_1","",3.466, 1 ,100);
  RooRealVar n_2("n_2","",70, 1, 100);

  
  RooDoubleCB DCB_x("DCB_x","DCB_x",x,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);

  RooRealVar a0("a0","a0",-0.02387, -1, 1);
  RooRealVar a1("a1","a1",-0.02991, -1 ,1);

  RooChebychev comb_x("comb_x","comb_x",x,RooArgList(a0,a1));
 
  RooRealVar gfrac_Jpsi("gfrac_Jpsi","",0.6,0.5,0.85);
  
  RooAddPdf model_mass("model_x","model_x",RooArgList(DCB_x,comb_x),RooArgList(gfrac_Jpsi));
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
  RooRealVar DPS_width2("DPS_width2", "DPS width2",0.00369649, 0, 0.01);
  
  RooGaussian DPS_gauss1("gauss1","gaussian1 PDF",x,DPS_mean1,DPS_width1);
  RooGaussian DPS_gauss2("gauss2","gaussian2 PDF",x,DPS_mean1,DPS_width2);
  RooRealVar gfrac_DPS_gauss("gfrac_DPS_gauss","",0.801,0.5,0.9);
  RooAddPdf DPS_model("DPS_model","DPS_model",RooArgList(DPS_gauss1,DPS_gauss2),RooArgList(gfrac_DPS_gauss));
  
  //SPS double gauss fit
  
  RooRealVar SPS_mean1("SPS_mean1", "SPS mean1", 0.000042, 0, 0.001); 
  RooRealVar SPS_width1("SPS_width1", "SPS width1",0.00116, 0, 0.004);

  //RooRealVar SPS_mean2("SPS_mean2", "SPS mean2", 0,-0.02,0.04); 
  RooRealVar SPS_width2("SPS_width2", "SPS width2",0.003367, 0, 0.01);
  
  RooGaussian SPS_gauss1("gauss1","gaussian1 PDF",x,SPS_mean1,SPS_width1);
  RooGaussian SPS_gauss2("gauss2","gaussian2 PDF",x,SPS_mean1,SPS_width2);
  RooRealVar gfrac_SPS_gauss("gfrac_SPS_gauss","",0.703,0.5,0.9);
  //RooAddPdf SPS_model("SPS_model","SPS_model",RooArgList(SPS_gauss1,SPS_gauss2),RooArgList(gfrac_SPS_gauss));
  //RooAddPdf prompt_model("prompt_model","prompt_model",RooArgList(SPS_gauss1,SPS_gauss2),RooArgList(gfrac_SPS_gauss));
  
  //SPS fit function parameter
  
  RooRealVar SPS_mean("SPS_mean", "SPS mean", 0.000042, 0, 0.001); 
  RooRealVar SPS_width("SPS_width", "SPS width",0.01, 0, 0.08);
  RooRealVar SPS_alpha_1("SPS_alpha_1", "",1.00068, 0.8, 1.1);
  RooRealVar SPS_alpha_2("SPS_alpha_2", "",0.9571);
  RooRealVar SPS_n_1("SPS_n_1","",7.52, 5.0, 8.0);
  RooRealVar SPS_n_2("SPS_n_2","",6.2);
  //RooDoubleCB SPS_DCB("DCB","DCB",x,SPS_mean,SPS_width,SPS_alpha_1,SPS_n_1,SPS_alpha_2,SPS_n_2);
  //RooDoubleCB prompt_model("prompt_model","prompt_model",x,SPS_mean,SPS_width,SPS_alpha_1,SPS_n_1,SPS_alpha_2,SPS_n_2);
  RooDoubleCB prompt_model("prompt_model","prompt_model",x,SPS_mean,SPS_width,SPS_alpha_1,SPS_n_1,SPS_alpha_1,SPS_n_1);
  

  //Step Function
  RooRealVar x_StepValue("x_StepValue", "x_StepValue", 1);
  RooRealVar x_LowerLimit("x_LowerLimit", "x_LowerLimit", 0);
  RooRealVar x_HigherLimit("x_HigherLimit", "x_HigherLimit", 1);
  RooAbsReal* Stepx = new RooStepFunction("Stepx", "Stepx", x, RooArgList(x_StepValue), RooArgList(x_LowerLimit, x_HigherLimit));

  //BBbarToJpsiJpsi fit function parameter
    //Generic Exponential 
  RooRealVar BBbar_Tau("BBbar_Tau", "BBbar_Tau", -25, -29, -22);      //1

  RooAbsPdf* Expx_BBbar = new RooGenericPdf("BBbar", "BBbar", "Stepx * exp(x * BBbar_Tau)", RooArgSet(*Stepx, x, BBbar_Tau));
  
    //Gauss
  RooRealVar BBbar_Mean("BBbar_Mean", "BBbar_Mean",0.00171);
  RooRealVar BBbar_Sigma("BBbar_Sigma", "BBbar_Sigma", 0.00235);
  RooAbsPdf* Gauss_BBbar = new RooGaussian("Gauss_BBbar", "Gauss_BBbar", x, BBbar_Mean, BBbar_Sigma);
   
    //Convolution
  RooFFTConvPdf non_prompt_model("ConvPdf_BBbar", "ConvPdf_BBbar", x, *Expx_BBbar, *Gauss_BBbar);

  RooRealVar gfrac1("gfrac1","",0.5,0.2,0.85);
  RooAddPdf model_ctau("model_ctau","prompt_model+non_prompt_model",RooArgList(prompt_model,non_prompt_model),RooArgList(gfrac1));

  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  //name.push_back("MuMuFit_ctauPV_MC"); bin.push_back(200); Min.push_back(-0.2); Max.push_back(0.6); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("fourMuFit_ups1_cTau_MC"); bin.push_back(200); Min.push_back(-0.02); Max.push_back(0.2); axis.push_back("J/#psi c#tau [cm]");
  //name.push_back("Jpsi1_cTau"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.1); axis.push_back("J/#psi1 c#tau [cm]");
  //name.push_back("Jpsi2_cTau"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.1); axis.push_back("J/#psi2 c#tau [cm]");
  name.push_back("Jpsi1_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{12}[GeV]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    char CUT_Jpsi1_SR[1000];
    char CUT_Jpsi1_CR[1000];
    char CUT_Jpsi2_SR[1000];
    char CUT_Jpsi2_CR[1000];
    char CUT_Jpsi[1000];
    char CUT_Jpsi_Jpsi[1000];
    char CUT_mumu[1000];

    sprintf(CUT_Jpsi1_SR, "(Jpsi1_mass>%f && Jpsi1_mass<%f)", Jpsi_mass_SR_do, Jpsi_mass_SR_up);
    sprintf(CUT_Jpsi1_CR, "((Jpsi1_mass>%f && Jpsi1_mass<%f)||(Jpsi1_mass>%f && Jpsi1_mass<%f))", Jpsi_mass_CR_do, Jpsi_mass_SR_do, Jpsi_mass_SR_up, Jpsi_mass_SR_up);

    sprintf(CUT_Jpsi2_SR, "(Jpsi2_mass>%f && Jpsi2_mass<%f)", Jpsi_mass_SR_do, Jpsi_mass_SR_up);
    sprintf(CUT_Jpsi2_CR, "((Jpsi2_mass>%f && Jpsi2_mass<%f)||(Jpsi2_mass>%f && Jpsi2_mass<%f))", Jpsi_mass_CR_do, Jpsi_mass_SR_do, Jpsi_mass_SR_up, Jpsi_mass_SR_up);

    //sprintf(CUT_Jpsi,"(trigger_matched && %s )", CUT_Jpsi1_SR);
    sprintf(CUT_Jpsi_Jpsi,"(trigger_matched && %s && %s)", CUT_Jpsi1_SR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu,"(trigger_matched && %s )", CUT_Jpsi1_CR);


    //GetHisto(CUT_Jpsi, Tree_Run2018, Run2018 ,plot,bin[i],Min[i],Max[i]);
    GetHisto(CUT_Jpsi_Jpsi, Tree_Run2018, Run2018 ,plot,bin[i],Min[i],Max[i]);
       
    //cout<<DPS->Integral()<<endl;; 
    cout<<Run2018->Integral()<<endl;; 
   
    for(int j=1; j<Run2018->GetNbinsX()+1; j++){
      Run2018->SetBinContent(j,Run2018->GetBinContent(j));
      Run2018->SetBinError(j,sqrt(Run2018->GetBinContent(j)));
    }

    SPS = (TH1F*)Run2018->Clone();
    BBbar = (TH1F*)Run2018->Clone();
    Total = (TH1F*)Run2018->Clone();

    //DPS->SetLineColor(kBlue); DPS->SetLineWidth(2); DPS->SetMarkerColor(kBlue); DPS->SetMarkerSize(1.3); DPS->SetMarkerStyle(20);
    SPS->SetLineColor(kBlue); SPS->SetLineWidth(1);  SPS->SetLineStyle(1);
    BBbar->SetLineColor(kGreen); BBbar->SetLineWidth(1);  BBbar->SetLineStyle(7);
    Total->SetLineColor(kRed); Total->SetLineWidth(1);  Total->SetLineStyle(1);

    RooDataHist roohist_Run2018("roohist_Run2018","roohist_Run2018",x,Run2018);


    RooPlot* frame_Run2018 = x.frame();
    //model_ctau.fitTo(roohist_Run2018,SumW2Error(kTRUE));
    //roohist_Run2018.plotOn(frame_Run2018);
    //model_ctau.plotOn(frame_Run2018,Components(prompt_model),LineColor(kBlue)); 
    //model_ctau.plotOn(frame_Run2018,Components(non_prompt_model),LineColor(kGreen));
    //model_ctau.plotOn(frame_Run2018,LineColor(kRed));

    model_mass.fitTo(roohist_Run2018,SumW2Error(kTRUE));
    roohist_Run2018.plotOn(frame_Run2018);
    model_mass.plotOn(frame_Run2018,LineColor(kRed));

    RooHist *hpull = frame_Run2018->pullHist();
    RooPlot* x_Pull_Frame = x.frame();
    x_Pull_Frame->addPlotable(hpull,"P");

    model_mass.plotOn(frame_Run2018,Components(DCB_x),LineColor(kBlue)); 
    model_mass.plotOn(frame_Run2018,Components(comb_x),LineColor(kGreen),LineStyle(7));

    canvas->Divide(1,2);
    canvas->cd(1)->SetPad(0.01, 0.23, 0.99, 0.98);
    gPad->SetLeftMargin(0.1);      
    gPad->SetLogy();      
    frame_Run2018->SetTitle(""); 
    frame_Run2018->GetXaxis()->SetTitle(axis[i]); 
    frame_Run2018->GetXaxis()->SetLabelSize(0);
    frame_Run2018->GetYaxis()->SetTitle("Events"); 
    //frame_Run2018->SetMaximum(10000);
    frame_Run2018->SetMaximum(50000);
    frame_Run2018->SetMinimum(10);
    frame_Run2018->Draw();

    canvas->cd(2)->SetPad(0.01, 0.03, 0.99, 0.29);
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.27);
    //gPad->SetLogy(0);
    x_Pull_Frame->SetTitle(""); 
    x_Pull_Frame->GetYaxis()->SetTitle("Pull");
    x_Pull_Frame->GetYaxis()->SetTitleOffset(0.4);
    x_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
    x_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
    x_Pull_Frame->GetXaxis()->SetTitleSize(0.13);
    x_Pull_Frame->GetXaxis()->SetLabelSize(0.13);
    x_Pull_Frame->GetXaxis()->SetTitle(axis[i]); 
    x_Pull_Frame->Draw();
    TLine* line = new TLine(Min[i],0,Max[i],0);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw("same");

    canvas->cd(1);
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

    TString NAME = name[i];
    cout<<"Finished "<<NAME<<endl;
    TLegendEntry *ple2 = pl2->AddEntry(&roohist_Run2018, "data",  "PE");
    pl2->AddEntry(Total, "total fit", "L");
    pl2->AddEntry(SPS, "J/#psi fit", "L");
    pl2->AddEntry(BBbar, "comb. fit", "L");
    pl2->SetBorderSize(0);
    pl2->Draw();

    //TString NAME = name[i];
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
