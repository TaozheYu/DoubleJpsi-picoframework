#include "RooFit_data_Jpsi.h"
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

void RooFit4Dfunction_combine_data_Jpsi(){
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
  canvas->SetLogy();
  TH2F *data_2D_mass; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH2F *data_2D_ctau; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_x; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_y; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_z; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_u; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *BBbar; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *SPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *DPS; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);

  float twomuon_cut = 0.05;
  float fourmuon_cut = 0.01;
  float Jpsi_mass_do = 2.8;
  float Jpsi_mass_up = 3.4;

  RooRealVar x("x", "x", 2.80, 3.40);
  RooRealVar y("y", "y", 2.80, 3.40);
  RooRealVar z("z", "z", -0.02, 0.1);
  RooRealVar u("u", "u", -0.02, 0.1);
  z.setBins(1000,"cache");
  u.setBins(1000,"cache");
  z.setRange("R1", -0.02, 0.1);  
  u.setRange("R2", -0.02, 0.1);  

  //Construct Jpsi Jpsi pdf
  RooRealVar PDGmZ("PDGmZ","PDGmZ",3.092);
  RooRealVar PDGwZ("PDGwZ","PDGwZ",0.05,0.01,0.08);
  RooBreitWigner PDGBW("PDGBW","PDGBW",x,PDGmZ,PDGwZ);
  PDGmZ.setConstant(kTRUE);
  PDGwZ.setConstant(kTRUE);
  //Jpsi+Jpsi

  RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi mean", 3.092); 
  RooRealVar Jpsi_width("Jpsi_width", "Jpsi width",0.032786, 0.01, 0.05);
  RooRealVar alpha_1("alpha_1", "",1.28, 0.3, 3);
  RooRealVar alpha_2("alpha_2", "",1.2424, 0.3, 3);
  RooRealVar n_1("n_1","",3.466, 1 ,100);
  RooRealVar n_2("n_2","",70, 1, 100);

  
  RooDoubleCB DCB_x("DCB_x","DCB_x",x,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);
  RooDoubleCB DCB_y("DCB_y","DCB_y",y,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);

  //RooCBShape DCB_x("DCB_x","DCB_x",x,Jpsi_mean,Jpsi_width,alpha_1,n_1);
  //RooCBShape DCB_y("DCB_y","DCB_y",y,Jpsi_mean,Jpsi_width,alpha_2,n_2);

  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,DCB);
  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,CB_1);
 
 
  //Construct combinatorial pdf
  RooRealVar a0("a0","a0",-0.02387, -1, 1);
  RooRealVar a1("a1","a1",-0.02991, -1 ,1);

  RooChebychev comb_x("comb_x","comb_x",x,RooArgList(a0,a1));
  RooChebychev comb_y("comb_y","comb_y",y,RooArgList(a0,a1));
  
  //construct composite pdf
  RooRealVar gfrac_Jpsi("gfrac_Jpsi","",0.6,0.5,0.85);
  
  RooAddPdf model_x("model_x","model_x",RooArgList(DCB_x,comb_x),RooArgList(gfrac_Jpsi));
  RooAddPdf model_y("model_y","model_y",RooArgList(DCB_y,comb_y),RooArgList(gfrac_Jpsi));
  //RooAddPdf model_x("model_x","model_x",RooArgList(CB_x,comb_x),RooArgList(gfrac_Jpsi));
  //RooAddPdf model_y("model_y","model_y",RooArgList(CB_y,comb_y),RooArgList(gfrac_Jpsi));
  RooProdPdf model_mass("","model_x*model_y",RooArgList(model_x,model_y));

    //Step Function
  RooRealVar StepValue("StepValue", "StepValue", 1);
  RooRealVar LowerLimit("LowerLimit", "LowerLimit", 0);
  RooRealVar HigherLimit("HigherLimit", "HigherLimit", 10);
  RooAbsReal* Step_z = new RooStepFunction("Step_z", "Step_z", z, RooArgList(StepValue), RooArgList(LowerLimit, HigherLimit));
  RooAbsReal* Step_u = new RooStepFunction("Step_u", "Step_u", u, RooArgList(StepValue), RooArgList(LowerLimit, HigherLimit));

  //Non-prompt and prompt pdf
    //prompt fit
  /* 
  RooRealVar prompt_mean1("prompt_mean1", " mean1", 0.00007376); //0.00007376
  RooRealVar prompt_width1("prompt_width1", " width1", 0.00225);//0.00225

  RooRealVar prompt_width2("prompt_width2", " width2", 0.006728);//0.006728
    //z 
  RooGaussian prompt_gauss1_z("gauss1","gaussian1 PDF",z,prompt_mean1,prompt_width1);
  RooGaussian prompt_gauss2_z("gauss2","gaussian2 PDF",z,prompt_mean1,prompt_width2);
    //u
  RooGaussian prompt_gauss1_u("gauss1","gaussian1 PDF",u,prompt_mean1,prompt_width1);
  RooGaussian prompt_gauss2_u("gauss2","gaussian2 PDF",u,prompt_mean1,prompt_width2);
   //Add
  RooRealVar gfrac_gauss1("gfrac_gauss1","", 0.676);
  RooAddPdf prompt_model_z("prompt_model_z","prompt_model_z",RooArgList(prompt_gauss1_z,prompt_gauss2_z),RooArgList(gfrac_gauss1));
  RooAddPdf prompt_model_u("prompt_model_u","prompt_model_u",RooArgList(prompt_gauss1_u,prompt_gauss2_u),RooArgList(gfrac_gauss1));
  */
  
  RooRealVar prompt_mean("prompt_mean", "prompt mean", 0.0001);
  RooRealVar prompt_width("prompt_width", "prompt width",0.00225);
  RooRealVar prompt_alpha_1("prompt_alpha_1", "", 1.006);
  RooRealVar prompt_alpha_2("prompt_alpha_2", "", 0.95);
  RooRealVar prompt_n_1("prompt_n_1","", 3.426);
  RooRealVar prompt_n_2("prompt_n_2","", 3.476);
  RooDoubleCB prompt_model_z("prompt_model_z","prompt_model_z",z,prompt_mean,prompt_width,prompt_alpha_1,prompt_n_1,prompt_alpha_2,prompt_n_2);
  RooDoubleCB prompt_model_u("prompt_model_u","prompt_model_u",u,prompt_mean,prompt_width,prompt_alpha_1,prompt_n_1,prompt_alpha_2,prompt_n_2);
  
   //SPS Exponential
     //Generic Exponential 
   /*
   RooRealVar SPS_Tau("SPS_Tau", "SPS_Tau", -1.133);      //1
   RooAbsPdf* Exp_SPS_z = new RooGenericPdf("SPS", "SPS", "Step_z * exp(z * SPS_Tau)", RooArgSet(*Step_z, z, SPS_Tau));
   RooAbsPdf* Exp_SPS_u = new RooGenericPdf("SPS", "SPS", "Step_u * exp(u * SPS_Tau)", RooArgSet(*Step_u, u, SPS_Tau));

     //Gauss
   RooRealVar SPS_Mean("SPS_Mean", "SPS_Mean",0.3127);
   RooRealVar SPS_Sigma("SPS_Sigma", "SPS_Sigma", 0.22215);
   RooAbsPdf* Gauss_SPS_z = new RooGaussian("Gauss_SPS", "Gauss_SPS", z, SPS_Mean, SPS_Sigma);
   RooAbsPdf* Gauss_SPS_u = new RooGaussian("Gauss_SPS", "Gauss_SPS", u, SPS_Mean, SPS_Sigma);
    */
   //Convolution
   //RooFFTConvPdf prompt_model_z("ConvPdf_SPS_z", "ConvPdf_SPS_z", z, *Exp_SPS_z, *Gauss_SPS_z);
   //RooFFTConvPdf prompt_model_u("ConvPdf_SPS_u", "ConvPdf_SPS_u", u, *Exp_SPS_u, *Gauss_SPS_u);
   //RooNumConvPdf prompt_model_z("ConvPdf_SPS_z", "ConvPdf_SPS_z", z, *Exp_SPS_z, *Gauss_SPS_z);
   //RooNumConvPdf prompt_model_u("ConvPdf_SPS_u", "ConvPdf_SPS_u", u, *Exp_SPS_u, *Gauss_SPS_u);
    //Non-prompt fit function parameter
  /* 
  RooRealVar BBbar_mean("BBbar_mean", "BBbar mean", 0.00231);
  RooRealVar BBbar_width("BBbar_width", "BBbar width",0.00233);
  RooRealVar BBbar_alpha_1("BBbar_alpha_1", "",0.464);
  RooRealVar BBbar_alpha_2("BBbar_alpha_2", "",0.566);
  RooRealVar BBbar_n_1("BBbar_n_1","",125);
  RooRealVar BBbar_n_2("BBbar_n_2","",100);
  RooDoubleCB non_prompt_model_z("BBbar_DCB","BBbar_DCB",z,BBbar_mean,BBbar_width,BBbar_alpha_1,BBbar_n_1,BBbar_alpha_2,BBbar_n_2);
  RooDoubleCB non_prompt_model_u("BBbar_DCB","BBbar_DCB",u,BBbar_mean,BBbar_width,BBbar_alpha_1,BBbar_n_1,BBbar_alpha_2,BBbar_n_2);
  RooRealVar gfrac2("gfrac2","",0.51, 0.5, 0.52);
  //RooAddPdf model_ctau("model_ctau","prompt_model+non_prompt_model",RooArgList(prompt_model,non_prompt_model),RooArgList(gfrac2));
  */

    //Generic Exponential
  
  RooRealVar BBbar_Tau("BBbar_Tau", "BBbar_Tau", -26.3);      //-2.63

  RooAbsPdf* Exp_BBbar_z = new RooGenericPdf("BBbar_z", "BBbar_z", "Step_z * exp(z * BBbar_Tau)", RooArgSet(*Step_z, z, BBbar_Tau));
  RooAbsPdf* Exp_BBbar_u = new RooGenericPdf("BBbar_u", "BBbar_u", "Step_u * exp(u * BBbar_Tau)", RooArgSet(*Step_u, u, BBbar_Tau));

    //Gauss
  RooRealVar BBbar_Mean("BBbar_Mean", "BBbar_Mean",0.00171);  //0.00171
  RooRealVar BBbar_Sigma("BBbar_Sigma", "BBbar_Sigma", 0.00235); //0.00235
  RooAbsPdf* Gauss_BBbar_z = new RooGaussian("Gauss_BBbar_z", "Gauss_BBbar_z", z, BBbar_Mean, BBbar_Sigma);
  RooAbsPdf* Gauss_BBbar_u = new RooGaussian("Gauss_BBbar_u", "Gauss_BBbar_u", u, BBbar_Mean, BBbar_Sigma);
  
    //Convolution
  RooFFTConvPdf non_prompt_model_z("ConvPdf_BBbar_z", "ConvPdf_BBbar_z", z, *Exp_BBbar_z, *Gauss_BBbar_z);
  RooFFTConvPdf non_prompt_model_u("ConvPdf_BBbar_u", "ConvPdf_BBbar_u", u, *Exp_BBbar_u, *Gauss_BBbar_u);

  //RooNumConvPdf non_prompt_model_z("ConvPdf_BBbar_z", "ConvPdf_BBbar_z", z, *Exp_BBbar_z, *Gauss_BBbar_z);
  //RooNumConvPdf non_prompt_model_u("ConvPdf_BBbar_u", "ConvPdf_BBbar_u", u, *Exp_BBbar_u, *Gauss_BBbar_u);
  
  RooRealVar gfrac_prompt("gfrac_prompt","",0.518, 0.49, 0.6);
  //RooRealVar gfrac_SPS("gfrac_SPS","",0.0001);
  RooAddPdf model_z("model_z","prompt_model_z+non_prompt_model_z",RooArgList(prompt_model_z,non_prompt_model_z),RooArgList(gfrac_prompt));
  RooAddPdf model_u("model_u","prompt_model_u+non_prompt_model_u",RooArgList(prompt_model_u,non_prompt_model_u),RooArgList(gfrac_prompt));
  RooProdPdf model_ctau("model_ctau","model_z*model_u",RooArgList(model_z,model_u));
  //RooProdPdf model_test("","prompt_model_z*prompt_model_u",RooArgList(prompt_model_z,prompt_model_u));
  
  // sim
  //RooDataSet* data_ctau_sim = model_ctau.generate(RooArgSet(z,u), 12475);
   
  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  name.push_back("Jpsi1_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{12}[GeV]");
  name.push_back("Jpsi2_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{34}[GeV]");
  name.push_back("Jpsi1_cTau"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.1); axis.push_back("ctau_{J/#psi1} [cm]");
  name.push_back("Jpsi2_cTau"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.1); axis.push_back("ctau_{J/#psi2} [cm]");

  char CUTpre[1000];
  char CUT[1000];
  sprintf(CUTpre,"(trigger_matched && Jpsi1_mass>%f && Jpsi1_mass<%f && Jpsi2_mass>%f && Jpsi2_mass<%f)",  Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
  sprintf(CUT,"1 *%s ", CUTpre);
  /*
  const char *plot = name[2];
  GetHisto(CUT, Tree_BBbarToJpsi2018, BBbar ,plot,bin[2],Min[2],Max[2]);
  RooDataHist roohist_BBbar_z("roohist_BBbar_z","roohist_BBbar_z",z,BBbar);
  RooDataHist roohist_BBbar_u("roohist_BBbar_u","roohist_BBbar_u",u,BBbar);
  Int_t N = 0;
  RooHistPdf non_prompt_model_z("non_prompt_model_z","non_prompt_model_z",z, roohist_BBbar_z,N);
  RooHistPdf non_prompt_model_u("non_prompt_model_u","non_prompt_model_u",u, roohist_BBbar_u,N);

  GetHisto(CUT, Tree_SPS2018, SPS ,plot,bin[2],Min[2],Max[2]);
  RooDataHist roohist_SPS_z("roohist_SPS_z","roohist_SPS_z",z,SPS);
  RooDataHist roohist_SPS_u("roohist_SPS_u","roohist_SPS_u",u,SPS);
  RooHistPdf SPS_z("SPS_z","SPS_z",z, roohist_SPS_z,N);
  RooHistPdf SPS_u("SPS_u","SPS_u",u, roohist_SPS_u,N);

  GetHisto(CUT, Tree_DPS2018, DPS ,plot,bin[2],Min[2],Max[2]);
  RooDataHist roohist_DPS_z("roohist_DPS_z","roohist_DPS_z",z,DPS);
  RooDataHist roohist_DPS_u("roohist_DPS_u","roohist_DPS_u",u,DPS);
  RooHistPdf DPS_z("DPS_z","DPS_z",z, roohist_DPS_z,N);
  RooHistPdf DPS_u("DPS_u","DPS_u",u, roohist_DPS_u,N);

  RooAddPdf prompt_model_z("model_z","SPS_z+DPS_z",RooArgList(SPS_z,DPS_z),RooArgList(gfrac_SPS));
  RooAddPdf prompt_model_u("model_u","SPS_u+DPS_u",RooArgList(SPS_u,DPS_u),RooArgList(gfrac_SPS));
  */
  //RooAddPdf model_z("model_z","prompt_model_z+non_prompt_model_z",RooArgList(prompt_model_z,non_prompt_model_z),RooArgList(gfrac2));
  //RooAddPdf model_u("model_u","prompt_model_u+non_prompt_model_u",RooArgList(prompt_model_u,non_prompt_model_u),RooArgList(gfrac2));
  //RooProdPdf model_ctau("model_ctau","model_z*model_u",RooArgList(model_z,model_u));

  
  RooProdPdf model_total1("model_total1","model_x*model_y*model_z*model_u",RooArgList(model_x,model_y,model_z,model_u));

  //same pdf function but different component
  
  RooProdPdf model_prompt_prompt("prompt_prompt","prompt_z*prompt_u",RooArgList(prompt_model_z,prompt_model_u));
  RooProdPdf model_prompt_nonprompt("prompt_nonprompt","prompt_z*nonprompt_u",RooArgList(prompt_model_z,non_prompt_model_u));
  RooProdPdf model_nonprompt_prompt("nonprompt_prompt","nonprompt_z*prompt_u",RooArgList(non_prompt_model_z,prompt_model_u));
  RooProdPdf model_nonprompt_nonprompt("nonprompt_nonprompt","nonprompt_z*nonprompt_u",RooArgList(non_prompt_model_z,non_prompt_model_u));

  RooRealVar gfrac_prompt_prompt("gfrac_prompt_prompt","gfrac_prompt_prompt", 0.25, 0.2, 0.3);  
  RooRealVar gfrac_nonprompt_prompt("gfrac_nonprompt_prompt","gfrac_nonprompt_prompt", 0.25, 0.2, 0.3);  
  RooRealVar gfrac_PmixNP("PmixNP","PmixNP", 0.5, 0.4, 0.6);  
  RooRealVar gfrac_nonprompt_nonprompt("gfrac_nonprompt_nonprompt","gfrac_nonprompt_nonprompt", 0.25, 0.2, 0.3);  
  //RooAddPdf model_ctau2("model_ctau2","model_ctau2", RooArgList(model_prompt_prompt, model_prompt_nonprompt, model_nonprompt_prompt, model_nonprompt_nonprompt),
   //RooArgList(gfrac_prompt_prompt, gfrac_nonprompt_prompt, gfrac_nonprompt_prompt));



  //RooDataSet* data_nonprompt_nonprompt_sim = model_nonprompt_nonprompt.generate(RooArgSet(z,u), 12475);
  //RooDataSet* data_prompt_nonprompt_sim = model_prompt_nonprompt.generate(RooArgSet(z,u), 12475);
  //RooDataSet* data_nonprompt_prompt_sim = model_nonprompt_prompt.generate(RooArgSet(z,u), 12475);
  //RooDataSet* data_prompt_prompt_sim = model_prompt_prompt.generate(RooArgSet(z,u), 12475);
  
  RooProdPdf model_signal_prompt_prompt("model_signal_prompt_prompt","DCB_x*DCB_y*prompt_z*prompt_u",RooArgList(DCB_x,DCB_y,prompt_model_z,prompt_model_u));

  RooProdPdf model_signal_prompt_nonprompt("model_signal_prompt_nonprompt","DCB_x*DCB_y*prompt_z*nonprompt_u",RooArgList(DCB_x,DCB_y,prompt_model_z,non_prompt_model_u));
  RooProdPdf model_signal_nonprompt_prompt("model_signal_nonprompt_prompt","DCB_x*DCB_y*nonprompt_z*prompt_u",RooArgList(DCB_x,DCB_y,non_prompt_model_z,prompt_model_u));
  RooRealVar gfrac_signal_PmixNP("PmixNP","PmixNP", 0.5);  
  RooAddPdf model_signal_PmixNP("model_signal_PmixNP","model_signal_prompt_nonprompt+model_signal_nonprompt_prompt",RooArgList(model_signal_prompt_nonprompt,model_signal_nonprompt_prompt), RooArgList(gfrac_signal_PmixNP));  

  RooProdPdf model_signal_nonprompt_nonprompt("model_signal_nonprompt_nonprompt","DCB_x*DCB_y*nonprompt_z*nonprompt_u",RooArgList(DCB_x,DCB_y,non_prompt_model_z,non_prompt_model_u));

  RooProdPdf model_DCB1_comb2("model_bkg1","DCB_x*comb_y*model_ctau",RooArgList(DCB_x,comb_y,model_ctau));
  RooProdPdf model_comb1_DCB2("model_bkg2","comb_x*DCB_y*model_ctau",RooArgList(comb_x,DCB_y,model_ctau));
  RooRealVar gfrac_DCBmixComb("","", 0.5);  
  RooAddPdf model_DCBmixComb("DCBmixComb","model_DCB1_comb2+model_comb1_DCB2",RooArgList(model_DCB1_comb2,model_comb1_DCB2), RooArgList(gfrac_DCBmixComb));  
  RooProdPdf model_comb1_comb2("model_bkg3","comb_x*comb_y*model_ctau",RooArgList(comb_x,comb_y,model_ctau));
  
  for(int i=0; i<name.size(); i++){
   for(int j=i+1; j<name.size(); j++){
    for(int k=j+1; k<name.size(); k++){
     for(int l=k+1; l<name.size(); l++){
        const char *plotX = name[i];
	const char *plotY = name[j];
        const char *plotZ = name[k];
        const char *plotU = name[l];
	int   XBIN = bin[i];   int YBIN = bin[j];   int ZBIN = bin[k];   int UBIN = bin[l];
	float XMIN = Min[i]; float YMIN = Min[j]; float ZMIN = Min[k]; float UMIN = Min[l];
	float XMAX = Max[i]; float YMAX = Max[j]; float ZMAX = Max[k]; float UMAX = Max[l];
	//GetHisto(CUT, Tree, data_x ,plotX,bin[i],Min[i],Max[i]);
	//GetHisto(CUT, Tree, data_y ,plotY,bin[j],Min[j],Max[j]);
        GetHisto(CUT, Tree, data_z ,plotZ,bin[k],Min[k],Max[k]);
        GetHisto(CUT, Tree, data_u ,plotU,bin[l],Min[l],Max[l]);
	GetHisto(CUT, Tree, data ,plotX,bin[i],Min[i],Max[i]);

        GetHisto2D(CUT, Tree, data_2D_mass ,plotX,plotY,XBIN,XMIN,XMAX,YBIN,YMIN,YMAX);
        GetHisto2D(CUT, Tree, data_2D_ctau ,plotZ,plotU,ZBIN,ZMIN,ZMAX,UBIN,UMIN,UMAX);
        cout<<data_2D_mass->Integral()<<endl;
        cout<<data_2D_ctau ->Integral()<<endl;
        /*
	for(int n=1; n<data_2D_mass->GetNbinsX()+1; n++){
          for(int m=1; m<data_2D_mass->GetNbinsY()+1; m++){
              data_2D_mass->SetBinContent(n,m,data_2D_mass->GetBinContent(n,m));
              data_2D_mass->SetBinError(n,m,sqrt(data_2D_mass->GetBinContent(n,m)));

            }
          }
        
	for(int p=1; p<data_2D_ctau->GetNbinsX()+1; p++){
          for(int q=1; q<data_2D_ctau->GetNbinsY()+1; q++){
              data_2D_ctau->SetBinContent(p,q,data_2D_ctau->GetBinContent(p,q));
              data_2D_ctau->SetBinError(p,q,sqrt(data_2D_ctau->GetBinContent(p,q)));

            }
          }
       */ 
     }
    }
   }
  }
      //data->Draw("LEGO");
      //canvas->SaveAs("data.pdf") ;
  data->SetMarkerColor(1); data->SetMarkerStyle(20); data->SetMarkerSize(1.3);
  TH1F *h_DCB_DCB_prompt_prompt; h_DCB_DCB_prompt_prompt = (TH1F*)data->Clone(); h_DCB_DCB_prompt_prompt->SetLineColor(kBlue); h_DCB_DCB_prompt_prompt->SetLineWidth(1); 
  TH1F *h_DCB_DCB_prompt_nonprompt; h_DCB_DCB_prompt_nonprompt = (TH1F*)data->Clone(); h_DCB_DCB_prompt_nonprompt->SetMarkerSize(0); h_DCB_DCB_prompt_nonprompt->SetLineColor(kGreen); h_DCB_DCB_prompt_nonprompt->SetLineWidth(2);
  TH1F *h_DCB_DCB_nonprompt_prompt; h_DCB_DCB_nonprompt_prompt = (TH1F*)data->Clone(); h_DCB_DCB_nonprompt_prompt->SetLineColor(kViolet-2); h_DCB_DCB_nonprompt_prompt->SetLineWidth(1);
  TH1F *h_DCB_DCB_nonprompt_nonprompt; h_DCB_DCB_nonprompt_nonprompt = (TH1F*)data->Clone(); h_DCB_DCB_nonprompt_nonprompt->SetLineColor(kOrange); h_DCB_DCB_nonprompt_nonprompt->SetLineWidth(1);  

  TH1F *h_DCB_comb; h_DCB_comb = (TH1F*)data->Clone();  h_DCB_comb->SetLineColor(1); h_DCB_comb->SetLineWidth(2); h_DCB_comb->SetLineStyle(7);
  TH1F *h_comb_DCB; h_comb_DCB = (TH1F*)data->Clone(); h_comb_DCB->SetLineColor(kViolet-2); h_comb_DCB->SetLineWidth(1); h_comb_DCB->SetLineStyle(7);
  TH1F *h_comb_comb; h_comb_comb = (TH1F*)data->Clone(); h_comb_comb->SetLineColor(kOrange); h_comb_comb->SetLineWidth(1); h_comb_comb->SetLineStyle(7);
  //TH1F *h_prompt; h_prompt = (TH1F*)data->Clone(); h_prompt->SetLineColor(1); h_prompt->SetLineStyle(1);
  //TH1F *h_nonprompt; h_nonprompt = (TH1F*)data->Clone(); h_nonprompt->SetLineColor(1); h_nonprompt->SetLineStyle(7);
  TH1F *h_total; h_total = (TH1F*)data->Clone(); h_total->SetLineColor(kRed); h_total->SetLineWidth(1);


  RooDataHist roohist_data_2D_mass("roohist_data","roohist_data",RooArgSet(x,y),data_2D_mass);
  RooDataHist roohist_data_2D_ctau("roohist_data_ctau","roohist_data_ctau",RooArgSet(z,u),data_2D_ctau);
  RooDataHist roohist_data_z("roohist_data_z","roohist_data_z",RooArgSet(z),data_z);
  RooDataHist roohist_data_u("roohist_data_u","roohist_data_u",RooArgSet(u),data_u);
  RooPlot* xframe_data = x.frame();
  RooPlot* yframe_data = y.frame();
  RooPlot* zframe_data = z.frame();
  RooPlot* uframe_data = u.frame();
  //data_2D_ctau->Draw("LEGO");
  //data_2D_ctau->Draw("SURF");
  //data_2D_mass->Draw("LEGO");
  //data_2D_ctau->Draw("COLZ");
  //data_2D_mass->Draw("COLZ");
  //canvas->SaveAs("data2D_ctau.pdf");
  //canvas->SaveAs("data2D_mass.pdf");
  //data_z->Draw("PE");
  //data_z->Draw("COLZ");
  //canvas->SaveAs("data_ctau.pdf");
  //2D+2DD fit
  model_total1.fitTo(roohist_data_2D_mass,SumW2Error(kTRUE));
  model_total1.fitTo(roohist_data_2D_ctau,SumW2Error(kTRUE));
  //roohist_data_2D_mass.plotOn(xframe_data);
  //model_total1.plotOn(xframe_data, LineColor(kRed), LineWidth(1)); //Normalization(12475, RooAbsReal::Relative),  
  //xframe_data->Draw();
  //canvas->SaveAs("data_fit_mass.pdf");
  //roohist_data_2D_ctau.plotOn(zframe_data);
  //model_total1.plotOn(zframe_data, LineColor(kRed), LineWidth(1)); //Normalization(12475, RooAbsReal::Relative),  
  //zframe_data->Draw();
  //canvas->SaveAs("data_fit_ctau.pdf");
  //model_ctau.fitTo(roohist_data_2D_ctau);//,SumW2Error(kTRUE)
  //model_ctau.fitTo(*data_ctau_sim ,SumW2Error(kTRUE));//,SumW2Error(kTRUE)
  //model_ctau.fitTo(*data_ctau_sim ,Extended(kTRUE));//,SumW2Error(kTRUE)
  //model_ctau.fitTo(roohist_data_2D_ctau);//,SumW2Error(kTRUE)
  //roohist_data_2D_ctau.plotOn(zframe_data);
  //data_ctau_sim->plotOn(zframe_data);
  //model_ctau2.plotOn(zframe_data, LineColor(kRed), LineWidth(1)); //Normalization(12475, RooAbsReal::Relative),  
  //model_ctau.plotOn(zframe_data, ProjectionRange("R2") ,LineColor(kRed), LineWidth(1)); 
  //zframe_data->Draw();
  //model_u.fitTo(roohist_data_u,SumW2Error(kTRUE));
  //roohist_data_u.plotOn(uframe_data);
  //model_u.plotOn(uframe_data,LineColor(kRed),LineWidth(1)); 
  //uframe_data->Draw();
  //canvas->SaveAs("data_fit_ctau.pdf");
  //canvas->SaveAs("data_sim_ctau.pdf");
  //model_mass.fitTo(roohist_data_2D_mass,SumW2Error(kTRUE));
  //roohist_data_2D_mass.plotOn(xframe_data);
  //model_mass.plotOn(xframe_data, LineColor(kRed), LineWidth(1)); //, Normalization(12475, RooAbsReal::NumEvent)
  //xframe_data->Draw();
  //canvas->SaveAs("data_fit_mass.pdf");
  /*
  model_nonprompt_nonprompt.fitTo(*data_nonprompt_nonprompt_sim ,Extended(kTRUE));//,SumW2Error(kTRUE)
  data_nonprompt_nonprompt_sim->plotOn(zframe_data);
  model_nonprompt_nonprompt.plotOn(zframe_data, LineColor(kRed), LineWidth(1)); 
  zframe_data->Draw();
  canvas->SaveAs("data_sim_nonprompt_nonprompt.pdf");
  */
  /*  
  model_prompt_nonprompt.fitTo(*data_prompt_nonprompt_sim ,Extended(kTRUE));//,SumW2Error(kTRUE)
  data_prompt_nonprompt_sim->plotOn(zframe_data);
  model_prompt_nonprompt.plotOn(zframe_data, LineColor(kRed), LineWidth(1)); 
  zframe_data->Draw();
  canvas->SaveAs("data_sim_prompt_nonprompt.pdf");
  */
  /* 
  model_nonprompt_prompt.fitTo(*data_nonprompt_prompt_sim ,Extended(kTRUE));//,SumW2Error(kTRUE)
  data_nonprompt_prompt_sim->plotOn(zframe_data);
  model_nonprompt_prompt.plotOn(zframe_data, LineColor(kRed), LineWidth(1)); 
  zframe_data->Draw();
  canvas->SaveAs("data_sim_nonprompt_prompt.pdf");
  */
  /*  
  model_prompt_prompt.fitTo(*data_prompt_prompt_sim ,Extended(kTRUE));//,SumW2Error(kTRUE)
  data_prompt_prompt_sim->plotOn(zframe_data);
  model_prompt_prompt.plotOn(zframe_data,  LineColor(kRed), LineWidth(1)); 
  zframe_data->Draw();
  canvas->SaveAs("data_sim_prompt_prompt.pdf");
  */
   
  Double_t frac_Jpsi = gfrac_Jpsi.getVal();
  Double_t frac_Jpsi_err = gfrac_Jpsi.getError(); 
  Double_t frac_prompt = gfrac_prompt.getVal();
  Double_t frac_prompt_err = gfrac_prompt.getError(); 
  cout<<"frac Jpsi is :"<<frac_Jpsi<<endl; 
  cout<<"frac Jpsi err is :"<<frac_Jpsi_err<<endl; 
  cout<<"frac prompt is :"<<frac_prompt<<endl; 
  cout<<"frac prompt err is :"<<frac_prompt_err<<endl; 
  //parameter (g1+(1-g1))*(g1+(1-g1))*(g2+(1-g2))*(g2+(1-g2))

  Double_t frac1111 = frac_Jpsi*frac_Jpsi*frac_prompt*frac_prompt;
  Double_t frac1112 = frac_Jpsi*frac_Jpsi*frac_prompt*(1-frac_prompt); 
  Double_t frac1121 = frac_Jpsi*frac_Jpsi*(1-frac_prompt)*frac_prompt; 
  Double_t frac1122 = frac_Jpsi*frac_Jpsi*(1-frac_prompt)*(1-frac_prompt); 
  Double_t frac11xx = frac1112 + frac1121;

  Double_t frac12 = frac_Jpsi*(1-frac_Jpsi); 
  Double_t frac21= (1-frac_Jpsi)*frac_Jpsi; 
  Double_t frac22 = (1-frac_Jpsi)*(1-frac_Jpsi); 
  Double_t fracxx = frac12 + frac21;

  RooRealVar gfrac1111("gfrac1111","",frac1111);  
  RooRealVar gfrac1112("gfrac1112","",frac1112);  
  RooRealVar gfrac1121("gfrac1121","",frac1121);  
  RooRealVar gfrac1122("gfrac1122","",frac1122);  
  RooRealVar gfrac11xx("gfrac11xx","",frac11xx);  

  RooRealVar gfrac12("gfrac12","",frac12);  
  RooRealVar gfrac21("gfrac21","",frac21);  
  RooRealVar gfrac22("gfrac22","",frac22);  
  RooRealVar gfracxx("gfracxx","",fracxx);  

  //RooAddPdf model_total2("model_total2","model_total2", RooArgList(model_signal_prompt_prompt, model_signal_prompt_nonprompt, model_signal_nonprompt_prompt, model_signal_nonprompt_nonprompt,
  //                       model_DCB1_comb2, model_comb1_DCB2, model_comb1_comb2),
  // RooArgList(gfrac1111, gfrac1112, gfrac1121, gfrac1122, gfrac12, gfrac21));
  
  RooAddPdf model_total2("model_total2","model_total2", RooArgList(model_signal_prompt_prompt, model_signal_PmixNP, model_signal_nonprompt_nonprompt,
                         model_DCBmixComb, model_comb1_comb2),
   RooArgList(gfrac1111, gfrac11xx, gfrac1122, gfracxx));

  //RooAddPdf model_total2("model_total2","model_total2", RooArgList(model_prompt_prompt, model_prompt_nonprompt, model_nonprompt_prompt, model_nonprompt_nonprompt), RooArgList(gfrac11, gfrac12, gfrac21));

  //model_total2.fitTo(roohist_data_2D_ctau,SumW2Error(kTRUE));
  //model_total2.fitTo(roohist_data_2D_mass,SumW2Error(kTRUE));
  //gStyle->SetPadTickY(1); 
  canvas->Divide(1,2);
  for(int k=0; k<name.size(); k++){
  //for(int k=0; k<1; k++){
    if(k==0){
      roohist_data_2D_mass.plotOn(xframe_data);
      model_total2.plotOn(xframe_data,LineColor(kRed),LineWidth(1));
      RooHist *hpull = xframe_data->pullHist();
      RooPlot* x_Pull_Frame = x.frame();
      x_Pull_Frame->addPlotable(hpull,"P");
      model_total2.plotOn(xframe_data,LineColor(kRed),LineWidth(1)); 
 
      model_total2.plotOn(xframe_data,Components(model_signal_prompt_prompt),LineColor(kBlue),LineWidth(1)); 

      //model_total2.plotOn(xframe_data,Components(model_signal_prompt_nonprompt),LineColor(kGreen),LineWidth(2)); 

      //model_total2.plotOn(xframe_data,Components(model_signal_nonprompt_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(xframe_data,Components(model_signal_PmixNP),LineColor(kViolet-2),LineWidth(1)); 

      model_total2.plotOn(xframe_data,Components(model_signal_nonprompt_nonprompt),LineColor(kOrange),LineWidth(1)); 

      //model_total2.plotOn(xframe_data,Components(model_DCB1_comb2),LineColor(1),LineStyle(7),LineWidth(2)); 

      //model_total2.plotOn(xframe_data,Components(model_comb1_DCB2),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 
      model_total2.plotOn(xframe_data,Components(model_DCBmixComb),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(xframe_data,Components(model_comb1_comb2),LineColor(kOrange),LineStyle(7),LineWidth(1)); 

      canvas->cd(1)->SetPad(0.01, 0.23, 0.99, 0.98);
      gPad->SetLeftMargin(0.1);      
      gPad->SetLogy();      
      xframe_data->SetTitle(""); 
      xframe_data->GetXaxis()->SetTitle(axis[k]); 
      xframe_data->GetXaxis()->SetLabelSize(0);
      xframe_data->GetYaxis()->SetTitle("Events"); 
      xframe_data->GetYaxis()->SetTitleOffset(1.30);
      xframe_data->SetMaximum(50000);
      xframe_data->SetMinimum(10);
      xframe_data->Draw();
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
      x_Pull_Frame->GetXaxis()->SetTitle(axis[k]); 
      x_Pull_Frame->Draw();
      TLine* line = new TLine(Min[k],0,Max[k],0);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw("same");
      } 

    if(k==1){
      roohist_data_2D_mass.plotOn(yframe_data);
      model_total2.plotOn(yframe_data,LineColor(kRed),LineWidth(1));
      RooHist *hpull = yframe_data->pullHist();
      RooPlot* y_Pull_Frame = y.frame();
      y_Pull_Frame->addPlotable(hpull,"P");
      model_total2.plotOn(yframe_data,LineColor(kRed),LineWidth(1)); 
 
      model_total2.plotOn(yframe_data,Components(model_signal_prompt_prompt),LineColor(kBlue),LineWidth(1)); 

      //model_total2.plotOn(yframe_data,Components(model_signal_prompt_nonprompt),LineColor(kGreen),LineWidth(2)); 

      //model_total2.plotOn(yframe_data,Components(model_signal_nonprompt_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(yframe_data,Components(model_signal_PmixNP),LineColor(kViolet-2),LineWidth(1)); 

      model_total2.plotOn(yframe_data,Components(model_signal_nonprompt_nonprompt),LineColor(kOrange),LineWidth(1)); 

      //model_total2.plotOn(yframe_data,Components(model_DCB1_comb2),LineColor(1),LineStyle(7),LineWidth(2)); 

      //model_total2.plotOn(yframe_data,Components(model_comb1_DCB2),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 
      model_total2.plotOn(yframe_data,Components(model_DCBmixComb),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(yframe_data,Components(model_comb1_comb2),LineColor(kOrange),LineStyle(7),LineWidth(1)); 

      canvas->cd(1)->SetPad(0.01, 0.23, 0.99, 0.98);
      gPad->SetLeftMargin(0.1);      
      gPad->SetLogy();      
      yframe_data->SetTitle(""); 
      yframe_data->GetXaxis()->SetTitle(axis[k]); 
      yframe_data->GetXaxis()->SetLabelSize(0);
      yframe_data->GetYaxis()->SetTitle("Events"); 
      yframe_data->GetYaxis()->SetTitleOffset(1.30);
      yframe_data->SetMaximum(50000);
      yframe_data->SetMinimum(10);
      yframe_data->Draw();
      canvas->cd(2)->SetPad(0.01, 0.03, 0.99, 0.29);
      gPad->SetLeftMargin(0.1);
      gPad->SetBottomMargin(0.27);
      gPad->SetLogy(0);
      y_Pull_Frame->SetTitle(""); 
      y_Pull_Frame->GetYaxis()->SetTitle("Pull");
      y_Pull_Frame->GetYaxis()->SetTitleOffset(0.4);
      y_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
      y_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
      y_Pull_Frame->GetXaxis()->SetTitleSize(0.13);
      y_Pull_Frame->GetXaxis()->SetLabelSize(0.13);
      y_Pull_Frame->GetXaxis()->SetTitle(axis[k]); 
      y_Pull_Frame->Draw();
      TLine* line = new TLine(Min[k],0,Max[k],0);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw("same");
    }
    if(k==2){
      roohist_data_2D_ctau.plotOn(zframe_data);
      model_total2.plotOn(zframe_data,LineColor(kRed),LineWidth(1));
      RooHist *hpull = zframe_data->pullHist();
      RooPlot* z_Pull_Frame = z.frame();
      z_Pull_Frame->addPlotable(hpull,"P");

      model_total2.plotOn(zframe_data,Components(model_signal_prompt_prompt),LineColor(kBlue),LineWidth(1)); 

      //model_total2.plotOn(zframe_data,Components(model_signal_prompt_nonprompt),LineColor(kGreen),LineWidth(2)); 

      //model_total2.plotOn(zframe_data,Components(model_signal_nonprompt_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(zframe_data,Components(model_signal_PmixNP),LineColor(kViolet-2),LineWidth(1)); 

      model_total2.plotOn(zframe_data,Components(model_signal_nonprompt_nonprompt),LineColor(kOrange),LineWidth(1)); 

      //model_total2.plotOn(zframe_data,Components(model_DCB1_comb2),LineColor(1),LineStyle(7),LineWidth(2)); 

      //model_total2.plotOn(zframe_data,Components(model_comb1_DCB2),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 
      model_total2.plotOn(zframe_data,Components(model_DCBmixComb),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(zframe_data,Components(model_comb1_comb2),LineColor(kOrange),LineStyle(7),LineWidth(1)); 
      canvas->cd(1)->SetPad(0.01, 0.23, 0.99, 0.98);
      gPad->SetLeftMargin(0.1);      
      gPad->SetLogy();      
      zframe_data->SetTitle(""); 
      zframe_data->GetXaxis()->SetTitle(axis[k]); 
      zframe_data->GetXaxis()->SetLabelSize(0);
      zframe_data->GetYaxis()->SetTitle("Events"); 
      zframe_data->GetYaxis()->SetTitleOffset(1.30);
      zframe_data->SetMaximum(10000);
      zframe_data->SetMinimum(1);
      zframe_data->Draw();
      canvas->cd(2)->SetPad(0.01, 0.03, 0.99, 0.29);
      gPad->SetLeftMargin(0.1);
      gPad->SetBottomMargin(0.27);
      gPad->SetLogy(0);
      z_Pull_Frame->SetTitle(""); 
      z_Pull_Frame->GetYaxis()->SetTitle("Pull");
      z_Pull_Frame->GetYaxis()->SetTitleOffset(0.4);
      z_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
      z_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
      z_Pull_Frame->GetXaxis()->SetTitleSize(0.13);
      z_Pull_Frame->GetXaxis()->SetLabelSize(0.13);
      z_Pull_Frame->GetXaxis()->SetTitle(axis[k]); 
      z_Pull_Frame->Draw();
      TLine* line = new TLine(Min[k],0,Max[k],0);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw("same");
    }
    if(k==3){
      roohist_data_2D_ctau.plotOn(uframe_data);
      model_total2.plotOn(uframe_data,LineColor(kRed),LineWidth(1));
      RooHist *hpull = uframe_data->pullHist();
      RooPlot* u_Pull_Frame = u.frame();
      u_Pull_Frame->addPlotable(hpull,"P");
      model_total2.plotOn(uframe_data, LineColor(kRed), LineWidth(1));

      model_total2.plotOn(uframe_data,Components(model_signal_prompt_prompt),LineColor(kBlue),LineWidth(1)); 

      //model_total2.plotOn(uframe_data,Components(model_signal_prompt_nonprompt),LineColor(kGreen),LineWidth(2)); 

      //model_total2.plotOn(uframe_data,Components(model_signal_nonprompt_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(uframe_data,Components(model_signal_PmixNP),LineColor(kViolet-2),LineWidth(1)); 

      model_total2.plotOn(uframe_data,Components(model_signal_nonprompt_nonprompt),LineColor(kOrange),LineWidth(1)); 

      //model_total2.plotOn(uframe_data,Components(model_DCB1_comb2),LineColor(1),LineStyle(7),LineWidth(2)); 

      //model_total2.plotOn(uframe_data,Components(model_comb1_DCB2),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 
      model_total2.plotOn(uframe_data,Components(model_DCBmixComb),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(uframe_data,Components(model_comb1_comb2),LineColor(kOrange),LineStyle(7),LineWidth(1)); 
      canvas->cd(1)->SetPad(0.01, 0.23, 0.99, 0.98);
      gPad->SetLeftMargin(0.1);      
      gPad->SetLogy();      
      uframe_data->SetTitle(""); 
      uframe_data->GetXaxis()->SetTitle(axis[k]); 
      uframe_data->GetXaxis()->SetLabelSize(0);
      uframe_data->GetYaxis()->SetTitle("Events"); 
      uframe_data->GetYaxis()->SetTitleOffset(1.30);
      uframe_data->SetMaximum(10000);
      uframe_data->SetMinimum(1);
      uframe_data->Draw();
      canvas->cd(2)->SetPad(0.01, 0.03, 0.99, 0.29);
      gPad->SetLeftMargin(0.1);
      gPad->SetBottomMargin(0.27);
      gPad->SetLogy(0);
      u_Pull_Frame->SetTitle(""); 
      u_Pull_Frame->GetYaxis()->SetTitle("Pull");
      u_Pull_Frame->GetYaxis()->SetTitleOffset(0.4);
      u_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
      u_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
      u_Pull_Frame->GetXaxis()->SetTitleSize(0.13);
      u_Pull_Frame->GetXaxis()->SetLabelSize(0.13);
      u_Pull_Frame->GetXaxis()->SetTitle(axis[k]); 
      u_Pull_Frame->Draw();
      TLine* line = new TLine(Min[k],0,Max[k],0);
      line->SetLineColor(2);
      line->SetLineWidth(2);
      line->Draw("same");
    }
    //if(k==1){ 
    //  Double_t data_num = data_3D->Integral();
    //  CalculateCS(frac1, frac1_err, data_num); 
    //}
    //TH3D* model3D = model_xy.createHistogram("model3D",x,YVar(y));
    //model3D->Draw(); 
    
    canvas->cd(1);
    TPad *pad = new TPad("pad","pad",0.01,0.01,0.99,0.99);
    gPad->RedrawAxis();
    TString channelText = "";
    float channelTextFont   = 42;
    float channelTextSize   = 0.060;
    //TString cmsText     = "CMS";
    TString cmsText     = "";
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
    //lumi_13TeV = "36.3 fb^{-1}";
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
    //latex.DrawLatex(l+0.08, 1-t+lumiTextOffset*t-0.19, extraText);
    //latex.SetTextFont(channelTextFont);
    //latex.SetTextSize(channelTextSize);
    //latex.DrawLatex(l+0.165, 1-t+lumiTextOffset*t-0.19, channelText1);
    //latex.DrawLatex(l+0.165, 1-t+lumiTextOffset*t-0.28, channelText2);

    TLegend *pl2 = new TLegend(0.55,0.55,0.86,0.89);
    pl2->SetTextSize(0.028);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(data, "Observed",  "PE");
    pl2->AddEntry(h_total, "total fit model", "L");
    pl2->AddEntry(h_DCB_DCB_prompt_prompt, "J/#psi1 prompt+J/#psi2 prompt", "L");
    //pl2->AddEntry(h_DCB_DCB_prompt_nonprompt, "J/#psi1 prompt+J/#psi2 nonprompt", "L");
    //pl2->AddEntry(h_DCB_DCB_nonprompt_prompt, "J/#psi1 nonprompt+J/#psi2 prompt", "L");
    pl2->AddEntry(h_DCB_DCB_nonprompt_prompt, "J/#psi1(2) nonprompt+J/#psi2(1) prompt", "L");
    pl2->AddEntry(h_DCB_DCB_nonprompt_nonprompt, "J/#psi1 nonprompt+J/#psi2 nonprompt", "L");
    //pl2->AddEntry(h_DCB_comb, "J/#psi1+comb ", "L");
    //pl2->AddEntry(h_comb_DCB, "comb.+J/#psi2 ", "L");  
    pl2->AddEntry(h_comb_DCB, "J/#psi1(comb.)+comb.(J/#psi2) ", "L");  
    pl2->AddEntry(h_comb_comb, "comb.+comb. ", "L");  
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[k];
    cout<<"Finished "<<NAME<<endl;
    canvas->SaveAs("Fit_plots/data4D_fit_"+NAME+".pdf");
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

void GetHisto2D(char CUT[1000], TTree *Tree, TH2F* & histo, const char *plotX,const char *plotY, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX){
  char input[1000]; sprintf(input,"%s:%s>>h(%i,%f,%f,%i,%f,%f)",plotY,plotX,XBIN,XMIN,XMAX,YBIN,YMIN,YMAX);
  TH2F *temp = new TH2F("", "", XBIN,XMIN,XMAX,YBIN,YMIN,YMAX);
  if(Tree->Draw(input,CUT)){ Tree->Draw(input,CUT); TH2F* h=(TH2F*)gDirectory->Get("h"); histo = (TH2F*)h->Clone(); delete h; } else { histo = (TH2F*)temp->Clone(); }
  histo->SetDirectory(0);
  histo->SetName("histo");
  delete temp;
}
 
void GetHisto3D(char CUT[1000], TTree *Tree, TH3F* & histo, const char *plotX, const char *plotY, const char *plotZ, int XBIN, float XMIN, float XMAX, int YBIN, float YMIN, float YMAX, int ZBIN, float ZMIN, float ZMAX){
  char input[1000]; sprintf(input,"%s:%s:%s>>h(%i,%f,%f,%i,%f,%f,%i,%f,%f)",plotZ,plotY,plotX,XBIN,XMIN,XMAX,YBIN,YMIN,YMAX,ZBIN,ZMIN,ZMAX);
  TH3F *temp = new TH3F("", "", XBIN,XMIN,XMAX,YBIN,YMIN,YMAX,ZBIN,ZMIN,ZMAX);
  if(Tree->Draw(input,CUT)){ Tree->Draw(input,CUT); TH3F* h=(TH3F*)gDirectory->Get("h"); histo = (TH3F*)h->Clone(); delete h; } else { histo = (TH3F*)temp->Clone(); }
  histo->SetDirectory(0);
  histo->SetName("histo");
  delete temp;
}
