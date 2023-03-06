#include "plotterCRliftime_data_Jpsi.h"
//#include "CalculateCS.h"
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

void plotterCRliftime_data_Jpsi(){
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
 
  bool TriMatched = true;
  bool noMatch = false;
  bool fourMucut = false;
 
  TH1F *data;
  TH1F *data_Jpsi_mumu;
  TH1F *data_mumu_Jpsi;
  TH1F *data_mumu_mumu;
  TH1F *data_mumu;
  TH1F *data_Jpsi;


  float twomuon_cut = 0.05;
  float fourmuon_cut = 0.01;
  float Jpsi_mass_CR_do = 2.7;
  float Jpsi_mass_SR_do = 3.0;
  float Jpsi_mass_SR_up = 3.2;
  float Jpsi_mass_CR_up = 3.5;

  RooRealVar x("x", "x", 2.80, 3.40);
  RooRealVar y("y", "y", 2.80, 3.40);
  RooRealVar z("z", "z", 0, 5);
  RooRealVar u("u", "u", 0, 5);
  x.setBins(1000,"cache");

  //Jpsi+Jpsi

  RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi mean", 3.092); 
  RooRealVar Jpsi_width("Jpsi_width", "Jpsi width",0.03,0.02,0.04);
  RooRealVar alpha_1("alpha_1", "",1,0,10);
  RooRealVar alpha_2("alpha_2", "",1,0,10);
  RooRealVar n_1("n_1","",24,0,70);
  RooRealVar n_2("n_2","",10,0,70);
  
  RooDoubleCB DCB_x("DCB_x","DCB_x",x,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);
  RooDoubleCB DCB_y("DCB_y","DCB_y",y,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);
 
  //Construct combinatorial pdf
  RooRealVar a0("a0","a0",-0.365,-1,-0.001);
  RooRealVar a1("a1","a1",-0.0734,-1,1);

  RooChebychev comb_x("comb_x","comb_x",x,RooArgList(a0,a1));
  RooChebychev comb_y("comb_y","comb_y",y,RooArgList(a0,a1));

    //Step Function
  RooRealVar StepValue("StepValue", "StepValue", 1);
  RooRealVar LowerLimit("LowerLimit", "LowerLimit", 0);
  RooRealVar HigherLimit("HigherLimit", "HigherLimit", 10);
  RooAbsReal* Step_z = new RooStepFunction("Step_z", "Step_z", z, RooArgList(StepValue), RooArgList(LowerLimit, HigherLimit));
  RooAbsReal* Step_u = new RooStepFunction("Step_u", "Step_u", u, RooArgList(StepValue), RooArgList(LowerLimit, HigherLimit));

  //Non-prompt and prompt pdf
    //prompt fit
  /*
  RooRealVar SPS_mean("SPS_mean", "SPS mean", 0);
  RooRealVar SPS_width("SPS_width", "SPS width",0.0017,0.0015,0.0025);
  RooRealVar SPS_alpha_1("SPS_alpha_1", "",1.2,0.7,2);
  RooRealVar SPS_alpha_2("SPS_alpha_2", "",1.2,0.7,2);
  RooRealVar SPS_n_1("SPS_n_1","",5.8,2,10);
  RooRealVar SPS_n_2("SPS_n_2","",4.99,4,6);
  RooDoubleCB prompt_model("DCB","DCB",z,SPS_mean,SPS_width,SPS_alpha_1,SPS_n_1,SPS_alpha_2,SPS_n_2);
  */
   //SPS Exponential
     //Generic Exponential 
   RooRealVar SPS_Tau("SPS_Tau", "SPS_Tau", -1.133);      //1
   RooAbsPdf* Exp_SPS_z = new RooGenericPdf("SPS", "SPS", "Step_z * exp(z * SPS_Tau)", RooArgSet(*Step_z, z, SPS_Tau));
   RooAbsPdf* Exp_SPS_u = new RooGenericPdf("SPS", "SPS", "Step_u * exp(u * SPS_Tau)", RooArgSet(*Step_u, u, SPS_Tau));

     //Gauss
   RooRealVar SPS_Mean("SPS_Mean", "SPS_Mean",0.3127);
   RooRealVar SPS_Sigma("SPS_Sigma", "SPS_Sigma", 0.22215);
   RooAbsPdf* Gauss_SPS_z = new RooGaussian("Gauss_SPS", "Gauss_SPS", z, SPS_Mean, SPS_Sigma);
   RooAbsPdf* Gauss_SPS_u = new RooGaussian("Gauss_SPS", "Gauss_SPS", u, SPS_Mean, SPS_Sigma);
   
   //Convolution
   //RooFFTConvPdf prompt_model_z("ConvPdf_SPS_z", "ConvPdf_SPS_z", z, *Exp_SPS_z, *Gauss_SPS_z);
   //RooFFTConvPdf prompt_model_u("ConvPdf_SPS_u", "ConvPdf_SPS_u", u, *Exp_SPS_u, *Gauss_SPS_u);
   //RooNumConvPdf prompt_model_z("ConvPdf_SPS_z", "ConvPdf_SPS_z", z, *Exp_SPS_z, *Gauss_SPS_z);
   //RooNumConvPdf prompt_model_u("ConvPdf_SPS_u", "ConvPdf_SPS_u", u, *Exp_SPS_u, *Gauss_SPS_u);
    //Non-prompt fit function parameter
  /*
  RooRealVar BBbar_mean("BBbar_mean", "BBbar mean", 0.002323);
  RooRealVar BBbar_width("BBbar_width", "BBbar width",0.00233);
  RooRealVar BBbar_alpha_1("BBbar_alpha_1", "",1.72262);
  RooRealVar BBbar_alpha_2("BBbar_alpha_2", "",1.22250);
  RooRealVar BBbar_n_1("BBbar_n_1","",31.663, 1,100);
  RooRealVar BBbar_n_2("BBbar_n_2","",8.7886, 1,100);
  RooDoubleCB non_prompt_model("BBbar_DCB","BBbar_DCB",z,BBbar_mean,BBbar_width,BBbar_alpha_1,BBbar_n_1,BBbar_alpha_2,BBbar_n_2);
  RooRealVar gfrac2("gfrac2","",0.7,0.5,0.90);
  RooAddPdf model_ctau("model_ctau","prompt_model+non_prompt_model",RooArgList(prompt_model,non_prompt_model),RooArgList(gfrac2));
  */

    //Generic Exponential
  RooRealVar BBbar_Tau("BBbar_Tau", "BBbar_Tau", -0.8418);      //1

  RooAbsPdf* Exp_BBbar_z = new RooGenericPdf("BBbar_z", "BBbar_z", "Step_z * exp(z * BBbar_Tau)", RooArgSet(*Step_z, z, BBbar_Tau));
  RooAbsPdf* Exp_BBbar_u = new RooGenericPdf("BBbar_u", "BBbar_u", "Step_u * exp(u * BBbar_Tau)", RooArgSet(*Step_u, u, BBbar_Tau));

    //Gauss
  RooRealVar BBbar_Mean("BBbar_Mean", "BBbar_Mean",0.5374);
  RooRealVar BBbar_Sigma("BBbar_Sigma", "BBbar_Sigma", 0.40001);
  RooAbsPdf* Gauss_BBbar_z = new RooGaussian("Gauss_BBbar_z", "Gauss_BBbar_z", z, BBbar_Mean, BBbar_Sigma);
  RooAbsPdf* Gauss_BBbar_u = new RooGaussian("Gauss_BBbar_u", "Gauss_BBbar_u", u, BBbar_Mean, BBbar_Sigma);
    //Convolution
  //RooFFTConvPdf non_prompt_model_z("ConvPdf_BBbar_z", "ConvPdf_BBbar_z", z, *Exp_BBbar_z, *Gauss_BBbar_z);
  //RooFFTConvPdf non_prompt_model_u("ConvPdf_BBbar_u", "ConvPdf_BBbar_u", u, *Exp_BBbar_u, *Gauss_BBbar_u);

  RooRealVar gfrac2("gfrac2","",0.7, 0.4, 0.9);
  RooRealVar gfrac_SPS("gfrac_SPS","",0.0001);
  //RooAddPdf model_z("model_z","prompt_model_z+non_prompt_model_z",RooArgList(prompt_model_z,non_prompt_model_z),RooArgList(gfrac2));
  //RooAddPdf model_u("model_u","prompt_model_u+non_prompt_model_u",RooArgList(prompt_model_u,non_prompt_model_u),RooArgList(gfrac2));
  
  //construct composite pdf
  RooRealVar gfrac1("gfrac1","",0.6,0.5,0.85);
  
  RooAddPdf model_x("model_x","model_x",RooArgList(DCB_x,comb_x),RooArgList(gfrac1));
  RooAddPdf model_y("model_y","model_y",RooArgList(DCB_y,comb_y),RooArgList(gfrac1));
  //RooProdPdf model_SigLxyPV("","model_z*model_u",RooArgList(model_z,model_u));
  //RooProdPdf model_test("","prompt_model_z*prompt_model_u",RooArgList(prompt_model_z,prompt_model_u));


  vector<TString> name;  vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis; vector<string> variable;
  //name.push_back("Jpsi1_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{12}[GeV]");
  //name.push_back("Jpsi2_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{34}[GeV]");
  name.push_back("Jpsi1_cTau"); bin.push_back(100); Min.push_back(-0.02); Max.push_back(0.1); axis.push_back("ctau_{J/#psi1} [cm]"); variable.push_back("Jpsi1_cTau");
  //name.push_back("Jpsi2_cTau"); bin.push_back(30); Min.push_back(-0.02); Max.push_back(0.02); axis.push_back("ctau_{J/#psi2} [cm]"); variable.push_back("Jpsi2_cTau");
  //name.push_back("Jpsi1_LxyPVSig"); bin.push_back(50); Min.push_back(0); Max.push_back(5); axis.push_back("Significance of L_{xy}(J/#psi1) [cm]"); variable.push_back("Jpsi1_LxyPVSig");
  //name.push_back("Jpsi2_LxyPVSig"); bin.push_back(50); Min.push_back(0); Max.push_back(5); axis.push_back("Significance of L_{xy}(J/#psi2) [cm]"); variable.push_back("Jpsi2_LxyPVSig");

  char CUT_Jpsi_mumu[1000];
  char CUT_mumu_Jpsi[1000];
  char CUT_mumu_mumu[1000];
  char CUT_Jpsi[1000];
  char CUT_mumu[1000];
  char CUT_Jpsi1_SR[1000];
  char CUT_Jpsi1_CR[1000];
  char CUT_Jpsi2_SR[1000];
  char CUT_Jpsi2_CR[1000];

  sprintf(CUT_Jpsi1_SR, "(Jpsi1_mass>%f && Jpsi1_mass<%f)", Jpsi_mass_SR_do, Jpsi_mass_SR_up);
  sprintf(CUT_Jpsi1_CR, "((Jpsi1_mass>%f && Jpsi1_mass<%f)||(Jpsi1_mass>%f && Jpsi1_mass<%f))", Jpsi_mass_CR_do, Jpsi_mass_SR_do, Jpsi_mass_SR_up, Jpsi_mass_SR_up);

  sprintf(CUT_Jpsi2_SR, "(Jpsi2_mass>%f && Jpsi2_mass<%f)", Jpsi_mass_SR_do, Jpsi_mass_SR_up);
  sprintf(CUT_Jpsi2_CR, "((Jpsi2_mass>%f && Jpsi2_mass<%f)||(Jpsi2_mass>%f && Jpsi2_mass<%f))", Jpsi_mass_CR_do, Jpsi_mass_SR_do, Jpsi_mass_SR_up, Jpsi_mass_SR_up);


 for(int i=0; i<name.size(); i++){
  if(fourMucut){
   if((!TriMatched) && (!noMatch)){
    sprintf(CUT_Jpsi_mumu,"((!trigger_matched) && fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_SR, CUT_Jpsi2_CR);
    sprintf(CUT_mumu_Jpsi,"((!trigger_matched) && fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_CR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu_mumu,"((!trigger_matched) && fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_CR, CUT_Jpsi2_CR);

    if(variable[i].find("Jpsi1")!=string::npos){
      sprintf(CUT_Jpsi,"((!trigger_matched) && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi1_SR);
      sprintf(CUT_mumu,"((!trigger_matched) && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi1_CR);
    }

    if(variable[i].find("Jpsi2")!=string::npos){
      sprintf(CUT_Jpsi,"((!trigger_matched) && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi2_SR);
      sprintf(CUT_mumu,"((!trigger_matched) && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi2_CR);
    }
   }

   if(TriMatched && (!noMatch)){
    sprintf(CUT_Jpsi_mumu,"(trigger_matched && fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_SR, CUT_Jpsi2_CR);
    sprintf(CUT_mumu_Jpsi,"(trigger_matched && fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_CR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu_mumu,"(trigger_matched && fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_CR, CUT_Jpsi2_CR);

    if(variable[i].find("Jpsi1")!=string::npos){
      sprintf(CUT_Jpsi,"(trigger_matched && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi1_SR);
      sprintf(CUT_mumu,"(trigger_matched && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi1_CR);
    }
    if(variable[i].find("Jpsi2")!=string::npos){
      sprintf(CUT_Jpsi,"(trigger_matched && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi2_SR);
      sprintf(CUT_mumu,"(trigger_matched && fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi2_CR);
    }
   }

   if(noMatch){
    sprintf(CUT_Jpsi_mumu,"(fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_SR, CUT_Jpsi2_CR);
    sprintf(CUT_mumu_Jpsi,"(fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_CR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu_mumu,"(fourMu_VtxProb>%f && %s && %s)", fourmuon_cut, CUT_Jpsi1_CR, CUT_Jpsi2_CR);

    if(variable[i].find("Jpsi1")!=string::npos){
      sprintf(CUT_Jpsi,"(fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi1_SR);
      sprintf(CUT_mumu,"(fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi1_CR);
    }
    if(variable[i].find("Jpsi2")!=string::npos){
      sprintf(CUT_Jpsi,"(fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi2_SR);
      sprintf(CUT_mumu,"(fourMu_VtxProb>%f && %s )", fourmuon_cut, CUT_Jpsi2_CR);
    }
   }
  }  else {
   if((!TriMatched) && (!noMatch)){
    sprintf(CUT_Jpsi_mumu,"((!trigger_matched) && %s && %s)", CUT_Jpsi1_SR, CUT_Jpsi2_CR);
    sprintf(CUT_mumu_Jpsi,"((!trigger_matched) && %s && %s)", CUT_Jpsi1_CR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu_mumu,"((!trigger_matched) && %s && %s)", CUT_Jpsi1_CR, CUT_Jpsi2_CR);

    if(variable[i].find("Jpsi1")!=string::npos){
      sprintf(CUT_Jpsi,"((!trigger_matched) && %s )", CUT_Jpsi1_SR);
      sprintf(CUT_mumu,"((!trigger_matched) && %s )", CUT_Jpsi1_CR);
    }

    if(variable[i].find("Jpsi2")!=string::npos){
      sprintf(CUT_Jpsi,"((!trigger_matched) && %s )", CUT_Jpsi2_SR);
      sprintf(CUT_mumu,"((!trigger_matched) && %s )", CUT_Jpsi2_CR);
    }
   }

   if(TriMatched && (!noMatch)){
    sprintf(CUT_Jpsi_mumu,"(trigger_matched && %s && %s)", CUT_Jpsi1_SR, CUT_Jpsi2_CR);
    sprintf(CUT_mumu_Jpsi,"(trigger_matched && %s && %s)", CUT_Jpsi1_CR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu_mumu,"(trigger_matched && %s && %s)", CUT_Jpsi1_CR, CUT_Jpsi2_CR);

    if(variable[i].find("Jpsi1")!=string::npos){
      sprintf(CUT_Jpsi,"(trigger_matched && %s )", CUT_Jpsi1_SR);
      sprintf(CUT_mumu,"(trigger_matched && %s )", CUT_Jpsi1_CR);
    }
    if(variable[i].find("Jpsi2")!=string::npos){
      sprintf(CUT_Jpsi,"(trigger_matched && %s )", CUT_Jpsi2_SR);
      sprintf(CUT_mumu,"(trigger_matched && %s )", CUT_Jpsi2_CR);
    }
   }

   if(noMatch){
    sprintf(CUT_Jpsi_mumu,"(%s && %s)", CUT_Jpsi1_SR, CUT_Jpsi2_CR);
    sprintf(CUT_mumu_Jpsi,"(%s && %s)", CUT_Jpsi1_CR, CUT_Jpsi2_SR);
    sprintf(CUT_mumu_mumu,"(%s && %s)", CUT_Jpsi1_CR, CUT_Jpsi2_CR);

    if(variable[i].find("Jpsi1")!=string::npos){
      sprintf(CUT_Jpsi,"( %s )", CUT_Jpsi1_SR);
      sprintf(CUT_mumu,"( %s )", CUT_Jpsi1_CR);
    }
    if(variable[i].find("Jpsi2")!=string::npos){
      sprintf(CUT_Jpsi,"( %s )", CUT_Jpsi2_SR);
      sprintf(CUT_mumu,"( %s )", CUT_Jpsi2_CR);
    }
   }
  }
        const char *plotX = name[i];
	int   XBIN = bin[i];
	float XMIN = Min[i];
	float XMAX = Max[i];
	GetHisto(CUT_Jpsi_mumu, Tree, data ,plotX,bin[i],Min[i],Max[i]);
	GetHisto(CUT_Jpsi_mumu, Tree, data_Jpsi_mumu ,plotX,bin[i],Min[i],Max[i]);
	GetHisto(CUT_mumu_Jpsi, Tree, data_mumu_Jpsi ,plotX,bin[i],Min[i],Max[i]);
	GetHisto(CUT_mumu_mumu, Tree, data_mumu_mumu ,plotX,bin[i],Min[i],Max[i]);

	GetHisto(CUT_Jpsi, Tree, data_Jpsi ,plotX,bin[i],Min[i],Max[i]);
	GetHisto(CUT_mumu, Tree, data_mumu ,plotX,bin[i],Min[i],Max[i]);
        cout<<data_Jpsi_mumu ->Integral()<<endl;
        cout<<data_mumu_Jpsi ->Integral()<<endl;
        cout<<data_mumu_mumu ->Integral()<<endl;

        TString NAME = name[i];

        if (fourMucut) {
          if ((!TriMatched) &&(!noMatch)) NAME = NAME + "_triNotMatched";
          if (TriMatched &&(!noMatch)) NAME = NAME + "_triMatched";
          if (noMatch) NAME = NAME + "_noMatch";
        } else {
          if ((!TriMatched) &&(!noMatch)) NAME = NAME + "_No4Mucut_triNotMatched";
          if (TriMatched &&(!noMatch)) NAME = NAME + "_No4Mucut_triMatched";
          if (noMatch) NAME = NAME + "_No4Mucut_noMatch";
        }
        /*
        data_Jpsi_mumu->SetTitle("");
        data_Jpsi_mumu->GetXaxis()->SetTitle(axis[i]);
        data_Jpsi_mumu->GetYaxis()->SetTitle("Events");
        data_Jpsi_mumu->GetYaxis()->SetTitleOffset(1.30);
        data_Jpsi_mumu->SetMarkerColor(1); 
        data_Jpsi_mumu->SetMarkerStyle(20); 
        data_Jpsi_mumu->SetMarkerSize(1.3);
        data_Jpsi_mumu->Draw("PE"); 
        canvas->SaveAs(NAME+"_Jpsi_mumu.pdf");

        data_mumu_Jpsi->SetTitle("");
        data_mumu_Jpsi->GetXaxis()->SetTitle(axis[i]);
        data_mumu_Jpsi->GetYaxis()->SetTitle("Events");
        data_mumu_Jpsi->GetYaxis()->SetTitleOffset(1.30);
        data_mumu_Jpsi->Draw();
        canvas->SaveAs(NAME+"_mumu_Jpsi.pdf");

        data_mumu_mumu->SetTitle("");
        data_mumu_mumu->GetXaxis()->SetTitle(axis[i]);
        data_mumu_mumu->GetYaxis()->SetTitle("Events");
        data_mumu_mumu->GetYaxis()->SetTitleOffset(1.30);
        data_mumu_mumu->Draw();
        canvas->SaveAs(NAME+"_mumu_mumu.pdf");
        */
        data_Jpsi->SetTitle("");
        data_Jpsi->GetXaxis()->SetTitle(axis[i]);
        data_Jpsi->GetYaxis()->SetTitle("Events");
        data_Jpsi->GetYaxis()->SetTitleOffset(1.30);
        data_Jpsi->Draw(); 
        canvas->SaveAs(NAME+"_Jpsi.pdf");
        /*
        data_mumu->SetTitle("");
        data_mumu->GetXaxis()->SetTitle(axis[i]);
        data_mumu->GetYaxis()->SetTitle("Events");
        data_mumu->GetYaxis()->SetTitleOffset(1.30);
        data_mumu->Draw(); 
        canvas->SaveAs(NAME+"_mumu.pdf");
        */
        cout<<"Finished "<<NAME<<endl;
  
      //data->Draw("LEGO");
      //canvas->SaveAs("data.pdf") ;

  TH1F *h_DCB_comb; h_DCB_comb = (TH1F*)data->Clone();  h_DCB_comb->SetLineColor(kGreen); h_DCB_comb->SetLineWidth(1); h_DCB_comb->SetLineStyle(7);
  TH1F *h_comb_DCB; h_comb_DCB = (TH1F*)data->Clone(); h_comb_DCB->SetLineColor(kViolet-2); h_comb_DCB->SetLineWidth(1); h_comb_DCB->SetLineStyle(7);
  TH1F *h_comb_comb; h_comb_comb = (TH1F*)data->Clone(); h_comb_comb->SetLineColor(kOrange); h_comb_comb->SetLineWidth(1); h_comb_comb->SetLineStyle(7);
  //TH1F *h_prompt; h_prompt = (TH1F*)data->Clone(); h_prompt->SetLineColor(1); h_prompt->SetLineStyle(1);
  //TH1F *h_nonprompt; h_nonprompt = (TH1F*)data->Clone(); h_nonprompt->SetLineColor(1); h_nonprompt->SetLineStyle(7);
  TH1F *h_total; h_total = (TH1F*)data->Clone(); h_total->SetLineColor(kRed); h_total->SetLineWidth(1);


  //data_2D_SigLxyPV->Draw("LEGO");
  //data_2D_SigLxyPV->Draw("COLZ");
  //canvas->SaveAs("data2D_SigLxyPV.pdf");
  
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
