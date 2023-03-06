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

void RooFit3D_data_Jpsi(){
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
  TH3F *data_3D; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_x; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_y; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);
  TH1F *data_z; //= new TH1F("DPS",  "DPS",  30,2.8,3.4);

  float twomuon_cut = 0.05;
  float fourmuon_cut = 0.01;
  float Jpsi_mass_do = 2.8;
  float Jpsi_mass_up = 3.4;

  RooRealVar x("x", "x", 2.80, 3.40);
  RooRealVar y("y", "y", 2.80, 3.40);
  RooRealVar z("z", "z", 0, 5);
  x.setBins(1000,"cache");

  //Construct Jpsi Jpsi pdf
  RooRealVar PDGmZ("PDGmZ","PDGmZ",3.092);
  RooRealVar PDGwZ("PDGwZ","PDGwZ",0.05,0.01,0.08);
  RooBreitWigner PDGBW("PDGBW","PDGBW",x,PDGmZ,PDGwZ);
  PDGmZ.setConstant(kTRUE);
  PDGwZ.setConstant(kTRUE);
  //Jpsi+Jpsi

  RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi mean", 3.092); 
  RooRealVar Jpsi_width("Jpsi_width", "Jpsi width",0.03,0.02,0.04);
  RooRealVar alpha_1("alpha_1", "",1,0,10);
  RooRealVar alpha_2("alpha_2", "",1,0,10);
  RooRealVar n_1("n_1","",24,0,70);
  RooRealVar n_2("n_2","",10,0,70);

  
  RooDoubleCB DCB_x("DCB_x","DCB_x",x,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);
  RooDoubleCB DCB_y("DCB_y","DCB_y",y,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);

  //RooCBShape CB_1("CB_1","CB_1",x,Jpsi_mean,Jpsi_width,alpha_1,n_1);
  //RooCBShape CB_2("CB_2","CB_2",x,Jpsi_mean,Jpsi_width,alpha_2,n_2);

  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,DCB);
  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,CB_1);
 
 
  //Construct combinatorial pdf
  RooRealVar a0("a0","a0",-0.365,-1,-0.001);
  RooRealVar a1("a1","a1",-0.0734,-1,1);

  RooChebychev comb_x("comb_x","comb_x",x,RooArgList(a0,a1));
  RooChebychev comb_y("comb_y","comb_y",y,RooArgList(a0,a1));

    //Step Function
  RooRealVar x_StepValue("x_StepValue", "x_StepValue", 1);
  RooRealVar x_LowerLimit("x_LowerLimit", "x_LowerLimit", 0);
  RooRealVar x_HigherLimit("x_HigherLimit", "x_HigherLimit", 1000000);
  RooAbsReal* Stepx = new RooStepFunction("Stepx", "Stepx", z, RooArgList(x_StepValue), RooArgList(x_LowerLimit, x_HigherLimit));

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
   RooAbsPdf* Expx_SPS = new RooGenericPdf("SPS", "SPS", "Stepx * exp(z * SPS_Tau)", RooArgSet(*Stepx, z, SPS_Tau));
   //Gauss
   RooRealVar SPS_Mean("SPS_Mean", "SPS_Mean",0.3127);
   RooRealVar SPS_Sigma("SPS_Sigma", "SPS_Sigma", 0.22215);
   RooAbsPdf* Gauss_SPS = new RooGaussian("Gauss_SPS", "Gauss_SPS", z, SPS_Mean, SPS_Sigma);
   
   //Convolution
   RooFFTConvPdf prompt_model("ConvPdf_SPS", "ConvPdf_SPS", z, *Expx_SPS, *Gauss_SPS);
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

  RooAbsPdf* Expx_BBbar = new RooGenericPdf("BBbar", "BBbar", "Stepx * exp(z * BBbar_Tau)", RooArgSet(*Stepx, z, BBbar_Tau));

    //Gauss
  RooRealVar BBbar_Mean("BBbar_Mean", "BBbar_Mean",0.5374);
  RooRealVar BBbar_Sigma("BBbar_Sigma", "BBbar_Sigma", 0.40001);
  RooAbsPdf* Gauss_BBbar = new RooGaussian("Gauss_BBbar", "Gauss_BBbar", z, BBbar_Mean, BBbar_Sigma);
    //Convolution
  RooFFTConvPdf non_prompt_model("ConvPdf_BBbar", "ConvPdf_BBbar", z, *Expx_BBbar, *Gauss_BBbar);

  RooRealVar gfrac2("gfrac2","",0.5,0.2,0.85);
  RooAddPdf model_ctau("model_ctau","prompt_model+non_prompt_model",RooArgList(prompt_model,non_prompt_model),RooArgList(gfrac2));
  
  //construct composite pdf
  RooRealVar gfrac1("gfrac1","",0.6,0.5,0.85);
  
  RooAddPdf model_x("model_x","model_x",RooArgList(DCB_x,comb_x),RooArgList(gfrac1));
  RooAddPdf model_y("model_y","model_y",RooArgList(DCB_y,comb_y),RooArgList(gfrac1));
  RooProdPdf model_total1("model_total1","model_x*model_y*model_z",RooArgList(model_x,model_y,model_ctau));
  
  //same pdf function but different component
  RooProdPdf model_signal_prompt("model_signal_prompt","DCB_x*DCB_y*prompt_z",RooArgList(DCB_x,DCB_y,prompt_model));
  RooProdPdf model_signal_nonprompt("model_signal_nonprompt","DCB_x*DCB_y*nonprompt_z",RooArgList(DCB_x,DCB_y,non_prompt_model));

  RooProdPdf model_DCB1_comb2_prompt("model_bkg1_prompt","DCB_x*comb_y*prompt_z",RooArgList(DCB_x,comb_y,prompt_model));
  RooProdPdf model_DCB1_comb2_nonprompt("model_bkg1_nonprompt","DCB_x*comb_y*nonprompt_z",RooArgList(DCB_x,comb_y,non_prompt_model));

  RooProdPdf model_comb1_DCB2_prompt("model_bkg2_prompt_nonprompt","comb_x*DCB_y*prompt_z",RooArgList(comb_x,DCB_y,prompt_model));
  RooProdPdf model_comb1_DCB2_nonprompt("model_bkg2_nonprompt","comb_x*DCB_y*nonprompt_z",RooArgList(comb_x,DCB_y,non_prompt_model));

  RooProdPdf model_comb1_comb2_prompt("model_bkg3_prompt","comb_x*comb_y*prompt_z",RooArgList(comb_x,comb_y,prompt_model));
  RooProdPdf model_comb1_comb2_nonprompt("model_bkg3_nonprompt","comb_x*comb_y*nonprompt_z",RooArgList(comb_x,comb_y,non_prompt_model));




  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  name.push_back("fourMuFit_ups1_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{12}[GeV]");
  name.push_back("fourMuFit_ups2_mass"); bin.push_back(60); Min.push_back(2.8); Max.push_back(3.4); axis.push_back("m_{34}[GeV]");
  //name.push_back("fourMuFit_ups1_cTau_MC"); bin.push_back(200); Min.push_back(-0.02); Max.push_back(0.2); axis.push_back("J/#psi1 c#tau [cm]");
  //name.push_back("fourMuFit_ups1_cTau_MC"); bin.push_back(200); Min.push_back(-0.02); Max.push_back(0.2); axis.push_back("J/#psi1 c#tau [cm]");
  //name.push_back("(abs(fourMuFit_ups1_LxyPVSig_MC+fourMuFit_ups2_LxyPVSig_MC)+abs(fourMuFit_ups1_LxyPVSig_MC-fourMuFit_ups2_LxyPVSig_MC))/2"); bin.push_back(100); Min.push_back(0); Max.push_back(5); axis.push_back("Significance of L_{xy}(J/#psi) [cm]");
  //name.push_back("abs( fourMuFit_ups1_LxyPVSig_MC + fourMuFit_ups2_LxyPVSig_MC )"); bin.push_back(100); Min.push_back(0); Max.push_back(5); axis.push_back("Significance of L_{xy}(J/#psi) [cm]");
  name.push_back("abs(fourMuFit_ups1_LxyPVSig_MC)"); bin.push_back(100); Min.push_back(0); Max.push_back(5); axis.push_back("Significance of L_{xy}(J/#psi) [cm]");


  for(int i=0; i<name.size(); i++){
    for(int j=i+1; j<name.size(); j++){
      for(int k=j+1; k<name.size(); k++){
        const char *plotX = name[i];
	const char *plotY = name[j];
        const char *plotZ = name[k];
	int XBIN = bin[i];   int YBIN = bin[j];     int ZBIN = bin[k];
	float XMIN = Min[i]; float YMIN = Min[j]; float ZMIN = Min[k];
	float XMAX = Max[i]; float YMAX = Max[j]; float ZMAX = Max[k];
	char CUTpre[1000];
	char CUT[1000];
	//sprintf(CUTpre,"(((Jpsi1_pt>6 && abs(Jpsi1_eta)>1)||(Jpsi1_pt>7 && abs(Jpsi1_eta)<1)) && ((Jpsi2_pt>6 && abs(Jpsi2_eta)>1)||(Jpsi2_pt>7 && abs(Jpsi2_eta)<1)))");
	sprintf(CUTpre,"fourMuFit_VtxProb>%f && fourMuFit_ups1_mass>%f && fourMuFit_ups1_mass<%f && fourMuFit_ups2_mass>%f && fourMuFit_ups2_mass<%f", fourmuon_cut, Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
	//sprintf(CUTpre,"fourMuFit_ups1_mass>%f && fourMuFit_ups1_mass<%f && fourMuFit_ups2_mass>%f && fourMuFit_ups2_mass<%f", Jpsi_mass_do, Jpsi_mass_up, Jpsi_mass_do, Jpsi_mass_up );
	sprintf(CUT,"1 *%s ", CUTpre);
	//sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12)) * *%s ", CUTpre);
	GetHisto(CUT, Tree, data_x ,plotX,bin[i],Min[i],Max[i]);
	GetHisto(CUT, Tree, data_y ,plotY,bin[j],Min[j],Max[j]);
        GetHisto(CUT, Tree, data_z ,plotZ,bin[k],Min[k],Max[k]);
	GetHisto(CUT, Tree, data ,plotX,bin[i],Min[i],Max[i]);

        GetHisto3D(CUT, Tree, data_3D ,plotX,plotY,plotZ,XBIN,XMIN,XMAX,YBIN,YMIN,YMAX,ZBIN,ZMIN,ZMAX);
        cout<<data_x->Integral()<<endl;
        cout<<data_y->Integral()<<endl;
        cout<<data_z->Integral()<<endl;
        cout<<data_3D->Integral()<<endl;

	for(int n=1; n<data_3D->GetNbinsX()+1; n++){
          for(int m=1; m<data_3D->GetNbinsY()+1; m++){
            for(int l=1; l<data_3D->GetNbinsZ()+1; l++){
              data_3D->SetBinContent(n,m,l,data_3D->GetBinContent(n,m,l));
              data_3D->SetBinError(n,m,l,sqrt(data_3D->GetBinContent(n,m,l)));

            }
          }
        }
      }
    }
  }
      //data->Draw("LEGO");
      //canvas->SaveAs("data.pdf") ;
  data->SetMarkerColor(1); data->SetMarkerStyle(20); data->SetMarkerSize(1.3);
  TH1F *h_DCB_DCB_prompt; h_DCB_DCB_prompt = (TH1F*)data->Clone(); h_DCB_DCB_prompt->SetLineColor(kBlue); h_DCB_DCB_prompt->SetLineWidth(1); 
  TH1F *h_DCB_comb_prompt; h_DCB_comb_prompt = (TH1F*)data->Clone(); h_DCB_comb_prompt->SetMarkerSize(0); h_DCB_comb_prompt->SetLineColor(kGreen); h_DCB_comb_prompt->SetLineWidth(1);
  TH1F *h_comb_DCB_prompt; h_comb_DCB_prompt = (TH1F*)data->Clone(); h_comb_DCB_prompt->SetLineColor(kViolet-2); h_comb_DCB_prompt->SetLineWidth(1);
  TH1F *h_comb_comb_prompt; h_comb_comb_prompt = (TH1F*)data->Clone(); h_comb_comb_prompt->SetLineColor(kOrange); h_comb_comb_prompt->SetLineWidth(1);
  TH1F *h_DCB_DCB_nonprompt; h_DCB_DCB_nonprompt = (TH1F*)data->Clone(); h_DCB_DCB_nonprompt->SetLineColor(kBlue); h_DCB_DCB_nonprompt->SetLineWidth(1); h_DCB_DCB_nonprompt->SetLineStyle(7); 
  TH1F *h_DCB_comb_nonprompt; h_DCB_comb_nonprompt = (TH1F*)data->Clone(); h_DCB_comb_nonprompt->SetMarkerSize(0); h_DCB_comb_nonprompt->SetLineColor(kGreen); h_DCB_comb_nonprompt->SetLineWidth(1); h_DCB_comb_nonprompt->SetLineStyle(7);
  TH1F *h_comb_DCB_nonprompt; h_comb_DCB_nonprompt = (TH1F*)data->Clone(); h_comb_DCB_nonprompt->SetLineColor(kViolet-2); h_comb_DCB_nonprompt->SetLineWidth(1); h_comb_DCB_nonprompt->SetLineStyle(7);
  TH1F *h_comb_comb_nonprompt; h_comb_comb_nonprompt = (TH1F*)data->Clone(); h_comb_comb_nonprompt->SetLineColor(kOrange); h_comb_comb_nonprompt->SetLineWidth(1); h_comb_comb_nonprompt->SetLineStyle(7);
  //TH1F *h_prompt; h_prompt = (TH1F*)data->Clone(); h_prompt->SetLineColor(1); h_prompt->SetLineStyle(1);
  //TH1F *h_nonprompt; h_nonprompt = (TH1F*)data->Clone(); h_nonprompt->SetLineColor(1); h_nonprompt->SetLineStyle(7);
  TH1F *h_total; h_total = (TH1F*)data->Clone(); h_total->SetLineColor(kRed); h_total->SetLineWidth(1);


  RooDataHist roohist_data_3D("roohist_data","roohist_data",RooArgSet(x,y,z),data_3D);
  RooPlot* xframe_data = x.frame();
  RooPlot* yframe_data = y.frame();
  RooPlot* zframe_data = z.frame();

  //3D fit
  model_total1.fitTo(roohist_data_3D,SumW2Error(kTRUE));
  Double_t frac1 = gfrac1.getVal();
  Double_t frac2 = gfrac2.getVal();
  Double_t frac1_err = gfrac1.getError(); 
  Double_t frac2_err = gfrac2.getError(); 

  //parameter (g1+(1-g1))*(g1+(1-g1))*(g2+(1-g2))
  Double_t frac111 = frac1*frac1*frac2;
  Double_t frac112 = frac1*frac1*(1-frac2); 
  Double_t frac121 = frac1*(1-frac1)*frac2; 
  Double_t frac122 = frac1*(1-frac1)*(1-frac2); 

  Double_t frac211 = (1-frac1)*frac1*frac2; 
  Double_t frac212 = (1-frac1)*frac1*(1-frac2); 
  Double_t frac221 = (1-frac1)*(1-frac1)*frac2; 
  Double_t frac222 = (1-frac1)*(1-frac1)*(1-frac2);

  RooRealVar gfrac111("gfrac111","",frac111);  
  RooRealVar gfrac112("gfrac112","",frac112);  
  RooRealVar gfrac121("gfrac121","",frac121);  
  RooRealVar gfrac122("gfrac122","",frac122);  

  RooRealVar gfrac211("gfrac211","",frac211);  
  RooRealVar gfrac212("gfrac212","",frac212);  
  RooRealVar gfrac221("gfrac221","",frac221);  
  RooRealVar gfrac222("gfrac222","",frac222);  

  RooAddPdf model_total2("model_total2","model_total2", RooArgList(model_signal_prompt, model_signal_nonprompt,
                                                       model_DCB1_comb2_prompt, model_DCB1_comb2_nonprompt,
                                                       model_comb1_DCB2_prompt, model_comb1_DCB2_nonprompt,
                                                       model_comb1_comb2_prompt,model_comb1_comb2_nonprompt),
   RooArgList(gfrac111, gfrac112, gfrac121, gfrac122, gfrac211, gfrac212, gfrac221));
  
  for(int k=0; k<name.size(); k++){
    if(k==0){
      roohist_data_3D.plotOn(xframe_data);
      model_total2.plotOn(xframe_data,LineColor(kRed),LineWidth(1)); 
 
      model_total2.plotOn(xframe_data,Components(model_signal_prompt),LineColor(kBlue),LineWidth(1)); 
      model_total2.plotOn(xframe_data,Components(model_signal_nonprompt),LineColor(kBlue),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(xframe_data,Components(model_DCB1_comb2_prompt),LineColor(kGreen),LineWidth(1)); 
      model_total2.plotOn(xframe_data,Components(model_DCB1_comb2_nonprompt),LineColor(kGreen),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(xframe_data,Components(model_comb1_DCB2_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(xframe_data,Components(model_comb1_DCB2_nonprompt),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(xframe_data,Components(model_comb1_comb2_prompt),LineColor(kOrange+9),LineWidth(1)); 
      model_total2.plotOn(xframe_data,Components(model_comb1_comb2_nonprompt),LineColor(kOrange+9),LineStyle(7),LineWidth(1)); 

      xframe_data->SetTitle(""); 
      xframe_data->GetXaxis()->SetTitle(axis[k]); 
      xframe_data->GetYaxis()->SetTitle("Events"); 
      xframe_data->GetYaxis()->SetTitleOffset(1.30);
      xframe_data->Draw();
      } 

    if(k==1){
      roohist_data_3D.plotOn(yframe_data);
      model_total2.plotOn(yframe_data,LineColor(kRed),LineWidth(1)); 
 
      model_total2.plotOn(yframe_data,Components(model_signal_prompt),LineColor(kBlue),LineWidth(1)); 
      model_total2.plotOn(yframe_data,Components(model_signal_nonprompt),LineColor(kBlue),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(yframe_data,Components(model_DCB1_comb2_prompt),LineColor(kGreen),LineWidth(1)); 
      model_total2.plotOn(yframe_data,Components(model_DCB1_comb2_nonprompt),LineColor(kGreen),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(yframe_data,Components(model_comb1_DCB2_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(yframe_data,Components(model_comb1_DCB2_nonprompt),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(yframe_data,Components(model_comb1_comb2_prompt),LineColor(kOrange+9),LineWidth(1)); 
      model_total2.plotOn(yframe_data,Components(model_comb1_comb2_nonprompt),LineColor(kOrange+9),LineStyle(7),LineWidth(1)); 

      yframe_data->SetTitle(""); 
      yframe_data->GetXaxis()->SetTitle(axis[k]); 
      yframe_data->GetYaxis()->SetTitle("Events"); 
      yframe_data->GetYaxis()->SetTitleOffset(1.30);
      yframe_data->Draw();
    }
    
    if(k==2){
      roohist_data_3D.plotOn(zframe_data);
      model_total2.plotOn(zframe_data,LineColor(kRed),LineWidth(1)); 
 
      model_total2.plotOn(zframe_data,Components(model_signal_prompt),LineColor(kBlue),LineWidth(1)); 
      model_total2.plotOn(zframe_data,Components(model_signal_nonprompt),LineColor(kBlue),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(zframe_data,Components(model_DCB1_comb2_prompt),LineColor(kGreen),LineWidth(1)); 
      model_total2.plotOn(zframe_data,Components(model_DCB1_comb2_nonprompt),LineColor(kGreen),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(zframe_data,Components(model_comb1_DCB2_prompt),LineColor(kViolet-2),LineWidth(1)); 
      model_total2.plotOn(zframe_data,Components(model_comb1_DCB2_nonprompt),LineColor(kViolet-2),LineStyle(7),LineWidth(1)); 

      model_total2.plotOn(zframe_data,Components(model_comb1_comb2_prompt),LineColor(kOrange+9),LineWidth(1)); 
      model_total2.plotOn(zframe_data,Components(model_comb1_comb2_nonprompt),LineColor(kOrange+9),LineStyle(7),LineWidth(1)); 

      zframe_data->SetTitle(""); 
      zframe_data->GetXaxis()->SetTitle(axis[k]); 
      zframe_data->GetYaxis()->SetTitle("Events"); 
      zframe_data->GetYaxis()->SetTitleOffset(1.30);
      zframe_data->Draw();
    }
    if(k==1){ 
      Double_t data_num = data_3D->Integral();
      CalculateCS(frac1, frac1_err, data_num); 
    }
    
    //TH3D* model3D = model_xy.createHistogram("model3D",x,YVar(y));
    //model3D->Draw(); 

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

    TLegend *pl2 = new TLegend(0.55,0.55,0.85,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(data, "Observed",  "PE");
    pl2->AddEntry(h_total, "total fit model", "L");
    pl2->AddEntry(h_DCB_DCB_prompt, "J/#psi1+J/#psi2 prompt", "L");
    pl2->AddEntry(h_DCB_comb_prompt, "J/#psi1+comb prompt", "L");
    pl2->AddEntry(h_comb_DCB_prompt, "comb.+J/#psi2 prompt", "L");  
    pl2->AddEntry(h_comb_comb_prompt, "comb.+comb. prompt", "L");  
    pl2->AddEntry(h_DCB_DCB_nonprompt, "J/#psi1+J/#psi2 non-prompt", "L");
    pl2->AddEntry(h_DCB_comb_nonprompt, "J/#psi1+comb non-prompt", "L");
    pl2->AddEntry(h_comb_DCB_nonprompt, "comb.+J/#psi2 non-prompt", "L");  
    pl2->AddEntry(h_comb_comb_nonprompt, "comb.+comb. non-prompt", "L");  
    //pl2->AddEntry(h_prompt, "Prompt", "L");  
    //pl2->AddEntry(h_nonprompt, "Non-prompt", "L");  
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[k];
    cout<<"Finished "<<NAME<<endl;
    canvas->SaveAs("data3D_fit_"+NAME+".pdf");
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
  char input[50]; sprintf(input,"%s:%s>>h(%i,%f,%f,%i,%f,%f)",plotY,plotX,XBIN,XMIN,XMAX,YBIN,YMIN,YMAX);
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
