#include "RooFit_DPS_Jpsi.h"
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

void RooFit_DPS_Jpsi(){
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
  RooRealVar x("x", "x", -0.2 , 0.2);
  x.setBins(1000,"cache");
  //RooRealVar PDGmZ("PDGmZ","PDGmZ",3.092);
  //RooRealVar PDGwZ("PDGwZ","PDGwZ",0.05,0.01,0.08);
  //RooBreitWigner PDGBW("PDGBW","PDGBW",x,PDGmZ,PDGwZ);
  //PDGmZ.setConstant(kTRUE);
  //PDGwZ.setConstant(kTRUE);
  //RooRealVar mean("mean", "mean",0.00065,-0.002,0.002); 
  //RooRealVar sigma("sigma", "sigma",0.0037,0,0.01);
  //RooRealVar rlife("rlife", "rlife",0.016,0,0.5);
  RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi mean", 0.005,0,0.02); 
  RooRealVar Jpsi_width("Jpsi_width", "Jpsi width",0.01,0,0.06);
  RooRealVar alpha_1("alpha_1", "",1,0,10);
  RooRealVar alpha_2("alpha_2", "",1,0,10);
  RooRealVar n_1("n_1","",10,0,100);
  RooRealVar n_2("n_2","",10,0,500);
  //RooRealVar g1frac("g1frac","",0.5);
  //RooRealVar g2frac("g2frac","",0,1.0);
  RooDoubleCB DCB("DCB","DCB",x,Jpsi_mean,Jpsi_width,alpha_1,n_1,alpha_2,n_2);
  //RooCBShape CB_1("CB_1","CB_1",x,Jpsi_mean,Jpsi_width,alpha_1,n_1);
  //RooCBShape CB_2("CB_2","CB_2",x,Jpsi_mean,Jpsi_width,alpha_2,n_2);
  //RooAddPdf  DCB("DCB","CB_1+CB_2",RooArgList(CB_1,CB_2),RooArgList(g1frac));
  //RooAddPdf  BW_DCB("BW_DCB","BW+CB_1+CB_2",RooArgList(PDGBW,DCB),RooArgList(g2frac));

  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,DCB);
  //RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",x,PDGBW,CB_1);
  //RooGExpModel CtauModel("ctauModel","ctauModel",x,mean,sigma,rlife);

  vector<TString> name;         vector<int> bin;   vector<float> Min;  vector<float> Max;  vector<TString> axis;
  //name.push_back("MuMuFit_ctauPV_MC"); bin.push_back(50); Min.push_back(-0.01); Max.push_back(0.02); axis.push_back("J/#psi c#tau [cm]");
  name.push_back("MuMuFit_lxyPV_MC"); bin.push_back(50); Min.push_back(-0.05); Max.push_back(0.08); axis.push_back("L_{xy}(J/#psi) [cm]");


  for(int i=0; i<name.size(); i++){
    const char *plot = name[i];
    char CUT[1000];
    //sprintf(CUT,"(1./((w_acc_Jpsi1*w_acc_Jpsi2) * (w_reco_Jpsi1*w_reco_Jpsi2) * (w_eff_Jpsi1*w_eff_Jpsi2) *w_trig_Jpsi12))");
    //sprintf(CUT,"");
    sprintf(CUT,"fourMuFit_VtxProb>0.01");
    GetHisto(CUT, Tree, DPS ,plot,bin[i],Min[i],Max[i]);
       
    cout<<DPS->Integral()<<endl;; 
   
    for(int j=1; j<DPS->GetNbinsX()+1; j++){
      //data->SetBinContent(j,data_pre->GetBinContent(j));

      DPS->SetBinContent(j,DPS->GetBinContent(j));
      DPS->SetBinError(j,sqrt(DPS->GetBinContent(j)));
    }

    TH1F *h_DCB; h_DCB = (TH1F*)DPS->Clone(); h_DCB->SetLineColor(kBlue);

    RooDataHist roohist_DPS("roohist_DPS","roohist_DPS",x,DPS);
    RooPlot* frame_DPS = x.frame();
    DCB.fitTo(roohist_DPS,SumW2Error(kTRUE));
    roohist_DPS.plotOn(frame_DPS);
    DCB.plotOn(frame_DPS,LineColor(kBlue)); 
    frame_DPS->SetTitle(""); 
    frame_DPS->GetXaxis()->SetTitle(axis[i]); 
    frame_DPS->GetYaxis()->SetTitle("Events"); 
    frame_DPS->SetMaximum(800);;
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
    lumi_13TeV = "36.3 fb^{-1}";
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




    TLegend *pl2 = new TLegend(0.63,0.70,0.88,0.87);
    pl2->SetTextSize(0.030);
    pl2->SetTextFont(62);
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(&roohist_DPS, "BBar to Jpsi Jpsi",  "PE");
    pl2->AddEntry(h_DCB, "fit", "L");
    pl2->SetBorderSize(0);
    pl2->Draw();

    TString NAME = name[i];
    cout<<"Finished "<<NAME<<endl;

    canvas->SaveAs("BBar_fit_"+NAME+".pdf");
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
