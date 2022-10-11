void Calculate_corrector(){
  //gROOT->Reset();
  //gStyle->SetCanvasColor(0);
  //gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat("M");
  //gStyle->SetPalette(1,0);
  //gStyle->SetPaintTextFormat(".2f");
  //gErrorIgnoreLevel = kError;
  using namespace std;
 
  float Cal_stat_err(float N_num, float N_den);
  const int N = 11;  Double_t xbins[N] = {5.,6.,7.,8.,9.,10.,12.,15.,20.,25.,40.};
  const int M = 11;   Double_t ybins[M] = {-2.,-1.75,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,1.75,2.};
 
  TH2F *h_w_acc   = new TH2F("h_w_acc",  "h_w_acc", N-1, xbins, M-1, ybins);
  TH2F *h_w_reco   = new TH2F("h_w_reco",  "h_w_reco", N-1, xbins, M-1, ybins);
  TH2F *h_w_eff   = new TH2F("h_w_eff",  "h_w_eff", N-1, xbins, M-1, ybins);
  TH2F *h_w_vtx   = new TH2F("h_w_vtx",  "h_w_vtx", N-1, xbins, M-1, ybins);

  const int pt_bin = 7; Double_t pt_bins[pt_bin] = {5.,7.,9.,12.,15.,25.,40.}; 
  TH2F *h_w_HLT   = new TH2F("h_w_HLT",  "h_w_HLT", pt_bin-1, pt_bins, pt_bin-1, pt_bins);
  TH2F *h_w_evt   = new TH2F("h_w_evt",  "h_w_evt", pt_bin-1, pt_bins, pt_bin-1, pt_bins);



  TFile *fileJpsi_gen = new TFile("Gen_Event.root");
  
  TH2F *h_total_Jpsi = (TH2F*) fileJpsi_gen->Get("h_total_Jpsi");
  TH2F *h_acc_Jpsi = (TH2F*) fileJpsi_gen->Get("h_acc_Jpsi");
  TH2F *h_reco_Jpsi = (TH2F*) fileJpsi_gen->Get("h_reco_Jpsi");
  TH2F *h_eff_Jpsi = (TH2F*) fileJpsi_gen->Get("h_eff_Jpsi");
  TH2F *h_vtx_Jpsi = (TH2F*) fileJpsi_gen->Get("h_vtx_Jpsi");
  TH2F *h_JpsiPair = (TH2F*) fileJpsi_gen->Get("h_JpsiPair");
  TH2F *h_HLT_JpsiPair = (TH2F*) fileJpsi_gen->Get("h_HLT_JpsiPair");
  TH2F *h_evt_JpsiPair = (TH2F*) fileJpsi_gen->Get("h_evt_JpsiPair");

  cout<<h_total_Jpsi->GetXaxis()->GetNbins()<<endl;
  cout<<h_total_Jpsi->GetYaxis()->GetNbins()<<endl;


  for(int i=1; i<h_total_Jpsi->GetXaxis()->GetNbins()+1; i++){
    for(int j=1; j<h_total_Jpsi->GetYaxis()->GetNbins()+1; j++){
    
    float N_total = h_total_Jpsi->GetBinContent(i,j);
    float N_acc   = h_acc_Jpsi  ->GetBinContent(i,j);
    float N_reco  = h_reco_Jpsi ->GetBinContent(i,j);
    float N_eff   = h_eff_Jpsi  ->GetBinContent(i,j);
    float N_vtx   = h_vtx_Jpsi  ->GetBinContent(i,j);

    if (N_total!=0&&N_acc!=0&&N_reco!=0&&N_eff!=0&&N_vtx!=0){ 
      float N_total_err = sqrt(N_total);
      float N_acc_err = sqrt(N_acc);
      float N_reco_err = sqrt(N_reco);
      float N_eff_err = sqrt(N_eff);

      float w_acc = N_acc/N_total;
      float w_reco= N_reco/N_acc;
      float w_eff = N_eff/N_reco;
      float w_vtx = N_vtx/N_eff;

      float w_acc_err = sqrt((N_acc_err*N_acc_err*N_total*N_total + N_total_err*N_total_err*N_acc*N_acc)/(N_total*N_total*N_total*N_total));
      float w_reco_err= sqrt((N_reco*N_reco*N_acc_err*N_acc_err + N_reco_err*N_reco_err*N_acc*N_acc)/(N_acc*N_acc*N_acc*N_acc));
      float w_eff_err = sqrt((N_eff*N_eff*N_reco_err*N_reco_err + N_eff_err*N_eff_err*N_reco*N_reco)/(N_reco*N_reco*N_reco*N_reco));

      float w_vtx_err = Cal_stat_err(N_vtx,N_eff);

      h_w_acc->SetBinContent(i,j,w_acc);
      h_w_reco->SetBinContent(i,j,w_reco);
      h_w_eff->SetBinContent(i,j,w_eff);
      h_w_vtx->SetBinContent(i,j,w_vtx);

      h_w_acc->SetBinError(i,j,w_acc_err);
      h_w_reco->SetBinError(i,j,w_reco_err);
      h_w_eff->SetBinError(i,j,w_eff_err);
      h_w_vtx->SetBinError(i,j,w_vtx_err);
      }
    }
  }

  for(int i=1; i<h_HLT_JpsiPair->GetXaxis()->GetNbins()+1; i++){
    for(int j=1; j<h_HLT_JpsiPair->GetYaxis()->GetNbins()+1; j++){
    
    float N_JpsiPair = (h_JpsiPair ->GetBinContent(i,j)+h_JpsiPair ->GetBinContent(j,i))/2;
    float N_HLT_JpsiPair = (h_HLT_JpsiPair ->GetBinContent(i,j)+h_HLT_JpsiPair ->GetBinContent(j,i))/2;
    float N_evt_JpsiPair = (h_evt_JpsiPair ->GetBinContent(i,j)+h_evt_JpsiPair ->GetBinContent(j,i))/2;
 
    //float N_betri = h_betri_Jpsi ->GetBinContent(i,j);
    //float N_patri = h_patri_Jpsi ->GetBinContent(i,j);
 
    if(N_JpsiPair!=0&&N_HLT_JpsiPair!=0&&N_evt_JpsiPair!=0){

     
      float w_HLT = N_HLT_JpsiPair/N_JpsiPair;
      float w_evt = N_evt_JpsiPair/N_HLT_JpsiPair;
      
      float w_HLT_err = Cal_stat_err(N_HLT_JpsiPair,N_JpsiPair);
      float w_evt_err = Cal_stat_err(N_evt_JpsiPair,N_HLT_JpsiPair);    
 
      h_w_HLT->SetBinContent(i,j,w_HLT);
      h_w_evt->SetBinContent(i,j,w_evt);

      h_w_HLT->SetBinError(i,j,w_HLT_err);
      h_w_evt->SetBinError(i,j,w_evt_err);

      }
    } 
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",0,0,900,600);
  h_w_acc->GetYaxis()->SetTitleSize(0.040);
  h_w_acc->GetXaxis()->SetTitleSize(0.040);
  h_w_acc->GetYaxis()->SetLabelSize(0.040);
  h_w_acc->GetXaxis()->SetLabelSize(0.040);
  h_w_acc->SetMinimum(0.8);
  h_w_acc->SetMaximum(1.0);
  h_w_acc->SetTitle("acceptance probability");
  h_w_acc->GetYaxis()->SetTitle("Jpsi eta");
  h_w_acc->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_acc->Draw("COLZ");
  //h_w_acc->Draw("COLZTEXTE");
  //h_w_acc->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("acc.png");
  canvas->SaveAs("acc.pdf");
    
  h_w_reco->GetYaxis()->SetTitleSize(0.040);
  h_w_reco->GetXaxis()->SetTitleSize(0.040);
  h_w_reco->GetYaxis()->SetLabelSize(0.040);
  h_w_reco->GetXaxis()->SetLabelSize(0.040);
  h_w_reco->SetMinimum(0.0);
  h_w_reco->SetMaximum(0.7);
  h_w_reco->SetTitle("reconstruct probability");
  h_w_reco->GetYaxis()->SetTitle("Jpsi eta");
  h_w_reco->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_reco->Draw("COLZ");
  //h_w_reco->Draw("COLZTEXTE");
  //h_w_reco->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("reco.png");
  canvas->SaveAs("reco.pdf");

  h_w_eff->GetYaxis()->SetTitleSize(0.040);
  h_w_eff->GetXaxis()->SetTitleSize(0.040);
  h_w_eff->GetYaxis()->SetLabelSize(0.040);
  h_w_eff->GetXaxis()->SetLabelSize(0.040);
  h_w_eff->SetMinimum(0.4);
  h_w_eff->SetMaximum(1.0);
  h_w_eff->SetTitle("selection probability");
  h_w_eff->GetYaxis()->SetTitle("Jpsi eta");
  h_w_eff->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_eff->Draw("COLZ");
  //h_w_eff->Draw("COLZTEXTE");
  //h_w_eff->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("eff.png");
  canvas->SaveAs("eff.pdf");

  h_w_vtx->GetYaxis()->SetTitleSize(0.040);
  h_w_vtx->GetXaxis()->SetTitleSize(0.040);
  h_w_vtx->GetYaxis()->SetLabelSize(0.040);
  h_w_vtx->GetXaxis()->SetLabelSize(0.040);
  h_w_vtx->SetMinimum(0.4);
  h_w_vtx->SetMaximum(1.0);
  h_w_vtx->SetTitle("vertex probability");
  h_w_vtx->GetYaxis()->SetTitle("Jpsi eta");
  h_w_vtx->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_vtx->Draw("COLZ");
  //h_w_vtx->Draw("COLZTEXTE");
  //h_w_vtx->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("vtx.png");
  canvas->SaveAs("vtx.pdf");

  TCanvas* canvas2 = new TCanvas("canvas","canvas",0,0,600,600);
  h_w_HLT->GetYaxis()->SetTitleSize(0.040);
  h_w_HLT->GetXaxis()->SetTitleSize(0.040);
  h_w_HLT->GetYaxis()->SetLabelSize(0.040);
  h_w_HLT->GetXaxis()->SetLabelSize(0.040);
  h_w_HLT->SetMinimum(0.);
  h_w_HLT->SetMaximum(1.0);
  h_w_HLT->SetTitle(" pass HLT probability");
  h_w_HLT->GetYaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_HLT->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_HLT->Draw("COLZ");
  //h_w_HLT->Draw("COLZTEXTE");
  //h_w_HLT->Draw("COLZTEXT");
  canvas2->Draw();
  canvas2->SaveAs("HLT.png");
  canvas2->SaveAs("HLT.pdf");

  h_w_evt->GetYaxis()->SetTitleSize(0.040);
  h_w_evt->GetXaxis()->SetTitleSize(0.040);
  h_w_evt->GetYaxis()->SetLabelSize(0.040);
  h_w_evt->GetXaxis()->SetLabelSize(0.040);
  h_w_evt->SetMinimum(0.);
  h_w_evt->SetMaximum(1.0);
  h_w_evt->SetTitle(" pass evt probability");
  h_w_evt->GetYaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_evt->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  h_w_evt->Draw("COLZ");
  //h_w_evt->Draw("COLZTEXTE");
  //h_w_evt->Draw("COLZTEXT");
  canvas2->Draw();
  canvas2->SaveAs("evt.png");
  canvas2->SaveAs("evt.pdf");

  char NewFileName[500]; sprintf(NewFileName, "corrector.root");
  TFile f(NewFileName,"new");
  h_w_acc->Write();
  h_w_reco->Write();
  h_w_eff->Write();
  h_w_vtx->Write();
  h_w_HLT->Write();
  h_w_evt->Write();
  f.Close();

}

float Cal_stat_err(float N_num, float N_den){
  //value=num/den;
  if (N_num==0||N_den==0) return 0;
  float N_num_err = sqrt(N_num);
  float N_den_err = sqrt(N_den);
  float err = sqrt((N_num_err*N_num_err*N_den*N_den + N_den_err*N_den_err*N_num*N_num)/(N_den*N_den*N_den*N_den));
  return err;
}





