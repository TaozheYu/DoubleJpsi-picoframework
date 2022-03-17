void Calculate_corrector(){
  

  const int N = 10;  Double_t xbins[N] = {5.,6.,7.,8.,10.,12.,15.,20.,30.,40.};
  const int M = 19;   Double_t ybins[M] = {-2.2,-2.,-1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.,2.2};

  TH2F *h_w_acc   = new TH2F("h_w_acc",  "h_w_acc", N-1, xbins, M-1, ybins);
  TH2F *h_w_reco   = new TH2F("h_w_reco",  "h_w_reco", N-1, xbins, M-1, ybins);
  TH2F *h_w_eff   = new TH2F("h_w_eff",  "h_w_eff", N-1, xbins, M-1, ybins);



  TFile *fileJpsi_gen = new TFile("Gen_Events.root");
  
  TH2F *h_total_Jpsi = (TH2F*) fileJpsi_gen->Get("plots/h_total_Jpsi");
  TH2F *h_acc_Jpsi = (TH2F*) fileJpsi_gen->Get("plots/h_acc_Jpsi");
  TH2F *h_reco_Jpsi = (TH2F*) fileJpsi_gen->Get("plots/h_reco_Jpsi");
  TH2F *h_eff_Jpsi = (TH2F*) fileJpsi_gen->Get("plots/h_eff_Jpsi");

  cout<<h_total_Jpsi->GetXaxis()->GetNbins()<<endl;
  cout<<h_total_Jpsi->GetYaxis()->GetNbins()<<endl;


  for(int i=1; i<h_total_Jpsi->GetXaxis()->GetNbins()+1; i++){
    for(int j=1; j<h_total_Jpsi->GetYaxis()->GetNbins()+1; j++){
    
    float N_total = h_total_Jpsi->GetBinContent(i,j);
    float N_acc   = h_acc_Jpsi  ->GetBinContent(i,j);
    float N_reco  = h_reco_Jpsi ->GetBinContent(i,j);
    float N_eff   = h_eff_Jpsi  ->GetBinContent(i,j);

    if (N_total!=0&&N_acc!=0&&N_reco!=0&&N_eff!=0){ 
      float N_total_err = sqrt(N_total);
      float N_acc_err = sqrt(N_acc);
      float N_reco_err = sqrt(N_reco);
      float N_eff_err = sqrt(N_eff);

      float w_acc = N_acc/N_total;
      float w_reco= N_reco/N_acc;
      float w_eff = N_eff/N_reco;

      float w_acc_err = sqrt((N_acc_err*N_acc_err*N_total*N_total + N_total_err*N_total_err*N_acc*N_acc)/(N_total*N_total*N_total*N_total));
      float w_reco_err= sqrt((N_reco*N_reco*N_acc_err*N_acc_err + N_reco_err*N_reco_err*N_acc*N_acc)/(N_acc*N_acc*N_acc*N_acc));
      float w_eff_err = sqrt((N_eff*N_eff*N_reco_err*N_reco_err + N_eff_err*N_eff_err*N_reco*N_reco)/(N_reco*N_reco*N_reco*N_reco));

      h_w_acc->SetBinContent(i,j,w_acc);
      h_w_reco->SetBinContent(i,j,w_reco);
      h_w_eff->SetBinContent(i,j,w_eff);

      h_w_acc->SetBinError(i,j,w_acc_err);
      h_w_reco->SetBinError(i,j,w_reco_err);
      h_w_eff->SetBinError(i,j,w_eff_err);
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
  //h_w_acc->Draw("COLZ");
  //h_w_acc->Draw("COLZTEXTE");
  h_w_acc->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("acc.png");
  canvas->SaveAs("acc.pdf");
    
  h_w_reco->GetYaxis()->SetTitleSize(0.040);
  h_w_reco->GetXaxis()->SetTitleSize(0.040);
  h_w_reco->GetYaxis()->SetLabelSize(0.040);
  h_w_reco->GetXaxis()->SetLabelSize(0.040);
  h_w_reco->SetMinimum(0.0);
  h_w_reco->SetMaximum(0.7);
  h_w_reco->SetTitle("rconstruct probability");
  h_w_reco->GetYaxis()->SetTitle("Jpsi eta");
  h_w_reco->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  //h_w_reco->Draw("COLZ");
  //h_w_reco->Draw("COLZTEXTE");
  h_w_reco->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("reco.png");
  canvas->SaveAs("reco.pdf");

  h_w_eff->GetYaxis()->SetTitleSize(0.040);
  h_w_eff->GetXaxis()->SetTitleSize(0.040);
  h_w_eff->GetYaxis()->SetLabelSize(0.040);
  h_w_eff->GetXaxis()->SetLabelSize(0.040);
  h_w_eff->SetMinimum(0.8);
  h_w_eff->SetMaximum(1.0);
  h_w_eff->SetTitle("selection eff probability");
  h_w_eff->GetYaxis()->SetTitle("Jpsi eta");
  h_w_eff->GetXaxis()->SetTitle("Jpsi pt [GeV]");
  //h_w_eff->Draw("COLZ");
  //h_w_eff->Draw("COLZTEXTE");
  h_w_eff->Draw("COLZTEXT");
  canvas->Draw();
  canvas->SaveAs("eff.png");
  canvas->SaveAs("eff.pdf");
}





