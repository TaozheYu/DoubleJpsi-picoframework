void Calculate(){
    /*
    Double_t data_num = 22550; 
    Double_t frac1 = 0.796;
    Double_t frac1_err = 0.009;

    Double_t frac2 = 0.798;
    Double_t frac2_err = 0.013;
    */
    
    Double_t data_num = 401615; 
    Double_t frac1 = 0.750;
    Double_t frac1_err = 0.002;

    Double_t frac2 = 0.765;
    Double_t frac2_err = 0.005;
       
    Double_t frac12 = frac1*frac2;
    Double_t frac12_err = sqrt(frac1*frac1*frac2_err*frac2_err + frac2*frac2*frac1_err*frac1_err);
 
   cout<<"double Jpsi events is "<<frac12*data_num<<endl;
   cout<<"double Jpsi events error is "<<frac12_err*data_num<<endl;

   float N_corr = frac12*data_num;
   float N_corr_err = frac12_err*data_num;
 
   float Lumi = 36.3;

   float Br = 0.0593;
   float Br_err = 0.0006;

   float cross_section = N_corr/(Lumi*Br*Br);
   float Br_2_err = sqrt(2*Br*Br*Br_err*Br_err);

   float cross_section_err = sqrt((N_corr*N_corr*Br_2_err*Br_2_err + N_corr_err*N_corr_err*Br*Br*Br*Br)/(Br*Br*Br*Br*Br*Br*Br*Br))/Lumi;

   cout<<"the Jpsi pair production is "<<cross_section/1000000<<"+-"<<cross_section_err/1000000; 


}
