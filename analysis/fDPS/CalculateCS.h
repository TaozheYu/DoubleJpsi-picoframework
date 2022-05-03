void CalculateCS(const char* sampletype, float Frac1, float Frac1_err, float Data_num){
   /*
   Double_t data_num = 22550; 
   Double_t frac1 = 0.796;
   Double_t frac1_err = 0.009;

   Double_t frac2 = 0.798;
   Double_t frac2_err = 0.013;
   */
   string Sampletype = sampletype;

   Double_t data_num = Data_num;//401615; 
   Double_t frac1 = Frac1;//0.750;
   Double_t frac1_err = Frac1_err;//0.002;

   Double_t frac2 = Frac1;//0.765;
   Double_t frac2_err = Frac1_err;//0.005;
       
   Double_t frac12 = frac1*frac2;
   Double_t frac12_err = sqrt(frac1*frac1*frac2_err*frac2_err + frac2*frac2*frac1_err*frac1_err);
 
   float N_corr; 
   float N_corr_err; 
   
   if(Sampletype.find("data")!=string::npos){  
     N_corr = frac12*data_num;
     N_corr_err = frac12_err*data_num;                                                
   } else{
     N_corr =data_num;
     N_corr_err = sqrt(data_num);
   }
                                   
   cout<<"double Jpsi events is "<<N_corr<<endl;
   cout<<"double Jpsi events error is "<<N_corr_err<<endl;

 
   float Lumi = 36.3;

   float Br = 0.0593;
   float Br_err = 0.0006;

   float cross_section = N_corr/(Lumi*Br*Br);
   float Br_2_err = sqrt(2*Br*Br*Br_err*Br_err);

   float cross_section_err = sqrt((N_corr*N_corr*Br_2_err*Br_2_err + N_corr_err*N_corr_err*Br*Br*Br*Br)/(Br*Br*Br*Br*Br*Br*Br*Br))/Lumi;

   cout<<"the Jpsi pair production cross section is "<<cross_section/1000000<<endl;//unit is nb 
   cout<<"the Jpsi pair production cross section error is "<<cross_section_err/1000000<<endl;//unit is nb


}
