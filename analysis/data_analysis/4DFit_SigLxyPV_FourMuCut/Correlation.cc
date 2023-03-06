void Correlation(){
  float E_a = 1.52;   float s_a = 1.406;
  float E_b = 1.504;  float s_b = 1.364;
  float E_ab = 2.764; float V_ab= 8.762*8.762;

  float correlation = abs((E_ab - E_a*E_b)/(s_a*s_b));
  cout<<correlation<<endl;
}
