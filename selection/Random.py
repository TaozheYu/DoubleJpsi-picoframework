import ROOT

for i in range(0,100):
  gRandom = ROOT.TRandom2(0)
  k =gRandom.Uniform(0,1)
  print k 
