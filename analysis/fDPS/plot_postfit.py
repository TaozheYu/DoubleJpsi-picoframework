import ROOT
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import gPad
import math

#ROOT.gROOT.Reset()
ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetFrameBorderMode(0)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPalette(1,0,1.)
ROOT.gStyle.SetTitleX(0.5) #title X location
ROOT.gStyle.SetTitleY(0.96) #title Y location
ROOT.gStyle.SetPaintTextFormat(".2f")

c1 = ROOT.TCanvas("c1","c1",0,0,900,600)

f = ROOT.TFile.Open("fitDiagnosticsTest.root","READ")
f1 = ROOT.TFile.Open("cross_section_JpsiPair_mass.root","READ")

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-c', '--cutcategory', dest='cutcategory', action='store', type=str, default='JpsiPair_mass')

cutcategory = args.cutcategory

if cutcategory=='JpsiPair_mass':
   bins = [5.,10.,15.,20.,25.,30.]
elif cutcategory=='JpsiPair_DeltaEta':
   bins = [0.,0.5,1.0,1.5,2.0,2.5]
else:
   raise ValueError('"cutcategory" must be"JpsiPair_mass" or "JpsiPair_DetaEta" ')

#DPS = ROOT.TH1F("","",5,5,30)
DPS = ROOT.TH1F("","",len(bins),bin[0],bins[len(bins)-1])
SPS = ROOT.TH1F("","",len(bins),bin[0],bins[len(bins)-1])
total = ROOT.TH1F("","",len(bins),bin[0],bins[len(bins)-1])
data = ROOT.TH1F("","",len(bins),bin[0],bins[len(bins)-1])


DPS_pre = f.Get("shapes_fit_s/cat0/DPS")
SPS_pre = f.Get("shapes_fit_s/cat0/SPS")
total_pre = f.Get("shapes_fit_s/cat0/total")
data_pre = f1.Get("data_obs")

for i in range(data_pre.GetNbinsX()):
 DPS.SetBinContent(i+1,DPS_pre.GetBinContent(i+1))
 DPS.SetBinError(i+1,DPS_pre.GetBinError(i+1))

 SPS.SetBinContent(i+1,SPS_pre.GetBinContent(i+1))
 SPS.SetBinError(i+1,SPS_pre.GetBinError(i+1))

 total.SetBinContent(i+1,total_pre.GetBinContent(i+1))
 total.SetBinError(i+1,total_pre.GetBinError(i+1))

 data.SetBinContent(i+1,data_pre.GetBinContent(i+1))
 data.SetBinError(i+1,data_pre.GetBinError(i+1))

 print data_pre.GetBinError(i+1)

DPS.SetLineWidth(1)
DPS.SetLineStyle(1)
DPS.SetLineColor(ROOT.kGreen)

SPS.SetLineWidth(1)
SPS.SetLineStyle(1)
SPS.SetLineColor(ROOT.kMagenta)

total.SetLineWidth(1)
total.SetLineStyle(1)
total.SetLineColor(ROOT.kRed)

data.SetLineWidth(2) 
data.SetLineColor(1) 
data.SetMarkerColor(1) 
data.SetMarkerStyle(20) 
data.SetMarkerSize(1.3)

DPS.GetYaxis().SetTitle("cross section per bin(pb)")
DPS.GetXaxis().SetTitle("m_{JpsiJpsi} [GeV]")
DPS.SetMaximum(500)
DPS.SetMinimum(0)

DPS.Draw("HIST")
SPS.Draw("HIST same")
total.Draw("HIST same")
data_pre.Draw("same epx0e0")

pad = ROOT.TPad("pad","pad",0.01,0.01,0.99,0.99)
ROOT.gPad.RedrawAxis()
channelText = ""
channelTextFont   = 42
channelTextSize   = 0.060
cmsText     = "CMS"
cmsTextFont   = 61  # default is helvetic-bold
#writeExtraText = true
extraText   = "Simulation"
extraTextFont = 52  # default is helvetica-italics
#text sizes and text offsets with respect to the top frame in unit of the top margin size
lumiTextSize     = 0.7
lumiTextOffset   = 0.2
cmsTextSize      = 0.75
cmsTextOffset    = 0.1  # only used in outOfFrame version
relPosX    = 0.045
relPosY    = 0.035
relExtraDY = 1.2
# ratio of "CMS" and extra text size
extraOverCmsTextSize  = 0.65
lumi_13TeV = "36.3 fb^{-1}"
#lumiText += lumi_13TeV
lumiText ="36.3 fb^{-1}"+" (13 TeV)"
t = pad.GetTopMargin()
b = pad.GetBottomMargin()
r = pad.GetRightMargin()
l = pad.GetLeftMargin()
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
extraTextSize = extraOverCmsTextSize*cmsTextSize
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(lumiTextSize*t)
latex.DrawLatex(1-r,0.92,lumiText)
latex.SetTextFont(cmsTextFont)
latex.SetTextAlign(11)
latex.SetTextSize(cmsTextSize*t)
latex.DrawLatex(l+0.08, 1-t+lumiTextOffset*t-0.12,cmsText)#0.075
latex.SetTextFont(extraTextFont)
latex.SetTextSize(extraTextSize*t)

pl2 = ROOT.TLegend(0.60,0.60,0.88,0.87)
pl2.SetTextSize(0.030)
pl2.SetTextFont(62)
pl2.SetFillColor(0)
pl2.AddEntry(data, "Observed",  "PE")
pl2.AddEntry(total, "Sum of DPS and SPS ", "L")
pl2.AddEntry(DPS, "DPS", "L")
pl2.AddEntry(SPS, "SPS", "L")
pl2.SetBorderSize(0)
pl2.Draw()

print 'Finish'
c1.SaveAs("c1.pdf")

