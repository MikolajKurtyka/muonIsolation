import ROOT
import sys
import os

if len(sys.argv) != 3:
	print "USAGE: %s <input file> <output file>"%(sys.argv[0])
	sys.exit(1)

inFileName = sys.argv[1] 
outFileName = sys.argv[2]

print "Reading from", inFileName, " and writing to ", outFileName

dir_name = "jpt_cut_03/"
#dir_name = "img/"


inFile = ROOT.TFile.Open(inFileName, "READ")

def draw_projection(hist_name, dr):
	h = inFile.Get(hist_name)
	h.SetDirectory(0)
	h_y = h.ProjectionY(hist_name , dr , dr)
	c = ROOT.TCanvas("c")
	c.cd()
	h_y.Rebin(2)
	h_y.GetXaxis().SetRangeUser(0., 1.)
	h_y.Draw("h")
	c.Print(dir_name + fileName + hist_name + "plot3.pdf")

def draw_projection_eff(hist_name, dr):
	h = inFile.Get(hist_name)
	if not h:
		print "Failed to get " + hist_name
		sys.exit(1)
	h.SetDirectory(0)
	c = ROOT.TCanvas("c")
	c.cd()
	h_y_eff = h.ProjectionY("Pt_Sum", dr , dr)
	h_y_eff.Rebin(2)
	h_y_eff.GetXaxis().SetRangeUser(0., 1.)
	h_y_eff.Scale(1./h_y_eff.GetEntries())
	for i in range(1, h_y_eff.GetNbinsX()):
		value = h_y_eff.Integral(i-1, i)
		h_y_eff.SetBinContent(i, value)
	h_y_eff.GetXaxis().SetTitle("PtSumInCone/MuonPt")	
	h_y_eff.GetYaxis().SetTitle("Eff")
	h_y_eff.Draw("h")
	c.Print(dir_name + fileName + hist_name + "plotEff.pdf")



h_dr_Pt_Sum = inFile.Get("h_dr_Pt_Sum")
h_dz_drho = inFile.Get("h_dz_drho")
h_dz = inFile.Get("h_dz")

if not h_dr_Pt_Sum:
	print "Failed to get hist h_dr_Pt_Sum"
	sys.exit(1)
if not h_dz_drho:
	print "Failed to get hist h_dz_drho"

h_dr_Pt_Sum.SetDirectory(0)
h_dz_drho.SetDirectory(0)
h_dz.SetDirectory(0)



fileName = inFileName.replace(".root", "")

c1 = ROOT.TCanvas("c1")
c1.cd()
h_dr_Pt_Sum.Draw("h")
c1.Print(dir_name + fileName + "plot1.pdf")

c2 = ROOT.TCanvas("c2")
c2.cd()
h_dz_drho.Draw("h")
c2.Print(dir_name + fileName + "plot2.pdf")

h_dr_Pt_Sum_y = h_dr_Pt_Sum.ProjectionY("Pt_Sum", 30 , 30)
c3 = ROOT.TCanvas("c3")
c3.cd()
h_dr_Pt_Sum_y.SetTitle(fileName)
h_dr_Pt_Sum_y.Rebin(2)
h_dr_Pt_Sum_y.GetXaxis().SetRangeUser(0., 1.)
h_dr_Pt_Sum_y.Draw("h")
c3.Print(dir_name + fileName + "plot3.pdf")




c4 = ROOT.TCanvas("c4")
c4.cd()
h_dr_Pt_Sum_y_eff = h_dr_Pt_Sum.ProjectionY("Pt_Sum", 30 , 30)
h_dr_Pt_Sum_y_eff.SetTitle(fileName)
h_dr_Pt_Sum_y_eff.Rebin(2)
h_dr_Pt_Sum_y_eff.GetXaxis().SetRangeUser(0., 1.)
h_dr_Pt_Sum_y_eff.Scale(1./h_dr_Pt_Sum_y.GetEntries())
for i in range(1, h_dr_Pt_Sum_y_eff.GetNbinsX()):
	value = h_dr_Pt_Sum_y_eff.Integral(i-1, i)
	h_dr_Pt_Sum_y_eff.SetBinContent(i, value)
h_dr_Pt_Sum_y_eff.GetXaxis().SetTitle("PtSumInCone/MuonPt")
h_dr_Pt_Sum_y_eff.GetYaxis().SetTitle("Eff")
h_dr_Pt_Sum_y_eff.Draw("h")
c4.Print(dir_name + fileName + "plotEff.pdf")

draw_projection_eff("h_dr_Pt_Sum_Ch", 30)
#draw_projection("h_dr_Pt_Sum_Ch", 30)
draw_projection_eff("h_dr_Pt_Sum_dz05", 30)
#draw_projection("h_dr_Pt_Sum_dz05", 30)
draw_projection_eff("h_dr_Pt_Sum_dz05_Ch", 30)
#draw_projection("h_dr_Pt_Sum_dz05_Ch", 30)

c5 = ROOT.TCanvas("c5")
c5.cd()
h_dz.SetTitle(fileName)
h_dz.Rebin(5)
h_dz.GetXaxis().SetRangeUser(-0.3,0.3)
h_dz.Draw("h")
c5.Print(dir_name + fileName + "plot4.pdf")

inFile.Close()
outHistFile = ROOT.TFile.Open(outFileName, "RECREATE")
outHistFile.cd()

outHistFile.Close()
