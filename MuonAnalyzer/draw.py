import ROOT
import sys

if len(sys.argv) != 3:
	print "USAGE: %s <input file> <output file>"%(sys.argv[0])
	sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]

print "Reading from", inFileName, " and writing to ", outFileName


inFile = ROOT.TFile.Open(inFileName, "READ")

h_dr_Pt_Sum = inFile.Get("h_dr_Pt_Sum")
h_dz_drho = inFile.Get("h_dz_drho")

if not h_dr_Pt_Sum:
	print "Failed to get hist h_dr_Pt_Sum"
	sys.exit(1)
if not h_dz_drho:
	print "Failed to get hist h_dz_drho"

h_dr_Pt_Sum.SetDirectory(0)
h_dz_drho.SetDirectory(0)

inFile.Close()

fileName = inFileName.replace(".root", "")

c1 = ROOT.TCanvas("c1")
c1.cd()
h_dr_Pt_Sum.Draw("h")
c1.Print("img/" + fileName + "plot1.pdf")

c2 = ROOT.TCanvas("c2")
c2.cd()
h_dz_drho.Draw("h")
c2.Print("img/" + fileName + "plot2.pdf")

h_dr_Pt_Sum_y = h_dr_Pt_Sum.ProjectionY("Pt_Sum", 10 , 10)
c3 = ROOT.TCanvas("c3")
c3.cd()
#h_dr_Pt_Sum_y.Rebin(5)
h_dr_Pt_Sum_y.GetXaxis().SetRangeUser(0., 1.)

h_dr_Pt_Sum_y.Draw("h")
c3.Print("img/" + fileName + "plot3.pdf")


outHistFile = ROOT.TFile.Open(outFileName, "RECREATE")
outHistFile.cd()

outHistFile.Close()
