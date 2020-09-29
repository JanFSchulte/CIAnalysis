import argparse	
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT
import ratios
from setTDRStyle import setTDRStyle
gROOT.SetBatch(True)
from helpers import *
from defs import getPlot, Backgrounds, Backgrounds2016, Backgrounds2018, Signals, Signals2016, Signals2016ADD, SignalsADD, Signals2018ADD, Signals2018, Data, Data2016, Data2018, path, plotList, zScale, zScale2016, zScale2018
import math
import os
from copy import copy


def plotDataMC(datahist,mchist,usedata, label1, label2,name,filename):
	

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	if usedata==True:
                plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
                ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
                setTDRStyle()
                plotPad.UseCurrentStyle()
                ratioPad.UseCurrentStyle()
                plotPad.Draw("hist")
                ratioPad.Draw("hist")
                plotPad.cd()
        else:
                plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
                setTDRStyle()
                plotPad.UseCurrentStyle()
                plotPad.Draw()
                plotPad.cd()	
	colors = createMyColors()		
	
	legend = TLegend(0.55, 0.6, 0.925, 0.925)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextFont(42)
	legendEta = TLegend(0.45, 0.75, 0.925, 0.925)
	legendEta.SetFillStyle(0)
	legendEta.SetBorderSize(0)
	legendEta.SetTextFont(42)
	legendEta.SetNColumns(2)
	latex = ROOT.TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = ROOT.TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.06)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.045)
	latexCMSExtra.SetNDC(True)	
	legendHists = []
	legendHistData = ROOT.TH1F()
	category = ROOT.TLatex()
	category.SetNDC(True)
	category.SetTextSize(0.04)
	if usedata==True:	
		legend.AddEntry(legendHistData,label1,"pe")	
		legendEta.AddEntry(legendHistData,label1,"pe")	

	#process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+^{*}}e^{-}")
	temphist = ROOT.TH1F()
	temphist.SetFillColor(3)
	legendHists.append(temphist.Clone)
	#legend.AddEntry(temphist,label2,"l")
	#legendEta.AddEntry(temphist,process.label,"f")
	
	# Modify plot pad information	
	nEvents=-1

	ROOT.gStyle.SetOptStat(0)
	
	intlumi = ROOT.TLatex()
	intlumi.SetTextAlign(12)
	intlumi.SetTextSize(0.045)
	intlumi.SetNDC(True)
	intlumi2 = ROOT.TLatex()
	intlumi2.SetTextAlign(12)
	intlumi2.SetTextSize(0.07)
	intlumi2.SetNDC(True)
	scalelabel = ROOT.TLatex()
	scalelabel.SetTextAlign(12)
	scalelabel.SetTextSize(0.03)
	scalelabel.SetNDC(True)
	metDiffLabel = ROOT.TLatex()
	metDiffLabel.SetTextAlign(12)
	metDiffLabel.SetTextSize(0.03)
	metDiffLabel.SetNDC(True)
	chi2Label = ROOT.TLatex()
	chi2Label.SetTextAlign(12)
	chi2Label.SetTextSize(0.03)
	chi2Label.SetNDC(True)
	hCanvas.SetLogy()


	# Luminosity information	
	plotPad.cd()
	plotPad.SetLogy(0)
	plotPad.SetLogy()
	if usedata==True:
		yMax = datahist.GetBinContent(datahist.GetMaximumBin())*1000
		yMin = 0.00000001
		xMax = datahist.GetXaxis().GetXmax()
		xMin = datahist.GetXaxis().GetXmin()
	else:	
		yMax = mchist.GetBinContent(datahist.GetMaximumBin())
		yMin = 0.00000001
		xMax = mchist.GetXaxis().GetXmax()
		xMin = mchist.GetXaxis().GetXmin()	
		yMax = yMax*10000
	if name.find("dimuon")!=-1:
		plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; m_{#mu#mu}[GeV] ;Events/GeV")
	else:
		plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; m_{ee}[GeV] ;Events/GeV")
	# Draw signal information
	
	# Draw background from stack
	mchist.SetFillColor(0)
	mchist.SetLineColor(2)
	legend.AddEntry(mchist,label2,"l")
	mchist.Draw("samehist")		

	# Draw data
	datahist.SetMinimum(0.0001)
	if usedata==True:
		datahist.SetMarkerStyle(8)
		datahist.Draw("samepehist")	

	# Draw legend
	legend.Draw()
	plotPad.SetLogx()
	latex.DrawLatex(0.95,0.96,"13 TeV")
	yLabelPos = 0.85
	cmsExtra = "Preliminary"
	if not usedata==True:
		cmsExtra = "#splitline{Preliminary}{Simulation}"
		yLabelPos = 0.82	
	latexCMS.DrawLatex(0.19,0.89,"CMS")
	category.DrawLatex(0.3,0.7,name)
	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
	#~ print datahist.Integral()
	if usedata==True:
		try:
			ratioPad.cd()
			ratioPad.SetLogx()
		except AttributeError:
			print ("Plot fails. Look up in errs/failedPlots.txt")
			outFile =open("errs/failedPlots.txt","a")
			outFile.write('%s\n'%plot.filename%("_"+run.label+"_"+dilepton))
			outFile.close()
			#plot.cuts=baseCut
			return 1
		ratioGraphs =  ratios.RatioGraph(datahist,mchist, xMin=xMin, xMax=xMax,title="After/Before",yMin=0,yMax=2.0,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000000000000,labelSize=0.125,pull=False)
		ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
					

	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	if usedata==True:

		ratioPad.RedrawAxis()
	if not os.path.exists("plots"):
		os.makedirs("plots")	
	hCanvas.Print("plots/%s.jpg"%filename)

					

						
f=ROOT.TFile.Open("RooUnfold/invVsiter.root")

recBB_18_mu=f.Get("datamubb2018")
unfoldBB_18_mu=f.Get("BB2018mu_invert")
unfoldBB_18_muiter1=f.Get("BB2018muiter_1")
unfoldBB_18_muiter3=f.Get("BB2018muiter_3")
unfoldBB_18_muiter5=f.Get("BB2018muiter_5")
unfoldBB_18_muiter10=f.Get("BB2018muiter_10")
unfoldBB_18_muiter100=f.Get("BB2018muiter_100")
unfoldBB_18_muiter1000=f.Get("BB2018muiter_1000")
unfoldBB_18_muiter10000=f.Get("BB2018muiter_10000")

plotDataMC(unfoldBB_18_mu,recBB_18_mu,True,"Unfolded","Reconsructed","2018 BB dimuon","data_BB_18_mu_invt" )
plotDataMC(unfoldBB_18_muiter1,recBB_18_mu,True,"Unfolded_iter1","Reconstructed","2018 BB dimuon","data_BE_18_muiter1" )
plotDataMC(unfoldBB_18_muiter3,recBB_18_mu,True,"Unfolded_iter3","Reconstructed","2018 BB dimuon","data_BE_18_muiter3" )
plotDataMC(unfoldBB_18_muiter5,recBB_18_mu,True,"Unfolded_iter5","Reconstructed","2018 BB dimuon","data_BE_18_muiter5" )
plotDataMC(unfoldBB_18_muiter10,recBB_18_mu,True,"Unfolded_iter10","Reconstructed","2018 BB dimuon","data_BE_18_muiter10" )
plotDataMC(unfoldBB_18_muiter100,recBB_18_mu,True,"Unfolded_iter100","Reconstructed","2018 BB dimuon","data_BE_18_muiter100" )
plotDataMC(unfoldBB_18_muiter1000,recBB_18_mu,True,"Unfolded_iter1000","Reconstructed","2018 BB dimuon","data_BE_18_muiter1000" )
plotDataMC(unfoldBB_18_muiter10000,recBB_18_mu,True,"Unfolded_iter10000","Reconstructed","2018 BB dimuon","data_BE_18_muiter10000" )

f2=ROOT.TFile.Open("RooUnfold/invVsiter_obkg.root")
bkgBB_18_mu=f2.Get("bkgmubb2018")
recBB_18_mu1=f2.Get("datamubb2018")
unfolditer=f2.Get("BB2018muiter_10000")
unfoldInv=f2.Get("BB2018mu_invert")
recBB_18_mu1.Add(bkgBB_18_mu,-1)
for i in range(recBB_18_mu.GetNbinsX()+2):
	val=recBB_18_mu1.GetBinContent(i)
	err=recBB_18_mu1.GetBinError(i)
	wid=recBB_18_mu1.GetBinWidth(i)
	val=val/wid
	err=err/wid
	recBB_18_mu1.SetBinContent(i,val)
	recBB_18_mu1.SetBinError(i,err)
plotDataMC(unfoldInv,recBB_18_mu1,True,"Unfolded","Reconsructed","2018 BB dimuon","BB_18_mu_invt" )
plotDataMC(unfolditer,recBB_18_mu1,True,"Unfoldediter10000","Reconstructed","2018 BB dimuon","BB_18_muiter" )
	

