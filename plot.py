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


def plotDataMC(args,plot):
	

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	if args.ratio:
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
		setTDRStyle()		
		plotPad.UseCurrentStyle()
		ratioPad.UseCurrentStyle()
		plotPad.Draw()	
		ratioPad.Draw()	
		plotPad.cd()
	else:
		plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
		setTDRStyle()
		plotPad.UseCurrentStyle()
		plotPad.Draw()	
		plotPad.cd()	
		
	# Data load processes
	colors = createMyColors()		
	if args.use2016:
		data = Process(Data2016, normalized=True)
	elif args.use2018:
		data = Process(Data2018, normalized=True)
	elif args.useRun2:
		data2016 = Process(Data2016, normalized=True)
		data2017 = Process(Data, normalized=True)
		data2018 = Process(Data2018, normalized=True)
	else:	
		data = Process(Data, normalized=True)

	eventCounts = totalNumberOfGeneratedEvents(path,plot.muon)	
	negWeights = negWeightFractions(path,plot.muon)

	# Background load processes	
	backgrounds = copy(args.backgrounds)
	if plot.useJets:
		if "Wjets" in backgrounds:
			backgrounds.remove("Wjets")
		if not "Jets" in backgrounds:
			backgrounds.insert(0,"Jets")
	if args.use2016 and not plot.muon and "Other" in backgrounds:
		backgrounds.remove("Other")
		backgrounds.insert(1,"OtherEle")
	processes = []
	processes2016 = []
	processes2017 = []
	processes2018 = []
	for background in backgrounds:
		if args.use2016:
			if background == "Jets":
				processes.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
			elif background == "Other" and not plot.muon:
				processes.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts,negWeights))				
			else:	
				processes.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
		elif args.use2018:
			if background == "Jets":
				processes.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights,normalized=True))
			else:	
				processes.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))
		elif args.useRun2:

			if background == "Jets":
				processes.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
				processes2016.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
				processes2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
				processes2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights,normalized=True))
			elif background == "Other" and not plot.muon:
				processes.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts,negWeights))
				processes2016.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts,negWeights))
				processes2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
				processes2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))
			else:	
				processes.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
				processes2016.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
				processes2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
				processes2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))
		else:
			if background == "Jets":
				processes.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
			else:	
				processes.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
	
	# Signal load processes
	signals = []
	signals2016 = []
	signals2017 = []
	signals2018 = []	
	for signal in args.signals:
		if args.use2016:
			if args.ADD: signals.append(Process(getattr(Signals2016ADD, signal), eventCounts, negWeights))
			else: signals.append(Process(getattr(Signals2016,signal),eventCounts,negWeights))
		elif args.use2018:
			if args.ADD: signals.append(Process(getattr(Signals2018ADD, signal), eventCounts, negWeights))
			else: signals.append(Process(getattr(Signals2018,signal),eventCounts,negWeights))
		elif args.useRun2:
			if args.ADD: 
					signals.append(Process(getattr(Signals2016ADD, signal), eventCounts, negWeights))
					signals2016.append(Process(getattr(Signals2016ADD, signal), eventCounts, negWeights))
					signals2017.append(Process(getattr(SignalsADD, signal), eventCounts, negWeights))
					signals2018.append(Process(getattr(Signals2018ADD, signal), eventCounts, negWeights))
			else: 
					signals.append(Process(getattr(Signals2016,signal),eventCounts,negWeights))
					signals2016.append(Process(getattr(Signals2016,signal),eventCounts,negWeights))
					signals2017.append(Process(getattr(Signals,signal),eventCounts,negWeights))
					signals2018.append(Process(getattr(Signals2018,signal),eventCounts,negWeights))
		else:	
			if args.ADD: signals.append(Process(getattr(SignalsADD, signal), eventCounts, negWeights))
			else: signals.append(Process(getattr(Signals,signal),eventCounts,negWeights))
		
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
	
	# Modify legend information
	legendHistData = ROOT.TH1F()
	if args.data:	
		legend.AddEntry(legendHistData,"Data","pe")	
		legendEta.AddEntry(legendHistData,"Data","pe")	

	for process in reversed(processes):
		if not plot.muon and "#mu^{+}#mu^{-}" in process.label:
			process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+^{*}}e^{-}")
		temphist = ROOT.TH1F()
		temphist.SetFillColor(process.theColor)
		legendHists.append(temphist.Clone)
		legend.AddEntry(temphist,process.label,"f")
		legendEta.AddEntry(temphist,process.label,"f")
	
	if args.signals !=0:
		processesWithSignal = []
		for process in processes:
			processesWithSignal.append(process)
		for Signal in signals:
			processesWithSignal.append(Signal)
			temphist = ROOT.TH1F()
			temphist.SetFillColor(Signal.theColor)
			temphist.SetLineColor(Signal.theLineColor)
			legendHists.append(temphist.Clone)		
			legend.AddEntry(temphist,Signal.label,"l")
			legendEta.AddEntry(temphist,Signal.label,"l")
	

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
	logScale = plot.log
	
	if logScale == True:
		plotPad.SetLogy()
	if args.use2016:	
		lumi = 35.9*1000
		if plot.muon:
			lumi = 36294.6
	elif args.use2018:	
		lumi = 59401
		if plot.muon:
			lumi = 61298.775231718995
	elif args.useRun2:	
		lumi2016 = 35.9*1000
		lumi2017 = 41.529*1000
		lumi2018 = 59401
		if plot.muon:
			lumi2016 = 36294.6
			lumi2018 = 61298.775231718995
			lumi2017 = 42079.880396
	else:
		lumi = 41.529*1000
		if plot.muon:
			lumi = 42079.880396
	if args.use2016:		
		zScaleFac = zScale2016["muons"]
		if not plot.muon:
			if "bbbe" in plot.histName:
				zScaleFac = zScale2016["electrons"][0]
			elif "bb" in plot.histName:
				zScaleFac = zScale2016["electrons"][1]
			elif "be" in plot.histName:
				zScaleFac = zScale2016["electrons"][2]
			else:
				zScaleFac = zScale2016["electrons"][0]
	elif args.use2018:		
		zScaleFac = zScale2018["muons"]
		if not plot.muon:
			if "bbbe" in plot.histName:
				zScaleFac = zScale2018["electrons"][0]
			elif "bb" in plot.histName:
				zScaleFac = zScale2018["electrons"][1]
			elif "be" in plot.histName:
				zScaleFac = zScale2018["electrons"][2]
			else:
				zScaleFac = zScale2018["electrons"][0]

	elif args.useRun2:		
		zScaleFac2016 = zScale2016["muons"]
		zScaleFac2017 = zScale["muons"]
		zScaleFac2018 = zScale2018["muons"]
		if not plot.muon:
			if "bbbe" in plot.histName:
				zScaleFac2016 = zScale2016["electrons"][0]
				zScaleFac2017 = zScale["electrons"][0]
				zScaleFac2018 = zScale2018["electrons"][0]
			elif "bb" in plot.histName:
				zScaleFac2016 = zScale2016["electrons"][1]
				zScaleFac2017 = zScale["electrons"][1]
				zScaleFac2018 = zScale2018["electrons"][1]
			elif "be" in plot.histName:
				zScaleFac2016 = zScale2016["electrons"][2]
				zScaleFac2017 = zScale["electrons"][2]
				zScaleFac2018 = zScale2018["electrons"][2]
			else:
				zScaleFac2016 = zScale2016["electrons"][0]			
				zScaleFac2017 = zScale["electrons"][0]			
				zScaleFac2018 = zScale2018["electrons"][0]			

	else:
		zScaleFac = zScale["muons"]
		if not plot.muon:
			if "bbbe" in plot.histName:
				zScaleFac = zScale["electrons"][0]
			elif "bb" in plot.histName:
				zScaleFac = zScale["electrons"][1]
			elif "be" in plot.histName:
				zScaleFac = zScale["electrons"][2]
			else:
				zScaleFac = zScale["electrons"][0]			
	CIBins = args.CIBins		
	# Data and background loading	
	if plot.plot2D:	
		if args.useRun2:
			
			datahist = data2016.loadHistogramProjected(plot,lumi2016,zScaleFac2016)	
			datahist.Add(data2017.loadHistogramProjected(plot,lumi2017,zScaleFac2017))
			datahist.Add(data2018.loadHistogramProjected(plot,lumi2018,zScaleFac2018))
			
			stack = TheStack2DRun2([processes2016,processes2017,processes2018],[lumi2016,lumi2017,lumi2018],plot,[zScaleFac2016,zScaleFac2017,zScaleFac2018])

		else:
			datahist = data.loadHistogramProjected(plot,lumi,zScaleFac)	

			stack = TheStack2D(processes,lumi,plot,zScaleFac)
	else:
		
		if args.useRun2:
			datahist = data2016.loadHistogram(plot,lumi2016,zScaleFac2016,CIBins=args.CIBins)	
			datahist.Add(data2017.loadHistogram(plot,lumi2017,zScaleFac2017,CIBins=args.CIBins))
			datahist.Add(data2018.loadHistogram(plot,lumi2018,zScaleFac2018,CIBins=args.CIBins))	
			stack = TheStackRun2([processes2016,processes2017,processes2018],[lumi2016,lumi2017,lumi2018],plot,[zScaleFac2016,zScaleFac2017,zScaleFac2018],CIBins=args.CIBins)
		else:
			datahist = data.loadHistogram(plot,lumi,zScaleFac,CIBins=args.CIBins)	
			stack = TheStack(processes,lumi,plot,zScaleFac,CIBins=CIBins)
	if args.data:
		yMax = datahist.GetBinContent(datahist.GetMaximumBin())
		if "Mass" in plot.fileName:
			yMin = 0.00001
		else:
			yMin = 0.01
		xMax = datahist.GetXaxis().GetXmax()
		xMin = datahist.GetXaxis().GetXmin()
	else:	
		yMax = stack.theHistogram.GetBinContent(datahist.GetMaximumBin())
		yMin = 0.01
		xMax = stack.theHistogram.GetXaxis().GetXmax()
		xMin = stack.theHistogram.GetXaxis().GetXmin()
		

	print (datahist.Integral(datahist.FindBin(60),datahist.FindBin(120-0.001)), datahist.Integral(datahist.FindBin(60),datahist.FindBin(120-0.001))**0.5, datahist.Integral(datahist.FindBin(60),datahist.FindBin(120-0.001))**0.5/datahist.Integral(datahist.FindBin(60),datahist.FindBin(120-0.001)))

			
	if plot.yMax == None:
		if logScale:
			yMax = yMax*10000
		else:
			yMax = yMax*1.5
	else: yMax = plot.yMax
	
	if "Mass" in plot.fileName:
		yMax = 20000000	
	
	if not plot.yMin == None:
		yMin = plot.yMin
	if not plot.xMin == None:
		xMin = plot.xMin
	if not plot.xMax == None:
		xMax = plot.xMax
		
	if args.CIBins:
		xMin = 400
		xMax = 6000			
		
	#if args.ADD and args.use2016: 
	#	xMin = 1700
	#	xMax = 4000
	#	yMax = 1.0
	if "CosThetaStarBBM1800" in plot.fileName:
		yMax = 3
	plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %(plot.xaxis,plot.yaxis))

	
	drawStack = stack

					

	# Draw signal information
	if len(args.signals) != 0:
		if args.useRun2:
			signalhists = []
			for index, Signal in enumerate(signals2016):
				if plot.plot2D: # plot collins-soper angle
					signalhist = Signal.loadHistogramProjected(plot,lumi2016, zScaleFac2016)
					signalhist.Add(signals2017[index].loadHistogramProjected(plot,lumi2017, zScaleFac2017))
					signalhist.Add(signals2018[index].loadHistogramProjected(plot,lumi2018, zScaleFac2018))
					signalhist.SetLineWidth(2)
					signalBackgrounds = deepcopy(backgrounds)
					signalBackgrounds.remove("DrellYan")
					signalProcesses2016 = []
					signalProcesses2017 = []
					signalProcesses2018 = []
					for background in signalBackgrounds:
						if background == "Jets":
							signalProcesses2016.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
							signalProcesses2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
							signalProcesses2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights,normalized=True))
						elif background == "Other" and not plot.muon:
							signalProcesses2016.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts,negWeights))
							signalProcesses2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
							signalProcesses2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))							
						else:	
							signalProcesses2016.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
							signalProcesses2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
							signalProcesses2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))
					signalStack = TheStack2DRun2([signalProcesses2016,signalProcesses2017,signalProcesses2018],[lumi2016,lumi2017,lumi2018],plot, [zScaleFac2016,zScaleFac2017,zScaleFac2018],CIBins=args.CIBins)
					signalhist.Add(signalStack.theHistogram)
					signalhist.SetMinimum(0.1)
					signalhist.Draw("samehist")
					signalhists.append(signalhist)	
				else:
					signalhist = Signal.loadHistogram(plot,lumi2016, zScaleFac2016)
					signalhist.Add(signals2017[index].loadHistogram(plot,lumi2017, zScaleFac2017))
					signalhist.Add(signals2018[index].loadHistogram(plot,lumi2018, zScaleFac2018))
					signalhist.SetLineWidth(2)
					signalBackgrounds = deepcopy(backgrounds)
					signalBackgrounds.remove("DrellYan") # signalBackgrounds = ["Jets", "Other"]
					signalProcesses2016 = []
					signalProcesses2017 = []
					signalProcesses2018 = []
					for background in signalBackgrounds:
						if background == "Jets":
							signalProcesses2016.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
							signalProcesses2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
							if plot.muon:
								signalProcesses2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights,normalized=True))
							else:
								signalProcesses2018.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
						else:	
							if background == "Other" and not plot.muon:
								signalProcesses2016.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts,negWeights))
							else:	
								signalProcesses2016.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
							signalProcesses2017.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
							signalProcesses2018.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))
					signalStack = TheStackRun2([signalProcesses2016,signalProcesses2017,signalProcesses2018],[lumi2016,lumi2017,lumi2018],plot, [zScaleFac2016,zScaleFac2017,zScaleFac2018],CIBins=args.CIBins)
					signalhist.Add(signalStack.theHistogram)
					signalhist.SetMinimum(0.0001)
					signalhist.Draw("samehist")
					signalhists.append(signalhist)	
		else:
			signalhists = []
			for Signal in signals:
				if plot.plot2D: # plot collins-soper angle
					signalhist = Signal.loadHistogramProjected(plot,lumi, zScaleFac)
					signalhist.SetLineWidth(2)
					signalBackgrounds = deepcopy(backgrounds)
					signalBackgrounds.remove("DrellYan")
					signalProcesses = []
					for background in signalBackgrounds:
						if args.use2016:
							if background == "Jets":
								processes.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
							else:	
								processes.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
						elif args.use2018:
							if background == "Jets":
								processes.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights,normalized=True))
							else:	
								processes.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))						
						else:		
							if background == "Jets":
								signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
							else:	
								signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
							
					signalStack = TheStack2D(signalProcesses,lumi,plot, zScaleFac,CIBins=args.CIBins)
					signalhist.Add(signalStack.theHistogram)
					signalhist.SetMinimum(0.1)
					signalhist.Draw("samehist")
					signalhists.append(signalhist)	
				else:
					signalhist = Signal.loadHistogram(plot,lumi,zScaleFac)
					signalhist.SetLineWidth(2)
					signalBackgrounds = deepcopy(backgrounds)
					signalBackgrounds.remove("DrellYan") # signalBackgrounds = ["Jets", "Other"]
					signalProcesses = []
					for background in signalBackgrounds:
						if args.use2016:
							if background == "Jets":
								signalProcesses.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights,normalized=True))
							else:	
								signalProcesses.append(Process(getattr(Backgrounds2016,background),eventCounts,negWeights))
						elif args.use2018:
							if background == "Jets":
								signalProcesses.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights,normalized=True))
							else:	
								signalProcesses.append(Process(getattr(Backgrounds2018,background),eventCounts,negWeights))	
						else:							
							if background == "Jets":
								signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
							else:	
								signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
					signalStack = TheStack(signalProcesses,lumi,plot,zScaleFac,CIBins=args.CIBins)
					signalhist.Add(signalStack.theHistogram)
					signalhist.SetMinimum(0.0001)
					signalhist.Draw("samehist")
					signalhists.append(signalhist)	

	

	# Draw background from stack
	drawStack.theStack.Draw("samehist")	

	# Draw data
	datahist.SetMinimum(0.0001)

	if "Mass" in plot.fileName:
		histForIntegral = datahist.Clone("blubb")
		histForIntegral2 = drawStack.theHistogram.Clone("blubb")
		for i in range(0,histForIntegral.GetNbinsX()):
			histForIntegral.SetBinContent(i,histForIntegral.GetBinContent(i)*histForIntegral.GetBinWidth(i))		
			histForIntegral2.SetBinContent(i,histForIntegral2.GetBinContent(i)*histForIntegral2.GetBinWidth(i))		

	if args.data:
		datahist.Draw("samep")	

	# Draw legend
	if "Eta" in plot.fileName or "CosTheta" in plot.fileName:
		legendEta.Draw()
	else:
		legend.Draw()

	plotPad.SetLogx(plot.logX)
	if args.useRun2:
		if plot.muon:
			latex.DrawLatex(0.95, 0.96, "139.7 fb^{-1} (13 TeV)")
		else:	
			latex.DrawLatex(0.95, 0.96, "136.8 fb^{-1} (13 TeV)")
	else:	
		latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV)"%(float(lumi)/1000,))
	yLabelPos = 0.85
	cmsExtra = "Preliminary"
	if not args.data:
		cmsExtra = "#splitline{Preliminary}{Simulation}"
		yLabelPos = 0.82	
	latexCMS.DrawLatex(0.19,0.89,"CMS")
	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
	if args.ratio:
		try:
			ratioPad.cd()
			ratioPad.SetLogx(plot.logX)
		except AttributeError:
			print ("Plot fails. Look up in errs/failedPlots.txt")
			outFile =open("errs/failedPlots.txt","a")
			outFile.write('%s\n'%plot.filename%("_"+run.label+"_"+dilepton))
			outFile.close()
			plot.cuts=baseCut
			return 1
		ratioGraphs =  ratios.RatioGraph(datahist,drawStack.theHistogram, xMin=xMin, xMax=xMax,title="(Data - Bkg) / Bkg",titleX="",yMin=-1.0,yMax=1.0,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000000000000,labelSize=0.125,pull=True)
		ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
					

	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	if args.ratio:

		ratioPad.RedrawAxis()
	if not os.path.exists("plots"):
		os.makedirs("plots")	
	if args.use2016:
		hCanvas.Print("plots/"+plot.fileName+"_2016.pdf")
	elif args.use2018:
		hCanvas.Print("plots/"+plot.fileName+"_2018.pdf")
	elif args.useRun2:
		hCanvas.Print("plots/"+plot.fileName+"_Run2.pdf")
	else:	
		hCanvas.Print("plots/"+plot.fileName+"_2017.pdf")

					
if __name__ == "__main__":
	
	
	parser = argparse.ArgumentParser(description='Process some integers.')
	
	parser.add_argument("-d", "--data", action="store_true", dest="data", default=False,
						  help="plot data points.")
	parser.add_argument("-m", "--mc", action="store_true", dest="mc", default=False,
						  help="plot mc backgrounds.")
	parser.add_argument("-p", "--plot", dest="plot", nargs=1, default="",
						  help="plot to plot.")
	parser.add_argument("-n", "--norm", action="store_true", dest="norm", default=False,
						  help="normalize to data.")
	parser.add_argument("-2016", "--2016", action="store_true", dest="use2016", default=False,
						  help="use 2016 data and MC.")
	parser.add_argument("-2018", "--2018", action="store_true", dest="use2018", default=False,
						  help="use 2018 data with 2017 MC.")
	parser.add_argument("-Run2", "--Run2", action="store_true", dest="useRun2", default=False,
						  help="use 2018 data with 2017 MC.")
	parser.add_argument("-r", "--ratio", action="store_true", dest="ratio", default=False,
						  help="plot ratio plot")
	parser.add_argument("-l", "--log", action="store_true", dest="log", default=False,
						  help="plot with log scale for y axis")
	parser.add_argument("-s", "--signal", dest="signals", action="append", default=[],
						  help="signals to plot.")
	parser.add_argument("-b", "--backgrounds", dest="backgrounds", action="append", default=[],
						  help="backgrounds to plot.")
	parser.add_argument("-a", "--ADD", action="store_true", dest="ADD", default=False, help="plot add signals")
	parser.add_argument("--CIBins", action="store_true", dest="CIBins", default=False, help="use CI limit binning for mass bins")


	args = parser.parse_args()
	if len(args.backgrounds) == 0:
		args.backgrounds = ["Wjets","Other","DrellYan"]
		#~ args.backgrounds = ["Diboson","DrellYan"]

	if len(args.signals) != 0:
		args.plotSignal = True

	if args.plot == "":
		args.plot = plotList
	
	signals = args.signals
	for plot in args.plot:
		args.signals = signals
		plotObject = getPlot(plot)
		if len(args.signals) > 0:
			#~ if ("To2E" in args.signals[0] and plotObject.muon) or ("To2Mu" in args.signals[0] and not plotObject.muon):
			args.signals = []
			if plotObject.muon:
				for signal in signals:
					if args.ADD: args.signals.append("ADDGravTo2Mu_"+signal)
					else: args.signals.append("CITo2Mu_"+signal)
			else:
				for signal in signals:
					if args.ADD: args.signals.append("ADDGravTo2E_"+signal)
					else: args.signals.append("CITo2E_"+signal)
		#~ print args.plotSignal	
		plotDataMC(args,plotObject)
	
