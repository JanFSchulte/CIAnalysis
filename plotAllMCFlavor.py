import argparse	
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT, TGaxis
import ratios
from setTDRStyle import setTDRStyle
gROOT.SetBatch(True)
from helpers import *
from defs import getPlot, Backgrounds, Backgrounds2016, Backgrounds2018, Signals, Signals2016, Signals2016ADD, Data, Data2016, Data2018, path, plotList, zScale, zScale2016, zScale2018
import math
import os
from copy import copy
import numpy as np
import root_numpy
from copy import deepcopy
# Muon sys uncertainty (%) 
# as a function of mass
def znorm(hist):
	lowLimit=hist.FindBin(60)
	upLimit=hist.FindBin(120)
	znorm=0
	for i in range(lowLimit+1,upLimit+1):
		znorm+=hist.GetBinContent(i)
	return znorm

def scale(hist1,hist2):
	lowLimit=hist1.FindBin(200)
	upLimit=hist1.FindBin(400)
	sFac1=hist1.Integral(lowLimit,upLimit)
	sFac2=hist2.Integral(lowLimit,upLimit)
	print(sFac2)
	return sFac1/sFac2
def getMuErr(mass, chann, norm=False):
	lumi = 0.0
	znorm = 0.0
	pileup = 0.0   # we don't use pileup 0.046
	dybkg = 0.0    # 0.07
	#pdf = 0.01*(0.433+0.003291*mass-2.159e-6*mass**2+9.044e-10*mass**3-1.807e-13*mass**4+1.51e-17*mass**5)
	pdf = 0.0
	
	# muons only next
	if chann: mutrig = 0.003
	else: mutrig = 0.007
	resolution = 0.01
	muid = 0.05

	if norm: 
		lumi = 0.0
		znorm = 0.0
		dybkg = 0
	return math.sqrt(lumi**2+znorm**2+pileup**2+dybkg**2+pdf**2+mutrig**2+resolution**2+muid**2)


# chann = True if BB
# chann = False if BE
def getElErr(mass, chann, norm=False):

	lumi = 0.0
	znorm = 0.0
	pileup = 0.046  # we don't use pileup 0.046
	dybkg = 0.0   # 0.07
	
	# poly values are in %
	#pdf = 0.01*(0.433 + 0.003291*mass - 2.159e-6*mass**2 + 9.044e-10*mass**3 - 1.807e-13*mass**4 + 1.51e-17*mass**5)
	pdf = 0.0
	
	# the following two are electrons only
	if chann: energyscale = 0.02
	else: energyscale = 0.01
	
	if chann: 
		if mass < 90: idscale = 0.01
		elif mass < 1000: idscale = 0.00002198 * mass + 0.008
		else: idscale = 0.03
	else:
		if mass < 90: idscale = 0.01
		elif mass < 300: idscale = 0.00014286 * mass - 0.00285
		else: idscale = 0.04
	
	if chann: scalefac = 0.03
	else: scalefac = 0.05
	
	if norm:
		lumi = 0.0
		znorm = 0.0
		dybkg = 0
	return math.sqrt(lumi**2+znorm**2+ pileup**2 + dybkg**2 + pdf**2 + energyscale**2 + idscale**2 + scalefac**2)

def getErrors(default, others, keys):
	dfarr=root_numpy.hist2array(default)
	errs={}
	for other, key in zip(others,keys):
		if type(other)==list:
			err1=root_numpy.hist2array(other[0])-dfarr
			err1=abs(err1)
			err2=root_numpy.hist2array(other[1])-dfarr
			err2=abs(err2)
			err=np.maximum(err1,err2)
			errs[key]=err
		else:           
			err=root_numpy.hist2array(other)-dfarr
			err=abs(err)
			errs[key]=err
	return errs                           

# multiply hist by 1/(Acceptance x Efficiency)
def Rebin(hist):
	bng=[50, 120,150,200,300,400,500,690,900,1250,1610, 3490, 6070]
	#bng=([j for j in range(50, 120, 5)] +
        #                        [j for j in range(120, 150, 5)] +
        #                        [j for j in range(150, 200, 10)] +
        #                        [j for j in range(200, 600, 20)]+
        #                        [j for j in range(600, 900, 30) ]+
        #                        [j for j in range(900, 1250, 50)] +
        #                        [j for j in range(1250, 1610, 60) ] +
        #                        [j for j in range(1610, 1890, 70) ] +
        #                        [j for j in range(1890, 3970, 80) ] +
        #                        [j for j in range(3970, 6070, 100) ] +
        #                        [6070])
	#bng=[50,200,1500,6070]
	hist.Sumw2()
	for i in range(0,hist.GetNbinsX()):
		hist.SetBinContent(i,hist.GetBinContent(i)*hist.GetBinWidth(i))
		hist.SetBinError(i,hist.GetBinError(i)*hist.GetBinWidth(i))
	hist=hist.Rebin(len(bng) - 1, 'hist_' + uuid.uuid4().hex, array('d', bng))
	for i in range(0,hist.GetNbinsX()):
		hist.SetBinContent(i,hist.GetBinContent(i)/hist.GetBinWidth(i))
		hist.SetBinError(i,hist.GetBinError(i)/hist.GetBinWidth(i))
	return deepcopy(hist)

def inverseAE(hist, plotObj, year):
		# muon and electron
		# BB and BE
	
	if year == 2016:
		if plotObj.muon:
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 600:
						ae = 2.129-0.1268*math.exp(-(mass-119.2)/22.35)-2.386*mass**(-0.03619)
					else:
						ae = 2.891-2.291e+04/(mass+8294.)-0.0001247*mass
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
					hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
					#print(mass)
					#print("2016 BB mu ae = %f" %(ae))
			elif "BE" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 450:
									ae = 13.56-6.672*math.exp((mass+4.879e+06)/7.233e+06)-826*mass**(-1.567)
					else:
									ae =  0.2529+0.06511*mass**0.8755*math.exp(-(mass+4601.)/1147)
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
					hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
					#print(mass)
					#print("2016 BE mu ae = %f" %(ae))
		else: # is electron
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					ae = 0.6386-497.7/(mass+348.7) + 69570.0/(mass**2+115400.0)
					#print mass, ae
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
					hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
					#print(mass)
					#print("2016 BB e ae = %f" %(ae))
					
			elif "BE" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					ae = -0.03377+735.1/(mass+1318)-88890.0/(mass**2+75720)+14240000.0/(mass**3+23420000)
					#print mass, ae
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
					hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
					#print(mass)
					#print("2016 BE e ae = %f" %(ae))

	elif year == 2017:
		if plotObj.muon:
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 600:
						ae = 2.13-0.1313*math.exp(-(mass-110.9)/20.31)-2.387*mass**(-0.03616)
					else:
						ae = 4.931-55500.0/(mass+11570.0)-0.0002108*mass
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
					hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
					#print(mass)
					#print("2017 BB mu ae = %f" %(ae))
			elif "BE" in plotObj.fileName:
					for i in range(1, hist.GetSize()-1):
						mass = hist.GetBinCenter(i)
						if mass < 450:
							ae = 13.39-6.696*math.exp((mass+4855000.0)/7431000.0)-108.8*mass**(-1.138)
						else:
							ae = 0.3148+0.04447*mass**1.42*math.exp(-(mass+5108.0)/713.5)
						#print mass, ae
						if mass < 120: ae = float("inf")
						hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
						hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
						#print(mass)
						#print("2017 BE mu ae = %f" %(ae))
		else: # is electron
				if "BB" in plotObj.fileName:
						for i in range(1, hist.GetSize()-1):
								mass = hist.GetBinCenter(i)
								ae = 0.585-404.6/(mass+279.5) + 56180.0/(mass**2+91430.0)
								#print mass, ae
								hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
								hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
								#print(mass)
								#print("2017 BB e ae = %f" %(ae))
				elif "BE" in plotObj.fileName:
						for i in range(1, hist.GetSize()-1):
								mass = hist.GetBinCenter(i)
								ae = 0.02066+429.7/(mass+729)-90960.0/(mass**2+71900)+13780000.0/(mass**3+22050000)
								#print mass, ae
								hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
								hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
								#print(mass)
								#print("2017 BE e ae = %f" %(ae))
	elif year == 2018:
			if plotObj.muon:
					if "BB" in plotObj.fileName:
							for i in range(1, hist.GetSize()-1):
									mass = hist.GetBinCenter(i)
									if mass < 600:
											ae = 2.14-0.1286*math.exp(-(mass-110.6)/22.44)-2.366*mass**(-0.03382)
									else:
											ae = 5.18-58450.0/(mass+11570.0)-0.0002255*mass
									#print mass, ae
									if mass < 120: ae = float("inf")
									hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
									hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
									#print(mass)
									#print("2018 BB mu ae = %f" %(ae))
					elif "BE" in plotObj.fileName:
							for i in range(1, hist.GetSize()-1):
									mass = hist.GetBinCenter(i)
									if mass < 450:
											ae = 13.4-6.693*math.exp((mass+4852000.0)/7437000.0)-81.43*mass**(-1.068)
									else:
											ae = 0.3154+0.04561*mass**1.362*math.exp(-(mass+4927.0)/727.5)
									#print mass, ae
									if mass < 120: ae = float("inf")
									hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
									hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
									#print(mass)
									#print("2018 BE mu ae = %f" %(ae))
			else: # is electron
					if "BB" in plotObj.fileName:
							for i in range(1, hist.GetSize()-1):
									mass = hist.GetBinCenter(i)
									ae = 0.576-417.7/(mass+381.8) + 46070.0/(mass**2+107200)
									#print mass, ae
									hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
									hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
									#print(mass)
									#print("2018 BB e ae = %f" %(ae))
					elif "BE" in plotObj.fileName:
							for i in range(1, hist.GetSize()-1):
									mass = hist.GetBinCenter(i)
									ae = 0.01443+475.7/(mass+639.1)-105600.0/(mass**2+82810)+12890000.0/(mass**3+23170000)
									#print mass, ae
									hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
									hist.SetBinError(i, hist.GetBinError(i)*1.0/ae)
									#print(mass)
									#print("2018 BE e ae = %f" %(ae))																			
def Stacks(processes,lumi,plot,zScale):
	stacks=[]
	for i in range(3):
		stacks.append(TheStack(processes[i],lumi[i],plot,zScale[i]))
	return stacks
def Addhist(histlist):
	tempHist=histlist[0]
	for i in range(1,3):
		tempHist.Add(histlist[i])
	return tempHist         
def Addstack(Stacklist):
	tempStack=Stacklist[0]
	for i in range(1,3):
		tempStack.Add(Stacklist[i])
	return tempStack                                                                                                                                                                          
def plotDataMC(args,plot_mu,plot_el, category):
	

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	if args.ratio:
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.5,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.5)
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
		data_mu = Process(Data2016, normalized=True)
		data_el = Process(Data2016, normalized=True)
	elif args.use2018:
		data_mu = Process(Data2018, normalized=True)
		data_el = Process(Data2018, normalized=True)
	elif args.useall:
		data_all=[Process(Data2016, normalized=True),Process(Data, normalized=True),Process(Data2018, normalized=True)]              
	else:	
		data_mu = Process(Data, normalized=True)
		data_el = Process(Data, normalized=True)
	
	eventCounts_mu = totalNumberOfGeneratedEvents(path,plot_mu["default"].muon)	
	eventCounts_el = totalNumberOfGeneratedEvents(path,plot_el["default"].muon)
	negWeights_mu = negWeightFractions(path,plot_mu["default"].muon)
	negWeights_el = negWeightFractions(path,plot_el["default"].muon)

	# Background load processes	
	backgrounds = copy(args.backgrounds)
	if plot_mu["default"].useJets:
		if "Wjets" in backgrounds:
			backgrounds.remove("Wjets")
			backgrounds.insert(0,"Jets")
	processes_mu = []
	processes_el = []
	processes_mu2016=[]
	processes_mu2017=[]
	processes_mu2018=[]
	processes_el2016=[]
	processes_el2017=[]
	processes_el2018=[]
	for background in backgrounds:
		if args.use2016:
			if background == "Jets":
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu,normalized=True))
				processes_el.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el,normalized=True))
			elif background == "Other":
				processes_el.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts_el,negWeights_el))
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
			else:	
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
				processes_el.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el))
		elif args.use2018:
				if background == "Jets":
						processes_mu.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el,normalized=True))
				else:
						processes_mu.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu))
						processes_el.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el))
		elif args.useall:
				if background == "Jets":
						processes_mu2016.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el2016.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el,normalized=True))
						processes_mu2017.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el2017.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el,normalized=True))
						processes_mu2018.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el2018.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el,normalized=True))
						processes_mu=[processes_mu2016,processes_mu2017,processes_mu2018]
						processes_el=[processes_mu2016,processes_mu2017,processes_mu2018]
				elif background == "Other":
						processes_el2016.append(Process(getattr(Backgrounds2016,"OtherEle"),eventCounts_el,negWeights_el))
                                		processes_mu2016.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
						processes_mu2017.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu))
                                        	processes_el2017.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el))
                                        	processes_mu2018.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu))
                                        	processes_el2018.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el))
                                        	processes_mu=[processes_mu2016,processes_mu2017,processes_mu2018]
                                        	processes_el=[processes_el2016,processes_el2017,processes_el2018]

				else:
						processes_mu2016.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
						processes_el2016.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el))
						processes_mu2017.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu))
						processes_el2017.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el))
						processes_mu2018.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu))
						processes_el2018.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el))
						processes_mu=[processes_mu2016,processes_mu2017,processes_mu2018]
						processes_el=[processes_el2016,processes_el2017,processes_el2018]
		else:             
				if background == "Jets":
						processes_mu.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el,normalized=True))
				else:
						processes_mu.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu))
						processes_el.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el))

	
	legend = TLegend(0.55, 0.75, 0.925, 0.925)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextFont(42)
	
	
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
	legendHistData_mu = ROOT.TH1F()
	legendHistData_el = ROOT.TH1F()
	dy_mu = ROOT.TH1F()
	dy_el = ROOT.TH1F()
	if args.data:	
		legendHistData_mu.SetMarkerColor(ROOT.kViolet)
		legendHistData_el.SetMarkerColor(ROOT.kOrange)
		dy_mu.SetLineColor(ROOT.kBlue-3)
		dy_el.SetLineColor(ROOT.kRed-3)
		legend.AddEntry(legendHistData_mu,"Data #rightarrow #mu^{+}#mu^{-}","pe")
		legend.AddEntry(legendHistData_el,"Data #rightarrow e^{+}e^{-}", "pe")
		legend.AddEntry(dy_mu, "MC Inclusive #rightarrow #mu^{+}#mu^{-}", "l")
		legend.AddEntry(dy_el, "MC Inclusive #rightarrow e^{+}e^{-}", "l")
		 	
	if args.useall:
		for i in range(3):
			for process in reversed(processes_mu[i]):
				if not plot_mu["default"].muon and "#mu^{+}#mu^{-}" in process.label:
					process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
				process.theColor = ROOT.kBlue
				process.theLineColor = ROOT.kBlue
				temphist = ROOT.TH1F()
				temphist.SetFillColor(process.theColor)

			for process in reversed(processes_el[i]):
				if not plot_el["default"].muon and "#mu^{+}#mu^{-}" in process.label:
					process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
				process.theColor = ROOT.kRed
				process.theLineColor = ROOT.kRed
				temphist = ROOT.TH1F()
				temphist.SetFillColor(process.theColor)
	else:
		for process in reversed(processes_mu):
			if not plot_mu["default"].muon and "#mu^{+}#mu^{-}" in process.label:
				process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
			process.theColor = ROOT.kBlue
			process.theLineColor = ROOT.kBlue
			temphist = ROOT.TH1F()
			temphist.SetFillColor(process.theColor)
	
		for process in reversed(processes_el):
			if not plot_el["default"].muon and "#mu^{+}#mu^{-}" in process.label:
				process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
			process.theColor = ROOT.kRed
			process.theLineColor = ROOT.kRed
			temphist = ROOT.TH1F()
			temphist.SetFillColor(process.theColor)

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
	logScale = plot_mu["default"].log
	
	if logScale == True:
		plotPad.SetLogy()

	if args.use2016:	
		lumi_el = 35.9*1000
		lumi_mu = 36.3*1000
	elif args.use2018:	
		lumi_el = 59.97*1000
		lumi_mu = 61.608*1000
	elif args.useall:
		lumi_el = [35.9*1000,41.529*1000,59.97*1000]
		lumi_mu = [36.3*1000,42.135*1000,61.608*1000]
	else:
		lumi_el = 41.529*1000
		lumi_mu = 42.135*1000

	if args.use2016:		
		zScaleFac_mu = zScale2016["muons"]
		if "BE" in plot_el["default"].fileName:
			zScaleFac_el = zScale2016["electrons"][0]
		elif "BB" in plot_el["default"].fileName:
			zScaleFac_el = zScale2016["electrons"][1]
		elif "BE" in plot_el["default"].fileName:
			zScaleFac_el = zScale2016["electrons"][2]
		else:
                        zScaleFac_el = zScale2016["electrons"][0]
	elif args.use2018:		
		zScaleFac_mu = zScale2018["muons"]
		if "BBBE" in plot_el["default"].fileName:
			zScaleFac_el = zScale2018["electrons"][0]
		elif "BB" in plot_el["default"].fileName:
			zScaleFac_el = zScale2018["electrons"][1]
		elif "BE" in plot_el["default"].fileName:
			zScaleFac_el = zScale2018["electrons"][2]
		else:
			zScaleFac_el = zScale2018["electrons"][0]

	elif args.useall:
		zScaleFac_mu = [zScale2016["muons"],zScale["muons"],zScale2018["muons"]]
		if "BBBE" in plot_el["default"].fileName:
			zScaleFac_el = [zScale2016["electrons"][0],zScale["electrons"][0],zScale2018["electrons"][0]]
		elif "BB" in plot_el["default"].fileName:
			zScaleFac_el = [zScale2016["electrons"][1],zScale["electrons"][1],zScale2018["electrons"][1]]
		elif "BE" in plot_el["default"].fileName:
			zScaleFac_el = [zScale2016["electrons"][2],zScale["electrons"][2],zScale2018["electrons"][2]]
		else:
			zScaleFac_el = [zScale2016["electrons"][0],zScale["electrons"][0],zScale2018["electrons"][0]]
	else:	
		zScaleFac_mu = zScale["muons"]
		if "BBBE" in plot_el["default"].fileName:
			zScaleFac_el = zScale["electrons"][0]
		elif "BB" in plot_el["default"].fileName:
			zScaleFac_el = zScale["electrons"][1]
		elif "BE" in plot_el["default"].fileName:
			zScaleFac_el = zScale["electrons"][2]
		else:
			zScaleFac_el = zScale["electrons"][0]

	
			
	# Data and background loading
	if args.useall:
		datamu=[]
		datael=[]
		for i in range(3):
			datamu.append(data_all[i].loadHistogram(plot_mu["default"],lumi_mu[i],zScaleFac_mu[i]))
			datael.append(data_all[i].loadHistogram(plot_el["default"],lumi_el[i],zScaleFac_el[i]))
		stackmu = Stacks(processes_mu,lumi_mu,plot_mu["default"],zScaleFac_mu)
		mu_scaleup=Stacks(processes_mu,lumi_mu,plot_mu["scale_up"],zScaleFac_mu)
		mu_scaledown=Stacks(processes_mu,lumi_mu,plot_mu["scale_down"],zScaleFac_mu)
		mu_ID=Stacks(processes_mu,lumi_mu,plot_mu["ID"],zScaleFac_mu)
		mu_reso=Stacks(processes_mu,lumi_mu,plot_mu["reso"],zScaleFac_mu)
		stackel = Stacks(processes_el,lumi_el,plot_el["default"],zScaleFac_el)
		el_scaleup=Stacks(processes_el,lumi_el,plot_el["scale_up"],zScaleFac_el)
		el_scaledown=Stacks(processes_el,lumi_el,plot_el["scale_down"],zScaleFac_el)
		el_PUup=Stacks(processes_el,lumi_el,plot_el["PU_up"],zScaleFac_el)
		el_PUdown=Stacks(processes_el,lumi_el,plot_el["PU_down"],zScaleFac_el)
	else:	
		datamu = data_mu.loadHistogram(plot_mu["default"],lumi_mu,zScaleFac_mu)
		datael = data_el.loadHistogram(plot_el["default"],lumi_el,zScaleFac_el)
		stackmu = TheStack(processes_mu,lumi_mu,plot_mu["default"],zScaleFac_mu)
		mu_scaleup=TheStack(processes_mu,lumi_mu,plot_mu["scale_up"],zScaleFac_mu)
		mu_scaledown=TheStack(processes_mu,lumi_mu,plot_mu["scale_down"],zScaleFac_mu)
		mu_ID=TheStack(processes_mu,lumi_mu,plot_mu["ID"],zScaleFac_mu)
		mu_reso=TheStack(processes_mu,lumi_mu,plot_mu["reso"],zScaleFac_mu)
			
		stackel = TheStack(processes_el,lumi_el,plot_el["default"],zScaleFac_el)
		el_scaleup=TheStack(processes_el,lumi_el,plot_el["scale_up"],zScaleFac_el)
		el_scaledown=TheStack(processes_el,lumi_el,plot_el["scale_down"],zScaleFac_el)
		el_PUup=TheStack(processes_el,lumi_el,plot_el["PU_up"],zScaleFac_el)
		el_PUdown=TheStack(processes_el,lumi_el,plot_el["PU_down"],zScaleFac_el)
	
	if args.alter:
		if args.useall:
			for i in range(3):
				sFac=scale(datael[i],datamu[i])	
				datamu[i].Scale(sFac)
				sFac=scale(stackel[i].theHistogram,stackmu[i].theHistogram)
				stackmu[i].theHistogram.Scale(sFac)
				mu_scaleup[i].theHistogram.Scale(sFac)
				mu_scaledown[i].theHistogram.Scale(sFac)
				mu_ID[i].theHistogram.Scale(sFac)
				mu_reso[i].theHistogram.Scale(sFac)
			func=ROOT.TF1("func","[0]+[1]*log(x)+[2]*x^2",200,3500)
			
			i=0
                	Errs_mu_fit=[]
                	Errs_el_fit=[]
                	for year in range(2016,2019):
                        	lis_mu_fit=[[mu_scaleup[i].theHistogram,mu_scaledown[i].theHistogram],mu_ID[i].theHistogram,mu_reso[i].theHistogram]
                        	key_mu_fit=["massScale","ID","resolution"]
                        	lis_el_fit=[[el_scaleup[i].theHistogram,el_scaledown[i].theHistogram],[el_PUup[i].theHistogram,el_PUdown[i].theHistogram]]
                        	key_el_fit=["massScale","pileUp"]
                        	Errs_mu_fit.append(getErrors(stackmu[i].theHistogram,lis_mu_fit,key_mu_fit))
                        	Errs_el_fit.append(getErrors(stackel[i].theHistogram,lis_el_fit,key_el_fit))
                        	i+=1
                	if category=="bb":
                        	errmu_fit=Errs_mu_fit[0]["massScale"]**2+(Errs_mu_fit[1]["massScale"]+Errs_mu_fit[2]["massScale"])**2+Errs_mu_fit[0]["resolution"]**2+Errs_mu_fit[1]["resolution"]**2+Errs_mu_fit[2]["resolution"]**2+Errs_mu_fit[0]["ID"]**2+(Errs_mu_fit[1]["ID"]+Errs_mu_fit[2]["ID"])**2
                	else:
                        	errmu_fit=Errs_mu_fit[0]["massScale"]**2+(Errs_mu_fit[1]["massScale"]+Errs_mu_fit[2]["massScale"])**2+Errs_mu_fit[0]["resolution"]**2+Errs_mu_fit[1]["resolution"]**2+Errs_mu_fit[2]["resolution"]**2+Errs_mu_fit[0]["ID"]**2+Errs_mu_fit[1]["ID"]**2+Errs_mu_fit[2]["ID"]**2
                	errel_fit=Errs_el_fit[0]["massScale"]**2+Errs_el_fit[0]["pileUp"]**2+Errs_el_fit[1]["massScale"]**2+Errs_el_fit[1]["pileUp"]**2+Errs_el_fit[2]["massScale"]**2+Errs_el_fit[2]["pileUp"]**2
                	scale_el_fit=Errs_el_fit[0]["massScale"]**2+Errs_el_fit[1]["massScale"]**2+Errs_el_fit[2]["massScale"]**2
                	scale_mu_fit=Errs_mu_fit[0]["massScale"]**2+(Errs_mu_fit[1]["massScale"]+Errs_mu_fit[2]["massScale"])**2
                	stackmu_fit=Addstack(stackmu)
                	stackel_fit=Addstack(stackel)
                	mu_scaleup_fit=Addstack(mu_scaleup)
                	mu_scaledown_fit=Addstack(mu_scaledown)
                	mu_ID_fit=Addstack(mu_ID)
                	mu_reso_fit=Addstack(mu_reso)
                	el_scaleup_fit=Addstack(el_scaleup)
                	el_scaledown_fit=Addstack(el_scaledown)
                	el_PUup_fit=Addstack(el_PUup)
                	el_PUdown_fit=Addstack(el_PUdown)
                	datamu_fit=Addhist(datamu)
                	datael_fit=Addhist(datael)
       			if args.ratio:
				hhmu_fit = stackmu_fit.theHistogram
                		hhel_fit = stackel_fit.theHistogram
                		h_mu_scaleup_fit=mu_scaleup_fit.theHistogram
                		h_mu_scaledown_fit=mu_scaledown_fit.theHistogram
                		h_mu_ID_fit=mu_ID_fit.theHistogram
                		h_mu_reso_fit=mu_reso_fit.theHistogram
                		h_el_scaleup_fit=el_scaleup_fit.theHistogram
                		h_el_scaledown_fit=el_scaledown_fit.theHistogram
                		h_el_PUup_fit=el_PUup_fit.theHistogram
                		h_el_PUdown_fit=el_PUdown_fit.theHistogram
                		ratioGraphs_fit = ROOT.TGraphAsymmErrors(hhmu_fit.GetSize()-2)
                		for i in range(1, hhmu_fit.GetSize()-1):
                        		xval_fit = hhmu_fit.GetBinCenter(i)
                        		xerr_fit = hhmu_fit.GetBinWidth(i)/2
                        		if hhel_fit.GetBinContent(i) == 0: continue
                        		if hhmu_fit.GetBinContent(i) == 0: continue
                        #print(hhmu.GetBinContent(i))
                        #print(hhel.GetBinContent(i))
                        		yval_fit = hhmu_fit.GetBinContent(i)*1.0/hhel_fit.GetBinContent(i)
                        		hhmu_fit.SetBinError(i,math.sqrt(errmu_fit[i-1]+hhmu_fit.GetBinError(i)**2))
                        		hhel_fit.SetBinError(i,math.sqrt(errel_fit[i-1]+hhel_fit.GetBinError(i)**2))
                        		yerr_fit = yval_fit*math.sqrt((hhmu_fit.GetBinError(i)/hhmu_fit.GetBinContent(i))**2+(hhel_fit.GetBinError(i)/hhel_fit.GetBinContent(i))**2)
                        		ratioGraphs_fit.SetPoint(i, xval_fit, yval_fit)
                        		ratioGraphs_fit.SetPointError(i, xerr_fit, xerr_fit, yerr_fit, yerr_fit)
				#func.SetParLimits(2,0.8,1.2)
				#func.SetParLimits(1,1.5e-1,1.9e-1)
                        	ratioGraphs_fit.Fit(func,"","",200,3500)
				
				#ratioGraphs_fit.Fit(func,"","",200,3500)
				#ratioGraphs_fit.Fit(func,"","",200,3500)
				#ratioGraphs_fit.Fit(func,"","",200,3500)
				#ratioGraphs_fit.Fit(func,"","",200,3500)
				print(func.GetNDF())
				print(func.GetChisquare())
				print(func.GetChisquare()/func.GetNDF())	
				
	if args.ae: 
				
		if args.useall:
			i=0
			for year in range(2016,2019):
				print(year)
				for h in stackmu[i].theStack.GetHists(): 
					inverseAE(h, plot_mu["default"], year)
					h=Rebin(h)
				for h in stackel[i].theStack.GetHists(): 
					inverseAE(h, plot_el["default"], year)
					h=Rebin(h)
				inverseAE(stackmu[i].theHistogram, plot_mu["default"], year)
				stackmu[i].theHistogram=Rebin(stackmu[i].theHistogram)
				inverseAE(stackel[i].theHistogram, plot_el["default"], year)
				stackel[i].theHistogram=Rebin(stackel[i].theHistogram)
				inverseAE(mu_scaleup[i].theHistogram, plot_mu["scale_up"], year)
				mu_scaleup[i].theHistogram=Rebin(mu_scaleup[i].theHistogram)
				inverseAE(mu_scaledown[i].theHistogram, plot_mu["scale_down"], year)
				mu_scaledown[i].theHistogram=Rebin(mu_scaledown[i].theHistogram)
				inverseAE(mu_ID[i].theHistogram, plot_mu["ID"], year)
				mu_ID[i].theHistogram=Rebin(mu_ID[i].theHistogram)
				inverseAE(mu_reso[i].theHistogram, plot_mu["reso"], year)
				mu_reso[i].theHistogram=Rebin(mu_reso[i].theHistogram)
				inverseAE(el_scaleup[i].theHistogram, plot_el["scale_up"], year)
				el_scaleup[i].theHistogram=Rebin(el_scaleup[i].theHistogram)
				inverseAE(el_scaledown[i].theHistogram, plot_el["scale_down"], year)
				el_scaledown[i].theHistogram=Rebin(el_scaledown[i].theHistogram)
				inverseAE(el_PUup[i].theHistogram, plot_el["PU_up"], year)
				el_PUup[i].theHistogram=Rebin(el_PUup[i].theHistogram)
				inverseAE(el_PUdown[i].theHistogram, plot_el["PU_down"], year)
				el_PUdown[i].theHistogram=Rebin(el_PUdown[i].theHistogram)
				inverseAE(datamu[i], plot_mu["default"], year)
				datamu[i]=Rebin(datamu[i])
				inverseAE(datael[i], plot_el["default"], year)
				datael[i]=Rebin(datael[i])
				i+=1             
		else:
			if args.use2016: year = 2016
			elif args.use2018: year = 2018
			else: year =2017
			for h in stackmu.theStack.GetHists():
				inverseAE(h, plot_mu["default"], year)
				h=Rebin(h)
			for h in stackel.theStack.GetHists(): 
				inverseAE(h, plot_el["default"], year)
				h=Rebin(h)
			inverseAE(stackmu.theHistogram, plot_mu["default"], year)
			stackmu.theHistogram=Rebin(stackmu.theHistogram)
			inverseAE(stackel.theHistogram, plot_el["default"], year)
			stackel.theHistogram=Rebin(stackel.theHistogram)
			inverseAE(mu_scaleup.theHistogram, plot_mu["scale_up"], year)
			mu_scaleup.theHistogram=Rebin(mu_scaleup.theHistogram)
			inverseAE(mu_scaledown.theHistogram, plot_mu["scale_down"], year)
			mu_scaledown.theHistogram=Rebin(mu_scaledown.theHistogram)
			inverseAE(mu_ID.theHistogram, plot_mu["ID"], year)
			mu_ID.theHistogram=Rebin(mu_ID.theHistogram)
			inverseAE(mu_reso.theHistogram, plot_mu["reso"], year)
			mu_reso.theHistogram=Rebin(mu_reso.theHistogram)
			inverseAE(el_scaleup.theHistogram, plot_el["scale_up"], year)
			el_scaleup.theHistogram=Rebin(el_scaleup.theHistogram)
			inverseAE(el_scaledown.theHistogram, plot_el["scale_down"], year)
			el_scaledown.theHistogram=Rebin(el_scaledown.theHistogram)
			inverseAE(el_PUup.theHistogram, plot_el["PU_up"], year)
			el_PUup.theHistogram=Rebin(el_PUup.theHistogram)
			inverseAE(el_PUdown.theHistogram, plot_el["PU_down"], year)
			el_PUdown.theHistogram=Rebin(el_PUdown.theHistogram)
			inverseAE(datamu, plot_mu["default"], year)
			datamu=Rebin(datamu)
			inverseAE(datael, plot_el["default"], year)
			datael=Rebin(datael)
	#else:
		#if args.useall:
                        #i=0
                        #for year in range(2016,2019):
                                #for h in stackmu[i].theStack.GetHists():
                                #        h=Rebin(h)
                                #for h in stackel[i].theStack.GetHists():
                                #        h=Rebin(h)
                                #stackmu[i].theHistogram=Rebin(stackmu[i].theHistogram)
                                #stackel[i].theHistogram=Rebin(stackel[i].theHistogram)
                                #mu_scaleup[i].theHistogram=Rebin(mu_scaleup[i].theHistogram)
                                #mu_scaledown[i].theHistogram=Rebin(mu_scaledown[i].theHistogram)
                                #mu_ID[i].theHistogram=Rebin(mu_ID[i].theHistogram)
                                #mu_reso[i].theHistogram=Rebin(mu_reso[i].theHistogram)
                                #el_scaleup[i].theHistogram=Rebin(el_scaleup[i].theHistogram)
                                #el_scaledown[i].theHistogram=Rebin(el_scaledown[i].theHistogram)
                                #el_PUup[i].theHistogram=Rebin(el_PUup[i].theHistogram)
                                #el_PUdown[i].theHistogram=Rebin(el_PUdown[i].theHistogram)
                                #datamu[i]=Rebin(datamu[i])
                                #datael[i]=Rebin(datael[i])
                                #i+=1
                #else:
                        #if args.use2016: year = 2016
                        #elif args.use2018: year = 2018
                        #else: year =2017
                        #for h in stackmu.theStack.GetHists():
                        #        h=Rebin(h)
                        #for h in stackel.theStack.GetHists():
                        #        h=Rebin(h)
                        #stackmu.theHistogram=Rebin(stackmu.theHistogram)
                        #stackel.theHistogram=Rebin(stackel.theHistogram)
                        #mu_scaleup.theHistogram=Rebin(mu_scaleup.theHistogram)
                        #mu_scaledown.theHistogram=Rebin(mu_scaledown.theHistogram)
                        #mu_ID.theHistogram=Rebin(mu_ID.theHistogram)
                        #mu_reso.theHistogram=Rebin(mu_reso.theHistogram)
                        #el_scaleup.theHistogram=Rebin(el_scaleup.theHistogram)
                        #el_scaledown.theHistogram=Rebin(el_scaledown.theHistogram)
                        #el_PUup.theHistogram=Rebin(el_PUup.theHistogram)
                        #el_PUdown.theHistogram=Rebin(el_PUdown.theHistogram)
                        #datamu=Rebin(datamu)
                        #datael=Rebin(datael)
	#if args.znorm:
	#	if args.useall:
	#		for i in range(3):
	#			znum=znorm(stackmu[i].theHistogram)
         #               	stackmu[i].theHistogram.Scale(1./znum)
        #                	znum=znorm(stackel[i].theHistogram)
         #               	stackel[i].theHistogram.Scale(1./znum)
          #              	znum=znorm(mu_scaleup[i].theHistogram)
         #               	mu_scaleup[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(mu_scaledown[i].theHistogram)
         #               	mu_scaledown[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(mu_ID[i].theHistogram)
         #               	mu_ID[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(mu_reso[i].theHistogram)
         #               	mu_reso[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(el_scaleup[i].theHistogram)
         #               	el_scaleup[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(el_scaledown[i].theHistogram)
         #               	el_scaledown[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(el_PUup[i].theHistogram)
         #               	el_PUup[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(el_PUdown[i].theHistogram)
         #               	el_PUdown[i].theHistogram.Scale(1./znum)
         #               	znum=znorm(datamu[i])
         #               	datamu[i].Scale(1./znum)
         #               	znum=znorm(datael[i])
         #               	datael[i].Scale(1./znum)			
	#	else:
	#		znum=znorm(stackmu.theHistogram)
	#		stackmu.theHistogram.Scale(1./znum)
	#		znum=znorm(stackel.theHistogram)
         #               stackel.theHistogram.Scale(1./znum)	
         #       	znum=znorm(mu_scaleup.theHistogram)
          #              mu_scaleup.theHistogram.Scale(1./znum)
	#		znum=znorm(mu_scaledown.theHistogram)
         #               mu_scaledown.theHistogram.Scale(1./znum)
	#		znum=znorm(mu_ID.theHistogram)
        #                mu_ID.theHistogram.Scale(1./znum)
	#		znum=znorm(mu_reso.theHistogram)
        #                mu_reso.theHistogram.Scale(1./znum)
	#		znum=znorm(el_scaleup.theHistogram)
        #                el_scaleup.theHistogram.Scale(1./znum)
	#		znum=znorm(el_scaledown.theHistogram)
        #                el_scaledown.theHistogram.Scale(1./znum)
	#		znum=znorm(el_PUup.theHistogram)
        #                el_PUup.theHistogram.Scale(1./znum)
	#		znum=znorm(el_PUdown.theHistogram)
        #                el_PUdown.theHistogram.Scale(1./znum)
	#		znum=znorm(datamu)
	#		datamu.Scale(1./znum)
	#		znum=znorm(datael)
	#		datael.Scale(1./znum)
	if args.useall:
		i=0
		Errs_mu=[]
		Errs_el=[]
		for year in range(2016,2019):
			lis_mu=[[mu_scaleup[i].theHistogram,mu_scaledown[i].theHistogram],mu_ID[i].theHistogram,mu_reso[i].theHistogram]
			key_mu=["massScale","ID","resolution"]
			lis_el=[[el_scaleup[i].theHistogram,el_scaledown[i].theHistogram],[el_PUup[i].theHistogram,el_PUdown[i].theHistogram]]
			key_el=["massScale","pileUp"]
			Errs_mu.append(getErrors(stackmu[i].theHistogram,lis_mu,key_mu))
			Errs_el.append(getErrors(stackel[i].theHistogram,lis_el,key_el))
			i+=1
		if category=="bb":
			errmu=Errs_mu[0]["massScale"]**2+(Errs_mu[1]["massScale"]+Errs_mu[2]["massScale"])**2+Errs_mu[0]["resolution"]**2+Errs_mu[1]["resolution"]**2+Errs_mu[2]["resolution"]**2+Errs_mu[0]["ID"]**2+(Errs_mu[1]["ID"]+Errs_mu[2]["ID"])**2
			ErrMu={}
			ErrMu["massScale"]=Errs_mu[0]["massScale"]**2+(Errs_mu[1]["massScale"]+Errs_mu[2]["massScale"])**2
			ErrMu["resolution"]=Errs_mu[0]["resolution"]**2+Errs_mu[1]["resolution"]**2+Errs_mu[2]["resolution"]**2
			ErrMu["ID"]=Errs_mu[0]["ID"]**2+(Errs_mu[1]["ID"]+Errs_mu[2]["ID"])**2
			
		else:
			errmu=Errs_mu[0]["massScale"]**2+(Errs_mu[1]["massScale"]+Errs_mu[2]["massScale"])**2+Errs_mu[0]["resolution"]**2+Errs_mu[1]["resolution"]**2+Errs_mu[2]["resolution"]**2+Errs_mu[0]["ID"]**2+Errs_mu[1]["ID"]**2+Errs_mu[2]["ID"]**2
			ErrMu={}
			ErrMu["massScale"]=Errs_mu[0]["massScale"]**2+(Errs_mu[1]["massScale"]+Errs_mu[2]["massScale"])**2
			ErrMu["resolution"]=Errs_mu[0]["resolution"]**2+Errs_mu[1]["resolution"]**2+Errs_mu[2]["resolution"]**2
			ErrMu["ID"]=Errs_mu[0]["ID"]**2+Errs_mu[1]["ID"]**2+Errs_mu[2]["ID"]**2
		errel=Errs_el[0]["massScale"]**2+Errs_el[0]["pileUp"]**2+Errs_el[1]["massScale"]**2+Errs_el[1]["pileUp"]**2+Errs_el[2]["massScale"]**2+Errs_el[2]["pileUp"]**2
		ErrEl={}
		ErrEl["massScale"]=Errs_el[0]["massScale"]**2+Errs_el[1]["massScale"]**2+Errs_el[2]["massScale"]**2
		ErrEl["pileUp"]=Errs_el[0]["pileUp"]**2+Errs_el[1]["pileUp"]**2+Errs_el[2]["pileUp"]**2			
		scale_el=Errs_el[0]["massScale"]**2+Errs_el[1]["massScale"]**2+Errs_el[2]["massScale"]**2
		scale_mu=Errs_mu[0]["massScale"]**2+(Errs_mu[1]["massScale"]+Errs_mu[2]["massScale"])**2         
		stackmu=Addstack(stackmu)
		stackel=Addstack(stackel)
		mu_scaleup=Addstack(mu_scaleup)
		mu_scaledown=Addstack(mu_scaledown)
		mu_ID=Addstack(mu_ID)
		mu_reso=Addstack(mu_reso)
		el_scaleup=Addstack(el_scaleup)
		el_scaledown=Addstack(el_scaledown)
		el_PUup=Addstack(el_PUup)
		el_PUdown=Addstack(el_PUdown)
		datamu=Addhist(datamu)
		datael=Addhist(datael)               
	else:
		lis_mu=[[mu_scaleup.theHistogram,mu_scaledown.theHistogram],mu_ID.theHistogram,mu_reso.theHistogram]
		lis_el=[[el_scaleup.theHistogram,el_scaledown.theHistogram],[el_PUup.theHistogram,el_PUdown.theHistogram]]
		key_mu=["massScale","ID","resolution"]
		Errs_mu=getErrors(stackmu.theHistogram,lis_mu,key_mu)
		errmu=Errs_mu["massScale"]**2+Errs_mu["resolution"]**2+Errs_mu["ID"]**2
		key_el=["massScale","pileUp"]
		Errs_el=getErrors(stackel.theHistogram,lis_el,key_el)
		errel=Errs_el["massScale"]**2+Errs_el["pileUp"]**2
		scale_el=Errs_el["massScale"]**2
                scale_mu=Errs_mu["massScale"]**2
	if args.data:
		yMax = datamu.GetBinContent(datamu.GetMaximumBin())
		if "Mass" in plot_mu["default"].fileName:
			yMin = 0.00001
		else:
			yMin = 0.01
		xMax = datamu.GetXaxis().GetXmax()
		xMin = datael.GetXaxis().GetXmin()
	else:	
		yMax = stackmu.theHistogram.GetBinContent(datamu.GetMaximumBin())
		yMin = 0.01
		xMax = stackmu.theHistogram.GetXaxis().GetXmax()
		xMin = stackmu.theHistogram.GetXaxis().GetXmin()	
	if plot_mu["default"].yMax == None:
		if logScale:
			yMax = yMax*10000
		else:
			yMax = yMax*1.5
	else: yMax = plot_mu["default"].yMax
	
	if "Mass" in plot_mu["default"].fileName:
		yMax = 20000000	
	
	if not plot_mu["default"].yMin == None:
		yMin = plot_mu.yMin
	if not plot_mu["default"].xMin == None:
		xMin = plot_mu["default"].xMin
	if not plot_mu["default"].xMax == None:
		xMax = plot_mu["default"].xMax


	xMin = 200
	xMax = 3490
	yMin = 1e-3
	yMax = 1e4

	if args.ae: 
		yMin = 0.00001 / 40
		yMax = 200000000.0 / 40
		xMin = 200
		xMax = 3490
		yMin *= 10000
		yMax /= 10
	if args.ae:
		# ~ vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","3 years data"))
		vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","Events / GeV * 1/(acc. x eff)"))
	else:
		# ~ vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","Lumi #times d#sigma(pp#rightarrow ll)"))
		vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","Events / GeV"))
	vh.GetXaxis().SetMoreLogLabels()
	
	drawStack_mu = stackmu
	drawStack_el = stackel

	
	# Draw background from stack
	drawStack_mu.theHistogram.SetLineColor(ROOT.kBlue-3)
	drawStack_el.theHistogram.SetLineColor(ROOT.kRed-3)
	drawStack_mu.theHistogram.SetMarkerStyle(6)
        drawStack_el.theHistogram.SetMarkerStyle(6)
	drawStack_mu.theHistogram.SetFillStyle(0)
        drawStack_el.theHistogram.SetFillStyle(0)
	#for i in range(1, drawStack_mu.theHistogram.GetSize()-1):
	#	val=drawStack_mu.theHistogram.GetBinContent(i)
	#	drawStack_mu.theHistogram.SetBinError(i,math.sqrt(errmu[i-1]+drawStack_mu.theHistogram.GetBinError(i)**2))
	#	drawStack_mu.theHistogram.SetBinContent(i,val)
	#	val=drawStack_el.theHistogram.GetBinContent(i)
	#	drawStack_el.theHistogram.SetBinError(i,math.sqrt(errel[i-1]+drawStack_el.theHistogram.GetBinError(i)**2))
	#	drawStack_el.theHistogram.SetBinContent(i,val)
	#	val=datamu.GetBinContent(i)
	#	datamu.SetBinError(i,math.sqrt(scale_mu[i-1]+datamu.GetBinError(i)**2))
	#	datamu.SetBinContent(i,val)
	#	val=datael.GetBinContent(i)
	#	datael.SetBinError(i,math.sqrt(scale_el[i-1]+datael.GetBinError(i)**2))
	#	datael.SetBinContent(i,val)
	#drawStack_mu.theHistogram.SetLineWidth(2)
	#drawStack_el.theHistogram.SetLineWidth(2)
	drawStack_mu.theHistogram.Draw("same histe ")
	#for i in range(1, drawStack_mu.theHistogram.GetSize()-1):
		#print(drawStack_mu.theHistogram.GetBinContent(i))
		#print(drawStack_el.theHistogram.GetBinContent(i))
	drawStack_el.theHistogram.Draw("same histe ")
	drawStack_mu.theHistogram.SetLineColor(ROOT.kBlue)
        drawStack_el.theHistogram.SetLineColor(ROOT.kRed)
        drawStack_mu.theHistogram.SetFillStyle(0)
        drawStack_el.theHistogram.SetFillStyle(0)


	# Draw data
	datamu.SetMinimum(0.0001)
	Box = ROOT.TBox(0.39,0.79,0.5,0.9)
        Box.SetLineColor(ROOT.kRed)
        Box.SetFillStyle(0)
        #Box.SetFillAlphaColor(0,0)
        Box.SetLineWidth(2)
	if args.data:
		datamu.SetMarkerColor(ROOT.kViolet)
		datael.SetMarkerColor(ROOT.kOrange)
		#for i in range(1, datamu.GetSize()-1):
                	#print(datamu.GetBinContent(i))
                	#print(datael.GetBinContent(i))
		datamu.Draw("ESAMEp")	
		datael.Draw("ESAMEp")
		Box.Draw("same")
		#a1=datael.Integral(60,120)
		#a2=drawStack_el.theHistogram.Integral(60,120)
		#print(a1/a2)
	# Draw legend
	if "Eta" in plot_mu["default"].fileName or "CosTheta" in plot_mu["default"].fileName:
		legendEta.Draw()
	else:
		legend.Draw()

	plotPad.SetLogx(plot_mu["default"].logX)
	if args.useall:
		latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV, #mu#mu), %.1f fb^{-1} (13 TeV, ee)"%((lumi_mu[0]+lumi_mu[1]+lumi_mu[2])*0.001, (lumi_el[0]+lumi_el[1]+lumi_el[2])*0.001))
	else:	
		latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV, #mu#mu), %.1f fb^{-1} (13 TeV, ee)"%(lumi_mu*0.001, lumi_el*0.001))
	yLabelPos = 0.85
	cmsExtra = "Preliminary"
	#if not args.data:
		#cmsExtra = "#splitline{Private Work}{Simulation}"
		#yLabelPos = 0.82
	if category=="bb":
		if args.useall:
			dataLabel = "Full Run 2   BB"
		elif args.use2016:
			dataLabel = "2016   BB"	
		elif args.use2018:
			dataLabel = "2018   BB"
		else:
			dataLabel = "2017   BB"
	else:
                if args.useall:
                        dataLabel = "Full Run 2   BE"
                elif args.use2016:
                        dataLabel = "2016   BE"      
                elif args.use2018:
                        dataLabel = "2018   BE"
                else:
                        dataLabel = "2017   BE"	
	pave=ROOT.TPaveText(0.39,0.79,0.5,0.9)
	pave.AddBox(0.39,0.79,0.5,0.9)
	pave.Draw()
	latexCMS.DrawLatex(0.19,0.89,"CMS")
	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
	latexCMSExtra.DrawLatex(0.39,0.79,"%s"%(dataLabel))	
	if args.ratio:
		ratioPad.cd()
		ratioPad.SetLogx(plot_mu["default"].logX)

		hhmu = drawStack_mu.theHistogram
		hhel = drawStack_el.theHistogram
		h_mu_scaleup=mu_scaleup.theHistogram
		h_mu_scaledown=mu_scaledown.theHistogram
		h_mu_ID=mu_ID.theHistogram
		h_mu_reso=mu_reso.theHistogram
		h_el_scaleup=el_scaleup.theHistogram
		h_el_scaledown=el_scaledown.theHistogram
		h_el_PUup=el_PUup.theHistogram
		h_el_PUdown=el_PUdown.theHistogram
		ratioGraphs = ROOT.TGraphAsymmErrors(hhmu.GetSize()-2)
		chann = True if "BB" in plot_mu["default"].fileName else False
		#print(errel)
		#print(errmu)
		#print(hhmu.GetSize())
		for i in range(1, hhmu.GetSize()-1):
			xval = hhmu.GetBinCenter(i)
			xerr = hhmu.GetBinWidth(i)/2
			if hhel.GetBinContent(i) == 0: continue
			if hhmu.GetBinContent(i) == 0: continue
			muScale=math.sqrt(ErrMu["massScale"][i-1])/hhmu.GetBinContent(i)
			muReso=math.sqrt(ErrMu["resolution"][i-1])/hhmu.GetBinContent(i)
			muID=math.sqrt(ErrMu["ID"][i-1])/hhmu.GetBinContent(i)
			muStat=hhmu.GetBinError(i)/hhmu.GetBinContent(i)
			elScale=math.sqrt(ErrEl["massScale"][i-1])/hhel.GetBinContent(i)
			elPile=math.sqrt(ErrEl["pileUp"][i-1])/hhel.GetBinContent(i)
			elStat=hhel.GetBinError(i)/hhel.GetBinContent(i)
			print "mass scale %f, resolution %f, ID %f, stat %f, mass %f, dimuon channel" %(muScale, muReso, muID, muStat, xval)
			print "mass scale %f, pile up %f, stat %f, mass %f, dielectron channel" %(elScale, elPile, elStat, xval)
			#print(hhmu.GetBinContent(i))
			#print(hhel.GetBinContent(i))
			yval = hhmu.GetBinContent(i)*1.0/hhel.GetBinContent(i)
			hhmu.SetBinError(i,math.sqrt(errmu[i-1]+hhmu.GetBinError(i)**2))
			hhel.SetBinError(i,math.sqrt(errel[i-1]+hhel.GetBinError(i)**2))
			
			yerr = yval*math.sqrt((hhmu.GetBinError(i)/hhmu.GetBinContent(i))**2+(hhel.GetBinError(i)/hhel.GetBinContent(i))**2)
			ratioGraphs.SetPoint(i, xval, yval)
                        ratioGraphs.SetPointError(i, xerr, xerr, yerr, yerr)
		if args.alter:
			#h_fit=deepcopy(hhmu)
			#hhmu.Sumw2()
			#hhel.Sumw2()
			#h_fit.Divide(hhmu,hhel)	
			#func1=ROOT.TF1("func1","[0]*log(x)+[1]",200,3500)
			#ratioGraphs.Fit(func1,"","",200,3500)
			#print(func1.GetChisquare()/func1.GetNDF())
			print("M = %f, r+-e = %f +- %f"%(xval, yval, yerr/yval))
		ratioData = ROOT.TGraphAsymmErrors(datamu.GetSize()-2)
		for i in range(1, datamu.GetSize()-1):
			xval = datamu.GetBinCenter(i)
			xerr = datamu.GetBinWidth(i)/2
			if datael.GetBinContent(i) == 0: continue
			if datamu.GetBinContent(i) == 0: continue
			#print(datael.GetBinContent(i))
			#print(datael.GetBinError(i))
			yval = datamu.GetBinContent(i)*1.0/datael.GetBinContent(i)
			yerr = yval*math.sqrt((datamu.GetBinError(i)/datamu.GetBinContent(i))**2+(datael.GetBinError(i)/datael.GetBinContent(i))**2+scale_mu[i-1]/(datamu.GetBinContent(i)**2)+scale_el[i-1]/(datael.GetBinContent(i)**2))
			print(yval)
			print(yerr)
			ratioData.SetPoint(i, xval, yval)
			ratioData.SetPointError(i, xerr, xerr, yerr, yerr)
			print("M = %f, r+-e = %f +- %f"%(xval, yval, yerr/yval))
		nBinsX = 20
		nBinsY = 10
		if args.ae: 
			hAxis = ROOT.TH2F("hAxis", "", nBinsX, xMin, xMax, nBinsY, 0.5, 2.5)
		else:	
			hAxis = ROOT.TH2F("hAxis", "", nBinsX, xMin, xMax, nBinsY, 0.5, 2.5)
		hAxis.Draw("AXIS")

		hAxis.GetYaxis().SetNdivisions(408)
		hAxis.SetTitleOffset(0.4, "Y")
		hAxis.SetTitleSize(0.09, "Y")
		hAxis.SetTitleSize(0.06, "X")
		hAxis.SetYTitle("R_{#mu#mu/ee}")
		hAxis.SetXTitle("m(l^{+}l^{-}) [GeV]")
		hAxis.GetXaxis().SetLabelSize(0.048)
		hAxis.GetYaxis().SetLabelSize(0.048)
		#hAxis.GetXaxis().SetTicks("+")
		#hAxis.SetTitleSize(0.15, "Y")
		hAxis.GetXaxis().SetMoreLogLabels()	
		oneLine = ROOT.TLine(xMin, 1.0, xMax, 1.0)
		oneLine.SetLineStyle(2)
		oneLine.Draw()
	
		ratioGraphs.SetFillColor(ROOT.kBlue-3)
		ratioGraphs.SetMarkerColor(ROOT.kBlue-3)
		ratioGraphs.GetXaxis().SetLabelSize(0.0)
		ratioGraphs.SetFillStyle(3002)	
		ratioGraphs.Draw("same p")
		#func.Draw("same")
		ratioData.SetMarkerColor(ROOT.kViolet)
		ratioData.Draw("same p")
		if args.alter:
			Expression="f(x)=p_{0}+p_{1}lnx+p_{2}x^{2}"
			latexCMSExtra.DrawLatex(0.59,yLabelPos,"%s"%(Expression))
			func.Draw("same")
		rlegend = TLegend(0.2, 0.65, 0.5, 0.925)
		rlegend.SetFillStyle(0)
		rlegend.SetBorderSize(1)
		rlegend.SetTextFont(42)
		rlegend.AddEntry(ratioGraphs, "mumu/ee in MC inclusive", "pe")
		rlegend.AddEntry(ratioData, "mumu/ee in data", "pe")
		rlegend.Draw("same")
		
		ratioPad.Update()
					
	
	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	if args.ratio:
		ratioPad.RedrawAxis()
	if not os.path.exists("lepFlavor"):
		os.makedirs("lepFlavor")	

	if args.use2016: year = "2016"
	elif args.useall: year = "2016_to_2018"
	elif args.use2018: year = "2018"
	else: year = "2017"

	if args.ae: year += "_inverseAE"
	#if args.znorm: year += "_znorm"
	
	hCanvas.Print("lepFlavor/%s_%s_datamc.pdf"%(plot_mu["default"].fileName, year))
	if args.datadiviMC:
		ers=0
		dmCanvas = TCanvas("hCanvas", "Distribution", 800,800)
		ratioDataMC = ROOT.TGraphAsymmErrors(datamu.GetSize()-2)
		for i in range(1, datamu.GetSize()-1):
			xval = datamu.GetBinCenter(i)
			xerr = datamu.GetBinWidth(i)/2
			if datael.GetBinContent(i) == 0: continue
			if datamu.GetBinContent(i) == 0: continue
			dataval = datamu.GetBinContent(i)*1.0/datael.GetBinContent(i)
			dataerr = math.sqrt((datamu.GetBinError(i)/datamu.GetBinContent(i))**2+(datael.GetBinError(i)/datael.GetBinContent(i))**2+scale_mu[i-1]/(datamu.GetBinContent(i)**2)+scale_el[i-1]/(datael.GetBinContent(i)**2))
			if hhel.GetBinContent(i) == 0: continue
			if hhmu.GetBinContent(i) == 0: continue
			MCval = hhmu.GetBinContent(i)*1.0/hhel.GetBinContent(i)
			MCerr = math.sqrt((errel[i-1]**0.5/hhel.GetBinContent(i))**2+(errmu[i-1]**0.5/hhmu.GetBinContent(i))**2+(hhmu.GetBinError(i)/hhmu.GetBinContent(i))**2+(hhel.GetBinError(i)/hhel.GetBinContent(i))**2)
			yval = dataval/MCval
			yerr=yval*math.sqrt((MCerr)**2+(dataerr)**2)
			ratioDataMC.SetPointError(i, xerr, xerr, yerr, yerr)
                        print ("M = %f, r+-e = %f +- %f"%(xval, yval, yerr))
	if args.alter:
		ers=0
                dmCanvas = TCanvas("hCanvas", "Distribution", 800,800)
                ratioDataMC = ROOT.TGraphAsymmErrors(datamu.GetSize()-2)
		for i in range(1, datamu.GetSize()-1):
			xval=datamu.GetBinCenter(i)
			fval=func.Eval(xval)
			print "At the bin center with the mass %f, scaling function value is %f, dimuon event yield is %f, dielectron event yield is %f"%(datamu.GetBinCenter(i), fval, datamu.GetBinContent(i)*datamu.GetBinWidth(i), datael.GetBinContent(i)*datael.GetBinWidth(i))
			muval=datamu.GetBinContent(i)/func.Eval(xval)
			ParErr0=func.GetParError(0)
			ParErr1=func.GetParError(1)
			ParErr2=func.GetParError(2)
			fErr=math.sqrt(ParErr0**2+(ParErr1*math.log(xval))**2+(ParErr2*xval**2)**2)
			#fErr=0
			#err=math.sqrt(datamu.GetBinError(i)**2+scale_mu[i-1]+muval**2*fErr**2)/fval
			err=muval*fErr/fval
			datamu.SetBinContent(i,muval)
			datamu.SetBinError(i,err)
			#elerr=math.sqrt(datael.GetBinError(i)**2+scale_el[i-1])
			#datael.SetBinError(i,elerr)
		datamu.Sumw2()
		datamu=Rebin(datamu)
		datael.Sumw2()
		datael=Rebin(datael)
		for i in range(1, datamu.GetSize()-1):
			xval=datamu.GetBinCenter(i)
			xerr=datamu.GetBinWidth(i)
			if datael.GetBinContent(i) == 0: continue
                        if datamu.GetBinContent(i) == 0: continue
			yval=datamu.GetBinContent(i)/datael.GetBinContent(i)
			yerr=yval*math.sqrt((datamu.GetBinError(i)/datamu.GetBinContent(i))**2+(datael.GetBinError(i)/datael.GetBinContent(i))**2)
			ratioDataMC.SetPoint(i, xval, yval)
			ratioDataMC.SetPointError(i, xerr/2, xerr/2, yerr, yerr)
			print ("M = %f, r+-e = %f +- %f"%(xval, yval, yerr))
			

		#print("max err")
		#print(ers)
		dmPad = ROOT.TPad("dmPad","dmPad",0,0,1,1)
		dmPad.SetLogx(True)
                #setTDRStyle()
                #dmPad.UseCurrentStyle()
		Box.Draw("same")
                dmPad.Draw()
                dmPad.cd()
		dmAxis = ROOT.TH2F("dmAxis", "", nBinsX, xMin, xMax, nBinsY, 0.5, 2.5)
		dmAxis.Draw("AXIS")
                dmAxis.GetYaxis().SetNdivisions(408)
                dmAxis.SetTitleOffset(1.1, "Y")
		dmAxis.SetTitleOffset(1.2, "X")
                dmAxis.SetTitleSize(0.06, "Y")
                dmAxis.SetTitleSize(0.04, "X")
		if args.datadiviMC:
                	dmAxis.SetYTitle("R^{data}_{#mu#mu/ee}/R^{mc}_{#mu#mu/ee}")
		else:
                        dmAxis.SetYTitle("R^{data}_{#mu#mu/ee}")
                dmAxis.SetXTitle("m(l^{+}l^{-}) [GeV]")
                dmAxis.GetXaxis().SetLabelSize(0.048)
                dmAxis.GetYaxis().SetLabelSize(0.048)
                #hAxis.GetXaxis().SetTicks("+")
                #hAxis.SetTitleSize(0.15, "Y")
                dmAxis.GetXaxis().SetMoreLogLabels()
                dmLine = ROOT.TLine(xMin, 1.0, xMax, 1.0)
                dmLine.SetLineStyle(2)
                dmLine.Draw()

                ratioDataMC.SetFillColor(ROOT.kBlue-3)
                ratioDataMC.SetMarkerColor(ROOT.kBlue-3)
                ratioDataMC.GetXaxis().SetLabelSize(0.0)
                ratioDataMC.SetFillStyle(3002)
                ratioDataMC.Draw("SAME p")
		dmPad.Update()
		dmPad.RedrawAxis()
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
		latexCMSExtra.DrawLatex(0.59,yLabelPos,"%s"%(dataLabel))
		
		if args.useall:
                	latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV, #mu#mu), %.1f fb^{-1} (13 TeV, ee)"%((lumi_mu[0]+lumi_mu[1]+lumi_mu[2])*0.001, (lumi_el[0]+lumi_el[1]+lumi_el[2])*0.001))
        	else:
                	latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV, #mu#mu), %.1f fb^{-1} (13 TeV, ee)"%(lumi_mu*0.001, lumi_el*0.001))
		dmCanvas.Print("lepFlavor/%s_%s_datadivimc.pdf"%(plot_mu["default"].fileName, year))
					
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
						  help="use 2018 data with 2018 MC.")
	parser.add_argument("-r", "--ratio", action="store_true", dest="ratio", default=False,
						  help="plot ratio plot")
	parser.add_argument("-l", "--log", action="store_true", dest="log", default=False,
						  help="plot with log scale for y axis")
	parser.add_argument("-b", "--backgrounds", dest="backgrounds", action="append", default=[],
						  help="backgrounds to plot.")
	parser.add_argument("--ae", action="store_true", dest="ae", default=False,help="times inverse Acceptance x Efficiency")
	parser.add_argument("--alter", action="store_true", dest="alter", default=False, help="Another approach")
	parser.add_argument("--all", action="store_true", dest="useall", default=False, help="add the data from 2016 to 2018")
	parser.add_argument("-c", action="store_true", dest="datadiviMC", default=False, help="get data/MC")
	args = parser.parse_args()
	if len(args.backgrounds) == 0:
		#args.backgrounds = ["Wjets","Other","DrellYan"]
		args.backgrounds = ["DrellYan","Other"]
		
	muplots_bb={"default":"massPlotBB","scale_up":"massPlotBBScaleUpNoLog","scale_down":"massPlotBBScaleDownNoLog","reso":"massPlotBBSmearNoLog","ID":"massPlotBBMuonIDNoLog"}
	muplots_be={"default":"massPlotBE","scale_up":"massPlotBEScaleUpNoLog","scale_down":"massPlotBEScaleDownNoLog","reso":"massPlotBESmearNoLog","ID":"massPlotBEMuonIDNoLog"}
	eleplots_bb={"default":"massPlotEleBB","scale_up":"massPlotEleBBScaleUpNoLog","scale_down":"massPlotEleBBScaleDownNoLog","PU_up":"massPlotEleBBPUScaleUpNoLog","PU_down":"massPlotEleBBPUScaleDownNoLog"}
	eleplots_be={"default":"massPlotEleBE","scale_up":"massPlotEleBEScaleUpNoLog","scale_down":"massPlotEleBEScaleDownNoLog","PU_up":"massPlotEleBEPUScaleUpNoLog","PU_down":"massPlotEleBEPUScaleDownNoLog"}
	plot_mu_bb={}
	plot_el_bb={}
	plot_mu_be={}
	plot_el_be={}
	for key in muplots_bb.keys():
		plot_mu_bb[key] = getPlot(muplots_bb[key])
		plot_mu_bb[key].logX=True
	for key in eleplots_bb.keys():
		plot_el_bb[key] = getPlot(eleplots_bb[key])
		plot_el_bb[key].logX=True
	for key in muplots_be.keys():
		plot_mu_be[key] = getPlot(muplots_be[key])
		plot_mu_be[key].logX=True
	for key in eleplots_be.keys():
		plot_el_be[key] = getPlot(eleplots_be[key])
		plot_el_be[key].logX=True

	plotDataMC(args,plot_mu_bb,plot_el_bb,"bb")
	plotDataMC(args,plot_mu_be,plot_el_be,"be")
