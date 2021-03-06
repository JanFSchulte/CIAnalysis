#!/bin/env python

import sys, os
import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("-inFile", help="Input file", type=str)
parser.add_argument("-flav", help="Lepton flavor", type=str)
parser.add_argument("-unc",  help="Uncertainty: 'nominal'*, 'scaleup', 'scaledown', 'muonid', 'smeared', 'pdfWeightsUp', 'pdfWeightsDown', 'prefireUp', 'prefireDown'", type=str, default="nominal")
parser.add_argument("-cs",   help="CS bin: 'inc', 'cspos', 'csneg'", type=str, default="inc")
parser.add_argument("-d",    help="debug", action='store_true')
parser.add_argument("-do2016", help="do 2016", action='store_true')
parser.add_argument("-do2018", help="do 2018", action='store_true')
parser.add_argument("-add", help="add", action="store_true")
parser.add_argument("-truncation", help="truncation", action="store_true")

args = parser.parse_args()


outDir = "fitPlots/"
if args.truncation: 
	outDir = "fitPlotsTruncation/"


import ROOT as r
import numpy as np
#~ from nesteddict import nesteddict as ndict
import json
from copy import copy

from defs import getPlot, Backgrounds, Signals, Data, path, Signals2016, Signals2016ADD, SignalsADD, Signals2018, Signals2018ADD, zScale2016, zScale2018
from helpers import *

from setTDRStyle import setTDRStyle
setTDRStyle()


def applyPDFCorrection(hist):
	
	for i in range(0,hist.GetNbinsX()+1):
		binCenter = hist.GetBinCenter(i)
		scaleFac = 0.86 - 3.72e-05 * binCenter + 2.72e-08 * binCenter **2
		hist.SetBinContent(i,hist.GetBinContent(i)*scaleFac)
		
	return copy(hist)
	
	
def truncateADD(signalhist, signalhist100,l):
	
	returnHist = copy(signalhist)
	for i in range(0,returnHist.GetNbinsX()+1):
		binCenter = returnHist.GetBinCenter(i)
		if binCenter > l:
			returnHist.SetBinContent(i,signalhist100.GetBinContent(i))
		
	return copy(returnHist)
antypes=[
	# ["E","e","Ele","/store/user/sturdy/ZprimeAnalysis/histosHLTWeighted"],
	# ["Mu","mu","Mu","/store/user/sturdy/ZprimeAnalysis/histosCutHLT"]
	["E","e","Ele"],
	["Mu","mu","Mu"]
	]

r.gROOT.SetBatch(True)

lvals=["1","10","16", "24", "32", "40", "100k"]
if args.do2016:
	lvals=["1", "10", "16", "22", "28", "34", "100k"]
helis=["LL","LR","RL","RR"]
intfs=["Con","Des"]
# ~ helis=["LL"]
# ~ intfs=["Con"]
supers = [400,500,700,1100,1900,3500,10000]
extrabins = [1000+i for i in range(0, 2500, 200)]

if not os.path.exists(outDir):
	os.mkdir(outDir)

if args.add:
	if args.do2016:
		lvals = ["%.1f"%(4.0+i*0.5) for i in range(11)]
		lvals.remove("6.5")
		lvals.append("10")
		lvals.append("15")
		lvals.append("30")
		lvals.append("50")
		lvals.append("100")
		lvals.append("1000")
		lvals.append("100k")
		helis = [""]
		intfs = [""]
		supers = [1900, 2200, 2600, 3000, 3400, 10000]
		extrabins = [1900+i for i in range(0, 1500, 200)]
	else:
		lvals = ["%.1f"%(4+i*1) for i in range(9)]
		lvals.append("100")
		helis = [""]
		intfs = [""]
		supers = [1800, 2200, 2600, 3000, 3400, 10000]
		extrabins = [1900+i for i in range(0, 1500, 200)]
#	else:
#		lvals = ["%.1f"%(4+i*1) for i in range(9)]
#		lvals.append("100")
#		helis = [""]
#		intfs = [""]
#		supers = [400, 700, 1500, 2500, 3500, 10000]
#		extrabins = [1900+i for i in range(0, 1500, 200)]

uncertainties = [
	"nominal",
	"scaleup",
	"scaledown",
	"pdfWeightsUp",
	"pdfWeightsDown",
	## ele only
	"pileup",
	"piledown",
	"prefireup",
	"prefiredown",
	## muon only
	"smeared",
	"muonid",
	]
	
plots = {

	"Mubbnominal": "massPlotBBNoLog",
	"Mubenominal": "massPlotBENoLog",
	"Elebbnominal": "massPlotEleBBNoLog",
	"Elebenominal": "massPlotEleBENoLog",
	"MubbpdfWeightsUp": "massPlotBBNoLog",
	"MubepdfWeightsUp": "massPlotBENoLog",
	"ElebbpdfWeightsUp": "massPlotEleBBNoLog",
	"ElebepdfWeightsUp": "massPlotEleBENoLog",
	"MubbpdfWeightsDown": "massPlotBBNoLog",
	"MubepdfWeightsDown": "massPlotBENoLog",
	"ElebbpdfWeightsDown": "massPlotEleBBNoLog",
	"ElebepdfWeightsDown": "massPlotEleBENoLog",
	"Mubbscaleup": "massPlotBBScaleUpNoLog",
	"Mubescaleup": "massPlotBEScaleUpNoLog",
	"Elebbscaleup": "massPlotEleBBScaleUpNoLog",
	"Elebescaleup": "massPlotEleBEScaleUpNoLog",
	"Mubbscaledown": "massPlotBBScaleDownNoLog",
	"Mubescaledown": "massPlotBEScaleDownNoLog",
	"Elebbscaledown": "massPlotEleBBScaleDownNoLog",
	"Elebescaledown": "massPlotEleBEScaleDownNoLog",

	"Elebbpileup": "massPlotEleBBPUScaleUpNoLog",
	"Elebepileup": "massPlotEleBEPUScaleUpNoLog",
	"Elebbpiledown": "massPlotEleBBPUScaleDownNoLog",
	"Elebepiledown": "massPlotEleBEPUScaleDownNoLog",
	"Elebbprefireup": "massPlotEleBBPrefireUpNoLog",
	"Elebeprefireup": "massPlotEleBEPrefireUpNoLog",
	"Elebbprefiredown": "massPlotEleBBPrefireDownNoLog",
	"Elebeprefiredown": "massPlotEleBEPrefireDownNoLog",
	
	"Mubbsmeared": "massPlotBBSmearNoLog",
	"Mubesmeared": "massPlotBESmearNoLog",	
	"Mubbmuonid": "massPlotBBMuonIDNoLog",
	"Mubemuonid": "massPlotBEMuonIDNoLog",	

	"Mubbnominalcspos": "massPlotBBCSPosNoLog",
	"Mubenominalcspos": "massPlotBECSPosNoLog",
	"Elebbnominalcspos": "massPlotEleBBCSPosNoLog",
	"Elebenominalcspos": "massPlotEleBECSPosNoLog",

	"MubbpdfWeightsUpcspos": "massPlotBBCSPosNoLog",
	"MubepdfWeightsUpcspos": "massPlotBECSPosNoLog",
	"ElebbpdfWeightsUpcspos": "massPlotEleBBCSPosNoLog",
	"ElebepdfWeightsUpcspos": "massPlotEleBECSPosNoLog",

	"MubbpdfWeightsDowncspos": "massPlotBBCSPosNoLog",
	"MubepdfWeightsDowncspos": "massPlotBECSPosNoLog",
	"ElebbpdfWeightsDowncspos": "massPlotEleBBCSPosNoLog",
	"ElebepdfWeightsDowncspos": "massPlotEleBECSPosNoLog",
	
	"Mubbscaleupcspos": "massPlotBBScaleUpCSPosNoLog",
	"Mubescaleupcspos": "massPlotBEScaleUpCSPosNoLog",
	"Elebbscaleupcspos": "massPlotEleBBScaleUpCSPosNoLog",
	"Elebescaleupcspos": "massPlotEleBEScaleUpCSPosNoLog",
	"Mubbscaledowncspos": "massPlotBBScaleDownCSPosNoLog",
	"Mubescaledowncspos": "massPlotBEScaleDownCSPosNoLog",
	"Elebbscaledowncspos": "massPlotEleBBScaleDownCSPosNoLog",
	"Elebescaledowncspos": "massPlotEleBEScaleDownCSPosNoLog",

	"Elebbpileupcspos": "massPlotEleBBPUScaleUpCSPosNoLog",
	"Elebepileupcspos": "massPlotEleBEPUScaleUpCSPosNoLog",
	"Elebbpiledowncspos": "massPlotEleBBPUScaleDownCSPosNoLog",
	"Elebepiledowncspos": "massPlotEleBEPUScaleDownCSPosNoLog",

	"Elebbprefireupcspos": "massPlotEleBBPrefireUpCSPosNoLog",
	"Elebeprefireupcspos": "massPlotEleBEPrefireUpCSPosNoLog",
	"Elebbprefiredowncspos": "massPlotEleBBPrefireDownCSPosNoLog",
	"Elebeprefiredowncspos": "massPlotEleBEPrefireDownCSPosNoLog",
	
	"Mubbsmearedcspos": "massPlotBBSmearCSPosNoLog",
	"Mubesmearedcspos": "massPlotBESmearCSPosNoLog",	
	"Mubbmuonidcspos": "massPlotBBMuonIDCSPosNoLog",
	"Mubemuonidcspos": "massPlotBEMuonIDCSPosNoLog",	
	
	"Mubbnominalcsneg": "massPlotBBCSNegNoLog",
	"Mubenominalcsneg": "massPlotBECSNegNoLog",
	"Elebbnominalcsneg": "massPlotEleBBCSNegNoLog",
	"Elebenominalcsneg": "massPlotEleBECSNegNoLog",
	"MubbpdfWeightsUpcsneg": "massPlotBBCSNegNoLog",
	"MubepdfWeightsUpcsneg": "massPlotBECSNegNoLog",
	"ElebbpdfWeightsUpcsneg": "massPlotEleBBCSNegNoLog",
	"ElebepdfWeightsUpcsneg": "massPlotEleBECSNegNoLog",
	"MubbpdfWeightsDowncsneg": "massPlotBBCSNegNoLog",
	"MubepdfWeightsDowncsneg": "massPlotBECSNegNoLog",
	"ElebbpdfWeightsDowncsneg": "massPlotEleBBCSNegNoLog",
	"ElebepdfWeightsDowncsneg": "massPlotEleBECSNegNoLog",
	"Mubbscaleupcsneg": "massPlotBBScaleUpCSNegNoLog",
	"Mubescaleupcsneg": "massPlotBEScaleUpCSNegNoLog",
	"Elebbscaleupcsneg": "massPlotEleBBScaleUpCSNegNoLog",
	"Elebescaleupcsneg": "massPlotEleBEScaleUpCSNegNoLog",
	"Mubbscaledowncsneg": "massPlotBBScaleDownCSNegNoLog",
	"Mubescaledowncsneg": "massPlotBEScaleDownCSNegNoLog",
	"Elebbscaledowncsneg": "massPlotEleBBScaleDownCSNegNoLog",
	"Elebescaledowncsneg": "massPlotEleBEScaleDownCSNegNoLog",

	"Elebbpileupcsneg": "massPlotEleBBPUScaleUpCSNegNoLog",
	"Elebepileupcsneg": "massPlotEleBEPUScaleUpCSNegNoLog",
	"Elebbpiledowncsneg": "massPlotEleBBPUScaleDownCSNegNoLog",
	"Elebepiledowncsneg": "massPlotEleBEPUScaleDownCSNegNoLog",
	"Elebbprefireupcsneg": "massPlotEleBBPrefireUpCSNegNoLog",
	"Elebeprefireupcsneg": "massPlotEleBEPrefireUpCSNegNoLog",
	"Elebbprefiredowncsneg": "massPlotEleBBPrefireDownCSNegNoLog",
	"Elebeprefiredowncsneg": "massPlotEleBEPrefireDownCSNegNoLog",
	
	"Mubbsmearedcsneg": "massPlotBBSmearCSNegNoLog",
	"Mubesmearedcsneg": "massPlotBESmearCSNegNoLog",	
	"Mubbmuonidcsneg": "massPlotBBMuonIDCSNegNoLog",
	"Mubemuonidcsneg": "massPlotBEMuonIDCSNegNoLog",	
	
}

lumis = {

	"Mu": 42.1,
	"Ele": 41.5

} 
if args.do2016:
	lumis = {"Mu":36.3, "Ele":35.9}   
	
etabins = ["bb","be"]
#             0    1    2    3
csbins = ["inc","csneg","cspos"]
#            0     1     2

csbin  = args.cs
unc    = args.unc
debug = args.d

if csbin not in csbins:
	print("CS bin '{0}' not in:".format(csbin),csbins)
	exit(1)
if unc not in uncertainties:
	print("Plot type '{0}' not in:".format(unc),uncertainties)
	exit(1)


outf = "parametrizations/"
if args.truncation:
	outf = "parametrizationsTruncation/"
if not os.path.exists(outf):
    os.makedirs(outf)
for etabin in etabins:

	for antype in antypes:
		muonlyuncs = ["muonid", "smeared"]
		eleonlyuncs = ["piledown", "pileup",'prefireup','prefiredown']
		if unc in muonlyuncs and antype[2] == "Ele":
			print("Not processing uncertainty '{0:s}' for lepton flavour '{1:s}'".format(unc,antype[2]))
			continue
		if unc in eleonlyuncs and antype[2] == "Mu":
			print("Not processing uncertainty '{0:s}' for lepton flavour '{1:s}'".format(unc,antype[2]))
			continue
		params = {}
		addLabel = ""
		if args.do2016:
			addLabel="_2016"
		elif args.do2018:
			addLabel="_2018"
		model = "CI"
		if args.add: model = "ADD"
		with open(outf+"{0:s}parametrization_2{1:s}_{2:s}_{3:s}_{4:s}{5:s}.json".format(model,antype[1],unc,etabin,csbin,addLabel),"w") as js:
			with open(outf+"{0:s}counts_2{1:s}_{2:s}_{3:s}_{4:s}{5:s}.txt".format(model,antype[1],unc,etabin,csbin,addLabel),"w") as out:
				for intf in intfs:
					for heli in helis:
						files=[]
						for point in supers[:-1]:
							# params["{0:s}{1:s}_{2:d}GeV".format(intf,heli,point)] = np.zeros(len(lvals),'float64')
							params["{0:s}{1:s}_{2:d}GeV".format(intf,heli,point)]     = [0. for j in range(len(lvals))]
							params["{0:s}{1:s}_{2:d}GeV_err".format(intf,heli,point)] = [0. for j in range(len(lvals))]
							pass
						for point in extrabins:
							params["singleBin_{0:s}{1:s}_{2:d}GeV".format(intf,heli,point)]     = [0. for j in range(len(lvals))]
							params["singleBin_{0:s}{1:s}_{2:d}GeV_err".format(intf,heli,point)] = [0. for j in range(len(lvals))]
							pass
						for i,lval in enumerate(lvals):
							hist = None
							can   = r.TCanvas("can","",800,800)
							stack = r.THStack("stack","")
							plotName = antype[2]+etabin+unc
							if not csbin == "inc":
								plotName = antype[2]+etabin+unc+csbin
							plot = getPlot(plots[plotName])

	 
							eventCounts = totalNumberOfGeneratedEvents(path,plot.muon)  
							negWeights = negWeightFractions(path,plot.muon)
							if args.do2016:	
								lumi = 35.9*1000
								if plot.muon:
									lumi = 36294.6
							elif args.do2018:	
								lumi = 59.4*1000
								if plot.muon:
									lumi = 61298.775231718995
							else:
								lumi = 41.529*1000
								if plot.muon:
									lumi = 42079.880396
							if args.do2016:		
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
							elif args.do2018:		
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
							signal = "%sTo2%s_Lam%sTeV%s%s"%(model,antype[0],lval,intf,heli)
							if args.add:
									if lval == "100k":
										signal = "ADDGravTo2%s_Lam100k"%(antype[0])
									else:	
										signal = "ADDGravTo2%s_Lam%s"%(antype[0],str(int(float(lval)*1000)))


							
								
							if signal == "CITo2E_Lam40TeVDesRR" or signal == "CITo2Mu_Lam40TeVConRR":
								# list for 2017
								continue
							if signal == 'CITo2E_Lam100kTeVDesRR':
								signal = 'CITo2E_Lam100kTeVDesLL'	
							if args.do2016 and (signal == "CITo2Mu_Lam10TeVConRL" or signal == "CITo2Mu_Lam10TeVConLR"):
								# list for 2016
								continue
								
							if '100k' in signal and not args.add:			
								if "E" in signal:	
									signalConLL = "CITo2E_Lam100kTeVConLL"	
									signalConLR = "CITo2E_Lam100kTeVConLR"	
									signalConRL = "CITo2E_Lam100kTeVConRL"	
									signalConRR = "CITo2E_Lam100kTeVConRR"	
									signalDesLL = "CITo2E_Lam100kTeVDesLL"	
									signalDesLR = "CITo2E_Lam100kTeVDesLR"	
									signalDesRL = "CITo2E_Lam100kTeVDesRL"	
									signalDesRR = "CITo2E_Lam100kTeVDesRR"	
								else:	
									signalConLL = "CITo2Mu_Lam100kTeVConLL"	
									signalConLR = "CITo2Mu_Lam100kTeVConLR"	
									signalConRL = "CITo2Mu_Lam100kTeVConRL"	
									signalConRR = "CITo2Mu_Lam100kTeVConRR"	
									signalDesLL = "CITo2Mu_Lam100kTeVDesLL"	
									signalDesLR = "CITo2Mu_Lam100kTeVDesLR"	
									signalDesRL = "CITo2Mu_Lam100kTeVDesRL"	
									signalDesRR = "CITo2Mu_Lam100kTeVDesRR"	
									
								if args.do2016:	
									SignalConLL = Process(getattr(Signals2016,signalConLL),eventCounts,negWeights) 
									SignalConLR = Process(getattr(Signals2016,signalConLR),eventCounts,negWeights) 
									SignalConRL = Process(getattr(Signals2016,signalConRL),eventCounts,negWeights) 
									SignalConRR = Process(getattr(Signals2016,signalConRR),eventCounts,negWeights) 
									SignalDesLL = Process(getattr(Signals2016,signalDesLL),eventCounts,negWeights) 
									SignalDesLR = Process(getattr(Signals2016,signalDesLR),eventCounts,negWeights) 
									SignalDesRL = Process(getattr(Signals2016,signalDesRL),eventCounts,negWeights) 
									SignalDesRR = Process(getattr(Signals2016,signalDesRR),eventCounts,negWeights) 
								elif args.do2018:	
									SignalConLL = Process(getattr(Signals2018,signalConLL),eventCounts,negWeights) 
									SignalConLR = Process(getattr(Signals2018,signalConLR),eventCounts,negWeights) 
									SignalConRL = Process(getattr(Signals2018,signalConRL),eventCounts,negWeights) 
									SignalConRR = Process(getattr(Signals2018,signalConRR),eventCounts,negWeights) 
									SignalDesLL = Process(getattr(Signals2018,signalDesLL),eventCounts,negWeights) 
									SignalDesLR = Process(getattr(Signals2018,signalDesLR),eventCounts,negWeights) 
									SignalDesRL = Process(getattr(Signals2018,signalDesRL),eventCounts,negWeights) 
								else:
									SignalConLL = Process(getattr(Signals,signalConLL),eventCounts,negWeights) 
									SignalConLR = Process(getattr(Signals,signalConLR),eventCounts,negWeights) 
									SignalConRL = Process(getattr(Signals,signalConRL),eventCounts,negWeights) 
									SignalConRR = Process(getattr(Signals,signalConRR),eventCounts,negWeights) 
									SignalDesLL = Process(getattr(Signals,signalDesLL),eventCounts,negWeights) 
									SignalDesLR = Process(getattr(Signals,signalDesLR),eventCounts,negWeights) 
									SignalDesRL = Process(getattr(Signals,signalDesRL),eventCounts,negWeights) 
								
								signalhist = SignalConLL.loadHistogram(plot,lumi,zScaleFac)
								signalhist.Add(SignalConLR.loadHistogram(plot,lumi,zScaleFac))
								signalhist.Add(SignalConRL.loadHistogram(plot,lumi,zScaleFac))
								signalhist.Add(SignalConRR.loadHistogram(plot,lumi,zScaleFac))
								signalhist.Add(SignalDesLL.loadHistogram(plot,lumi,zScaleFac))
								signalhist.Add(SignalDesLR.loadHistogram(plot,lumi,zScaleFac))
								signalhist.Add(SignalDesRL.loadHistogram(plot,lumi,zScaleFac))
								if args.do2016:
									signalhist.Add(SignalDesRR.loadHistogram(plot,lumi,zScaleFac))
									signalhist.Scale(1./8)
								else:	
									signalhist.Scale(1./7)
								signalhist.Rebin(20)
								signalHistBefore = signalhist.Clone()
								if not args.do2016:
									if unc == "pdfWeightsUp":
										signalhist = applyPDFCorrection(signalhist)
										signalhist = applyPDFCorrection(signalhist)
									elif unc == "pdfWeightsDown":
										print ("not applying weights for down uncertainty")
									else:
										signalhist = applyPDFCorrection(signalhist)
							else:		
								if args.do2016 and not args.add:	
									Signal = Process(getattr(Signals2016,signal),eventCounts,negWeights) 
								elif args.do2018 and not args.add:	
									Signal = Process(getattr(Signals2018,signal),eventCounts,negWeights) 
								elif args.add and args.do2016:
										Signal = Process(getattr(Signals2016ADD, signal),eventCounts,negWeights)
								elif args.add and args.do2018:
									Signal = Process(getattr(Signals2018ADD, signal),eventCounts,negWeights)
								elif args.add:
									Signal = Process(getattr(SignalsADD, signal),eventCounts,negWeights)
								else:	
									Signal = Process(getattr(Signals,signal),eventCounts,negWeights)                        
								
								signalhist = Signal.loadHistogram(plot,lumi,zScaleFac)
								signalhist.Rebin(20)
								signalHistBefore = signalhist.Clone()
								if not args.do2016:
									if unc == "pdfWeightsUp":
										signalhist = applyPDFCorrection(signalhist)
										signalhist = applyPDFCorrection(signalhist)
									elif unc == "pdfWeightsDown":
										print ("not applying weights for down uncertainty")
									else:
										signalhist = applyPDFCorrection(signalhist)
							

							signalhist.SetMinimum(0.8*signalhist.GetMinimum(0.001))
							signalhist.SetMaximum(1.25*signalhist.GetMaximum())
							signalhist.Draw("hist")
							signalHistBefore.SetLineStyle(ROOT.kDashed)
							signalHistBefore.SetLineColor(ROOT.kRed)
							signalHistBefore.Draw("samehist")

							signalhist.GetXaxis().SetRangeUser(0,5000)
							signalhist.GetXaxis().SetNdivisions(505)
							r.gPad.SetLogy(True)

							### truncate ADD samples
							
							if args.add:
								if args.do2016:
									if not lval == "100k":
										signal100k = "ADDGravTo2%s_Lam100k"%(antype[0])
										Signal100k = Process(getattr(Signals2016ADD, signal100k),eventCounts,negWeights)
										signalhist100k = Signal100k.loadHistogram(plot,lumi,zScaleFac)
										signalhist100k.Rebin(20)
										signalhistTruncated = truncateADD(signalhist, signalhist100k,int(float(lval)*1000))
								elif args.do2018:
									if not lval == "100":
										signal100 = "ADDGravTo2%s_Lam%s"%(antype[0],100000)
										Signal100 = Process(getattr(Signals2018ADD, signal100),eventCounts,negWeights)
										signalhist100 = Signal100.loadHistogram(plot,lumi,zScaleFac)
										signalhist100.Rebin(20)
										signalhistTruncated = truncateADD(signalhist, signalhist100,int(float(lval)*1000))
								else:
									if not lval == "100":
										signal100 = "ADDGravTo2%s_Lam%s"%(antype[0],100000)
										Signal100 = Process(getattr(SignalsADD, signal100),eventCounts,negWeights)
										signalhist100 = Signal100.loadHistogram(plot,lumi,zScaleFac)
										signalhist100.Rebin(20)
										signalhistTruncated = truncateADD(signalhist, signalhist100,int(float(lval)*1000))
	 	
								signalhistTruncated.SetLineStyle(ROOT.kDotted)
								signalhistTruncated.Draw("samehist")
								if not lval == "100k":
									signalhist.GetXaxis().SetRangeUser(0,max(5000,int(float(lval)*1000)+1000))
									signalhist.SetMinimum(1e-5)
								
							latex = ROOT.TLatex()
							latex.SetTextFont(42)
							latex.SetTextAlign(31)
							latex.SetTextSize(0.04)
							latex.SetNDC(True)						

							latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV)"%(lumis[antype[2]],))
							
							
							latex.DrawLatex(0.85, 0.8, "%s"%(signal))
							
							can.Update()
							# raw_input()
							for ftype in ["png","C","pdf","eps"]:
								if args.do2016:
									can.SaveAs(outDir+"/{8:s}to2{1:s}_{2:s}{3:s}{4:s}_{5:s}_{6:s}_{7:s}_2016.{0:s}".format(ftype,antype[1],lval,intf,heli,unc,etabin,csbin,model))									
								elif args.do2018:
									can.SaveAs(outDir+"/{8:s}to2{1:s}_{2:s}{3:s}{4:s}_{5:s}_{6:s}_{7:s}_2018.{0:s}".format(ftype,antype[1],lval,intf,heli,unc,etabin,csbin,model))									
								else:	
									can.SaveAs(outDir+"/{8:s}to2{1:s}_{2:s}{3:s}{4:s}_{5:s}_{6:s}_{7:s}.{0:s}".format(ftype,antype[1],lval,intf,heli,unc,etabin,csbin,model))
							# raw_input("enter to continue")
							
							
							if args.add and args.truncation:
								signalhist = signalhistTruncated						
							
							
							for p,point in enumerate(supers[:-1]):
								bval  = signalhist.FindBin(point)
								upval = signalhist.FindBin(supers[p+1]-0.05)
								val   = signalhist.Integral(bval,upval)
								err   = r.Double(0)
								val   = signalhist.Integral(bval,upval)
								val2  = signalhist.IntegralAndError(bval,upval,err)
								out.write("{0:s} {1:d} {2:d} {3:d} {4:2.4f} {5:2.4f}\n".format(lval,point,bval,upval,val,err))
								params["{0:s}{1:s}_{2:d}GeV".format(intf,heli,point)][i] = val
								params["{0:s}{1:s}_{2:d}GeV_err".format(intf,heli,point)][i] = err
								pass

							# Mass bin scan above 1 TeV
							for point in extrabins:
								bval  = signalhist.FindBin(point)
								upval = signalhist.FindBin(100000000)
								val   = signalhist.Integral(bval,upval)
								err   = r.Double(0)
								val   = signalhist.Integral(bval,upval)
								val2  = signalhist.IntegralAndError(bval,upval,err)
								out.write("SingleBin {0:s} {1:d} {2:d} {3:d} {4:2.4f} {5:2.4f}\n".format(lval,point,bval,upval,val,err))
								params["singleBin_{0:s}{1:s}_{2:d}GeV".format(intf,heli,point)][i]     = val
								params["singleBin_{0:s}{1:s}_{2:d}GeV_err".format(intf,heli,point)][i] = err
								pass
							pass
						pass
					pass
				pass
			json.dump(params,js)
			pass
		pass
	pass
