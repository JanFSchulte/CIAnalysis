from ROOT import gROOT, TFile, TGraphErrors, TCanvas
from numpy import array as ar
from array import array
import argparse
from copy import deepcopy, copy
import pickle

from defs import getPlot, Backgrounds, Signals, Data, path, Signals2016, Signals2016ADD, SignalsADD, Signals2018, Signals2018ADD, zScale2016, zScale2018, Backgrounds, Backgrounds2016, Backgrounds2018
from helpers import *


def applyPDFCorrection(hist):
	
	for i in range(0,hist.GetNbinsX()+1):
		binCenter = hist.GetBinCenter(i)
		scaleFac = 0.86 - 3.72e-05 * binCenter + 2.72e-08 * binCenter **2
		hist.SetBinContent(i,hist.GetBinContent(i)*scaleFac)
	return copy(hist)


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
	"MubbpdfUncertainty": "massPlotBBNoLog",
	"MubepdfUncertainty": "massPlotBENoLog",
	"ElebbpdfUncertainty": "massPlotEleBBNoLog",
	"ElebepdfUncertainty": "massPlotEleBENoLog",
	
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
	
	"MubbpdfUncertaintycspos": "massPlotBBCSPosNoLog",
	"MubepdfUncertaintycspos": "massPlotBECSPosNoLog",
	"ElebbpdfUncertaintycspos": "massPlotEleBBCSPosNoLog",
	"ElebepdfUncertaintycspos": "massPlotEleBECSPosNoLog",

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
	"MubbpdfUncertaintycsneg": "massPlotBBCSNegNoLog",
	"MubepdfUncertaintycsneg": "massPlotBECSNegNoLog",
	"ElebbpdfUncertaintycsneg": "massPlotEleBBCSNegNoLog",
	"ElebepdfUncertaintycsneg": "massPlotEleBECSNegNoLog",
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


def main():
	gROOT.SetBatch(True)
	
	parser = argparse.ArgumentParser(description='Process some integers.')
	
	parser.add_argument("-add", "--add", action="store_true", dest="useADD", default=False,
						  help="use ADD instead of CI.")
	parser.add_argument("-s", "--suffix", dest="suffix", default='nominal',
						  help="name of systematic to use")
	args = parser.parse_args()					  
	useADD = args.useADD					  
	histos = ["BB","BE"]
	labels = ["dielectron_2016","dimuon_2016","dimuon_2017","dielectron_2017","dimuon_2018","dielectron_2018"]
	suffixesMu = ["nominal","scaleup","scaledown","smeared","muonid","pdfWeightsUp","pdfWeightsDown","pdfUncertainty"]
	suffixesEle = ["nominal","scaledown","scaleup","pileup","piledown","pdfWeightsUp","pdfWeightsDown",'prefireup','prefiredown',"pdfUncertainty"]
	css = ["inc","cspos","csneg"]
	lambdas = [10,16,22,28,34,40,46]
	interferences = ["Con","Des"]
	hels = ["LL","RL","LR","RR"]
	massBins = [400,500,700,1100,1900,3500]
	if useADD:
		labels = ["dielectron_2016","dimuon_2016","dimuon_2017","dielectron_2017","dimuon_2018","dielectron_2018"]
		lambdas = [3500+i*500 for i in range(12)]; lambdas.append(10000)
		interferences = [""]
		hels = [""]
		massBins = [1800, 2200, 2600, 3000, 3400]
		
	uncertainties = {}		
	for label in labels:
		print (label)
		if "dimuon" in label:
			suffixes = suffixesMu
		else:
			suffixes = suffixesEle
		if not args.suffix in suffixes: continue	
		suffix = args.suffix
		for cs in css:
			for histo in histos:
				for hel in hels:
					for interference in interferences:			
						model = interference+hel
						addci = "CI"
						if useADD: addci = "ADD"
						if "dimuon" in label:
							name = "%sto2mu"%addci
						else:			#~ print signalYields

							name = "%sto2e"	%addci
						
						if "dimuon" in label:
							plotName = "Mu"+histo.lower()+suffix
							plotNameDefault = "Mu"+histo.lower()+"nominal"
							antype = ["Mu"]
						else:	
							plotName = "Ele"+histo.lower()+suffix
							plotNameDefault = "Ele"+histo.lower()+"nominal"
							antype = ["E"]

						if not cs == "inc":
							plotName = plotName + cs	
							plotNameDefault = plotNameDefault + cs	
						plot = getPlot(plots[plotName])
						plotDefault = getPlot(plots[plotNameDefault])

 
						eventCounts = totalNumberOfGeneratedEvents(path,plot.muon)  
						negWeights = negWeightFractions(path,plot.muon)
						if "2016" in label:	
							lumi = 35.9*1000
							if plot.muon:
								lumi = 36294.6
						elif "2018" in label:	
							lumi = 59.4*1000
							if plot.muon:
								lumi = 61298.775231718995
						else:
							lumi = 41.529*1000
							if plot.muon:
								lumi = 42079.880396
						if "2016" in label:		
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
						elif "2018" in label:		
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
						background = "DrellYan"
						if "2016" in label:	
							DY = Process(getattr(Backgrounds2016,background),eventCounts,negWeights) 
						elif "2018" in label:	
							DY = Process(getattr(Backgrounds2018,background),eventCounts,negWeights) 
						else:	
							DY = Process(getattr(Backgrounds,background),eventCounts,negWeights)						
						
						backgroundHist = DY.loadHistogram(plot,lumi,zScaleFac)	
						backgroundHist = backgroundHist.Rebin(len(massBins), name, array('d', massBins+[10000]))
						backgroundHistDefault = DY.loadHistogram(plotDefault,lumi,zScaleFac)	
						backgroundHistDefault = backgroundHistDefault.Rebin(len(massBins), name, array('d', massBins+[10000]))


						if "2016" in label:
							
							
							if args.useADD:
							
								signal3500TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"3500")
								signal4000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"4000")							
								signal4500TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"4500")
								signal5000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"5000")
								signal5500TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"5500")
								signal6000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"6000")
								signal6500TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"6500")
								signal7000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"7000")
								signal7500TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"7500")
								signal8000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"8000")
								signal8500TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"8500")
								signal9000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"9000")
								signal10000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"10000")

								Signal3500TeV = Process(getattr(Signals2016ADD,signal3500TeV),eventCounts,negWeights) 
								Signal4000TeV = Process(getattr(Signals2016ADD,signal4000TeV),eventCounts,negWeights) 
								Signal4500TeV = Process(getattr(Signals2016ADD,signal4500TeV),eventCounts,negWeights) 
								Signal5000TeV = Process(getattr(Signals2016ADD,signal5000TeV),eventCounts,negWeights) 
								Signal5500TeV = Process(getattr(Signals2016ADD,signal5500TeV),eventCounts,negWeights) 
								Signal6000TeV = Process(getattr(Signals2016ADD,signal6000TeV),eventCounts,negWeights) 
								Signal6500TeV = Process(getattr(Signals2016ADD,signal6500TeV),eventCounts,negWeights) 
								Signal7000TeV = Process(getattr(Signals2016ADD,signal7000TeV),eventCounts,negWeights) 
								Signal7500TeV = Process(getattr(Signals2016ADD,signal7500TeV),eventCounts,negWeights) 
								Signal8000TeV = Process(getattr(Signals2016ADD,signal8000TeV),eventCounts,negWeights) 
								Signal8500TeV = Process(getattr(Signals2016ADD,signal8500TeV),eventCounts,negWeights) 
								Signal9000TeV = Process(getattr(Signals2016ADD,signal9000TeV),eventCounts,negWeights) 
								Signal10000TeV = Process(getattr(Signals2016ADD,signal10000TeV),eventCounts,negWeights) 

								signalhist3500TeV = Signal3500TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist4000TeV = Signal4000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist4500TeV = Signal4500TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist5000TeV = Signal5000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist5500TeV = Signal5500TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist6000TeV = Signal6000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist6500TeV = Signal6500TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist7000TeV = Signal7000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist7500TeV = Signal7500TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist8000TeV = Signal8000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist8500TeV = Signal8500TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist9000TeV = Signal9000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist10000TeV = Signal10000TeV.loadHistogram(plot,lumi,zScaleFac)
						
								
								signalhist3500TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist4000TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist4500TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist5000TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist5500TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist6000TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist6500TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist7000TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist7500TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist8000TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist8500TeV.Add(signalhist10000TeV.Clone(),-1)
								signalhist9000TeV.Add(signalhist10000TeV.Clone(),-1)

								
								signalhistDefault3500TeV = Signal3500TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault4000TeV = Signal4000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault4500TeV = Signal4500TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault5000TeV = Signal5000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault5500TeV = Signal5500TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault6000TeV = Signal6000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault6500TeV = Signal6500TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault7000TeV = Signal7000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault7500TeV = Signal7500TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault8000TeV = Signal8000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault8500TeV = Signal8500TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault9000TeV = Signal9000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault10000TeV = Signal10000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
						
								
								signalhistDefault3500TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault4000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault4500TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault5000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault5500TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault6000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault6500TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault7000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault7500TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault8000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault8500TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault9000TeV.Add(signalhistDefault10000TeV.Clone(),-1)


								for i in range(0,signalhist3500TeV.GetNbinsX()+1):
									signalhist3500TeV.SetBinContent(i,(signalhist3500TeV.GetBinContent(i) + signalhist4000TeV.GetBinContent(i) + signalhist4500TeV.GetBinContent(i) + signalhist5000TeV.GetBinContent(i) + signalhist5500TeV.GetBinContent(i) + signalhist6000TeV.GetBinContent(i) + signalhist6500TeV.GetBinContent(i) + signalhist7000TeV.GetBinContent(i) + signalhist7500TeV.GetBinContent(i) + signalhist8000TeV.GetBinContent(i) + signalhist8500TeV.GetBinContent(i) + signalhist9000TeV.GetBinContent(i))  / 12)
									signalhistDefault3500TeV.SetBinContent(i,(signalhistDefault3500TeV.GetBinContent(i) + signalhistDefault4000TeV.GetBinContent(i) + signalhistDefault4500TeV.GetBinContent(i) + signalhistDefault5000TeV.GetBinContent(i) + signalhistDefault5500TeV.GetBinContent(i) + signalhistDefault6000TeV.GetBinContent(i) + signalhistDefault6500TeV.GetBinContent(i) + signalhistDefault7000TeV.GetBinContent(i) + signalhistDefault7500TeV.GetBinContent(i) + signalhistDefault8000TeV.GetBinContent(i) + signalhistDefault8500TeV.GetBinContent(i) + signalhistDefault9000TeV.GetBinContent(i)) / 12)
										
									
									
									
								signalHist = signalhist3500TeV
								signalHistDefault = signalhistDefault3500TeV
							else:	
								signal1TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"1",interference,hel)
								signal10TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"10",interference,hel)							
								signal16TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"16",interference,hel)
								signal22TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"22",interference,hel)
								signal28TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"28",interference,hel)
								signal34TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"34",interference,hel)
								signal100kTeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"100k",interference,hel)

								Signal1TeV = Process(getattr(Signals2016,signal1TeV),eventCounts,negWeights) 
								Signal16TeV = Process(getattr(Signals2016,signal16TeV),eventCounts,negWeights) 
								Signal22TeV = Process(getattr(Signals2016,signal22TeV),eventCounts,negWeights) 
								Signal28TeV = Process(getattr(Signals2016,signal28TeV),eventCounts,negWeights) 
								Signal34TeV = Process(getattr(Signals2016,signal34TeV),eventCounts,negWeights) 
								Signal100kTeV = Process(getattr(Signals2016,signal100kTeV),eventCounts,negWeights) 

								signalhist1TeV = Signal1TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist16TeV = Signal16TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist22TeV = Signal22TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist28TeV = Signal28TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist34TeV = Signal34TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist100kTeV = Signal100kTeV.loadHistogram(plot,lumi,zScaleFac)
						
								
								signalhist1TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist16TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist22TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist28TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist34TeV.Add(signalhist100kTeV.Clone(),-1)

								
								signalhistDefault1TeV = Signal1TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault16TeV = Signal16TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault22TeV = Signal22TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault28TeV = Signal28TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault34TeV = Signal34TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault100kTeV = Signal100kTeV.loadHistogram(plotDefault,lumi,zScaleFac)
						
								
								signalhistDefault1TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault16TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault22TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault28TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault34TeV.Add(signalhistDefault100kTeV.Clone(),-1)




								if not "ConRL" in signal10TeV and not "ConLR" in signal10TeV:
									Signal10TeV = Process(getattr(Signals2016,signal10TeV),eventCounts,negWeights) 
									signalhist10TeV = Signal10TeV.loadHistogram(plot,lumi,zScaleFac)
									signalhist10TeV.Add(signalhist100kTeV.Clone(),-1)
									signalhist1TeV.Add(signalhist10TeV.Clone())
			
									signalhistDefault10TeV = Signal10TeV.loadHistogram(plotDefault,lumi,zScaleFac)
									signalhistDefault10TeV.Add(signalhistDefault100kTeV.Clone(),-1)
									signalhistDefault1TeV.Add(signalhistDefault10TeV.Clone())
									
									for i in range(0,signalhist1TeV.GetNbinsX()+1):
										signalhist1TeV.SetBinContent(i,(signalhist1TeV.GetBinContent(i) + signalhist10TeV.GetBinContent(i) + signalhist16TeV.GetBinContent(i) + signalhist22TeV.GetBinContent(i) + signalhist28TeV.GetBinContent(i) + signalhist34TeV.GetBinContent(i)) / 6)
										signalhistDefault1TeV.SetBinContent(i,(signalhistDefault1TeV.GetBinContent(i) + signalhistDefault10TeV.GetBinContent(i) + signalhistDefault16TeV.GetBinContent(i) + signalhistDefault22TeV.GetBinContent(i) + signalhistDefault28TeV.GetBinContent(i) + signalhistDefault34TeV.GetBinContent(i)) / 6)

								else:

									for i in range(0,signalhist1TeV.GetNbinsX()+1):
										signalhist1TeV.SetBinContent(i,(signalhist1TeV.GetBinContent(i) + signalhist16TeV.GetBinContent(i) + signalhist22TeV.GetBinContent(i) + signalhist28TeV.GetBinContent(i) + signalhist34TeV.GetBinContent(i)) / 5)
										signalhistDefault1TeV.SetBinContent(i,(signalhistDefault1TeV.GetBinContent(i) + signalhistDefault16TeV.GetBinContent(i) + signalhistDefault22TeV.GetBinContent(i) + signalhistDefault28TeV.GetBinContent(i) + signalhistDefault34TeV.GetBinContent(i)) / 5)
										
									
									
									
								signalHist = signalhist1TeV
								signalHistDefault = signalhistDefault1TeV

						if "2017" in label:
							
							if args.useADD:
							
								signal4000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"4000")
								signal5000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"5000")							
								signal6000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"6000")
								signal7000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"7000")
								signal8000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"8000")
								signal9000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"9000")
								signal10000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"10000")
								signal11000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"11000")
								signal12000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"12000")
								signal100000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"100000")

								Signal4000TeV = Process(getattr(SignalsADD,signal4000TeV),eventCounts,negWeights) 
								Signal5000TeV = Process(getattr(SignalsADD,signal5000TeV),eventCounts,negWeights) 
								Signal6000TeV = Process(getattr(SignalsADD,signal6000TeV),eventCounts,negWeights) 
								Signal7000TeV = Process(getattr(SignalsADD,signal7000TeV),eventCounts,negWeights) 
								Signal8000TeV = Process(getattr(SignalsADD,signal8000TeV),eventCounts,negWeights) 
								Signal9000TeV = Process(getattr(SignalsADD,signal9000TeV),eventCounts,negWeights) 
								Signal10000TeV = Process(getattr(SignalsADD,signal10000TeV),eventCounts,negWeights) 
								Signal11000TeV = Process(getattr(SignalsADD,signal11000TeV),eventCounts,negWeights) 
								Signal12000TeV = Process(getattr(SignalsADD,signal12000TeV),eventCounts,negWeights) 
								Signal100000TeV = Process(getattr(SignalsADD,signal100000TeV),eventCounts,negWeights) 

								signalhist4000TeV = Signal4000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist5000TeV = Signal5000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist6000TeV = Signal6000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist7000TeV = Signal7000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist8000TeV = Signal8000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist9000TeV = Signal9000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist10000TeV = Signal10000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist11000TeV = Signal11000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist12000TeV = Signal12000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist100000TeV = Signal100000TeV.loadHistogram(plot,lumi,zScaleFac)
						
								
								signalhist4000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist5000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist6000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist7000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist8000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist9000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist10000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist11000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist12000TeV.Add(signalhist100000TeV.Clone(),-1)

								
								signalhistDefault4000TeV = Signal4000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault5000TeV = Signal5000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault6000TeV = Signal6000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault7000TeV = Signal7000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault8000TeV = Signal8000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault9000TeV = Signal9000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault10000TeV = Signal10000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault11000TeV = Signal11000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault12000TeV = Signal12000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault100000TeV = Signal100000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
						
								
								signalhistDefault4000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault5000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault6000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault7000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault8000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault9000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault10000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault11000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault12000TeV.Add(signalhistDefault10000TeV.Clone(),-1)


								for i in range(0,signalhist4000TeV.GetNbinsX()+1):
									signalhist4000TeV.SetBinContent(i,(signalhist4000TeV.GetBinContent(i) + signalhist5000TeV.GetBinContent(i) + signalhist6000TeV.GetBinContent(i) + signalhist7000TeV.GetBinContent(i) + signalhist8000TeV.GetBinContent(i) + signalhist9000TeV.GetBinContent(i) + signalhist10000TeV.GetBinContent(i) + signalhist11000TeV.GetBinContent(i) + signalhist12000TeV.GetBinContent(i) + signalhist100000TeV.GetBinContent(i)) / 9)
									signalhistDefault4000TeV.SetBinContent(i,(signalhistDefault4000TeV.GetBinContent(i) + signalhistDefault5000TeV.GetBinContent(i) + signalhistDefault6000TeV.GetBinContent(i) + signalhistDefault7000TeV.GetBinContent(i) + signalhistDefault8000TeV.GetBinContent(i) + signalhistDefault9000TeV.GetBinContent(i) + signalhistDefault10000TeV.GetBinContent(i) + signalhistDefault11000TeV.GetBinContent(i) + signalhistDefault12000TeV.GetBinContent(i) + signalhistDefault100000TeV.GetBinContent(i)) / 9)
										
									
									
									
								signalHist = signalhist4000TeV
								signalHistDefault = signalhistDefault4000TeV							
							
							else:
								
								signal16TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"16",interference,hel)
								signal24TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"24",interference,hel)
								signal32TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"32",interference,hel)
								signal40TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"40",interference,hel)
								signal100kTeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"100k",interference,hel)

								if signal100kTeV == 'CITo2E_Lam100kTeVDesRR':
									signal100kTeV = 'CITo2E_Lam100kTeVDesLL'	


								Signal16TeV = Process(getattr(Signals,signal16TeV),eventCounts,negWeights) 
								Signal24TeV = Process(getattr(Signals,signal24TeV),eventCounts,negWeights) 
								Signal32TeV = Process(getattr(Signals,signal32TeV),eventCounts,negWeights) 
								Signal100kTeV = Process(getattr(Signals,signal100kTeV),eventCounts,negWeights) 

								signalhist16TeV = Signal16TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist24TeV = Signal24TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist32TeV = Signal32TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist100kTeV = Signal100kTeV.loadHistogram(plot,lumi,zScaleFac)
													
								signalhist16TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist24TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist32TeV.Add(signalhist100kTeV.Clone(),-1)
								
								signalhistDefault16TeV = Signal16TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault24TeV = Signal24TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault32TeV = Signal32TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault100kTeV = Signal100kTeV.loadHistogram(plotDefault,lumi,zScaleFac)
						
								signalhistDefault16TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault24TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault32TeV.Add(signalhistDefault100kTeV.Clone(),-1)



								if not "DesRR" in signal40TeV and not "ConRR" in signal40TeV:
									Signal40TeV = Process(getattr(Signals,signal40TeV),eventCounts,negWeights) 
									signalhist40TeV = Signal40TeV.loadHistogram(plot,lumi,zScaleFac)
									signalhist40TeV.Add(signalhist100kTeV.Clone(),-1)

			
									signalhistDefault40TeV = Signal40TeV.loadHistogram(plotDefault,lumi,zScaleFac)
									signalhistDefault40TeV.Add(signalhistDefault100kTeV.Clone(),-1)

									
									
									for i in range(0,signalhist1TeV.GetNbinsX()+1):
										signalhist16TeV.SetBinContent(i,(signalhist16TeV.GetBinContent(i) + signalhist24TeV.GetBinContent(i) + signalhist32TeV.GetBinContent(i) + signalhist40TeV.GetBinContent(i)) / 4)
										signalhistDefault16TeV.SetBinContent(i,(signalhistDefault16TeV.GetBinContent(i) + signalhistDefault24TeV.GetBinContent(i) + signalhistDefault32TeV.GetBinContent(i) + signalhistDefault40TeV.GetBinContent(i)) / 4)

								else:

									for i in range(0,signalhist1TeV.GetNbinsX()+1):
										signalhist16TeV.SetBinContent(i,(signalhist16TeV.GetBinContent(i) + signalhist24TeV.GetBinContent(i) + signalhist32TeV.GetBinContent(i)) / 3)
										signalhistDefault16TeV.SetBinContent(i,(signalhistDefault16TeV.GetBinContent(i) + signalhistDefault24TeV.GetBinContent(i) + signalhistDefault32TeV.GetBinContent(i)) / 3)
										
									
									
								signalHist = signalhist16TeV
								signalHistDefault = signalhistDefault16TeV								
						if "2018" in label:
							
							
							if args.useADD:
							
								signal4000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"4000")
								signal5000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"5000")							
								signal6000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"6000")
								signal7000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"7000")
								signal8000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"8000")
								signal9000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"9000")
								signal10000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"10000")
								signal11000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"11000")
								signal12000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"12000")
								signal100000TeV = "%sGravTo2%s_Lam%s"%(addci,antype[0],"100000")

								Signal4000TeV = Process(getattr(Signals2018ADD,signal4000TeV),eventCounts,negWeights) 
								Signal5000TeV = Process(getattr(Signals2018ADD,signal5000TeV),eventCounts,negWeights) 
								Signal6000TeV = Process(getattr(Signals2018ADD,signal6000TeV),eventCounts,negWeights) 
								Signal7000TeV = Process(getattr(Signals2018ADD,signal7000TeV),eventCounts,negWeights) 
								Signal8000TeV = Process(getattr(Signals2018ADD,signal8000TeV),eventCounts,negWeights) 
								Signal9000TeV = Process(getattr(Signals2018ADD,signal9000TeV),eventCounts,negWeights) 
								Signal10000TeV = Process(getattr(Signals2018ADD,signal10000TeV),eventCounts,negWeights) 
								Signal11000TeV = Process(getattr(Signals2018ADD,signal11000TeV),eventCounts,negWeights) 
								Signal12000TeV = Process(getattr(Signals2018ADD,signal12000TeV),eventCounts,negWeights) 
								Signal100000TeV = Process(getattr(Signals2018ADD,signal100000TeV),eventCounts,negWeights) 

								signalhist4000TeV = Signal4000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist5000TeV = Signal5000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist6000TeV = Signal6000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist7000TeV = Signal7000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist8000TeV = Signal8000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist9000TeV = Signal9000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist10000TeV = Signal10000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist11000TeV = Signal11000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist12000TeV = Signal12000TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist100000TeV = Signal100000TeV.loadHistogram(plot,lumi,zScaleFac)
						
								
								signalhist4000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist5000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist6000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist7000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist8000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist9000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist10000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist11000TeV.Add(signalhist100000TeV.Clone(),-1)
								signalhist12000TeV.Add(signalhist100000TeV.Clone(),-1)

								
								signalhistDefault4000TeV = Signal4000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault5000TeV = Signal5000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault6000TeV = Signal6000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault7000TeV = Signal7000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault8000TeV = Signal8000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault9000TeV = Signal9000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault10000TeV = Signal10000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault11000TeV = Signal11000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault12000TeV = Signal12000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault100000TeV = Signal100000TeV.loadHistogram(plotDefault,lumi,zScaleFac)
						
								
								signalhistDefault4000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault5000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault6000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault7000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault8000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault9000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault10000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault11000TeV.Add(signalhistDefault10000TeV.Clone(),-1)
								signalhistDefault12000TeV.Add(signalhistDefault10000TeV.Clone(),-1)


								for i in range(0,signalhist4000TeV.GetNbinsX()+1):
									signalhist4000TeV.SetBinContent(i,(signalhist4000TeV.GetBinContent(i) + signalhist5000TeV.GetBinContent(i) + signalhist6000TeV.GetBinContent(i) + signalhist7000TeV.GetBinContent(i) + signalhist8000TeV.GetBinContent(i) + signalhist9000TeV.GetBinContent(i) + signalhist10000TeV.GetBinContent(i) + signalhist11000TeV.GetBinContent(i) + signalhist12000TeV.GetBinContent(i) + signalhist100000TeV.GetBinContent(i)) / 9)
									signalhistDefault4000TeV.SetBinContent(i,(signalhistDefault4000TeV.GetBinContent(i) + signalhistDefault5000TeV.GetBinContent(i) + signalhistDefault6000TeV.GetBinContent(i) + signalhistDefault7000TeV.GetBinContent(i) + signalhistDefault8000TeV.GetBinContent(i) + signalhistDefault9000TeV.GetBinContent(i) + signalhistDefault10000TeV.GetBinContent(i) + signalhistDefault11000TeV.GetBinContent(i) + signalhistDefault12000TeV.GetBinContent(i) + signalhistDefault100000TeV.GetBinContent(i)) / 9)
										
									
									
									
								signalHist = signalhist4000TeV
								signalHistDefault = signalhistDefault4000TeV							
							else:
								signal16TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"16",interference,hel)
								signal24TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"24",interference,hel)
								signal32TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"32",interference,hel)
								signal40TeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"40",interference,hel)
								signal100kTeV = "%sTo2%s_Lam%sTeV%s%s"%(addci,antype[0],"100k",interference,hel)

								if signal100kTeV == 'CITo2E_Lam100kTeVDesRR':
									signal100kTeV = 'CITo2E_Lam100kTeVDesLL'	


								Signal16TeV = Process(getattr(Signals2018,signal16TeV),eventCounts,negWeights) 
								Signal24TeV = Process(getattr(Signals2018,signal24TeV),eventCounts,negWeights) 
								Signal32TeV = Process(getattr(Signals2018,signal32TeV),eventCounts,negWeights) 
								Signal100kTeV = Process(getattr(Signals2018,signal100kTeV),eventCounts,negWeights) 

								signalhist16TeV = Signal16TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist24TeV = Signal24TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist32TeV = Signal32TeV.loadHistogram(plot,lumi,zScaleFac)
								signalhist100kTeV = Signal100kTeV.loadHistogram(plot,lumi,zScaleFac)
						
								
								signalhist16TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist24TeV.Add(signalhist100kTeV.Clone(),-1)
								signalhist32TeV.Add(signalhist100kTeV.Clone(),-1)



								
								signalhistDefault16TeV = Signal16TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault24TeV = Signal24TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault32TeV = Signal32TeV.loadHistogram(plotDefault,lumi,zScaleFac)
								signalhistDefault100kTeV = Signal100kTeV.loadHistogram(plotDefault,lumi,zScaleFac)
						
								signalhistDefault16TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault24TeV.Add(signalhistDefault100kTeV.Clone(),-1)
								signalhistDefault32TeV.Add(signalhistDefault100kTeV.Clone(),-1)




								if not "DesRR" in signal40TeV and not "ConRR" in signal40TeV:
									Signal40TeV = Process(getattr(Signals2018,signal40TeV),eventCounts,negWeights) 
									signalhist40TeV = Signal40TeV.loadHistogram(plot,lumi,zScaleFac)
									signalhist40TeV.Add(signalhist100kTeV.Clone(),-1)
			
									signalhistDefault40TeV = Signal40TeV.loadHistogram(plotDefault,lumi,zScaleFac)
									signalhistDefault40TeV.Add(signalhistDefault100kTeV.Clone(),-1)

									for i in range(0,signalhist1TeV.GetNbinsX()+1):
										signalhist16TeV.SetBinContent(i,(signalhist16TeV.GetBinContent(i) + signalhist24TeV.GetBinContent(i) + signalhist32TeV.GetBinContent(i) + signalhist40TeV.GetBinContent(i)) / 4)
										signalhistDefault16TeV.SetBinContent(i,(signalhistDefault16TeV.GetBinContent(i) + signalhistDefault24TeV.GetBinContent(i) + signalhistDefault32TeV.GetBinContent(i) + signalhistDefault40TeV.GetBinContent(i)) / 4)

								else:

									for i in range(0,signalhist1TeV.GetNbinsX()+1):
										signalhist16TeV.SetBinContent(i,(signalhist16TeV.GetBinContent(i) + signalhist24TeV.GetBinContent(i) + signalhist32TeV.GetBinContent(i)) / 3)
										signalhistDefault16TeV.SetBinContent(i,(signalhistDefault16TeV.GetBinContent(i) + signalhistDefault24TeV.GetBinContent(i) + signalhistDefault32TeV.GetBinContent(i)) / 3)



								signalHist = signalhist16TeV
								signalHistDefault = signalhistDefault16TeV	

			
						signalHist = signalHist.Rebin(len(massBins), name, array('d', massBins+[10000]))
						signalHistDefault = signalHistDefault.Rebin(len(massBins), name, array('d', massBins+[10000]))
						
						
						if suffix == "pdfUncertainty":
								pdfUncert = [0.029808675319316736,0.032492944178882606,0.046693632068081733,0.06985339137281024,0.11520346316811705,0.17965523337854952]
								pdfUncertDY = [0.0133126577202494, 0.0147788328624159, 0.01842209757115, 0.0243644300786998, 0.039572839847093, 0.1144248733154136]
								for i in range(0,signalHist.GetNbinsX()):
									signalHist.SetBinContent(i,signalHist.GetBinContent(i)*(1.+pdfUncert[i]))
									backgroundHist.SetBinContent(i,backgroundHist.GetBinContent(i)*(1.+pdfUncertDY[i]))
							
						
						if suffix == "pdfWeightsUp":
							signalHist = applyPDFCorrection(signalHist)
							signalHist = applyPDFCorrection(signalHist)
							signalHistDefault = applyPDFCorrection(signalHistDefault)
						elif suffix == "pdfWeightsDown":
							signalHistDefault = applyPDFCorrection(signalHistDefault)							
							print ("not applying weights for down uncertainty")
						else:
							signalHist = applyPDFCorrection(signalHist)						
							signalHistDefault = applyPDFCorrection(signalHistDefault)
					
						signalHist.Add(backgroundHist)
						signalHistDefault.Add(backgroundHistDefault)
						
						if not cs == "inc":
							uncertainties["%s_%s_%s"%(cs,label,histo)] = {}
						else:	
							uncertainties["%s_%s"%(label,histo)] = {}
						
						
						
						for index, massBin in enumerate(massBins):
								# ~ print (abs(signalHist.GetBinContent(index+1) / signalHistDefault.GetBinContent(index+1) -1.))
								if not cs == "inc":
									uncertainties["%s_%s_%s"%(cs,label,histo)][str(index)] = 1. + abs(signalHist.GetBinContent(index+1) / signalHistDefault.GetBinContent(index+1) -1.)								
								else:
									uncertainties["%s_%s"%(label,histo)][str(index)] = 1. + abs(signalHist.GetBinContent(index+1) / signalHistDefault.GetBinContent(index+1) -1.)								
								

	fileName = "%ssystematicsForPriors"%addci
	
	if suffix == "nominal":
		otherSuffix = "default"
	elif suffix == "scaledown":
		otherSuffix = "scaleDown"
	elif suffix == "scaleup":
		otherSuffix = "scaleUp"
	elif suffix == "smeared":
		otherSuffix = "resolution"
	elif suffix == "muonid":
		otherSuffix = "ID"
	elif suffix == "pileup":
		otherSuffix = "pileup"
	elif suffix == "piledown":
		otherSuffix = "piledown"
	elif suffix == "prefireup":
		otherSuffix = "prefireup"
	elif suffix == "prefiredown":
		otherSuffix = "prefiredown"
	elif suffix == "pdfWeightsUp":
		otherSuffix = "pdfWeightsUp"
	elif suffix == "pdfWeightsDown":
		otherSuffix = "pdfWeightsDown"
	elif suffix == "pdfUncertainty":
		otherSuffix = "pdfUncertainty"
	else:
		print (suffix)
	outFilePkl = open("%s_%s.pkl"%(fileName,otherSuffix),"wb")
	pickle.dump(uncertainties, outFilePkl, protocol=2)
	outFilePkl.close()									
main()
