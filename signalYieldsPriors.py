from ROOT import gROOT, TFile
from numpy import array as ar
from array import array
import argparse
from copy import deepcopy
import pickle



def main():
	gROOT.SetBatch(True)
	
	parser = argparse.ArgumentParser(description='Process some integers.')
	
	parser.add_argument("-add", "--add", action="store_true", dest="useADD", default=False,
						  help="use ADD instead of CI.")
	parser.add_argument("-truncation", "--truncation", action="store_true", dest="truncation", default=False,
						  help="use ADD instead of CI.")
	parser.add_argument("-s", "--suffix", dest="suffix", default='nominal',
						  help="name of systematic to use")
	args = parser.parse_args()					  
	useADD = args.useADD					  
	histos = ["BB","BE"]
	labels = ["dielectron_2016","dimuon_2016","dimuon_2017","dielectron_2017","dimuon_2018","dielectron_2018"]
	suffixesMu = ["nominal","scaledown","smeared","muonid","pdfWeightsUp","pdfWeightsDown"]
	suffixesEle = ["nominal","scaledown","scaleup","pileup","piledown","pdfWeightsUp","pdfWeightsDown",'prefireup','prefiredown']
	lambdas = [10,16,22,28,34,40,46]
	interferences = ["Con","Des"]
	hels = ["LL","RL","LR","RR"]
	css = ["inc","cspos","csneg"]



	
	massBins = [400,500,700,1100,1900,3500]
	if useADD:
		labels = ["dielectron_2016","dimuon_2016","dimuon_2017","dielectron_2017","dimuon_2018","dielectron_2018"]
		lambdas = [3500+i*500 for i in range(12)]; lambdas.append(10000)
		interferences = [""]
		hels = [""]
		massBins = [1800, 2200, 2600, 3000, 3400]


	outDir = "parametrizations"
	if args.truncation:
		outDir = "parametrizationsTruncation"
		
	signalYields = {}	
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
						if args.useADD:
							if "2016" in label:		
								massBins = [1900, 2200, 2600, 3000, 3400]
							else:	
								massBins = [1800, 2200, 2600, 3000, 3400]
		
						model = interference+hel
						addci = "CI"
						if useADD: addci = "ADD"
						if "dimuon" in label:
							name = "%sto2mu"%addci
						else:			#~ print signalYields

							name = "%sto2e"	%addci
						if useADD:

							if "2016" in label:
								fitFile = TFile(outDir+"/%s_%s_%s_%s_parametrizationForPriors_fixinf_2016.root"%(name,suffix,histo.lower(),cs),"READ")
							elif "2018" in label:
								fitFile = TFile(outDir+"/%s_%s_%s_%s_parametrizationForPriors_fixinf_2018.root"%(name,suffix,histo.lower(),cs),"READ")
							else:
								fitFile = TFile(outDir+"/%s_%s_%s_%s_parametrizationForPriors_fixinf.root"%(name,suffix,histo.lower(),cs),"READ")
						else:	
							if "2016" in label:
								fitFile = TFile(outDir+"/%s_%s_%s_%s_parametrizationForPriors_fixinf_limitp0_limitp1_limitp2_2016.root"%(name,suffix,histo.lower(),cs),"READ")
							elif "2018" in label:
								fitFile = TFile(outDir+"/%s_%s_%s_%s_parametrizationForPriors_fixinf_limitp0_limitp1_limitp2_2018.root"%(name,suffix,histo.lower(),cs),"READ")
							else:	
								fitFile = TFile(outDir+"/%s_%s_%s_%s_parametrizationForPriors_fixinf_limitp0_limitp1_limitp2.root"%(name,suffix,histo.lower(),cs),"READ")
						# ~ print (fitFile.ls())		
						# ~ print ("%s_%s_%s_%s_parametrization_fixinf.root"%(name,suffix,histo.lower(),cs))
						for l in lambdas:
							if "dimuon" in label:
								name = "CITo2Mu_Lam%dTeV%s"%(l,model)
							else:	
								name = "CITo2E_Lam%dTeV%s"%(l,model)
							if useADD:
								name = "ADDGravTo2Mu_Lam%d"%l
								if "dielectron" in label: name = "ADDGravTo2E_Lam%d"%l
								l = l * 0.001
							if not cs == "inc":
								name += "_%s"%cs
							signalYields["%s_%s_%s"%(name,label,histo)] = {}
							for index, massBin in enumerate(massBins):
								function = fitFile.Get("fn_m%d_%s"%(massBin,model))
								fitR = fitFile.Get("fitR_m%d_%s"%(massBin,model))
								pars = fitR.GetParams()
								errs = fitR.Errors()
								function.SetParameter(0,pars[0])
								function.SetParameter(1,pars[1])
								function.SetParameter(2,pars[2])
								function.SetParError(0,errs[0])
								function.SetParError(1,errs[1])
								function.SetParError(2,errs[2])
								# ~ if useADD:
									# ~ function.SetParameter(3, pars[3])
									# ~ function.SetParError(3, errs[3])
								functionUnc = fitFile.Get("fn_unc_m%d_%s"%(massBin,model))
								uncert = (abs((functionUnc.Eval(l)/function.Eval(l))**2 + (functionUnc.Eval(100000)/function.Eval(100000))))**0.5	
								signalYields["%s_%s_%s"%(name,label,histo)][str(index)] = [(function.Eval(l)-function.Eval(100000)),uncert]
					if useADD: continue
					for l in lambdas:	
						if "dimuon" in label:
							name = "CITo2Mu_Lam%dTeV%s"%(l,hel)
							nameCon = "CITo2Mu_Lam%dTeV%s"%(l,"Con"+hel)
							nameDes = "CITo2Mu_Lam%dTeV%s"%(l,"Des"+hel)
						else:	
							name = "CITo2E_Lam%dTeV%s"%(l,hel)
							nameCon = "CITo2E_Lam%dTeV%s"%(l,"Con"+hel)
							nameDes = "CITo2E_Lam%dTeV%s"%(l,"Des"+hel)
						if not cs == "inc":
							name += "_%s"%cs							
							nameCon += "_%s"%cs							
							nameDes += "_%s"%cs							
						signalYields["%s_%s_%s"%(name,label,histo)] = {}	
						for index, massBin in enumerate(massBins):
							signalYields["%s_%s_%s"%(name,label,histo)][str(index)] = [(signalYields["%s_%s_%s"%(nameCon,label,histo)][str(index)][0] + signalYields["%s_%s_%s"%(nameDes,label,histo)][str(index)][0])/2,((signalYields["%s_%s_%s"%(nameCon,label,histo)][str(index)][1]**2 + signalYields["%s_%s_%s"%(nameDes,label,histo)][str(index)][1]**2))**0.5]
					

			addci = "CI"
			if useADD: addci = "ADD"
			if "dimuon" in label:
				fileName = "%ssignalYieldsPriors"%addci
			else:
				fileName = "%ssignalYieldsPriorsEle"%addci
			
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
			else:
				print (suffix)
			outFilePkl = open("%s_%s.pkl"%(fileName,otherSuffix),"wb")
			pickle.dump(signalYields, outFilePkl, protocol=2)
			outFilePkl.close()		
	
	
			
							
main()
