from ROOT import gROOT, TFile
from numpy import array as ar
from array import array
import argparse
from copy import deepcopy
import pickle

from sys import argv

def main():
	gROOT.SetBatch(True)

	parser = argparse.ArgumentParser(description='Process some integers.')
	
	parser.add_argument("-add", "--add", action="store_true", dest="useADD", default=False,
						  help="use ADD instead of CI.")
	parser.add_argument("-truncation", "--truncatio", action="store_true", dest="truncation", default=False,
						  help="use ADD instead of CI.")
	parser.add_argument("-s", "--suffix", dest="suffix", default='nominal',
						  help="name of systematic to use")
	args = parser.parse_args()					  
	useADD = args.useADD	

	suffix = args.suffix
	
	histos = ["BB","BE"]
	labels = ["dimuon_2016","dielectron_2016","dimuon_2017","dielectron_2017","dimuon_2018","dielectron_2018"]
	#~ channels = ["cito2mu","cito2e"]
	suffixesMu = ["nominal","scaledown","smeared","muonid","pdfWeightsUp","pdfWeightsDown"]
	suffixesEle = ["nominal","scaledown","scaleup","pileup","piledown","pdfWeightsUp","pdfWeightsDown",'prefireup','prefiredown']
	css = ["inc","cspos","csneg"]	
	#~ suffixes = ["smeared"]
	lambdas = [10,16,22,28,34,40,46,52,58]
	interferences = ["Con","Des"]
	# ~ hels = ["LL","LR","RR"]
	# ~ interferences = ["Con"]
	hels = ["LL","LR","RL","RR"]

	outf = "parametrizations/"
	if args.truncation:
			outf = "parametrizationsTruncation/"


	if useADD:
		labels = ["dielectron_2016","dimuon_2016","dimuon_2017","dielectron_2017","dimuon_2018","dielectron_2018"]
		lambdas = [3500+i*500 for i in range(12)]; lambdas.append(10000)
		interferences = [""]
		hels = [""]

	#~ massBins = [1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400]
	massBins = [1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400]
	# ~ massBins = [1800]
	signalYields = {}
	
	#~ names = ["default","resolution","scale","ID"]
	
	for label in labels:
			print (label)
		# ~ if "dimuon" in label:
			# ~ suffixes = suffixesMu
		# ~ else:
			# ~ suffixes = suffixesEle
		# ~ for suffix in suffixes:
			if "dimuon" in label and not suffix in suffixesMu:
				continue
			if "dielectron" in label and not suffix in suffixesEle:
				continue
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
							if useADD:
								if "2016" in label:
									# ~ fitFile = TFile("%s_%s_%s_%s_parametrization_fixinf_limitp0_limitp1_limitp2_2016.root"%(name,suffix,histo.lower(),cs),"READ")
									fitFile = TFile(outf+"%s_%s_%s_%s_scanmass_fixinf_2016.root"%(name,suffix,histo.lower(),cs),"READ")
								elif "2018" in label:
									fitFile = TFile(outf+"%s_%s_%s_%s_scanmass_fixinf_2018.root"%(name,suffix,histo.lower(),cs),"READ")
								else:
									fitFile = TFile(outf+"%s_%s_%s_%s_scanmass_fixinf.root"%(name,suffix,histo.lower(),cs),"READ")
							else:																	
								if "2016" in label:	
									fitFile = TFile(outf+"%s_%s_%s_%s_scanmass_fixinf_limitp0_limitp1_limitp2_2016.root"%(name,suffix,histo.lower(),cs),"READ")
								elif "2018" in label:	
									fitFile = TFile(outf+"%s_%s_%s_%s_scanmass_fixinf_limitp0_limitp1_limitp2_2018.root"%(name,suffix,histo.lower(),cs),"READ")
								else:	
									fitFile = TFile(outf+"%s_%s_%s_%s_scanmass_fixinf_limitp0_limitp1_limitp2.root"%(name,suffix,histo.lower(),cs),"READ")
							# ~ print (fitFile.ls())		
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
									functionUnc = fitFile.Get("fn_unc_m%d_%s"%(massBin,model))
									uncert = ((functionUnc.Eval(l)/function.Eval(l))**2 + (functionUnc.Eval(100000)/function.Eval(100000)))**0.5
									signalYields["%s_%s_%s"%(name,label,histo)][str(massBin)] = [(function.Eval(l)-function.Eval(100000)),uncert]
						if useADD: continue
			
						if "dimuon" in label:
							name = "CITo2Mu_Lam%dTeV%s"%(l,hel)
							nameCon = "CITo2Mu_Lam%dTeV%s"%(l,"Con"+hel)
							nameDes = "CITo2Mu_Lam%dTeV%s"%(l,"Des"+hel)
						else:	
							name = "CITo2E_Lam%dTeV%s"%(l,hel)
							nameCon = "CITo2E_Lam%dTeV%s"%(l,"Con"+hel)
							nameDes = "CITo2E_Lam%dTeV%s"%(l,"Des"+hel)
						signalYields["%s_%s_%s"%(name,label,histo)] = {}	
						for index, massBin in enumerate(massBins):
							signalYields["%s_%s_%s"%(name,label,histo)][str(massBin)] = [(signalYields["%s_%s_%s"%(nameCon,label,histo)][str(massBin)][0] + signalYields["%s_%s_%s"%(nameDes,label,histo)][str(massBin)][0])/2,((signalYields["%s_%s_%s"%(nameCon,label,histo)][str(massBin)][1]**2 + signalYields["%s_%s_%s"%(nameDes,label,histo)][str(massBin)][1]**2))**0.5]
					
			


			addci = "CI"
			if useADD: addci = "ADD"
			if "dimuon" in label:
				fileName = "%ssignalYieldsSingleBin"%addci
			else:
				fileName = "%ssignalYieldsSingleBinEle"%addci
			
			if args.truncation:
				fileName+="Truncation"

			
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
