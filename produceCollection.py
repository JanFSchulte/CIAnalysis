from defs import getPlot, Backgrounds, Backgrounds2016, Backgrounds2018, Signals, Signals2016, Signals2016ADD, Data, Data2016, Data2016el, Data2018, path, plotList, zScale, zScale2016, zScale2018
import ROOT
from helpers import *




#savePath="/depot/cms/users/minxi/lepratio/LimitSetting/data/"
#f=ROOT.TFile.Open(savePath+"NewCollection.root","RECREATE")

def getTag(flavor, year, uncertainty,category):

	print(flavor+"_"+str(year)+"_"+uncertainty+"_"+category)
	flavors={"mu":1000,"el":2000}
	years={"2016":100,"2017":200,"2018":300}
	categories={"bb":10,"be":20}
	uncertainties={"scaleUp":1,"scaleDown":1,"reso":2,"ID":3,"PUUp":4,"PUDown":4,"prefireUp":5,"prefireDown":5}
	if flavor=="mu":
		if uncertainty=="scaleUp" or uncertainty=="scaleDown":
			category="bb"
			if year==2018: year=2017
		elif uncertainty=="ID" or uncertainty=="reso":
			if category=="bb":
				if year==2018: 
					year=2017

	elif flavor=="el":
		if uncertainty.find("PU")!=-1 or uncertainty.find("prefire")!=-1:
			year=2016
			category="bb"
	
	tag=flavors[flavor]+years[str(year)]+uncertainties[uncertainty]+categories[category]		 
	print(str(tag))		
	return str(tag)





if __name__ == "__main__":
	
	savePath="/depot/cms/users/minxi/lepratio/LimitSetting/data/"
	#savePath="/depot/cms/users/minxi/study_novfit/mass_spectrum/"
	bins=binning()	
	lumi_el = [35.9*1000,41.529*1000,59.97*1000]
	lumi_mu = [36.3*1000,42.135*1000,61.608*1000]
	muplots_bb={"default":"massPlotBB","scaleUp":"massPlotBBScaleUpNoLog","scaleDown":"massPlotBBScaleDownNoLog","reso":"massPlotBBSmearNoLog","ID":"massPlotBBMuonIDNoLog"}
	muplots_be={"default":"massPlotBE","scaleUp":"massPlotBEScaleUpNoLog","scaleDown":"massPlotBEScaleDownNoLog","reso":"massPlotBESmearNoLog","ID":"massPlotBEMuonIDNoLog"}
	eleplots_bb={"default":"massPlotEleBB","scaleUp":"massPlotEleBBScaleUpNoLog","scaleDown":"massPlotEleBBScaleDownNoLog","PUUp":"massPlotEleBBPUScaleUpNoLog","PUDown":"massPlotEleBBPUScaleDownNoLog","prefireUp":"massPlotEleBBPrefireUpNoLog","prefireDown":"massPlotEleBBPrefireDownNoLog"}
 	eleplots_be={"default":"massPlotEleBE","scaleUp":"massPlotEleBEScaleUpNoLog","scaleDown":"massPlotEleBEScaleDownNoLog","PUUp":"massPlotEleBEPUScaleUpNoLog","PUDown":"massPlotEleBEPUScaleDownNoLog","prefireUp":"massPlotEleBEPrefireUpNoLog","prefireDown":"massPlotEleBEPrefireDownNoLog"}
	plot_mu_bb={}
	plot_el_bb={}
	plot_mu_be={}
	plot_el_be={}
	plotGenbb_mu = getPlot("massPlotBBGenNoLog")
	plotGenbe_mu = getPlot("massPlotBEGenNoLog")
	plotGenbb_el = getPlot("massPlotEleBBGenNoLog")
        plotGenbe_el = getPlot("massPlotEleBEGenNoLog")
	plotGenee_el = getPlot("massPlotEleEEGenNoLog")
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
	eventCounts_mu = totalNumberOfGeneratedEvents(path,plot_mu_bb["default"].muon)
	eventCounts_el = totalNumberOfGeneratedEvents(path,plot_el_bb["default"].muon)
	negWeights_mu = negWeightFractions(path,plot_mu_bb["default"].muon)
	negWeights_el = negWeightFractions(path,plot_el_bb["default"].muon)
	f=ROOT.TFile.Open(savePath+"NewCollection.root","RECREATE")
	#f=ROOT.TFile.Open(savePath+"vfitCollection.root","RECREATE")
	#f_jet=ROOT.TFile.Open(".root","RECREATE")
	yieldFile=open(savePath+"processYield.txt",'w')
	#for year in [2016,2017,2018]:
	for year in [2016,2017,2018]:
		print year
		if year==2016:
			processes_other_el=Process(getattr(Backgrounds2016,"OtherEle"),eventCounts_el,negWeights_el)
        		processes_other_mu=Process(getattr(Backgrounds2016,"Other"),eventCounts_mu,negWeights_mu)
			processes_dy_el=Process(getattr(Backgrounds2016,"DrellYan"),eventCounts_el,negWeights_el)
        		processes_dy_mu=Process(getattr(Backgrounds2016,"DrellYan"),eventCounts_mu,negWeights_mu)
			data = Process(Data2016, normalized=True)
			data_el = Process(Data2016el, normalized=True)
			i=0
			Znorm=zScale2016
		elif year==2018:
			processes_other_el=Process(getattr(Backgrounds2018,"Other"),eventCounts_el,negWeights_el)
                        processes_other_mu=Process(getattr(Backgrounds2018,"Other"),eventCounts_mu,negWeights_mu)
                        processes_dy_el=Process(getattr(Backgrounds2018,"DrellYan"),eventCounts_el,negWeights_el)
                        processes_dy_mu=Process(getattr(Backgrounds2018,"DrellYan"),eventCounts_mu,negWeights_mu)
                        data = Process(Data2018, normalized=True)
                        i=2
                        Znorm=zScale2018
		else:
			processes_other_el=Process(getattr(Backgrounds,"Other"),eventCounts_el,negWeights_el)
                        processes_other_mu=Process(getattr(Backgrounds,"Other"),eventCounts_mu,negWeights_mu)
                        processes_dy_el=Process(getattr(Backgrounds,"DrellYan"),eventCounts_el,negWeights_el)
                        processes_dy_mu=Process(getattr(Backgrounds,"DrellYan"),eventCounts_mu,negWeights_mu)
                        data = Process(Data, normalized=True)
                        i=1
                        Znorm=zScale
	
	
			
		hist=data.loadHistogram(plot_mu_bb["default"],lumi_mu[i],Znorm["muons"])
		hist.SetName("muon_%s_bb_data_obs"%str(year))
		Ls=[]
		for j in range(hist.GetNbinsX()+1):
			Ls.append(hist.GetBinContent(j))
		print(Ls)
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

		f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral()
		yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		hist=data.loadHistogram(plot_mu_be["default"],lumi_mu[i],Znorm["muons"])
        	hist.SetName("muon_%s_be_data_obs"%str(year))
		Ls=[]
		for j in range(hist.GetNbinsX()+1):
                        Ls.append(hist.GetBinContent(j))
                print(Ls)

		#for j in range(hist.GetNbinsX()+1):
                #	print(hist.GetBinContent(j))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

        	f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral()
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		if year==2016:
			hist=data_el.loadHistogram(plot_el_bb["default"],lumi_mu[i],Znorm["muons"])
		else:
			hist=data.loadHistogram(plot_el_bb["default"],lumi_mu[i],Znorm["muons"])
        	hist.SetName("electron_%s_bb_data_obs"%str(year))
		Ls=[]
		for j in range(hist.GetNbinsX()+1):
                        Ls.append(hist.GetBinContent(j))
                print(Ls)

		#for j in range(hist.GetNbinsX()+1):
                #        print(hist.GetBinContent(j))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

        	f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral()
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		if year==2016:
			hist=data_el.loadHistogram(plot_el_be["default"],lumi_mu[i],Znorm["muons"])
		else:
			hist=data.loadHistogram(plot_el_be["default"],lumi_mu[i],Znorm["muons"])
        	hist.SetName("electron_%s_be_data_obs"%str(year))
		Ls=[]
		for j in range(hist.GetNbinsX()+1):
                        Ls.append(hist.GetBinContent(j))
                print(Ls)

		#for j in range(hist.GetNbinsX()+1):
                #        print(hist.GetBinContent(j))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

        	f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral()
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		hist=processes_dy_mu.loadHistogram(plotGenbb_mu,lumi_mu[i],Znorm["muons"])
		hist.SetName("muon_%s_bb_genMass_DrellYan"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

		f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		hist=processes_dy_mu.loadHistogram(plotGenbe_mu,lumi_mu[i],Znorm["muons"])
        	hist.SetName("muon_%s_be_genMass_DrellYan"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

		f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		hist=processes_other_mu.loadHistogram(plotGenbb_mu,lumi_mu[i],Znorm["muons"])
        	hist.SetName("muon_%s_bb_genMass_Other"%str(year)) 
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

        	f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
        	hist=processes_other_mu.loadHistogram(plotGenbe_mu,lumi_mu[i],Znorm["muons"])
        	hist.SetName("muon_%s_be_genMass_Other"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

        	f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		hist=processes_dy_el.loadHistogram(plotGenbb_el,lumi_el[i],Znorm["electrons"][0])
        	hist.SetName("electron_%s_bb_genMass_DrellYan"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

		f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
        	hist=processes_dy_el.loadHistogram(plotGenbe_el,lumi_el[i],Znorm["electrons"][1])
        	hist.SetName("electron_%s_be_genMass_DrellYan"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)
	 
        	f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
        	hist=processes_other_el.loadHistogram(plotGenbb_el,lumi_el[i],Znorm["electrons"][0])
        	hist.SetName("electron_%s_bb_genMass_Other"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

		f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")

        	hist=processes_other_el.loadHistogram(plotGenbe_el,lumi_el[i],Znorm["electrons"][1])
        	hist.SetName("electron_%s_be_genMass_Other"%str(year))
                for j in range(hist.GetNbinsX()+1):
                        if hist.GetBinContent(j)<0:
                                hist.SetBinContent(j,0)

 		f.WriteObject(hist,hist.GetName())
		Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")

		hist=processes_other_el.loadHistogram(plotGenee_el,lumi_el[i],Znorm["electrons"][2])
        	hist.SetName("electron_%s_ee_genMass_Other"%str(year))
                for j in range(hist.GetNbinsX()+1):
                	if hist.GetBinContent(j)<0:
                        	hist.SetBinContent(j,0)
        	f.WriteObject(hist,hist.GetName())
        	Yield=hist.Integral(bins[0],bins[-1])
        	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
		norm_el_bb=processes_dy_el.loadHistogram(plot_el_bb["default"],lumi_el[i],Znorm["electrons"][0]).Integral()
		norm_el_be=processes_dy_el.loadHistogram(plot_el_be["default"],lumi_el[i],Znorm["electrons"][1]).Integral()
		norm_mu_bb=processes_dy_mu.loadHistogram(plot_mu_bb["default"],lumi_mu[i],Znorm["muons"]).Integral()
		norm_mu_be=processes_dy_mu.loadHistogram(plot_mu_be["default"],lumi_mu[i],Znorm["muons"]).Integral()
		print((norm_el_bb,norm_el_be,norm_mu_bb,norm_mu_be))
		for key in eleplots_bb.keys():
                	hist=processes_dy_el.loadHistogram(plot_el_bb[key],lumi_el[i],Znorm["electrons"][0])
                	if key == "default":
                        	hist.SetName("electron_%s_bb_recoMass_DrellYan1"%str(year))
			#norm_el_bb=hist.Integral()
                	else:
				tag = getTag("el", year, key, "bb")
                        	hist.SetName(("electron_%s_bb_recoMass_DrellYan1_"%str(year))+tag+key)
		
			Yield=hist.Integral()
                	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
                        for j in range(hist.GetNbinsX()+1):
                                if hist.GetBinContent(j)<0:
                                        hist.SetBinContent(j,0)
                	f.WriteObject(hist,hist.GetName())
                	hist=processes_dy_el.loadHistogram(plot_el_be[key],lumi_el[i],Znorm["electrons"][1])
                	if key == "default":
                        	hist.SetName("electron_%s_be_recoMass_DrellYan1"%str(year))
			#norm_el_be=hist.Integral()
                	else:
				tag = getTag("el", year, key, "be")
                        	hist.SetName(("electron_%s_be_recoMass_DrellYan1_"%str(year))+tag+key)

			Yield=hist.Integral()
        		yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
                        for j in range(hist.GetNbinsX()+1):
                                if hist.GetBinContent(j)<0:
                                        hist.SetBinContent(j,0)
                	f.WriteObject(hist,hist.GetName())
                	hist=processes_other_el.loadHistogram(plot_el_bb[key],lumi_el[i],Znorm["electrons"][0])
                	if key == "default":
                        	hist.SetName("electron_%s_bb_recoMass_Other"%str(year))
				ee_bb_list=[]
				print("ee bb mc")
                		for j in range(hist.GetNbinsX()+1):
					ee_bb_list.append(hist.GetBinContent(j))
				print(ee_bb_list)
                	else:
				tag = getTag("el", year, key, "bb")
                        	hist.SetName(("electron_%s_bb_recoMass_Other_"%str(year))+tag+key)
			Yield=hist.Integral()
        		yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
                        for j in range(hist.GetNbinsX()+1):
                                if hist.GetBinContent(j)<0:
                                        hist.SetBinContent(j,0)
                	f.WriteObject(hist,hist.GetName())
                	hist=processes_other_el.loadHistogram(plot_el_be[key],lumi_el[i],Znorm["electrons"][1])
                	if key == "default":
                        	hist.SetName("electron_%s_be_recoMass_Other"%str(year))
				ee_be_list=[]
				print("ee be mc")
                                for j in range(hist.GetNbinsX()+1):
					ee_be_list.append(hist.GetBinContent(j))
                                print(ee_be_list)
   
                	else:
				tag = getTag("el", year, key, "be")
                        	hist.SetName(("electron_%s_be_recoMass_Other_"%str(year))+tag+key)
			Yield=hist.Integral()
                	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
                        for j in range(hist.GetNbinsX()+1):
                        	if hist.GetBinContent(j)<0:
                                	hist.SetBinContent(j,0)
                	f.WriteObject(hist,hist.GetName())

	
		for key in muplots_bb.keys():
			hist=processes_dy_mu.loadHistogram(plot_mu_bb[key],lumi_mu[i],Znorm["muons"])
			if key == "default":
				hist.SetName("muon_%s_bb_recoMass_DrellYan"%str(year))
			#norm_mu_bb=hist.Integral()
			else:
				tag = getTag("mu", year, key, "bb")
				hist.SetName(("muon_%s_bb_recoMass_DrellYan_"%str(year))+tag+key)
			Yield=hist.Integral()
                	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
			if key=="ID" or key=="reso":
				tag = getTag("mu", year, key, "bb")
				hist.SetName(("muon_%s_bb_recoMass_DrellYan_"%str(year))+tag+key+"Up")
				default=processes_dy_mu.loadHistogram(plot_mu_bb["default"],lumi_mu[i],Znorm["muons"])
				default.Scale(2)
				default.Add(hist,-1)
				default.SetName(("muon_%s_bb_recoMass_DrellYan_"%str(year))+tag+key+"Down")
				hist.Scale(norm_el_bb/norm_mu_bb)
				default.Scale(norm_el_bb/norm_mu_bb)
                                for j in range(hist.GetNbinsX()+1):
                                        if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)
                                        if default.GetBinContent(j)<0:
                                                default.SetBinContent(j,0)

				f.WriteObject(hist,hist.GetName())
				f.WriteObject(default,default.GetName())
			else:
				print(hist.Integral())
				hist.Scale(norm_el_bb/norm_mu_bb)
				print(hist.Integral()) 
                                for j in range(hist.GetNbinsX()+1):
                                        if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)
				f.WriteObject(hist,hist.GetName())

			hist=processes_dy_mu.loadHistogram(plot_mu_be[key],lumi_mu[i],Znorm["muons"])
			if key == "default":
                        	hist.SetName("muon_%s_be_recoMass_DrellYan"%str(year))
			#norm_mu_be=hist.Integral()
                	else:
				tag = getTag("mu", year, key, "be")
                        	hist.SetName(("muon_%s_be_recoMass_DrellYan_"%str(year))+tag+key)
			Yield=hist.Integral()
                	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
			if key=="ID" or key=="reso":
				tag = getTag("mu", year, key, "be")
                        	hist.SetName(("muon_%s_be_recoMass_DrellYan_"%str(year))+tag+key+"Up")
                        	default=processes_dy_mu.loadHistogram(plot_mu_be["default"],lumi_mu[i],Znorm["muons"])
                        	default.Scale(2)
                        	default.Add(hist,-1)
                        	default.SetName(("muon_%s_be_recoMass_DrellYan_"%str(year))+tag+key+"Down")
                        	hist.Scale(norm_el_be/norm_mu_be)
                        	default.Scale(norm_el_be/norm_mu_be)
                                for j in range(hist.GetNbinsX()+1):
                                        if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)
                                        if default.GetBinContent(j)<0:
                                                default.SetBinContent(j,0)
                        	f.WriteObject(hist,hist.GetName())
                        	f.WriteObject(default,default.GetName())
                	else:
                        	hist.Scale(norm_el_be/norm_mu_be) 
				for j in range(hist.GetNbinsX()+1):
                                        if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)
                        	f.WriteObject(hist,hist.GetName())
                
			hist=processes_other_mu.loadHistogram(plot_mu_bb[key],lumi_mu[i],Znorm["muons"])

			if key == "default":
                        	hist.SetName("muon_%s_bb_recoMass_Other"%str(year))
				print("mu bb mc")
				mu_bb_list=[]
                                for j in range(hist.GetNbinsX()+1):
					mu_bb_list.append(hist.GetBinContent(j))
                                print(mu_bb_list)
                	else:
				tag = getTag("mu", year, key, "bb")
                        	hist.SetName(("muon_%s_bb_recoMass_Other_"%str(year))+tag+key)
			Yield=hist.Integral()
                	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n")
			if key=="ID" or key=="reso":
				tag = getTag("mu", year, key, "bb")
                        	hist.SetName(("muon_%s_bb_recoMass_Other_"%str(year))+tag+key+"Up")
                        	default=processes_other_mu.loadHistogram(plot_mu_bb["default"],lumi_mu[i],Znorm["muons"])
                        	default.Scale(2)
                        	default.Add(hist,-1)
                        	default.SetName(("muon_%s_bb_recoMass_Other_"%str(year))+tag+key+"Down")
				for j in range(hist.GetNbinsX()+1):
                                        if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)
                                        if default.GetBinContent(j)<0:
                                                default.SetBinContent(j,0)

                        	f.WriteObject(hist,hist.GetName())
                        	f.WriteObject(default,default.GetName())
                	else:
				for j in range(hist.GetNbinsX()+1):
                                	if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)
                        	f.WriteObject(hist,hist.GetName()) 
              
                	hist=processes_other_mu.loadHistogram(plot_mu_be[key],lumi_mu[i],Znorm["muons"])
                	if key == "default":
                        	hist.SetName("muon_%s_be_recoMass_Other"%str(year))
				print("mu be mc")
				mu_be_list=[]
                                for j in range(hist.GetNbinsX()+1):
					mu_be_list.append(hist.GetBinContent(j))
                                print(mu_be_list)
                	else:
				tag = getTag("mu", year, key, "be")
                        	hist.SetName(("muon_%s_be_recoMass_Other_"%str(year))+tag+key)
			Yield=hist.Integral()
                	yieldFile.write(hist.GetName()+": "+str(Yield)+"\n") 
			if key=="ID" or key=="reso":
				tag = getTag("mu", year, key, "be")
                        	hist.SetName(("muon_%s_be_recoMass_Other_"%str(year))+tag+key+"Up")
                        	default=processes_other_mu.loadHistogram(plot_mu_be["default"],lumi_mu[i],Znorm["muons"])
                        	default.Scale(2)
                        	default.Add(hist,-1)
                        	default.SetName(("muon_%s_be_recoMass_Other_"%str(year))+tag+key+"Down")
				for j in range(hist.GetNbinsX()+1):
					if hist.GetBinContent(j)<0:
						hist.SetBinContent(j,0)
					if default.GetBinContent(j)<0:
						default.SetBinContent(j,0)
                        	f.WriteObject(hist,hist.GetName())
                        	f.WriteObject(default,default.GetName())
                	else:
				for j in range(hist.GetNbinsX()+1):
                                        if hist.GetBinContent(j)<0:
                                                hist.SetBinContent(j,0)   
                        	f.WriteObject(hist,hist.GetName())
	f.Save()
	#print(processes2016_other_el_bb)
	
	
