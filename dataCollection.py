import ROOT
from defs import Backgrounds, Backgrounds2016, Backgrounds2018, zScale, zScale2016, zScale2018, fileNamesEle, fileNames, crossSections  
import math
import os
from copy import copy
import numpy 
import root_numpy
from copy import deepcopy
from histHandle import histHandle


def reAverage(hist):
        for i in range(hist.GetNbinsX()+2):
                val=hist.GetBinContent(i)
                err=hist.GetBinError(i)
                wid=hist.GetBinWidth(i)
		print wid
		print val
                val=val*wid
                err=err*wid
                hist.SetBinContent(i,val)
                hist.SetBinError(i,err)

path="/depot/cms/users/minxi/"
class processHandle:
	def __init__(self, name):
		self.name=name
		self.filename=""
		self.crossSection=0
		self.year="2017"
		self.flavor=""


processList=[]
bng=[200,300,400,500,690,900,1250,1610, 2000, 3000, 4000]
bng=numpy.asarray(bng,dtype=numpy.float64)
tempBng=[200,300,400,500,690,900,1250,1610, 2000,2500 ,3000,3500, 4000]
tempBng=numpy.asarray(tempBng,dtype=numpy.float64)

for filename in os.listdir(path):
	flag=False
	for key, value in fileNamesEle.items():
		if value==filename and filename.find("Double")==-1 and filename.find("Single")==-1:
			process=processHandle(key)
			process.filename=filename
			print filename
			process.flavor="electron"
			process.crossSection=crossSections[key]
			if filename.find("2018")!=-1:
				process.year="2018"
			elif filename.find("2016")!=-1:
				process.year="2016"
			flag=True
	for key, value in fileNames.items(): 
		if value==filename and filename.find("Double")==-1 and filename.find("Single")==-1:
                        process=processHandle(key)
                        process.filename=filename
			print filename
                        process.flavor="muon"
                        process.crossSection=crossSections[key]
			if filename.find("2018")!=-1:
                                process.year="2018"
                        elif filename.find("2016")!=-1:
                                process.year="2016"
			flag=True
	if flag:
		processList.append(process)


#temphist1d=ROOT.TH1D("temp1d","temphist",len(tempBng)-1,tempBng)
#temphist2d=ROOT.TH2D("temp2d","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
histContainer=histHandle("dataCollection",["flavor","year","category","type","source","uncertainty","rebin"])
struct={"flavor":["muon","electron"],"year":["2016"],"category":["bb","be"],"type":["genMass","recoMass","response"],"source":["DrellYan"],"uncertainty":["resolution","scaleUp","scaleDown","pileUp","pileDown","ID","prefireUp","prefireDown","kFac","PDFUp","PDFDown"]}


directoryEle="ElectronSelectionElectronsAllSignsHistos/"
directoryMu2016="Our2016MuonsPlusMuonsMinusHistos/"
directoryMu2017="Our2017MuonsPlusMuonsMinusHistos/"
directoryMu2018="Our2018MuonsPlusMuonsMinusHistos/"
histosListEle=["DielectronMass_bb","DielectronMass_be","DielectronMassPDFUp_bb","DielectronMassPDFUp_be","DielectronMassPDFDown_bb","DielectronMassPDFDown_be","DielectronMass_gen_bb","DielectronMass_gen_be","DielectronMassPDFUp_gen_bb","DielectronMassPDFUp_gen_be","DielectronMassPDFDown_gen_bb","DielectronMassPDFDown_gen_be","DielectronMass_bb_kFac","DielectronMass_be_kFac","DielectronMass_gen_bb_kFac","DielectronMass_gen_be_kFac","DielectronMassPrefireUp_gen_bb","DielectronMassPrefireUp_gen_be","DielectronMassPrefireDown_gen_bb","DielectronMassPrefireDown_gen_be","DielectronMassPUScaleUp_gen_bb","DielectronMassPUScaleUp_gen_be","DielectronMassPUScaleDown_gen_bb","DielectronMassPUScaleDown_gen_be","DielectronResponse_bb","DielectronResponse_be","DielectronResponse_bb_kFac","DielectronResponse_be_kFac","DielectronResponseMassScaleUp_bb","DielectronResponseMassScaleUp_be","DielectronResponseMassScaleDown_bb","DielectronResponseMassScaleDown_be","DielectronResponsePDFUp_be","DielectronResponsePDFDown_be","DielectronResponsePDFUp_bb","DielectronResponsePDFDown_bb"]
histosListMu=["DimuonMassVertexConstrained_bb","DimuonMassVertexConstrained_be","DimuonMassVertexConstrainedPDFUp_bb","DimuonMassVertexConstrainedPDFUp_be","DimuonMassVertexConstrainedPDFDown_bb","DimuonMassVertexConstrainedPDFDown_be","DimuonMassVertexConstrained_gen_bb","DimuonMassVertexConstrained_gen_be","DimuonMassVertexConstrainedPDFUp_gen_bb","DimuonMassVertexConstrainedPDFUp_gen_be","DimuonMassVertexConstrainedPDFDown_gen_bb","DimuonMassVertexConstrainedPDFDown_gen_be","DimuonMassVertexConstrained_bb_kFac","DimuonMassVertexConstrained_be_kFac","DimuonMassVertexConstrained_gen_bb_kFac","DimuonMassVertexConstrained_gen_be_kFac","DimuonMassVertexConstrainedMuonID_gen_bb","DimuonMassVertexConstrainedMuonID_gen_be","DimuonResponse_bb","DimuonResponse_be","DimuonResponsePDFUp_bb","DimuonResponsePDFUp_be","DimuonResponsePDFDown_bb","DimuonResponsePDFDown_be","DimuonResponseMassScaleUp_bb","DimuonResponseMassScaleUp_be","DimuonResponseMassScaleDown_bb","DimuonResponseMassScaleDown_be","DimuonResponse_bb_kFac","DimuonResponse_be_kFac","DimuonResponseSmear_bb","DimuonResponseSmear_be"]
for process in processList:
	filename=path+process.filename
	#print filename
	f=ROOT.TFile.Open(filename)
	nEvents=f.FindObjectAny("Events").GetBinContent(1)
	neg=1-2*f.FindObjectAny("weights").GetBinContent(1)/(f.FindObjectAny("weights").GetBinContent(1)+f.FindObjectAny("weights").GetBinContent(2))
	nEvents=nEvents/neg
	xSec=process.crossSection
	if process.name in Backgrounds2018.DrellYan.subprocesses:
		if process.flavor=="electron":
			for histname in histosListEle:
				treePath={}
				treePath["source"]="DrellYan"
        			treePath["year"]=process.year
        			treePath["flavor"]=process.flavor
				if histname.find("_bb")!=-1:
					treePath["category"]="bb"
				else:
					treePath["category"]="be"
				if histname.find("Response")!=-1:
					treePath["type"]="response"
				elif histname.find("_gen")!=-1:
					treePath["type"]="genMass"
				else:
					treePath["type"]="recoMass"
				if histname.find("MassScaleUp")!=-1:
					treePath["uncertainty"]="scaleUp"
				elif histname.find("MassScaleDown")!=-1:
					treePath["uncertainty"]="scaleDown"
				elif histname.find("_kFac")!=-1:
					treePath["uncertainty"]="kFac"
				elif histname.find("PrefireUp")!=-1:        
                                        treePath["uncertainty"]="prefireUp"
				elif histname.find("PrefireDown")!=-1:
                                        treePath["uncertainty"]="prefireDown"
				elif histname.find("PUScaleUp")!=-1:
                                        treePath["uncertainty"]="pileUp"
                                elif histname.find("PUScaleDown")!=-1:
                                        treePath["uncertainty"]="pileDown"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
				hist=f.Get(directoryEle+histname)
				hist.Scale(xSec/nEvents)
				histContainer.addHist(hist,treePath)			
			#	print "f-1"
		elif process.flavor=="muon":
			for histname in histosListMu:
				treePath={}
				treePath["source"]="DrellYan"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor

                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("Smear")!=-1:
                                        treePath["uncertainty"]="resolution"
                                elif histname.find("MuonID")!=-1:
                                        treePath["uncertainty"]="ID"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"        
                                hist=f.Get(directoryMu2018+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f0"
	if process.name in Backgrounds2018.Other.subprocesses:
		if process.flavor=="electron":
			for histname in histosListEle:
				treePath={}
				treePath["source"]="Other"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
                        	if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("PrefireUp")!=-1:
                                        treePath["uncertainty"]="prefireUp"
                                elif histname.find("PrefireDown")!=-1:
                                        treePath["uncertainty"]="prefireDown"
                                elif histname.find("PUScaleUp")!=-1:
                                        treePath["uncertainty"]="pileUp"
                                elif histname.find("PUScaleDown")!=-1:
                                        treePath["uncertainty"]="pileDown"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryEle+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f1"
                elif process.flavor=="muon":
			for histname in histosListMu:
				treePath={}
				treePath["source"]="Other"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("Smear")!=-1:
                                        treePath["uncertainty"]="resolution"
                                elif histname.find("MuonID")!=-1:
                                        treePath["uncertainty"]="ID"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
				#print histname
                                hist=f.Get(directoryMu2018+histname)
                                hist.Scale(xSec/nEvents)
				#print histname
				#print treePath
				#print filename
                                histContainer.addHist(hist,treePath)
			#	print "f2"

	if process.name in Backgrounds.DrellYan.subprocesses:
                if process.flavor=="electron":
			for histname in histosListEle:
				treePath={}
				treePath["source"]="DrellYan"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("PrefireUp")!=-1:
                                        treePath["uncertainty"]="prefireUp"
                                elif histname.find("PrefireDown")!=-1:
                                        treePath["uncertainty"]="prefireDown"
                                elif histname.find("PUScaleUp")!=-1:
                                        treePath["uncertainty"]="pileUp"
                                elif histname.find("PUScaleDown")!=-1:
                                        treePath["uncertainty"]="pileDown"
                                hist=f.Get(directoryEle+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f3"
                elif process.flavor=="muon":
			for histname in histosListMu:
				treePath={}
				treePath["source"]="DrellYan"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("Smear")!=-1:
                                        treePath["uncertainty"]="resolution"
                                elif histname.find("MuonID")!=-1:
                                        treePath["uncertainty"]="ID"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryMu2017+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f4"
        if process.name in Backgrounds.Other.subprocesses:
                if process.flavor=="electron":
			for histname in histosListEle:
				treePath={}
				treePath["source"]="Other"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("PrefireUp")!=-1:
                                        treePath["uncertainty"]="prefireUp"
                                elif histname.find("PrefireDown")!=-1:
                                        treePath["uncertainty"]="prefireDown"
                                elif histname.find("PUScaleUp")!=-1:
                                        treePath["uncertainty"]="pileUp"
                                elif histname.find("PUScaleDown")!=-1:
                                        treePath["uncertainty"]="pileDown"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryEle+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f5"
                elif process.flavor=="muon":
			for histname in histosListMu:
				treePath={}
				treePath["source"]="Other"
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("Smear")!=-1:
                                        treePath["uncertainty"]="resolution"
                                elif histname.find("MuonID")!=-1:
                                        treePath["uncertainty"]="ID"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryMu2017+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f6"
	if process.name in Backgrounds2016.DrellYan.subprocesses:
                if process.flavor=="electron":
			for histname in histosListEle:
				treePath={}
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
				treePath["source"]="DrellYan"
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("PrefireUp")!=-1:
                                        treePath["uncertainty"]="prefireUp"
                                elif histname.find("PrefireDown")!=-1:
                                        treePath["uncertainty"]="prefireDown"
                                elif histname.find("PUScaleUp")!=-1:
                                        treePath["uncertainty"]="pileUp"
                                elif histname.find("PUScaleDown")!=-1:
                                        treePath["uncertainty"]="pileDown"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
					#if histname.find("Response")!=-1:
					#	.replace("Mass","")
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
					#if directoryEle.find("Response")!=-1:
                                        #        directoryEle.replace("Mass","")
				print histname
                                hist=f.Get(directoryEle+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
			#	print "f7"
                elif process.flavor=="muon":
			for histname in histosListMu:
				treePath={}
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
				treePath["source"]="DrellYan"
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("Smear")!=-1:
                                        treePath["uncertainty"]="resolution"
                                elif histname.find("MuonID")!=-1:
                                        treePath["uncertainty"]="ID"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryMu2016+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
				#print "f8"
        if process.name in Backgrounds2016.Other.subprocesses:
                if process.flavor=="muon":
			for histname in histosListMu:
				treePath={}
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
				treePath["source"]="Other"
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("Smear")!=-1:
                                        treePath["uncertainty"]="resolution"
                                elif histname.find("MuonID")!=-1:
                                        treePath["uncertainty"]="ID"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryMu2016+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
				#print "f9"
	if process.name in Backgrounds2016.OtherEle.subprocesses:
		if process.flavor=="electron":
			for histname in histosListEle:
				treePath={}
                                treePath["year"]=process.year
                                treePath["flavor"]=process.flavor
				treePath["source"]="Other"
                                if histname.find("_bb")!=-1:
                                        treePath["category"]="bb"
                                else:
                                        treePath["category"]="be"
                                if histname.find("Response")!=-1:
                                        treePath["type"]="response"
                                elif histname.find("_gen")!=-1:
                                        treePath["type"]="genMass"
                                else:
                                        treePath["type"]="recoMass"
                                if histname.find("MassScaleUp")!=-1:
                                        treePath["uncertainty"]="scaleUp"
                                elif histname.find("MassScaleDown")!=-1:
                                        treePath["uncertainty"]="scaleDown"
                                elif histname.find("_kFac")!=-1:
                                        treePath["uncertainty"]="kFac"
                                elif histname.find("PrefireUp")!=-1:
                                        treePath["uncertainty"]="prefireUp"
                                elif histname.find("PrefireDown")!=-1:
                                        treePath["uncertainty"]="prefireDown"
                                elif histname.find("PUScaleUp")!=-1:
                                        treePath["uncertainty"]="pileUp"
                                elif histname.find("PUScaleDown")!=-1:
                                        treePath["uncertainty"]="pileDown"
				elif histname.find("PDFUp")!=-1:
                                        treePath["uncertainty"]="PDFUp"
                                elif histname.find("PDFDown")!=-1:
                                        treePath["uncertainty"]="PDFDown"
                                hist=f.Get(directoryEle+histname)
                                hist.Scale(xSec/nEvents)
                                histContainer.addHist(hist,treePath)
				print filename


print "produce combine"


lumi_e ={"2016":35.9*1000,"2017":41.529*1000,"2018":59.97*1000}
zScale_e={}
zScale_e["bb"] = {"2016":zScale2016["electrons"][1],"2017":zScale["electrons"][1],"2018":zScale2018["electrons"][1]}
zScale_e["be"] = {"2016":zScale2016["electrons"][2],"2017":zScale["electrons"][2],"2018":zScale2018["electrons"][2]}
lumi_mu = {"2016":36.3*1000,"2017":42.135*1000,"2018":61.608*1000}
zScale_mu = {"2016":zScale2016["muons"],"2017":zScale["muons"],"2018":zScale2018["muons"]}

for flavor in struct["flavor"]:
	for Type in struct["type"]:
		if Type=="response" and flavor=="muon":
                	for uncertainty in ["resolution","scaleUp","scaleDown","kFac","none"]:
				temphist1d=ROOT.TH1F("temp1d","temphist",20000,0,20000)
                		temphist2d=ROOT.TH2F("temp2d","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                		for category in struct["category"]:
					temphist1dex=ROOT.TH1F("temp1dex","temphist",20000,0,20000)
                                	temphist2dex=ROOT.TH2F("temp2dex","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexA=ROOT.TH1F("temp1dexA","temphist",20000,0,20000)
                                        temphist2dexA=ROOT.TH2F("temp2dexA","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexB=ROOT.TH1F("temp1dexB","temphist",20000,0,20000)
                                        temphist2dexB=ROOT.TH2F("temp2dexB","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                        		for year in struct["year"]:
						if uncertainty !="none":
                                        		path1={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
                                        		#path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"Other"}
                                        		path3={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"combine"}
						else:
							path1={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
                                                        #path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"Other"}
                                                        path3={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"combine"}
                                        	h1=histContainer.getHist(path1)
                                        	h2=histContainer.getHist(path2)
						h1.Scale(lumi_mu[year]*zScale_mu[year])
						h2.Scale(lumi_mu[year]*zScale_mu[year])
						#temphist2dexB.Add(h2)
						temphist2dexA.Add(h1)
                                        	h1.Add(h2)
						h3=h1.Clone()
						temphist2d.Add(h3)
						temphist2dex.Add(h3)
						#temphist2dexA.Add(h1)
						temphist2dexB.Add(h2)
                                        	histContainer.addHist(h1,path3)
					if uncertainty !="none":
						path4={"flavor":flavor,"category":"combine","uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
						#path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"Other"}
					else:
						path4={"flavor":flavor,"category":"combine","type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"Other"}
					histContainer.addHist(temphist2dex,path5)
					histContainer.addHist(temphist2dexA,path6)
					histContainer.addHist(temphist2dexB,path7)
				histContainer.addHist(temphist2d,path4)        
		elif Type=="response":
                	for uncertainty in ["scaleUp","scaleDown","kFac","none"]:
				temphist1d=ROOT.TH1F("temp1d","temphist",20000,0,20000)
                                temphist2d=ROOT.TH2F("temp2d","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
				for category in struct["category"]:
					temphist1dex=ROOT.TH1F("temp1dex","temphist",20000,0,20000)
                                        temphist2dex=ROOT.TH2F("temp2dex","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexA=ROOT.TH1F("temp1dexA","temphist",20000,0,20000)
                                        temphist2dexA=ROOT.TH2F("temp2dexA","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                        temphist1dexB=ROOT.TH1F("temp1dexB","temphist",20000,0,20000)
                                        temphist2dexB=ROOT.TH2F("temp2dexB","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                        for year in struct["year"]:
						if uncertainty!="none":
                                        		path1={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
                                                	#path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"Other"}
                                                	path3={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"combine"}
						else:
							path1={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
                                                        #path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"Other"}
							path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
                                                        path3={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"combine"}
                                                h1=histContainer.getHist(path1)
                                                h2=histContainer.getHist(path2)
						h1.Scale(lumi_e[year]*zScale_e[category][year])
						h2.Scale(lumi_e[year]*zScale_e[category][year])
						temphist2dexA.Add(h1)
                                                h1.Add(h2)
						h3=h1.Clone()
                                                temphist2d.Add(h3)
						temphist2dex.Add(h3)
						#temphist2dexA.Add(h1)
						temphist2dexB.Add(h2)
                                                histContainer.addHist(h1,path3)
					if uncertainty!="none":
                                		path4={"flavor":flavor,"category":"combine","uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"Other"}
					else:
						path4={"flavor":flavor,"category":"combine","type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"Other"}
					histContainer.addHist(temphist2dex,path5)
					histContainer.addHist(temphist2dexA,path6)
					histContainer.addHist(temphist2dexB,path7)
                                histContainer.addHist(temphist2d,path4)
		elif flavor=="electron" and Type=="genMass": 
                	for uncertainty in ["prefireUp","prefireDown","pileUp","pileDown","kFac","none"]:
				temphist1d=ROOT.TH1F("temp1d","temphist",20000,0,20000)
                                temphist2d=ROOT.TH2F("temp2d","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
				for category in struct["category"]:
					temphist1dex=ROOT.TH1F("temp1dex","temphist",20000,0,20000)
                                        temphist2dex=ROOT.TH2F("temp2dex","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexA=ROOT.TH1F("temp1dexA","temphist",20000,0,20000)
                                        temphist2dexA=ROOT.TH2F("temp2dexA","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexB=ROOT.TH1F("temp1dexB","temphist",20000,0,20000)
                                        temphist2dexB=ROOT.TH2F("temp2dexB","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                        for year in struct["year"]:
						if uncertainty !="none":
                                        		path1={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
                                                	#path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"Other"}
                                                	path3={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"combine"}
						else:
							path1={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
                                                        #path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"Other"}
                                                        path3={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"combine"}
                                                h1=histContainer.getHist(path1)
						h1.Scale(lumi_e[year]*zScale_e[category][year])
                                                h2=histContainer.getHist(path2)
						h2.Scale(lumi_e[year]*zScale_e[category][year])
						temphist1dexA.Add(h1)
                                                h1.Add(h2)
                                                histContainer.addHist(h1,path3)
						h3=h1.Clone()
                                                temphist1d.Add(h3)
						temphist1dex.Add(h3)
						temphist1dexB.Add(h2)
                                                #histContainer.addHist(h1,path3)
					if uncertainty !="none":
						path4={"flavor":flavor,"category":"combine","uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"Other"}
						path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
					else:
						path4={"flavor":flavor,"category":"combine","type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"Other"}
						path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
					histContainer.addHist(temphist1dex,path5)
					histContainer.addHist(temphist1dexA,path6)
					histContainer.addHist(temphist1dexB,path7)
                                histContainer.addHist(temphist1d,path4)
		elif flavor=="muon" and Type=="genMass":
                  	for uncertainty in ["ID","kFac","none"]:
				temphist1d=ROOT.TH1F("temp1d","temphist",20000,0,20000)
                                temphist2d=ROOT.TH2F("temp2d","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                for category in struct["category"]:
					temphist1dex=ROOT.TH1F("temp1dex","temphist",20000,0,20000)
                                	temphist2dex=ROOT.TH2F("temp2dex","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexA=ROOT.TH1F("temp1dexA","temphist",20000,0,20000)
                                        temphist2dexA=ROOT.TH2F("temp2dexA","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexB=ROOT.TH1F("temp1dexB","temphist",20000,0,20000)
                                        temphist2dexB=ROOT.TH2F("temp2dexB","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                        for year in struct["year"]:
						if uncertainty!="none":
							path1={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
							#path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"Other"}
							path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
							path3={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"combine"}
						else:
							path1={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
                                                        #path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"Other"}
                                                        path3={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"combine"}
						h1=histContainer.getHist(path1)
						h2=histContainer.getHist(path2)
						h1.Scale(lumi_mu[year]*zScale_mu[year])
						h2.Scale(lumi_mu[year]*zScale_mu[year])
						temphist1dexA.Add(h1)
						h1.Add(h2)
						histContainer.addHist(h1,path3)
						h3=h1.Clone()
                                                temphist1d.Add(h3)
						temphist1dex.Add(h3)
						temphist1dexB.Add(h2)
					if uncertainty!="none":
						path4={"flavor":flavor,"category":"combine","uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
                                                path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
						#path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"Other"}
					else:
						path4={"flavor":flavor,"category":"combine","type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"Other"}
					histContainer.addHist(temphist1dex,path5)
					histContainer.addHist(temphist1dexA,path6)
					histContainer.addHist(temphist1dexB,path7)
                                histContainer.addHist(temphist1d,path4)
		else:
                	for uncertainty in ["kFac","none"]:
				temphist1d=ROOT.TH1F("temp1d","temphist",20000,0,20000)
                                temphist2d=ROOT.TH2F("temp2d","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                for category in struct["category"]:
					temphist1dex=ROOT.TH1F("temp1dex","temphist",20000,0,20000)
                                	temphist2dex=ROOT.TH2F("temp2dex","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexA=ROOT.TH1F("temp1dexA","temphistA",20000,0,20000)
                                        temphist2dexA=ROOT.TH2F("temp2dexA","temphistA",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
					temphist1dexB=ROOT.TH1F("temp1dexB","temphist",20000,0,20000)
                                        temphist2dexB=ROOT.TH2F("temp2dexB","temphist",len(tempBng)-1,tempBng,len(tempBng)-1,tempBng)
                                        for year in struct["year"]:
						if uncertainty!="none":
                                                	path1={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"DrellYan"}
                                                	#path2={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"Other"}
                                                	path3={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":year,"source":"combine"}
						else:
							path1={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
							path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"DrellYan"}
                                                        #path2={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"Other"}
                                                        path3={"flavor":flavor,"category":category,"type":Type,"year":year,"source":"combine"}
                                                h1=histContainer.getHist(path1)
                                                h2=histContainer.getHist(path2)
						if flavor=="muon":
                                                        h1.Scale(lumi_mu[year]*zScale_mu[year])
							h2.Scale(lumi_mu[year]*zScale_mu[year])
                                                else:
                                                        h1.Scale(lumi_e[year]*zScale_e[category][year])
							h2.Scale(lumi_e[year]*zScale_e[category][year])
                                                #h1.Add(h2)
						#print path3
                                                #histContainer.addHist(h1,path3)
						h3=h1.Clone()	
                                                temphist1dexA.Add(h1)
						temphist1dexB.Add(h2)
						temphist1d.Add(h1)
						temphist1d.Add(h2)
                                                temphist1dex.Add(h1)
						temphist1dex.Add(h2)
						histContainer.addHist(h1,path3)
						histContainer.addHist(h2,path3)
					if uncertainty!="none":
                                		path4={"flavor":flavor,"category":"combine","uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"uncertainty":uncertainty,"type":Type,"year":"combine","source":"Other"}
					else:
						path4={"flavor":flavor,"category":"combine","type":Type,"year":"combine","source":"combine"}
						path5={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"combine"}
						path6={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
						path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"DrellYan"}
                                                #path7={"flavor":flavor,"category":category,"type":Type,"year":"combine","source":"Other"}
					histContainer.addHist(temphist1dex,path5)
					histContainer.addHist(temphist1dexA,path6)
					print path5
					histContainer.addHist(temphist1dexB,path7)
                                histContainer.addHist(temphist1d,path4)

					

print "add data and bkg"
fData=ROOT.TFile.Open("unfoldingData.root")
datamubb2018=fData.Get("datamubb2018")
histContainer.addHist(datamubb2018,{"flavor":"muon","category":"bb","year":"2018","type":"data"})
datamube2018=fData.Get("datamube2018")
histContainer.addHist(datamube2018,{"flavor":"muon","category":"be","year":"2018","type":"data"})
datamubb2017=fData.Get("datamubb2017")
histContainer.addHist(datamubb2017,{"flavor":"muon","category":"bb","year":"2017","type":"data"})
datamube2017=fData.Get("datamube2017")
histContainer.addHist(datamube2017,{"flavor":"muon","category":"be","year":"2017","type":"data"})
datamubb2016=fData.Get("datamubb2016")
histContainer.addHist(datamubb2016,{"flavor":"muon","category":"bb","year":"2016","type":"data"})
datamube2016=fData.Get("datamube2016")
histContainer.addHist(datamube2016,{"flavor":"muon","category":"be","year":"2016","type":"data"})
dataebb2018=fData.Get("dataebb2018")
histContainer.addHist(dataebb2018,{"flavor":"electron","category":"bb","year":"2018","type":"data"})
dataebe2018=fData.Get("dataebe2018")
histContainer.addHist(dataebe2018,{"flavor":"electron","category":"be","year":"2018","type":"data"})
dataebb2017=fData.Get("dataebb2017")
histContainer.addHist(dataebb2017,{"flavor":"electron","category":"bb","year":"2017","type":"data"})
dataebe2017=fData.Get("dataebe2017")
histContainer.addHist(dataebe2017,{"flavor":"electron","category":"be","year":"2017","type":"data"})
dataebb2016=fData.Get("dataebb2016")
histContainer.addHist(dataebb2016,{"flavor":"electron","category":"bb","year":"2016","type":"data"})
dataebe2016=fData.Get("dataebe2016")
histContainer.addHist(dataebe2016,{"flavor":"electron","category":"be","year":"2016","type":"data"})
datamu=datamubb2018.Clone()
datamu.Add(datamube2018)
datamu.Add(datamubb2017)
datamu.Add(datamube2017)
datamu.Add(datamubb2016)
datamu.Add(datamube2016)
histContainer.addHist(datamu,{"flavor":"muon","category":"combine","year":"combine","type":"data"})
datae=dataebb2018.Clone()
datae.Add(dataebe2018)
datae.Add(dataebb2017)
datae.Add(dataebe2017)
datae.Add(dataebb2016)
datae.Add(dataebe2016)
histContainer.addHist(datae,{"flavor":"electron","category":"combine","year":"combine","type":"data"})

datamubb=datamubb2018.Clone()
datamubb.Add(datamubb2017)
datamubb.Add(datamubb2016)
histContainer.addHist(datamubb,{"flavor":"muon","category":"bb","year":"combine","type":"data"})
dataebb=dataebb2018.Clone()
dataebb.Add(dataebb2017)
dataebb.Add(dataebb2016)
histContainer.addHist(dataebb,{"flavor":"electron","category":"bb","year":"combine","type":"data"})
datamube=datamube2018.Clone()
datamube.Add(datamube2017)
datamube.Add(datamube2016)
histContainer.addHist(datamube,{"flavor":"muon","category":"be","year":"combine","type":"data"})
dataebe=dataebe2018.Clone()
dataebe.Add(dataebe2017)
dataebe.Add(dataebe2016)
histContainer.addHist(dataebe,{"flavor":"electron","category":"be","year":"combine","type":"data"})


fMC=ROOT.TFile.Open("unfoldingMC_V2.root")
fMC1=ROOT.TFile.Open("unfoldingMC_V3.root")
bkgmubb2018=fMC.Get("bkgmubb2018")
histContainer.addHist(bkgmubb2018,{"flavor":"muon","category":"bb","year":"2018","type":"bkg"})
bkgmube2018=fMC.Get("bkgmube2018")
histContainer.addHist(bkgmube2018,{"flavor":"muon","category":"be","year":"2018","type":"bkg"})
bkgmubb2017=fMC.Get("bkgmubb2017")
histContainer.addHist(bkgmubb2017,{"flavor":"muon","category":"bb","year":"2017","type":"bkg"})
bkgmube2017=fMC.Get("bkgmube2017")
histContainer.addHist(bkgmube2017,{"flavor":"muon","category":"be","year":"2017","type":"bkg"})
bkgmubb2016=fMC.Get("bkgmubb2016")
histContainer.addHist(bkgmubb2016,{"flavor":"muon","category":"bb","year":"2016","type":"bkg"})
bkgmube2016=fMC.Get("bkgmube2016")
histContainer.addHist(bkgmube2016,{"flavor":"muon","category":"be","year":"2016","type":"bkg"})
bkgebb2018=fMC1.Get("bkgebb2018")
histContainer.addHist(bkgebb2018,{"flavor":"electron","category":"bb","year":"2018","type":"bkg"})
bkgebe2018=fMC1.Get("bkgebe2018")
histContainer.addHist(bkgebe2018,{"flavor":"electron","category":"be","year":"2018","type":"bkg"})
bkgebb2017=fMC1.Get("bkgebb2017")
histContainer.addHist(bkgebb2017,{"flavor":"electron","category":"bb","year":"2017","type":"bkg"})
bkgebe2017=fMC1.Get("bkgebe2017")
histContainer.addHist(bkgebe2018,{"flavor":"electron","category":"be","year":"2017","type":"bkg"})
bkgebb2016=fMC1.Get("bkgebb2016")
histContainer.addHist(bkgebb2017,{"flavor":"electron","category":"bb","year":"2016","type":"bkg"})
bkgebe2016=fMC1.Get("bkgebe2016")
histContainer.addHist(bkgebb2017,{"flavor":"electron","category":"be","year":"2016","type":"bkg"})

bkgmu=bkgmubb2018.Clone()
bkgmu.Add(bkgmube2018)
bkgmu.Add(bkgmubb2017)
bkgmu.Add(bkgmube2017)
bkgmu.Add(bkgmubb2016)
bkgmu.Add(bkgmube2016)
histContainer.addHist(bkgmu,{"flavor":"muon","category":"combine","year":"combine","type":"bkg"})
bkge=bkgebb2018.Clone()
bkge.Add(bkgebe2018)
bkge.Add(bkgebb2017)
bkge.Add(bkgebe2017)
bkge.Add(bkgebb2016)
bkge.Add(bkgebe2016)
histContainer.addHist(bkge,{"flavor":"electron","category":"combine","year":"combine","type":"bkg"})

bkgmubb=bkgmubb2018.Clone()
bkgmubb.Add(bkgmubb2017)
bkgmubb.Add(bkgmubb2016)
histContainer.addHist(bkgmubb,{"flavor":"muon","category":"bb","year":"combine","type":"bkg"})
bkgmube=bkgmube2018.Clone()
bkgmube.Add(bkgmube2017)
bkgmube.Add(bkgmube2016)
histContainer.addHist(bkgmube,{"flavor":"muon","category":"be","year":"combine","type":"bkg"})
bkgebb=bkgebb2018.Clone()
bkgebb.Add(bkgebb2017)
bkgebb.Add(bkgebb2016)
histContainer.addHist(bkgebb,{"flavor":"electron","category":"bb","year":"combine","type":"bkg"})
bkgebe=bkgebe2018.Clone()
bkgebe.Add(bkgebe2017)
bkgebe.Add(bkgebe2016)
histContainer.addHist(bkgebe,{"flavor":"electron","category":"be","year":"combine","type":"bkg"})
#bng=[60,120,400,600,900,1300,1800,6000]
#bng=numpy.asarray(bng,dtype=numpy.float64)
databe=histContainer.getHist({"flavor":"electron","category":"be","year":"combine","type":"data"})
databe=databe.Rebin(len(bng)-1,"",bng)
reAverage(databe)
databb=histContainer.getHist({"flavor":"electron","category":"bb","year":"combine","type":"data"})
databb=databb.Rebin(len(bng)-1,"",bng)
reAverage(databb)
histContainer.saveHists("dataCollectionDY2016.root")
mcbb=histContainer.getHist({"flavor":"electron","category":"bb","year":"combine","type":"recoMass","source":"combine"})
mcbb=mcbb.Rebin(len(bng)-1,"",bng)
mcbe=histContainer.getHist({"flavor":"electron","category":"be","year":"combine","type":"recoMass","source":"combine"})
mcbe=mcbe.Rebin(len(bng)-1,"",bng)
mcbe_DY=histContainer.getHist({"flavor":"electron","category":"be","year":"combine","type":"recoMass","source":"DrellYan"})
mcbe_DY=mcbe_DY.Rebin(len(bng)-1,"",bng)
mcbe_Other=histContainer.getHist({"flavor":"electron","category":"be","year":"combine","type":"recoMass","source":"Other"})
mcbe_Other=mcbe_Other.Rebin(len(bng)-1,"",bng)
mcbb_DY=histContainer.getHist({"flavor":"electron","category":"bb","year":"combine","type":"recoMass","source":"DrellYan"})
mcbb_DY=mcbb_DY.Rebin(len(bng)-1,"",bng)
mcbb_Other=histContainer.getHist({"flavor":"electron","category":"bb","year":"combine","type":"recoMass","source":"Other"})
mcbb_Other=mcbb_Other.Rebin(len(bng)-1,"",bng)
genbb=histContainer.getHist({"flavor":"electron","category":"bb","year":"combine","type":"genMass","source":"combine"})
genbb=genbb.Rebin(len(bng)-1,"",bng)
genbe=histContainer.getHist({"flavor":"electron","category":"be","year":"combine","type":"genMass","source":"combine"})
genbe=genbe.Rebin(len(bng)-1,"",bng)
for i in range(1,len(bng)):
	val1=databb.GetBinContent(i)
	val2=databe.GetBinContent(i)
	val3=mcbb.GetBinContent(i)
	val4=mcbe.GetBinContent(i)
	val5=genbb.GetBinContent(i)
	val6=genbe.GetBinContent(i)
	val7=mcbe_DY.GetBinContent(i)
	val8=mcbe_Other.GetBinContent(i)
	val9=mcbb_DY.GetBinContent(i)
        val10=mcbb_Other.GetBinContent(i)
	print "data for bb in mass region %f and %f is %f" %(bng[i-1],bng[i],val1)
	print "data for be in mass region %f and %f is %f" %(bng[i-1],bng[i],val2)
	print "reconstructed mass for bb in region %f and %f is %f" %(bng[i-1],bng[i],val3)
	print "reconstructed mass for be in region %f and %f is %f" %(bng[i-1],bng[i],val4)
	print "DY reconstructed mass for be in region %f and %f is %f" %(bng[i-1],bng[i],val7)
        print "Other reconstructed mass for be in region %f and %f is %f" %(bng[i-1],bng[i],val8)
	print "DY reconstructed mass for bb in region %f and %f is %f" %(bng[i-1],bng[i],val9)
        print "Other reconstructed mass for bb in region %f and %f is %f" %(bng[i-1],bng[i],val10)
	print "genrated mass for bb in region %f and %f is %f" %(bng[i-1],bng[i],val5)
        print "genrated mass for be in region %f and %f is %f" %(bng[i-1],bng[i],val6)
#print "reBin"
#bng1=[200,300,400,500,690,900,1250,1610, 2000, 3000]
#bng1=numpy.asarray(bng1,dtype=numpy.float64)
#def reBin1D(hist):
#	name=hist.GetName()
#	hist1=hist.Rebin(len(bng)-1,name,bng)
#	hist1.SetBinContent(len(bng)+1,0)
#	hist1.SetBinError(len(bng)+1,0)
#	hist2=hist1.Rebin(len(bng1)-1,name,bng1)
#	hist=hist2.Clone()
#def reBin2D(hist):
#	name=hist.GetName()
#	title=hist.GetTitle()
#	tmp=ROOT.TH2F(name,title,len(bng1)-1,bng1,len(bng1)-1,bng1)
#	for i in range(hist.GetNbinsX()+1):
#		for j in range(hist.GetNbinsY()+1):
#			if i>10:
#				k=11
#			else:
#				k=i
#			if j>10:
#				l=11
#			else:
#				l=j
#				val1=hist.GetBinContent(i,j)
#				val2=tmp.GetBinContent(k,l)
#				tmp.SetBinContent(k,l,val1+val2)
#
#	hist=tmp.Clone()
#histContainer.processes(reBin2D,["response"])
#histContainer.processes(reBin1D,["recMass","data","bkg"])			
#histContainer.saveHists("dataCollectionDY2016.root")











