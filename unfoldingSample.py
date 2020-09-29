import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy
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
from ROOT import TUnfold
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
def Stacks(processes,lumi,plot,zScale):
        stacks=[]
        for i in range(3):
                stacks.append(TheStack(processes[i],lumi[i],plot,zScale[i]))
        return stacks

datamubb=[]
datamube=[]
data_mu=[Process(Data2016, normalized=True),Process(Data, normalized=True),Process(Data2018, normalized=True)]
plot_mu_bb = getPlot("massPlotBB")
plot_mu_be = getPlot("massPlotBE")
lumi_mu = [36.3*1000,42.135*1000,61.608*1000]
zScaleFac_mu = [zScale2016["muons"],zScale["muons"],zScale2018["muons"]]
for i in range(3):
	datamubb.append(data_mu[i].loadHistogram(plot_mu_bb,lumi_mu[i],zScaleFac_mu[i]))
	datamube.append(data_mu[i].loadHistogram(plot_mu_be,lumi_mu[i],zScaleFac_mu[i]))

bng=[50, 120,150,200,300,400,500,690,900,1250,1610, 2000, 4000, 6070]
#bng=[60,120,400,600,900,1300,1800,6000]
bng=numpy.asarray(bng,dtype=numpy.float64)
datamubb2016=ROOT.TH1D("datamubb2016","datamubb2016",len(bng)-1,bng)
datamubb2017=ROOT.TH1D("datamubb2017","datamubb2017",len(bng)-1,bng)
datamubb2018=ROOT.TH1D("datamubb2018","datamubb2018",len(bng)-1,bng)
datamube2016=ROOT.TH1D("datamube2016","datamube2016",len(bng)-1,bng)
datamube2017=ROOT.TH1D("datamube2017","datamube2017",len(bng)-1,bng)
datamube2018=ROOT.TH1D("datamube2018","datamube2018",len(bng)-1,bng)
for i in range(len(bng)+1):
	val=datamubb[0].GetBinContent(i)
	err=datamubb[0].GetBinError(i)
	datamubb2016.SetBinContent(i,val)
	datamubb2016.SetBinError(i,err)
	val=datamubb[1].GetBinContent(i)
        err=datamubb[1].GetBinError(i)
        datamubb2017.SetBinContent(i,val)
        datamubb2017.SetBinError(i,err)
	val=datamubb[2].GetBinContent(i)
        err=datamubb[2].GetBinError(i)
        datamubb2018.SetBinContent(i,val)
        datamubb2018.SetBinError(i,err)
	val=datamube[0].GetBinContent(i)
        err=datamube[0].GetBinError(i)
        datamube2016.SetBinContent(i,val)
        datamube2016.SetBinError(i,err)
        val=datamube[1].GetBinContent(i)
        err=datamube[1].GetBinError(i)
        datamube2017.SetBinContent(i,val)
        datamube2017.SetBinError(i,err)
        val=datamube[2].GetBinContent(i)
        err=datamube[2].GetBinError(i)
        datamube2018.SetBinContent(i,val)
        datamube2018.SetBinError(i,err)

dataebb=[]
dataebe=[]
data_e=[Process(Data2016, normalized=True),Process(Data, normalized=True),Process(Data2018, normalized=True)]
plot_e_bb = getPlot("massPlotEleBB")
plot_e_be = getPlot("massPlotEleBE")
lumi_e = [35.9*1000,41.529*1000,59.97*1000]
zScaleFac_e = [zScale2016["electrons"],zScale["electrons"],zScale2018["electrons"]]
for i in range(3):
        dataebb.append(data_e[i].loadHistogram(plot_e_bb,lumi_e[i],zScaleFac_e[i][0]))
        dataebe.append(data_e[i].loadHistogram(plot_e_be,lumi_e[i],zScaleFac_e[i][0]))

dataebb2016=ROOT.TH1D("dataebb2016","dataebb2016",len(bng)-1,bng)
dataebb2017=ROOT.TH1D("dataebb2017","dataebb2017",len(bng)-1,bng)
dataebb2018=ROOT.TH1D("dataebb2018","dataebb2018",len(bng)-1,bng)
dataebe2016=ROOT.TH1D("dataebe2016","dataebe2016",len(bng)-1,bng)
dataebe2017=ROOT.TH1D("dataebe2017","dataebe2017",len(bng)-1,bng)
dataebe2018=ROOT.TH1D("dataebe2018","dataebe2018",len(bng)-1,bng)
for i in range(len(bng)+1):
        val=dataebb[0].GetBinContent(i)
        err=dataebb[0].GetBinError(i)
        dataebb2016.SetBinContent(i,val)
        dataebb2016.SetBinError(i,err)
        val=dataebb[1].GetBinContent(i)
        err=dataebb[1].GetBinError(i)
        dataebb2017.SetBinContent(i,val)
	dataebb2017.SetBinError(i,err)
	val=dataebb[2].GetBinContent(i)
        err=dataebb[2].GetBinError(i)
        dataebb2018.SetBinContent(i,val)
	dataebb2018.SetBinError(i,err)
	val=dataebe[0].GetBinContent(i)
        err=dataebe[0].GetBinError(i)
        dataebe2016.SetBinContent(i,val)
        dataebe2016.SetBinError(i,err)
        val=dataebe[1].GetBinContent(i)
        err=dataebe[1].GetBinError(i)
        dataebe2017.SetBinContent(i,val)
	dataebe2017.SetBinError(i,err)
        val=dataebe[2].GetBinContent(i)
        err=dataebe[2].GetBinError(i)
        dataebe2018.SetBinContent(i,val)
	dataebe2018.SetBinError(i,err)

backgrounds=["Wjets","Other"]

f=ROOT.TFile("unfoldingData.root","RECREATE")
datamubb2016.Write()
datamubb2017.Write()
datamubb2018.Write()
datamube2016.Write()
datamube2017.Write()
datamube2018.Write()
dataebb2016.Write()
dataebb2017.Write()
dataebb2018.Write()
dataebe2016.Write()
dataebe2017.Write()
dataebe2018.Write()
f.Close()

print plot_mu_be.useJets
print plot_e_be.useJets
print plot_e_bb.fileName
print plot_e_be.fileName
eventCounts_e = totalNumberOfGeneratedEvents(path,False)
eventCounts_mu = totalNumberOfGeneratedEvents(path,True)
print eventCounts_e
negWeights_mu = negWeightFractions(path,True)
negWeights_e = negWeightFractions(path,False)
processes_mu2016=[Process(getattr(Backgrounds2016,"Jets"),eventCounts_mu,negWeights_mu,normalized=True),Process(getattr(Backgrounds2016,"Other"),eventCounts_mu,negWeights_mu)]
processes_e2016=[Process(getattr(Backgrounds2016,"Jets"),eventCounts_e,negWeights_e,normalized=True),Process(getattr(Backgrounds2016,"OtherEle"),eventCounts_e,negWeights_e)]
processes_mu2017=[Process(getattr(Backgrounds,"Jets"),eventCounts_mu,negWeights_mu,normalized=True),Process(getattr(Backgrounds,"Other"),eventCounts_mu,negWeights_mu)]
processes_e2017=[Process(getattr(Backgrounds,"Jets"),eventCounts_e,negWeights_e,normalized=True),Process(getattr(Backgrounds,"Other"),eventCounts_e,negWeights_e)]
processes_mu2018=[Process(getattr(Backgrounds2018,"Jets"),eventCounts_mu,negWeights_mu,normalized=True),Process(getattr(Backgrounds2018,"Other"),eventCounts_mu,negWeights_mu)]
processes_e2018=[Process(getattr(Backgrounds2018,"Jets"),eventCounts_e,negWeights_e,normalized=True),Process(getattr(Backgrounds2018,"Other"),eventCounts_e,negWeights_e)]
stackmu_bb2016 = TheStack(processes_mu2016,lumi_mu[0],plot_mu_bb,zScaleFac_mu[0])
stackmu_bb2017 = TheStack(processes_mu2017,lumi_mu[1],plot_mu_bb,zScaleFac_mu[1])
stackmu_bb2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_bb,zScaleFac_mu[2])
stackmu_be2016 = TheStack(processes_mu2016,lumi_mu[0],plot_mu_be,zScaleFac_mu[0])
stackmu_be2017 = TheStack(processes_mu2017,lumi_mu[1],plot_mu_be,zScaleFac_mu[1])
stackmu_be2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_be,zScaleFac_mu[2])
stacke_bb2016 = TheStack(processes_mu2016,lumi_mu[0],plot_mu_bb,zScaleFac_e[0][1])
stacke_bb2017 = TheStack(processes_mu2017,lumi_mu[1],plot_mu_bb,zScaleFac_e[1][1])
stacke_bb2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_bb,zScaleFac_e[2][1])
stacke_be2016 = TheStack(processes_mu2016,lumi_mu[0],plot_mu_be,zScaleFac_e[0][2])
stacke_be2017 = TheStack(processes_mu2017,lumi_mu[1],plot_mu_be,zScaleFac_e[1][2])
stacke_be2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_be,zScaleFac_e[2][2])
bkgmubb=[stackmu_bb2016.theHistogram,stackmu_bb2017.theHistogram,stackmu_bb2018.theHistogram]
bkgmube=[stackmu_be2016.theHistogram,stackmu_be2017.theHistogram,stackmu_be2018.theHistogram]
bkgebb=[stacke_bb2016.theHistogram,stacke_bb2017.theHistogram,stacke_bb2018.theHistogram]
bkgebe=[stacke_be2016.theHistogram,stacke_be2017.theHistogram,stacke_be2018.theHistogram]
bkgebb2016=ROOT.TH1D("bkgebb2016","bkgebb2016",len(bng)-1,bng)
bkgebb2017=ROOT.TH1D("bkgebb2017","bkgebb2017",len(bng)-1,bng)
bkgebb2018=ROOT.TH1D("bkgebb2018","bkgebb2018",len(bng)-1,bng)
bkgebe2016=ROOT.TH1D("bkgebe2016","bkgebe2016",len(bng)-1,bng)
bkgebe2017=ROOT.TH1D("bkgebe2017","bkgebe2017",len(bng)-1,bng)
bkgebe2018=ROOT.TH1D("bkgebe2018","bkgebe2018",len(bng)-1,bng)
bkgmubb2016=ROOT.TH1D("bkgmubb2016","bkgmubb2016",len(bng)-1,bng)
bkgmubb2017=ROOT.TH1D("bkgmubb2017","bkgmubb2017",len(bng)-1,bng)
bkgmubb2018=ROOT.TH1D("bkgmubb2018","bkgmubb2018",len(bng)-1,bng)
bkgmube2016=ROOT.TH1D("bkgmube2016","bkgmube2016",len(bng)-1,bng)
bkgmube2017=ROOT.TH1D("bkgmube2017","bkgmube2017",len(bng)-1,bng)
bkgmube2018=ROOT.TH1D("bkgmube2018","bkgmube2018",len(bng)-1,bng)
for i in range(len(bng)+1):
	val1=0
	val2=0
        val=bkgmubb[0].GetBinContent(i)
        err=bkgmubb[0].GetBinError(i)
        bkgmubb2016.SetBinContent(i,val)
        bkgmubb2016.SetBinError(i,err)
        val=bkgmubb[1].GetBinContent(i)
        err=bkgmubb[1].GetBinError(i)
        bkgmubb2017.SetBinContent(i,val)
        bkgmubb2017.SetBinError(i,err)
        val=bkgmubb[2].GetBinContent(i)
        err=bkgmubb[2].GetBinError(i)
        bkgmubb2018.SetBinContent(i,val)
        bkgmubb2018.SetBinError(i,err)
        val=bkgmube[0].GetBinContent(i)
        err=bkgmube[0].GetBinError(i)
        bkgmube2016.SetBinContent(i,val)
        bkgmube2016.SetBinError(i,err)
        val=bkgmube[1].GetBinContent(i)
        err=bkgmube[1].GetBinError(i)
        bkgmube2017.SetBinContent(i,val)
        bkgmube2017.SetBinError(i,err)
        val=bkgmube[2].GetBinContent(i)
        err=bkgmube[2].GetBinError(i)
        bkgmube2018.SetBinContent(i,val)
        bkgmube2018.SetBinError(i,err)
	val=bkgebb[0].GetBinContent(i)
	val1+=val
        err=bkgebb[0].GetBinError(i)
        bkgebb2016.SetBinContent(i,val)
        bkgebb2016.SetBinError(i,err)
        val=bkgebb[1].GetBinContent(i)
	val1+=val
        err=bkgebb[1].GetBinError(i)
        bkgebb2017.SetBinContent(i,val)
        bkgebb2017.SetBinError(i,err)
        val=bkgebb[2].GetBinContent(i)
        err=bkgebb[2].GetBinError(i)
	val1+=val
        bkgebb2018.SetBinContent(i,val)
        bkgebb2018.SetBinError(i,err)
        val=bkgebe[0].GetBinContent(i)
        err=bkgebe[0].GetBinError(i)
	val2+=val
        bkgebe2016.SetBinContent(i,val)
        bkgebe2016.SetBinError(i,err)
        val=bkgebe[1].GetBinContent(i)
        err=bkgebe[1].GetBinError(i)
	val2+=val
        bkgebe2017.SetBinContent(i,val)
        bkgebe2017.SetBinError(i,err)
        val=bkgebe[2].GetBinContent(i)
        err=bkgebe[2].GetBinError(i)
	val2+=val
        bkgebe2018.SetBinContent(i,val)
        bkgebe2018.SetBinError(i,err)
	print val1
	print val2
f=ROOT.TFile("unfoldingMC.root","RECREATE")
bkgmubb2016.Write()
bkgmubb2017.Write()
bkgmubb2018.Write()
bkgmube2016.Write()
bkgmube2017.Write()
bkgmube2018.Write()
bkgebb2016.Write()
bkgebb2017.Write()
bkgebb2018.Write()
bkgebe2016.Write()
bkgebe2017.Write()
bkgebe2018.Write()
f.Close()





