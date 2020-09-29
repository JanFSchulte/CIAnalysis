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
def average(hist):
        for i in range(hist.GetNbinsX()+2):
                val=hist.GetBinContent(i)
                err=hist.GetBinError(i)
                wid=hist.GetBinWidth(i)
                val=val/wid
                err=err/wid
                if i==hist.GetNbinsX()+1:
                        print wid
                hist.SetBinContent(i,val)
                hist.SetBinError(i,err)
def reAverage(hist):
        for i in range(hist.GetNbinsX()+2):
                val=hist.GetBinContent(i)
                err=hist.GetBinError(i)
                wid=hist.GetBinWidth(i)
                val=val*wid
                err=err*wid
                hist.SetBinContent(i,val)
                hist.SetBinError(i,err)

plot_mu_bb = getPlot("massPlotBB")
plot_mu_be = getPlot("massPlotBE")
lumi_mu = [36.3*1000,42.135*1000,61.608*1000]
zScaleFac_mu = [zScale2016["muons"],zScale["muons"],zScale2018["muons"]]
bng=[50, 120,150,200,300,400,500,690,900,1250,1610, 2000, 3000,4000, 6070]
bng1=[50, 120,150,200,300,400,500,690,900,1250,1610, 1970, 3010,3970, 6070]
bng=numpy.asarray(bng,dtype=numpy.float64)
bng1=numpy.asarray(bng1,dtype=numpy.float64)
plot_e_bb = getPlot("massPlotEleBB")
plot_e_be = getPlot("massPlotEleBE")
lumi_e = [35.9*1000,41.529*1000,59.97*1000]
zScaleFac_e = [zScale2016["electrons"],zScale["electrons"],zScale2018["electrons"]]
eventCounts_e = totalNumberOfGeneratedEvents(path,False)
eventCounts_mu = totalNumberOfGeneratedEvents(path,True)
negWeights_mu = negWeightFractions(path,True)
negWeights_e = negWeightFractions(path,False)
processes_mu2016=[Process(getattr(Backgrounds2016,"Jets"),eventCounts_mu,negWeights_mu,normalized=True)]
processes_e2016=[Process(getattr(Backgrounds2016,"Jets"),eventCounts_e,negWeights_e,normalized=True)]
processes_mu2017=[Process(getattr(Backgrounds,"Jets"),eventCounts_mu,negWeights_mu,normalized=True)]
processes_e2017=[Process(getattr(Backgrounds,"Jets"),eventCounts_e,negWeights_e,normalized=True)]
processes_mu2018=[Process(getattr(Backgrounds2018,"Jets"),eventCounts_mu,negWeights_mu,normalized=True)]
processes_e2018=[Process(getattr(Backgrounds2018,"Jets"),eventCounts_e,negWeights_e,normalized=True)]
stackmu_bb2016 = TheStack(processes_mu2016,lumi_mu[0],plot_mu_bb,zScaleFac_mu[0])
stackmu_bb2017 = TheStack(processes_mu2017,lumi_mu[1],plot_mu_bb,zScaleFac_mu[1])
stackmu_bb2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_bb,zScaleFac_mu[2])
stackmu_be2016 = TheStack(processes_mu2016,lumi_mu[0],plot_mu_be,zScaleFac_mu[0])
stackmu_be2017 = TheStack(processes_mu2017,lumi_mu[1],plot_mu_be,zScaleFac_mu[1])
stackmu_be2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_be,zScaleFac_mu[2])
stacke_bb2016 = TheStack(processes_e2016,lumi_e[0],plot_e_bb,zScaleFac_e[0][1])
stacke_bb2017 = TheStack(processes_e2017,lumi_e[1],plot_e_bb,zScaleFac_e[1][1])
stacke_bb2018 = TheStack(processes_e2018,lumi_e[2],plot_e_bb,zScaleFac_e[2][1])
stacke_be2016 = TheStack(processes_e2016,lumi_e[0],plot_e_be,zScaleFac_e[0][2])
stacke_be2017 = TheStack(processes_e2017,lumi_e[1],plot_e_be,zScaleFac_e[1][2])
stackmu_be2018 = TheStack(processes_mu2018,lumi_mu[2],plot_mu_be,zScaleFac_mu[2])
stacke_bb2016 = TheStack(processes_e2016,lumi_e[0],plot_e_bb,zScaleFac_e[0][1])
stacke_bb2017 = TheStack(processes_e2017,lumi_e[1],plot_e_bb,zScaleFac_e[1][1])
stacke_bb2018 = TheStack(processes_e2018,lumi_e[2],plot_e_bb,zScaleFac_e[2][1])
stacke_be2016 = TheStack(processes_e2016,lumi_e[0],plot_e_be,zScaleFac_e[0][2])
stacke_be2017 = TheStack(processes_e2017,lumi_e[1],plot_e_be,zScaleFac_e[1][2])
stacke_be2018 = TheStack(processes_e2018,lumi_e[2],plot_e_be,zScaleFac_e[2][2])
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
reAverage(bkgebb[0])
bkgebb0=bkgebb[0].Rebin(len(bng1)-1,"",bng1)
reAverage(bkgebb[1])
bkgebb1=bkgebb[1].Rebin(len(bng1)-1,"",bng1)
reAverage(bkgebb[2])
bkgebb2=bkgebb[2].Rebin(len(bng1)-1,"",bng1)
reAverage(bkgebe[0])
bkgebe0=bkgebe[0].Rebin(len(bng1)-1,"",bng1)
reAverage(bkgebe[1])
bkgebe1=bkgebe[1].Rebin(len(bng1)-1,"",bng1)
reAverage(bkgebe[2])
bkgebe2=bkgebe[2].Rebin(len(bng1)-1,"",bng1)
for i in range(len(bng)+1):
        val1=0
        val2=0
        #val=bkgmubb[0].GetBinContent(i)
        #err=bkgmubb[0].GetBinError(i)
        #bkgmubb2016.SetBinContent(i,val)
        #bkgmubb2016.SetBinError(i,err)
        #val=bkgmubb[1].GetBinContent(i)
        #err=bkgmubb[1].GetBinError(i)
        #bkgmubb2017.SetBinContent(i,val)
        #bkgmubb2017.SetBinError(i,err)
        #val=bkgmubb[2].GetBinContent(i)
        #err=bkgmubb[2].GetBinError(i)
        #bkgmubb2018.SetBinContent(i,val)
        #bkgmubb2018.SetBinError(i,err)
        #val=bkgmube[0].GetBinContent(i)
        #err=bkgmube[0].GetBinError(i)
        #bkgmube2016.SetBinContent(i,val)
        #bkgmube2016.SetBinError(i,err)
        #val=bkgmube[1].GetBinContent(i)
        #err=bkgmube[1].GetBinError(i)
        #bkgmube2017.SetBinContent(i,val)
        #bkgmube2017.SetBinError(i,err)
        #val=bkgmube[2].GetBinContent(i)
        #err=bkgmube[2].GetBinError(i)
        #bkgmube2018.SetBinContent(i,val)
        #bkgmube2018.SetBinError(i,err)
        val=bkgebb0.GetBinContent(i)
        err=bkgebb0.GetBinError(i)
        val1+=val
        bkgebb2016.SetBinContent(i,val)
	bkgebb2016.SetBinError(i,err)
        val=bkgebb1.GetBinContent(i)
        err=bkgebb1.GetBinError(i)
        val1+=val
        bkgebb2017.SetBinContent(i,val)
        bkgebb2017.SetBinError(i,err)
        val=bkgebb2.GetBinContent(i)
        err=bkgebb2.GetBinError(i)
        val1+=val
        bkgebb2018.SetBinContent(i,val)
        bkgebb2018.SetBinError(i,err)
        val=bkgebe0.GetBinContent(i)
        err=bkgebe0.GetBinError(i)
        val2+=val
        bkgebe2016.SetBinContent(i,val)
        bkgebe2016.SetBinError(i,err)
        val=bkgebe1.GetBinContent(i)
        err=bkgebe1.GetBinError(i)
        val2+=val
        bkgebe2017.SetBinContent(i,val)
        bkgebe2017.SetBinError(i,err)
        val=bkgebe2.GetBinContent(i)
        err=bkgebe2.GetBinError(i)
        val2+=val
        bkgebe2018.SetBinContent(i,val)
        bkgebe2018.SetBinError(i,err)
        print val1
        print val2
average(bkgebb2016)
average(bkgebb2017)
average(bkgebb2018)
average(bkgebe2016)
average(bkgebe2017)
average(bkgebe2018)
f=ROOT.TFile("unfoldingMC_V3.root","RECREATE")
bkgebb2016.Write()
bkgebb2017.Write()
bkgebb2018.Write()
bkgebe2016.Write()
bkgebe2017.Write()
bkgebe2018.Write()
f.Close()

