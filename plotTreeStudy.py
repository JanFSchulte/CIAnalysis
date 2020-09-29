import ROOT
import os
import numpy
import subprocess
import root_numpy
import math
import copy


def loopOver(path):
	pathList=[]
	for filename in os.listdir(path):
		pathNew=path+"/"+filename
		if os.path.isfile(pathNew) and pathNew.endswith(".root"):
			pathList.append(pathNew)
		elif os.path.isdir(pathNew): 
			return loopOver(pathNew)	
		return pathList

def getBkg(path,name):
	pathList=[]
	for filename in os.listdir(path):
		pathNew=path+"/"+filename
		if os.path.isdir(pathNew) and filename.startswith(name):
			pathList+=loopOver(pathNew)
	return pathList


ROOT.gStyle.SetOptStat(0)
path="/mnt/hadoop/store/user/minxi"
pathList=getBkg(path,"ZToMuMu")
pathList+=getBkg(path,"DY")
pathList2016=[]
pathList2017=[]
for path in pathList:
	if(path.find("2016")==-1):
		pathList2016.append(path)
	else:
		pathList2017.append(path)

#chain=ROOT.TChain("SimpleNtupler")
mass2016=y2016=pt2016=CS2016=CS_label2016=numpy.ones(0)
mass2017=y2017=pt2017=CS2017=CS_label2017=numpy.ones(0)
for path in pathList2016:
	#path=path.replace("/mnt/hadoop/","root://xrootd.rcac.purdue.edu//")
	splitPath=path.split("/")
	pathMap="/home/yang1452/lepratio/CIAnalysis/crabFile/"+splitPath[-3]+"/"+splitPath[-2]
	if not os.path.exists(pathMap):
		os.makedirs(pathMap)
	cmd="cp "+path+" "+pathMap
	os.system(cmd)
	pathMap=pathMap+"/"+splitPath[-1]
	f2016=ROOT.TFile.Open(pathMap)
	tree2016=f2016.Get("SimpleNtupler/t")
	mass2016=numpy.concatenate([mass2016,root_numpy.tree2array(tree2016,"dil_mass")])
	y2016=numpy.concatenate([y2016,root_numpy.tree2array(tree2016,"dil_rap")])
	pt2016=numpy.concatenate([pt2016,root_numpy.tree2array(tree2016,"dil_pt")])
	CS2016=numpy.concatenate([CS2016,root_numpy.tree2array(tree2016,"cos_cs")])
	CS_label2016=numpy.concatenate([CS_label2016,root_numpy.tree2array(tree2016,"CS_label")])



for path in pathList2017:
        #path=path.replace("/mnt/hadoop/","root://xrootd.rcac.purdue.edu//")
        splitPath=path.split("/")
        pathMap="/home/yang1452/lepratio/CIAnalysis/crabFile/"+splitPath[-3]+"/"+splitPath[-2]
        if not os.path.exists(pathMap):
                os.makedirs(pathMap)
        cmd="cp "+path+" "+pathMap
        os.system(cmd)
        pathMap=pathMap+"/"+splitPath[-1]
        f2017=ROOT.TFile.Open(pathMap)
        tree2017=f2017.Get("SimpleNtupler/t")
        mass2017=numpy.concatenate([mass2017,root_numpy.tree2array(tree2017,"dil_mass")])
        y2017=numpy.concatenate([y2017,root_numpy.tree2array(tree2017,"dil_rap")])
        pt2017=numpy.concatenate([pt2017,root_numpy.tree2array(tree2017,"dil_pt")])
        CS2017=numpy.concatenate([CS2017,root_numpy.tree2array(tree2017,"cos_cs")])
        CS_label2017=numpy.concatenate([CS_label2017,root_numpy.tree2array(tree2017,"CS_label")])





#h_CS800to1000=ROOT.TH1D("h_CS","CS angle",100,-1.,1.)
#h_CS800to1000.SetLineColor(2)
#h_CS800to1000.GetXaxis().SetTitle("cos #theta")
#h_CS5500to6000=ROOT.TH1D("h_CS","CS angle",100,-1.,1.)
#h_CS5500to6000.SetLineColor(3)
#h_CS5500to6000T=ROOT.TH1D("h_CS","CS angle",100,-1.,1.)
#h_CS5500to6000T.SetLineColor(4)
bng=[i for i in range(200,400,20)]+[i for i in range(400,800,40)]+[i for i in range(800,1400,60)]+[i for i in range(1400,2300,90)]+[i for i in range(2300,3500,120)]+[i for i in range(3500,4500,100)]+[i for i in range(4500,6000,150)]
bng=numpy.asarray(bng,dtype=numpy.float64)
CS_2016true=copy.deepcopy(CS2016)
CS_2016true[CS_label2016==0]=-CS_2016true[CS_label2016==0]
CS_2017true=copy.deepcopy(CS2017)
CS_2017true[CS_label2017==0]=-CS_2017true[CS_label2017==0]
h_prob_yVSmassn2016=ROOT.TH2D("h_prob_yVSmassn2016","y_vs_mass",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmassd2016=ROOT.TH2D("h_prob_yVSmassd2016","y_vs_mass",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmass2016=ROOT.TH2D("h_prob_yVSmass2016","Probability to identify quark with a correct direction",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmassErr2016=ROOT.TH2D("h_prob_yVSmassErr2016","Relative Error",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmassn2017=ROOT.TH2D("h_prob_yVSmassn2017","y_vs_mass",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmassd2017=ROOT.TH2D("h_prob_yVSmassd2017","y_vs_mass",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmass2017=ROOT.TH2D("h_prob_yVSmass2017","Probability to identify quark with a correct direction",len(bng)-1,bng,30,0,2.4)
h_prob_yVSmassErr2017=ROOT.TH2D("h_prob_yVSmassErr2017","Relative Error",len(bng)-1,bng,30,0,2.4)

h_prob_yVSmassn2016.FillN(len(y2016[CS_label2016==1]),mass2016[CS_label2016==1],numpy.absolute(y2016[CS_label2016==1]),numpy.ones(len(y2016[CS_label2016==1])))
h_prob_yVSmassn2016.Sumw2()
h_prob_yVSmassd2016.FillN(len(y2016),mass2016,numpy.absolute(y2016),numpy.ones(len(y2016)))
h_prob_yVSmassd2016.Sumw2()
h_prob_yVSmass2016.Divide(h_prob_yVSmassn2016,h_prob_yVSmassd2016)
h_prob_yVSmassn2017.FillN(len(y2017[CS_label2017==1]),mass2017[CS_label2017==1],numpy.absolute(y2017[CS_label2017==1]),numpy.ones(len(y2017[CS_label2017==1])))
h_prob_yVSmassn2017.Sumw2()
h_prob_yVSmassd2017.FillN(len(y2017),mass2017,numpy.absolute(y2017),numpy.ones(len(y2017)))
h_prob_yVSmassd2017.Sumw2()
h_prob_yVSmass2017.Divide(h_prob_yVSmassn2017,h_prob_yVSmassd2017)

for i in range(0,h_prob_yVSmass2016.GetNbinsX()):
	for j in range(0,h_prob_yVSmass2016.GetNbinsY()):
		if h_prob_yVSmassn2016.GetBinContent(i+1,j+1)!=0:
			err=h_prob_yVSmass2016.GetBinError(i+1,j+1)/h_prob_yVSmass2016.GetBinContent(i+1,j+1)
		else:
			err=0		
		h_prob_yVSmassErr2016.SetBinContent(i+1,j+1,err)

for i in range(0,h_prob_yVSmass2017.GetNbinsX()):
        for j in range(0,h_prob_yVSmass2017.GetNbinsY()):
                if h_prob_yVSmassn2017.GetBinContent(i+1,j+1)!=0:
                        err=h_prob_yVSmass2017.GetBinError(i+1,j+1)/h_prob_yVSmass2017.GetBinContent(i+1,j+1)
                else:
                        err=0   
                h_prob_yVSmassErr2017.SetBinContent(i+1,j+1,err)

c1=ROOT.TCanvas("c1","c1",800,800)
c1.Divide(2,2)
c1.cd(1)
h_prob_yVSmass2016.GetXaxis().SetTitle("mll[GeV]")
h_prob_yVSmass2016.GetYaxis().SetTitle("y")
h_prob_yVSmass2016.Draw("COLZ")
c1.cd(2)
h_prob_yVSmassErr2016.GetXaxis().SetTitle("mll[GeV]")
h_prob_yVSmassErr2016.GetYaxis().SetTitle("y")
h_prob_yVSmassErr2016.Draw("COLZ")
c1.cd(3)
h_prob_yVSmass2017.GetXaxis().SetTitle("mll[GeV]")
h_prob_yVSmass2017.GetYaxis().SetTitle("y")
h_prob_yVSmass2017.Draw("COLZ")
c1.cd(4)
h_prob_yVSmassErr2017.GetXaxis().SetTitle("mll[GeV]")
h_prob_yVSmassErr2017.GetYaxis().SetTitle("y")
h_prob_yVSmassErr2017.Draw("COLZ")
c1.Print("Study/Result1.pdf")









h_CS1_16=ROOT.TH1D("h_CS1_16","CS angle",100,-1.,1.)
h_CS2_16=ROOT.TH1D("h_CS2_16","CS angle",100,-1.,1.)
h_CS3_16=ROOT.TH1D("h_CS3_16","CS angle",100,-1.,1.)
h_CS4_16=ROOT.TH1D("h_CS4_16","CS angle",100,-1.,1.)
h_CS5_16=ROOT.TH1D("h_CS5_16","CS angle",100,-1.,1.)
h_CS1_16t=ROOT.TH1D("h_CS1_16t","CS angle",100,-1.,1.)
h_CS2_16t=ROOT.TH1D("h_CS2_16t","CS angle",100,-1.,1.)
h_CS3_16t=ROOT.TH1D("h_CS3_16t","CS angle",100,-1.,1.)
h_CS4_16t=ROOT.TH1D("h_CS4_16t","CS angle",100,-1.,1.)
h_CS5_16t=ROOT.TH1D("h_CS5_16t","CS angle",100,-1.,1.)
h_CS1_17=ROOT.TH1D("h_CS1_17","CS angle",100,-1.,1.)
h_CS2_17=ROOT.TH1D("h_CS2_17","CS angle",100,-1.,1.)
h_CS3_17=ROOT.TH1D("h_CS3_17","CS angle",100,-1.,1.)
h_CS4_17=ROOT.TH1D("h_CS4_17","CS angle",100,-1.,1.)
h_CS5_17=ROOT.TH1D("h_CS5_17","CS angle",100,-1.,1.)
h_CS1_17t=ROOT.TH1D("h_CS1_17t","CS angle",100,-1.,1.)
h_CS2_17t=ROOT.TH1D("h_CS2_17t","CS angle",100,-1.,1.)
h_CS3_17t=ROOT.TH1D("h_CS3_17t","CS angle",100,-1.,1.)
h_CS4_17t=ROOT.TH1D("h_CS4_17t","CS angle",100,-1.,1.)
h_CS5_17t=ROOT.TH1D("h_CS5_17t","CS angle",100,-1.,1.)
#h_CS800to1000.SetLineColor(2)
#h_CS800to1000.GetXaxis().SetTitle("cos #theta")
#h_CS5500to6000=ROOT.TH1D("h_CS","CS angle",100,-1.,1.)
#h_CS5500to6000.SetLineColor(3)
#h_CS5500to6000T=ROOT.TH1D("h_CS","CS angle",100,-1.,1.)
#h_CS5500to6000T.SetLineColor(4)
h_16=[h_CS1_16,h_CS2_16,h_CS3_16,h_CS4_16,h_CS5_16]
h_16t=[h_CS1_16t,h_CS2_16t,h_CS3_16t,h_CS4_16t,h_CS5_16t]
h_17=[h_CS1_17,h_CS2_17,h_CS3_17,h_CS4_17,h_CS5_17]
h_17t=[h_CS1_17t,h_CS2_17t,h_CS3_17t,h_CS4_17t,h_CS5_17t]
massLimit=[1000.,2000.,3000.,4000.,5000.,6000.]
h_CS1_16t.FillN(len(CS2016[(mass2016<massLimit[1]) & (mass2016>massLimit[0]) & (y2016<1.)]),CS2016[(mass2016<massLimit[1]) & (mass2016>massLimit[0]) & (y2016<1.)],numpy.ones(len(CS2016[(mass2016<massLimit[1]) & (mass2016>massLimit[0]) & (y2016<1.)])))
c4=ROOT.TCanvas("c4","c4",800,800)
h_CS1_16t.Draw()
c4.Print("Study/result4.pdf")
i=0
for h in h_16t:
	h_16[i].FillN(len(CS2016[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)]),CS2016[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)],numpy.ones(len(CS2016[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)])))
	#print len(CS2016[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)])
	h_16[i].Scale(1./len(CS2016[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)]))
	h.FillN(len(CS_2016true[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)]),CS_2016true[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)],numpy.ones(len(CS_2016true[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)])))
        h.Scale(1./len(CS_2016true[(mass2016<massLimit[i+1]) & (mass2016>massLimit[i]) & (y2016<1.)]))
 	h_17[i].FillN(len(CS2017[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)]),CS2017[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)],numpy.ones(len(CS2017[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)])))
        h_17[i].Scale(1./len(CS2017[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)]))
        h_17t[i].FillN(len(CS_2017true[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)]),CS_2017true[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)],numpy.ones(len(CS_2017true[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)])))
        h_17t[i].Scale(1./len(CS_2017true[(mass2017<massLimit[i+1]) & (mass2017>massLimit[i]) & (y2017<1.)]))
	i+=1


c2=ROOT.TCanvas("c2","c2",800,800)
c2.Divide(2,3)
c2.cd(1)
lyt1=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
lyt1.SetTextFont(28)
leg1=["1000 to 2000 GeV","2000 to 3000 GeV"," 3000 to 4000 GeV"," 4000 to 5000 GeV"," 5000 to 6000 GeV"]
#h_plott1=ROOT.TH1D("h_plott1","2016 truth of the CS angle",100,-1.,1.)
#h_plott1.GetXaxis().SetTitle("cos #theta")
#h_plott1.Draw()
h_16t[0].SetTitle("2016 truth of the CS angle")
for i in range(len(h_16)):
	h_16t[i].SetLineColor(i+2)
	lyt1.AddEntry(h_16t[i],leg1[i])
	if i>0:
		h_16t[i].Draw("samehist")
	else:
		h_16t[i].Draw("hist")
lyt1.Draw()

c2.cd(2)
lyt2=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
lyt2.SetTextFont(28)
#h_plott2=ROOT.TH1D("h_plott2","2017 truth of the CS angle",100,-1.,1.)
#h_plott2.GetXaxis().SetTitle("cos #theta")
#h_plott2.Draw()
h_17t[0].SetTitle("2017 truth of the CS angle")
for i in range(len(h_17)):
        h_17t[i].SetLineColor(i+2)
        lyt2.AddEntry(h_17t[i],leg1[i])
	if i>0:
        	h_17t[i].Draw("samehist")
	else:
		h_17t[i].Draw("hist")
lyt2.Draw()

c2.cd(3)
lyt3=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
lyt3.SetTextFont(28)
#h_plott3=ROOT.TH1D("h_plott3","Comparsion of the CS truth between 2016 and 2017 at 2 TeV to 3 TeV",100,-1.,1.)
#h_plott3.GetXaxis().SetTitle("cos #theta")
h_16t[1].SetTitle("Comparsion of the CS truth between 2016 and 2017 at 2 TeV to 3 TeV")
#h_plott3.Draw()
h_16t[1].SetLineColor(2)
h_17t[1].SetLineColor(3)
lyt3.AddEntry(h_16t[1],"2016")
lyt3.AddEntry(h_17t[1],"2017")
h_16t[1].Draw("hist")
h_17t[1].Draw("samehist")
lyt3.Draw()


c2.cd(4)
lyt4=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
lyt4.SetTextFont(28)
#h_plott4=ROOT.TH1D("h_plott4","Comparsion of the CS truth between 2016 and 2017 at 3 TeV to 4 TeV",100,-1.,1.)
#h_plott4.GetXaxis().SetTitle("cos #theta")
#h_plott4.Draw()
h_16t[2].SetTitle("Comparsion of the CS truth between 2016 and 2017 at 3 TeV to 4 TeV")
h_16t[2].SetLineColor(2)
h_17t[2].SetLineColor(3)
lyt4.AddEntry(h_16t[2],"2016")
lyt4.AddEntry(h_17t[2],"2017")
h_16t[2].Draw("hist")
h_17t[2].Draw("samehist")
lyt4.Draw()

c2.cd(5)
lyt5=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
lyt5.SetTextFont(28)
#h_plott5=ROOT.TH1D("h_plott5","Comparsion of the CS truth between 2016 and 2017 at 4 TeV to 5 TeV",100,-1.,1.)
#h_plott5.GetXaxis().SetTitle("cos #theta")
h_16t[3].SetTitle("Comparsion of the CS truth between 2016 and 2017 at 4 TeV to 5 TeV")
#h_plott5.Draw()
h_16t[3].SetLineColor(2)
h_17t[3].SetLineColor(3)
lyt5.AddEntry(h_16t[3],"2016")
lyt5.AddEntry(h_17t[3],"2017")
h_16t[3].Draw("hist")
h_17t[3].Draw("samehist")
lyt5.Draw()

c2.cd(6)
lyt6=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
lyt6.SetTextFont(28)
#h_plott6=ROOT.TH1D("h_plott6","Comparsion of the CS truth between 2016 and 2017 at 5 TeV to 6 TeV",100,-1.,1.)
#h_plott6.GetXaxis().SetTitle("cos #theta")
h_16t[4].SetTitle("Comparsion of the CS truth between 2016 and 2017 at 5 TeV to 6 TeV")
#h_plott6.Draw()
h_16t[4].SetLineColor(2)
h_17t[4].SetLineColor(3)
#lyt6.AddEntry(h_plott6,"")
lyt6.AddEntry(h_16t[4],"2016")
lyt6.AddEntry(h_17t[4],"2017")
h_16t[4].Draw("hist")
h_17t[4].Draw("samehist")
lyt6.Draw()
c2.Update()
c2.Print("Study/Result2.pdf")




c3=ROOT.TCanvas("c3","c3",800,800)
c3.Divide(2,3)
c3.cd(1)
ly1=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
ly1.SetTextFont(28)
#h_plot1=ROOT.TH1D("h_plot1","2016 CS angle",100,-1.,1.)
#h_plot1.GetXaxis().SetTitle("cos #theta")
#h_plot1.Draw()
h_16[0].SetTitle("2016 CS angle")
for i in range(len(h_16)):
        h_16[i].SetLineColor(i+2)
        ly1.AddEntry(h_16[i],leg1[i])
	if i>0:
        	h_16[i].Draw("samehist")
	else:
		h_16[i].Draw("hist")
ly1.Draw()

c3.cd(2)
ly2=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
ly2.SetTextFont(28)
#h_plot2=ROOT.TH1D("h_plot2","2017 CS angle",100,-1.,1.)
#h_plot2.GetXaxis().SetTitle("cos #theta")
#h_plot2.Draw()
h_17[0].SetTitle("2017 CS angle")
for i in range(len(h_17)):
        h_17[i].SetLineColor(i+2)
        ly2.AddEntry(h_17[i],leg1[i])
	if i>0:
        	h_17[i].Draw("samehist")
	else:
		h_17[i].Draw("hist")
ly2.Draw()

c3.cd(3)
ly3=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
ly3.SetTextFont(28)
#h_plot3=ROOT.TH1D("h_plot3","Comparsion of the CS between 2016 and 2017 at 2 TeV to 3 TeV",100,-1.,1.)
#h_plot3.GetXaxis().SetTitle("cos #theta")
#h_plot3.Draw()
h_16[1].SetTitle("Comparsion of the CS between 2016 and 2017 at 2 TeV to 3 TeV")
h_16[1].SetLineColor(2)
h_17[1].SetLineColor(3)
ly3.AddEntry(h_16[1],"2016")
ly3.AddEntry(h_17[1],"2017")
h_16[1].Draw("hist")
h_17[1].Draw("samehist")
ly3.Draw()


c3.cd(4)
ly4=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
ly4.SetTextFont(28)
#h_plot4=ROOT.TH1D("h_plot4","Comparsion of the CS between 2016 and 2017 at 3 TeV to 4 TeV",100,-1.,1.)
#h_plot4.GetXaxis().SetTitle("cos #theta")
#h_plot4.Draw()
h_16[2].SetTitle("Comparsion of the CS between 2016 and 2017 at 3 TeV to 4 TeV")
h_16[2].SetLineColor(2)
h_17[2].SetLineColor(3)
ly4.AddEntry(h_16[2],"2016")
ly4.AddEntry(h_17[2],"2017")
h_16[2].Draw("hist")
h_17[2].Draw("samehist")
ly4.Draw()

c3.cd(5)
ly5=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
ly5.SetTextFont(28)
#h_plot5=ROOT.TH1D("h_plot5","Comparsion of the CS between 2016 and 2017 at 4 TeV to 5 TeV",100,-1.,1.)
#h_plot5.GetXaxis().SetTitle("cos #theta")
#h_plot5.Draw()
h_16[3].SetTitle("Comparsion of the CS between 2016 and 2017 at 4 TeV to 5 TeV")
h_16[3].SetLineColor(2)
h_17[3].SetLineColor(3)
ly5.AddEntry(h_16[3],"2016")
ly5.AddEntry(h_17[3],"2017")
h_16[3].Draw("hist")
h_17[3].Draw("samehist")
ly5.Draw()


c3.cd(6)
ly6=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
ly6.SetTextFont(28)
#h_plot6=ROOT.TH1D("h_plot6","Comparsion of the CS between 2016 and 2017 at 5 TeV to 6 TeV",100,-1.,1.)
#h_plot6.GetXaxis().SetTitle("cos #theta")
#h_plot6.Draw()
h_16[4].SetTitle("Comparsion of the CS between 2016 and 2017 at 5 TeV to 6 TeV")
h_16[4].SetLineColor(2)
h_17[4].SetLineColor(3)
ly6.AddEntry(h_16[4],"2016")
ly6.AddEntry(h_17[4],"2017")
h_16[4].Draw("hist")
h_17[4].Draw("samehist")
ly6.Draw()

c3.Print("Study/Result3.pdf")

#masses
#bng1=[i for i in range(200,1000,5)]+[i for i in range(1000,2000,10)]+[i for i in range(2000,3500,20)]
#bng1=numpy.asarray(bng1,dtype=numpy.float64)
#h_Corr1=ROOT.TH1D("h_Corr1", "Pearson Correlation Coefficient",len(bng1)-1,bng1)
#h_Corr1.SetLineColor(2)
#h_Corr2=ROOT.TH1D("h_Corr2", "Pearson Correlation Coefficient",len(bng1)-1,bng1)
#h_Corr2.SetLineColor(3)
#h_Corr3=ROOT.TH1D("h_Corr3", "Pearson Correlation Coefficient",len(bng1)-1,bng1)
#h_Corr3.SetLineColor(4)
#h_Corr4=ROOT.TH1D("h_Corr4", "Pearson Correlation Coefficient",len(bng1)-1,bng1)
#h_Corr4.SetLineColor(5)
#h_Corr4.GetXaxis().SetTitle("mll[GeV]")
#h_Corr4.GetYaxis().SetRangeUser(-0.15,0.15)
#h_CorrT=ROOT.TH1D("h_CorrT", "Pearson Correlation Coefficient",len(bng1)-1,bng1)
#h_CorrT.SetLineColor(2)
#h_Corr=ROOT.TH1D("h_Corr", "Pearson Correlation Coefficient",len(bng1)-1,bng1)
#h_Corr.SetLineColor(3)
#h_Corr.GetXaxis().SetTitle("mll[GeV]")
#h_Corr.GetYaxis().SetRangeUser(-0.02,0.06)
#ly2=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)
#ly3=ROOT.TLegend(0.1, 0.75, 0.4, 0.925)


#for i in range(1,len(bng1)):
#	if bng1[i]<1000:
#		lowerLimit=0.9*bng1[i]
#		upperLimit=1.1*bng1[i]
#	elif bng1[i]<2000:
#		lowerLimit=0.85*bng1[i]
#		upperLimit=1.15*bng1[i]
#	else:
#		lowerLimit=0.8*bng1[i]
#		upperLimit=1.2*bng1[i]
#	obsmass1=mass[(mass<upperLimit) & (mass>lowerLimit) &(numpy.absolute(y)<0.25)]
#	obsCS1=CS[(mass<upperLimit)&(mass>lowerLimit)&(numpy.absolute(y)<0.25)]
#	corr1=numpy.corrcoef(obsmass1,obsCS1)
#	h_Corr1.SetBinContent(i,corr1[0,1])
#	obsmass2=mass[(mass<upperLimit) & (mass>lowerLimit) &(numpy.absolute(y)>0.25)&(numpy.absolute(y)<0.5)]
#        obsCS2=CS[(mass<upperLimit)&(mass>lowerLimit)&(numpy.absolute(y)>0.25)&(numpy.absolute(y)<0.5)]
#        corr2=numpy.corrcoef(obsmass2,obsCS2)
#        h_Corr2.SetBinContent(i,corr2[0,1])
#        obsmass3=mass[(mass<upperLimit) & (mass>lowerLimit) &(numpy.absolute(y)>0.5)&(numpy.absolute(y)<0.75)]
#        obsCS3=CS[(mass<upperLimit)&(mass>lowerLimit)&(numpy.absolute(y)>0.5)&(numpy.absolute(y)<0.75)]
#        corr3=numpy.corrcoef(obsmass3,obsCS3)
#        h_Corr3.SetBinContent(i,corr3[0,1])
#	obsmass4=mass[(mass<upperLimit) & (mass>lowerLimit) &(numpy.absolute(y)>0.75)&(numpy.absolute(y)<1.)]
#        obsCS4=CS[(mass<upperLimit)&(mass>lowerLimit)&(numpy.absolute(y)>0.75)&(numpy.absolute(y)<1.)]
#        corr4=numpy.corrcoef(obsmass4,obsCS4)
#        h_Corr4.SetBinContent(i,corr4[0,1])
#	obsmass=mass[(mass<upperLimit) & (mass>lowerLimit)]
#        obsCST=CS_true[(mass<upperLimit)&(mass>lowerLimit)]
#	obsCS=CS[(mass<upperLimit)&(mass>lowerLimit)]
#        corrT=numpy.corrcoef(obsmass,obsCST)
#        h_CorrT.SetBinContent(i,corrT[0,1])
#	corr=numpy.corrcoef(obsmass,obsCS)
#	h_Corr.SetBinContent(i,corr[0,1])
#plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
#plotPad.Draw()
#plotPad.cd()
#c2=ROOT.TCanvas("c2","c2",800,800)
#c2.Divide(2,1)
#c2.cd(1)
#ly2.AddEntry(h_Corr1,"0<|y|<0.25")
#ly2.AddEntry(h_Corr2,"0.25<|y|<0.5")
#ly2.AddEntry(h_Corr3,"0.5<|y|<0.75")
#ly2.AddEntry(h_Corr4,"0.75<|y|<1")
#h_Corr4.Draw()
#h_Corr3.Draw("same")
#h_Corr2.Draw("same")
#h_Corr1.Draw("same")
#ly2.SetTextFont(32)
#ly2.Draw()
#ROOT.gPad.SetLogx()
#c2.cd(2)
#ly3.AddEntry(h_Corr,"CS")
#ly3.AddEntry(h_CorrT,"True CS")
#h_Corr.Draw()
#h_CorrT.Draw("same")
#ly3.SetTextFont(32)
#ly3.Draw()
#plotPad.SetLogx()
#ROOT.gPad.SetLogx()
#plotPad.Update()
#c2.Print("Study/Result2.pdf")
#plotPad.Update()

#h_CS1000to1200=ROOT.TH1D("h_CS1","CS angle",100,-1.,1.)
#h_CS1000to1200.SetLineColor(2)
#h_CS1000to1200.GetXaxis().SetTitle("cos #theta")
#h_CS3000to3400=ROOT.TH1D("h_CS2","CS angle",100,-1.,1.)
#h_CS3000to3400.SetLineColor(3)

#h_CS1000to1200T=ROOT.TH1D("h_CST1","CS angle truth",100,-1.,1.)
#h_CS1000to1200T.SetLineColor(2)
#h_CS3000to3400T=ROOT.TH1D("h_CST2","CS angle truth",100,-1.,1.)
#h_CS3000to3400T.SetLineColor(3)
#h_CS3000to3400T.GetXaxis().SetTitle("cos #theta")

#h_CS1000to1200r1=ROOT.TH1D("h_CS1r1","CS angle |y|<0.5",100,-1.,1.)
#h_CS1000to1200r1.SetLineColor(2)
#h_CS3000to3400r1=ROOT.TH1D("h_CS2r1","CS angle |y|<0.5",100,-1.,1.)
#h_CS3000to3400r1.SetLineColor(3)
#h_CS3000to3400r1.GetXaxis().SetTitle("cos #theta")

#h_CS1000to1200r2=ROOT.TH1D("h_CS1r2","CS angle",100,-1.,1.)
#h_CS1000to1200r2.SetLineColor(2)
#h_CS3000to3400r2=ROOT.TH1D("h_CS2r2","CS angle 0.5<|y|<1",100,-1.,1.)
#h_CS3000to3400r2.SetLineColor(3)
#h_CS3000to3400r2.GetXaxis().SetTitle("cos #theta")

#h_CS1000to1200Tr=ROOT.TH1D("h_CS1tr","CS angle truth |y|<1",100,-1.,1.)
#h_CS1000to1200Tr.SetLineColor(2)
#h_CS3000to3400Tr=ROOT.TH1D("h_CS2tr","CS angle truth |y|<1",100,-1.,1.)
#h_CS3000to3400Tr.SetLineColor(3)
#h_CS3000to3400Tr.GetXaxis().SetTitle("cos #theta")

#h_CS1000to1200.FillN(len(CS[(mass<1200) & (mass>1000)]),CS[(mass<1200) & (mass>1000)],numpy.ones(len(CS[(mass<1200) & (mass>1000)])))
#h_CS1000to1200.Scale(1./len(CS[(mass<1200) & (mass>1000)]))
#h_CS3000to3400.FillN(len(CS[(mass<3400) & (mass>3000)]),CS[(mass<3400) & (mass>3000)],numpy.ones(len(CS[(mass<3400) & (mass>3000)])))
#h_CS3000to3400.Scale(1./len(CS[(mass<3400) & (mass>3000)]))
#h_CS1000to1200T.FillN(len(CS_true[(mass<1200) & (mass>1000)]),CS_true[(mass<1200) & (mass>1000)],numpy.ones(len(CS_true[(mass<1200) & (mass>1000)])))
#h_CS1000to1200T.Scale(1./len(CS_true[(mass<1200) & (mass>1000)]))
#h_CS3000to3400T.FillN(len(CS_true[(mass<3400) & (mass>3000)]),CS_true[(mass<3400) & (mass>3000)],numpy.ones(len(CS_true[(mass<3400) & (mass>3000)])))
#h_CS3000to3400T.Scale(1./len(CS_true[(mass<3400) & (mass>3000)]))
#h_CS1000to1200r1.FillN(len(CS[(mass<1200) & (mass>1000) & (numpy.absolute(y)<0.5)]),CS[(mass<1200) & (mass>1000)&(numpy.absolute(y)<0.5)],numpy.ones(len(CS[(mass<1200) & (mass>1000) & (numpy.absolute(y)<0.5)])))
#h_CS1000to1200r1.Scale(1./len(CS[(mass<1200) & (mass>1000) & (numpy.absolute(y)<0.5)]))
#h_CS3000to3400r1.FillN(len(CS[(mass<3400) & (mass>3000) & (numpy.absolute(y)<0.5)]),CS[(mass<3400) & (mass>3000)& (numpy.absolute(y)<0.5)],numpy.ones(len(CS[(mass<3400) & (mass>3000)& (numpy.absolute(y)<0.5)])))
#h_CS3000to3400r1.Scale(1./len(CS[(mass<3400) & (mass>3000)& (numpy.absolute(y)<0.5)]))
#h_CS1000to1200r2.FillN(len(CS[(mass<1200) & (mass>1000) & (numpy.absolute(y)>0.5)&(numpy.absolute(y)<1)]),CS[(mass<1200) & (mass>1000)&(numpy.absolute(y)>0.5)&(numpy.absolute(y)<1)],numpy.ones(len(CS[(mass<1200) & (mass>1000) & (numpy.absolute(y)>0.5)&(numpy.absolute(y)<1)])))
#h_CS1000to1200r2.Scale(1./len(CS[(mass<1200) & (mass>1000) & (numpy.absolute(y)>0.5)&(numpy.absolute(y)<1)]))
#h_CS3000to3400r2.FillN(len(CS[(mass<3400) & (mass>3000) & (numpy.absolute(y)>0.5)& (numpy.absolute(y)<1)]),CS[(mass<3400) & (mass>3000)& (numpy.absolute(y)>0.5) & (numpy.absolute(y)<1)],numpy.ones(len(CS[(mass<3400) & (mass>3000)& (numpy.absolute(y)>0.5) & (numpy.absolute(y)<1)])))
#h_CS3000to3400r2.Scale(1./len(CS[(mass<3400) & (mass>3000)& (numpy.absolute(y)>0.5) & (numpy.absolute(y)<1)]))
#h_CS1000to1200Tr.FillN(len(CS_true[(mass<1200) & (mass>1000) &(numpy.absolute(y)<1)]),CS_true[(mass<1200) & (mass>1000)&(numpy.absolute(y)<1)],numpy.ones(len(CS_true[(mass<1200) & (mass>1000) &(numpy.absolute(y)<1)])))
#h_CS1000to1200Tr.Scale(1./len(CS_true[(mass<1200) & (mass>1000) &(numpy.absolute(y)<1)]))
#h_CS3000to3400Tr.FillN(len(CS_true[(mass<3400) & (mass>3000) & (numpy.absolute(y)<1)]),CS_true[(mass<3400) & (mass>3000)& (numpy.absolute(y)<1)],numpy.ones(len(CS_true[(mass<3400) & (mass>3000)& (numpy.absolute(y)<1)])))
#h_CS3000to3400Tr.Scale(1./len(CS_true[(mass<3400) & (mass>3000)&  (numpy.absolute(y)<1)]))

#c3=ROOT.TCanvas("c3","c3",800,800)
#c3.Divide(2,3)
#c3.cd(1)
#h_CS3000to3400.Draw()
#h_CS1000to1200.Draw("same")
#c3.cd(2)
#h_CS3000to3400T.Draw()
#h_CS1000to1200T.Draw("same")
#c3.cd(3)
#h_CS3000to3400r1.Draw()
#h_CS1000to1200r1.Draw("same")
#c3.cd(4)
#h_CS3000to3400r2.Draw()
#h_CS1000to1200r2.Draw("same")
#c3.cd(5)
#h_CS3000to3400Tr.Draw()
#h_CS1000to1200Tr.Draw("same")

#c3.Print("Study/Result3.pdf")


