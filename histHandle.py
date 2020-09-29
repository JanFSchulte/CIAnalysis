import ROOT
from copy import deepcopy


class histHandle(object):

	def __init__(self,name,layers):
		self.name=name
		self.layers=layers
		self.nodes={}
		self.__tree={}
		ROOT.TH1.AddDirectory(0)
		ROOT.TH2.AddDirectory(0)
		for layer in layers:
			self.nodes[layer]=set([])

	def __addLeaf(self,histName,nodes,h,tree="none"):
		if tree=="none":
			tree=self.__tree
		if len(nodes)>1:
			if nodes[0] not in tree.keys():
				#print nodes[0]
				tree[nodes[0]]={}
			self.__addLeaf(histName,nodes[1:],h,tree[nodes[0]])	
        	else:
			if nodes[0] in tree.keys():
				#print nodes[0]
				if "hist" in tree[nodes[0]].keys():
					tree[nodes[0]]["hist"].Add(h.Clone())
				else:
					tree[nodes[0]]["hist"]=h.Clone(histName)
			else:
                		tree[nodes[0]]={"hist":h.Clone(histName)}

	def __loopOver(self,node,func,selection=[]):
		flag=True
		for key in node.keys():
			#print key
			if key =="hist":
				func(node[key])	
			elif key in selection:
				self.__loopOver(node[key],func,selection)
				flag=False
		if flag:
			for key in node.keys():
				#print key
				if key !="hist":
					self.__loopOver(node[key],func,selection)

	def loadHist(self,filename,histName,leafPath):
		f=ROOT.TFile.Open(filename)
		hist=f.Get(histName)
		layerNum=len(leafPath)
		for key in leafPath.keys():
			self.nodes[key].add(leafPath[key])
		subLayers=self.layers[0:layerNum]
		nodes=[]
		for layer in subLayers:
			nodes.append(leafPath[layer])
		s="_"
                histName=s.join(nodes)
		#print histName
		self.__addLeaf(histName,nodes,histm)
		f.Close()
	
	def addHist(self,hist,leafPath):
                layerNum=len(leafPath)
                for key in leafPath.keys():
                        self.nodes[key].add(leafPath[key])
                subLayers=self.layers[0:layerNum]
                nodes=[]
                for layer in subLayers:
                        nodes.append(leafPath[layer])
		s="_"
		histName=s.join(nodes)
		#print histName
                self.__addLeaf(histName,nodes,hist)
	
	def loadHists(self,filename,histList):
		for dual in histList:
			self.loadHist(filename,dual[0],dual[1])
		
	def getHist(self,leafPath):
		tree=self.__tree.copy()
		nodes=[]
		num=len(leafPath)
		layers=self.layers[0:num]
		for layer in layers:
			nodes.append(leafPath[layer])
		#print nodes
		for node in nodes:
			#print tree
			tree=tree[node].copy()
		h=tree["hist"].Clone()
		return h
	def process(self,func,leafPath): 	
		tree=self.__tree.copy()
		num=len(leafPath)
		layers=self.layers[0:num]
                nodes=[]
                for layer in layers:
                        nodes.append(leafPath[layer])
                for node in nodes:
                        tree=tree[node].copy()
		func(tree["hist"])
	def processes(self,func,selection=""):
		tree=self.__tree.copy()
                if selection=="":
			self.__loopOver(tree,func)
		else:
			self.__loopOver(tree,func,selection)
	def saveHists(self,filename,selection=""):
		f=ROOT.TFile(filename,"RECREATE")
		def saveHist(hist):
			hist.Write()
		self.processes(saveHist,selection)	
