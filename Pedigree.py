######################################################################################################
# This class interfaces with the SQL database to allow for searching through the pedigree information
# the functions are
# __init__(SQLdb): connects to the SQL database and creates the pedigree tree.
# getOriginators(node,depth):  Get "start" of the pedigree for a specific node, i.e. parents/grandparents 
#     Searches up the pedigree from a set node and retrieves all nodes at a specific depth.  
#			If the search research a parent at a lower depth which has no parents, then this parent node is returned as an originator.
#     If the search cannot progree upto the depth, all originators at the maximum possible depth are returned as well as the maximum possible depth
# getProgeny(parents,depth): Get all the children from a set if parents.
#     Searches down the pedigree from a set of parent nodes and returns all children
# getExtendedFamily(node,depth): Gets all family members related to node, to a specified depth
#			First calls getOriginators to get the "start" of the pedigree at the specified depth
#     Then calls getProgeny to get all members within that pedigree starting from the originators.
######################################################################################################

import sqlite3
import collections

class PedigreeNode:
	index = -1
	ID = None
	mother = None
	father = None
	sex = None
	pedigree = None
	children = []
	
	def __init__(self,INDEX,ID,MA,FA,SEX,PEDIGREE):
		self.index =	INDEX
		self.ID = ID
		self.mother = MA
		self.father = FA
		self.sex = SEX
		self.pedigree = PEDIGREE
		self.children = []
	
	def addChild(self,child):
		self.children.append(child)
		self.children = list(set(self.children))


class Pedigree:
	PED = collections.defaultdict(list)
	
	def __init__(self, SQLdb):
		self.PED = collections.defaultdict(list)

		# Make all the nodes
		queryres = SQLdb.execute("select id, name, maternal, paternal, male, family from ind;")
		for row in queryres:
			self[str(row[1])] = PedigreeNode( row[0],str(row[1]),str(row[2]),str(row[3]),row[4],row[5] )		

		# add the mother/father edges
		for node in self.PED:
			if self[node].mother != "0":
				self[ str(self[node].mother) ].addChild(self[node].ID)
			if self[node].father != "0":
				self[ str(self[node].father) ].addChild(self[node].ID)
				

	def __setitem__(self,key,value):
		self.PED[key] = value

	def keys(self):
		return self.PED.keys()

	def __getitem__(self,node):
		return self.PED[node]

	def getParents(self,node):
		pars = []
		if self[node].mother != "0":
			pars.append(self.PED[node].mother)
		if self[node].father != "0":
			pars.append(self.PED[node].father)		
		return pars

	def getChildren(self,node):
		return self[node].children

	# search back through the generations to a specified depth
	# to find the originator couples.  0 = parents, 1 = grandparents etc
	def getOriginators(self,node,depth):
		d = 0
		originators = []
		z = self.getOriginatorsProxy(node,originators,d,depth)
		return {'originators': originators,'maxdepth': d}
		
	def getOriginatorsProxy(self,node,originators,depth,maxdepth):
 		zpars = self.getParents(node)
		if len(zpars) == 0:
			return False	
	
		if len(zpars) > 0:
			depth = depth + 1
			# get the parents at the maximum depth
			if depth == maxdepth: 
				originators += zpars
				return True
				
			# get any non-maximum depth parents
			for x in zpars:
				xpars = self.getParents(x)		
				if len(xpars) == 0:
					originators += [x]

		for x in zpars:
			self.getOriginatorsProxy(x,originators,depth,maxdepth)
	
	# get all the children from a set of parents upto a specified depth
	# 0 = children, 1 = grandchildren etc..
	def getProgeny(self,parents,depth):
		d = 0
		progeny = []
		if isinstance(parents,str):
			z = self.getProgenyProxy([parents],progeny,d,depth)
		else: 
			z = self.getProgenyProxy(parents,progeny,d,depth) 
		
		zret = list(set(progeny))
		return {"progeny": zret,"maxdepth": d}

	def getProgenyProxy(self,parents,progeny,depth,maxdepth):		
		# get all kids
		allkids = []
		for x in parents:
			allkids += self.getChildren(x)
		
		allkids = list(set(allkids))

		# if there are no kids
		if len(allkids) == 0:
			return False
		
		# otherwise add the kids
		progeny += allkids
		if depth == maxdepth:
			return True			

		# keep going if not at maxdepth
		depth = depth + 1	
		for x in allkids:
			self.getProgenyProxy([x],progeny,depth,maxdepth)

	# gets all the members of an extended family
	# 0 = parents
 	# 1 = grand parents ... etc	
 	def getExtendedFamily(self,node,depth):
		zoriginators = self.getOriginators(node,depth)
		zfamily = self.getProgeny(zoriginators["originators"],zoriginators["maxdepth"])
		zret = list(set(zfamily["progeny"] + zoriginators["originators"]))
		return {"members": zret,"maxdepth":zoriginators["maxdepth"]}
	





