##
## This python script converts the genos files which I
## received from Baker IDI to a ped file for plink
##
import tarfile
import re
import sqlite3
import collections
from bitarray import bitarray
import os

chromosomes = map(str,range(1,23))

# pedigree initialization
indiv = {"ID": [],"SEX": [], "MO": [], "FA": [], "PEDIGREE": [],"HASOUTCOME": [],"HASGENOTYPE": []}
pedfile = open("pedigree/SAFHSpedigree.csv","r")
#pedfile = open("/home/thancock/Data/test/indiv.txt")

# outcomes initialization
VARS = {'INDEX': [], 'CLASS': [], 'SHORTNAME': [], 'NAME': [],'TYPE': []}
outfile = open("/home/thancock/Data/SAFHS/lipidcsv/SAFHS_LIPID_NICTA.csv","r")
#outfile = open("/home/thancock/Data/test/pheno.txt")

# probe location initialization
mapfile = {'SNPID': [], 'CHROMOSOME': [], 'POSITION': [], \
				 'FMaj': [], 'FMin': [], 'SNPBLOB': [], \
				 'Gmaj': [],'Amaj': [],'Tmaj': [], 'Cmaj': [],'0maj': [], 'Imaj': [], 'Dmaj': [], \
				 'Gmin': [],'Amin': [],'Tmin': [], 'Cmin': [],'0min': [], 'Imin': [], 'Dmin': []}
mapbasedir = "/home/thancock/Data/SAFHS_23052013/chr"
#mapbasedir = "/home/thancock/Data/test/chr"

# genotype initialization 
genobasedir = "/home/thancock/Data/SAFHS_23052013/chr"
#genobasedir = "/home/thancock/Data/test/chr"

###
## Process the pedigree 
###
lines = pedfile.readlines()
lines.remove(lines[0])
for line in lines:
	fields = line.strip().split(",")
	fields = ["0" if x == "" else x for x in fields]
	indiv["ID"].append(fields[0])
	indiv["SEX"].append(fields[3])
	indiv["MO"].append(fields[2])
	indiv["FA"].append(fields[1])
	indiv["PEDIGREE"].append(fields[5])
	indiv["HASOUTCOME"].append(False)
	indiv["HASGENOTYPE"].append(False)

maxpedigree = max(map(int,indiv["PEDIGREE"]))

###
## Process the outcome information 
###
# create new phenotype names because the actual ones have invalid characters
lines = outfile.readlines()
VARS["SHORTNAME"] = ["X" + str(x) for x in range(0,len(lines[0].strip().split(";")))]
VARS["NAME"] = lines[0].strip().split(";")
VARS["CLASS"] = lines[1].strip().split(";")
VARS["INDEX"] = range(1,len(VARS["SHORTNAME"]))
VARS["TYPE"] = lines[2].strip().split(";")

outcomes = collections.defaultdict(list)
for vname in VARS["SHORTNAME"]:
	outcomes[vname] = ["" for _ in range(len(indiv["ID"]))]
	
outcomes["X0"] = indiv["ID"][:]
missingindiv = []
for i in range(3,len(lines)):
	line = lines[i].strip().split(";")
	if line[0] in indiv["ID"]:
		zid = indiv["ID"].index(line[0])
		for j in range(0,len(VARS["SHORTNAME"])):
			outcomes[ VARS["SHORTNAME"][j] ][zid] = line[j]
		indiv["HASOUTCOME"][zid] = True
	else:
		print "Couldn't find " + line[0] + " in pedigree"
		print "... appending information and patching the pedigree"
		missingindiv.append(line[0])
		for j in range(0,len(VARS["SHORTNAME"])): #patch the outcomes table
			outcomes[ VARS["SHORTNAME"][j] ].append( line[j] ) 
		outid = outcomes["X0"].index(line[0])
		indiv["ID"].append(line[0])
		indiv["SEX"].append( outcomes[VARS["SHORTNAME"][VARS["NAME"].index("SEX")]][outid] )
		indiv["MO"].append("")
		indiv["FA"].append("")
		indiv["PEDIGREE"].append(maxpedigree + 1 + i)
		indiv["HASOUTCOME"].append(True)
		indiv["HASGENOTYPE"].append(False)

###
## Process the SNP bim information
###
nindiv = len(indiv["ID"])
for chromosomeid in chromosomes:
	ztar = tarfile.open(mapbasedir + chromosomeid + ".maps")
	mapsnps = 0
	ztarmembers = ztar.getmembers()
	ztarnames = ztar.getnames()
	for zid in range(0,len(ztarmembers)):
		fname = "map." + chromosomeid + "_" + str(zid+1)
		print "Processing " + fname
		tarpart = ztarmembers[ztarnames.index(fname)]
		mappart = ztar.extractfile(tarpart).readlines()
		mappart.remove(mappart[0]) # remove header row
		mapsnps += len(mappart)
		print "  Samples in file: " + str(len(mappart)) + " total samples: " + str(mapsnps)
		for zid in range(0,len(mappart)):
			maprow = re.split(" +",mappart[zid].strip())
			mapfile['CHROMOSOME'].append(chromosomeid)
			mapfile['SNPID'].append(maprow[0])
			mapfile['POSITION'].append(maprow[1])
			mapfile['FMaj'].append(0)
			mapfile['FMin'].append(0)
			mapfile['SNPBLOB'].append(bytearray(nindiv))
			mapfile['Gmaj'].append(0) # all this stuff is needed to get the
			mapfile['Tmaj'].append(0) # consensus major and minor allele for
			mapfile['Amaj'].append(0) # each sample .... yay
			mapfile['Cmaj'].append(0)
			mapfile['0maj'].append(0)
			mapfile['Imaj'].append(0)
			mapfile['Dmaj'].append(0)
			mapfile['Gmin'].append(0)
			mapfile['Tmin'].append(0)
			mapfile['Amin'].append(0)
			mapfile['Cmin'].append(0)
			mapfile['0min'].append(0)
			mapfile['Imin'].append(0)
			mapfile['Dmin'].append(0)
	ztar.close()

###
## Process the genotype information
###
def multi_delete(list_, delindex):
    indexes = sorted(delindex, reverse=True)
    for index in indexes:
        del list_[index]
    return list_

# 1 byte per SNP
# 4 bits for the major and 4 bits for the minor allele
# missing is 0
def snptoint(x):	
	ret = 0	
	if x == "": return ret                       # missing values are mapped to 0
	elif x[0] == "G": ret = 1 << 4
	elif x[0] == "A": ret = 2 << 4
	elif x[0] == "T": ret = 3 << 4
	elif x[0] == "C": ret = 4 << 4
	elif x[0] == "I": ret = 5 << 4
	elif x[0] == "D": ret = 6 << 4
	
	if x[2] == "G": ret = ret | 1
	elif x[2] == "A": ret = ret | 2
	elif x[2] == "T":	ret = ret | 3
	elif x[2] == "C":	ret = ret | 4
	elif x[2] == "I": ret = ret | 5
	elif x[2] == "D": ret = ret | 6
	
	return ret	

probeids = []
missingprobes = []
pedsnps = 0
alleles = ["G","T","A","C","I","D","0"] 		
for chromosomeid in chromosomes:
	ztar = tarfile.open(genobasedir + chromosomeid + ".genos.gz","r:gz")
	ztarmembers = ztar.getmembers()
	ztarnames = ztar.getnames()
	badsamples = collections.defaultdict(list)
	for zid in range(0,len(ztarmembers)):
		fname = "geno." + chromosomeid + "_" + str(zid+1)
		print "Processing " + fname
		tarpart = ztarmembers[ztarnames.index(fname)]
		pedpart = ztar.extractfile(tarpart).readlines()
		pedpart = [line.strip().split(",") for line in pedpart]
		fileprobes = pedpart[0][1:len(pedpart[0])] 
		
		# delete probes we don't have mapped
		pdiff = list(set(fileprobes) - set(mapfile["SNPID"])) # find any probes in the genos file not in the map files
		nonmatchedprobesidx = [] # index of the ID column (always remove this one)
		if len(pdiff) > 0:  # get rid of probes that are not in the map file - i dont have their positions...
			missingprobes += pdiff
			print "  File: " + fname + " has the following unmatched probes: " + str(pdiff)
			print "  ... deleting them ..."
			nonmatchedprobesidx += [fileprobes.index(zi) for zi in pdiff] # get the index of the unmatched probes
			fileprobes = multi_delete(fileprobes,nonmatchedprobesidx) # delete the probes not in the map file.
		
		probeids += fileprobes # append to the list of all probes
		pedsnps += len(fileprobes) 
		print "  Samples: " + str(len(pedpart) -1) + " SNPs: " + str( len(fileprobes) ) + " Cumulative SNPs: " + str(pedsnps)
		pedpart.remove(pedpart[0]) # remove the sample and probe id row
		probeidx = [mapfile["SNPID"].index(x) for x in fileprobes] # get the position of the probes in the map file
		rowids = []
		
		# count the major and minor allele for each SNP over all samples
		# and remove probes that are not in the mapfile.
		for sample in pedpart:
			idx = indiv['ID'].index(sample[0]) # find the id of the sample in the databse
			indiv["HASGENOTYPE"][idx] = True			
			rowids.append(idx)			
			sample.remove(sample[0]) # remove the ID column
			
			# remove non-matched probes
			if (len(missingprobes) > 0): sample = multi_delete(sample,nonmatchedprobesidx) 
					
			# count the major and minor allele for each SNP over all samples
			for x in zip(sample,probeidx):  
				if x[0] == "": 
					mapfile["0maj"][ x[1] ] += 1
					mapfile["0min"][ x[1] ] += 1
				else:
					mapfile[x[0][0] + "maj"][ x[1] ] += 1
					mapfile[x[0][2] + "min"][ x[1] ] += 1
					
		# create the binary representation of the genotype data 
		# (assume probes to not duplicate across files and maintaining the order of individuals)
		for i in range(0,len(probeidx)):
			# get the major and minor allele for each SNP 
			pidx = probeidx[i]
			maxmaj = [mapfile[x][pidx] for x in ["Gmaj","Tmaj","Amaj","Cmaj","Imaj","Dmaj","0maj"]]
			maxallele = alleles[maxmaj.index(max(maxmaj))]
			maxmin = [mapfile[x][pidx] for x in ["Gmin","Tmin","Amin","Cmin","Imin","Dmin","0min"]]
			minallele = alleles[maxmin.index(max(maxmin))]
			mapfile["FMaj"][pidx] = maxallele
			mapfile["FMin"][pidx] = minallele
			
			# at the moment there is 1 byte per SNP, (this will be packed to 4 bytes per SNP later).
			# pidx = row position of the id in the mapfile dictionary
			# rowids[k] = the position of the person in the SNP bytearray.
			for k in range(0,len(rowids)):
				mapfile["SNPBLOB"][pidx][ rowids[k] ] = snptoint(pedpart[k][i]) 
	ztar.close()

################################################################################################
################################################################################################
## Begin the database construction
################################################################################################
################################################################################################

genotable = "create table snp(id integer primary key," + \
                "name text," + \
                "chr integer," + \
                "position integer," + \
                "min text," + \
                "maj text," + \
                "bit_stride integer," + \
                "genotype blob)" 
indivtable =  "create table ind(id integer primary key," + \
                "name text," + \
                "family text," + \
                "paternal text," + \
                "maternal text," + \
                "male bool," + \
					 	"hasphenotype bool," + \
					 	"hasgenotype bool)"
phenomapping = "create table phenotypeMap(id text primary key," + \
						"mappedname text," + \
					   "fullname text," + \
					   "type text," + \
					   "class text)" 
freestr = ""
for i in range(0,len(VARS["SHORTNAME"])):
	freestr += "," + VARS["SHORTNAME"][i] + " "
	if VARS["TYPE"][i] == "STRING": freestr += "text"
	if VARS["TYPE"][i] == "BIN": freestr += "bool"
	if VARS["TYPE"][i] == "REAL": freestr += "real"
	if VARS["TYPE"][i] == "CAT": freestr += "text"

phenotable = "create table phenotype(id integer primary key"  + freestr + ")"

os.system("rm SQLdb/safhs.db")
conn = sqlite3.connect("SQLdb/safhs.db")

# create and populate the individuals table
conn.execute(indivtable)
tablerows = zip( range(0,len(indiv["ID"])), \
	   						 indiv["ID"], \
		 						 indiv["PEDIGREE"], \
	               indiv["FA"], \
		             indiv["MO"], \
		             map(lambda x: x == "M", indiv["SEX"]), \
							   indiv["HASOUTCOME"], \
						     indiv["HASGENOTYPE"])
conn.executemany("insert into ind values(?,?,?,?,?,?,?,?)", tablerows)
conn.commit()

# create and populate the phenotype and phenotype mapping table
rows = zip(range(0,len(VARS["SHORTNAME"])),VARS["SHORTNAME"],VARS["NAME"],VARS["TYPE"],VARS["CLASS"])
conn.execute(phenomapping)
conn.executemany("insert into phenotypeMap values(?,?,?,?,?)",rows)

rows = [() for _ in range(len(indiv["ID"]))]
for i in range(0, len(indiv["ID"]) ):
	rows[i] += (i,)
	for vname in VARS["SHORTNAME"]:
		rows[i] = rows[i] + (outcomes[vname][i],)

conn.execute(phenotable)
conn.executemany("insert into phenotype values(" + ",".join(["?"]*(len(rows[0]))) + ")", rows)
conn.commit()

# create and populate the genotype files

for i in range(0,len(mapfile["SNPBLOB"])):
	mapfile["SNPBLOB"][i] = buffer(mapfile["SNPBLOB"][i])
			 
rows = zip(range(0,len(mapfile["SNPID"])), \
				mapfile["SNPID"], \
				mapfile["CHROMOSOME"], \
				mapfile["POSITION"], \
				mapfile["FMin"], \
				mapfile["FMaj"], \
				[4]*len(mapfile["SNPID"]), \
				mapfile["SNPBLOB"] )

conn.execute(genotable)
conn.executemany("insert into snp values(?,?,?,?,?,?,?,?)", rows)
conn.commit()
conn.close()






