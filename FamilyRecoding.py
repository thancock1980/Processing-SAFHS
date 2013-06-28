#######################################################################################################
# This code performs the compares the phenotype and genotypes of individuals in the same family and
# outputs a plink like file which can be used for GWAS or GWIS.
# the key functions are:
# __init__(self,SQLdb,plinkstem,phenotype,recodingMethod,depth): which takes a pointer to an SQL database,
#          the location and name (plinkstem) of the returned plink database 
#          (e.g. /test/test will write the /test/test.bed, /test/test.bim, /test/test.fam)
#          a string with the phenotype label, the recoding method (1,2,3), and the depth which to consider the family.
#
# There are 3 recoding methods which consider the comparision of family members
# 1) label change of phenotype between family members
# 2) if one person has a disease then label as affected
# 3) if progeny has a SNP which causes or protects (label by progeny)
# 
# The genotypes are recoded in the following way: (sticking with plink conventions)
# 1) if one individuals SNP is missing, the interaction is labelled as missing
# 2) if both individuals SNPs are identical, then label as homozygous.
# 3) if one allele is different between individuals, then label as heterozygous
# 4) if both alleles are different between individuals, then label as non-reference homozygous
#
# recode() writes out the database.  Each individual is treated as a parent, and using
# the functions in Pedigree.py, all progeny are extracted, and the comparison is perform with respect
# to parent and all progeny.

import sqlite3

class FamilyRecoding:	
	plinkstem = "out"
	SQLdb = None
	phenotype = ""
	recodingMethod = None
	PED = None
	comparisons = None
	depth = None
	
	def __init__(self,SQLdb,plinkstem,phenotype,recodingMethod,depth):
		self.plinkstem = plinkstem
		self.SQLdb = SQLdb
		self.phenotype = phenotype
		self.recodingMethod = recodingMethod
		self.PED = Pedigree(SQLdb)
		self.depth = depth
		comparisons = []

	def recode(self):
		self.recodeFAM()
		self.recodeBIM()
		self.recodeBED()

	# scheme 1: label change of phenotype between family members
	def recodeOne(self,p1,p2):
		if p1 == p2: return 1
		else: return 2

	# scheme 2: if one person has a disease then label as affected
	def recodeTwo(self,p1,p2):
		if p1 == 0 or p2 == 0: 	return 1
		else: return 2
	
	# scheme 3: progeny has a SNP which causes or protects (label by progeny)
	def recodeTwo(self,p1,p2):
		if p2 == 0: return 1
		else: return 2


	def recodeFAM(self):
		query = "select ind.id, family, name, paternal, maternal, male, " + self.phenotype + " " + \
				    "from ind join phenotype on ind.id = phenotype.id " + \
				    "where ind.hasphenotype = '1' and ind.hasgenotype = '1' order by ind.id;"
		indiv = self.SQLdb.execute( query ).fetchall()
		
		print "Writing FAM file to: " + self.plinkstem + ".fam"
		famfile = open(self.plinkstem + ".fam","w")
		
		depth = 2
		self.comparisons = []
		for sample in indiv:
			zfam = self.PED.getProgeny([sample[2]],depth)["progeny"]

			if len(zfam) > 0:	
				for fid in zfam:
					# get the comparisions information
					query = "select ind.id, name, " + self.phenotype + " " + \
				    			"from ind join phenotype on ind.id = phenotype.id where ind.name = ?"
					compsample = self.SQLdb.execute(query, (fid,) ).fetchone()
					
					# Update the comparsions list so that the BED file is populated correctly
					self.comparisons.append([sample[0],compsample[0]])

					# recode the label
					pheno = 1
					if self.recodingMethod == 1: pheno = self.recodeOne(sample[6],compsample[2])
					if self.recodingMethod == 2: pheno = self.recodeTwo(sample[6],compsample[2])
					if self.recodingMethod == 3: pheno = self.recodeThree(sample[6],compsample[2])

					# family, name, paternal, maternal, male, phenotype
					famline = "1" + "\t" + \
									  ":".join([sample[2],compsample[1]]) + "\t" + \
									  "0" + "\t" + \
									  "0" + "\t" + \
									  "1" + "\t" + \
									  str(pheno)
					famfile.write(famline + "\n")
		print "Done"

	def recodeBIM(self):
		print "Writing BIM file to: " + self.plinkstem + ".bim"
		rows = self.SQLdb.execute("select chr, name, position, min, maj from snp order by id").fetchall()
	
		bimfile = open(self.plinkstem + ".bim","w")
		bimfile.write("\n".join("%s\t%s\t0\t%d\t%s\t%s" % x for x in rows) + "\n")
		bimfile.close()
		print "Done"

	def recodeBED(self):
		print "Writing BED file to: " + self.plinkstem + ".bed"

		## indiv ped packing
		# get the number of SNPs, individuals and all the major/minor SNP associations
		nsnps = self.SQLdb.execute("select count(*) from snp;").fetchone()[0]
		majmin = self.SQLdb.execute("select maj, min from snp;").fetchall()

		# convert the sample ids to a list
		nid = len(self.comparisons)

		# get the size of the buffer to write
		bsize = nid / 4
		if (nid % 4) > 0: bsize = bsize + 1

		# plink mapping scheme
		# Homozygote "1"/"1": 3
    # Heterozygote: 2
    # Homozygote "2"/"2": 0
    # Missing genotype: 1
		# bitmasks 240 = 11110000 and 15 = 00001111 
		def plinkMAP(snp1, snp2):
			if snp1 == 0 or snp2 == 0: 
				return 1 # missing values (plink stores as 2)
			elif snp1 == snp2: 
				return 3  # both individual have the same alleles for both major and minor (store as plink homozyogous: 0)
			elif ((snp1 & 240) == (snp2 & 240)) or ((snp1 & 15) == (snp2 & 15)): 		
				return 2 # one allele change between individuals (store as plink heterozygous: 1) (240 = 11110000, 15 = 00001111) 
			else: 
				return 0 # complete change of both alleles(store as plink non-reference homozygous: 3)
		
		# open the file
		bedfile = open(self.plinkstem + ".bed","wb")
		
		# define and write the magic number
		magic = bytearray(3)
		magic[0] = 108
		magic[1] = 27
		magic[2] = 1
		bedfile.write(magic)
		
		# write out in SNP dominant form
		shift = [0,2,4,6]
		gtypes = ["G","A","T","C","0"]
		for snpid in range(0,nsnps):
			buf = bytearray(bsize)
			# get one SNP
			snps = self.SQLdb.execute("select genotype from snp where id = ?;",(snpid,)).fetchone()[0]
			j = 0
			idx = 0
			for zid in range(0,len(self.comparisons)):
				mapped = plinkMAP( ord(snps[self.comparisons[zid][0]]), ord(snps[self.comparisons[zid][1]]) )		
				buf[j] = buf[j] | (mapped << shift[idx])
				if idx == 3: 
					j = j + 1 # go to the next byte
					idx = 0
				else: idx = idx + 1		
			
			#write the buffer
			bedfile.write(buf)		
	
		# all done!
		bedfile.close()
		print "Done"
			
			



