#######################################################################################################
# This code interfaces with the SQL database to convert it into plink format
# The key function is:
# 	__init__(SQLdb,plinkstem,phenotype): which takes a SQLdb object (already connected); the location and name (plinkstem) of the
#   returned plink database (e.g. /test/test will write the /test/test.bed, /test/test.bim, /test/test.fam) and the phenotype name
#   as string
#  writePlink(): writes out the plink database in the correct format.
#
# I could probably do some more checking here, but i will leave that to the future.
#######################################################################################################

import sqlite3

class SqliteToPlink:
	plinkstem = "out"
	SQLdb = None
	phenotype = ""

	def __init__(self,SQLdb,plinkstem,phenotype):
		self.plinkstem = plinkstem
		self.SQLdb = SQLdb
		self.phenotype = phenotype

	def writePlink(self):
		self.writeFam()		
		self.writeBim()
		self.writeBed()
		

	def writeBim(self):
		print "Writing BIM file to: " + self.plinkstem + ".bim"
		rows = self.SQLdb.execute("select chr, name, position, min, maj from snp order by id").fetchall()
		
		bimfile = open(self.plinkstem + ".bim","w")
		bimfile.write("\n".join("%s\t%s\t0\t%d\t%s\t%s" % x for x in rows) + "\n")
		bimfile.close()
		print "Done"

	def writeFam(self):
		print "Writing FAM file to: " + self.plinkstem + ".fam"

		query = "select family, name, paternal, maternal, male, " + self.phenotype \
						+ " from ind join phenotype on ind.id = phenotype.id where ind.hasphenotype = '1' and ind.hasgenotype = '1' order by ind.id ;"

		rows = self.SQLdb.execute(query).fetchall()

		famfile = open(self.plinkstem + ".fam","w")
		# strictly limit to binary phenotypes
		for row in rows:
			sex = 2 - row[4]	
			pheno = 1
			if row[5] == 1: 
				pheno = 2
			if row[5] == -1:
				pheno = 0 # missing phenotype value.
	
			famfile.write("\t".join([row[0],row[1],row[2],row[3],str(sex),str(pheno)]) + "\n")

		famfile.close()

		print "Done"

	def writeBed(self):
		print "Writing BED file to: " + self.plinkstem + ".bed"

		## indiv ped packing
		# get the number of SNPs, individuals and all the major/minor SNP associations
		nsnps = self.SQLdb.execute("select count(*) from snp;").fetchone()[0]
		ids =  self.SQLdb.execute("select id from ind where hasphenotype = '1' and hasgenotype = '1';").fetchall()
		majmin = self.SQLdb.execute("select maj, min from snp;").fetchall()

		# convert the sample ids to a list
		nid = len(ids)

		# get the size of the buffer to write
		bsize = nid / 4
		if (nid % 4) > 0: bsize = bsize + 1

		# plink mapping scheme
		# Homozygote "1"/"1": 3
    # Heterozygote: 2
    # Homozygote "2"/"2": 0
    # Missing genotype: 1
		# bitmasks 240 = 11110000 and 15 = 00001111 
		def plinkMAP(x,gmaj,gmin):
			if ((x & 240) == 0) or ((x & 15) == 0): 
				return 1 # missing values (plink stores as 1)
			elif x == ( gmaj << 4 | gmaj ): 
				return 3 # reference homozygote (plink stores as 3)
			elif x == ( gmin << 4 | gmin ): 
				return 0 # non-reference homozygote (plink stores as 0)
			else:
				return 2 # heterzygote (plink stores as 2)

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
		alleles = ["0","G","A","T","C","I","D"]
		for snpid in range(0,nsnps):
			buf = bytearray(bsize)
			# get one SNP
			snps = self.SQLdb.execute("select genotype from snp where id = ?;",(snpid,)).fetchone()[0]
			gmaj = alleles.index(majmin[snpid][0])
			gmin = alleles.index(majmin[snpid][1])
			j = 0
			idx = 0
			for zid in ids:
				mapped = plinkMAP( ord(snps[zid[0]]),gmaj,gmin)		
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

		

