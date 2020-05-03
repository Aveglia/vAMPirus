#!/usr/local/python/bin/python

#    Copyright 2002,2003,2004,2005 Rasmus Wernersson, Technical University of Denmark
#
#    This file is part of VirtualRibosome.
#
#    VirtualRibosome is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    VirtualRibosome is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RevTrans; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# $Id: mod_translate.py,v 1.7 2005/06/09 09:58:54 raz Exp $
#
# $Log: mod_translate.py,v $

"""
mod_translate.py

The translation library behind the Virtual Ribosome and the RevTrans server.
"""

import sys,string

class TransTableRec:
	def __init__(self):
		self.description  = ""
		self.d_all        = {}
		self.d_first      = {}
		
	def toString(self):
		return "Description: %s\nd_all: %s\nd_first: %s" % (self.description,str(self.d_all),str(self.d_first))

# Standard genetic code:
d = {}
# nTn
d["TTT"] = d["TTC"] = "F"
d["TTA"] = d["TTG"] = d["CTT"] = d["CTC"] = d["CTA"] = d["CTG"] = "L"
d["ATT"] = d["ATC"] = d["ATA"] = "I"
d["ATG"] = "M"
d["GTT"] = d["GTC"] = d["GTA"] = d["GTG"] = "V"

# nCn
d["TCT"] = d["TCC"] = d["TCA"] = d["TCG"] = "S"
d["CCT"] = d["CCC"] = d["CCA"] = d["CCG"] = "P"
d["ACT"] = d["ACC"] = d["ACA"] = d["ACG"] = "T"
d["GCT"] = d["GCC"] = d["GCA"] = d["GCG"] = "A"

# nAn
d["TAT"] = d["TAC"] = "Y"
d["TAA"] = d["TAG"] = "*" 	#Stop
d["CAT"] = d["CAC"] = "H"
d["CAA"] = d["CAG"] = "Q"
d["AAT"] = d["AAC"] = "N"
d["AAA"] = d["AAG"] = "K"
d["GAT"] = d["GAC"] = "D"
d["GAA"] = d["GAG"] = "E"

# nGn
d["TGT"] = d["TGC"] = "C"
d["TGA"] = "*"  			#Stop
d["TGG"] = "W"
d["CGT"] = d["CGC"] = d["CGA"] = d["CGG"] = "R"
d["AGT"] = d["AGC"] = "S"
d["AGA"] = d["AGG"] = "R"
d["GGT"] = d["GGC"] = d["GGA"] = d["GGG"] = "G"

dStdRec = TransTableRec()
dStdRec.d_all = dStdRec.d_first = d

iupac = {}

iupac["A"] = "A"
iupac["C"] = "C"
iupac["G"] = "G"
iupac["T"] = "T"

iupac["R"] = "AG"	#puRine
iupac["Y"] = "CT"	#pYrimidine
iupac["M"] = "CA"
iupac["K"] = "TG"	#Keto
iupac["W"] = "AT"	#Weak
iupac["S"] = "GC"	#Strong
iupac["B"] = "CGT"	#not A
iupac["D"] = "AGT"	#not C
iupac["H"] = "ACT"	#not G
iupac["V"] = "ACG"	#not T ("U" is Uracile)

iupac["N"] = "ACGT"	#aNy

alphaDNA       = "ACGTRYMKWSBDHVN" 
alphaDNAStrict = "ACGT"
alphaPep       ="*ACDEFGHIKLMNPQRSTVWY"

DEBUG = 1;

d_ncbi_table = {}

def parseNcbiTable(lines):
	result = {}
	tab_id = ""
	desc = aa_all = aa_first = ""
	dRec = TransTableRec()
	for line in lines.split("\n"):
		line = line.strip()
		#print "!!!"+line
		
		if   line.startswith("name ") and (desc == ""):
			dRec.description += line.split('"')[1]+" "
			
		elif line.startswith("id "):
			tab_id = line.split()[1]
			
		elif line.startswith("ncbieaa"):
			aa_all = line.split('"')[1]
			
		elif line.startswith("sncbieaa"):
			aa_first = line.split('"')[1]
		
			c = 0
			for b1 in "TCAG":
				for b2 in "TCAG":
					for b3 in "TCAG":
						codon = b1+b2+b3
						dRec.d_all[codon] = aa_all[c]
						aaf = aa_first[c]
						if aaf == "-": aaf = aa_all[c]
						dRec.d_first[codon] = aaf
						
						c += 1
			result[tab_id] = dRec
			dRec = TransTableRec()
			tab_id = ""
	return result

def parseMatrixLines(iterator):
	result = {}
	for line in iterator:
		line = line.strip()
		
		if not line: 
			continue                # Ignore blank lines
		if line.startswith("#"): 
			continue                # Ignore comment lines
			
		tokens = line.split()
		try:
			codon, aa = tokens
			
			# Skip invalid entries
			if len(codon) <> 3: 
				badCodon = 1
			else:
				codon = codon.upper().replace("U","T")
				for c in codon: 
					if not c in alphaDNAStrict: badCodon = 1
				badCodon = 0
				
			if badCodon:
				raise "Bad codon: %s [%s]" % (codon,line)
					
				
			if len(aa) <> 1: #or (not aa in alphaPep):
				raise "Bad aa: %s [%s]" % (aa,line)
				
			result[codon] = aa
				
		except Exception, e:
			if DEBUG:
				sys.stderr.write("Matrix Error - %s\n" % e)
	
	if len(d) != 64 and DEBUG:
		sys.stderr.write("Matrix Error - size of matrix differs from 64 [%i]\n" % len(d))
	
	return result
				

def parseMatrixFile(filename):
	if d_ncbi_table.has_key(filename):
		dRec = d_ncbi_table[filename]
		return dRec
		
	dRec = TransTableRec()
	dRec.d_all = parseMatrixLines(open(filename,"r").xreadlines())
	dRec.d_first = dRec.d_all
	dRec.description = "Custom translation table '%s'" % filename
	return dRec

def trim(seq):
	return seq.upper().replace("U","T")

def trim_old(seq):
	result = []
	for c in seq.upper():
		if c in alphaDNA: result.append(c)
		elif c == "U":    result.append("T")
		else:             result.append("N")
	return result

# Assumption: Degenerate codons are rare - speed is not an issue
def decode(codon,dRec,isFirst):
	if len(codon) != 3: return []
	
	#Use the relevant translation table
	if isFirst: d_gc = dRec.d_first
	else:       d_gc = dRec.d_all
	
	#Check for the simple case: this is a standard non-degenerate codon
	if d_gc.has_key(codon):
		return [ d_gc[codon]]
	
	#The codon is to some degree degenerate - start the whole recursive scheme
	result = []
	
	for i in range(0,3):
		p = iupac[codon[i]]
		if (len(p) > 1):
			for c in p:
				result += (decode(codon[0:i]+c+codon[i+1:3],dRec,isFirst))
			return result
		
		if len(p) == 0:
			return [] # Unknown/illegal char
	#return [d[codon]]
	
def condense(lst):
	result = []
	for e in lst:
		if not e in result: result.append(e)
	return result

def translate(seq,transRec):
	return translate(seq,transRec,True,True)
			
def translate(seq,transRec,firstIsStartCodon,readThroughStopCodon):
	debug = False
	if not transRec:
		transRec = dStdRec # Use the Standard Genetic Code if no custom matrix is supplied

	result = []
	seq = trim(seq)
	if firstIsStartCodon: isFirst = True
	else:                 isFirst = False 
	for i in range(0,len(seq),3):
		aa = condense(decode(seq[i:i+3],transRec,isFirst))
		if debug: print seq[i:i+3], aa
		if aa:
			if aa[0] == "*" and not readThroughStopCodon: 
				break
			
			if len(aa) == 1: result.append(aa[0])
			else:
				s = string.join(aa,"")
				if   s == "DN" or s == "ND": result.append("B") #N or D
				elif s == "EQ" or s == "QE": result.append("Z") #E or Q
				else :                       result.append("X") #Any
				#print seq[i:i+3]
		isFirst = False
		
	pepseq = "".join(result)		
	if debug:
		print pepseq
	return pepseq
	
# Annotate all possible start and stop codons
def annotate(seq,transRec):
	debug = False

	result = []
	seq = trim(seq)
	for i in range(0,len(seq),3):
		codon = seq[i:i+3]
		aa = condense(decode(codon,transRec,True))
		if debug: print seq[i:i+3], aa
		if aa:
			if   "*" in aa: 
				result.append("*")
			elif ["M"] == aa:
				aa_int = condense(decode(codon,transRec,False))
				if aa_int == ["M"]: result.append("M")
				else:               result.append("m")
			else:           
				result.append(".")	
		
	return "".join(result)			

#def translate(seq):
#	translate(seq,None)

try:
	import ncbi_genetic_codes
	d_ncbi_table = parseNcbiTable(ncbi_genetic_codes.ncbi_gc_table)
			
except:
	pass
			
if __name__ == "__main__":
	for line in sys.stdin.readlines():
		print translate(line,None,True,False)
