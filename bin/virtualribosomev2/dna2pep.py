#!/usr/bin/env python3

# Script was modified form original to use python3

#    Copyright 2006 Rasmus Wernersson, Technical University of Denmark
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
NAME
	dna2pep - full featured computational translation of DNA to peptide.
	(The program behind the "Virtual Ribosome" webserver.

SYNOPSIS
	dna2pep [options] [input files] [-f outfile]

DESCRIPTION
	TRANSLATION: The translation engine of dna2pep has full support for handling
	degenerate nucleotides (IUPAC definition, e.g. W = A or T, S = G or C).
	All translation table defined by the NCBI taxonomy group is included,
	and a number of options determining the behaviour of STOP and START
	codons is avialable.

	INTRON and EXONS: dna2pep natively understands TAB files containing
	Intron/Exon annotation (gb2tab / FeatureExtract). When translating
	files containing Intron/Exon structure, dna2pep will annotate the
	underlying gene-structure in the annotation of the translated
	sequence.

	Input files can be in FASTA (no Intron/Exon annotation) RAW (single
	sequence with no header - all non-letters are discarded) or TAB
	(incluing annotation) FORMAT. The output format will by default be FASTA
	for files without annotation and TAB for files including annotation.
	The file format is autodetected by investigating the first line
	of the input.

	If no input files are specified, dna2pep will read from STDIN.

OPTIONS
	-F, --outfile
		Optional - specify an output file. If no output file is
		specified the output will go to STDOUT.

	-O, --outformat
		Specify output format (see also the --fasta, --tab,
		--report options below):

		FASTA:  Fasta format (plain DNA, no sequence annotation)

		TAB:    Tab format. Each line contains the following four
		        fields, separated by tabs:
			name, seq, ann, comment

			See gb2tab (FeatureExtract) for details.

		REPORT:	A nice visualization of the results.

		AUTO:  [Default] Generate a both a report and sequence output
			(use the same format as the one detected from the for
			the input files).

	--fasta filename
		Write output sequences in FASTA format to the specified file.
		Use '-' to indicate STDOUT.

	--tab filename
		Write output sequences in TAB format to the specified file.
		Use '-' to indicate STDOUT.

	--report filename
		Write report to the specified file.
		Use '-' to indicate STDOUT.

        -m, --matrix tablename/file
                Use alternative translation matrix instead of the build-in
                Standard Genetic Code for translation.

                If "tablename" is 1-6,9-16 or 21-23 one of the alternative
                translation tables defined by the NCBI taxonomy group will be
                used.

                Briefly, the following tables are defined:
                -----------------------------------------
                 1: The Standard Code
                 2: The Vertebrate Mitochondrial Code
                 3: The Yeast Mitochondrial Code
                 4: The Mold, Protozoan, and Coelenterate Mitochondrial Code
                    and the Mycoplasma/Spiroplasma Code
                 5: The Invertebrate Mitochondrial Code
                 6: The Ciliate, Dasycladacean and Hexamita Nuclear Code
                 9: The Echinoderm and Flatworm Mitochondrial Code
                10: The Euplotid Nuclear Code
                11: The Bacterial and Plant Plastid Code
                12: The Alternative Yeast Nuclear Code
                13: The Ascidian Mitochondrial Code
                14: The Alternative Flatworm Mitochondrial Code
                15: Blepharisma Nuclear Code
                16: Chlorophycean Mitochondrial Code
                21: Trematode Mitochondrial Code
                22: Scenedesmus obliquus mitochondrial Code
                23: Thraustochytrium Mitochondrial Code

                See http://www.ncbi.nlm.nih.gov/Taxonomy [Genetic Codes]
                for a detailed description. Please notice that the table
                of start codons is also used (see the --allinternal option
                below for details).

                If a filename is supplied the translation table is read from
                file instead.

                The file should contain one line per codon in the format:

                codon<whitespace>aa-single letter code

                All 64 codons must be included. Stop codons is specified
                by "*". T and U is interchangeable. Blank lines and lines
                starting with "#" are ignored.

                See the "gcMitVertebrate.mtx" file in the dna2pep source
                distribution for a well documented example.

	-r x, --readingframe=x
		Specify the reading frame. For input files in TAB format this
		options is ignored, and the reading frame is build from the
		annotated Intron/Exon structure.

		 1: Reading frame 1 (e.g. ATGxxxxxx). DEFAULT.
		 2: Reading frame 2 (e.g. xATGxxxxx).
		 3: Reading frame 3 (e.g. xxATGxxxx).

		-1: Reading frame 1 on the minus strand.
		-2: Reading frame 2 on the minus strand.
		-3: Reading frame 3 on the minus strand.

		all: 	Try all reading frames.
		     	This option also implies the -x option.

		plus: 	All positive reading frames.
		     	This option also implies the -x option.

		minus: 	All negative reading frames.
			This option also implies the -x option.

	-o mode, --orf mode
		Report longest ORF in the reading frame(s) specified with the
		-r option.

		Mode governs which criterias are used to allow the opening of
		an ORF. "Strict start codons" => codons _always_ coding for
		methione (e.g. ATG in the standard code), "Minor start codons"
		=> codon only coding for methionine at the start positon
		(e.g. TTG in the standard genetic code).

		Mode can be:
		------------
		strict:		Open an ORF at "strict start codons" only.
		any:		Open an ORF at any start codon.
		none:		Do not use start codons - look for the longest
				fragment before a STOP codon.

		The DNA fragment usedfor encoding the ORF will be added to the
		comment field (TAB format only).

        -a, --allinternal
                By default the very first codon in each sequences is assumed
                to be the initial codon on the transcript. This means certain
                non-methionine codons actually codes for metionine at this
                position. For example "TTG" in the standard genetic code (see
                above).

                Selecting this option treats all codons as internal codons.

        -x, --readthroughstop
                Allow the translation to continue after a stop codon is reached.
                The stop codon will be marked as "*".

	-p, --plain, --ignoreannotation
		Ignore annotation for TAB files. If this options is selected
		TAB files will be treated in same way as FASTA files.

	-c, --comment
		Preserve the comment field in TAB files. Normally the comment
		field is silently dropped, since it makes no sense for FASTA
		files.

	-C, --processcomment
		Works as the -c option described above, except a bit of intelligent
		parsing is done on the comment field: If a "/spliced_product"
		sub-field is found (from TAB files create by FeatureExtract / gb2tab)
		only the part of the comment field before the DNA specific information
		is kept in the comment field.

	-e, --exonstructure
		Default for TAB files. Annotate the underlying exons structure
		of the translated sequence the following way. Positions that
		are fully or partially encoded within the first exon get the
		annotation character "1", positions in the secon exon get the
		character "2" etc.

		The hex-decimal system is used, which means up to 15 exons can
		be uniquely annotated, before the numbering wraps around to "0".

	-i, --intronphase
		Annotate where an intron interrupted the DNA sequences, and how
		the intron did cut the readingframe.

		0 : phase-0 intron (inbetween the previous and current position).
		1 : phase-1 intron.
		2 : phase-2 intron.


AUTHOR
	Rasmus Wernersson, raz@cbs.dtu.dk
	Feb-Mar 2006

FILES
	dna2pep.py, mod_translate.py, ncbi_genetic_codes.py

WEB PAGE
	http://www.cbs.dtu.dk/services/VirtualRibosome/

REFERENCE
	Rasmus Wernersson
	Virtual Ribosome - Comprehensive DNA translation tool.
	Submitted to Nucleic Acids Research, 2006
"""

import sys, re, mod_translate,string
from optparse import OptionParser

validDNA = "ATUGCYRSWKMBDHVN"
complDNA = "TAACGRYWSMKVHDBN"

allValid = validDNA+validDNA.lower()
transTable = str.maketrans(validDNA+validDNA.lower(),complDNA+complDNA.lower())

pwidth = 90

REPORTHEADER = """
VIRTUAL RIBOSOME
----------------
"""

ORF_ANNOTATION = """
>>>   :   START codon (strict)
)))   :   START codon (alternative)
***   :   STOP
"""

def makePretty(title,vals,labels,max_len):
	l = []
	l.append(">%s\n\n" % title)
	for i in range(0,max_len,pwidth):
		for j in range(0,len(vals)):
			val = vals[j]
			lab = labels[j]
			pos = min(max_len,i+pwidth)
			if lab:
				spos = "%d" % pos
			else:
				spos = ""

			s = "%-2s %s %s\n" % (lab,val[i:pos],spos)
			l.append(s)
		l.append("\n")
	return l

def explodePep(s):
	l = list(s)
	return " "+"  ".join(l)+" "

def revCom(dna):
	l = list(dna)
	l.reverse()
	rev_dna = "".join(l)
	return rev_dna.translate(transTable)

def dnaComplement(dna):
	return dna.translate(transTable)

def revStr(s):
	l = list(s)
	l.reverse()
	return "".join(l)


def isDNAValid(dna):
	for c in dna:
		if not c in allValid: return False

	return True

def combineToTab(name, seql):
	seq = "".join(seql)
	ann= "."*len(seq)
	com = ""
	return (name, seq, ann, com)

# Read RAW (single entry, no name etc).
def readRaw(lines):
	l = []
	seql = []
	for line in lines:
		sl = []
		for c in line.strip().upper():
			if c.isalpha(): sl.append(c)

		seql.append("".join(sl))

	l.append(combineToTab("Seq1",seql))
	return l

# Read FASTA format
def readFasta(lines):
	l = []

	name = ""
	seql = []

	for line in lines:
		line = line.strip()
		if line.startswith(">"):
			if name: l.append(combineToTab(name,seql))
			name = line[1:]
			if not name: name ="NoNameSeq"
			seql = []
		else:
			seql.append(line)

	if name: l.append(combineToTab(name,seql))

	return l

# read TAB format
def readTab(lines):
	l = []

	for line in lines:
		line = line.strip()
		tokens = line.split("\t")
		if len(tokens) < 2: continue

		name = tokens[0]
		seq  = tokens[1]

		if len(tokens) > 2:
			ann = tokens[2]
		else:
			ann = "."*len(seq)

		if len(tokens) > 3:
			com = tokens[3]
		else:
			com = ""

		l.append( (name,seq,ann,com) )

	return l


def readInput(lines):
	if not lines: return ([], True)

	line = lines[0]
	tokens = line.split("\t")


	if line.startswith(">") and len(tokens) < 3:
		l = readFasta(lines)
		isFasta = True

	elif len(tokens) == 1:
		l = readRaw(lines)
		isFasta = True

	else:
		l = readTab(lines)
		isFasta = False

	return (l, isFasta)

def writeFasta(seqs, outstream):
	for (name, seq, ann, com) in seqs:
		print(">"+name, file=outstream)
		for i in range(0,len(seq),60):
			print(seq[i:i+60], file=outstream)

def writeTab(seqs, outstream):
	for tokens in seqs:
		print("\t".join(tokens), file=outstream)

def openForWriteOrDie(outfile):
	try:
		outstream = open(outfile,"w")

	except IOError as xxx_todo_changeme:
		(strerror) = xxx_todo_changeme
		print("ERROR - cannot write to the specified file %s [%s]" % (outfile,strerror), file=sys.stderr)
		sys.exit(-1)

	return outstream

def parseOpts():
	# Quick hack to overrule the -h and --help feature
	# build into the optpase module
	if "-h" in sys.argv or "--help" in sys.argv:
		print(__doc__)
		sys.exit(0)

	parser = OptionParser()

	# File handling
	parser.add_option("-F","--outfile",   type="string", dest="outfile",    default="")
	parser.add_option("-O","--outformat", type="string", dest="outformat",  default="AUTO")
	parser.add_option("--tab",            type="string", dest="tabfile",    default="")
	parser.add_option("--fasta",          type="string", dest="fastafile",  default="")
	parser.add_option("--report",         type="string", dest="reportfile", default="")

	# Matrix
	parser.add_option("-m","--matrix", type="string", dest="matrix", default="1")

	# Reading frame
	parser.add_option("-r","--readingframe", type="string", dest="readingframe", default="1")

	# ORF finding
	parser.add_option("-o","--orf", type="string", dest="orf", default="")

	# Stop and Start codons
	parser.add_option("-a","--allinternal",     action="store_true", dest="allinternal",     default = False)
	parser.add_option("-x","--readthroughstop", action="store_true", dest="readthroughstop", default = False)

	# Annotation
	parser.add_option("-p","--plain",              action="store_true", dest="ignoreann", default=False)
	parser.add_option("-e","--exonstructure",      action="store_true", dest="exonann",   default=True)
	parser.add_option("-i","--intronphase",        action="store_true", dest="intronrf",  default=False)

	# Comments
	parser.add_option("-c","--comment",        action = "store_true", dest = "keepcomment",    default=False)
	parser.add_option("-C","--processcomment", action = "store_true", dest = "processcomment", default=False)

	# Debug
	parser.add_option("-d","--debug", action="store_true", dest="debug", default=False)

	(opt, args) = parser.parse_args()

	# Check reading frame
	if not opt.readingframe in ["1", "2", "3", "-1", "-2", "-3", "all", "plus","minus"]:
		sys.stderr.write("Invalid reading frame [%s]\n" % opt.readingframe)
		sys.exit(-1)

	if opt.readingframe in ["all","plus","minus"]:
		opt.readthroughstop = True

	# Chech ORF mode
	if not opt.orf in ["","strict","any","none"]:
		print("Invalid ORF mode [%s]\n" % opt.orf, file=sys.stderr)
		sys.exit(-1)

	# Chech output format
	opt.outformat = opt.outformat.upper()
	if not opt.outformat in ["AUTO","TAB","FASTA","REPORT"]:
		print(file=sys.stderr("Invalid output format [%s]\n" % opt.outformat))
		sys.exit(-1)

	# Check mutually exclusive options
	if opt.intronrf: opt.exonann=False

	return (opt, args)

if __name__ == "__main__":
	reports = []
	pepseqs = []
	opt, args = parseOpts()

	# Initialize translation matrix
	mtx = mod_translate.parseMatrixFile(opt.matrix)
	if not mtx:
		sys.stderr.write("Cannot initialize matrix [%s]\n" % matrix)
		sys.exit(-1)

	# Start the translation report
	reports.append( [REPORTHEADER,"Translation table: %s\n\n" % mtx.description] )

	# Read input data
	if not args:
		sys.stderr.write("Reading from STDIN. -h for help\n")
		lines = sys.stdin.readlines()
	else:
		lines = []
		for fn in args:
			try:
				lines += (open(fn,"r").readlines())
			except:
				sys.stderr.write("ERROR: Cannot read from file '%s'\n" % fn)
				sys.exit(-1)

	try:
		(seq_list, isFasta) = readInput(lines)
	except Exception as msg:
		print("ERROR parsing input files. Please verify the format (FASTA, RAW or TAB)", file=sys.stderr)
		print("[%s]" % (str(msg)), file=sys.stderr)
		sys.exit(-1)

	# Ignore annotation?
	if opt.ignoreann:
		isFasta = True

	# Process the sequences
	rx = re.compile("(\(E+\))")

	for (name, dna, ann, com) in seq_list:
		seq_proc = ann_proc = exnum_proc = ""

		# Test if the DNA sequence is valid
		if not isDNAValid(dna):
			sys.stderr.write("Non IUPAC characters is detected in sequence '%s' - skipping this entry\n" %name)
			#sys.stderr.write("seq: %s\n" % seq)
			continue

		# Files without intron/exon annotation -------------------------------------------
		if isFasta:
#			reports.append([ORF_ANNOTATION])

			d_collect = {}
			if   opt.readingframe in ["all"]:
				rf_list = ["1","2","3","-1","-2","-3"]
				echo_rf = True

			elif opt.readingframe == "plus":
				rf_list = ["1","2","3"]
				echo_rf = True

			elif opt.readingframe == "minus":
				rf_list = ["-1","-2","-3"]
				echo_rf = True

			else:
				rf_list = [opt.readingframe]
				echo_rf = False

			for rf in rf_list:

				# Find current reading frame
				if   rf == "1":
					qseq = dna
				elif rf == "2":
					qseq = dna[1:]
				elif rf == "3":
					qseq = dna[2:]
				elif rf == "-1":
					qseq = revCom(dna)
				elif rf == "-2":
					qseq = revCom(dna)[1:]
				elif rf == "-3":
					qseq = revCom(dna)[2:]
				else:
					qseq = dna

				# Do the actual translation
				pep = mod_translate.translate(qseq,mtx,not opt.allinternal,opt.readthroughstop)
				pa  = mod_translate.annotate(qseq,mtx)

				# The annotation string may be longer that the peptide, if the -x more is not used
				pa = pa[:len(pep)]

				# Store translated sequence
				if echo_rf:
					cname = name+"_rframe"+rf
				else:
					cname = name

				data = ( cname,pep,pa,qseq )
				d_collect[rf] = data

			# Do ORF finding?
			if opt.orf:
				#Find longest ORF
				bestlen  = -1
				bestspan = (0,0)
				bestdata = None
				bestrf   = ""
				for key in list(d_collect.keys()):
					data = d_collect[key]
					seq = data[1]
					ann = data[2]
					j   = -1
					m   = -1
					while(True):
						#Step 1 - find start
						if opt.orf == "none":
							j = m+1
							if j >= len(ann): break
						else:
							j_strict = ann.find("M",j+1)
							j_any    = ann.find("m",j+1)
							if opt.orf == "strict":
								j = j_strict
							else:
								if   j_strict == -1:
									j = j_any
								elif j_any == -1:
									j = j_strict
								else:
									j = min(j_strict,j_any)

							if j == -1: break

						#Step 2 - find stop
						m = ann.find("*",j)
						if m == -1: m = len(ann)

						#print j,m

						if (m - j) > bestlen:
							bestlen  = m - j
							bestspan = (j, m)
							bestdata = data
							bestrf   = key

							#print rf, bestlen
							#print seq[j:m+1]
							#print ann[j:m+1]

						j = m
				# Format the best hit
				if not bestdata:
					bestdata = list(d_collect.values())[0]
					bestspan = (0,0)
					bestrf   = list(d_collect.keys())[0]

					msg      = "NO ORF FOUND (given the criteria '%s') for sequence '%s'\n\n" % (opt.orf,name)
					reports.append([msg])

#				clist = list(bestdata[1].lower())
				name       = bestdata[0]
				pep        = bestdata[1]
				ann        = bestdata[2]
				dna_work   = bestdata[3]

				bpos, epos = bestspan
				orf_dna    = dna_work[bpos*3:epos*3]
				orf        = mod_translate.translate(orf_dna,mtx,True,False)
				new_pep    = " "*bpos + orf + " "*(len(pep)-epos)

				name += "_ORF"

				d_collect = {}
				d_collect[bestrf] = (name,new_pep,ann,orf_dna)
				#print d_collect

			# Processing and Pretty printing...
			if opt.readingframe in ["all","plus","minus"] and (not opt.orf):
				if   opt.readingframe == "all":
					doPlus  = True
					doMinus = True
				elif opt.readingframe == "plus":
					doPlus  = True
					doMinus = False
				else:
					doPlus  = False
					doMinus = True

				dna_plus  = dna
				dna_minus = revCom(dna)

				#dnapl = list(dna_plus)
				#dnaml = list(dna_minus)
#				dnapann = [" "]*len(dna_plus)
#				dnamann = [" "]*len(dna_minus)
				dnapann = ["."]*len(dna_plus)
				dnamann = ["."]*len(dna_minus)

				vals   = []
				labels = []

				if doPlus:
					pep1  =  d_collect["1"][1]
					ann1  =  d_collect["1"][2]
					pep1e = explodePep(pep1)

					pep2  =  d_collect["2"][1]
					ann2  =  d_collect["2"][2]
					pep2e = " "+explodePep(pep2)

					pep3  =  d_collect["3"][1]
					ann3  =  d_collect["3"][2]
					pep3e =  "  "+explodePep(pep3)

					anns = [ann1,ann2,ann3]
					peps = [pep1,pep2,pep3]
					for i in range(0,len(pep1)):
						for j in range(0,3):
							ann = anns[j]
							if i >= len(ann): continue
							ac = ann[i]
							if ac in ["M","m","*"]:
								dnapos = (i*3) + j
								if   ac == "M": c = ">"
								elif ac == "m": c = ")"
								else:           c = "*"
								for k in range(dnapos,dnapos+3):
									#dnapl[k] = dnapl[k].upper()
									dnapann[k] = c

					#dna_plus = "".join(dnapl)
					vals   += [pep3e,pep2e,pep1e,dna_plus,"".join(dnapann)]
					labels += ["","","","5'",""]

				if doMinus:
					pepm1  = d_collect["-1"][1]
					annm1  = d_collect["-1"][2]
					pepm1e = explodePep(pepm1)+"  "+"   "
					pepm1e = pepm1e[:len(dna_minus)]

					pepm2  = d_collect["-2"][1]
					annm2  = d_collect["-2"][2]
					pepm2e = " "+explodePep(pepm2)+" "+"   "
					pepm2e = pepm2e[:len(dna_minus)]

					pepm3  = d_collect["-3"][1]
					annm3  = d_collect["-3"][2]
					pepm3e = "  "+explodePep(pepm3)+"   "
					pepm3e = pepm3e[:len(dna_minus)]

					anns = [annm1,annm2,annm3]
					peps = [pepm1,pepm2,pepm3]
					maxlen = max(len(pepm1),len(pepm2))
					maxlen = max(maxlen,len(pepm3))

					for i in range(0,maxlen):
						for j in range(0,3):
							ann = anns[j]
							if i >= len(ann): continue
							ac = ann[i]
							if ac in ["M","m","*"]:
								dnapos = (i*3) + j
								if   ac == "M": c = "<"
								elif ac == "m": c = "("
								else:           c = "*"
								for k in range(dnapos,dnapos+3):
									#dnaml[k] = dnaml[k].upper()
									dnamann[k] = c

					#dna_minus = "".join(dnaml)
					dna_ann_str = "".join(dnamann)
					vals   += [revStr(dna_ann_str),revStr(dna_minus),revStr(pepm1e),revStr(pepm2e),revStr(pepm3e)]
					labels += ["","3'","","",""]
#					vals   += [revStr(dna_minus),revStr(dna_ann_str),revStr(pepm1),revStr(pepm2),revStr(pepm3)]
#					labels += ["3'","","","",""]

#				vals = [pep3,pep2,pep1,dna_plus,revStr(dna_minus),revStr(pepm1),revStr(pepm2),revStr(pepm3)]
#				labels = ["","","","5'","3'","","",""]
				title = "%s - reading frame(s): %s" % (name,opt.readingframe)

				l = makePretty(title,vals,labels,len(dna_plus))

				###print "".join(l)
				reports.append(l)

			else:
				rf_list = list(d_collect.keys())
				rf      = rf_list[0]

				vals   = []
				labels = []
				title  = "%s\nReading frame: %s" % (name,rf)

				#rf_int = abs(int(rf)) - 1
				pep    = d_collect[rf][1]
				pepx   = explodePep(pep)

				vals   = []
				labels = []

				vals.append(pepx)
				labels.append("")

				pa = d_collect[rf][2]
				pax = ["."]*len(pa)*3
				if rf.startswith("-"):
					dna_minus = revCom(dna)
					#print len(dna_minus)
					dnamann   = ["."] * len(dna_minus)
					if   rf == "-1":
						pepm  = d_collect["-1"][1]
						annm  = d_collect["-1"][2]
						pepme = explodePep(pepm)+"  "+"   "
						j     = 0
					elif rf == "-2":
						pepm  = d_collect["-2"][1]
						annm  = d_collect["-2"][2]
						pepme = " "+explodePep(pepm)+"  "+"   "
						j     = 1
					else:
						pepm  = d_collect["-3"][1]
						annm  = d_collect["-3"][2]
						pepme = "  "+explodePep(pepm)+"  "+"   "
						j     = 2

					pepme = pepme[:len(dna_minus)]

					for i in range(0,len(pepm)):
						ac = annm[i]
						if ac in ["M","m","*"]:
							dnapos = (i*3) + j
							if   ac == "M": c = "<"
							elif ac == "m": c = "("
							else:           c = "*"
							for k in range(dnapos,dnapos+3):
								dnamann[k] = c

					dna_ann_str = "".join(dnamann)
					#print len(dna_ann_str)
					#print revStr(dna_ann_str)+"'"
					#print revStr(dna_minus)  +"'"
					#print revStr(pepme)      +"'"
					vals   = [revStr(dna_ann_str),revStr(dna_minus),revStr(pepme)]
					labels = ["","3'",""]

				else:
					dna_plus = dna
					if   rf == "1":
						pepx = explodePep(pep)
						j    = 0
					elif rf == "2":
						pepx = " "+explodePep(pep)
						pepx = pepx[:len(dna_plus)]
						j    = 1
					else:
						pepx = "  "+explodePep(pep)
						pepx = pepx[:len(dna_plus)]
						j    = 2
					dnapann = ["."]*len(dna_plus)
					for i in range(0,len(pep)):
						c = "."
						ac = pa[i]
						if ac in ["M","m","*"]:
							dnapos = (i*3) + j
							if   ac == "M": c = ">"
							elif ac == "m": c = ")"
							else:           c = "*"
							for k in range(dnapos,dnapos+3):
								dnapann[k] = c

					ann = "".join(dnapann)
					vals   = [pepx,dna_plus,ann]
					labels = ["","5'",""]


				l = makePretty(title,vals,labels,len(dna))
				reports.append(l)
				###print "".join(l)


			for rf in rf_list:
				name,seq,ann,com = d_collect[rf]
				new_com = ""

				# Seq may contain leding and trailing spaces if we are in ORF finde mode
				if opt.orf:
					new_seq = []
					new_ann = []
					for i in range(0,len(seq)):
						if seq[i] != " ":
							new_seq.append(seq[i])
							new_ann.append(ann[i])

					seq = "".join(new_seq)
					ann = "".join(new_ann)

					new_com = '/orf_mode="%s"; /dna="%s";' % (opt.orf,com)

				pepseqs.append( (name,seq,ann,new_com) )
#				pepseqs.append( (cname,pep,pa,com.strip()) )

		# Files with intron/exon annotation ---------------------------------------------
		else:
			exon_count = 1
			for mo in rx.finditer(ann):
				start, end = mo.span()
				seq_proc += dna[start:end]
				ann_proc += ann[start:end]

				exnum_chr = "%x" % exon_count
				exnum_proc += exnum_chr * (end-start)
				exon_count = (exon_count + 1) % 0x10

			# DEBUG
			if opt.debug:
				print(seq_proc)
				print(ann_proc)
				print(rf_proc)
				print(exnum_proc)

			# Process comments
			if opt.keepcomment or opt.processcomment:
				if opt.processcomment:
					com = com[0:com.find("/spliced_product")]
			else:
				com = ""

			# Do the actual translation
			pep = mod_translate.translate(seq_proc,mtx,not opt.allinternal,opt.readthroughstop)

			# Calculate the pep annotation
			pa = ["."]*len(pep)

			if opt.exonann:
				for i in range(0,len(pep)):
					pa[i] = exnum_proc[i*3]

			if opt.intronrf:
				rf_proc = "012" * (len(ann_proc) / 3)
				rf_proc += "012"[:len(ann_proc) % 3]

				for i in range(1,len(pa)*3): # Skip first position
					if ann_proc[i] == "(":
						pa[i/3] = rf_proc[i]

			# Store translated sequence
			prot_ann = "".join(pa)
			pepseqs.append( (name,pep,prot_ann,com.strip()) )

			# Pretty printing for the report
			title = "%s - " % (name)
			if   opt.exonann:
				title += "translation and annotation of the exonic structure"

			elif opt.intronrf:
				title += "translation and annotation of the position and phase of the introns"

			vals   = [pep,prot_ann]
			labels = ["pep:","ann:"]
			l = makePretty(title,vals,labels,len(pep))
			reports.append(l)


	# Output the results -------------------------------------------------------
	# Step 1) Combined results.
	if opt.outfile: outstrem = openForWriteOrDie(opt.outfile)
	else:		outstream = sys.stdout

	if opt.outformat in ["REPORT","AUTO"]:
		for l in reports: outstream.writelines(l)
		print("//", file=outstream)

	if   "TAB" == opt.outformat:
		outFasta = False
	elif "FASTA" == opt.outformat:
		outFasta = True
	else:
		outFasta = isFasta

	if outFasta:
		writeFasta(pepseqs,outstream)
	else:
		writeTab(pepseqs,outstream)

	#outstream.close()

	# Step 2) Write specifik sub-result if requested
	if opt.reportfile:
		if opt.reportfile == "-": 	outstream = sys.stdout
		else:                	 	outstream = openForWriteOrDie(opt.reportfile)

		for l in reports: outstream.writelines(l)

	if opt.fastafile:
		if opt.fastafile == "-": 	outstream = sys.stdout
		else:				outstream = openForWriteOrDie(opt.fastafile)

		writeFasta(pepseqs,outstream)

	if opt.tabfile:
		if opt.tabfile == "-":		outstream = sys.stdout
		else:				outstream = openForWriteOrDie(opt.tabfile)

		writeTab(pepseqs,outstream)
