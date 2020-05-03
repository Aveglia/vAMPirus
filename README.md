

                                    ,---.,-.-.,---.o
                              .    ,|---|| | ||---'.,---..   .,---.
                               \  / |   || | ||    ||    |   |`---.
                                `'  `   `` ` ``    ``    `---``---`
        A light-weight, automated virus amplicon sequencing analysis pipeline

## Overview
The end goal of the vAMPirus project is to generate a robust workflow that can be used consistently across various
virus amplicon projects.

```
 ___________________________              ________________________________              ________________________________________              _____________________
|                           |            |                                |            |                                        |            |                     |
| Read processing/filtering |  ------->  | Read merging and dereplication |  ------->  | ASV/OTU generation and chimera removal |  ------->  | Taxonomy assignment | --
|___________________________|            |________________________________|            |________________________________________|            |_____________________|

                  ________________________________                __________________________________________________________________________________
                 |                                |   Optional   | - Translate nucleotide fastas to protein  - Diamond blastp with translated seqs  |
        -------> | Generate ASV/OTU counts tables |  --------->  | -- Nucleotide/protein model testing       - Phylogenetic tree generation         |
                 |________________________________|    steps     |__________________________________________________________________________________|

```

### Contact/support:

Please contact Alex Veglia at ajv5@rice.edu with any feedback or questions. This is a **early stage** of the automated pipeline and any kind of input from the
community is welcomed and encouraged!

### How to cite:

If you do use vAMPirus for your analyses, please cite using the following ->

Veglia, A.J. *et.al.,* 2020


### So, what can you do with vAMPirus?

**Features in v0.1.0 ->**  
1. Capable of performing all read processing including adapter removal, primer removal, read correction, read merging and length/quality trimming/filtering  using the programs `fastp` (Chen et al. 2018) and `bbduk.sh` (DOE Joint Genome Institute)
2. FastQC analysis before and after read processing  
3. Performs denoising with the UNOISE version 3 algorithm by Robert Edgar (2010) in `vsearch` (Rognes et al. 2016) to generate ASV's representing biologically true "haplotypes" within the sample set which are then clustered by user-specified percent similarity into Operational Taxonomic Units (OTUs). Both the clustered and non-clustered ASV fastas are available to be used in OTU counts table generation.
4. Assigns taxonomy to OTU sequences using the diamond `blastx` with a specified virus protein database (note: vAMPirus was tested with the protein translated Reference Viral DataBase v18 (RVDBv18; Goodacre et al. 2018).  
5. Phyloseq-formatted summary tables are generated with the diamond output for easy loading and manipulation in `R`
6. Generate OTU counts tables in both .tsv and .biome format for downstream manipulation in R or `mothur` for visualization/statistical analyses.  
7. Translates nucleotide sequences for all OTU fastas generated during analysis  
8. Performs nucleotide and protein model testing using `Modeltest-NG` (Darriba et al. 2020)  
9. Performs nucleotide/protein alignment with `MAFFT` (Rozewicki et al. 2019) with subsequent phylogenetic analysis using `RAxML-NG` (Kozlov et al. 2019)  

**Features in development ->**
1. Integration into NextFlow workflow manager for use of pipeline across all operating systems  
2. More flexibility in RAxML commands  
3. Add the use of RefSeq header style for taxonomy assignment  
4. Add capability of producing preliminary plots of results  


**Dependencies (version of program vAMPirus-v0.1.0 was tested with) ->**    
1. [FastQC (v0.11.1)](https://github.com/s-andrews/FastQC)
2. [fastp (v0.20.1)](https://github.com/OpenGene/fastp)
3. [BBtools](https://jgi.doe.gov/data-and-tools/bbtools/) -> vAMPirus uses the `bbduk.sh` and the `reformat.sh` scripts within BBtools
4. [vsearch (v2.14.2)](https://github.com/torognes/vsearch)
5. [diamond (v0.9.30)](https://github.com/bbuchfink/diamond)
6. [ModelTest-NG (vx.y.z)](https://github.com/ddarriba/modeltest)  
7. [MAFFT (v7.450)](https://mafft.cbrc.jp/alignment/software/)  
8. [RAxML-NG (v0.9.0)](https://github.com/amkozlov/raxml-ng)
9. [Virtual Ribosome-dna2pep (v.1.1)](http://www.cbs.dtu.dk/services/VirtualRibosome/download.php)
10. Python -> python2.4(and up) and python3(and up) are needed for full analysis.  
NOTE: if not interested in protein analysis, then only python3 is needed.  
