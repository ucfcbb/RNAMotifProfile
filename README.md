## RNAMotifProfile

**A graph-based approach to generate profiles from RNA structural motifs and search new instances**

* Md Mahfuzur Rahaman<sup>†</sup>, mahfuz at ucf dot edu
* Shaojie Zhang*<sup>†</sup>, shaojie dot zhang at ucf dot edu

<sup>†</sup>Department Computer Science, University of Central Florida, Orlando, FL 32816 USA \
*To whom correspondence should be addressed.

---

**Abstract**
RNA structural motifs are the recurrent segments in RNA 3D structures that play a crucial role in the functional diversity of RNAs. Understanding the similarities and variations within these recurrent motif groups is essential for gaining insights into RNA structure and function. While recurrent structural motifs are generally assumed to be composed of the same isosteric base interactions, this consistent pattern is not observed across all examples of these motifs. Existing methods for analyzing and comparing RNA structural motifs may overlook variations in base interactions and associated nucleotides. RNAMotifProfile is a novel profile-to-profile alignment algorithm that generates a comprehensive profile from a group of structural motifs, incorporating all base interactions and associated nucleotides at each position. By structurally aligning input motif instances using a guide-tree-based approach, RNAMotifProfile captures the similarities and variations within recurrent motif groups. Additionally, RNAMotifProfile can function as a motif search tool, enabling the identification of instances of a specific motif family by searching with the corresponding profile. The ability to generate accurate and comprehensive profiles for RNA structural motif families, and to search for these motifs, facilitates a deeper understanding of RNA structure-function relationships and potential applications in RNA engineering and therapeutic design.

### 1. Installation

#### Install prerequisites
RNAMotifProfile is implemented using Python in a 64-bit Linux environment. Additionally, it uses several Python libraries that are required to run RNAMotifProfile. These libraries are included in the [requirements.txt](requirements.txt) file. To install all required Python libraries, please navigate to the RNAMotifProfile home directory in the terminal and execute the following command:

```
pip install -r requirements.txt
```

### 2. Input Specifications

RNAMotifProfile takes input from a comma-separated file. Each line in the input file represents a motif family. The motif family starts with a name, followed by a comma-separated list of motifs (the indices for motifs are expected to be in the PBD index, but it can be changed to the FASTA index by setting *input_index_type* parameter in the configuration file). To see examples of input formats, please check the two sample input files ([sample1.in](sample1.in) and [sample2.in](sample2.in)) provided in the root directory.

While used in regular mode, RNAMotifProfile works as a profile building tool. However, it can be used as a motif search tool by using an *-m* parameter. To search for a specific type of RNA structures through an RNA chain, a profile file of that type in *\*.pfl* format needs to be provided along with an RNA chain. Comma-separated RNA chains can also be provided to search through multiple RNA chains.

By default, all input files are put in the *input* directory, all output files will be created under the *output* directory.

### 3. Parameter options

Provide input file containing motif locations:
```
-i <input_file>
```
Choose a subdirectory name under *output* directory:
```
-o <output_subdirectory>
```
Filter less occurring interactions and nucleotides while generating mean and standard deviation to use during *search* mode (will not be used in profile building process; providing a threshold value is optional - default is 0.1; value range 0.0-1.0):
```
-f1 <threshold_value>
```
Filter more occurring gaps to reduce noise generated due to extra gaps (providing a threshold value is optional - default is 75.0; value range 100.0-0.0):
```
-f2 <threshold_value>
```
Use *branch-and-bound* while building profile:
```
-bnb
```
Activate *search* mode to search the instances close to a profile:
```
-m search
```
Provide the profile file in \*pfl format to search similar instances in search mode:
```
-p <profile_file>
```
Provide the RNA chain information in which user wants to search for similar instances in *search* mode:
```
-c <PDB_ID><underscore><Chain_ID>
```
Calculate z-score for each query motif using mean and std of the used profile with the motifs used to generate this profile (in search mode):
```
-z
```

### 3. Usage examples

**Profile generation mode**
To build profile from *sample1.in* file and create output files in *sample* directory under *output* directory::

```
python3 main.py -i input/sample1.in -o sample1 -f1 -f2 -bnb
```

**Motif search mode**
Search *Kink-turn* instances using (previously build) *Kink-turn* profile in a 23S rRNA (PDB 1S72, chain 0). The z-score needs to be generated for the search results and the output should be created in the *search-1S72_0* subdirectory under the *output* directory:

```
python3 main.py -m search -p output/known_families_IL/Kink-turn_20240404-140856.pfl -c 1S72_0 -z -o search-1S72_0 -f1 -f2 -bnb
```

### 4. Output Specification

**Profile generation mode**
The generated profile will be provided in a text file in the "output_subdirectory" (if provided by user) under the [output](output/) directory. The name of the output file is generated by adding the timestamp with the family name provided in the input file to uniquely identify each outcome.

There will be five sections in the output profile - sequence, basepair, stacking, breakpoints, and total loops. The sequence part will have the frequency of one or more nucleotides in each particular position. The basepair part contains the frequency of each type of base-pairing interaction in all occurring nucleotide pairs. The base-stacking interactions will be in the stacking part following the same format. Finally, the breakpoints section will have the position of the breakpoints and the total loops section will contain the number of loops used to generate this profile. An example profile created from the sample input [sample1.in](sample1.in) is given below:

```
#info=sequence
U: 2, G: 2, A: 1
G: 4, C: 1
G: 4, U: 1
A: 5
G: 3, C: 1, A: 1
C: 3, G: 1, U: 1
G: 5
A: 4, U: 1
A: 3, U: 1, C: 1
G: 5
A: 5
A: 3, G: 2
C: 2, A: 1, U: 1, G: 1
#info=basepair
1-11	tSH - GG: 2, tSH - GA: 2, tSH - CA: 1
2-10	tSH - GA: 4, tWH - UA: 1
3-6	tWS - AG: 2, tSS - AG: 2
3-9	tHS - AG: 5
4-10	tSS - GA: 2, tSS - CA: 1
5-10	cSS - GA: 1
6-7	cSW - GA: 1
#info=stacking
0-1	upward: 5
1-2	upward: 5
1-12	outward: 3
2-3	upward: 2
3-7	inward: 3, upward: 2
3-10	outward: 5
4-6	outward: 5
5-6	upward: 5
7-9	outward: 2, upward: 2
10-11	upward: 5
11-12	upward: 5
2-9	inward: 4
#info=break_points
5, 13
#info=total_loops
5
```

**Motif search mode**

In the search result of the motif search mode, all the input motif instances with their corresponding alignment score against the given profile will be provided while the result will be sorted in decreasing order based on their alignment score.

A part of the search result with Kink-turn motif profile through a 23S rRNA (PDB 1S72, chain 0):
```
Motif	Profile aln score	Matching bp count	Z-score	
1S72_0:1147-1154_1213-1216	122.27	5	2.56
1S72_0:77-81_93-100	116.44	5	2.34
1S72_0:1587-1592_1602-1608	86.32	4	1.2
1S72_0:1312-1319_1338-1342	86.19	4	1.2
1S72_0:937-940_1026-1033	78.68	3	0.91
```

### ACKNOWLEDGEMENTS

RNAMotifProfile is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

For bug reports or comments please contact shaojie.zhang@ucf.edu
