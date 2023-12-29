## RNAMotifProfile

**A graph-based approach to generate profiles from RNA structural motifs and search new instances**

* Md Mahfuzur Rahaman<sup>†</sup>, mahfuz at ucf dot edu
* Shaojie Zhang*<sup>†</sup>, shaojie dot zhang at ucf dot edu

<sup>†</sup>Department Computer Science, University of Central Florida, Orlando, FL 32816 USA \
*To whom correspondence should be addressed.

---

**Abstract**
The knowledge of RNA 3D structures is essential in understanding the significance of their various biological functions in cell biology. The recurrent segments of these RNA 3D structures are mostly responsible for the functional diversity that can be used in drug discovery and RNA-based therapeutics. Ideally, these recurrent segments are formed by the same set of isosteric base interactions but that is not required to be true for all cases. Having an interaction profile will definitely provide a better understanding of the similarities and variations in such recurrent motif groups. However, it is very challenging to achieve such structural profiles due to high computational complexity. In this work, we present RNAMotifProfile, a profile-to-profile alignment algorithm that generates a profile from a 3D structural motif group. It considers each motif as an individual profile and structurally aligns the input motif instances by following a guide-tree-based approach. RNAMotifProfile incorporates all the base interactions along with their associated nucleotides in each position. All related information of a profile is stored in a strategically designed data structure. Additionally, RNAMotifProfile can perform as a motif search tool to identify the instances of a motif family when searched with a profile of that particular family.

### 1. Installation

#### Install prerequisites
RNAMotifProfile is implemented using Python in a 64-bit Linux environment. Additionally, it uses several Python libraries that are required to run RNAMotifProfile. These libraries are included in the [requirements.txt](requirements.txt) file. To install all required Python libraries, please navigate to the RNAMotifProfile home directory in the terminal and execute the following command:

```
pip install -r requirements.txt
```

### 2. Input Specifications

RNAMotifProfile takes input from a comma-separated file. Each line in the input file represents a motif family. The motif family starts with a name, followed by a comma-separated list of motifs (the indices for motifs are expected to be in the PBD index, but it can be changed to the FASTA index by setting a parameter in the configuration file). To see examples of input formats, please check the two sample input files ([sample1.in](sample1.in) and [sample2.in](sample2.in)) provided in the root directory.

While used in normal mode, RNAMotifProfile works as a profile generation tool. However, it can be used as a motif search tool by using an extra parameter. In the profile generation mode, the input file represents the motifs for which profiles will be generated. On the other hand, in motif search mode, that input file will represent the query motif instances. A profile file (\*.pfl) is required to be provided in the search mode against which the input instances will be aligned. A manually prepared consensus structure (\*.struct) file needs to be provided if a user wants to compare the alignment score with the RNAMotifScanX alignment score. For example, the structure files for C-loop, E-loop, Kink-turn, reverse Kink-turn, and Sarcin-ricin motif families are provided in the [models](my_lib/RNAMotifScanX-release/models/) directory which was included as examples in RNAMotifScanX.

### 3. Usage commands

**Profile generation mode**
Generate profiles for the motif families provided in the input file:

```
python3 main.py -i <input_file> -o <output_subdirectory>
```

Generate profiles for the motif families provided in the input file by filtering extra gaps and using branch-and-bound (recommended):

```
python3 main.py -i <input_file> -o <output_subdirectory> -f -bnb
```

**Motif search mode**
Search similar structures in the motif instances of the input file against the provided profile data:

```
python3 main.py -i <input_file> -o <output_subdirectory> -f -bnb -m search -p <profile_file>
```

Search similar structures in the motif instances of the input file against the provided profile data and compare the alignment score with RNAMotifScanX using a manually generated consensus structure:

```
python3 main.py -i <input_file> -o <output_subdirectory> -f -bnb -m search -p <profile_file> -x -xs <consensus_file>
```

**Example**

```
python3 main.py -i ./sample1.in -o sample
python3 main.py -i ./sample1.in -o sample -f -bnb
python3 main.py -i ./known_families_IL.in -o known_families -f -bnb

python3 main.py -i ./sample1.in -o sample -f -bnb -m search -p ./output/sample/Kink-turn-Sub-1av2_20231117-125132.pfl
python3 main.py -i ./sample1.in -o sample -f -bnb -m search -p ./output/sample/Kink-turn-Sub-1av2_20231117-125132.pfl -x -xs ./my_lib/RNAMotifScanX-release/models/k-turn_consensus.struct
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

### ACKNOWLEDGEMENTS

RNAMotifProfile is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

For bug reports or comments please contact shaojie.zhang@ucf.edu
