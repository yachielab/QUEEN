# dna.py Installation and User Manual

dna.py is a Python module that enables users to universally handle
annotated double-stranded DNA objects and design molecular cloning of
new plasmids and DNA parts using existing annotated DNA files. In
dna.py, all of the manipulations required in DNA engineering were
covered by four simple operational functions: "cut," "flip,"
"end-modification," and "join." Using dna.py, a new plasmid construction
process can be designed by programming a Python script or by an
interactive interpreter under Jupyter Notebook. The designed DNA can be
output as a GenBank format file and used to design other molecular
cloning processes.

## Software dependency
- python 3.8.0 or later
- biopython 1.7.8 later
- regex
- jupyterlab

## Installation
1.  Download the software by
```git clone https://github.com/yachielab/dna.py```
2.  Install the following Python packages by
```
pip install matplotlib
pip install biopython
pip install regex
pip install jupyterlab
```
3.  Set PYTHONPATH to the directory where you cloned the repository.

## Usage
dna.py provides a "DNA class" to define double-stranded DNA objects with annotations. The usage of DNA class objects and four functions to operate DNA class objects are described in the following sections. More usage examples are shown by Jupyter Notebook "example/tutorial/  tutorial_ex01_ex16.ipynb." & "example/tutorial/  tutorial_ex17_ex18.ipynb."

### DNA class
A DNA class object defines double-stranded DNA with sequence annotations. It can be created by specifying a DNA sequence or by incorporating a sequence data file in GenBank or FASTA format.

**Example code 1: Create a blunt-end DNA object by specifying its sequence**  
```python
#Soruce code#
from dna import
brick = DNA(seq="CCGGTATGCGTCGA")
```

**Example** **code 2: Create a sticky-end** **DNA object by specifying its structure**  
```python
#Soruce code#
from dna import *
brick = DNA(seq="CCGGTATGCG----/----ATACGCAGCT")
```
The two values separated by \"/\" show the top and bottom strand sequences of the generating DNA object, respectively. The top strand sequence is given in 5'-\>3' direction from left to right and the bottom strand sequence is given by 3'-\>5' direction from left to right. \"-\" indicates gap. A:T and G:C base pairing rule is required between the two strings except for the positions with gaps on one strand.  

**Example code 3: Create an annotated DNA object from GenBank format file**
```python
#Soruce code#
from dna import *
plasmid = DNA(record="pUC19.gbk")
```

#### Properties of DNA class objects
- **.project**:`str`
  Name of DNA object (project). If a DNA object is created from a GenBank or FASTA format file, its sequence ID will be inherited. 
  
- **.seq**:`str`
  The top strand sequence of`DNA`object. If the bottom strand sequence represents a sticky end, the gap region of the top strand sequence will be complemented according to the bottom strand sequence. The double-stranded DNA structure of `DNA`object can be obtained by`.getdnaseq`described in the following section.
  
- **.topology**:`str`(`"linear"` or `"circular"`)
  Defining sequence topology of`DNA`object.

- **.dnafeature**:`list`of`Bio.SeqFeature`objects
  List of Biopython SeqFeature (`Bio.SeqFeature`) objects. Each`Bio.SeqFeature`object provides a feature annotation for a a certain range of sequence in the`DNA`object. For details, please see https://biopython.org/docs/dev/api/Bio.SeqFeature.html.  

- **.record:**`Bio.SeqRecord`object
  A Biopython SeqRecord(`Bio.SeqRecord`) object that represents an optional description of the DNA object. For details, please see https://biopython.org/docs/dev/api/Bio.SeqRecord.html. 

- **.subject:**`None`or`DNA`object.
  If the DNA object is created as a part of the parental DNA object, the parental DNA object will be assigned to .subject.
  - **.start:**`int`
    Start position of the DNA in the parental DNA object.
  - **.end:**`int`
    End position of the DNA object in the parental DNA object.
  - **.strand:**`int`(-1,+1)
    Direction of the DNA object in the parental DNA object.
  - **.dna:**`DNA`object
    Parental DNA object.

### Analytical functions
dna.py module provides the following print and search functions to analyze DNA class objects.

- **`.getdnaseq`**_`(region=list, display=bool, whole=bool, end_length=int, linebreak=int)`_  
  
  ##### Parameters
  - **region**: `list` composed of [`start`, `end`, `strand`] or `Bio.SeqFeture` object; default: [0, length of the `DNA` object, 0]  
  The sequence region returned by the function. If the region specified by`Bio.SeqFeature` object, `start`, `end` and `strand` will be assigned by`Bio.SeqFeature.location.start`, `Bio.SeqFeature.location.end`, `Bio.SeqFeature.location.strand`, respectively. If `strand` is 0, the function will return the sequence of both strands. If `strand` is 1, it will return only the top strand sequence. If `strand` is -1, it will return only the bottom strand sequence.
  - **display**: `bool`(`True` or `False`; default: `True`)  
    If `True`, the function will print double-stranded DNA sequence structure of the `DNA` object and return None.if False, the following parameters will be ignored.
  - **whole**: `bool`(`True` or `False`; default: `True`)  
  If `False`, it will display the partial end sequence structures of `DNA` object of which lengths are given by `end_length`.
  - **end_length**: `int` (default: 10)
  - **linebreak**:`int` or `None` (default: `None`)  
   Length of sequence for line break.
  
  ##### Return
  `str`： DNA sequence in the specified region (from 5' to 3') , `[str:`*`top strand sequence (from 5' to 3')`*,`str：` *`bottom  strand sequence (from 5\' to 3\')`*\] or `None`
  
  ##### Example code 4: Print a double-strand DNA sequence with sticky ends
  ```python
  #Soruce codea#
  from dna import *
  brick = DNA(seq="CCGGTATGCG----/----ATACGCAGCT")
  brick.printdnaseq()
  ```
  ```
  #output#
  5' CCGGTATGCG---- 3'
  3' ----ATACGCAGCT 5'
  ```

- **`.finddna`**_`(query=str, attribute=str, min_match=int, max_mismatch=int)`_  
  **Parameters**
  
  - **query**: `str` (default: `None`)
    Query for the search. It allows fuzzy matching with regular expression. For details, see https://pypi.org/project/regex/. If query is None, “.+” will applied for the query.
  - **attribute**: `str` (default: `None`). 
    Attribute type to search target features or sequence (`"feature_id"`, `"feature_type"`, `“qualifier:*”`, `"strand"` or `“sequence:*”`). If the argument is `None` or not given, it will be applied to all of the possible attributes excluding sequence.  However, if query string is composed of standard [IUB/IUPAC nucleic acid codes](https://en.wikipedia.org/wiki/Nucleic_acid), attribute is set as “sequence”. If `attribute` is `"sequence"` and `query` is `None`, the entire sequence region will be target for the search. The partial sequences can be specified by using absolute positions like `“sequence:|`*`int`*`..`*`int`*`|”`. If you also want to specify the strand, please use  `“sequence:|`*`int`*`..`*`int`*`|+”` or  `“sequence:|`*`int`*`..`*`int`*`|-”`. If the attribute is not `"sequence"`, the function returns the existing features in the *`DNA`* object retrieved based on query and attribute. However, if the attribute is `"sequence:*"` and the features with the query sequence do not exist in the *`DNA`* object, it will return empty *`SeqFeature`* objects that include only positional information about the sequence  searched by the query. 
  - **min_match**: `int` (default: length of query)
    Minimal letters required to match to the query (only when the query is not provided with regular expression).
  - **max_mismatch**: `int` (default: 0)
  Maximum number of letters allowed to mismatch to the query (only when the query is not provided with regular expression).
  

  **Return**
    *`list`* (list of *`SeqFeature`* objects)

  ##### Example code 5: Search for a DNA sequence with a regular expression
  ```python
  #Source code (continued from the previous example)# 
  feature_list = plasmid.finddna(attribute="sequence",query="[ATGC]{20}[ATGC]GG")
for feature in feature_list:
print(feature.location.start, feature.location.end, feature.location.strand, plasmid.getdnaseq(feature), sep="\t")
  ```
  ```
  #Output#
  25	48	1	GAAAGCGCCACGCTTCCCGAAGG
55	78	1	GCGGACAGGTATCCGGTAAGCGG
83	106	1	TCGGAACAGGAGAGCGCACGAGG
107	130	1	AGCTTCCAGGGGGAAACGCCTGG
175	198	1	CGATTTTTGTGATGCTCGTCAGG
210	233	1	ATGGAAAAACGCCAGCAACGCGG
239	262	1	TACGGTTCCTGGCCTTTTGCTGG
︙
  ```

  ##### Example code 6: Search CDS features

  ````python
  #Source code (continued from the previous example)# 
  feature_list = plasmid.finddna(query="CDS")
  for feature in feature_list:
      print(feature)
  ````

  ```
  #Output#
  type: CDS
  location: [614:938](+)
  qualifiers:
      Key: codon_start, Value: ['1']
      Key: gene, Value: ['lacZ fragment']
      Key: label, Value: ['lacZ-alpha']
      Key: product, Value: ['LacZ-alpha fragment of beta-galactosidase']
      Key: translation, Value: ['MTMITPSLHACRSTLEDPRVPSSNSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRLMRYFLLTHLCGISHRIWCTLSTICSDAA']
  
  type: CDS
  location: [1283:2144](+)
  qualifiers:
      Key: codon_start, Value: ['1']
      Key: gene, Value: ['bla']
      Key: label, Value: ['AmpR']
      Key: note, Value: ['confers resistance to ampicillin, carbenicillin, and related antibiotics']
      Key: product, Value: ['beta-lactamase']
      Key: translation, Value: ['MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW']
  ```
- **`.printfeature`**_`(feature_list=None, attribute=list, detail=bool, separation=str, output=str, zero_based_index=bool)`_  
  Print a tidy data table of annotation features/attributes in DNA object. Default printing attributes are `"feature ID"`, `"feature type"`, `"qualifier:label"`, `"start"`, `"end"`, and `"strand"`. Unique `"feature ID"` is automatically assigned to each feature when DNA object is created by loading a GenBank format file.

  ##### Parameters
  - **­feature_list**: `list` (default: `None`)
    Feature list displayed into the output table. if the argument value is `None` or not given, all features of the `DNA` object will be displayed in the output table. 
  - **attribute**: `list` of feature IDs (default: `["feature_ID", "feature_type", "qualifier:label", "start position", "end position", "strand"]`)
    Selected attribute types for printing feature information. Default attributes can be specified by using `"$DEFAULT"` in the list. If detail is `True`, the value will be ignored.
  - **detail**: `bool` (`True` or `False`; default: `False`)
    If this is `True`, all attributes excluding `"sequence"` will be output to the table.
  - **­separation**: `str` (default: `None`)
    String to separate each line values. If the argument value is None, a well-formatted table will be generated with multiple spaces.
  - **output**: `str` (default: `stdout`)
    Output file name or file object. If the argument value is stdout or None, the table will be output to stdout.
  - **zero_based_index**: *bool* (`True` or `False`; default: `True`)
    As a default, positions of features are given in zero-based indexing (same as Python indexing). If the argument value is `False`, 1-based indexing will be applied (as seen in the GenBank format).
  
  **Return**
    `None`
  
  ##### Example code 7: Display all of the "primer_bind" features in a GenBank file
  
  ```python
  #Source code#
  features = plasmid.findna(“primer_bind”)
  brick.printfeature(featurs, attribute=["$DEFAULT", "sequence"])
  
  ```
  ```
   #Output#
  feature_ID  label            type         start  end   strand  sequence                      
  100         pBR322ori-F      primer_bind  117    137   +       GGGAAACGCCTGGTATCTTT     
  200         L4440            primer_bind  370    388   +       AGCGAGTCAGTGAGCGAG       
  600         M13/pUC Reverse  primer_bind  583    606   +       AGCGGATAACAATTTCACACAGG  
  700         M13 rev          primer_bind  602    619   +       CAGGAAACAGCTATGAC        
  800         M13 Reverse      primer_bind  602    619   +       CAGGAAACAGCTATGAC        
  1100        M13 Forward      primer_bind  688    706   -       TGTAAAACGACGGCCAGT       
  1200        M13 fwd          primer_bind  688    705   -       GTAAAACGACGGCCAGT        
  1300        M13/pUC Forward  primer_bind  697    720   -       CCCAGTCACGACGTTGTAAAACG  
  1400        pRS-marker       primer_bind  913    933   -       CGGCATCAGAGCAGATTGTA     
  1500        pGEX 3'          primer_bind  1032   1055  +       CCGGGAGCTGCATGTGTCAGAGG  
  1600        pBRforEco        primer_bind  1092   1111  -       AATAGGCGTATCACGAGGC      
  1900        Amp-R            primer_bind  1501   1521  -       ATAATACCGCGCCACATAGC     
  ```



### Operational functions

dna.py provides the following five fundamental operational functions to manipulate DNA objects.