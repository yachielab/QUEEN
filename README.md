# dna.py Installation and User Manual
dna.py is a Python module that enables users to handle annotated double-stranded DNA objects universally and design molecular cloning of new plasmids and DNA parts using existing annotated DNA files. In dna.py, all of the manipulations required in DNA engineering were covered by four simple operational functions: “cut,” “flip,” “end-modification,” and “join.” Using dna.py, a new plasmid construction process can be designed by programming a Python script or by an interactive interpreter under Jupyter Notebook. The designed DNA can be output as a GenBank format file including its building procedure and used to design other molecular cloning processes. 

## Software dependency
- python 3.8.0 or later
- biopython 1.7.8 later
- regex
- jupyter

## Installation
1.  Download the software by  ````git clone https://github.com/yachielab/dna.py````
2.  Install the following Python packages by  
	```
	pip install matplotlib
	pip install numpy
	pip install python-levenshtein
	pip install biopython
	pip install regex
	pip install jupyter
	```
3.  Set PYTHONPATH to the directory where you cloned the repository.

## Usage
dna.py provides a “DNA class” to define double-stranded DNA objects with annotations. The usage of DNA class objects and five functions to operate DNA class objects are described in the following sections. The following usage examples can be executed on the Jupyter Notebooks in “demo/tutorial.” 

### DNA class
A DNA class object defines double-stranded DNA with sequence annotations. It can be created by specifying a DNA sequence or by incorporating a sequence data file in GenBank or FASTA format.

**Example code 1: Create a blunt-end DNA object by specifying its sequence**  
```python
#Soruce code#
from dna import
brick = DNA(seq="CCGGTATGCGTCGA")
```

**Example code 2: Create a sticky-end** **DNA object by specifying its structure**  
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
- **.project**: *`str`*
  Name of `DNA` object (project). If a DNA object is created from a GenBank or FASTA format file, its sequence ID will be inherited. 

- **.seq**: *`str`*
  The top strand sequence of `DNA` object. If the bottom strand sequence represents a sticky end, the gap region of the top strand sequence will be complemented according to the bottom strand sequence. The double-stranded DNA structure of `DNA` object can be obtained by `.getdnaseq` described in the following section.

- **.topology**: *`str`* (`"linear"` or `"circular"`)
  Defining sequence topology of `DNA` object.

- **.dnafeatures**: *`list`* of *`DNAFeature`* objects
  List of *`DNAfeature`* objects. It provides a feature annotation for a certain range of sequence in the *`DNA`* object.  The following attributes are assigned for *`DNAfeature`* object. 
	- **.feature\_id:** *`str`*
	Unique identifier for *DNAfeature* object. It is automatically assigned to each feature when *DNA* object is created by loading a GenBank or FASTA format file.
	- **.feature\_type:** *`str`*
	Defining biological nature of *DNAfeature* object. 
  	- **.start:** *`int`*
	Start position of *DNAfeature*.
 	- **.end:** *`int`* 
	End position of *DNAfeature*.
	- **.span:** *`(int,int)`*
	Base span of *DNAfeature.* It is composed of (.start, .end).
	- **.strand:** *`int`* `(-1 or 1)`
	Direction of *DNAfeature*. 
	- **.qualifiers:** *`dict`* 
	Qualifiers of the feature. it reflects the qualifiers on a GenBank feature. The keys of the dictionary are qualifier names.
	- **.sequence:** *`str`*
	Sequence from .start to .end on .strand in the *DNA* object.

	*`DNAfeature`* object is implemented as subclass of Biopython *`Seqfeature`* object. Therefore, the other attributes and functions are totally inherited from *SeqFeature* object. For details, please see https://biopython.org/docs/dev/api/Bio.SeqFeature.html .

### Analytical functions
dna.py module provides the following print and search functions to analyze DNA class objects.

- **`.getdnaseq`**_`(region=list, display=bool, whole=bool, end_length=int, linebreak=int)`_  
	**Parameters**
	- **start**: *`int`*  (default: 0)
	Start position of the sequence. 
  	- **end**: *`int`*  (default: 0)
	Start position of the sequence. 
  	- **strandt**: *`int`* 1, 0, or -1 (default: 0)
	Start position of the sequence. 
  	- **display**: *`bool`* (`True` or `False`; default: `True`)  
	If `True`, the function will print double-stranded DNA sequence structure of the `DNA` object and return None.  
	if `display` is `False`, the following parameters will be ignored.  
  	- **whole**: *`bool`*(`True` or `False`; default: `True`)  
	If `False`, it will display the partial end sequence structures of `DNA` object of which lengths are given by `end_length`.
  	- **end\_length**: *`int`* (default: 10)
  	- **linebreak**: *`int`* or `None` (default: `None`)  
	Length of sequence for line break.

  	**Return**  
	if `strand` is 1 or -1, `str: DNA sequence from start to end on the strand ( 5' to 3')`.  
	if `strand` is `None` or `0`, `[str: top strand sequence from start to end (5' to 3'), str：bottom strand sequence from start to end (5' to 3')`.  
	if `display` is `True`,`None`.   
	
	**Example code 4: Print a double-strand DNA sequence with sticky ends**  	
	```python
	#Soruce codea#
	from dna import *
	brick = DNA(seq="CCGGTATGCG----/----ATACGCAGCT")
	brick.getdnaseq(display=True)
	#output#
	5' CCGGTATGCG---- 3'
	3' ----ATACGCAGCT 5'
	```
- **`.finddna`**_`(query=str, key_attribute=str, min_match=int, max_mismatch=int)`_  
	Search specific features that hold query values at given attributes from the *DNA* object. However, if the attribute is `"sequence:*"` and there are no features on the query sequences in the *`DNA`* object, features including its location and a digestion format (the details described later) about the query sequence will be returned.  
  	**Parameters**
  	- **query**: *`str`* (default: `None`)  
	Query for the search. It allows fuzzy matching with regular expression. For details, see https://pypi.org/project/regex/. If query is None, “.+” will applied for the query. 
	- **key\_attribute**: *`str`* (default: `None`)  
	Attribute type to search target features or sequence (`feature_id`, `feature_type`, `“qualifier:*”`, `strand` or `"sequence:*"`). If the argument is `None` or not given, it will be applied to all of the possible attributes excluding `sequence`. However, if `query` string is composed of characters expression DNA sequence clearly, `key_attribute` will be automatically set by `sequence`. If attribute is `sequence` and `query` is `None`, the entire sequence region will be target for the search. The partial sequences can be specified by using absolute positions in the entire sequence like `“sequence:|int..int|”`. To specify target strand for search, please use`"sequence:|int..int|+"` or `"sequence:|int..int|-"`.  
	As an advanced method, if attribute is sequence or `"sequence:*"`, users can specify cutting site in the query sequence using `"^"`, `"_"` and `"(*/*)"`. For example, EcoRI cutting site can be given by `"G^AATT_C"`. `"^"` indicates a cut position on the top strand. “\_” indicates a cut position on the bottom strand. Cutting site for TypeIIS restriction enzyme such as BsaI can be given by GGTCTC(1/5) that is same meaning with `"GGTCTCN^NNN_N"`. If any cutting site is not given in a query sequence, default cutting site will be end of the sequence. It is default for all *`DNAFeature`* object.   
 	- **min\_match**: *`int`* (default: length of query - `max_mismatch`)
	Minimal letters required to match to the query (only when the query is not provided with regular expression).  
	- **max\_mismatch**: *`int`* (default: 0)  
	Maximum number of letters allowed to mismatch to the query (only when the query is not provided with regular expression).
	 
	**Return**
	*`list`* (list of *`DNAFeature`* objects)
	
	**Example code 5: Search for a DNA sequence with a regular expression** 
	```python
	#Source code (continued from the previous example)# 
	feature_list = plasmid.finddna(key_attribute="sequence",query="[ATGC]{20}[ATGC]GG")
	for feature in feature_list:
	    print(feature.start, feature.end, plasmid.getdnaseq(feature.start, feature.end, feature.strand), sep="\t")
	#Output#
	88	111	ACGACTCAACAGCTTAACGTTGG
	137	160	ACTCTCACTCTTACCGAACTTGG
	200	223	AAACGAATCGACCGATTGTTAGG
	237	260	AGGAAGGTTTAAACGCATTTAGG
	267	290	ATAGAAGTGTGTATCGCTCGAGG
	293	316	CCGAATTCGAAGACTTGGTACGG
	327	350	TTTCCAGATCTGATAACTTGTGG
	376	399	TTAAGACGTCAGAATTCTCGAGG
	422	445	AGTGAGTCGTATTAATTTCGCGG
	︙
	```  

	**Example code 6: Search CDS features**
	```python
	#Source code (continued from the previous example)# 
	feature_list = plasmid.finddna(query="CDS")
	for feature in feature_list:
	    print(feature)
	#Output#
	type: CDS
	location: [546:1206](-)
	qualifiers:
	    Key: codon_start, Value: ['1']
	    Key: gene, Value: ['cat']
	    Key: label, Value: ['CmR']
	    Key: note, Value: ['confers resistance to chloramphenicol']
	    Key: product, Value: ['chloramphenicol acetyltransferase']
	    Key: translation, Value: ['MEKKITGYTTVDISQWHRKEHFEAFQSVAQCTYNQTVQLDITAFLKTVKKNKHKFYPAFIHILARLMNAHPEFRMAMKDGELVIWDSVHPCYTVFHEQTETFSSLWSEYHDDFRQFLHIYSQDVACYGENLAYFPKGFIENMFFVSANPWVSFTSFDLNVANMDNFFAPVFTMGKYYTQGDKVLMPLAIQVHHAVCDGFHVGRMLNELQQYCDEWQGGA']
	```  

- **`.printfeature`**_`(feature_list=None, attribute=list, detail=bool, separation=str, output=str, zero_based_index=bool)`_ 
	Print a tidy data table of annotation features/attributes in DNA object. Default printing attributes are `"feature ID"`, `"feature type"`, `"qualifier:label"`, `"start"`, `"end"`, and `"strand"`. Unique `"feature ID"` is automatically assigned to each feature when DNA object is created by loading a GenBank format file.
	**Parameters**
	- **­feature\_list**: *`list`* of *`DNAfeaure`* objects (default: `.dnafeatures`)
	Feature list displayed into the output table. if the argument value is `None` or not given, all features of the `DNA` object will be displayed in the output table. 
	- **attribute**: *`list`* of feature IDs (default: `["feature_id", "feature_type", "qualifier:label", "start position", "end position", "strand"]`)
	Selected attribute types for printing feature information. Default attributes can be specified by using `"$DEFAULT"` in the list. If detail is `True`, the value will be ignored.
	- **detail**: `bool` (`True` or `False`; default: `False`)
	If this is `True`, all attributes excluding `"sequence"` will be output to the table.
	- **seq**: bool (`True` or `False`; default: `False`)
	If this is `True`, sequence of each feature is output to the table.
	- **­separation**: `str` (default: `None`)
	String to separate each line values. If the argument value is `None`, a well-formatted table will be generated with multiple spaces.
	- **output**: `str` (default: `stdout`)
	Output file name or file object. If the argument value is stdout or None, the table will be output to stdout.
	- **zero\_based\_index**: *bool* (`True` or `False`; default: `True`)
	As a default, positions of features are given in zero-based indexing (same as Python indexing). If the argument value is `False`, 1-based indexing will be applied (as seen in the GenBank format).  
    	
	**Return**  
	`None`

	**Example code 7: Display all of the "primer\_bind" features in a GenBank file**
	```python
	#Source code (continued from the previous example)#
	feature_list = plasmid.finddna("primer_bind")
	plasmid.printfeature(feature_list)
	
	#Output#
	feature_id  qualifier:label                  feature_type  start  end  strand  
	100         Forward (CW) Analysis Primer     primer_bind   233    260  +       
	600         Cloning Analysis Reverse Primer  primer_bind   410    436  -         
	```

### Operational functions
dna.py provides the following five fundamental operational functions to manipulate DNA objects.  

- **`.cutdna`**_`dna_object=*DNA* object, *cutsites=**list* of int, "int/int," or  DNAFeature,  process_discription=”str”)`_
	Cut a DNA object at queried positions and return a list of linear DNA fragment objects. The cut positions can be specified differently for top strand and bottom strand sequences (0-based position). All respective features are inherited to the fragmented DNA objects. DNA features on the cut boundaries are also carried over to the fragments by specifying such information described later. If the fragments are ligated back, the DNA features will be also restored. It gives the cropped region information to split features on the cut boundary as `"qualifier:crop_trail"`. The qualifier is composed of the following contents.  
`[label of the original feature]:[project of the original feature]:[length of the original feature]:[start and end positions in the original feature]`
This function also allows linearization of a circular DNA object by having a single cut. As a branch function, `cropdna(dna_object=DNA object, start, end)` is also prepared. It returns a fragmented DNA object in sequence positions between `start` to `end`. 
	**Parameters**
	- **dna\_object**: *`DNA`* object
	- **positions**: *`list`* of *`int`* , *`"int/int"`*, or *`DNAFeature`* object 
	List of cut positions. If users want to cut at different positions on the top and bottom strands, the cut site can be set by the pairs of positions at the top and bottom strands. In that case, the first element of a position combination must be smaller than the second one of the next position combinations. Please refer to the following examples.  
	Right case : cutdna(object, \*[“50/55”, “100/105”])  
	Wrong case : cutdna(object, \*[“50/105”, “100/55”])  
	If cut positions are specified by *`DNAFeature`*  object, each cut position on top and bottom strand will be set by following the cutting format of the *`DNAFeaure`* object.
	
	**Return**
	*`list`* (list of *`DNAFeature`* objects)
	
	**Example code 8: Cut pGGA plasmid at multiple positions**  
	Cut a circular plasmid pGGA at two positions and generate two fragments. The one of the two fragments also cut at a position.  
	```python
	#Source code (continued from the previous example)# 
	fragment1, fragment2 = cutdna(plasmid ,1000, 2000)
	print(fragment1)
	print(fragment2)
	fragment3, fragment4 = cutdna(fragment1, 500)
	print(fragment3)
	print(fragment4)
	
	#Output#
	<dna.DNA object; project='pGGA', length='1000 bp', topology='linear'> 
	<dna.DNA object; project='pGGA', length='1174 bp', topology='linear'> 
	<dna.DNA object; project='pGGA', length='500 bp', topology='linear'> 
	<dna.DNA object; project='pGGA', length='674 bp', topology='linear'>
	```
  
	**Example code 9: Crop a fragmented dna object in a specific region**  
	To extract only fragment2, please use cropdna function as follows.  
	```python
	#Source code (continued from the previous example)# 
	fragment2 = cropdna(plasmid ,2000, 1000)
	print(fragment2)

	#Output#
	<dna.DNA object; project='pGGA', length='1174 bp', topology='linear'>
	```  
  	
	**Example code 10: Digest pGGA plasmid by EcoRI**  
	The recognition sequence where EcoRI cut is "5′-G^AATT\_C-3′". There are three EcoRI site in the pGGA plasmid. First, search the EcoRI sequence from the pUC19 plasmid. Then, cut the EcoRI site according to its cutting format.   
  	
	```python
	#Source code (continued from the previous example)#
	sites = plasmid.finddna("G^AATT_C", key_attribute="sequence")
	fragments = cutdna(plasmid, *sites)
	for fragment in fragments:
  	    print(fragment)
  	    fragment.getdnaseq(display=True, whole=False, end_length=10) 
	
	#Output#
	<dna.DNA object; project='pGGA', length='96 bp', topology='linear'>
	5' AATTCGAAGA...CGTCAG---- 3'
	3' ----GCTTCT...GCAGTCTTAA 5'

	<dna.DNA object; project='pGGA', length='604 bp', topology='linear'>
	5' AATTCTCGAG...ATACGG---- 3'
	3' ----GAGCTC...TATGCCTTAA 5'

	<dna.DNA object; project='pGGA', length='1486 bp', topology='linear'>
	5' AATTCCGGAT...GATCCG---- 3'
	3' ----GGCCTA...CTAGGCTTAA 5'
	```
