# dnaquine.py Installation and User Manual
dnaquine.py is a Python module that enables universal editing and annotation of double-stranded DNA objects for designing molecular cloning and simulating genome editing, DNA recombination and genetic circuits based on such. DNA parts information can be imported from external annotated DNA files. Output .q.gb (“quine” GenBank format) file encodes the full information of the constructed DNA and a quine code that can self-reproduce the file itself. In dnaquine.py, all of the manipulations required in DNA engineering are covered by five operational functions “cut,” “end-modification,” “flip,” “join,” and more like a super function “edit” along with various analytical and visualization functions. A new plasmid can be designed by programming a Python script or by an interactive interpreter under Jupyter Notebook. The designed DNA can be output as a q.gb that encode a quine code to reproduce the file itself. Therefore, designing DNA using dnaquine.py autonomously generate a process description for its perfect reproduction.

## Software dependency
- python 3.8.0 or later

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
dnaquine.py provides a “DNA class” to define double-stranded DNA objects with annotations. The DNA class and the usage of its functions to operate DNA class objects are described below. All of the provided examples can be executed in the Jupyter Notebooks in “demo/tutorial.”

### DNA class
A DNA class defines double-stranded DNA objects with sequence annotations. It can be created by specifying a DNA sequence or incorporating a sequence data file in GenBank or FASTA file format. 

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
The two values separated by "/" show the top and bottom strand sequences of the generating DNA object, respectively. The top strand sequence is given in 5’→3’ direction from left to right and the bottom strand sequence is given by 3’→5’ direction from left to right. "-" indicates gap. A:T and G:C base pairing rule is required between the two strings except for the positions with gaps.

**Example code 3: Create an annotated DNA object from GenBank format file**
```python
#Soruce code#
from dna import *
plasmid = DNA(record="pGGA.gbk")
```

#### Properties of DNA class objects
- **.project**: *`str`*  
  Name of `DNA` object (project). If a DNA object is created from a GenBank or FASTA format file, its sequence ID will be inherited. 

- **.seq**: *`str`*  
  The top strand sequence of DNA object. Sticky end gaps on the top strand if any are complemented according to the bottom strand sequence. The double-stranded DNA structure of the DNA object can be obtained by `.getdnaseq` described below.

- **.topology**: *`str`* (`"linear"` or `"circular"`)  
  Defining sequence topology of `DNA` object.

- **.dnafeatures**: *`list`* of *`DNAFeature`* objects
  List of *`DNAfeature`* objects. It provides a feature annotation for a certain range of sequence in the *`DNA`* object.  The following attributes are assigned for *`DNAfeature`* object. 
	- **.feature\_id:** *`str`*  
	Unique identifier for *DNAfeature* object. It is automatically assigned to each feature when DNA object is created by loading a GenBank or FASTA format file.
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
	Qualifiers of the feature. When GenBank file is imported, the qualifiers of the features are imported here. The keys of the dictionary are qualifier names.
	- **.sequence:** *`str`*
	Sequence from .start to .end on .strand in the *DNA* object.

	*`DNAfeature`* class is implemented as a subclass of Biopython *`SeqFeature`* class, enabling dnaquine to access all the instance variables and methods of *`SeqFeature`*. For details, see https://biopython.org/docs/dev/api/Bio.SeqFeature.html

### Output functions
dnaquine.py provides functions to output information of DNA class object.

- **.getdnaseq**_`(region=list, display=bool, whole=bool, hide_middl=int, linebreak=int)`_  
	Returns and displays the sequence of a specified region.
	**Parameters**
	- **start**: *`int`*  (default: 0)  
	Start position of the sequence.   
  	- **end**: *`int`*  (default: 0)  
	Start position of the sequence.   
  	- **strand**: *`int`* 1, -1 (default: 2)  
	Start position of the sequence.   
  	- **display**: *`bool`* (`True` or `False`; default: `True`)  
	If `True`, the function will print double-stranded DNA sequence structure of the *`DNA`* object. 
	If `False`, the following parameters will be ignored. 
	- **hide\_middle**: *`int`* or `None` (default: `None`)  
  	Length of the end sequences to be displayed
	- **linebreak**: *`int`* or `None` (default: `None`)  
	Length of sequence for line break.  

  	**Return**  
	if `strand` is 1 or -1, `str: DNA sequence from start to end on the strand ( 5' to 3')`.  
	if `strand` is `None` or `2`, `[str: top strand sequence from start to end (5' to 3'), str：bottom strand sequence from start to end (5' to 3')`.  
	
	**Example code 4: Print a double-strand DNA sequence with sticky ends**  	
	```python
	#Soruce codea#
	from dnaquine import *
	fragment = DNA(seq="CCGGTATGCG----/----ATACGCAGCT") 
	fragment.getdnaseq(display=True)

	#output#
	5' CCGGTATGCG---- 3'
	3' ----ATACGCAGCT 5'
	```

- **.printfeature**_`(feature_list=None, attribute=list, detail=bool, separation=str, output=str, x_based_index=bool)`_ 
	Print a tidy data table of annotation features/attributes in DNA object. Default printing attributes are `"feature ID"`, `"feature type"`, `"qualifier:label"`, `"start"`, `"end"`, and `"strand"`. Unique `"feature ID"` is automatically assigned to each feature when DNA object is created by loading a GenBank format file.
	**Parameters**
	- **­feature\_list**: *`list`* of *`DNAfeaure`* objects (default: `.dnafeatures`)  
	- Feature list displayed into the output table. if the argument value is `None` or not given, all features of the *`DNA`* object will be displayed in the output table.  
	- **attribute**: *`list`* of feature IDs (default: `["feature_id", "feature_type", "qualifier:label", "start position", "end position", "strand"]`)  
	Selected attribute types for printing feature information. If attribute is `None`, it will make a table for all attributes except for `"sequence"`. `"sequence"` can be given to the attribute list.  
	- **seq**: bool (`True` or `False`; default: `False`)  
	If this is `True`, sequence of each feature is output to the table.  
	- **­separation**: `str` (default: `None`)  
	String to separate each line values. If the argument value is `None`, a well-formatted table will be generated with multiple spaces.  
	- **output**: `str` (default: `stdout`)  
	Output file name or file object. If the argument value is stdout or None, the table will be output to stdout.  
	- **x\_based\_index**: *bool* (`0` or `1`; default: `0`)  
	As a default, positions of features are given in zero-based indexing (same as Python indexing). If the argument value is `False`, 1-based indexing will be applied (as seen in the GenBank format).  
    	
	**Return**  
	`None`

	**Example code 5: Print DNAfeatures with formatted table**
	```python
	#Source code#
	from dnaquine import *
	plasmid = DNA(record="pGGA.gb")
	plasmid.printfeature()

	#Output#
	feature_id  qualifier:label                  feature_type  start  end   strand
	0           null                             source        0      2174  +
	1100        pGGA                             source        0      2174  +
	100         Forward (CW) Analysis Primer     primer_bind   233    260   +
	200         SP6 promoter                     promoter      253    272   +
	300         upstream MCS                     misc_feature  282    304   +
	400         BsaI insert                      misc_feature  314    359   +
	500         downstream MCS                   misc_feature  372    406   +
	600         Cloning Analysis Reverse Primer  primer_bind   410    436   -
	700         T7 promoter                      promoter      417    436   -
	800         CmR                              CDS           546    1206  -
	900         cat promoter                     promoter      1206   1309  -
	1000        ori                              rep_origin    1409   1998  +
	```
- **.outputdna**_`(output=str)`_  
	Output DNA object with GenBank format file. When the GenBank output is generated, A feature covering the entire sequence will be crated, and the history of the executed operational functions to construct the DNA objects will be described as a qualifier of the feature.
	**Parameters**
	- **output**: *`str`* (default: `None`)  
	Output file name. If the argument value is `None`, the table will be output to stdout. 

-----

- **.finddna**_`(query=str, key_attribute=str, min_match=int, max_mismatch=int)`_  
	Search specific features that hold query values at given attributes from the *DNA* object. However, if the attribute is `"sequence:*"` and there are no features on the query sequences in the *`DNA`* object, features including its location and a digestion format (the details described later) about the query sequence will be returned.  
  	**Parameters**
  	- **query**: *`str`* (default: `None`)  
	Query for the search. It allows fuzzy matching with regular expression. For details, see https://pypi.org/project/regex/. If query is None, “.+” will applied for the query. 
	- **key\_attribute**: *`str`* (default: `None`)  
	Attribute type to search target features or sequence (`"feature_id"`, `"feature_type"`, `"qualifier:*"`, `"strand"` or `"sequence:*"`). If the argument is `None` or not given, it will be applied to all of the possible attributes excluding `"sequence"`. However, if `query` string is composed of characters expression DNA sequence clearly, `key_attribute` will be automatically set by `"sequence"`. If the attribute is `"sequence"` and `query` is `None`, the entire sequence region will be target for the search. The partial sequences can be specified by using absolute positions in the entire sequence like `“sequence:|int..int|”`. To specify target strand for search, please use`"sequence:|int..int|+"` or `"sequence:|int..int|-"`.  
	As an advanced method, if attribute is sequence or `"sequence:*"`, users can specify cutting site in the query sequence using `"^"`, `"_"` and `"(*/*)"`. For example, EcoRI cutting site can be given by `"G^AATT_C"`. `"^"` indicates a cut position on the top strand. “\_” indicates a cut position on the bottom strand. Cutting site for TypeIIS restriction enzyme such as BsaI can be given by `"GGTCTC(1/5)"` that is same meaning with `"GGTCTCN^NNN_N"`. If any cutting site is not given in a query sequence, default cutting site will be end of the sequence. It is default for all *`DNAFeature`* object.   
 	- **min\_match**: *`int`* (default: `len(query) - max_mismatch`)
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

### Operational functions
dna.py provides the following five fundamental operational functions to manipulate DNA objects.  

- **.cutdna**_`(dna_object=*DNA* object, *cutsites=**list* of int, "int/int," or  DNAFeature,  process_discription=”str”)`_
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
	  
	**Example code 11: Digest pGGA plasmid by BsaI**
	The recognition sequence where BsaI cut is "5′-GGTCTC(N1)/(N5)-3′". There are two BsaI in the pGGA plasmid. First, search the BsaI sequence from the pUC19 plasmid. Then, cut the BsaI site according to its cutting format.  
	```python
	#Source code (continued from the previous example)#
	sites = plasmid.finddna("GGTCTC(1/5)", key_attribute="sequence")
	plasmid.printfeature(sites)
	fragments = cutdna(plasmid,*sites)
	for fragment in fragments:
	    print(fragment)
	    fragment.getdnaseq(display=True, whole=False, end_length=10)

	#Output#
	feature_id  qualifier:label  feature_type  start  end  strand  
	null        null             misc_feature  348    359  +       
	null        null             misc_feature  314    325  -       

	<dna.DNA object; project='pGGA', length='45 bp', sequence='GGAGCGAGACCGCTTTCCAGATCTGATAACTTGTGGTCTCACCAT', topology='linear'>
	5' GGAGCGAGAC...GTCTCA---- 3'
	3' ----GCTCTG...CAGAGTGGTA 5'

	<dna.DNA object; project='pGGA', length='2137 bp', topology='linear'>
	5' CCATTCCTGT...TGGTAC---- 3'
	3' ----AGGACA...ACCATGCCTC 5'
	```
  
- **modifyends**_`(dna_object=DNA object, left ="str/str", right=”str/str”, process_description="str")`_
	Set end sequence structures of a given DNA object.  
	**Parameters**
	- **dna_object**: *`DNA`* object 
	- **left**: *`"str/str"`* (default: “/”)
	Left end sequence structure of DNA object. The following examples show how to describe the parameter.  
	- **right**: *`"str/str"`* (default: “/”) 
	Right end sequence structure of DNA object. The following examples show how to describe the parameter.  
	
	**Return**  
	*`DNA`* object  
	
	**Example code 12: Trim single-stranded DNA on both ends to generate sticky ends**  
	Sticky ends can be generated by trimming single-stranded DNA sequences when their end structures are given by top and bottom strand strings with "*" and "-" separated by "/." The letters "-" show nucleotide letters being trimmed and “*” shod nucleotide letters being remained.  
	```python
	#Source code (continued from the previous example)#
	fragment = cropdna(plasmid, 100, 120)
	fragment.getdnaseq(display=True)
	fragment = modifyends(fragment, "----/****", "**/--")
	fragment.getdnaseq(display=True)

	#Output#
	5' CTTAACGTTGGCTTGCCACG 3'
	3' GAATTGCAACCGAACGGTGC 5'

	5' ----ACGTTGGCTTGCCACG 3'
	3' GAATTGCAACCGAACGGT-- 5'
	```

	The following code also can execute same process with the above one.   
	```python
	#Source code (continued from the previous example)#
	fragment = cropdna(plasmid, (104,100), (120,118))
	fragment.getdnaseq(display=True)

	#Output#
	5' ----ACGTTGGCTTGCCACG 3'
	3' GAATTGCAACCGAACGGT-- 5'
	```

	A regex-like format can be used for the end stracture specification.  
	```python
	#Source code (continued from the previous example)#
	fragment = cropdna(plasmid, 100, 120)
	fragment = modifyends(fragment, "-{5}/*{5}","*{5}/-{5}")
	fragment.getdnaseq(display=True)

	#Output#
	5' -----CGTTGGCTTGCCACG 3'
	3' GAATTGCAACCGAAC----- 5'
	```

	If an invalid end sequence structure is given, error message will be returned.  
	```python
	#Source code (continued from the previous example)#
	fragment = cropdna(plasmid, 100, 120)
	fragment = modifyends(fragment, "******/****--", "----**/******)

	#Error message#
	TypeError: Please sepcify a proper sequence pattern for the 'left' argument
	```
  
- **flipdna**_`(dna_object=DNA object, process_description=”str”)`_
	Invert a DNA object.  
	**Parameters**
	- **dnaobject**: *`DNA`* object
						
	**Return**  
	*`DNA`* object  
  
- **joindna** _`(*dna_objects=list of DNA objects, topology=str, process_description="str")`_
	Join multiple DNA objects and return an assembled DNA object. Sticky ends need to be compatible. If the sequences of the split features are restored by the joining process, the original feature will be restored. If the linear DNA object is processed by the function, a circularized DNA object will be created by joining both ends. 
	**Parameters**
	- **dna_objects**: *`list`* of *`DNA`* objects
	- **topology**: *`str`* (`"linear"` or `"circular"`; default: `"linear"`) 
	Topology of returned DNA object.  

	**Return**  
	*`DNA`* object  
	  
	**Example code 14:  Flip chloramphenicol resistant gene in pGGA plasmid** 
	```python
	#Source code (continued from the previous example)#
	site = plasmid.finddna("CmR")[0]
	fragment1, fragment2 = cutdna(plasmid, site.location.start.position, site.location.end.position)
	print(fragment1, fragment2)
	fragment2   = flipdna(fragment2)
	new_plasmid = joindna(fragment1, fragment2, topology="circular")
	plasmid.printfeature(plasmid.finddna("CmR"))
	new_plasmid.printfeature(new_plasmid.finddna("CmR"))

	#Output 
	<dna.DNA object; project='pGGA', length='660 bp', topology='linear'> <dna.DNA object; project='pGGA', length='1514 bp', topology='linear'>
	feature_id  qualifier:label  feature_type  start  end   strand  
	800         CmR              CDS           546    1206  -       

	feature_id  qualifier:label  feature_type  start  end  strand  
	100         CmR              CDS           0      660  -       
	```
  
- **editdna**_`(dna_object=DNA object, key_attribute=str, query=reg, min_match=int, max_mismatch=int, target_attribute=str, operation=function, new_copy=bool, process_description="str")`_
	Edit given attribute values of DNAfeature objects according to the one of three type operations:  removeattribute, replaceattribute, createattribute. editdna is parental function of finddna, if operation is None or not given, editdna works just like finddna.
	**Parameters**
	- **dna_objects**: *`list`* of *`DNA`* objects
	- **key_attribute**: *`str`* (default: `None`)
	Same parameter with key_attribute of `finddna()`.  
	- **query**: *`reg`* (regular expression; default: `None`)
	Same parameter with `query` of `finddna()`.  
	- **min_match**: *`int`* (default: `len(query) – max_mismatch`)
	Same parameter with `min_match` of `finddna()`.
	- **max_mismatch**:  *`int`* (default: 0) 
	Same parameter with `max_mismatch` of `finddna()`.  
	- **target_attribute**: *`str`* (default: `None`)  
	Attribute type of the target features to be operated (`"feature_id"`, `"feature_type"`, `"qualifier:*"`, `"sequence:*"`, `"strand"` or `"span"`). If the argument value is `None` or not given, it will be applied to all of the possible attributes of the target features except `"sequence"`. When the attribute type is given by `"sequence"`, the entire sequence region associated with each target feature will be an operational target. When the attribute type is given by `"sequence:!int..int!+"`, `"sequence:!int..int!-"` or `"sequence:!int..int!"`,  a partial sequence region associated with each target feature given by two relative positions (0-based position) will be the operational target (`"+"` and `"-"` for top and bottom strand sequences, respectively, If `"+"` or `"-"` is not given, the feature sequence strand will be a target). 
	- **operation**: *`removeattribute()`*, *`createattribute(value=“str”)`* or *`replaceattribute(oldvalue, newvalue)`*  (default: `None`)
	If `operation` value is `None` or not given, no operation will happen, and it will return a list composed of the searched *`DNAfeature`* objects.   
	**`removeattribute()`** will work if `target_attribute` is `"feature_id"`, `"feature_type"`, `"qualifier:*"`, `"sequence:*"`. `removeattribute()` removes `target_attribute` value. If `target_attribute` is `"feature_id"` and the function is `removeattribute()`, the entire feature will be removed.   
	**`createattribute(“str”)`** will work if the `target_attribute` is `"feature_id"`, `"feature_type"`, or `"qualifier:*"`. If `target_attribute` is `"feature_id"`, a new feature whose `feature_id` is `"str"` and `feature_type` is `"misc_feature"` will be created on the queried region. If `target_attribute` is `"feature_type"`, a new feature with an automatically generated ID and whose feature type is `"str"` will be created. If "target_attribute" is `"qualifier:*"`, the qualifier whose value is `"str"` will be added to the feature.  
	**`replaceattribute(oldvalue, newvalue)`** works for all `target_attribute` values except `"start"` and `"end"`. The operation replaces *`oldvalue`* in `target_attribute` with *`newvalue`*. When the entire `target_attribute` value is replaced with a new one, the *`oldvalue`* argument can be omitted. If `target_attribute` is feature_id, the replacement will be accepted unless there is no conflict with other existing `feature_id` values in the object. To move *`DNAfeature`* location, `target_attribute` is assigned by `"span"` and `newvalue` should be set by `(start, end)`, and *`oldvalue`* can be omitted. The replacement of `span` value moves *`DNAfeature`*. It has no effect for the sequence of the *`DNA`* object. Replacement of `strand` value also only flips the feature direction. To move a feature with its sequence, please use the other operational functions such as `cropdna`, `joindna`, and `flipdna`. If the `target_attribute` is `"sequence:*"`, the substitution operation for the sequence can be executed with `replaceattribute(oldvalue=reg, newvalue="str")` that replaces the region matched to reg (regular expression) with `"str"`. For any DNA sequence editing that confers change in sequence length, DNA coordinates of all affected features will be adjusted.   
	- **new_copy**: *`bool`* (default: `True`)
	If True, the method will return an edited DNA object as different object with original one. Otherwise, original DNA objects will be edited and return None.  
	  	
	**Return**  
	If `operation` is `False`, list of *`DNAFeature`* objects.  
	If `new_copy` is `False`, None.  
	Otherwise, DNA object
	  
	**Example code 15: Convert feature type from ‘CDS’ to ‘gene’**
	```python
	#Source code (continued from the previous example)#
	plasmid.printfeature(plasmid.finddna("CDS"))
	new_plasmid = editdna(plasmid, key_attribute="feature type", query="CDS", target_attribute="feature type", operation=replaceattribute("gene"))
	new_plasmid.printfeature(new_plasmid.finddna("gene"))

	#Output#
	feature_id  qualifier:label  feature_type  start  end   strand
	800         CmR              CDS           546    1206  -

	feature_id  qualifier:label  feature_type  start  end   strand
	800         CmR              gene          546    1206  -
	```

	**Example code 16: Break start codon of the CDS features**  
	Replace ATG start codon on CDS region to GTG. The editing targeted only the first 3 nucleotides on each CDS region.  
	```python
	#Source code (continued from the previous example)#
	print("CmR:", plasmid.getdnaseq(plasmid.finddna("CmR")[0])[0:10])
	new_plasmid = editdna(plasmid, key_attribute="feature type", query="CDS", target_attribute="sequence:!0..3!", operation=replaceattribute("ATG", "GTG"))
	print("CmR:", new_plasmid.getdnaseq(plasmid.finddna("CmR")[0])[0:10])

	#Output#
	CmR: ATGGAGAAAA
	CmR: GTGGAGAAAA
	```
	
	**Example code 17: Create new features where SpCas9 bind**  
	```python
	#Source code (continued from the previous example)#
	new_plasmid = editdna(plasmid, key_attribute="sequence", query="[ATGC]{20}[ATGC]GG", target_attribute="feature id", operation=createattribute("spCas9_target*"))
	new_plasmid = editdna(new_plasmid, key_attribute="feature id", query="spCas9_target[0-9]+", target_attribute="feature type", operation=replaceattribute("misc_bind"))
	features = new_plasmid.finddna("misc_bind")
	new_plasmid.printfeature(features, seq=True)

	#Output#
	feature_id        qualifier:label  feature_type  start  end   strand  sequence
	spCas9_target1    null             misc_bind     137    160   +       ACTCTCACTCTTACCGAACTTGG
	spCas9_target2    null             misc_bind     200    223   +       AAACGAATCGACCGATTGTTAGG
	spCas9_target3    null             misc_bind     237    260   +       AGGAAGGTTTAAACGCATTTAGG
	```

### History functions  
DNA class manages all editing procedures executed for each DNA object generated in a script. Each editing process by an operational function is inherited as a python command into the newly generated DNA object. Additionally, users can give `process_description="str"` arguments to all operational function to leave experimental notes relating to the operation.  

- **exporthistory**_`(dna_object=DNA object, ouput=str, description_only=bool)`_
	Export executable python script that generates the DNA object. 
	**Parameters**
	- **dna_object**: *`DNA`* object. 
	- **output**: *`str`*
	Output file name.
	- **description_only**: *`bool`*
	If True, python scripts will be omitted and output only the contents of 'process_description' arguments given to each operational function.  
	
	**Return**  
	None

**Example code 18: Simulate the molecular cloning process of pCMV-Target-AID**  
Simulate the construction process of pCMV-Target-AID plasmid. The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two fragments encoding the N- and C-terminus halves of Target-AID, which were both amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129 and HM128/RS046, respectively, with a backbone fragment amplified from pCMV-ABE7.10 using RS047/RS048.  

```python
#Source code#
from dna import * 
#Source code# 
#Read GenBank inputs
pCMV_ABE             = DNA(record="input/addgene_102919.gbk",project="pCMV_ABE")
pCMV_nCas_PmCDA1_ugi = DNA(record="input/addgene_79620.gbk", project="pCMV-nCas-PmCDA1-ugi")

#Simulate PCR to amlify N-terminus of Target-AID
description1 = "The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT)\
 (Addgene 79620) using primer pairs RS045/HM129"
RS045 = DNA(record="input/RS045.fasta", project="RS045", process_description=description1) 
HM129 = DNA(record="input/HM129.fasta", project="HM129", process_description=description1) 
FW    = pCMV_nCas_PmCDA1_ugi.finddna(RS045.seq[-15:],key_attribute="sequence") #Search primer binding region
RV    = pCMV_nCas_PmCDA1_ugi.finddna(HM129.seq[-15:],key_attribute="sequence") #Search primer binding region
frag1 = joindna(RS045, cropdna(pCMV_nCas_PmCDA1_ugi,*FW,*RV), flipdna(HM129), project="N-term_Target-AID", process_description=description1) #Simulate PCR

#Simulate PCR to amlify C-terminus of Target-AID
description2 = "The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT)\
 (Addgene 79620) using primer pairs HM128/RS046"
HM128 = DNA(record="input/HM128.fasta", project="HM128", process_description=description2) 
RS046 = DNA(record="input/RS046.fasta", project="RS046", process_description=description2)
FW    = pCMV_nCas_PmCDA1_ugi.finddna(HM128.seq[-15:],key_attribute="sequence")
RV    = pCMV_nCas_PmCDA1_ugi.finddna(RS046.seq[-15:],key_attribute="sequence")
frag2 = joindna(HM128, cropdna(pCMV_nCas_PmCDA1_ugi,*FW,*RV), flipdna(RS046), project="C-term_Target-AID", process_description=description2)

#Simulate PCR to amplify a backbone fragment
description3 = "The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048"
RS047 = DNA(record="input/RS047.fasta", project="RS047", process_description=description3) 
RS048 = DNA(record="input/RS048.fasta", project="RS048", process_description=description3)
FW    = pCMV_ABE.finddna(RS047.seq[-15:],key_attribute="sequence") 
RV    = pCMV_ABE.finddna(RS048.seq[-15:],key_attribute="sequence")
frag3 = joindna(RS047, cropdna(pCMV_ABE,*FW,*RV), flipdna(RS048), project="backbone", process_description=description3)

#Gibson Assembly
description4 = "The Target-AID plasmid (pCMV-Target-AID) was constructed\
 by assembling two insert fragments and a backbone fragments."
frag1 = modifyends(frag1, "*{25}/-{25}","-{28}/*{28}", process_description=description4)
frag2 = modifyends(frag2, "*{28}/-{28}","-{25}/*{25}", process_description=description4)
frag3 = modifyends(frag3, "*{25}/-{25}","-{25}/*{25}", process_description=description4) 
pCMV_Target_AID = joindna(frag1,frag2,frag3,topology="circular", project="pCMV-Target-AID", process_description=description4)

pCMV_Target_AID.writedna("pCMV-Target-AID.gbk")
pCMV_Target_AID.printfeature()

#Output#
feature_id  qualifier:label     feature_type  start  end   strand  
3100        pCMV-Target-AID     source        0      8752  +       
100         SV40 NLS            CDS           33     40    +       
200         null                source        33     5363  +       
300         Cas9(D10A)          CDS           64     4168  +       
500         SV40 NLS            CDS           4180   4201  +       
600         3xFLAG              CDS           4408   4474  +       
700         SV40 NLS            CDS           5107   5128  +       
800         UGI                 CDS           5134   5363  +       
1000        null                source        5409   8745  +       
1100        BGH-rev             primer_bind   5411   5429  -       
1200        bGH poly(A) signal  polyA_signal  5417   5642  +       
1300        M13 rev             primer_bind   5712   5729  -       
1400        M13 Reverse         primer_bind   5712   5729  -       
1500        M13/pUC Reverse     primer_bind   5725   5748  -       
1600        lac operator        protein_bind  5736   5753  +       
1700        lac promoter        promoter      5760   5791  -       
1800        CAP binding site    protein_bind  5805   5827  +       
1900        L4440               primer_bind   5943   5961  -       
2000        ori                 rep_origin    6114   6703  -       
2100        pBR322ori-F         primer_bind   6194   6214  -       
2200        AmpR                CDS           6873   7734  -       
2300        Amp-R               primer_bind   7496   7516  +       
2400        AmpR promoter       promoter      7734   7839  -       
2500        pRS-marker          primer_bind   7917   7937  -       
2600        CMV enhancer        enhancer      8108   8488  +       
2700        CMV promoter        promoter      8488   8692  +       
2800        CMV-F               primer_bind   8642   8663  +       
2900        T7                  primer_bind   8733   8745  +       
3000        T7 promoter         promoter      8733   8745  +  
```

**Example code 19: Export python script to reconstruct the DNA object**  
```pyton
#Source code (continued from the previous example)#
exporthistory(pCMV_Target_AID)

#Output#
from dna import *
DNA.dna_dict['pCMV_ABE'] = DNA(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'] = DNA(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129'
DNA.dna_dict['RS045'] = DNA(seq=None, record='input/RS045.fasta', project='RS045', topology='linear', format=None, process_description=description0)
DNA.dna_dict['HM129'] = DNA(seq=None, record='input/HM129.fasta', project='HM129', topology='linear', format=None, process_description=description0)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_0'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='3914/3914', end='6442/6442', project='pCMV-nCas-PmCDA1-ugi', process_description=description0)
DNA.dna_dict['HM129_0'] = flipdna(DNA.dna_dict['HM129'], project='HM129', process_description=description0)
DNA.dna_dict['N-term_Target-AID'] = joindna(*[DNA.dna_dict['RS045'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_0'], DNA.dna_dict['HM129_0']], topology='linear', project='N-term_Target-AID', process_description=description0)

description1 = 'The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046'
DNA.dna_dict['HM128'] = DNA(seq=None, record='input/HM128.fasta', project='HM128', topology='linear', format=None, process_description=description1)
DNA.dna_dict['RS046'] = DNA(seq=None, record='input/RS046.fasta', project='RS046', topology='linear', format=None, process_description=description1)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='6484/6484', end='9244/9244', project='pCMV-nCas-PmCDA1-ugi', process_description=description1)
DNA.dna_dict['RS046_0'] = flipdna(DNA.dna_dict['RS046'], project='RS046', process_description=description1)
DNA.dna_dict['C-term_Target-AID'] = joindna(*[DNA.dna_dict['HM128'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'], DNA.dna_dict['RS046_0']], topology='linear', project='C-term_Target-AID', process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048'
DNA.dna_dict['RS047'] = DNA(seq=None, record='input/RS047.fasta', project='RS047', topology='linear', format=None, process_description=description2)
DNA.dna_dict['RS048'] = DNA(seq=None, record='input/RS048.fasta', project='RS048', topology='linear', format=None, process_description=description2)
DNA.dna_dict['pCMV_ABE_0'] = cropdna(DNA.dna_dict['pCMV_ABE'], start='5277/5277', end='8613/8613', project='pCMV_ABE', process_description=description2)
DNA.dna_dict['RS048_0'] = flipdna(DNA.dna_dict['RS048'], project='RS048', process_description=description2)
DNA.dna_dict['backbone'] = joindna(*[DNA.dna_dict['RS047'], DNA.dna_dict['pCMV_ABE_0'], DNA.dna_dict['RS048_0']], topology='linear', project='backbone', process_description=description2)

description3 = 'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two insert fragments and a backbone fragments.'
DNA.dna_dict['N-term_Target-AID_0'] = modifyends(DNA.dna_dict['N-term_Target-AID'], left='*{25}/-{25}', right='-{28}/*{28}', project='N-term_Target-AID', process_description=description3)
DNA.dna_dict['C-term_Target-AID_0'] = modifyends(DNA.dna_dict['C-term_Target-AID'], left='*{28}/-{28}', right='-{25}/*{25}', project='C-term_Target-AID', process_description=description3)
DNA.dna_dict['backbone_0'] = modifyends(DNA.dna_dict['backbone'], left='*{25}/-{25}', right='-{25}/*{25}', project='backbone', process_description=description3)
DNA.dna_dict['pCMV-Target-AID'] = joindna(*[DNA.dna_dict['N-term_Target-AID_0'], DNA.dna_dict['C-term_Target-AID_0'], DNA.dna_dict['backbone_0']], topology='circular', project='pCMV-Target-AID', process_description=description3)
DNA.dna_dict['pCMV-Target-AID'].writedna('reconstructed_pCMV-Target-AID.gbk')
```

### Visualization function
dna.py provides a function to generate a simple sequence map of a DNA object, including genomic feature and their annotation labels. All feature and label locations are automatically adjusted to prevent objects on the sequence map from overlapping each other. Face color and edge color of feature objects are also automatically set by the default color set. However, if a feature holds `"qualifier:edgecolor_dna.py"` and `"qualifier:facecolor_dna.py"` in the qualifiers,  the qualifier values will be used for the face color and edge color, respectively.  

- **visualize**_`(dna_object=DNA object, map_view=str, feature_list=list, start=int, end=int, width_scale=float, height_scale=float, label_location=str, linebreak=int, seq=bool, diameter=float)`_  
	**Parameters**  
	- **dna_object**: *`DNA`* object
	- **map_view**: *`str`* (`"linear"` or `"circular"`; default: `"linear"`)
	Type of map_view for the given DNA object. If sequence topology of the DNA object is circular, circular view can be select for a sequence map.
	- **feature_list**: *`list`* of *`DNAfeaure`* objects (default: `.dnafeatures` excluding *`DNAFeature`* objects whose type is `"source"`)
	Features are displayed on a sequence map. 
	  
	Parameters available for only linear map view
	
	- **start**: *`int`* (default: 0)
	Start position of the sequence map. 
	- **end**: *`int`* (default: Length of the entire sequence)
	End position of the sequence map. 
	- **width_scale**: *`float`* (default: 1)
	Scaling factor for width. The default width is automatically adjusted based on the sequence length.
	- **height_scale**: *`float`* (default: 0)
	Scaling factor for height. The default height is automatically adjusted based on the width value. 
	- **label_location**: *`str`* (default: “both” if seq is False, otherwise “top”)
	Feature label location of feature objects on a sequence map. Feature labels are generally described in feature objects. However, label text width is larger than the feature object width, label text will be put at `label_location`. If the value is `"both"`, labels will be put below or above the feature objects, whichever is available. If the value is "`top`", labels will be put above the feature objects. If `seq` is `True`, the value will be `"top"`.
	- **linebreak**: *`int`* or `None` (default: Length of the entire sequence)
	Length of sequence for linebreak.
	- **seq**:  *`bool`* (default: `False`)
	If `True`, a colormap representing a nucleotide sequence will be viewed in a sequence map. 
	  
	Parameters available for only linear map view
	
	- **diameter_scale**: *`float`* (default: 1) 
	Scaling factor for diameter of circular map. 
	  
	**Return**  
	`matplolib.pyplot.figure` object
  
**Example code 20: Visualize DNA object**  
```python
#Source code (continued from the previous example)#
fig1 = visualize(frag_1)
fig2 = visualize(frag_2)
fig3 = visualize(frag_3, linebreak=100, seq=True)
fig4 = visualize(pCMV_Target_AID, map_view=”circular”)
```
