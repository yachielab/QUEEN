# QUEEN Installation and User Manual

QUEEN (a framework to generate quinable and efficiently editable nucleotide sequence resources) is a Python programming module designed to describe, share credit DNA building processes and resources. DNA parts information can be imported from external annotated DNA files (GenBank and FASTA format). Output file (GenBank format) encodes the complete information of the constructed DNA and its annotations and enables the production of a quine code that self-reproduces the output file itself. In QUEEN, all of the manipulations required in DNA construction are covered by four simple operational functions, "cutdna", "modifyends", "flipdna", and "joindna" that can collectively represent any of the standard molecular DNA cloning processes, two search functions, "searchsequence" and "searchfeature", and two super functions, "editsequence" and "editfeature".  A new DNA can be designed by programming a Python script or using Jupyter Notebook, an interactive Python programming interpreter. The designed DNA product can be output in the GenBank file format that involves the history of its building process. The "quinable" feature of a QUEEN-generated GenBank file certifies that the annotated DNA material information and its production process are fully transparent, reproducible, inheritable, and modifiable by the community.

[https://www.dropbox.com/s/ykip1wbvc6ixedt/ga.jpg?dl=0](https://www.dropbox.com/s/ykip1wbvc6ixedt/ga.jpg?dl=0)


### Software dependency

Python 3.7.0 or later


### Installation

1. Download the QUEEN package from the GitHub repository. 

    `git clone [https://github.com/yachielab/QUEEN](https://github.com/yachielab/QUEEN)`

2. Move to the QUEEN directory and install QUEEN using the following command.

    `python setup.py install` 

3. Install Graphviz (optional; required for visualizing flowcharts of DNA building processes using the `visualizeflow()` function described below). Graphviz package is available at the following link.

    - [Graphviz](https://graphviz.org/download/source/)


### Usage

QUEEN provides the QUEEN class to define a double-stranded (ds)DNA object with sequence annotations. The QUEEN class and its operational functions are described below. Jupyter Notebook files for all of the example codes are provided in `./demo/tutorial` of QUEEN ([https://github.com/yachielab/QUEEN](https://github.com/yachielab/QUEEN)) and also made executable in Google Colaboratory:

- https://colab.research.google.com/drive/1ubN0O8SKXUr2t0pecu3I6Co8ctjTp0PS?usp=sharing
- https://colab.research.google.com/drive/1dPcNhsOl2sne_wq7ZULXXFUxizR6JQrR?usp=sharing

### QUEEN class 

The QUEEN class can define a dsDNA object with sequence annotations. It can be created by specifying a DNA sequence or importing a sequence data file in GenBank or FASTA file format (single sequence entry). When a GenBank format file is imported, its NCBI accession number, Addgene plasmid ID, or Benchling share link can be provided instead of downloading the file to your local environment.


##### Example code 1: Create a QUEEN class object (blunt-ends) 

A `QUEEN_object` (blunt-end) is created by providing its top-stranded sequence (5’-to-3’). By default, the DNA topology will be linear.

Source code

```python
from QUEEN.queen import *
dna = QUEEN(seq="CCGGTATGCGTCGA") 
```


##### Example code 2: Create a QUEEN class object (sticky-end)

The left and right values separated by `"/"` show the top and bottom strand sequences of the generating `QUEEN_object`, respectively. The top strand sequence is provided in the 5’-to-3’ direction from left to right, whereas the bottom strand sequence is provided in the 3′-to-5′ direction from left to right. Single-stranded regions can be provided by `"-"` for the corresponding nucleotide positions on the opposite strands. A:T and G:C base-pairing rule is required between the two strings except for the single-stranded positions.

**Source code**

```python
from QUEEN.queen import *
dna = QUEEN(seq="CCGGTATGCG----/----ATACGCAGCT") 
```


##### Example code 3: Create a circular QUEEN class object

The sequence topology of generating `QUEEN_object` can be specified by `"linear"` or` "circular"`.

**Source code**

```python
from QUEEN.queen import *
dna = QUEEN(seq="CCGGTATGCGTCGA", topology="circular") 
```


##### Example code 4.1: Create a QUEEN class object from a GenBank file in a local directory 

GenBank file can be loaded by specifying its local file path.

**Source code**

```python
from QUEEN.queen import *
pUC19 = QUEEN(record="./input/pUC19.gbk")
```


##### Example code 4.2: Create a QUEEN class object using a NCBI accession number

QUEEN_object can be generated from a NCBI accession number with `dbtype="ncbi"`. 

**Source code**

```python
from QUEEN.queen import *
#"M77789.2" is NCBI accession number for pUC19 plasmid
pUC19 = QUEEN(record="M77789.2", dbtype="ncbi") 
```


##### Example code 4.3: Create a QUEEN class object using an Addgene plasmid ID

`QUEEN_object` can be generated from an Addgene plasmid ID with `dbtype="addgene"`.

**Source code**

```python
from QUEEN.queen import *
#"50005" is Addgene plasmid ID for pUC19 plasmid
pUC19 = QUEEN(record="50005", dbtype="addgene")
```

##### Example code 4.4: Create a QUEEN class object from a Benchling share link

`QUEEN_object` can be generated from a Benchling shared link with `dbtype="benchling"`.

**Source code**

```python
from QUEEN.queen import *
plasmid = QUEEN(record="https://benchling.com/s/seq-U4pePb09KHutQzjyOPQV", dbtype="benchling")
```

pX330 plasmid encoding a Cas9 gene and a gRNA expression unit is provided in the above example. The `QUEEN_object` generated here is used in the following example codes in this document.


#### Properties of QUEEN class objects 

* **.project**: `str`  
Project name of `QUEEN_object`_ _construction. In QUEEN, this property is also used as a dictionary key to access the `.productdict` described below. If a `QUEEN_object` is created from a GenBank or FASTA file, its sequence ID will be inherited here. Otherwise, the project name is automatically generated to be unique amongst the existing `.productdict` keys.

* **.seq**: `str`  
Top strand sequence (5′→3′). Sticky end gaps, if any, can be represented by `"-"`. This property cannot be directly edited; only the built-in operational functions of QUEEN described below can edit this property_._

* **.rcseq**: `str`  
Bottom strand sequence (5′→3′). Sticky end gaps, if any, are represented by `"-"`. This property cannot be directly edited; only the built-in operational functions of QUEEN described below can edit this property_._

* **.topology**: `str` (`"linear"` or `"circular"`)   
Sequence topology. When a `QUEEN_object`_ _is created by loading from a GenBank file, the topology is set according to the description in the GenBank file. Only the built-in operational functions of QUEEN described below can edit this property_._

* **.dnafeatures**: `list` of `DNAfeature_objects`  
When a `QUEEN_object` is loaded from a GenBank file, `.dnafeatures` will automatically be generated from the GenBank file's sequence features. Otherwise, `.dnafeatures` will be an empty list. Each `DNAfeature_object` with the following attributes provides an annotation for a given range of DNA sequence in a `QUEEN_object`.

    * **.feature_id: **`str`    
      Unique identifier. It is automatically determined to each feature when a `QUEEN_object` is loaded from a GenBank file.

    * **.feature_type: **`str`    
      Biological nature. Any value is acceptable. The GenBank format requires registering a biological nature for each feature.

    * **.start: **`int`  
      Start position of `DNAfeature_object` in `QUEEN_object`.

    * **.end: **`int`  
      Start position of `DNAfeature_object` in `QUEEN_object`.

    * **.strand: **`int (1 or -1)`  
      Direction of `DNAfeature_object` in `QUEEN_object`._ _Top strand (`1`) or bottom strand (`-1`).  
      
    * **.sequence:** `str`  
      Sequence of the `DNAfeature_object` for its encoded direction.
      
    * **.qualifiers: **`dict`  
Qualifiers. When a GenBank file is imported, qualifiers of each feature will be registered here. Qualifier names and values will serve as dictionary keys and values, respectively. 

`DNAfeature_object` can be edited only by the `editfeature()` function described below. DNAfeature class is implemented as a subclass of the Biopython_ _SeqFeature class. Therefore, apart from the above attributes, DNAfeature class inherits all the attributes and methods of SeqFeature class. For details about SeqFeature class, see (https://biopython.org/docs/dev/api/Bio.SeqFeature.html) 

* **.productdict**: `dict`  
Dictionary for all of the inherited `QUEEN_objects`_ _used to construct the present `QUEEN_object`. The `.project` of each `QUEEN_object` serves as a key of this dictionary.

### Output functions

`QUEEN_objects` hold a simple set of functions to output its information.

* **`.printsequence(start=int, end=int,  strand=int, display=bool, hide_middle=int, linebreak=int)`**  
  Returns and displays partial or the entire dsDNA sequence and sequence end structures of `QUEEN_object`.
  
  
    ##### Parameters
  
    * **start: **`int` (zero-based indexing; default: `0`)  
      Start position of the sequence. 
  
    * **end: **`int` (zero-based indexing; default: the last sequence position of `QUEEN_object`)  
      End position of the sequence.
  
    * **strand:** `int`: `1` (top strand only), `-1` (bottom strand only), or `2` (both strands) (default: `2`)  
      Sequence strand(s) to be returned.
  
    * **display:** `bool` (`True` or `False`; default: `True`)   
      If `True`, the output will be displayed in `STDOUT`.
  
    * **hide_middle: **`int` or `None` (default: `None`)  
      Length of both end sequences to be displayed.
  
    * **linebreak:** `int` (default: length of the `QUEEN_object` sequence)  
      Length of sequence for linebreak.

  ##### Return  
  > If **strand** is `1` or `-1`, sequence of the defined strand (5’→3’)
  > If **strand** is `2`, `"str/str"`: "top strand sequence (5’→3’)/bottom strand sequence (3’→5’)"

##### Example code 5: Print a dsDNA object
**Source code**

```python
from queen import *
fragment **=** QUEEN(seq**=**"CCGGTATGCG----/----ATACGCAGCT") 
fragment.printsequence(display=True)
```
**Output**
```
5′ CCGGTATGCG---- 3′
3′ ----ATACGCAGCT 5′
```



* **`.printfeature(feature_list=list, attribute=list, seq=bool, separation=str, output=str, x_based_index=int)`**  
    Print a tidy data table of annotation features/attributes of `QUEEN_object`. Default output attributes are `"feature_id"`, `"feature_type"`, `"qualifiers:label"`, `"start"`, `"end"`, and `"strand"`.
    
    
    ##### Parameters
    
    * **feature_list: **`list` of `DNAfeaure_objects` (default: `.dnafeatures`)  
      List of features to be displayed in the output table. If not given, all features held by the QUEEN_object will be the subject.
    
    * **attribute:** `list` of feature attributes (default: `["feature_id", "feature_type", "qualifiers:label", "start", "end", "strand"]`)   List of feature attributes to be displayed in the output table. If the value is `"all"`, it will generate a table for all the attributes held by the `QUEEN_object` except for `"sequence"`.
    
    * **seq**: `bool` (`True` or `False`; default: `False`)  
      If `True`, the sequence of each feature for its encoded direction will be displayed in the output table.
    
    * **separation:** `str` (default: space(s) to generate a well-formatted table)  
      String to separate values of each line.
    
    * **output:** `str` (default: `STDOUT`)  
      Output file name.
    
    * **x_based_index:** `0` or `1` (default: `0`)  
      As a default, positions of all features are given in the zero-based indexing in QUEEN (same as Python). If this parameter is set to `1`, they will be shown in the 1-based indexing (as seen in the GenBank format).
    
    ##### Return
    > `None`


##### Example code 6: Print DNA features in a well-formatted table

**Source code**
```python
from queen import *
plasmid = QUEEN(record="input/px330.gb") 
plasmid.printfeature()
```
**Output**

```
feature_id  feature_type   qualifier:label     start  end   strand  
1           source         source              0      8484  +       
100         primer_bind    hU6-F               0      21    +       
200         promoter       U6 promoter         0      241   +       
300         primer_bind    LKO.1 5'            171    191   +       
400         misc_RNA       gRNA scaffold       267    343   +       
500         enhancer       CMV enhancer        439    725   +       
600         intron         hybrid intron       983    1211  +       
700         regulatory     Kozak sequence      1222   1232  +       
800         CDS            3xFLAG              1231   1297  +       
900         CDS            SV40 NLS            1303   1324  +       
1000        CDS            Cas9                1348   5449  +       
1100        CDS            nucleoplasmin NLS   5449   5497  +       
1200        primer_bind    BGH-rev             5524   5542  -       
1300        polyA_signal   bGH poly(A) signal  5530   5738  +       
1400        repeat_region  AAV2 ITR            5746   5876  +       
1500        repeat_region  AAV2 ITR            5746   5887  +       
1600        rep_origin     f1 ori              5961   6417  +       
1700        primer_bind    F1ori-R             6048   6068  -       
1800        primer_bind    F1ori-F             6258   6280  +       
1900        primer_bind    pRS-marker          6433   6453  -       
2000        primer_bind    pGEX 3'             6552   6575  +       
2100        primer_bind    pBRforEco           6612   6631  -       
2200        promoter       AmpR promoter       6698   6803  +       
2300        CDS            AmpR                6803   7664  +       
2400        primer_bind    Amp-R               7021   7041  -       
2500        rep_origin     ori                 7834   8423  +       
2600        primer_bind    pBR322ori-F         8323   8343  +     
```



* **`.outputgbk(output=str)`**  
    Output `QUEEN_object` to a GenBank file. In addition to all of the `DNAfeature_objects` in the input `QUEEN_object`, a `DNAfeature_object` encoding the entire construction processes that generated the `QUEEN_object` in `qualifiers:building_history` will also be output to the GenBank file.
    
    ##### Parameter
    * **output:** str (default: STDOUT)
    
        Output file name.
    
    #####  Return
    > `None`

### Search Function

`QUEEN_objects` hold `.searchsequene()` and `.searchfeature()` functions that enables users to search for query sequences and values in `DNAfeature_objects`.

* **`.searchsequence (query=regex or str, start=int, end=int, strand=int, product=str, process_name=str, process_description="str")`**    
  Search for specific sequences from a user-defined region of a `QUEEN_object` and return a list of `DNAfeature_objects`. Start and end attributes of returned `DNAfeature_objects` represent the sequence regions of the `QUEEN_object` that matched the user's query. Note that the returned `DNAfeature_objects` will not be generated with `.feature_id` and reflected to the parental `QUEEN_object`**. **The returned `DNAfeature_objects` can be added to `QUEEN_object.dnafeatures` by `editfeature()` with the `createattribute` option as explained below.
  
  
  ##### Parameters
  
  * **query: **`regex` or `str` (default: `".+"`)  
    Search query sequence. If the value is not provided, the user-specified search region of the `QUEEN_object` sequence with `start` and `end` explained below will be returned. It allows fuzzy matching and regular expression. For details, see [https://pypi.org/project/regex/](https://pypi.org/project/regex/). All IUPAC nucleotide symbols can be used. Restriction enzyme cut motif representation can be used to define a query with `"^"` and `"_"` or `"(_int_/_int_)"`. For example, EcoRI cut motif can be provided by `"G^AATT_C"`, where `"^"` and `"_"` represent the cut positions on the top and bottom strands, respectively, or by `"GAATTC(-5/-1)"` or `"(-5/-1)GAATTC"`, where the left and right integers between the parentheses represent the cut positions on the top and bottom strands, respectively. Similarly, the cut motif of a Type-IIS restriction enzyme BsaI can be given by `"GGTCTCN^NNN_N"`,  `"N^NNN_NGAGACC"`, `"GGTCTC(1/5)"` or `"(5/1)GAGACC"`. The returned `DNAfeature_objects`_ _obtained for a query restriction enzyme cut motif will hold the cutting rule in the `"qualifiers:cutsite"` attribute, which can be added to `QUEEN_object.dnafeatures` by `editfeature()` with the `createattribute` option as explained below. Regular expression is disabled for restriction enzyme cut motifs.
  
  * **start: **`int` (zero-based indexing; default: `0`)  
    Start position of the target range of the `QUEEN_object` sequence for the search. 
  
  * **end: **`int` (zero-based indexing; default: the last sequence position of `QUEEN_object`)  
    End position of the target range of the `QUEEN_object` sequence for the search. 
  
  * **strand:** `int`: `1` (top strand only), `-1` (bottom strand only), or `2` (both strands) (default: `2`)  
    Sequence strand to be searched.
  
  ##### Return
  > `list` (`list` of `DNAfeature_objects`)



##### Example code 7: Search for a DNA sequence motif with regular expression

**Source code (continued from the previous code)**

```python
match_list = plasmid.searchsequence(query="G[ATGC]{19}GGG")
plasmid.printfeature(match_list, seq=True, attribute=["start", "end", "strand"])
```

**Output**

```
start  end   strand  sequence                 
115    138   +       GTAGAAAGTAATAATTTCTTGGG  
523    546   +       GACTTTCCATTGACGTCAATGGG  
816    839   +       GTGCAGCGATGGGGGCGGGGGGG  
1372   1395  +       GACATCGGCACCAACTCTGTGGG  
1818   1841  +       GGCCCACATGATCAAGTTCCGGG  
3097   3120  +       GATCGGTTCAACGCCTCCCTGGG  
3300   3323  +       GCGGCGGAGATACACCGGCTGGG  
3336   3359  +       GAAGCTGATCAACGGCATCCGGG  
3529   3552  +       GGCAGCCCCGCCATTAAGAAGGG  
3577   3600  +       GACGAGCTCGTGAAAGTGATGGG  
︙
493    516   -       GCGTTACTATTGACGTCAATGGG  
654    677   -       GTCCCATAAGGTCATGTACTGGG  
758    781   -       GGTGGGGAGGGGGGGGAGATGGG  
1014   1037  -       GCGCGAGGCGGCGGCGGAGCGGG  
1301   1324  -       GACCTTCCGCTTCTTCTTTGGGG  
1820   1843  -       GCCCCGGAACTTGATCATGTGGG  
2090   2113  -       GAAGTTGCTCTTGAAGTTGGGGG  
2183   2206  -       GGCGTACTGGTCGCCGATCTGGG  
2288   2311  -       GATCATAGAGGCGCTCAGGGGGG  
2689   2712  -       GCCAGAGGGCCCACGTAGTAGGG  
︙
```



##### Example code 8: Search for a DNA sequence motif with fuzzy matching

Search for `"AAAAAAAA"` sequence, permitting a single nucleotide mismatch**.**

**Source code (continued from the previous code)**

```python
match_list = plasmid.searchsequence(query="(?:AAAAAAAA){s&lt;=1}")
plasmid.printfeature(match_list, seq=True) 
```

**Output**

```
feature_id  feature_type  qualifiers:label  start  end   strand  sequence  
null        misc_feature  null              5484   5492  +       AAAAAAGA  
null        misc_feature  null              6369   6377  +       AACAAAAA  
null        misc_feature  null              7872   7880  +       AAACAAAA  
null        misc_feature  null              346    354   -       AAAACAAA  
null        misc_feature  null              799    807   -       AAAAAATA  
null        misc_feature  null              1201   1209  -       GAAAAAAA  
null        misc_feature  null              6716   6724  -       AAAAATAA  
null        misc_feature  null              7844   7852  -       AGAAAAAA 
```



##### Example code 9: Search for a DNA sequence with the IUPAC nucleotide code

**Source code (continued from the previous code)**

```python
match_list = plasmid.searchsequence(query="SWSWSWDSDSBHBRHH")
plasmid.printfeature(match_list, seq=True)
```

**Output**

```
feature_id  feature_type  qualifiers:label  start  end   strand  sequence          
null        misc_feature  null              4098   4114  +       GAGACAGCTGGTGGAA  
null        misc_feature  null              3550   3566  -       CTGTCTGCAGGATGCC  
null        misc_feature  null              5239   5255  -       CTCTGATGGGCTTATC  
null        misc_feature  null              6415   6431  -       GAGAGTGCACCATAAA  
null        misc_feature  null              8357   8373  -       GTCAGAGGTGGCGAAA  
```



* **`.searchfeature(key_attribute=str, query=regex or str, source=list of DNAfeature_objects, start=int, end=int, strand=int, product=str, process_name=str, process_description="str")`**  
  Search for `DNAfeature_objects` holding a queried value in a designated `key_attribute` in `QUEEN_object`.
  
  
    ##### Parameters
  
    * **key_attribute**: `str` (default: `"all"`)  
      Attribute type to be searched (`feature_id`, `feature_type`, `"qualifiers:*"`, or `sequence`). If the value is not provided, it will be applied to all of the attributes in the `QUEEN_object`, excluding `sequence`. However, if the `query` value is provided with only the four nucleotide letters (A, T, G, and C), this value will be automatically set to `sequence`.
  
    * **query**: `regex`_ _or `str`(default: `".+"`)  
      Query term. `DNAfeature_objects` that have a value matches to this query for `key_attribute` designated above will be returned. It allows fuzzy matching and regular expression. For details, see [https://pypi.org/project/regex/](https://pypi.org/project/regex/). If the `key_attribute` is `sequence`, all IUPAC nucleotide symbols can be used.
  
    * **source: **`list`_ _of_ `DNAfeature_objects` (default: `QUEEN_object.dnafeatures`)  
      Source `DNAfeature_objects` to be searched. `DNAfeature_objects` outside the search range defined by `start`, `end`, and `strand` will be removed from the source. Any `DNAfeature_objects` can be provided here. For example, a list of `DNAfeature_objects` _returned from another `searchsequence()` or `searchfeature()` operation can be used as the source to achieve an AND search with multiple queries.
    * **start: **`int` (zero-based indexing; default: `0`)  
      Start position of the target range of the `QUEEN_object` sequence for the search. 
  
    * **end: **`int` (zero-based indexing; default: the last sequence position of `QUEEN_object`)    
      End position of the target range of the `QUEEN_object` sequence for the search. 
    * **strand:** `int`: `1` (top strand only), `-1` (bottom strand only), or `2` (both strands) (default: `2`)  
      Sequence strand to be searched.
  
  
    ##### Return
    > `list` (`list` of `DNAfeature_objects`)

  ​               


##### Example code 10: Search for sequence features having specific attribute values  

Search for `DNAfeature_objects` with a feature type `"primer_bind"`, and then further screen ones holding a specific string in `"qualifiers:label"`.

**Source code (continued from the previous code)**

```python
feature_list = plasmid.searchfeature(key_attribute="feature_type", query="primer_bind")
plasmid.printfeature(feature_list)
sub_feature_list = plasmid.searchfeature(key_attribute="qualifier:label", query=".+-R$", source=feature_list)
plasmid.printfeature(sub_feature_list)
```

**Output**

```
feature_id  feature_type  qualifiers:label  start  end   strand  
200         primer_bind   hU6-F            0      21    +       
300         primer_bind   LKO.1 5'         171    191   +       
1200        primer_bind   BGH-rev          5524   5542  -       
1700        primer_bind   F1ori-R          6048   6068  -       
1800        primer_bind   F1ori-F          6258   6280  +       
1900        primer_bind   pRS-marker       6433   6453  -       
2000        primer_bind   pGEX 3'          6552   6575  +       
2100        primer_bind   pBRforEco        6612   6631  -       
2400        primer_bind   Amp-R            7021   7041  -       
2600        primer_bind   pBR322ori-F      8323   8343  +       
feature_id  feature_type  qualifiers:label  start  end   strand  
1700        primer_bind   F1ori-R          6048   6068  -       
2400        primer_bind   Amp-R            7021   7041  -   
```