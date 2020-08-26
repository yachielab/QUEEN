**README.md**

#dbrick.py Installation and User Manual

dbrick.py is a Python module that enables users to design favored DNA sequences in silico. This module is composed Dbrick class and four fundamental methods that handling a DNA sequence as a brick. By using the methods, you can take partial fragments with sequence annotaion frorm genomic data such as GenBank, then generate new GenBank by combining their information. 

## Software dependency

Python 3.7.0 or later

## Installation

1. Download the software by
   `git clone https://github.com/yachielab/base-editing-prediction`

2. Install the necessary Python packages
   ``pip install matplotlib``
   ``pip install biopython``



## Usage

Dbrick class provides object to handle DNA seqeunces with identifiers, genomic infomation and optional annotaions for features such as genes.  Dbrick object DNA sequence as python `str`, ``Bio.Seq.Seq`` or ``Bio.Seq.SeqRecord`` object  



**Create Dbrick object from DNA sequnece**

````
from dbrick import *
brick = Dbrick(seq="ATGC")
````
Dbrick can also be created from  ``Bio.Seq.Seq`` object as an argument.
```
from dbrick import *
from Bio import Seq
from Bio.Alphabet import generic_dna
seq   = Seq("AGTACACTGGT", generic_protein)
brick = Dbrick(seq=seq)
```



**Create Dbrick object from genomic Data**

````
from dbrick import *
brick = Dbrick(record="example.gbk")
````
You can also give  ``Bio.Seq.SeqRecord`` object as an argument. 
```
from dbrick import *
from Bio import SeqIO
record = SeqIO.read("example.gbk",format="genbank")   
brick  = Dbrick(record=record)
```



**Isolate Dbrick object in a specific region**
By using slice indexing, you can take a fragment in the specified region . The feature annotations in the region are carried over to the extracted fragment.   

````
from dbrick import *
brick     = dbrick_object(record="example.gbk")  #The topology of examle.gbk is circular
fragment1 = brick[1000:2000]                      #same as substr(apart, 1000, 2000)
fragment2 = brick[1000:2000:-1]                   #same as complement(substr(apart, 1000, 2000))
fragment3 = brick[2000:1000]                      #same as substr(apart, 2000, 1000)
print(fragment1.seq)
print(fragment2.seq) 

#Fetures linked with seqeunce are also carried. 
for feat in sub_brcik1.features:
	print(feat)
````

In slice indexing, a unique DNA sequence in the top strand can be used for start position specification, and a unique DNA sequence in the bottom strand can be used for end position specification, as follows. If the sequence is not unique or does't exist in the strand , the process will not work and return a error status.  

````
from dbrick import *
fragment1 = brick[XXXXX:YYYYY] #Y should be specified by 5'to 3' sequnence in bottom strand` 
````



**Join multiple Dbrick objects** 
Dbrick objects can be joined by using  `+` operand. However, the process doesn't check if the overhung sequences of these objects can be joined or not. If you want to link Dbrick objects with being aware of the overhang sequence structure, please use the `link_dricks` method or `^` operand. 

````
from dbrick import *
fragment1 = Dbrick(record="example1.gbk")
fragment2 = Dbrick(record="example2.gbk")
print(fragment1.seq) 
print(fragment2.seq) 
joined_fragment1 = fragment1 + fragment2
print(joined_fragment3.seq)
joined_fragment2 = fragment1 ^ fragment2
print(joined_fragment4.seq)
````



**API**

####class `Dbrick(seq=None, record=None)`

> **Properties**
> *record*
> Bio.Seq.SeqRecord object. 
>
> 
>
> *features* 
> List of Bio.SeqFeature.SeqFeature objects.
>
> 
>
> *name*
> Specific name for Dbrick object.
>
> 
>
> *topology*
> Topology of DNA seqeunce. The parameter is specified by "linear" or "cicular".
>
> 
>
> **Methods**
> *write(output_name)*
> Output Dbrick contents as GenBank format 
>
> 
>
> *circularize(ovhg_check =True, min_overlap=10, max_overlao=500)*
> Circularize Dbrick object. As a default, the method check if three prime end and five prime end have common overhang sequence for joining. If you don't need to care it, please specified *ovhg_check* as ``False``.

`

### 







