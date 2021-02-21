**README.md**

# dbrick.py Installation and User Manual
dbrick.py is a Python module that enables users to design favored DNA sequences in silico. This module provides six fundamental methods (`substr`,`join_dbricks`,`shell`,`reverse_complement`,`linearize`,`circularize`) to handle the Dbrick object that represents a DNA sequence with genetic information.
By using the methods, you can take partial fragments with sequence annotaion from genomic data such as GenBank, then generate new GenBank by combining their information. 

## Software dependency
Python 3.7.0 or later

## Installation
1. Download the software by  
   `git clone https://github.com/yachielab/base-editing-prediction`

2. Install the necessary Python packages  
   ``pip install matplotlib``  
   ``pip install biopython``

3. Set PYTHONPATH to the directory where you cloned the repository.

## Usage
The python module provides Dbrick class to handle a double-stranded DNA sequence with genomic annotations. Dbrick objects can be joined, split, and transformed by using the six fundamental methods. The usage of Dbrick objects and the six methods are explained in the following section. 

### Dbrick class
Dbrick objects can be created from Python string or Biopython Seq object specifying a DNA sequence. By giving a Fasta or GenBank format file, you can also create a Dbrick object with additional information for the DNA sequence.

**Examples** 
```Python
from dbrick import *
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

#Create Dbrick object from Python string
brick1 = Dbrick(seq="ATGC")
brick4 = Dbrick(record="example.gbk")

#Create Dbrick object from Bio.Seq.Seq object.  
seq   = Seq("ATGC", generic_dna)
brick2 = Dbrick(seq=seq)

#Create Dbrick object from GenBank or Fasta format file.  
brick3 = Dbrick(record="example.gbk")

#Create Dbrick object from GenBank or Fasta format file.  
brick4 = Dbrick(record="example.gbk")
```

Dbrick object also can take a DNA sequence with sticky ends as follows. 
**Examples** 
```Python
from dbrick import *
brick = Dbrick(seq="--ATGCGG/AATAGC--")
```



#### Properties and methods

Dbrick objects have some properties and methods to confirm and edit annotation information for the DNA sequences.  

##### Properties 

- `.name`   
An identifier of the Dbrick object. If a dbrick object is created from GenBank or Fasta format file, its identifier used as `.name`. 

- `.seq`  
A top strand DNA sequence of the Dbrick object. Even If the bottom strand sequence has a sticky end, the digested nucleotides on the top strand are complemented based on the bottom strand. When you confirm the dsDNA sequence structure with sticky ends, please use the `.view_seq` method described in the following section.
  
- `.topology`  
DNA sequence topology of the Dbrick object. Dbrick object can take `linear` or `circular` topology.

- `.features`  
It is a list of Biopython SeqFeature objects. Each SeqFeature gives annotation information within the specific region of the DNA sequence. Please refer to the following URL for the detailed explanation of the SeqFeature object. 
https://biopython.org/docs/dev/api/Bio.SeqFeature.html

- `.record`  
Biopython SeqRecord object that has an optional description and annotation for the DNA sequence. Please refer to the following URL for the detailed explanation of the SeqRecord object.



##### Methods

- :`.print_dsdna(self, whole=False, end_length=10)`  
Print double strand DNA sequence. It the `whole` argument is True or the sequence length is less than end length, print the full length of the DNA sequence. If it is `False`, display the `end_length` bp sequence at both ends.   
    **Example: Dbrick object of blunt end DNA sequence**  
  
    ```Python
    >>> from dbrick import *
    >>> brick = Dbrick(seq="CCGGTATGCGTCGA")
    >>> brick.print_dsdna()
    CCGGTATGCGTCGA
    GGCCATACGCAGCT
    ```

  
    **Example: Dbrick object of stickey end DNA sequence**    
    ```python
    >>> from dbrick import *
    >>> brick = Dbrick(seq="CCGGTATGCG----/----ATACGCAGCT")
    >>> brick.print_dsdna()
    CCGGTATGCG
        ATACGCAGCT
    ```

- `.print_feature_table(self)`

- `.add_features(self)`

- `.remove_features(self)`

  

#### Fudamental functions
The dbrick.py provides six fundamental functions to manuplate Dbrick objects. 

- `substr`  
 The `substr` method extracts a region within the specifcic range. The feature annotations in the region are carried over to the extracted fragment.  

    **Extraction of a specific region by position specification**  
    A region can be extraced by position specification. 

    ```Python
    from dbrick import *
    #Load GenBank format file as circular DNA seqeucne.
    brick     = dbrick_object(record="example.gbk", topology="circular")

    #Extract sub-region from 100 to 200.
    fragment1 = substr(brick, 100, 200)
    fragment2 = brick[100:200] #same as substr(brick, 1000, 2000)

    #Extract sub-region from 200 to 100.
    fragment3 = brick[200:100] #same as substr(brick, 2000, 1000)
    print(fragment1)
    print(fragment2)
    print(framment3)

    #Fetures in the region are also extracted with DNA seqeunce. 
    for feat in fragment1.features:
        print(feat)
    ```

    **Extraction of a specific region by sequence specification**  
    A region can be extraced by specification of the unique DNA sequences on the top strand and bottom strand respectively. If the sequence is not unique or does't exist in genome, the process will not work and return a error status.  

    ```Python
    from dbrick import *
    #X should be specified by 5'to 3' sequnence on top strand
    #Y should be specified by 5'to 3' sequnence on bottom strand
    fragment1 = brick[X:Y]  
    ```
- `linearize`

- `circularize`

- `reverse_complement`

- `shell`

- `join_dbricks`
**Examples : Join multiple Dbrick objects**  
Dbrick objects can be joined by using  `+` operand. However, the process doesn't check if the overhung sequences of these objects can be joined or not. If you want to link Dbrick objects with being aware of the overhang sequence structure, please use the `link_dricks` method or `^` operand. 

```Python
from dbrick import *
fragment1 = Dbrick(record="example1.gbk")
fragment2 = Dbrick(record="example2.gbk")
joined_object = join_dbricks(fragment1, fragment2)
```







