**README.md**
# dbrick.py Installation and User Manual
The dna.py is a Python module that enables users to handle double-stranded DNA sequences in silico. In the dna.py, all manipulation of DNA engineering was converged only four fundamental methods cut, join, end modification, and flip. By using the four methods, users can design and construct any DNA sequences. The module also gives some functions to read and edit the GenBank file, then constructed sequences are can be output as GenBank format with sequence annotation.　

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
The module provides “DNA class” to handle a double-stranded DNA sequence with sequence annotations. The usage of DNA objects and the four fundamental methods are described in the following section.

### Dbrick class
DNA class object can be created from Python string specifying a DNA sequence or sequence data in GenBank/FASTA format.

**Example 1: Create DNA object from python string (blunt end DNA sequence).**

```Python
#Soruce code#
from DNA import * 
brick  = DNA(seq="ATGC") 
```

**Example 2: Create DNA object from python string (sticky end sequence).**  
You can also specify a DNA sequence with sticky ends as follows. 

```Python
#Soruce code#
from DNA import * 
brick  = DNA(seq="CCGGTATGCG----/----ATACGCAGCT") 
```

The left side of "/" means the top strand of the sequence, and the right side means the bottom strand. The "-" indicated the region of the gap nucleotides. The bottom strand sequence should be a complementary sequence to the top strand, excluding the gap region in the format.

**Example 3: Create DNA object from GenBank format file.**  

```Python
from DNA import * 
brick  = Dbrick(seq="ATGC") 

```

#### Properties and methods
Dbrick objects have some properties and methods to confirm and edit annotation information for the DNA sequences.  

##### Properties 
- `.project`  
An identifier of the Dbrick object. If a dbrick object is created from GenBank or Fasta format file, its identifier used as `.name`. 

- `.seq`  
A top strand DNA sequence of the Dbrick object. Even If the bottom strand sequence has a sticky end, the digested nucleotides on the top strand are complemented based on the bottom strand. When you confirm the dsDNA sequence structure with sticky ends, please use the `.printdnaseq` method described in the following section.
  
- `.topology`  
DNA sequence topology of the Dbrick object. Dbrick object can take `linear` or `circular` topology.

- `.features`  
It is a list of Biopython SeqFeature objects. Each SeqFeature gives annotation information within the specific region of the DNA sequence. Please refer to the following URL for the detailed explanation of the SeqFeature object. 
https://biopython.org/docs/dev/api/Bio.SeqFeature.html

- `.record`  
Biopython SeqRecord object that has an optional description and annotation for the DNA sequence. Please refer to the following URL for the detailed explanation of the SeqRecord object.

##### Methods
- `.printdnaseq(self, whole=False, end_length=10)`  
  Print double strand DNA sequence. It the `whole` argument is True or the sequence length is less than end length, print the full length of the DNA sequence. If it is `False`, display the `end_length` bp sequence at both ends.   
  
  **Example 1 : Dbrick object of blunt end DNA sequence**  
  ```Python
  #Source code#
  from dbrick import *
  brick = Dbrick(seq="CCGGTATGCGTCGA")
  brick.printdnaseq()
  ````
  ````
  CCGGTATGCGTCGA
  GGCCATACGCAGCT
  ````
  
  **Example 2 : Dbrick object of stickey end DNA sequence**
  ````Python
  #Source code#
  from dbrick import *
  brick = Dbrick(seq="CCGGTATGCG----/----ATACGCAGCT")
  brick.print_dsdna()
  ````
  ````
  CCGGTATGCG
      ATACGCAGCT
  ````
  
- `.printdnafeature(self, sep=None, output=None, feature_types= None, detail=False, with_seq=False, zero_index=True)`  
  Print a table with the features information of the Dbrick object.  Default columns of the table are composed of `Feature ID`,  `Label`, `Type`, `Start`, `End`, `Strand`.  The `Feature ID` is the index of each feature in the .features.

  **Parameters**
  - `sep` : str  
      String for separating columns. if the value is None, the output is separated by multiple spaces to create a well-formatted table. 

  - `output` : str or file object  
      File name or file object for output. If the value is None   

  - `feature_types` : list of feature types  
      Kinds of feature types displayed in the output. If the value is None, the information of all features is written in the table. 

  - `detail` : bool  
    If the `detail` is `True`,  all qualifiers of each feature are written in the additional columns on the table.

  - with_seq : bool  
    If the `with_seq` is `True`,  the genomic dna sequence of each feature are written in the additional columns on the table.

  - zero_index : bool  
    As a default, positions of features are displayed by zero-based index (same as Python index ).  If the value is False, the positions are shown by 1-based counting. 

  **Example : Display the information of `primer_bind` features.**  
  ````python
  #Source code#
  from dbrick import *
  brick = DNA(record="pUC19.gbk")
  brick.print_features(with_seq=True, feature_types=["primer_bind"])
  ````
  ```
  Feature_ID  Label            Type         Start  End   Strand  Seq                      
  1           pBR322ori-F      primer_bind  117    137   +       GGGAAACGCCTGGTATCTTT     
  2           L4440            primer_bind  370    388   +       AGCGAGTCAGTGAGCGAG       
  6           M13/pUC Reverse  primer_bind  583    606   +       AGCGGATAACAATTTCACACAGG  
  7           M13 rev          primer_bind  602    619   +       CAGGAAACAGCTATGAC        
  8           M13 Reverse      primer_bind  602    619   +       CAGGAAACAGCTATGAC        
  11          M13 Forward      primer_bind  688    706   -       TGTAAAACGACGGCCAGT       
  12          M13 fwd          primer_bind  688    705   -       GTAAAACGACGGCCAGT        
  13          M13/pUC Forward  primer_bind  697    720   -       CCCAGTCACGACGTTGTAAAACG  
  14          pRS-marker       primer_bind  913    933   -       CGGCATCAGAGCAGATTGTA     
  15          pGEX 3'          primer_bind  1032   1055  +       CCGGGAGCTGCATGTGTCAGAGG  
  16          pBRforEco        primer_bind  1092   1111  -       AATAGGCGTATCACGAGGC      
  19          Amp-R            primer_bind  1501   1521  -       ATAATACCGCGCCACATAGC  
  ```

- `.adddnafeatures(self, start, end, strand=1, feature_type="misc_feature", qualifiers={})`  
  Add feature(s) on the specified location(s).

  **Parameters**
  - `start` (int or str)  
    Start position of the feature. 
 
  - `end` (int)  
    End position of the feature.
    
  - `strand` (1 or -1)  
    Direction of the feature.   

  - `feature_type` (str)  
    Type of the feature. 

  - `qualifiers` (dict)  
    Annotation infomation of the feature.

  **Example: Add new features** 
  ````Python
  #Source code#
  from dna import *
  brick = DNA(record="pUC19.gbk")
  brick.printdnafeature(feature_types=["misc_feature"], with_seq=True)
  ````
  ````
  Feature_ID  Label  Type          Start  End  Strand  Seq                                                        
  1000        MCS    misc_feature  631    688  +       AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC  
  ````
  
  ````Python 
  #Source code (continued from previous one)#
  brick.adddnafeature(100, 120, qualifiers={"label":"feat1"}) 
  brick.adddnafeature(2680, 10, qualifiers={"label":"feat2"}) #feat on origin
  brick.adddnafeature(*brick.finddna("GGATCC")[0].sspan, qualifiers={"label":"BamHI"})
  brick.printdnafeature(feature_types=["misc_feature"], with_seq=True)        
  ````
  
  ````
  Feature_ID  Label  Type          Start  End  Strand  Seq                                                        
  2001        feat1  misc_feature  100    120  +       ACGAGGGAGCTTCCAGGGGG                                       
  1000        MCS    misc_feature  631    688  +       AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC  
  801         BamHI  misc_feature  661    667  +       GGATCC                                                     
  1901        feat2  misc_feature  2680   10   +       CGAACTGAGATACCTA           
  ````
  
  
- `.remove_feature(self, feature_id)`  
  Remove a feature specified by the `feature_id`.

  **Example: Remove a feature** 
  ````python
  #Source code (continued from previous one)#
  brick.removednafeature("101")
  brick.printdnafeature(feature_types=["misc_feature"], with_seq=True)
  ````
  
  ````
  Feature_ID  Label  Type          Start  End  Strand  Seq                                                        
  1           feat1  misc_feature  100    120  +       ACGAGGGAGCTTCCAGGGGG                                       
  11          MCS    misc_feature  631    688  +       AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC  
  22          feat2  misc_feature  2680   10   +       CGAACTGAGATACCTA     
  ````
  
  

### Fudamental functions
The dbrick.py provides 6 fundamental functions to handle Dbrick objects. 

- `cropdna(brick, start, end)`   
   The `substr` method extracts sub-Dbrick object within the specific range. All feature information in the area, including features on the boundaries, is carried over to the extracted fragments. It adds the cropped part's report for the features on the edges. 
   
 **Parameters**
   - `brick ` : A DNA object  

   - `start` : int or str  (DNA sequence)  
     If the start is specified by an integer, it indicates the start position of the region you want to extract.
     
   - `end` : int or str (DNA sequence)  
     If the start is specified by an integer, it indicates the start position of the region you want to extract.
     
 **Example 1 : Extraction of a specific region**   
   ````Python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, 1000, 2000)
   brick.printdnafeature()
   print()
   sub_brick.printdnafeature()
   ````

   ````
    Feature_ID  Label             Type          Start  End   Strand  
    0           source            source        0      2686  +       
    1           feat1             misc_feature  100    120   +       
    2           pBR322ori-F       primer_bind   117    137   +       
    3           L4440             primer_bind   370    388   +       
    4           CAP binding site  protein_bind  504    526   +       
    5           lac promoter      promoter      540    571   +       
    6           lac operator      protein_bind  578    595   +       
    7           M13/pUC Reverse   primer_bind   583    606   +       
    8           M13 rev           primer_bind   602    619   +       
    9           M13 Reverse       primer_bind   602    619   +       
    10          lacZ-alpha        CDS           614    938   +       
    11          MCS               misc_feature  631    688   +       
    12          M13 Forward       primer_bind   688    706   -       
    13          M13 fwd           primer_bind   688    705   -       
    14          M13/pUC Forward   primer_bind   697    720   -       
    15          pRS-marker        primer_bind   913    933   -       
    16          pGEX 3'           primer_bind   1032   1055  +       
    17          pBRforEco         primer_bind   1092   1111  -       
    18          AmpR promoter     promoter      1178   1283  +       
    19          AmpR              CDS           1283   2144  +       
    20          Amp-R             primer_bind   1501   1521  -       
    21          ori               rep_origin    2314   217   +       
    22          feat2             misc_feature  2680   10    +       

    Feature_ID  Label                   Type         Start  End   Strand  
    0           source:1001..2000:2686  source       0      1000  +       
    1           pGEX 3'                 primer_bind  32     55    +       
    2           pBRforEco               primer_bind  92     111   -       
    3           AmpR promoter           promoter     178    283   +       
    4           AmpR:1..717:861         CDS          283    1000  +       
    5           Amp-R                   primer_bind  501    521   -      
   ````

   For the features on the boundary, the information about the cropped region is recorded according to the following format. 

   ```Feature identifier```:```Relative start position in the original feature```..```Relative end position in the original feature```:```The length of the original feature```
   

- `shell(brick, left_end=None, right_end=None)`

   The shell methods set sequence structures at both ends of the Dbrick object.  The end sequences are used when the Dbrick object is assembled with the other Dbrick object. 

   **Parameters**
   - `brick` : Dbrick object 

   - `left_end  ` : str  
     The left end sequence of the Dbrck object.   The way to specify the parameter is shown in the following examples.

   - `right_end` : str  
     The right end sequence of the Dbrck object.  The way to specify the parameter is shown in the following examples.

   **Example 1: Set sticky ends at both ends**  
   You can specify end structures by using "`*`" and "`-`," "`/`."  
   The left side of "/" means the top strand of the sequence, and the right side means the bottom strand. The "`-`" means the truncating nucleotides, and the " `*`" means the remaining nucleotides at the end. 

   ````Python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, 100, 120)
   sub_brick.print_dsdna()  
   sub_brick = shell(sub_brick,"----**/******","******/****--")
   sub_brick.print_dsdna() 
   ````
   ````
   ACGAGGGAGC...TTCCAGGGGG
   TGCTCCCTCG...AAGGTCCCCC
       GGGAGC...TTCCAGGGGG
   TGCTCCCTCG...AAGGTCCC  
   ````
   
   **Example 2: Add favor blunt ends at both ends.** 
   ````Python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, 100, 120)
   sub_brick = shell(sub_brick,"ATGTACG","ATGCTAC")
   sub_brick.print_dsdna() 
   ````
   
   ````
   ATGTACGACG...GGGATGCTAC
   TACATGCTGC...CCCTACGATG
   ````
   
   **Example 3: Add favor sticky ends at both ends.**
   ````python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, 100, 120)
   sub_brick = shell(sub_brick,"---ATGC/ATGTACG","TACG---/ATGCTAC")
   sub_brick.print_dsdna() 
   ````
   
   ````
      ACGAGGG...CAGGGGG   
   ATGTGCTCCC...GTCCCCCTAC
   ````
   
   In the above case,  the top strand sequence should be compliment with the bottom strand sequence. 
   
   **Example4: Simulate SacI/SalI digestion. **
   
   The SacI recognizes GAGCT^C site, and the SalI recognizes G^TCGAC site. Therefore, you can simulate SalI/SacI digestion of pUC19 as follows. 
   
   ````Python
   SacI = "GAGCTC"
   SalI = "GTCGAC"
   sub_brick = substr(brick, SacI, SalI) #Extract the sub-region between SacI and SalI.
   sub_brick.print_dsdna() 
   sub_brick = shell(sub_brick, "-----*/-*****", "*-----/*****-") #Set sticky end structures generated by SacI and SalI.
   sub_brick.print_dsdna() 
   ````
   
   ````
   The region from start 677 to end 655 was extracted.
   GAGCTCGAAT...GCAGGTCGAC
   CTCGAGCTTA...CGTCCAGCTG
       CGAATT...TGCAGG    
   TCGAGCTTAA...ACGTCCAGCT
   ````

- `join_dbricks(*bricks, topology="linear", name=None, ovhg_check=False, min_overlap=10, max_overlap=300)`

   The join_dbricks method joins given multiple Dbrick objects and returns an assembled Dbrick object. By specifying the ovhg_check parameter as  True, you can check if the assembly process can work or not.

   **Parameters**
   
   - `*bricks ` : Variable numbers of Dbrick objects.  
     Dbrick objects used for the joining. All Dbrick objects should be linear topology.

   - `topology` : "linear" or "circular"  
     The given dbrick ojects are assembled into a linear or circular constructs.  

   - `name` : str  
     The name of the assembled Dbrick objects. If the name is None, the name of the first element of the given Dbrick objects is used. 

   - `ovhg_check` : bool (True or False)    
     If ovhg_check is True, adjacent Dbrick objects should either have compatible sticky ends or homology sequences.   

   - `min_overlap` : int  
     Minimum required overlap length between adjacent fragments.  

   - `max_overlap` : int 
     Maximum required overlap length between adjacent fragments.  

     

   **Example1  : Join multiple dbrick objects based on compatible sticky ends**
   ````Python
   #Source code (continued from previous one)#
   EGFP = Dbrick(record="EGFP.fasta") 
   EGFP = shell(EGFP, SalI, SacI) 
   EGFP.print_dsdna()
   EGFP = shell(EGFP, "-*****/-----*", "*****-/*-----")
   EGFP.print_dsdna()
   product = join_dbricks(sub_brick, EGFP, topology="circular")
   ````

   ````
   GTCGACATGT...CTGAGAGCTC
   CAGCTGTACA...GACTCTCGAG
   TCGACATGTC...ACTGAGAGCT
       GTACAG...TGACTC    
   The Dbrick objects were joined based on complementary sticky end of each fragment. The sticky end is 'TCGA'
   The Dbrick object was circularized based on complementary sticky end between 3' end and 5' end. The sticky end is 'AGCT'
   ````

   **Example2  : Join multiple dbrick objects based on homology ends**
   ````Python
   sub_brick = substr(brick, "ACTGGCCGTCGTTTTACA", "GTCATAGCTGTTTCCTG")
   EGFP = shell(EGFP, "CAGGAAACAGCTATGAC", "ACTGGCCGTCGTTTTACA")
   product = join_dbricks(sub_brick, EGFP, topology="circular")
   ````

   `````
   The region from start 689 to end 619 was extracted.
   The Dbrick objects were joined based on sequence homology at the end of each fragment. The overhang sequence is 'CAGGAAACAGCTATGAC'
   The Dbrick object was circularized based on sequence homology between 3' end and 5' end. The overhang sequence is 'ACTGGCCGTCGTTTTACA'
   `````

- `linearize(brick, zero_position=0)`  
   The method linearizes circular Dbrick object at a specified position. 

   **Parameters**
   - `brick ` : Dbrick object

   - `zero positon` : int  
     The position will be the start point of the linearized Dbrick object.   

- `circularize(brick, ovhg_check=True)`  
   The method circularizes linear Dbrick object.  

   **Parameters**
   - `brick ` : Dbrick object
   
   - `ovhg_check` : bool  
     If the parameter is True, the method checks if both end sequence structures can be joined based on complementary sticky ends or homology sequences. 

- `reverse_complement(brick)`  
  The method generates the reverse complement of the specified Dbrick object.
  
  **Parameters**
  - `brick ` : Dbrick object
  
    
## Additional package

### Usage of visualize_dbrick.py  

#### Linear visualization

- 

#### Circular visualization 

- 

