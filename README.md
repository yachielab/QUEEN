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

```Python
#Source code#
from dbrick import *
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

#Create Dbrick object from Python string. 
brick1 = Dbrick(seq="ATGC")

#Create Dbrick object from Bio.Seq.Seq object.  
seq   = Seq("ATGC", generic_dna)
brick2 = Dbrick(seq=seq)

#Create Dbrick object from GenBank or Fasta format file.  
brick3 = Dbrick(record="example.gbk")

#Create Dbrick object from GenBank or Fasta format file.  
brick4 = Dbrick(record="example.gbk")
```
Dbrick object also can take a DNA sequence with sticky ends as follows. The ''-'' 
```Python
#Source code#
from dbrick import *
brick = Dbrick(seq="CCGGTATGCG----/----ATACGCAGCT")
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
- `.print_dsdna(self, whole=False, end_length=10)`  
  Print double strand DNA sequence. It the `whole` argument is True or the sequence length is less than end length, print the full length of the DNA sequence. If it is `False`, display the `end_length` bp sequence at both ends.   
  
  **Example 1 : Dbrick object of blunt end DNA sequence**  
  ```Python
  #Source code#
  from dbrick import *
  brick = Dbrick(seq="CCGGTATGCGTCGA")
  brick.print_dsdna()
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
  
- `.print_features(self, sep=None, output=None, feature_types= None, detail=False, with_seq=False, zero_index=True)`  
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
  brick = Dbrick(record="pUC19.gbk")
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

- `.add_features(self, start, end=None, strand=1, feature_type="misc_feature", qualifiers={})`  
  Add feature(s) on the specified location(s).

  **Parameters**
  - `start` (int or str)  
    If the start is specified by an integer, it indicates the start position of the feature. You can also set it by the DNA sequence on the strand specified by the strand parameter. In this case, the feature is set for all regions composed of this DNA sequence.   
    
  - `end` (int)  
    End position of the feature. If the `start` is specified by the DNA seqeunce, the `end` value is ignored.  
    
  - `strand` (1 or -1)  
    Direction of the feature.   

  - `feature_type` (str)  
    Type of the feature. 

  - `qualifiers` (dict)  
    Annotation infomation of the feature.

  **Example: Add new features** 
  ````Python
  #Source code#
  from dbrick import *
  brick = Dbrick(record="pUC19.gbk")
  brick.print_features(feature_types=["misc_feature"], with_seq=True)
  ````
  ````
  Feature_ID  Label  Type          Start  End  Strand  Seq                                                        
  10          MCS    misc_feature  631    688  +       AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC 
  ````
  ````Python 
  #Source code (continued from previous one)#
  brick.add_feature(100, 120, qualifiers={"label":"feat1"}) 
  brick.add_feature(2680, 10, qualifiers={"label":"feat2"}) #feat on origin
  brick.add_feature("GGATCC", qualifiers={"label":"BamHI"})
  brick.print_features(feature_types=["misc_feature"], with_seq=True)           
  ````
  
  ````
  New feature was added in the range of start 661 to end 667.
  Feature_ID  Label  Type          Start  End  Strand  Seq                                                        
  1           feat1  misc_feature  100    120  +       ACGAGGGAGCTTCCAGGGGG                                       
  11          MCS    misc_feature  631    688  +       AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC  
  12          BamHI  misc_feature  661    667  +       GGATCC                                                     
  23          feat2  misc_feature  2680   10   +       CGAACTGAGATACCTA        
  ````
  
- `.remove_feature(self, feature_id)`

  Remove a feature specified by the `feature_id`.

  **Example: Remove a feature** 
  ````python
  #Source code (continued from previous one)#
  brick.remove_feature("12")
  brick.print_features(feature_types=["misc_feature"], with_seq=True)
  ````
  
  ````
  Feature_ID  Label  Type          Start  End  Strand  Seq                                                        
  1           feat1  misc_feature  100    120  +       ACGAGGGAGCTTCCAGGGGG                                       
  11          MCS    misc_feature  631    688  +       AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC  
  22          feat2  misc_feature  2680   10   +       CGAACTGAGATACCTA     
  ````
  
  

### Fudamental functions
The dbrick.py provides 6 fundamental functions to handle Dbrick objects. 

- `substr(brick, start, end, target=None)`   
   The `substr` method extracts sub-Dbrick object within the specific range. All feature information in the area, including features on the boundaries, is carried over to the extracted fragments. It adds the cropped part's report for the features on the edges. 

   **Parameters**
   - `brick ` : A Dbrick object  

   - `start` : int or str  (DNA sequence)  
     If the start is specified by an integer, it indicates the start position of the region you want to extract. You can also set the DNA sequence on the forward strand of the Dbrick object. In this case, the start of the DNA sequence is the start position of the region, and the DNA sequence should be unique on the forward strand. 

   - `end` : int or str (DNA sequence)  
     If the start is specified by an integer, it indicates the start position of the region you want to extract. You can also set the DNA sequence on the reverse complement strand of the Dbrick object. In this case, the start of the DNA sequence is the start position of the region, and the DNA sequence should be unique on the reverse complement strand. 

   - `target`:  int (Feature ID), list (Feature IDs)  or str (DNA seqeunce)  
     If you want to extract the specific area with the flanking region, please use this parameter. You can set Feature_ID (s) or the DNA sequence on the forward strand to the parameter. In that case, it extracts the area, including the feature(s) or the sequence, with the flanking region. The flanking length at both ends are specified by the `start` and `end` parameters, respectively. 

   **Example 1 : Extraction of a region specified by the position**     
   ````Python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, 1000, 2000)
   brick.print_features
   print()
   sub_brick.print_features() 
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

   ```  Feature identifier```:```Relative start position in the original feature```..```Relative end position in the original feature```:```The length of the original feature```
   
    **Example 2 : Extraction of a region specified by the DNA sequence** 
    Extract the region between M13_Forward and M13_reverse. 
   
    ````Python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, "CAGGAAACAGCTATGAC", "TGTAAAACGACGGCCAGT")
   sub_brick.print_dsdna(whole=True)
   sub_brick.print_features() 
    ````
   
    ````
   The region from start 603 to end 706 was extracted.
   CAGGAAACAGCTATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACA
   GTCCTTTGTCGATACTGGTACTAATGCGGTTCGAACGTACGGACGTCCAGCTGAGATCTCCTAGGGGCCCATGGCTCGAGCTTAAGTGACCGGCAGCAAAATGT
   Feature_ID  Label                      Type          Start  End  Strand  
   0           M13/pUC Reverse:20..23:23  primer_bind   0      4    +       
   1           M13 rev                    primer_bind   0      17   +       
   2           M13 Reverse                primer_bind   0      17   +       
   3           source:603..706:2686       source        0      104  +       
   4           lacZ-alpha:1..92:324       CDS           12     104  +       
   5           MCS                        misc_feature  29     86   +       
   6           M13 fwd                    primer_bind   86     103  -       
   7           M13 Forward                primer_bind   86     104  -       
   8           M13/pUC Forward:9..1:23    primer_bind   95     104  -    
    ````
   
   The methods allow the specified sequences to include an extra sequence that doesn't exist in the sequence of the Dbrick object. However, the specified sequence should include unique seqeunce in the sequence of Dbrick object. 
   
   ````Python
   #Source code (continued from previous one)#
   sub_brick = substr(brick, "TTTTTCAGGAAACAGCTATGAC", "GGGGGTGTAAAACGACGGCCAGT")
   sub_brick.print_dsdna(whole=True)
   ````
   
   `````
   The region from start 602 to end 706 was extracted. Adapter sequneces were detected at both ends. Right redundant sequence is TTTTT, Left redundant sequence is GGGGG.
   TTTTTCAGGAAACAGCTATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACACCCCC
   AAAAAGTCCTTTGTCGATACTGGTACTAATGCGGTTCGAACGTACGGACGTCCAGCTGAGATCTCCTAGGGGCCCATGGCTCGAGCTTAAGTGACCGGCAGCAAAATGTGGGGG
   `````
   
    **Example 3 : Extraction of a specific area with flanking region**  
    Extract the region of the MCS with 100 bp flanking length. 
   
    ````Python
    #Source code (continued from previous one)#
    sub_brick = substr(brick, -100, +100, target="11")
    sub_brick.print_features()
    ````
   
    ````
    The region from start 531 to end 788 was extracted.
    Feature_ID  Label                  Type          Start  End  Strand  
    0           source:532..788:2686   source        0      257  +       
    1           lac promoter           promoter      9      40   +       
    2           lac operator           protein_bind  47     64   +       
    3           M13/pUC Reverse        primer_bind   52     75   +       
    4           M13 rev                primer_bind   71     88   +       
    5           M13 Reverse            primer_bind   71     88   +       
    6           lacZ-alpha:1..174:324  CDS           83     257  +       
    7           MCS                    misc_feature  100    157  +       
    8           M13 fwd                primer_bind   157    174  -       
    9           M13 Forward            primer_bind   157    175  -       
    10          M13/pUC Forward        primer_bind   166    189  -   
    ````
   
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

