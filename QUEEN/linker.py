import os 
import sys 
import copy 
import regex as re
from types import SimpleNamespace
sys.path.append("/".join(__file__.split("/")[:-1]))
from qobj import QUEEN 

_lib = { 
          "G4S": {
            "amino_acid": "GGGGS",
            "dna_human": "GGCGGTGGAGGGAGC",
            "dna_mouse": "GGCGGTGGAGGGAGC",
            "dna_yeast": "GGTGGAGGCGGGTCT",
            "dna_common": "GGTGGTGGTGGTTCT"
          },
          "GGGGS": {
            "amino_acid": "GGGGS",
            "dna_human": "GGCGGTGGAGGGAGC",
            "dna_mouse": "GGCGGTGGAGGGAGC",
            "dna_yeast": "GGTGGAGGCGGGTCT",
            "dna_common": "GGTGGTGGTGGTTCT"
          },
          "G4Sx2": {
            "amino_acid": "GGGGSGGGGS",
            "dna_human": "GGCGGTGGAGGGAGCGGCGGTGGAGGGAGC",
            "dna_mouse": "GGCGGTGGAGGGAGCGGCGGTGGAGGGAGC",
            "dna_yeast": "GGTGGAGGCGGGTCTGGTGGAGGCGGGTCT",
            "dna_common": "GGTGGTGGTGGTTCTGGTGGTGGTGGTTCT"
          },
          "EAAAK": {
            "amino_acid": "EAAAK",
            "dna_human": "GAGGCCGCTGCGAAG",
            "dna_mouse": "GAGGCCGCTGCGAAG",
            "dna_yeast": "GAAGCTGCCGCAAAA",
            "dna_common": "GAAGCTGCTGCTAAA"
          },
          "EAAAKx2": {
            "amino_acid": "EAAAKEAAAK",
            "dna_human": "GAGGCCGCTGCGAAGGAGGCCGCTGCGAAG",
            "dna_mouse": "GAGGCCGCTGCGAAGGAGGCCGCTGCGAAG",
            "dna_yeast": "GAAGCTGCCGCAAAAGAAGCTGCCGCAAAA",
            "dna_common": "GAAGCTGCTGCTAAAGAAGCTGCTGCTAAA"
          },
          "TEV site": {
            "amino_acid": "ENLYFQG",
            "dna_human": "GAGAACCTGTACTTCCAGGGC",
            "dna_mouse": "GAGAACCTGTACTTCCAGGGC",
            "dna_yeast": "GAAAATTTGTATTTTCAAGGT",
            "dna_common": "GAAAACCTGTATTTTCAGGGT"
          },
          "HRV3C site": {
            "amino_acid": "LEVLFQGP",
            "dna_human": "CTGGAGGTGCTGTTCCAGGGCCCC",
            "dna_mouse": "CTGGAGGTGCTGTTCCAGGGCCCC",
            "dna_yeast": "TTGGAAGTTTTGTTTCAAGGTCCT",
            "dna_common": "CTGGAGGTGCTGTTTCAGGGTCCA"
          },
          "thrombin site": {
            "amino_acid": "LVPRGS",
            "dna_human": "CTGGTGCCCCGCGGCAGC",
            "dna_mouse": "CTGGTGCCCCGCGGCAGC",
            "dna_yeast": "TTGGTTCCTAGAGGTTCT",
            "dna_common": "CTGGTGCCTCGTGGTTCT"
          },
          "Factor Xa site": {
            "amino_acid": "IEGR",
            "dna_human": "ATCGAGGGCCGC",
            "dna_mouse": "ATCGAGGGCCGC",
            "dna_yeast": "ATTGAAGGTAGA",
            "dna_common": "ATCGAGGGTCGT"
          }
    }
lib = copy.deepcopy(_lib) 
for key1 in _lib:
    for key2 in ["dna_human", "dna_mouse", "dna_yeast", "dna_common"]:
        lib[key1][key2] = QUEEN(seq=_lib[key1][key2], supfeature={"feature_type": "CDS", "qualifier:label":key1, "qualifier:translation":_lib[key1]["amino_acid"]}, project=f"{key1.replace(' ','_')}_{key2}")
    lib[key1] = SimpleNamespace(**lib[key1])
