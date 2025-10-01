import os 
import sys 
import copy 
import regex as re
from types import SimpleNamespace
sys.path.append("/".join(__file__.split("/")[:-1]))
from qobj import QUEEN 

_lib = {
      "FLAG": {
        "amino_acid": "DYKDDDDK",
        "dna_human": "GACTACAAGGACGATGACGATAAG",
        "dna_mouse": "GACTACAAGGACGATGACGATAAG",
        "dna_yeast": "GATTATAAAGATGACGATGACAAA",
        "dna_common": "GACTACAAAGACGATGACGACAAG"
      },
      "3xFLAG": {
        "amino_acid": "DYKDHDGDYKDHDIDYKDDDDK",
        "dna_human": "GACTACAAGGACCACGACGGCGACTACAAGGACCACGACATCGACTACAAGGACGATGACGATAAG",
        "dna_mouse": "GACTACAAGGACCACGACGGCGACTACAAGGACCACGACATCGACTACAAGGACGATGACGATAAG",
        "dna_yeast": "GATTATAAAGATCATGATGGTGATTATAAAGATCATGATATTGATTATAAAGATGACGATGACAAA",
        "dna_common": "GACTACAAGGACCACGACGGTGACTACAAGGACCACGACATCGACTACAAGGACGACGACGACAAG"
      },
      "6xHis": {
        "amino_acid": "HHHHHH",
        "dna_human": "CACCATCACCATCACCAT",
        "dna_mouse": "CACCATCACCATCACCAT",
        "dna_yeast": "CATCACCATCACCATCAC",
        "dna_common": "CATCACCATCACCATCAC"
      },
      "10xHis": {
        "amino_acid": "HHHHHHHHHH",
        "dna_human": "CACCATCACCATCACCATCACCATCACCAT",
        "dna_mouse": "CACCATCACCATCACCATCACCATCACCAT",
        "dna_yeast": "CATCACCATCACCATCACCATCACCATCAC",
        "dna_common": "CATCACCATCACCATCACCATCACCATCAC"
      },
      "HA": {
        "amino_acid": "YPYDVPDYA",
        "dna_human": "TACCCCTACGACGTGCCCGACTACGCC",
        "dna_mouse": "TACCCCTACGACGTGCCCGACTACGCC",
        "dna_yeast": "TATCCTTATGATGTTCCTGATTATGCT",
        "dna_common": "TACCCATACGATGTTCCAGATTACGCT"
      },
      "3xHA": {
        "amino_acid": "YPYDVPDYAYPYDVPDYAYPYDVPDYA",
        "dna_human": "TACCCCTACGACGTGCCCGACTACGCCTACCCCTACGACGTGCCCGACTACGCCTACCCCTACGACGTGCCCGACTACGCC",
        "dna_mouse": "TACCCCTACGACGTGCCCGACTACGCCTACCCCTACGACGTGCCCGACTACGCCTACCCCTACGACGTGCCCGACTACGCC",
        "dna_yeast": "TATCCTTATGATGTTCCTGATTATGCTTATCCTTATGATGTTCCTGATTATGCTTATCCTTATGATGTTCCTGATTATGCT",
        "dna_common": "TACCCATACGATGTTCCAGATTACGCTTACCCATACGATGTTCCAGATTACGCTTACCCATACGATGTTCCAGATTACGCT"
      },
      "Myc": {
        "amino_acid": "EQKLISEEDL",
        "dna_human": "GAGCAGAAGCTGATCAGCGAGGAAGACCTG",
        "dna_mouse": "GAGCAGAAGCTGATCAGCGAGGAAGACCTG",
        "dna_yeast": "GAACAAAAATTGATTTCTGAAGAGGATTTG",
        "dna_common": "GAACAAAAACTCATCTCAGAAGAGGATCTG"
      },
      "V5": {
        "amino_acid": "GKPIPNPLLGLDST",
        "dna_human": "GGCAAGCCCATCCCCAACCCCCTGCTCGGCCTGGACAGCACC",
        "dna_mouse": "GGCAAGCCCATCCCCAACCCCCTGCTCGGCCTGGACAGCACC",
        "dna_yeast": "GGTAAACCTATTCCTAATCCTTTGTTAGGTTTGGATTCTACT",
        "dna_common": "GGCAAGCCCATCCCCAACCCCCTGCTCGGCCTGGACAGCACC"
      },
      "StrepII": {
        "amino_acid": "WSHPQFEK",
        "dna_human": "TGGAGCCACCCCCAGTTCGAGAAG",
        "dna_mouse": "TGGAGCCACCCCCAGTTCGAGAAG",
        "dna_yeast": "TGGTCTCATCCTCAATTTGAAAAA",
        "dna_common": "TGGAGCCACCCCCAGTTCGAGAAG"
      },
      "TwinStrep": {
        "amino_acid": "WSHPQFEKGGGSGGGSWSHPQFEK",
        "dna_human": "TGGAGCCACCCCCAGTTCGAGAAGGGCGGTGGAAGCGGCGGTGGAAGCTGGAGCCACCCCCAGTTCGAGAAG",
        "dna_mouse": "TGGAGCCACCCCCAGTTCGAGAAGGGCGGTGGAAGCGGCGGTGGAAGCTGGAGCCACCCCCAGTTCGAGAAG",
        "dna_yeast": "TGGTCTCATCCTCAATTTGAAAAAGGTGGAGGCTCTGGTGGAGGCTCTTGGTCTCATCCTCAATTTGAAAAA",
        "dna_common": "TGGAGCCACCCCCAGTTCGAGAAGGGCGGTGGAAGCGGCGGTGGAAGCTGGAGCCACCCCCAGTTCGAGAAG"
      },
      "T7": {
        "amino_acid": "MASMTGGQQMG",
        "dna_human": "ATGGCCAGCATGACCGGCGGTCAGCAAATGGGC",
        "dna_mouse": "ATGGCCAGCATGACCGGCGGTCAGCAAATGGGC",
        "dna_yeast": "ATGGCTTCTATGACTGGTGGACAACAGATGGGT",
        "dna_common": "ATGGCCAGCATGACCGGCGGTCAGCAAATGGGC"
      },
      "S-tag": {
        "amino_acid": "KETAAAKFERQHMDS",
        "dna_human": "AAGGAGACCGCCGCTGCGAAGTTCGAGCGCCAGCACATGGACAGC",
        "dna_mouse": "AAGGAGACCGCCGCTGCGAAGTTCGAGCGCCAGCACATGGACAGC",
        "dna_yeast": "AAAGAAACTGCTGCCGCAAAATTTGAAAGACAACATATGGATTCT",
        "dna_common": "AAGGAGACCGCCGCTGCGAAGTTCGAGCGCCAGCACATGGACAGC"
      },
      "AviTag": {
        "amino_acid": "GLNDIFEAQKIEWHE",
        "dna_human": "GGCCTGAACGACATCTTCGAGGCCCAGAAGATCGAGTGGCACGAG",
        "dna_mouse": "GGCCTGAACGACATCTTCGAGGCCCAGAAGATCGAGTGGCACGAG",
        "dna_yeast": "GGTTTGAATGATATTTTTGAAGCTCAAAAAATTGAATGGCATGAA",
        "dna_common": "GGCCTGAACGACATCTTCGAGGCCCAGAAGATCGAGTGGCACGAG"
      }
    } 
lib = copy.deepcopy(_lib) 
for key1 in _lib:
    for key2 in ["dna_human", "dna_mouse", "dna_yeast", "dna_common"]:
        lib[key1][key2] = QUEEN(seq=_lib[key1][key2], supfeature={"feature_type": "CDS", "qualifier:label":key1, "qualifier:translation":_lib[key1]["amino_acid"]}, project=f"{key1.replace(' ','_')}_{key2}")
    lib[key1] = SimpleNamespace(**lib[key1])
