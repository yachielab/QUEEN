import sys 
sys.path.append("../../")
from dna import *

#pCMV-Target-AID construction
nCas9_AID_source = DNA(record="input/addgene_79620.gbk")
nCas9_AID_source.name = "addgene_79620"
pCMV_ABE = DNA(record="input/addgene_102919.gbk")
pCMV_ABE.name = "addgene_102919"

##N-terminus of Target-AID
RS045 = DNA(record="input/RS045.fasta")
HM129 = DNA(record="input/HM129.fasta")
fwdna = nCas9_AID_source.finddna(RS045)
rvdna = nCas9_AID_source.finddna(HM129)
print(fwdna[0].__dict__,len(fwdna))
print(rvdna[0].__dict__,len(rvdna))