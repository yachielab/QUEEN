import sys
sys.path.append("../../")
from dnaquine import *

genome = DNA(record="GCF_000005845.2_ASM584v2_genomic.gbff") 
print(genome)
sub_genome = cropdna(genome, 0, 1000000)
print(sub_genome) 
cds_list = genome.searchdna(key_attribute="feature_type", query="CDS")
genome.printfeature(feature_list=cds_list[0:10]) 
