#import sys
#sys.path.append("../../")
from dna import *
DNA.dna_dict['pCMV_ABE'] = DNA(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'] = DNA(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS052'
DNA.dna_dict['pCMV_ABE_3'] = cropdna(DNA.dna_dict['pCMV_ABE'], start='5277/5277', end='5178/5178', project='pCMV_ABE', process_description=description0)
DNA.dna_dict['pCMV_ABE_4'] = modifyends(DNA.dna_dict['pCMV_ABE_3'], left='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', right='TGTCACAGCTTGGGGGTGACGGTGGAGGAGGT', project='pCMV_ABE', process_description=description0)

description1 = 'The insert fragment encoding the C-terminus region of Target-AID was amplified from pcDNA-pCMV-nCas9 using RS051/RS046.'
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_6'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='8102/8102', end='9244/9244', project='pCMV-nCas-PmCDA1-ugi', process_description=description1)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_6'], left='GCTTGGGGGTGACGGTGGAGGAGGTACCGGCGG', right='GAGAGAACAAGATCAAGATGCTATAATGAGTTTAAA', project='pCMV-nCas-PmCDA1-ugi', process_description=description1)

description2 = 'The Target-ACE plasmid (pCMV-Target-ACE) was constructed by assembling the insert fragment and a backbone fragment.'
DNA.dna_dict['pCMV_ABE_5'] = modifyends(DNA.dna_dict['pCMV_ABE_4'], left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV_ABE', process_description=description2)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_8'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'], left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV-nCas-PmCDA1-ugi', process_description=description2)
DNA.dna_dict['pCMV-Target-ACE'] = joindna(*[DNA.dna_dict['pCMV_ABE_5'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_8']], topology='circular', project='pCMV-Target-ACE', process_description=description2)
DNA.dna_dict['pCMV-Target-ACE_0'], = cutdna(DNA.dna_dict['pCMV-Target-ACE'], '4582/4582', crop=False, project='pCMV-Target-ACE', process_description=description2)
DNA.dna_dict['pCMV-Target-ACE_0'] = joindna(*[DNA.dna_dict['pCMV-Target-ACE_0']], topology='circular', project='pCMV-Target-ACE', process_description=description2)
DNA.dna_dict['pCMV-Target-ACE_0'].writedna('reconstructed_pCMV-Target-ACE.gbk')
