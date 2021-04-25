import sys
sys.path.append("../../")
from dna import *
DNA.dna_dict['pCMV_ABE'] = DNA(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'] = DNA(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_0'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='3914/3914', end='6442/6442', project='pCMV-nCas-PmCDA1-ugi', process_description=description0)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_0'], left='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', right='TCTCCGACTACGACGTGGATCATATCGTGCCCCAG', project='pCMV-nCas-PmCDA1-ugi', process_description=description0)
editdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f1'), project=None, new_copy=False, process_description=description0)
editdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'], query='f1', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description0)
editdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'], query='f1', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-1'), project=None, new_copy=False, process_description=description0)

description1 = 'The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046'
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_2'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='6484/6484', end='9244/9244', project='pCMV-nCas-PmCDA1-ugi', process_description=description1)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_3'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_2'], left='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', right='GAGAGAACAAGATCAAGATGCTATAATGAGTTTAAA', project='pCMV-nCas-PmCDA1-ugi', process_description=description1)
editdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_3'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f2'), project=None, new_copy=False, process_description=description1)
editdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_3'], query='f2', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description1)
editdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_3'], query='f2', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-2'), project=None, new_copy=False, process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048'
DNA.dna_dict['pCMV_ABE_0'] = cropdna(DNA.dna_dict['pCMV_ABE'], start='5277/5277', end='8613/8613', project='pCMV_ABE', process_description=description2)
DNA.dna_dict['pCMV_ABE_1'] = modifyends(DNA.dna_dict['pCMV_ABE_0'], left='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', right='CTATAGGGAGAGCCGCCACCATGGCACCGAAG', project='pCMV_ABE', process_description=description2)
editdna(DNA.dna_dict['pCMV_ABE_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f3'), project=None, new_copy=False, process_description=description2)
editdna(DNA.dna_dict['pCMV_ABE_1'], query='f3', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description2)
editdna(DNA.dna_dict['pCMV_ABE_1'], query='f3', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-3'), project=None, new_copy=False, process_description=description2)

description3 = 'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling the two insert fragments and the backbone fragment.'
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_4'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'], left='*{25}/-{25}', right='-{28}/*{28}', project='pCMV-nCas-PmCDA1-ugi', process_description=description3)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_5'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_3'], left='*{28}/-{28}', right='-{25}/*{25}', project='pCMV-nCas-PmCDA1-ugi', process_description=description3)
DNA.dna_dict['pCMV_ABE_2'] = modifyends(DNA.dna_dict['pCMV_ABE_1'], left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV_ABE', process_description=description3)
DNA.dna_dict['pCMV-Target-AID'] = joindna(*[DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_4'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_5'], DNA.dna_dict['pCMV_ABE_2']], topology='circular', project='pCMV-Target-AID', process_description=description3)
DNA.dna_dict['pCMV-Target-AID_0'], = cutdna(DNA.dna_dict['pCMV-Target-AID'], '64/64', crop=False, project='pCMV-Target-AID', process_description=description3)
DNA.dna_dict['pCMV-Target-AID_0'] = joindna(*[DNA.dna_dict['pCMV-Target-AID_0']], topology='circular', project='pCMV-Target-AID', process_description=description3)
DNA.dna_dict['pCMV-Target-AID_0'].writedna('reconstructed_pCMV-Target-AID.gbk')
