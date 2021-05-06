import sys
sys.path.append("../../")
from dnaquine import *
DNA.dna_dict['pCMV_ABE'] = DNA(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'] = DNA(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'
DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_0'] = DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'].finddna(query='GGCACCGAAGAAGAAGCG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description0)
DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_1'] = DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'].finddna(query='CCACGTCGTAGTCGGAGA', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description0)

description1 = 'The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046'
DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_2'] = DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'].finddna(query='TCGTGCCCCAGTCTTTTC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description1)
DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_3'] = DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'].finddna(query='TCTTGATCTTGTTCTCTC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048'
DNA.queried_features_dict['pCMV_ABE_4'] = DNA.dna_dict['pCMV_ABE'].finddna(query='AGTTTAAACCCGCTGATC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description2)
DNA.queried_features_dict['pCMV_ABE_5'] = DNA.dna_dict['pCMV_ABE'].finddna(query='TGGCGGCTCTCCCTATAG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description2)

description3 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS052'
DNA.queried_features_dict['pCMV_ABE_7'] = DNA.dna_dict['pCMV_ABE'].finddna(query='AGTTTAAACCCGCTGATC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description3)
DNA.queried_features_dict['pCMV_ABE_8'] = DNA.dna_dict['pCMV_ABE'].finddna(query='CACCCCCAAGCTGTGACA', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description3)
DNA.queried_feature_dict['pCMV_ABE_7_q0'] = DNA.queried_features_dict['pCMV_ABE_7'][0]
DNA.queried_feature_dict['pCMV_ABE_8_q0'] = DNA.queried_features_dict['pCMV_ABE_8'][0]
DNA.dna_dict['pCMV_ABE_3'] = cropdna(DNA.dna_dict['pCMV_ABE'], start=DNA.queried_feature_dict['pCMV_ABE_7_q0'].end, end=DNA.queried_feature_dict['pCMV_ABE_8_q0'].start, project='pCMV_ABE', process_description=description3)
DNA.dna_dict['pCMV_ABE_4'] = modifyends(DNA.dna_dict['pCMV_ABE_3'], left='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', right='TGTCACAGCTTGGGGGTGACGGTGGAGGAGGT', project='pCMV_ABE', process_description=description3)
editfeature(DNA.dna_dict['pCMV_ABE_4'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f4'), project=None, new_copy=False, process_description=description3)
editfeature(DNA.dna_dict['pCMV_ABE_4'], query='f4', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description3)
editfeature(DNA.dna_dict['pCMV_ABE_4'], query='f4', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-4'), project=None, new_copy=False, process_description=description3)

description4 = 'The insert fragment encoding the C-terminus region of Target-AID was amplified from pcDNA-pCMV-nCas9 using RS051/RS046.'
DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_9'] = DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'].finddna(query='TGGAGGAGGTACCGGCGG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description4)
DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_10'] = DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'].finddna(query='TCTTGATCTTGTTCTCTC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description4)
DNA.queried_feature_dict['pCMV-nCas-PmCDA1-ugi_9_q0'] = DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_9'][0]
DNA.queried_feature_dict['pCMV-nCas-PmCDA1-ugi_10_q0'] = DNA.queried_features_dict['pCMV-nCas-PmCDA1-ugi_10'][0]
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_6'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start=DNA.queried_feature_dict['pCMV-nCas-PmCDA1-ugi_9_q0'].end, end=DNA.queried_feature_dict['pCMV-nCas-PmCDA1-ugi_10_q0'].start, project='pCMV-nCas-PmCDA1-ugi', process_description=description4)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_6'], left='GCTTGGGGGTGACGGTGGAGGAGGTACCGGCGG', right='GAGAGAACAAGATCAAGATGCTATAATGAGTTTAAA', project='pCMV-nCas-PmCDA1-ugi', process_description=description4)
editfeature(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f5'), project=None, new_copy=False, process_description=description4)
editfeature(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'], query='f5', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description4)
editfeature(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'], query='f5', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-5'), project=None, new_copy=False, process_description=description4)

description5 = 'The Target-ACE plasmid (pCMV-Target-ACE) was constructed by assembling the insert fragment and a backbone fragment.'
DNA.dna_dict['pCMV_ABE_5'] = modifyends(DNA.dna_dict['pCMV_ABE_4'], left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV_ABE', process_description=description5)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_8'] = modifyends(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_7'], left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV-nCas-PmCDA1-ugi', process_description=description5)
DNA.dna_dict['pCMV-Target-ACE'] = joindna(*[DNA.dna_dict['pCMV_ABE_5'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_8']], topology='circular', project='pCMV-Target-ACE', process_description=description5)
DNA.queried_features_dict['pCMV-Target-ACE_11'] = DNA.dna_dict['pCMV-Target-ACE'].finddna(query='Cas9(D10A)', key_attribute='qualifier:label', min_match=None, max_mismatch=0, process_description=description5)
DNA.dna_dict['pCMV-Target-ACE_0'], = cutdna(DNA.dna_dict['pCMV-Target-ACE'], '4582/4582', crop=False, project='pCMV-Target-ACE', process_description=description5)
DNA.dna_dict['pCMV-Target-ACE_0'] = joindna(*[DNA.dna_dict['pCMV-Target-ACE_0']], topology='circular', project='pCMV-Target-ACE', process_description=description5)
DNA.dna_dict['pCMV-Target-ACE_0'].writedna('reconstructed_pCMV-Target-ACE.gbk')
