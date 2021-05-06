import sys
sys.path.append("../../")
from dnaquine import *
DNA.dna_dict['pCMV-ABEmax'] = DNA(seq=None, record='input/addgene_112095.gbk', project='pCMV-ABEmax', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-BE4max'] = DNA(seq=None, record='input/addgene_112093.gbk', project='pCMV-BE4max', topology='linear', format=None, process_description=None)
DNA.dna_dict['opt-pmCDA1-ugi'] = DNA(seq=None, record='input/puc-optimized-pmcda1-ugi.gbk', project='opt-pmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The fragment encoding a nCas9 fragment was amplified from pCMV-BE4max using SI1308/SI1305.'
DNA.queried_features_dict['pCMV-BE4max_12'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='GAAGTACAGCATCGGCCT', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description0)
DNA.queried_features_dict['pCMV-BE4max_13'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='CCCAGCTGAGACAGGTCG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description0)

description1 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax was amplified with primer pair SI1304/SI1307.'
DNA.queried_features_dict['opt-pmCDA1-ugi_14'] = DNA.dna_dict['opt-pmCDA1-ugi'].finddna(query='AGGAGGAACTGGAGGAGG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description1)
DNA.queried_features_dict['opt-pmCDA1-ugi_15'] = DNA.dna_dict['opt-pmCDA1-ugi'].finddna(query='CATTTTGATTTTGTTCTC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description1)
DNA.queried_feature_dict['opt-pmCDA1-ugi_14_q0'] = DNA.queried_features_dict['opt-pmCDA1-ugi_14'][0]
DNA.queried_feature_dict['opt-pmCDA1-ugi_15_q0'] = DNA.queried_features_dict['opt-pmCDA1-ugi_15'][0]
DNA.dna_dict['opt-pmCDA1-ugi_0'] = cropdna(DNA.dna_dict['opt-pmCDA1-ugi'], start=DNA.queried_feature_dict['opt-pmCDA1-ugi_14_q0'].end, end=DNA.queried_feature_dict['opt-pmCDA1-ugi_15_q0'].start, project='opt-pmCDA1-ugi', process_description=description1)
DNA.dna_dict['opt-pmCDA1-ugi_1'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_0'], left='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', right='GAGAACAAAATCAAAATGCTGTGATAATGAGTTTAAACCCGCTGA', project='opt-pmCDA1-ugi', process_description=description1)
editfeature(DNA.dna_dict['opt-pmCDA1-ugi_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f7'), project=None, new_copy=False, process_description=description1)
editfeature(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description1)
editfeature(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-7'), project=None, new_copy=False, process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABEmax using SI1310/SI1309.'
DNA.queried_features_dict['pCMV-ABEmax_16'] = DNA.dna_dict['pCMV-ABEmax'].finddna(query='TGAGTTTAAACCCGCTGA', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description2)
DNA.queried_features_dict['pCMV-ABEmax_17'] = DNA.dna_dict['pCMV-ABEmax'].finddna(query='TTCTTCTTTGGTGACTCG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description2)
DNA.queried_feature_dict['pCMV-ABEmax_16_q0'] = DNA.queried_features_dict['pCMV-ABEmax_16'][0]
DNA.queried_feature_dict['pCMV-ABEmax_17_q0'] = DNA.queried_features_dict['pCMV-ABEmax_17'][0]
DNA.dna_dict['pCMV-ABEmax_0'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start=DNA.queried_feature_dict['pCMV-ABEmax_16_q0'].end, end=DNA.queried_feature_dict['pCMV-ABEmax_17_q0'].start, project='pCMV-ABEmax', process_description=description2)
DNA.dna_dict['pCMV-ABEmax_1'] = modifyends(DNA.dna_dict['pCMV-ABEmax_0'], left='TGATAATGAGTTTAAACCCGCTGA', right='CGAGTCACCAAAGAAGAAGCGGAAAGTC', project='pCMV-ABEmax', process_description=description2)
editfeature(DNA.dna_dict['pCMV-ABEmax_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f8'), project=None, new_copy=False, process_description=description2)
editfeature(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description2)
editfeature(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-8'), project=None, new_copy=False, process_description=description2)

description3 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling the two insert fragments and the backbone fragment.'
DNA.dna_dict['opt-pmCDA1-ugi_2'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_1'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description3)
DNA.dna_dict['pCMV-ABEmax_2'] = modifyends(DNA.dna_dict['pCMV-ABEmax_1'], left='*{24}/-{24}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description3)

description4 = 'The rAPOBEC1 fragment was obtained from BE4max using SI1352/SI1357.'
DNA.queried_features_dict['pCMV-BE4max_21'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='GAGACTGGGCCTGTCGCC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description4)
DNA.queried_features_dict['pCMV-BE4max_22'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='GCCTGTGGCCCACAGGAT', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description4)
DNA.queried_feature_dict['pCMV-BE4max_21_q0'] = DNA.queried_features_dict['pCMV-BE4max_21'][0]
DNA.queried_feature_dict['pCMV-BE4max_22_q0'] = DNA.queried_features_dict['pCMV-BE4max_22'][0]
DNA.dna_dict['pCMV-BE4max_3'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_feature_dict['pCMV-BE4max_21_q0'].end, end=DNA.queried_feature_dict['pCMV-BE4max_22_q0'].start, project='pCMV-BE4max', process_description=description4)
DNA.dna_dict['pCMV-BE4max_4'] = modifyends(DNA.dna_dict['pCMV-BE4max_3'], left='GATGACGATGACAAGTCTGGCTCCTCAGAGACTGGGCCTGTCGCC', right='ATCCTGTGGGCCACAGGCCTGAAG', project='pCMV-BE4max', process_description=description4)
editfeature(DNA.dna_dict['pCMV-BE4max_4'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f10'), project=None, new_copy=False, process_description=description4)
editfeature(DNA.dna_dict['pCMV-BE4max_4'], query='f10', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description4)
editfeature(DNA.dna_dict['pCMV-BE4max_4'], query='f10', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-10'), project=None, new_copy=False, process_description=description4)

description5 = 'The 2×UGI fragment was obtained from BE4max using SI1359/SI1350.'
DNA.queried_features_dict['pCMV-BE4max_23'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='AATCTGAGCGACATCATT', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description5)
DNA.queried_features_dict['pCMV-BE4max_24'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='ACTTTCCTCTTCTTCTTG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description5)
DNA.queried_feature_dict['pCMV-BE4max_23_q0'] = DNA.queried_features_dict['pCMV-BE4max_23'][0]
DNA.queried_feature_dict['pCMV-BE4max_24_q0'] = DNA.queried_features_dict['pCMV-BE4max_24'][0]
DNA.dna_dict['pCMV-BE4max_5'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_feature_dict['pCMV-BE4max_23_q0'].end, end=DNA.queried_feature_dict['pCMV-BE4max_24_q0'].start, project='pCMV-BE4max', process_description=description5)
DNA.dna_dict['pCMV-BE4max_6'] = modifyends(DNA.dna_dict['pCMV-BE4max_5'], left='ATCCTGTGGGCCACAGGCCTGAAGACTAATCTGAGCGACATCATT', right='CAAGAAGAAGAGGAAAGTCTAATAATGAGTTTAAACCCGCTGATC', project='pCMV-BE4max', process_description=description5)
editfeature(DNA.dna_dict['pCMV-BE4max_6'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f11'), project=None, new_copy=False, process_description=description5)
editfeature(DNA.dna_dict['pCMV-BE4max_6'], query='f11', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description5)
editfeature(DNA.dna_dict['pCMV-BE4max_6'], query='f11', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-11'), project=None, new_copy=False, process_description=description5)

description6 = 'The backbone fragment was obtained from BE4max using SI1351/SI448.'
DNA.queried_features_dict['pCMV-BE4max_25'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='AGTTTAAACCCGCTGATC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description6)
DNA.queried_features_dict['pCMV-BE4max_26'] = DNA.dna_dict['pCMV-BE4max'].finddna(query='AAAGTTCTGCTATGTGGC', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description6)
DNA.queried_feature_dict['pCMV-BE4max_25_q0'] = DNA.queried_features_dict['pCMV-BE4max_25'][0]
DNA.queried_feature_dict['pCMV-BE4max_26_q0'] = DNA.queried_features_dict['pCMV-BE4max_26'][0]
DNA.dna_dict['pCMV-BE4max_7'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_feature_dict['pCMV-BE4max_25_q0'].end, end=DNA.queried_feature_dict['pCMV-BE4max_26_q0'].start, project='pCMV-BE4max', process_description=description6)
DNA.dna_dict['pCMV-BE4max_8'] = modifyends(DNA.dna_dict['pCMV-BE4max_7'], left='TAATAATGAGTTTAAACCCGCTGATC', right='GCCACATAGCAGAACTTTAAAAGTG', project='pCMV-BE4max', process_description=description6)
editfeature(DNA.dna_dict['pCMV-BE4max_8'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f12'), project=None, new_copy=False, process_description=description6)
editfeature(DNA.dna_dict['pCMV-BE4max_8'], query='f12', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description6)
editfeature(DNA.dna_dict['pCMV-BE4max_8'], query='f12', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-12'), project=None, new_copy=False, process_description=description6)
DNA.dna_dict['pCMV-BE4max_9'] = modifyends(DNA.dna_dict['pCMV-BE4max_4'], left='*{15}/-{15}', right='-{24}/*{24}', project='pCMV-BE4max', process_description=description6)
DNA.dna_dict['pCMV-BE4max_10'] = modifyends(DNA.dna_dict['pCMV-BE4max_6'], left='*{24}/-{24}', right='-{26}/*{26}', project='pCMV-BE4max', process_description=description6)
DNA.dna_dict['pCMV-BE4max_11'] = modifyends(DNA.dna_dict['pCMV-BE4max_8'], left='*{26}/-{26}', right='-{25}/*{25}', project='pCMV-BE4max', process_description=description6)

description7 = 'An ABEmax fragment obtained from pCMV-ABEmax using SI945/SI1305'
DNA.queried_features_dict['pCMV-ABEmax_28'] = DNA.dna_dict['pCMV-ABEmax'].finddna(query='AATACGACTCACTATAGG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description7)
DNA.queried_features_dict['pCMV-ABEmax_29'] = DNA.dna_dict['pCMV-ABEmax'].finddna(query='CCCAGCTGAGACAGGTCG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description7)
DNA.queried_feature_dict['pCMV-ABEmax_28_q0'] = DNA.queried_features_dict['pCMV-ABEmax_28'][0]
DNA.queried_feature_dict['pCMV-ABEmax_29_q0'] = DNA.queried_features_dict['pCMV-ABEmax_29'][0]
DNA.dna_dict['pCMV-ABEmax_3'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start=DNA.queried_feature_dict['pCMV-ABEmax_28_q0'].end, end=DNA.queried_feature_dict['pCMV-ABEmax_29_q0'].start, project='pCMV-ABEmax', process_description=description7)
DNA.dna_dict['pCMV-ABEmax_4'] = modifyends(DNA.dna_dict['pCMV-ABEmax_3'], left='AGATCCGCGGCCGCTAATACGACTCACTATAGG', right='CGACCTGTCTCAGCTGGGAGGTGAC', project='pCMV-ABEmax', process_description=description7)
editfeature(DNA.dna_dict['pCMV-ABEmax_4'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f13'), project=None, new_copy=False, process_description=description7)
editfeature(DNA.dna_dict['pCMV-ABEmax_4'], query='f13', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description7)
editfeature(DNA.dna_dict['pCMV-ABEmax_4'], query='f13', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-13'), project=None, new_copy=False, process_description=description7)

description8 = 'The Target-ACEmax plasmid (pCMV-BE4max(C)) was constructed by assembling a insert fragment and two backbone fragments.'
DNA.dna_dict['opt-pmCDA1-ugi_3'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_2'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description8)
DNA.dna_dict['pCMV-ABEmax_5'] = modifyends(DNA.dna_dict['pCMV-ABEmax_2'], left='*{24}/-{24}', right='-{103}/*{103}', project='pCMV-ABEmax', process_description=description8)
DNA.dna_dict['pCMV-ABEmax_6'] = modifyends(DNA.dna_dict['pCMV-ABEmax_4'], left='*{103}/-{103}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description8)
DNA.dna_dict['pCMV-Target-ACEmax'] = joindna(*[DNA.dna_dict['opt-pmCDA1-ugi_3'], DNA.dna_dict['pCMV-ABEmax_5'], DNA.dna_dict['pCMV-ABEmax_6']], topology='circular', project='pCMV-Target-ACEmax', process_description=description8)
DNA.queried_features_dict['pCMV-Target-ACEmax_30'] = DNA.dna_dict['pCMV-Target-ACEmax'].finddna(query='Cas9(D10A)', key_attribute='qualifier:label', min_match=None, max_mismatch=0, process_description=description8)
DNA.dna_dict['pCMV-Target-ACEmax_0'], = cutdna(DNA.dna_dict['pCMV-Target-ACEmax'], '5834/5834', crop=False, project='pCMV-Target-ACEmax', process_description=description8)
DNA.dna_dict['pCMV-Target-ACEmax_0'] = joindna(*[DNA.dna_dict['pCMV-Target-ACEmax_0']], topology='circular', project='pCMV-Target-ACEmax', process_description=description8)

description9 = 'An ABEmax fragment obtained from pCMV-Target-ACEmax using SI447/SI1105.'
DNA.queried_features_dict['pCMV-Target-ACEmax_0_31'] = DNA.dna_dict['pCMV-Target-ACEmax_0'].finddna(query='AGCAGAACTTTAAAAGTG', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description9)
DNA.queried_features_dict['pCMV-Target-ACEmax_0_32'] = DNA.dna_dict['pCMV-Target-ACEmax_0'].finddna(query='GTCATCGTCATCCTTGTA', key_attribute='sequence', min_match=None, max_mismatch=0, process_description=description9)
DNA.queried_feature_dict['pCMV-Target-ACEmax_0_31_q0'] = DNA.queried_features_dict['pCMV-Target-ACEmax_0_31'][0]
DNA.queried_feature_dict['pCMV-Target-ACEmax_0_32_q0'] = DNA.queried_features_dict['pCMV-Target-ACEmax_0_32'][0]
DNA.dna_dict['pCMV-Target-ACEmax_1'] = cropdna(DNA.dna_dict['pCMV-Target-ACEmax_0'], start=DNA.queried_feature_dict['pCMV-Target-ACEmax_0_31_q0'].end, end=DNA.queried_feature_dict['pCMV-Target-ACEmax_0_32_q0'].start, project='pCMV-Target-ACEmax', process_description=description9)
DNA.dna_dict['pCMV-Target-ACEmax_2'] = modifyends(DNA.dna_dict['pCMV-Target-ACEmax_1'], left='GCCACATAGCAGAACTTTAAAAGTG', right='TACAAGGATGACGATGACAAG', project='pCMV-Target-ACEmax', process_description=description9)
editfeature(DNA.dna_dict['pCMV-Target-ACEmax_2'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f14'), project=None, new_copy=False, process_description=description9)
editfeature(DNA.dna_dict['pCMV-Target-ACEmax_2'], query='f14', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description9)
editfeature(DNA.dna_dict['pCMV-Target-ACEmax_2'], query='f14', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-14'), project=None, new_copy=False, process_description=description9)
DNA.dna_dict['pCMV-Target-ACEmax_3'] = modifyends(DNA.dna_dict['pCMV-Target-ACEmax_2'], left='*{25}/-{25}', right='-{15}/*{15}', project='pCMV-Target-ACEmax', process_description=description9)
DNA.dna_dict['pCMV-BE4max_12'] = modifyends(DNA.dna_dict['pCMV-BE4max_9'], left='*{15}/-{15}', right='-{24}/*{24}', project='pCMV-BE4max', process_description=description9)
DNA.dna_dict['pCMV-BE4max_13'] = modifyends(DNA.dna_dict['pCMV-BE4max_10'], left='*{24}/-{24}', right='-{26}/*{26}', project='pCMV-BE4max', process_description=description9)
DNA.dna_dict['pCMV-BE4max_14'] = modifyends(DNA.dna_dict['pCMV-BE4max_11'], left='*{26}/-{26}', right='-{25}/*{25}', project='pCMV-BE4max', process_description=description9)

description10 = 'The ACBEmax plasmid (pCMV-ACBEmax) was constructed by assembling the ABEmax fragment with the three fragments encoding the rAPOBEC1 domain, 2×UGI domain, and the backbone fragment.'
DNA.dna_dict['pCMV-ACBEmax'] = joindna(*[DNA.dna_dict['pCMV-Target-ACEmax_3'], DNA.dna_dict['pCMV-BE4max_12'], DNA.dna_dict['pCMV-BE4max_13'], DNA.dna_dict['pCMV-BE4max_14']], topology='circular', project='pCMV-ACBEmax', process_description=description10)
DNA.queried_features_dict['pCMV-ACBEmax_33'] = DNA.dna_dict['pCMV-ACBEmax'].finddna(query='Cas9(D10A)', key_attribute='qualifier:label', min_match=None, max_mismatch=0, process_description=description10)
DNA.dna_dict['pCMV-ACBEmax_0'], = cutdna(DNA.dna_dict['pCMV-ACBEmax'], '2504/2504', crop=False, project='pCMV-ACBEmax', process_description=description10)
DNA.dna_dict['pCMV-ACBEmax_0'] = joindna(*[DNA.dna_dict['pCMV-ACBEmax_0']], topology='circular', project='pCMV-ACBEmax', process_description=description10)
DNA.dna_dict['pCMV-ACBEmax_0'].writedna('reconstructed_pCMV-ACBEmax.gbk')
