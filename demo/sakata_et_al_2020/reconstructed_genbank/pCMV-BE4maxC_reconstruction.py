import sys
sys.path.append("../../")
from dnaquine import *
DNA.dna_dict['pCMV-ABEmax'] = DNA(seq=None, record='input/addgene_112095.gbk', project='pCMV-ABEmax', topology='linear', format=None, product=None, process_description=None)
DNA.dna_dict['pCMV-BE4max'] = DNA(seq=None, record='input/addgene_112093.gbk', project='pCMV-BE4max', topology='linear', format=None, product=None, process_description=None)
DNA.dna_dict['opt-pmCDA1-ugi'] = DNA(seq=None, record='input/puc-optimized-pmcda1-ugi.gbk', project='opt-pmCDA1-ugi', topology='linear', format=None, product=None, process_description=None)

description0 = 'The fragment encoding a nCas9 fragment was amplified from pCMV-BE4max using SI1308/SI1305.'
DNA.queried_features_dict['pCMV-BE4max_12'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='GAAGTACAGCATCGGCCT', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description0)
DNA.queried_features_dict['pCMV-BE4max_13'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='CCCAGCTGAGACAGGTCG', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description0)
DNA.dna_dict['pCMV-BE4max_0'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_features_dict['pCMV-BE4max_12'][0].end, end=DNA.queried_features_dict['pCMV-BE4max_13'][0].start, project='pCMV-BE4max', product=None, process_description=description0)
DNA.dna_dict['pCMV-BE4max_1'] = modifyends(DNA.dna_dict['pCMV-BE4max_0'], left='ACCAAAGAAGAAGCGGAAAGTCGACAAGAAGTACAGCATCGGCCT', right='CGACCTGTCTCAGCTGGGAGGTGAC', project='pCMV-BE4max', product=None, process_description=description0)
editfeature(DNA.dna_dict['pCMV-BE4max_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f6'), project=None, new_copy=False, product=None, process_description=description0)
editfeature(DNA.dna_dict['pCMV-BE4max_1'], query='f6', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description0)
editfeature(DNA.dna_dict['pCMV-BE4max_1'], query='f6', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-6'), project=None, new_copy=False, product=None, process_description=description0)

description1 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax was amplified with primer pair SI1304/SI1307.'
DNA.queried_features_dict['opt-pmCDA1-ugi_14'] = DNA.dna_dict['opt-pmCDA1-ugi'].searchdna(query='AGGAGGAACTGGAGGAGG', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description1)
DNA.queried_features_dict['opt-pmCDA1-ugi_15'] = DNA.dna_dict['opt-pmCDA1-ugi'].searchdna(query='CATTTTGATTTTGTTCTC', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description1)
DNA.dna_dict['opt-pmCDA1-ugi_0'] = cropdna(DNA.dna_dict['opt-pmCDA1-ugi'], start=DNA.queried_features_dict['opt-pmCDA1-ugi_14'][0].end, end=DNA.queried_features_dict['opt-pmCDA1-ugi_15'][0].start, project='opt-pmCDA1-ugi', product=None, process_description=description1)
DNA.dna_dict['opt-pmCDA1-ugi_1'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_0'], left='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', right='GAGAACAAAATCAAAATGCTGTGATAATGAGTTTAAACCCGCTGA', project='opt-pmCDA1-ugi', product=None, process_description=description1)
editfeature(DNA.dna_dict['opt-pmCDA1-ugi_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f7'), project=None, new_copy=False, product=None, process_description=description1)
editfeature(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description1)
editfeature(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-7'), project=None, new_copy=False, product=None, process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABEmax using SI1310/SI1309.'
DNA.queried_features_dict['pCMV-ABEmax_16'] = DNA.dna_dict['pCMV-ABEmax'].searchdna(query='TGAGTTTAAACCCGCTGA', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description2)
DNA.queried_features_dict['pCMV-ABEmax_17'] = DNA.dna_dict['pCMV-ABEmax'].searchdna(query='TTCTTCTTTGGTGACTCG', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description2)
DNA.dna_dict['pCMV-ABEmax_0'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start=DNA.queried_features_dict['pCMV-ABEmax_16'][0].end, end=DNA.queried_features_dict['pCMV-ABEmax_17'][0].start, project='pCMV-ABEmax', product=None, process_description=description2)
DNA.dna_dict['pCMV-ABEmax_1'] = modifyends(DNA.dna_dict['pCMV-ABEmax_0'], left='TGATAATGAGTTTAAACCCGCTGA', right='CGAGTCACCAAAGAAGAAGCGGAAAGTC', project='pCMV-ABEmax', product=None, process_description=description2)
editfeature(DNA.dna_dict['pCMV-ABEmax_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f8'), project=None, new_copy=False, product=None, process_description=description2)
editfeature(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description2)
editfeature(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-8'), project=None, new_copy=False, product=None, process_description=description2)

description3 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling the two insert fragments and the backbone fragment.'
DNA.dna_dict['pCMV-BE4max_2'] = modifyends(DNA.dna_dict['pCMV-BE4max_1'], left='*{22}/-{22}', right='-{22}/*{22}', project='pCMV-BE4max', product=None, process_description=description3)
DNA.dna_dict['opt-pmCDA1-ugi_2'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_1'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', product=None, process_description=description3)
DNA.dna_dict['pCMV-ABEmax_2'] = modifyends(DNA.dna_dict['pCMV-ABEmax_1'], left='*{24}/-{24}', right='-{22}/*{22}', project='pCMV-ABEmax', product=None, process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax'] = joindna(*[DNA.dna_dict['pCMV-BE4max_2'], DNA.dna_dict['opt-pmCDA1-ugi_2'], DNA.dna_dict['pCMV-ABEmax_2']], topology='circular', project='pCMV-Target-AIDmax', product=None, process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'], = cutdna(DNA.dna_dict['pCMV-Target-AIDmax'], '22/22', crop=False, project='pCMV-Target-AIDmax', product=None, process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_1'] = joindna(*[DNA.dna_dict['pCMV-Target-AIDmax_0']], topology='circular', project='pCMV-Target-AIDmax', product=None, process_description=description3)

description4 = 'An nCas9 fragment was obtained from pCMV-Target-AIDmax using SI447/SI1105.'
DNA.queried_features_dict['pCMV-Target-AIDmax_1_19'] = DNA.dna_dict['pCMV-Target-AIDmax_1'].searchdna(query='AGCAGAACTTTAAAAGTG', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description4)
DNA.queried_features_dict['pCMV-Target-AIDmax_1_20'] = DNA.dna_dict['pCMV-Target-AIDmax_1'].searchdna(query='GTCATCGTCATCCTTGTA', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description4)
DNA.dna_dict['pCMV-Target-AIDmax_2'] = cropdna(DNA.dna_dict['pCMV-Target-AIDmax_1'], start=DNA.queried_features_dict['pCMV-Target-AIDmax_1_19'][0].end, end=DNA.queried_features_dict['pCMV-Target-AIDmax_1_20'][0].start, project='pCMV-Target-AIDmax', product=None, process_description=description4)
DNA.dna_dict['pCMV-Target-AIDmax_3'] = modifyends(DNA.dna_dict['pCMV-Target-AIDmax_2'], left='GCCACATAGCAGAACTTTAAAAGTG', right='TACAAGGATGACGATGACAAG', project='pCMV-Target-AIDmax', product=None, process_description=description4)
editfeature(DNA.dna_dict['pCMV-Target-AIDmax_3'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f9'), project=None, new_copy=False, product=None, process_description=description4)
editfeature(DNA.dna_dict['pCMV-Target-AIDmax_3'], query='f9', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description4)
editfeature(DNA.dna_dict['pCMV-Target-AIDmax_3'], query='f9', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-9'), project=None, new_copy=False, product=None, process_description=description4)

description5 = 'The rAPOBEC1 fragment was obtained from BE4max using SI1352/SI1357.'
DNA.queried_features_dict['pCMV-BE4max_21'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='GAGACTGGGCCTGTCGCC', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description5)
DNA.queried_features_dict['pCMV-BE4max_22'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='GCCTGTGGCCCACAGGAT', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description5)
DNA.dna_dict['pCMV-BE4max_3'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_features_dict['pCMV-BE4max_21'][0].end, end=DNA.queried_features_dict['pCMV-BE4max_22'][0].start, project='pCMV-BE4max', product=None, process_description=description5)
DNA.dna_dict['pCMV-BE4max_4'] = modifyends(DNA.dna_dict['pCMV-BE4max_3'], left='GATGACGATGACAAGTCTGGCTCCTCAGAGACTGGGCCTGTCGCC', right='ATCCTGTGGGCCACAGGCCTGAAG', project='pCMV-BE4max', product=None, process_description=description5)
editfeature(DNA.dna_dict['pCMV-BE4max_4'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f10'), project=None, new_copy=False, product=None, process_description=description5)
editfeature(DNA.dna_dict['pCMV-BE4max_4'], query='f10', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description5)
editfeature(DNA.dna_dict['pCMV-BE4max_4'], query='f10', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-10'), project=None, new_copy=False, product=None, process_description=description5)

description6 = 'The 2×UGI fragment was obtained from BE4max using SI1359/SI1350.'
DNA.queried_features_dict['pCMV-BE4max_23'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='AATCTGAGCGACATCATT', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description6)
DNA.queried_features_dict['pCMV-BE4max_24'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='ACTTTCCTCTTCTTCTTG', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description6)
DNA.dna_dict['pCMV-BE4max_5'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_features_dict['pCMV-BE4max_23'][0].end, end=DNA.queried_features_dict['pCMV-BE4max_24'][0].start, project='pCMV-BE4max', product=None, process_description=description6)
DNA.dna_dict['pCMV-BE4max_6'] = modifyends(DNA.dna_dict['pCMV-BE4max_5'], left='ATCCTGTGGGCCACAGGCCTGAAGACTAATCTGAGCGACATCATT', right='CAAGAAGAAGAGGAAAGTCTAATAATGAGTTTAAACCCGCTGATC', project='pCMV-BE4max', product=None, process_description=description6)
editfeature(DNA.dna_dict['pCMV-BE4max_6'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f11'), project=None, new_copy=False, product=None, process_description=description6)
editfeature(DNA.dna_dict['pCMV-BE4max_6'], query='f11', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description6)
editfeature(DNA.dna_dict['pCMV-BE4max_6'], query='f11', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-11'), project=None, new_copy=False, product=None, process_description=description6)

description7 = 'The backbone fragment was obtained from BE4max using SI1351/SI448.'
DNA.queried_features_dict['pCMV-BE4max_25'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='AGTTTAAACCCGCTGATC', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description7)
DNA.queried_features_dict['pCMV-BE4max_26'] = DNA.dna_dict['pCMV-BE4max'].searchdna(query='AAAGTTCTGCTATGTGGC', key_attribute='sequence', min_match=None, max_mismatch=0, product=None, process_description=description7)
DNA.dna_dict['pCMV-BE4max_7'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start=DNA.queried_features_dict['pCMV-BE4max_25'][0].end, end=DNA.queried_features_dict['pCMV-BE4max_26'][0].start, project='pCMV-BE4max', product=None, process_description=description7)
DNA.dna_dict['pCMV-BE4max_8'] = modifyends(DNA.dna_dict['pCMV-BE4max_7'], left='TAATAATGAGTTTAAACCCGCTGATC', right='GCCACATAGCAGAACTTTAAAAGTG', project='pCMV-BE4max', product=None, process_description=description7)
editfeature(DNA.dna_dict['pCMV-BE4max_8'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f12'), project=None, new_copy=False, product=None, process_description=description7)
editfeature(DNA.dna_dict['pCMV-BE4max_8'], query='f12', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, product=None, process_description=description7)
editfeature(DNA.dna_dict['pCMV-BE4max_8'], query='f12', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-12'), project=None, new_copy=False, product=None, process_description=description7)
DNA.dna_dict['pCMV-Target-AIDmax_4'] = modifyends(DNA.dna_dict['pCMV-Target-AIDmax_3'], left='*{25}/-{25}', right='-{15}/*{15}', project='pCMV-Target-AIDmax', product=None, process_description=description7)
DNA.dna_dict['pCMV-BE4max_9'] = modifyends(DNA.dna_dict['pCMV-BE4max_4'], left='*{15}/-{15}', right='-{24}/*{24}', project='pCMV-BE4max', product=None, process_description=description7)
DNA.dna_dict['pCMV-BE4max_10'] = modifyends(DNA.dna_dict['pCMV-BE4max_6'], left='*{24}/-{24}', right='-{26}/*{26}', project='pCMV-BE4max', product=None, process_description=description7)
DNA.dna_dict['pCMV-BE4max_11'] = modifyends(DNA.dna_dict['pCMV-BE4max_8'], left='*{26}/-{26}', right='-{25}/*{25}', project='pCMV-BE4max', product=None, process_description=description7)

description8 = 'The BE4max(C) plasmid (pCMV-BE4max(C)) was constructedby assembling the three insert fragments and the backbone fragment.'
DNA.dna_dict['pCMV-BE4max(C)'] = joindna(*[DNA.dna_dict['pCMV-Target-AIDmax_4'], DNA.dna_dict['pCMV-BE4max_9'], DNA.dna_dict['pCMV-BE4max_10'], DNA.dna_dict['pCMV-BE4max_11']], topology='circular', project='pCMV-BE4max(C)', product=None, process_description=description8)
DNA.dna_dict['pCMV-BE4max(C)_0'], = cutdna(DNA.dna_dict['pCMV-BE4max(C)'], '1316/1316', crop=False, project='pCMV-BE4max(C)', product=None, process_description=description8)
DNA.dna_dict['pCMV-BE4max(C)_1'] = joindna(*[DNA.dna_dict['pCMV-BE4max(C)_0']], topology='circular', project='pCMV-BE4max(C)', product=None, process_description=description8)
DNA.dna_dict['pCMV-BE4max(C)_1'].writedna('reconstructed_pCMV-BE4max(C).gbk')
