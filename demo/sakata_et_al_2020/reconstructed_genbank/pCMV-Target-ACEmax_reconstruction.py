import sys
sys.path.append("../../")
from dna import *
DNA.dna_dict['pCMV-ABEmax'] = DNA(seq=None, record='input/addgene_112095.gbk', project='pCMV-ABEmax', topology='linear', format=None, process_description=None)
DNA.dna_dict['opt-pmCDA1-ugi'] = DNA(seq=None, record='input/puc-optimized-pmcda1-ugi.gbk', project='opt-pmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax was amplified with primer pair SI1304/SI1307.'
DNA.dna_dict['opt-pmCDA1-ugi_0'] = cropdna(DNA.dna_dict['opt-pmCDA1-ugi'], start='525/525', end='1666/1666', project='opt-pmCDA1-ugi', process_description=description0)
DNA.dna_dict['opt-pmCDA1-ugi_1'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_0'], left='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', right='GAGAACAAAATCAAAATGCTGTGATAATGAGTTTAAACCCGCTGA', project='opt-pmCDA1-ugi', process_description=description0)
editdna(DNA.dna_dict['opt-pmCDA1-ugi_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f7'), project=None, new_copy=False, process_description=description0)
editdna(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description0)
editdna(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-7'), project=None, new_copy=False, process_description=description0)

description1 = 'The backbone fragment was amplified from pCMV-ABEmax using SI1310/SI1309.'
DNA.dna_dict['pCMV-ABEmax_0'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start='5861/5861', end='437/437', project='pCMV-ABEmax', process_description=description1)
DNA.dna_dict['pCMV-ABEmax_1'] = modifyends(DNA.dna_dict['pCMV-ABEmax_0'], left='TGATAATGAGTTTAAACCCGCTGA', right='CGAGTCACCAAAGAAGAAGCGGAAAGTC', project='pCMV-ABEmax', process_description=description1)
editdna(DNA.dna_dict['pCMV-ABEmax_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f8'), project=None, new_copy=False, process_description=description1)
editdna(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description1)
editdna(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-8'), project=None, new_copy=False, process_description=description1)

description2 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling the two insert fragments and the backbone fragment.'
DNA.dna_dict['opt-pmCDA1-ugi_2'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_1'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description2)
DNA.dna_dict['pCMV-ABEmax_2'] = modifyends(DNA.dna_dict['pCMV-ABEmax_1'], left='*{24}/-{24}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description2)

description3 = 'The BE4max(C) plasmid (pCMV-BE4max(C)) was constructedby assembling the three insert fragments and the backbone fragment.'
DNA.dna_dict['pCMV-ABEmax_3'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start='395/395', end='5729/5729', project='pCMV-ABEmax', process_description=description3)

description4 = 'An ABEmax fragment obtained from pCMV-ABEmax using SI945/SI1305'
DNA.dna_dict['pCMV-ABEmax_4'] = modifyends(DNA.dna_dict['pCMV-ABEmax_3'], left='AGATCCGCGGCCGCTAATACGACTCACTATAGG', right='CGACCTGTCTCAGCTGGGAGGTGAC', project='pCMV-ABEmax', process_description=description4)
editdna(DNA.dna_dict['pCMV-ABEmax_4'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f13'), project=None, new_copy=False, process_description=description4)
editdna(DNA.dna_dict['pCMV-ABEmax_4'], query='f13', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description4)
editdna(DNA.dna_dict['pCMV-ABEmax_4'], query='f13', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-13'), project=None, new_copy=False, process_description=description4)

description5 = 'The Target-ACEmax plasmid (pCMV-BE4max(C)) was constructed by assembling a insert fragment and two backbone fragments.'
DNA.dna_dict['opt-pmCDA1-ugi_3'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_2'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description5)
DNA.dna_dict['pCMV-ABEmax_5'] = modifyends(DNA.dna_dict['pCMV-ABEmax_2'], left='*{24}/-{24}', right='-{103}/*{103}', project='pCMV-ABEmax', process_description=description5)
DNA.dna_dict['pCMV-ABEmax_6'] = modifyends(DNA.dna_dict['pCMV-ABEmax_4'], left='*{103}/-{103}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description5)
DNA.dna_dict['pCMV-Target-ACEmax'] = joindna(*[DNA.dna_dict['opt-pmCDA1-ugi_3'], DNA.dna_dict['pCMV-ABEmax_5'], DNA.dna_dict['pCMV-ABEmax_6']], topology='circular', project='pCMV-Target-ACEmax', process_description=description5)
DNA.dna_dict['pCMV-Target-ACEmax_0'], = cutdna(DNA.dna_dict['pCMV-Target-ACEmax'], '5834/5834', crop=False, project='pCMV-Target-ACEmax', process_description=description5)
DNA.dna_dict['pCMV-Target-ACEmax_0'] = joindna(*[DNA.dna_dict['pCMV-Target-ACEmax_0']], topology='circular', project='pCMV-Target-ACEmax', process_description=description5)
DNA.dna_dict['pCMV-Target-ACEmax_0'].writedna('reconstructed_pCMV-Target-ACEmax.gbk')
