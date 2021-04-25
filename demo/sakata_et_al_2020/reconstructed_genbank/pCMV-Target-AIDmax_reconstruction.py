import sys
sys.path.append("../../")
from dna import *
DNA.dna_dict['pCMV-ABEmax'] = DNA(seq=None, record='input/addgene_112095.gbk', project='pCMV-ABEmax', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-BE4max'] = DNA(seq=None, record='input/addgene_112093.gbk', project='pCMV-BE4max', topology='linear', format=None, process_description=None)
DNA.dna_dict['opt-pmCDA1-ugi'] = DNA(seq=None, record='input/puc-optimized-pmcda1-ugi.gbk', project='opt-pmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The fragment encoding a nCas9 fragment was amplified from pCMV-BE4max using SI1308/SI1305.'
DNA.dna_dict['pCMV-BE4max_0'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start='1268/1268', end='5321/5321', project='pCMV-BE4max', process_description=description0)
DNA.dna_dict['pCMV-BE4max_1'] = modifyends(DNA.dna_dict['pCMV-BE4max_0'], left='ACCAAAGAAGAAGCGGAAAGTCGACAAGAAGTACAGCATCGGCCT', right='CGACCTGTCTCAGCTGGGAGGTGAC', project='pCMV-BE4max', process_description=description0)
editdna(DNA.dna_dict['pCMV-BE4max_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f6'), project=None, new_copy=False, process_description=description0)
editdna(DNA.dna_dict['pCMV-BE4max_1'], query='f6', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description0)
editdna(DNA.dna_dict['pCMV-BE4max_1'], query='f6', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-6'), project=None, new_copy=False, process_description=description0)

description1 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax was amplified with primer pair SI1304/SI1307.'
DNA.dna_dict['opt-pmCDA1-ugi_0'] = cropdna(DNA.dna_dict['opt-pmCDA1-ugi'], start='525/525', end='1666/1666', project='opt-pmCDA1-ugi', process_description=description1)
DNA.dna_dict['opt-pmCDA1-ugi_1'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_0'], left='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', right='GAGAACAAAATCAAAATGCTGTGATAATGAGTTTAAACCCGCTGA', project='opt-pmCDA1-ugi', process_description=description1)
editdna(DNA.dna_dict['opt-pmCDA1-ugi_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f7'), project=None, new_copy=False, process_description=description1)
editdna(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description1)
editdna(DNA.dna_dict['opt-pmCDA1-ugi_1'], query='f7', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-7'), project=None, new_copy=False, process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABEmax using SI1310/SI1309.'
DNA.dna_dict['pCMV-ABEmax_0'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start='5861/5861', end='437/437', project='pCMV-ABEmax', process_description=description2)
DNA.dna_dict['pCMV-ABEmax_1'] = modifyends(DNA.dna_dict['pCMV-ABEmax_0'], left='TGATAATGAGTTTAAACCCGCTGA', right='CGAGTCACCAAAGAAGAAGCGGAAAGTC', project='pCMV-ABEmax', process_description=description2)
editdna(DNA.dna_dict['pCMV-ABEmax_1'], query=None, key_attribute='sequence', min_match=None, max_mismatch=0, target_attribute='feature_id', operation=createattribute(value='f8'), project=None, new_copy=False, process_description=description2)
editdna(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='strand', operation=replaceattribute(query_re=None,value=0), project=None, new_copy=False, process_description=description2)
editdna(DNA.dna_dict['pCMV-ABEmax_1'], query='f8', key_attribute='feature_id', min_match=None, max_mismatch=0, target_attribute='qualifier:label', operation=createattribute(value='fragment-8'), project=None, new_copy=False, process_description=description2)

description3 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling the two insert fragments and the backbone fragment.'
DNA.dna_dict['pCMV-BE4max_2'] = modifyends(DNA.dna_dict['pCMV-BE4max_1'], left='*{22}/-{22}', right='-{22}/*{22}', project='pCMV-BE4max', process_description=description3)
DNA.dna_dict['opt-pmCDA1-ugi_2'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_1'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description3)
DNA.dna_dict['pCMV-ABEmax_2'] = modifyends(DNA.dna_dict['pCMV-ABEmax_1'], left='*{24}/-{24}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax'] = joindna(*[DNA.dna_dict['pCMV-BE4max_2'], DNA.dna_dict['opt-pmCDA1-ugi_2'], DNA.dna_dict['pCMV-ABEmax_2']], topology='circular', project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'], = cutdna(DNA.dna_dict['pCMV-Target-AIDmax'], '45/45', crop=False, project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'] = joindna(*[DNA.dna_dict['pCMV-Target-AIDmax_0']], topology='circular', project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'].writedna('reconstructed_pCMV-Target-AIDmax.gbk')
