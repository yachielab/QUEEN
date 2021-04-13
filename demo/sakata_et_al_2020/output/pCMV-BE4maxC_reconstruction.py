#import sys
#sys.path.append("../../")
from dna import *
DNA.dna_dict['pCMV-ABEmax'] = DNA(seq=None, record='input/addgene_112095.gbk', project='pCMV-ABEmax', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-BE4max'] = DNA(seq=None, record='input/addgene_112093.gbk', project='pCMV-BE4max', topology='linear', format=None, process_description=None)
DNA.dna_dict['opt-pmCDA1-ugi'] = DNA(seq=None, record='input/puc-optimized-pmcda1-ugi.gb', project='opt-pmCDA1-ugi', topology='linear', format=None, process_description=None)

description0 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax was amplified with primer pair SI1304/SI1307.'
DNA.dna_dict['opt-pmCDA1-ugi_0'] = cropdna(DNA.dna_dict['opt-pmCDA1-ugi'], start='525/525', end='1666/1666', project='opt-pmCDA1-ugi', process_description=description0)
DNA.dna_dict['opt-pmCDA1-ugi_1'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_0'], left='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', right='GAGAACAAAATCAAAATGCTGTGATAATGAGTTTAAACCCGCTGA', project='opt-pmCDA1-ugi', process_description=description0)

description1 = 'The fragment encoding a nCas9 fragment was amplified from pCMV-BE4max using SI1308/SI1305.'
DNA.dna_dict['pCMV-BE4max_0'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start='1268/1268', end='5321/5321', project='pCMV-BE4max', process_description=description1)
DNA.dna_dict['pCMV-BE4max_1'] = modifyends(DNA.dna_dict['pCMV-BE4max_0'], left='ACCAAAGAAGAAGCGGAAAGTCGACAAGAAGTACAGCATCGGCCT', right='CGACCTGTCTCAGCTGGGAGGTGAC', project='pCMV-BE4max', process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABEmax using SI1310/SI1309.'
DNA.dna_dict['pCMV-ABEmax_0'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start='5861/5861', end='437/437', project='pCMV-ABEmax', process_description=description2)
DNA.dna_dict['pCMV-ABEmax_1'] = modifyends(DNA.dna_dict['pCMV-ABEmax_0'], left='TGATAATGAGTTTAAACCCGCTGA', right='CGAGTCACCAAAGAAGAAGCGGAAAGTC', project='pCMV-ABEmax', process_description=description2)

description3 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling the two insert fragments and the backbone fragment.'
DNA.dna_dict['opt-pmCDA1-ugi_2'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_1'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description3)
DNA.dna_dict['pCMV-BE4max_2'] = modifyends(DNA.dna_dict['pCMV-BE4max_1'], left='*{22}/-{22}', right='-{22}/*{22}', project='pCMV-BE4max', process_description=description3)
DNA.dna_dict['pCMV-ABEmax_2'] = modifyends(DNA.dna_dict['pCMV-ABEmax_1'], left='*{24}/-{24}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax'] = joindna(*[DNA.dna_dict['opt-pmCDA1-ugi_2'], DNA.dna_dict['pCMV-ABEmax_2'], DNA.dna_dict['pCMV-BE4max_2']], topology='circular', project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'], = cutdna(DNA.dna_dict['pCMV-Target-AIDmax'], '4669/4669', crop=False, project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'] = joindna(*[DNA.dna_dict['pCMV-Target-AIDmax_0']], topology='circular', project='pCMV-Target-AIDmax', process_description=description3)

description4 = 'An nCas9 fragment was obtained from pCMV-Target-AIDmax using SI447/SI1105.'
DNA.dna_dict['pCMV-Target-AIDmax_1'] = cropdna(DNA.dna_dict['pCMV-Target-AIDmax_0'], start='7411/7411', end='4330/4330', project='pCMV-Target-AIDmax', process_description=description4)
DNA.dna_dict['pCMV-Target-AIDmax_2'] = modifyends(DNA.dna_dict['pCMV-Target-AIDmax_1'], left='GCCACATAGCAGAACTTTAAAAGTG', right='TACAAGGATGACGATGACAAG', project='pCMV-Target-AIDmax', process_description=description4)

description5 = 'The rAPOBEC1 fragment was obtained from BE4max using SI1352/SI1357.'
DNA.dna_dict['pCMV-BE4max_3'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start='489/489', end='1125/1125', project='pCMV-BE4max', process_description=description5)
DNA.dna_dict['pCMV-BE4max_4'] = modifyends(DNA.dna_dict['pCMV-BE4max_3'], left='GATGACGATGACAAGTCTGGCTCCTCAGAGACTGGGCCTGTCGCC', right='ATCCTGTGGGCCACAGGCCTGAAG', project='pCMV-BE4max', process_description=description5)

description6 = 'The 2×UGI fragment was obtained from BE4max using SI1359/SI1350.'
DNA.dna_dict['pCMV-BE4max_5'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start='5397/5397', end='5948/5948', project='pCMV-BE4max', process_description=description6)
DNA.dna_dict['pCMV-BE4max_6'] = modifyends(DNA.dna_dict['pCMV-BE4max_5'], left='ATCCTGTGGGCCACAGGCCTGAAGACTAATCTGAGCGACATCATT', right='CAAGAAGAAGAGGAAAGTCTAATAATGAGTTTAAACCCGCTGATC', project='pCMV-BE4max', process_description=description6)

description7 = 'The backbone fragment was obtained from BE4max using SI1351/SI448.'
DNA.dna_dict['pCMV-BE4max_7'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start='6013/6013', end='8110/8110', project='pCMV-BE4max', process_description=description7)
DNA.dna_dict['pCMV-BE4max_8'] = modifyends(DNA.dna_dict['pCMV-BE4max_7'], left='TAATAATGAGTTTAAACCCGCTGATC', right='GCCACATAGCAGAACTTTAAAAGTG', project='pCMV-BE4max', process_description=description7)
DNA.dna_dict['pCMV-Target-AIDmax_3'] = modifyends(DNA.dna_dict['pCMV-Target-AIDmax_2'], left='*{25}/-{25}', right='-{15}/*{15}', project='pCMV-Target-AIDmax', process_description=description7)
DNA.dna_dict['pCMV-BE4max_9'] = modifyends(DNA.dna_dict['pCMV-BE4max_4'], left='*{15}/-{15}', right='-{24}/*{24}', project='pCMV-BE4max', process_description=description7)
DNA.dna_dict['pCMV-BE4max_10'] = modifyends(DNA.dna_dict['pCMV-BE4max_6'], left='*{24}/-{24}', right='-{26}/*{26}', project='pCMV-BE4max', process_description=description7)
DNA.dna_dict['pCMV-BE4max_11'] = modifyends(DNA.dna_dict['pCMV-BE4max_8'], left='*{26}/-{26}', right='-{25}/*{25}', project='pCMV-BE4max', process_description=description7)

description8 = 'The BE4max(C) plasmid (pCMV-BE4max(C)) was constructedby assembling the three insert fragments and the backbone fragment.'
DNA.dna_dict['pCMV-BE4max(C)'] = joindna(*[DNA.dna_dict['pCMV-Target-AIDmax_3'], DNA.dna_dict['pCMV-BE4max_9'], DNA.dna_dict['pCMV-BE4max_10'], DNA.dna_dict['pCMV-BE4max_11']], topology='circular', project='pCMV-BE4max(C)', process_description=description8)
DNA.dna_dict['pCMV-BE4max(C)_0'], = cutdna(DNA.dna_dict['pCMV-BE4max(C)'], '1339/1339', crop=False, project='pCMV-BE4max(C)', process_description=description8)
DNA.dna_dict['pCMV-BE4max(C)_0'] = joindna(*[DNA.dna_dict['pCMV-BE4max(C)_0']], topology='circular', project='pCMV-BE4max(C)', process_description=description8)
DNA.dna_dict['pCMV-BE4max(C)_0'].writedna('reconstructed_pCMV-BE4max(C).gbk')
