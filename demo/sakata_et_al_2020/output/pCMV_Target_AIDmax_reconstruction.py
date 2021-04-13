#import sys
#sys.path.append(/Volumes/Mac_mini_ext/Dropbox\ \(Yachie\ Lab\)/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github/demo/sakata_et_al_2020/../../)
from dna import *
DNA.dna_dict['pCMV-ABEmax'] = DNA(seq=None, record='input/addgene_112095.gbk', project='pCMV-ABEmax', topology='linear', format=None, process_description=description-1)
DNA.dna_dict['pCMV-BE4max'] = DNA(seq=None, record='input/addgene_112093.gbk', project='pCMV-BE4max', topology='linear', format=None, process_description=description-1)
DNA.dna_dict['opt-pmCDA1-ugi'] = DNA(seq=None, record='input/puc-optimized-pmcda1-ugi.gb', project='opt-pmCDA1-ugi', topology='linear', format=None, process_description=description-1)

description0 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax was amplified with primer pair SI1304/SI1307.'
DNA.dna_dict['opt-pmCDA1-ugi_0'] = cropdna(DNA.dna_dict['opt-pmCDA1-ugi'], start='525/525', end='1666/1666', project='opt-pmCDA1-ugi', process_description=description0)
DNA.dna_dict['opt-pmCDA1-ugi_1'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_0'], left='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', right='GAGAACAAAATCAAAATGCTGTGATAATGAGTTTAAACCCGCTGA', project='opt-pmCDA1-ugi', process_description=description0)

description1 = 'The fragment encoding a nCas9 fragment was amplified from pCMV-BE4max using SI1308/SI1305.'
DNA.dna_dict['pCMV-BE4max_0'] = cropdna(DNA.dna_dict['pCMV-BE4max'], start='1268/1268', end='5321/5321', project='pCMV-BE4max', process_description=description1)
DNA.dna_dict['pCMV-BE4max_1'] = modifyends(DNA.dna_dict['pCMV-BE4max_0'], left='ACCAAAGAAGAAGCGGAAAGTCGACAAGAAGTACAGCATCGGCCT', right='CGACCTGTCTCAGCTGGGAGGTGAC', project='pCMV-BE4max', process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABEmax using SI1310/SI1309.'
DNA.dna_dict['pCMV-ABEmax_1'] = cropdna(DNA.dna_dict['pCMV-ABEmax'], start='5861/5861', end='437/437', project='pCMV-ABEmax', process_description=description2)
DNA.dna_dict['pCMV-ABEmax_2'] = modifyends(DNA.dna_dict['pCMV-ABEmax_1'], left='TGATAATGAGTTTAAACCCGCTGA', right='CGAGTCACCAAAGAAGAAGCGGAAAGTC', project='pCMV-ABEmax', process_description=description2)

description3 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling two insert fragments and a backbone fragment.'
DNA.dna_dict['opt-pmCDA1-ugi_2'] = modifyends(DNA.dna_dict['opt-pmCDA1-ugi_1'], left='*{22}/-{22}', right='-{24}/*{24}', project='opt-pmCDA1-ugi', process_description=description3)
DNA.dna_dict['pCMV-BE4max_2'] = modifyends(DNA.dna_dict['pCMV-BE4max_1'], left='*{22}/-{22}', right='-{22}/*{22}', project='pCMV-BE4max', process_description=description3)
DNA.dna_dict['pCMV-ABEmax_3'] = modifyends(DNA.dna_dict['pCMV-ABEmax_2'], left='*{24}/-{24}', right='-{22}/*{22}', project='pCMV-ABEmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax'] = joindna(*[DNA.dna_dict['opt-pmCDA1-ugi_2'], DNA.dna_dict['pCMV-ABEmax_3'], DNA.dna_dict['pCMV-BE4max_2']], topology='circular', project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'] = cutdna(DNA.dna_dict['pCMV-Target-AIDmax'], '4669/4669', crop='False', project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'] = joindna(*[DNA.dna_dict['pCMV-Target-AIDmax_0']], topology='linear', project='pCMV-Target-AIDmax', process_description=description3)
DNA.dna_dict['pCMV-Target-AIDmax_0'].writedna('reconstructed_pCMV-Target-AIDmax_0.gbk')
