import sys
sys.path.append("/Users/hideto/dropbox/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github")
from QUEEN.queen import *
set_namespace(globals())
QUEEN(record='input/addgene_112093.gbk', product='pCMV_BE4max')
QUEEN(record='input/addgene_112095.gbk', product='pCMV_ABEmax')
QUEEN(record='input/puc-optimized-pmcda1-ugi.gbk', product='opt_CDA1_ugi')

description1 = 'The fragment encoding a nCas9 fragment (fragment6) was amplified from pCMV-BE4max using SI1308/SI1305.'
QUEEN(seq='ACCAAAGAAGAAGCGGAAAGTCGACAAGAAGTACAGCATCGGCCT', product='SI1308', process_description=description1)
QUEEN(seq='GTCACCTCCCAGCTGAGACAGGTCG', product='SI1305', process_description=description1)
pCMV_BE4max.searchsequence(query=SI1308.seq[-18:], product='FW6', process_description=description1)
pCMV_BE4max.searchsequence(query=SI1305.seq[-18:], product='RV6', process_description=description1)
cropdna(pCMV_BE4max, start=FW6[0].end, end=RV6[0].start, project='addgene-112093', product='fragment6', process_description=description1)
modifyends(fragment6, left=SI1308.seq, right=SI1305.rcseq, product='fragment6', process_description=description1)

description2 = 'The fragment encoding the codon-optimized C-terminal region of Target-AIDmax (fragment7) was amplified with primer pair SI1304/SI1307.'
QUEEN(seq='CCTGTCTCAGCTGGGAGGTGACGGCGGAGGAGGAACTGGAGGAGG', product='SI1304', process_description=description2)
QUEEN(seq='TCAGCGGGTTTAAACTCATTATCACAGCATTTTGATTTTGTTCTC', product='SI1307', process_description=description2)
opt_CDA1_ugi.searchsequence(query=SI1304.seq[-18:], product='FW7', process_description=description2)
opt_CDA1_ugi.searchsequence(query=SI1307.seq[-18:], product='RV7', process_description=description2)
cropdna(opt_CDA1_ugi, start=FW7[0].end, end=RV7[0].start, project='puc-optimized-pmcda1-ugi', product='fragment7', process_description=description2)
modifyends(fragment7, left=SI1304.seq, right=SI1307.rcseq, product='fragment7', process_description=description2)

description3 = 'The backbone fragment (fragment8) was amplified from pCMV-ABEmax using SI1310/SI1309.'
QUEEN(seq='TGATAATGAGTTTAAACCCGCTGA', product='SI1310', process_description=description3)
QUEEN(seq='GACTTTCCGCTTCTTCTTTGGTGACTCG', product='SI1309', process_description=description3)
pCMV_ABEmax.searchsequence(query=SI1310.seq[-18:], product='FW8', process_description=description3)
pCMV_ABEmax.searchsequence(query=SI1309.seq[-18:], product='RV8', process_description=description3)
cropdna(pCMV_ABEmax, start=FW8[0].end, end=RV8[0].start, project='addgene-112095', product='fragment8', process_description=description3)
modifyends(fragment8, left=SI1310.seq, right=SI1309.rcseq, product='fragment8', process_description=description3)

description4 = 'The Target-AIDmax plasmid (pCMV-Target-AIDmax) was constructed by assembling the two insert fragments and the backbone fragment.'
modifyends(fragment6, left='*{22}/-{22}', right='-{22}/*{22}', product='fragment6', process_description=description4)
modifyends(fragment7, left='*{22}/-{22}', right='-{24}/*{24}', product='fragment7', process_description=description4)
modifyends(fragment8, left='*{24}/-{24}', right='-{22}/*{22}', product='fragment8_1', process_description=description4)
joindna(*[fragment6, fragment7, fragment8_1], topology='circular', product='pCMV_Target_AIDmax', process_description=description4)

description5 = 'An nCas9 fragment (fragment9) was obtained from pCMV-Target-AIDmax using SI447/SI1105.'
QUEEN(seq='GCCACATAGCAGAACTTTAAAAGTG', product='SI447', process_description=description5)
QUEEN(seq='CTTGTCATCGTCATCCTTGTA', product='SI1105', process_description=description5)
pCMV_Target_AIDmax.searchsequence(query=SI447.seq[-18:], product='FW9', process_description=description5)
pCMV_Target_AIDmax.searchsequence(query=SI1105.seq[-18:], product='RV9', process_description=description5)
cropdna(pCMV_Target_AIDmax, start=FW9[0].end, end=RV9[0].start, project='addgene-112093', product='fragment9', process_description=description5)
modifyends(fragment9, left=SI447.seq, right=SI1105.rcseq, product='fragment9', process_description=description5)

description6 = 'The rAPOBEC1 fragment (fragment10) was obtained from BE4max using SI1352/SI1357.'
QUEEN(seq='GATGACGATGACAAGTCTGGCTCCTCAGAGACTGGGCCTGTCGCC', product='SI1352', process_description=description6)
QUEEN(seq='CTTCAGGCCTGTGGCCCACAGGAT', product='SI1357', process_description=description6)
pCMV_BE4max.searchsequence(query=SI1352.seq[-18:], product='FW10', process_description=description6)
pCMV_BE4max.searchsequence(query=SI1357.seq[-18:], product='RV10', process_description=description6)
cropdna(pCMV_BE4max, start=FW10[0].end, end=RV10[0].start, project='addgene-112093', product='fragment10', process_description=description6)
modifyends(fragment10, left=SI1352.seq, right=SI1357.rcseq, product='fragment10', process_description=description6)

description7 = 'The 2×UGI fragment (fragment11) was obtained from BE4max using SI1359/SI1350.'
QUEEN(seq='ATCCTGTGGGCCACAGGCCTGAAGACTAATCTGAGCGACATCATT', product='SI1359', process_description=description7)
QUEEN(seq='GATCAGCGGGTTTAAACTCATTATTAGACTTTCCTCTTCTTCTTG', product='SI1350', process_description=description7)
pCMV_BE4max.searchsequence(query=SI1359.seq[-18:], product='FW11', process_description=description7)
pCMV_BE4max.searchsequence(query=SI1350.seq[-18:], product='RV11', process_description=description7)
cropdna(pCMV_BE4max, start=FW11[0].end, end=RV11[0].start, project='addgene-112093', product='fragment11', process_description=description7)
modifyends(fragment11, left=SI1359.seq, right=SI1350.rcseq, product='fragment11', process_description=description7)

description8 = 'The backbone fragment (fragment12) was obtained from BE4max using SI1351/SI448.'
QUEEN(seq='TAATAATGAGTTTAAACCCGCTGATC', product='SI1351', process_description=description8)
QUEEN(seq='CACTTTTAAAGTTCTGCTATGTGGC', product='SI448', process_description=description8)
pCMV_BE4max.searchsequence(query=SI1351.seq[-18:], product='FW12', process_description=description8)
pCMV_BE4max.searchsequence(query=SI448.seq[-18:], product='RV12', process_description=description8)
cropdna(pCMV_BE4max, start=FW12[0].end, end=RV12[0].start, project='addgene-112093', product='fragment12', process_description=description8)
modifyends(fragment12, left=SI1351.seq, right=SI448.rcseq, product='fragment12', process_description=description8)

description9 = 'The BE4max(C) plasmid (pCMV-BE4max(C)) was constructedby assembling the three insert fragments and the backbone fragment.'
modifyends(fragment9, left='*{25}/-{25}', right='-{15}/*{15}', product='fragment9', process_description=description9)
modifyends(fragment10, left='*{15}/-{15}', right='-{24}/*{24}', product='fragment10', process_description=description9)
modifyends(fragment11, left='*{24}/-{24}', right='-{26}/*{26}', product='fragment11', process_description=description9)
modifyends(fragment12, left='*{26}/-{26}', right='-{25}/*{25}', product='fragment12', process_description=description9)
joindna(*[fragment9, fragment10, fragment11, fragment12], topology='circular', product='pCMV_BE4maxC', process_description=description9)
pCMV_BE4maxC.writedna()
