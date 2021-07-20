import sys
sys.path.append("/Users/hideto/dropbox/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github")
from QUEEN.queen import *
set_namespace(globals())
QUEEN(record='input/addgene_112095.gbk', product='pCMV_ABEmax')
QUEEN(record='input/puc-optimized-pmcda1-ugi.gbk', product='opt_CDA1_ugi')

description1 = 'The fragment encoding a nCas9 fragment (fragment6) was amplified from pCMV-BE4max using SI1308/SI1305.'
QUEEN(seq='GTCACCTCCCAGCTGAGACAGGTCG', product='SI1305', process_description=description1)

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
modifyends(fragment7, left='*{22}/-{22}', right='-{24}/*{24}', product='fragment7', process_description=description4)

description5 = 'An ABEmax fragment (fragment13) obtained from pCMV-ABEmax using SI945/SI1305'
QUEEN(seq='AGATCCGCGGCCGCTAATACGACTCACTATAGG', product='SI945', process_description=description5)
pCMV_ABEmax.searchsequence(query=SI945.seq[-18:], product='FW13', process_description=description5)
pCMV_ABEmax.searchsequence(query=SI1305.seq[-18:], product='RV13', process_description=description5)
cropdna(pCMV_ABEmax, start=FW13[0].end, end=RV13[0].start, project='addgene-112095', product='fragment13', process_description=description5)
modifyends(fragment13, left=SI945.seq, right=SI1305.rcseq, product='fragment13', process_description=description5)

description6 = 'The Target-ACEmax plasmid (pCMV-Target-ACEmax) was constructed by assembling a insert fragment and two backbone fragments.'
modifyends(fragment8, left='*{24}/-{24}', right='-{103}/*{103}', product='fragment8_2', process_description=description6)
modifyends(fragment13, left='*{103}/-{103}', right='-{22}/*{22}', product='fragment13', process_description=description6)
joindna(*[fragment7, fragment8_2, fragment13], topology='circular', product='pCMV_Target_ACEmax', process_description=description6)
pCMV_Target_ACEmax.writedna()
