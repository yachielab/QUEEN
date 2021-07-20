import sys
sys.path.append("/Users/hideto/dropbox/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github")
from QUEEN.queen import *
set_namespace(globals())
QUEEN(record='input/addgene_102919.gbk', product='pCMV_ABE')
QUEEN(record='input/addgene_79620.gbk', product='pCMV_nCas_CDA1_ugi')

description1 = 'The N-terminus of Target-AID (fragment1) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'
QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', product='RS045', process_description=description1)
QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', product='HM129', process_description=description1)
pCMV_nCas_CDA1_ugi.searchsequence(query=RS045.seq[-18:], product='FW1', process_description=description1)
pCMV_nCas_CDA1_ugi.searchsequence(query=HM129.seq[-18:], product='RV1', process_description=description1)
cropdna(pCMV_nCas_CDA1_ugi, start=FW1[0].end, end=RV1[0].start, project='addgene-79620', product='fragment1', process_description=description1)
modifyends(fragment1, left=RS045.seq, right=HM129.rcseq, product='fragment1', process_description=description1)

description2 = 'The C-terminus of Target-AID (fragment2) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046.'
QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', product='HM128', process_description=description2)
QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', product='RS046', process_description=description2)
pCMV_nCas_CDA1_ugi.searchsequence(query=HM128.seq[-18:], product='FW2', process_description=description2)
pCMV_nCas_CDA1_ugi.searchsequence(query=RS046.seq[-18:], product='RV2', process_description=description2)
cropdna(pCMV_nCas_CDA1_ugi, start=FW2[0].end, end=RV2[0].start, project='addgene-79620', product='fragment2', process_description=description2)
modifyends(fragment2, left=HM128.seq, right=RS046.rcseq, product='fragment2', process_description=description2)

description3 = 'The backbone fragment (fragment3) was amplified from pCMV-ABE7.10 using RS047/RS048'
QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', product='RS047', process_description=description3)
QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', product='RS048', process_description=description3)
pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_description=description3)
pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_description=description3)
cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, project='addgene-102919', product='fragment3', process_description=description3)
modifyends(fragment3, left=RS047.seq, right=RS048.rcseq, product='fragment3', process_description=description3)

description4 = 'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling the two insert fragments and the backbone fragment.'
modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', product='fragment1', process_description=description4)
modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', product='fragment2', process_description=description4)
modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment3', process_description=description4)
joindna(*[fragment1, fragment2, fragment3], topology='circular', product='pCMV_Target_AID', process_description=description4)
pCMV_Target_AID.writedna()
