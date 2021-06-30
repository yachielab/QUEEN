import sys
sys.path.append("/Volumes/Mac_mini_ext/Dropbox (Yachie Lab)/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github/demo/tutorial/../..")
from QUEEN.queen import *
set_namespace(globals())
pCMV_ABE = QUEEN(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, product='pCMV_ABE', process_description=None)
pCMV_nCas_PmCDA1_ugi = QUEEN(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, product='pCMV_nCas_PmCDA1_ugi', process_description=None)

description0 = 'The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'
RS045 = QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', record=None, project='dna', topology='linear', format=None, product='RS045', process_description=description0)
HM129 = QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', record=None, project='dna', topology='linear', format=None, product='HM129', process_description=description0)
FW1 = pCMV_nCas_PmCDA1_ugi.searchsequence(query=RS045.seq[-18:], product='FW1', process_description=description0)
RV1 = pCMV_nCas_PmCDA1_ugi.searchsequence(query=HM129.seq[-18:], product='RV1', process_description=description0)
fragment1 = cropdna(pCMV_nCas_PmCDA1_ugi, start=FW1[0].end, end=RV1[0].start, project='pCMV-nCas-PmCDA1-ugi', product='fragment1', process_description=description0)
fragment1 = modifyends(fragment1, left=RS045.seq, right=HM129.getdnaseq(strand=-1), project='pCMV-nCas-PmCDA1-ugi', product='fragment1', process_description=description0)

description1 = 'The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046.'
HM128 = QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', record=None, project='dna', topology='linear', format=None, product='HM128', process_description=description1)
RS046 = QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', record=None, project='dna', topology='linear', format=None, product='RS046', process_description=description1)
FW2 = pCMV_nCas_PmCDA1_ugi.searchsequence(query=HM128.seq[-18:], product='FW2', process_description=description1)
RV2 = pCMV_nCas_PmCDA1_ugi.searchsequence(query=RS046.seq[-18:], product='RV2', process_description=description1)
fragment2 = cropdna(pCMV_nCas_PmCDA1_ugi, start=FW2[0].end, end=RV2[0].start, project='pCMV-nCas-PmCDA1-ugi', product='fragment2', process_description=description1)
fragment2 = modifyends(fragment2, left=HM128.seq, right=RS046.getdnaseq(strand=-1), project='pCMV-nCas-PmCDA1-ugi', product='fragment2', process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048.'
RS047 = QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', record=None, project='dna', topology='linear', format=None, product='RS047', process_description=description2)
RS048 = QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', record=None, project='dna', topology='linear', format=None, product='RS048', process_description=description2)
FW3 = pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_description=description2)
RV3 = pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_description=description2)
fragment3 = cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, project='pCMV_ABE', product='fragment3', process_description=description2)
fragment3 = modifyends(fragment3, left=RS047.seq, right=RS048.getdnaseq(strand=-1), project='pCMV_ABE', product='fragment3', process_description=description2)

description3 = 'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two insert fragments and a backbone fragments.'
fragment1 = modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', project='pCMV-nCas-PmCDA1-ugi', product='fragment1', process_description=description3)
fragment2 = modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', project='pCMV-nCas-PmCDA1-ugi', product='fragment2', process_description=description3)
fragment3 = modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV_ABE', product='fragment3', process_description=description3)
pCMV_Target_AID = joindna(*[fragment1, fragment2, fragment3], topology='circular', project='pCMV-Target-AID', product='pCMV_Target_AID', process_description=description3)
pCMV_Target_AID.writedna()
