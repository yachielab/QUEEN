import sys
sys.path.append("/Users/hideto/dropbox/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github")
from QUEEN.queen import *
set_namespace(globals())
QUEEN(record='input/addgene_102919.gbk', product='pCMV_ABE')
QUEEN(record='input/addgene_79620.gbk', product='pCMV_nCas_CDA1_ugi')

description1 = 'The C-terminus of Target-AID (fragment2) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046.'
QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', product='RS046', process_description=description1)

description2 = 'The backbone fragment (fragment3) was amplified from pCMV-ABE7.10 using RS047/RS048'
QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', product='RS047', process_description=description2)
QUEEN(seq='ACCTCCTCCACCGTCACCCCCAAGCTGTGACA', product='RS052', process_description=description2)

description3 = 'The backbone fragment (fragment4) was amplified from pCMV-ABE7.10 using RS047/RS052'
pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW4', process_description=description3)
pCMV_ABE.searchsequence(query=RS052.seq[-18:], product='RV4', process_description=description3)
cropdna(pCMV_ABE, start=FW4[0].end, end=RV4[0].start, project='addgene-102919', product='fragment4', process_description=description3)
modifyends(fragment4, left=RS047.seq, right=RS052.rcseq, product='fragment4', process_description=description3)
QUEEN(seq='GCTTGGGGGTGACGGTGGAGGAGGTACCGGCGG', product='RS051', process_description=description3)

description4 = 'The insert fragment (fragment5) encoding the C-terminus region of Target-AID was amplified from pcDNA-pCMV-nCas9 using RS051/RS046.'
pCMV_nCas_CDA1_ugi.searchsequence(query=RS051.seq[-18:], product='FW5', process_description=description4)
pCMV_nCas_CDA1_ugi.searchsequence(query=RS046.seq[-18:], product='RV5', process_description=description4)
cropdna(pCMV_nCas_CDA1_ugi, start=FW5[0].end, end=RV5[0].start, project='addgene-79620', product='fragment5', process_description=description4)
modifyends(fragment5, left=RS051.seq, right=RS046.rcseq, product='fragment5', process_description=description4)

description5 = 'The Target-ACE plasmid (pCMV-Target-ACE) was constructed by assembling the insert fragment and a backbone fragment.'
modifyends(fragment4, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment4', process_description=description5)
modifyends(fragment5, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment5', process_description=description5)
joindna(*[fragment4, fragment5], topology='circular', product='pCMV_Target_ACE', process_description=description5)
pCMV_Target_ACE.writedna()
