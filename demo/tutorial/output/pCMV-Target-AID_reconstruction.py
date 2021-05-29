import sys
sys.path.append("../../")
from dnaquine import *
set_namespace(globals())
DNA(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, product='pCMV_ABE', process_description=None)
DNA(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, product='pCMV_nCas_PmCDA1_ugi', process_description=None)

description0 = 'The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'
DNA(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', record=None, project='dna', topology='linear', format=None, product='RS045', process_description=description0)
DNA(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', record=None, project='dna', topology='linear', format=None, product='HM129', process_description=description0)
pCMV_nCas_PmCDA1_ugi.searchdna(query=RS045.seq[-18:], key_attribute='sequence', min_match=None, max_mismatch=0, product='FW', process_description=description0)
pCMV_nCas_PmCDA1_ugi.searchdna(query=HM129.seq[-18:], key_attribute='sequence', min_match=None, max_mismatch=0, product='RV', process_description=description0)
cropdna(pCMV_nCas_PmCDA1_ugi, start=FW[0].end, end=RV[0].start, project='pCMV-nCas-PmCDA1-ugi', product='frag1', process_description=description0)
flipdna(HM129, project='dna', product='HM129comp', process_description=description0)
modifyends(frag1, left=RS045.seq, right=HM129comp.seq, project='pCMV-nCas-PmCDA1-ugi', product='frag1', process_description=description0)

description1 = 'The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046.'
DNA(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', record=None, project='dna', topology='linear', format=None, product='HM128', process_description=description1)
DNA(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', record=None, project='dna', topology='linear', format=None, product='RS046', process_description=description1)
pCMV_nCas_PmCDA1_ugi.searchdna(query=HM128.seq[-18:], key_attribute='sequence', min_match=None, max_mismatch=0, product='FW', process_description=description1)
pCMV_nCas_PmCDA1_ugi.searchdna(query=RS046.seq[-18:], key_attribute='sequence', min_match=None, max_mismatch=0, product='RV', process_description=description1)
cropdna(pCMV_nCas_PmCDA1_ugi, start=FW[0].end, end=RV[0].start, project='pCMV-nCas-PmCDA1-ugi', product='frag2', process_description=description1)
flipdna(RS046, project='dna', product='RS046comp', process_description=description1)
modifyends(frag2, left=HM128.seq, right=RS046comp.seq, project='pCMV-nCas-PmCDA1-ugi', product='frag2', process_description=description1)

description2 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048.'
DNA(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', record=None, project='dna', topology='linear', format=None, product='RS047', process_description=description2)
DNA(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', record=None, project='dna', topology='linear', format=None, product='RS048', process_description=description2)
pCMV_ABE.searchdna(query=RS047.seq[-18:], key_attribute='sequence', min_match=None, max_mismatch=0, product='FW', process_description=description2)
pCMV_ABE.searchdna(query=RS048.seq[-18:], key_attribute='sequence', min_match=None, max_mismatch=0, product='RV', process_description=description2)
cropdna(pCMV_ABE, start=FW[0].end, end=RV[0].start, project='pCMV_ABE', product='frag3', process_description=description2)
flipdna(RS048, project='dna', product='RS048comp', process_description=description2)
modifyends(frag3, left=RS047.seq, right=RS048comp.seq, project='pCMV_ABE', product='frag3', process_description=description2)

description3 = 'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two insert fragments and a backbone fragments.'
modifyends(frag1, left='*{25}/-{25}', right='-{28}/*{28}', project='pCMV-nCas-PmCDA1-ugi', product='frag1', process_description=description3)
modifyends(frag2, left='*{28}/-{28}', right='-{25}/*{25}', project='pCMV-nCas-PmCDA1-ugi', product='frag2', process_description=description3)
modifyends(frag3, left='*{25}/-{25}', right='-{25}/*{25}', project='pCMV_ABE', product='frag3', process_description=description3)
joindna(*[frag1, frag2, frag3], topology='circular', project='pCMV-Target-AID', product='pCMV_Target_AID', process_description=description3)
pCMV_Target_AID.writedna('reconstructed_pCMV-Target-AID.gbk')
