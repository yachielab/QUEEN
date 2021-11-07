project='pCMV_Target_AID'
import sys
sys.path.append("/usr/local/lib/python3.9/site-packages/QUEEN-1.0.0-py3.9.egg")
from QUEEN.queen import *
from QUEEN import cutsite as cs
set_namespace(globals())

QUEEN(record='https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf', dbtype='benchling', product='pCMV_ABE', process_id='pCMV_Target_AID-8VZ56URU5AVQ9VYZDJJUEMGH', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
QUEEN(record='https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF', dbtype='benchling', product='pcDNA31_Target_AID', process_id='pCMV_Target_AID-86JH1MQ2J28S2J1AKDNMTO9P', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', product='RS045', process_id='pCMV_Target_AID-7RM7JLKHG34A94JPI89T0L1E', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', product='HM129', process_id='pCMV_Target_AID-5EBTBWHTA48KBOPKTTB9JT6A', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')

process1={'name':'PCR', 'description':'1. The N-terminus half of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair RS045/HM129.'}
pcDNA31_Target_AID.searchsequence(query=RS045.seq[-18:], product='FW1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-5YB6M433X8LB6QPABJLD7ISH', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
pcDNA31_Target_AID.searchsequence(query=HM129.seq[-18:], product='RV1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-389AO7EEP8NR13M9YT8JSG3P', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
cropdna(pcDNA31_Target_AID, start=FW1[0].end, end=RV1[0].start, product='extract1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-EX2BR21FOQZ56MRDO04MJN5D', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
modifyends(extract1, left=RS045.seq, right=HM129.rcseq, product='fragment1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-E6R3VS0JFWX919ZQ6SFJEYCT', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')

QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', product='HM128', process_id='pCMV_Target_AID-D17S3ZCH066ODNM8284G2WP5', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', product='RS046', process_id='pCMV_Target_AID-1S3LEXX6EV9MAV2MDAB7ASCE', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')

process2={'name':'PCR', 'description':'2. The C-terminus half of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair HM128/RS046.'}
pcDNA31_Target_AID.searchsequence(query=HM128.seq[-18:], product='FW2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-AA7IVK1DD30L7W9DQYU0EN74', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
pcDNA31_Target_AID.searchsequence(query=RS046.seq[-18:], product='RV2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-1572PU73FBRVB4HX3GY338EA', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
cropdna(pcDNA31_Target_AID, start=FW2[0].end, end=RV2[0].start, product='extract2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-86ANG0RMGV7GA2VDROZRK94F', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
modifyends(extract2, left=HM128.seq, right=RS046.rcseq, product='fragment2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-D4LHRLTQYVWS7DDHHNZIEBT6', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')

QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', product='RS047', process_id='pCMV_Target_AID-5EE95K2TXRV3AWL24E8JAUSU', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', product='RS048', process_id='pCMV_Target_AID-9568WJMB57R89TX74YCOVHCJ', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')

process3={'name':'PCR', 'description':'3. A backbone fragment was amplified from pCMV-ABE7.10 using the primer pair RS047/RS048.'}
pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-7XB61E4ZQ7Y9R8SGLIXCNHF4', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-6KKHICD6LV912NUOUFKKWIKE', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, product='extract3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-5RX0MGOVFY413CM08BGT6WC7', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
modifyends(extract3, left=RS047.seq, right=RS048.rcseq, product='fragment3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-6EBK7PW1I1LB3U0IV1WAHX4Z', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')

process4={'name':'Gibson Assembly', 'description':'4. The three fragments were assembled by Gibson Assembly.'}
modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', product='fragment1_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-1GNEIO7LG6GNCGTRMES7HIIK', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', product='fragment2_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-EGCOHDIF7TTO1FOTMSGJOP0E', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment3_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-5U5C5B6FKDQ4BA20NSQMGS2I', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
joindna(*[fragment1_mod, fragment2_mod, fragment3_mod], topology='circular', product='pCMV_Target_AID', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-8IXV92W0BW934FILCZ1TOH7R', original_ids=[], _sourcefile='pCMV-Target-AID_reconstruction')
if __name__ == '__main__':
    pCMV_Target_AID.outputgbk()
