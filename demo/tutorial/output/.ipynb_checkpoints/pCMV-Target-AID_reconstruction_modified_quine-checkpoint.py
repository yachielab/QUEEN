QUEEN(record='https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF', dbtype='benchling', project='pcDNA_Target_AID', product='pcDNA_Target_AID', process_id='pCMV_Target_AID-EZbi1moI0C', originals=[])
QUEEN(record='https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF', dbtype='benchling', project='pcDNA_Target_AID', product='pcDNA_Target_AID', process_id='pCMV_Target_AID_mod-EZbi1moI0C', originals=['pCMV_Target_AID-EZbi1moI0C'])

QUEEN(record='https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf', dbtype='benchling', project='pCMV_ABE', product='pCMV_ABE', process_id='pCMV_Target_AID-rA2RRtBu3d', originals=[])
QUEEN(record='https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf', dbtype='benchling', project='pCMV_ABE', product='pCMV_ABE', process_id='pCMV_Target_AID_mod-rA2RRtBu3d', originals=['pCMV_Target_AID-rA2RRtBu3d'])

QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', project='RS045', product='RS045', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-WpgGb7XWl5', originals=[])
QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', project='RS045', product='RS045', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-WpgGb7XWl5', originals=['pCMV_Target_AID-WpgGb7XWl5'])

QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', project='HM129', product='HM129', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-q3VxR7sxc2', originals=[])
QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', project='HM129', product='HM129', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-q3VxR7sxc2', originals=['pCMV_Target_AID-q3VxR7sxc2'])

pcDNA_Target_AID.searchsequence(query=RS045.seq[-18:], product='FW1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-ncQKuKboBk', originals=[])
pcDNA_Target_AID.searchsequence(query=RS045.seq[-18:], product='FW1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-ncQKuKboBk', originals=['pCMV_Target_AID-ncQKuKboBk'])

pcDNA_Target_AID.searchsequence(query=HM129.seq[-18:], product='RV1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-Yhy8TvI76F', originals=[])
pcDNA_Target_AID.searchsequence(query=HM129.seq[-18:], product='RV1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-Yhy8TvI76F', originals=['pCMV_Target_AID-Yhy8TvI76F'])

cropdna(pcDNA_Target_AID, start=FW1[0].end, end=RV1[0].start, product='extract1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-3AjRaQxViH', originals=[])
cropdna(pcDNA_Target_AID, start=FW1[0].end, end=RV1[0].start, product='extract1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-3AjRaQxViH', originals=['pCMV_Target_AID-3AjRaQxViH'])

modifyends(extract1, left=RS045.seq, right=HM129.rcseq, product='fragment1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-tm34jxHzZx', originals=[])
modifyends(extract1, left=RS045.seq, right=HM129.rcseq, product='fragment1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-tm34jxHzZx', originals=['pCMV_Target_AID-tm34jxHzZx'])

QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', project='HM128', product='HM128', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-Pdu8Zk8Kcu', originals=[])
QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', project='HM128', product='HM128', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-Pdu8Zk8Kcu', originals=['pCMV_Target_AID-Pdu8Zk8Kcu'])

QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', project='RS046', product='RS046', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-a96XNwwKYV', originals=[])
QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', project='RS046', product='RS046', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-a96XNwwKYV', originals=['pCMV_Target_AID-a96XNwwKYV'])

pcDNA_Target_AID.searchsequence(query=HM128.seq[-18:], product='FW2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-oBPkpzeOwc', originals=[])
pcDNA_Target_AID.searchsequence(query=HM128.seq[-18:], product='FW2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-oBPkpzeOwc', originals=['pCMV_Target_AID-oBPkpzeOwc'])

pcDNA_Target_AID.searchsequence(query=RS046.seq[-18:], product='RV2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-4NvpQbAVb8', originals=[])
pcDNA_Target_AID.searchsequence(query=RS046.seq[-18:], product='RV2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-4NvpQbAVb8', originals=['pCMV_Target_AID-4NvpQbAVb8'])

cropdna(pcDNA_Target_AID, start=FW2[0].end, end=RV2[0].start, product='extract2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-r6NI2fSBK4', originals=[])
cropdna(pcDNA_Target_AID, start=FW2[0].end, end=RV2[0].start, product='extract2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-r6NI2fSBK4', originals=['pCMV_Target_AID-r6NI2fSBK4'])

modifyends(extract2, left=HM128.seq, right=RS046.rcseq, product='fragment2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-bjx9zDibMA', originals=[])
modifyends(extract2, left=HM128.seq, right=RS046.rcseq, product='fragment2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-bjx9zDibMA', originals=['pCMV_Target_AID-bjx9zDibMA'])

QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', project='RS047', product='RS047', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-xbYwpI6slr', originals=[])
QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', project='RS047', product='RS047', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-xbYwpI6slr', originals=['pCMV_Target_AID-xbYwpI6slr'])

QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', project='RS048', product='RS048', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-uANRoOl8Yj', originals=[])
QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', project='RS048', product='RS048', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-uANRoOl8Yj', originals=['pCMV_Target_AID-uANRoOl8Yj'])

pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-ZRzlemHnPL', originals=[])
pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-ZRzlemHnPL', originals=['pCMV_Target_AID-ZRzlemHnPL'])

pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-OZqgFZtJ07', originals=[])
pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-OZqgFZtJ07', originals=['pCMV_Target_AID-OZqgFZtJ07'])

cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, product='extract3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-4Syq4DDNMa', originals=[])
cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, product='extract3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-4Syq4DDNMa', originals=['pCMV_Target_AID-4Syq4DDNMa'])

modifyends(extract3, left=RS047.seq, right=RS048.rcseq, product='fragment3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-nxOQby3rWS', originals=[])
modifyends(extract3, left=RS047.seq, right=RS048.rcseq, product='fragment3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-nxOQby3rWS', originals=['pCMV_Target_AID-nxOQby3rWS'])

modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', product='fragment1_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-fNk9Pxad1N', originals=[])
modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', product='fragment1_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-fNk9Pxad1N', originals=['pCMV_Target_AID-fNk9Pxad1N'])

modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', product='fragment2_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-V9witA8v4C', originals=[])
modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', product='fragment2_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-V9witA8v4C', originals=['pCMV_Target_AID-V9witA8v4C'])

modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment3_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-CGjYNsICes', originals=[])
modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment3_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-CGjYNsICes', originals=['pCMV_Target_AID-CGjYNsICes'])

joindna(*[fragment1_mod, fragment2_mod, fragment3_mod], topology='circular', product='pCMV_Target_AID', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-D3mSzAAq0h', originals=[])
joindna(*[fragment1_mod, fragment2_mod, fragment3_mod], topology='circular', product='pCMV_Target_AID', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-D3mSzAAq0h', originals=['pCMV_Target_AID-D3mSzAAq0h'])

################################################################################
#This source code was auto-generated by 'quine' funtion of QUEEN 1.0.0.
#Project Name    :pCMV_Target_AID
#File Name       :None
#Creation Date   :2021-08-26
################################################################################
project='pCMV_Target_AID_mod'
import sys
sys.path.append("/Users/hideto/Dropbox (Yachie Lab)/HIDETO_MORI.LAB/Experiments/Project/Dbrick/github/demo/tutorial/../..")
from QUEEN.queen import *
from QUEEN import cutsite as cs
set_namespace(globals())

QUEEN(record='https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF', dbtype='benchling', project='pcDNA_Target_AID', product='pcDNA_Target_AID', process_id='pCMV_Target_AID_mod-EZbi1moI0C', originals=['pCMV_Target_AID-EZbi1moI0C'])
QUEEN(record='https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf', dbtype='benchling', project='pCMV_ABE', product='pCMV_ABE', process_id='pCMV_Target_AID_mod-rA2RRtBu3d', originals=['pCMV_Target_AID-rA2RRtBu3d'])

process1={'name':'PCR', 'description':'The N-terminus of Target-AID (fragment1) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'}
QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', project='RS045', product='RS045', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-WpgGb7XWl5', originals=['pCMV_Target_AID-WpgGb7XWl5'])
QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', project='HM129', product='HM129', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-q3VxR7sxc2', originals=['pCMV_Target_AID-q3VxR7sxc2'])
pcDNA_Target_AID.searchsequence(query=RS045.seq[-18:], product='FW1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-ncQKuKboBk', originals=['pCMV_Target_AID-ncQKuKboBk'])
pcDNA_Target_AID.searchsequence(query=HM129.seq[-18:], product='RV1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-Yhy8TvI76F', originals=['pCMV_Target_AID-Yhy8TvI76F'])
cropdna(pcDNA_Target_AID, start=FW1[0].end, end=RV1[0].start, product='extract1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-3AjRaQxViH', originals=['pCMV_Target_AID-3AjRaQxViH'])
modifyends(extract1, left=RS045.seq, right=HM129.rcseq, product='fragment1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID_mod-tm34jxHzZx', originals=['pCMV_Target_AID-tm34jxHzZx'])

process2={'name':'PCR', 'description':'The C-terminus of Target-AID (fragment2) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046.'}
QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', project='HM128', product='HM128', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-Pdu8Zk8Kcu', originals=['pCMV_Target_AID-Pdu8Zk8Kcu'])
QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', project='RS046', product='RS046', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-a96XNwwKYV', originals=['pCMV_Target_AID-a96XNwwKYV'])
pcDNA_Target_AID.searchsequence(query=HM128.seq[-18:], product='FW2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-oBPkpzeOwc', originals=['pCMV_Target_AID-oBPkpzeOwc'])
pcDNA_Target_AID.searchsequence(query=RS046.seq[-18:], product='RV2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-4NvpQbAVb8', originals=['pCMV_Target_AID-4NvpQbAVb8'])
cropdna(pcDNA_Target_AID, start=FW2[0].end, end=RV2[0].start, product='extract2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-r6NI2fSBK4', originals=['pCMV_Target_AID-r6NI2fSBK4'])
modifyends(extract2, left=HM128.seq, right=RS046.rcseq, product='fragment2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID_mod-bjx9zDibMA', originals=['pCMV_Target_AID-bjx9zDibMA'])

process3={'name':'PCR', 'description':'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048.'}
QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', project='RS047', product='RS047', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-xbYwpI6slr', originals=['pCMV_Target_AID-xbYwpI6slr'])
QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', project='RS048', product='RS048', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-uANRoOl8Yj', originals=['pCMV_Target_AID-uANRoOl8Yj'])
pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-ZRzlemHnPL', originals=['pCMV_Target_AID-ZRzlemHnPL'])
pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-OZqgFZtJ07', originals=['pCMV_Target_AID-OZqgFZtJ07'])
cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, product='extract3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-4Syq4DDNMa', originals=['pCMV_Target_AID-4Syq4DDNMa'])
modifyends(extract3, left=RS047.seq, right=RS048.rcseq, product='fragment3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID_mod-nxOQby3rWS', originals=['pCMV_Target_AID-nxOQby3rWS'])

process4={'name':'Gibson Assembly', 'description':'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two insert fragments and a backbone fragments.'}
modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', product='fragment1_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-fNk9Pxad1N', originals=['pCMV_Target_AID-fNk9Pxad1N'])
modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', product='fragment2_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-V9witA8v4C', originals=['pCMV_Target_AID-V9witA8v4C'])
modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment3_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-CGjYNsICes', originals=['pCMV_Target_AID-CGjYNsICes'])
joindna(*[fragment1_mod, fragment2_mod, fragment3_mod], topology='circular', product='pCMV_Target_AID', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID_mod-D3mSzAAq0h', originals=['pCMV_Target_AID-D3mSzAAq0h'])
if __name__ == '__main__':
    check = quine(pCMV_Target_AID, author=None, project=project, _check=True)
    if check == True:
        pCMV_Target_AID.writedna()
