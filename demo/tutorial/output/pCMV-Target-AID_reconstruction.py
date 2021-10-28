################################################################################
#This source code was auto-generated by 'quine' funtion of QUEEN 1.0.0.
#Project Name    :pCMV_Target_AID_0
#File Name       :output/pCMV-Target-AID_reconstruction.py
#Creation Date   :2021-10-08
################################################################################
project='pCMV_Target_AID'
import sys
sys.path.append("/usr/local/lib/python3.9/site-packages/QUEEN-1.0.0-py3.9.egg")
from QUEEN.queen import *
from QUEEN import cutsite as cs
set_namespace(globals())

QUEEN(record='https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF', dbtype='benchling', project='pcDNA31_Target_AID', product='pcDNA31_Target_AID', process_id='pCMV_Target_AID-5d57650744fe2125d58dbc 6b7afb5fb4Dj02', originals=[])
QUEEN(record='https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf', dbtype='benchling', project='pCMV_ABE', product='pCMV_ABE', process_id='pCMV_Target_AID-73a4416a2922b40d0fafe8 7a9fda17afrA2R', originals=[])

process1={'name':'PCR', 'description':'The N-terminus of Target-AID (fragment1) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129.'}
QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', project='RS045', product='RS045', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-ca7e159835614b8a5dad7a 5423b33674WpgG', originals=[])
QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', project='HM129', product='HM129', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-36a68173712e1050ac417a 023534e59cq3Vx', originals=[])
pcDNA31_Target_AID.searchsequence(query=RS045.seq[-18:], product='FW1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-e8d4c19db976971fddbe0d 2135a2e575ncQK', originals=[])
pcDNA31_Target_AID.searchsequence(query=HM129.seq[-18:], product='RV1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-b909740d176b2aeeb2834c 9fd0c83a35Yhy8', originals=[])
cropdna(pcDNA31_Target_AID, start=FW1[0].end, end=RV1[0].start, product='extract1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-069014332f922d046f6ec8 8d20bbab7c3AjR', originals=[])
modifyends(extract1, left=RS045.seq, right=HM129.rcseq, product='fragment1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-f6381ccacc14fbfe3827df 4fb1bacb1btm34', originals=[])

process2={'name':'PCR', 'description':'The C-terminus of Target-AID (fragment2) was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046.'}
QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', project='HM128', product='HM128', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-129f3256cfa215011fa93 b5350fe917bPdu8', originals=[])
QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', project='RS046', product='RS046', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-b44d5e5fce0e2f3427caa 1a6abc78c92a96X', originals=[])
pcDNA31_Target_AID.searchsequence(query=HM128.seq[-18:], product='FW2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-3f475063e0f0e06712825 d4071309ccaoBPk', originals=[])
pcDNA31_Target_AID.searchsequence(query=RS046.seq[-18:], product='RV2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-d00950e1ddeb5a39340d9 fe6c9ca60c04Nvp', originals=[])
cropdna(pcDNA31_Target_AID, start=FW2[0].end, end=RV2[0].start, product='extract2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-0f864041dfaf0992268be 7e5fa72fdc3r6NI', originals=[])
modifyends(extract2, left=HM128.seq, right=RS046.rcseq, product='fragment2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-d18f41fd7315881066972 5867b75366cbjx9', originals=[])

process3={'name':'PCR', 'description':'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048.'}
QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', project='RS047', product='RS047', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-36b34c02bbe227b542b50 0f9f1018102xbYw', originals=[])
QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', project='RS048', product='RS048', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-41ae0b32e66f4ea83b404 1fed64bfc5buANR', originals=[])
pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-52d627a83d069fa426db0 8fffa18caf6ZRzl', originals=[])
pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-829186acdc21227e57c96 ca2e79d971eOZqg', originals=[])
cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, product='extract3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-54685a74d7e9137742c6a 163a2b0aa124Syq', originals=[])
modifyends(extract3, left=RS047.seq, right=RS048.rcseq, product='fragment3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-3a015c0b57a9a18a6c664 62041f4378fnxOQ', originals=[])

process4={'name':'Gibson Assembly', 'description':'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two insert fragments and a backbone fragments.'}
modifyends(fragment1, left='*{25}/-{25}', right='-{28}/*{28}', product='fragment1_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-deece5a05d671ff482c4f 4d1a11d71c3fNk9', originals=[])
modifyends(fragment2, left='*{28}/-{28}', right='-{25}/*{25}', product='fragment2_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-984ce7f58f32295866372 bb65de36029V9wi', originals=[])
modifyends(fragment3, left='*{25}/-{25}', right='-{25}/*{25}', product='fragment3_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-1eb2fbdae1a20f0f2e7ae b8546f3f15cCGjY', originals=[])
joindna(*[fragment1_mod, fragment2_mod, fragment3_mod], topology='circular', product='pCMV_Target_AID', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-a769896addd3976e289a9 d8903d7e5d8D3mS', originals=[])
if __name__ == '__main__':
    check = quine(pCMV_Target_AID, author=None, project=project, _check=True)
    if check == True:
        pCMV_Target_AID.writedna()