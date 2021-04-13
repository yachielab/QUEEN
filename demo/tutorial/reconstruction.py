from dna import *
DNA.dna_dict['pCMV_ABE'] = DNA(seq=None, record='input/addgene_102919.gbk', project='pCMV_ABE', topology='linear', format=None, process_description=None)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'] = DNA(seq=None, record='input/addgene_79620.gbk', project='pCMV-nCas-PmCDA1-ugi', topology='linear', format=None, process_description=None)
description0 = 'The N-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs RS045/HM129'
DNA.dna_dict['RS045'] = DNA(seq=None, record='input/RS045.fasta', project='RS045', topology='linear', format=None, process_description=description0)
DNA.dna_dict['HM129'] = DNA(seq=None, record='input/HM129.fasta', project='HM129', topology='linear', format=None, process_description=description0)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_0'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='3914/3914', end='6442/6442', project='pCMV-nCas-PmCDA1-ugi', process_description=description0)
DNA.dna_dict['HM129_0'] = flipdna(DNA.dna_dict['HM129'], project='HM129', process_description=description0)
DNA.dna_dict['N-term_Target-AID'] = joindna(*[DNA.dna_dict['RS045'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_0'], DNA.dna_dict['HM129_0']], topology='linear', project='N-term_Target-AID', process_description=description0)
description1 = 'The C-terminus of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) (Addgene 79620) using primer pairs HM128/RS046'
DNA.dna_dict['HM128'] = DNA(seq=None, record='input/HM128.fasta', project='HM128', topology='linear', format=None, process_description=description1)
DNA.dna_dict['RS046'] = DNA(seq=None, record='input/RS046.fasta', project='RS046', topology='linear', format=None, process_description=description1)
DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'] = cropdna(DNA.dna_dict['pCMV-nCas-PmCDA1-ugi'], start='6484/6484', end='9244/9244', project='pCMV-nCas-PmCDA1-ugi', process_description=description1)
DNA.dna_dict['RS046_0'] = flipdna(DNA.dna_dict['RS046'], project='RS046', process_description=description1)
DNA.dna_dict['C-term_Target-AID'] = joindna(*[DNA.dna_dict['HM128'], DNA.dna_dict['pCMV-nCas-PmCDA1-ugi_1'], DNA.dna_dict['RS046_0']], topology='linear', project='C-term_Target-AID', process_description=description1)
description2 = 'The backbone fragment was amplified from pCMV-ABE7.10 using RS047/RS048'
DNA.dna_dict['RS047'] = DNA(seq=None, record='input/RS047.fasta', project='RS047', topology='linear', format=None, process_description=description2)
DNA.dna_dict['RS048'] = DNA(seq=None, record='input/RS048.fasta', project='RS048', topology='linear', format=None, process_description=description2)
DNA.dna_dict['pCMV_ABE_0'] = cropdna(DNA.dna_dict['pCMV_ABE'], start='5277/5277', end='8613/8613', project='pCMV_ABE', process_description=description2)
DNA.dna_dict['RS048_0'] = flipdna(DNA.dna_dict['RS048'], project='RS048', process_description=description2)
DNA.dna_dict['backbone'] = joindna(*[DNA.dna_dict['RS047'], DNA.dna_dict['pCMV_ABE_0'], DNA.dna_dict['RS048_0']], topology='linear', project='backbone', process_description=description2)
description3 = 'The Target-AID plasmid (pCMV-Target-AID) was constructed by assembling two insert fragments and a backbone fragments.'
DNA.dna_dict['N-term_Target-AID_0'] = modifyends(DNA.dna_dict['N-term_Target-AID'], left='*{25}/-{25}', right='-{28}/*{28}', project='N-term_Target-AID', process_description=description3)
DNA.dna_dict['C-term_Target-AID_0'] = modifyends(DNA.dna_dict['C-term_Target-AID'], left='*{28}/-{28}', right='-{25}/*{25}', project='C-term_Target-AID', process_description=description3)
DNA.dna_dict['backbone_0'] = modifyends(DNA.dna_dict['backbone'], left='*{25}/-{25}', right='-{25}/*{25}', project='backbone', process_description=description3)
DNA.dna_dict['pCMV-Target-AID'] = joindna(*[DNA.dna_dict['N-term_Target-AID_0'], DNA.dna_dict['C-term_Target-AID_0'], DNA.dna_dict['backbone_0']], topology='circular', project='pCMV-Target-AID', process_description=description3)
DNA.dna_dict['pCMV-Target-AID'].writedna('reconstructed_pCMV-Target-AID.gbk')
