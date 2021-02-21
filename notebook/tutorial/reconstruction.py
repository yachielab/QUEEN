from dna import *
DNA.dna_dict['pUC19'] = DNA(seq=None, record='pUC19.gbk', project='pUC19', topology='linear', format=None)
DNA.dna_dict['EGFP'] = DNA(seq=None, record='EGFP.fasta', project='EGFP', topology='linear', format=None)
DNA.dna_dict['pUC19_0'] = cropdna(DNA.dna_dict['pUC19'], start=(681, 677), end=(632, 636), project='pUC19')
DNA.dna_dict['EGFP_0'] = editdna(DNA.dna_dict['EGFP'], key_attribute='sequence', query=None, target_attribute='feature ID', operation=createdna(value='Ins_1'), comment='', destructive=False, project='EGFP')
DNA.dna_dict['EGFP_1'] = editdna(DNA.dna_dict['EGFP_0'], key_attribute='feature ID', query='Ins_1', target_attribute='qualifier:label', operation=createdna(value='EGFP'), comment='', destructive=False, project='EGFP')
DNA.dna_dict['EGFP_2'] = modifyends(DNA.dna_dict['EGFP_1'], left='AGCTT/----A', right='GAGCT/C----', project='EGFP')
DNA.dna_dict['pUC19-EGFP'] = joindna(*[DNA.dna_dict['EGFP_2'],DNA.dna_dict['pUC19_0']], topology='circular', project='pUC19-EGFP')
DNA.dna_dict['pUC19-EGFP'].writedna('reconstructed_pUC19-EGFP.gbk')
