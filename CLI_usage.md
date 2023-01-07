## CLI usage
```
SYNOPSIS
QUEEN [--help] 
      (--protocol_description | --script_description | --feature_description | --dnamap_visualization | --protocolflow_visualization | --cutdna | --cropdna | --flipdna | --joindna) 
      [--input INPUT [INPUT ...]] [--output OUTPUT] [--positions POSITIONS [POSITIONS ...]] [--start START] [--end END] [--attribute ATTRIBUTE] [--query QUERY]
      [--columns COLUMNS [COLUMNS ...]] [--sequence] [--rcseq] [--linebreak LINEBREAK] [--map_view MAP_VIEW]

OPTIONS:
  -h, --help            show this help message and exit

  QUEEN function options: Please choose only one option from the following ones. 
  --get_gbk, -gg        
  			Download a GenBank file with a queried sequence ID or URL of a specified database. 
			When --get_gbk (-gg) is specified, --database (-db), --seqid (-si), and --output (-o) options are valid.
  --protocol_description, -pd
                        Describe the 'Materials and Methods' of the DNA construct in a QUEEN-generated GenBank input. 
			When --protocol_description (-pd) is specified, --input (-i), and --output (-o) options are valid.
  --script_description, -sd
                        Describe the python script to simulate the DNA construction process of a QUEEN-generated GenBank input. 
			When --script_description (-sd) is specified, --input (-i), and --output (-o) options are valid.
  --feature_description, -fd
                        Print a table of the sequence features in a GenBank input. 
			When --feauture_description (-fd) is specified, --input (-i), --attribute (-a), --query(-q), 
			--columns (-c) --sequence (-seq), and --output (-o) options are valid.
  --dnamap_visualization, -dv
                        Generate the annotated DNA sequence map in a GenBank input. 
			When --dnamap_visualization (-dv) is specified, --input (-i), --map_view (-m), --attribute (-a),
			--query(-q), --sequence(-seq), --rcseq(-rs), --linebreak (-lb), and --output (-o) options are valid.
  --protocolflow_visualization, -pv
                        Generate the flow chart representing the DNA construction processes of a QUEEN-generated GenBank input. 
			When --protocolflow_visualization (-pv) is specified, --input (-i) and --output (-o).
  --cutdna,  -cu        
  			Cut the DNA construct of a given GenBank/Fasta input. 
  			When --cutdna (-cu) is specified, --input (-i), --positions (-pos), and --output (-o) options are valid.
  --cropdna, -cr        
  			Extract a partial DNA fragment from a GenBank/Fasta input. 
			When --cropdna (-cr) is specified, --input (-i), --start (-s), --end (-e), and --output (-o) options are valid.
  --flipdna, -fl        
  			Generate the revese complement of a GenBank/Fasta input. 
  			When --flipdna (-f) is specified, --input (-i) and --output (-o) options are valid.
  --joindna, -jo        
  			Join multiple GenBank inputs and generate the single assembled GenBank/Fasta output. 
			When --joindna (-j) is specified, --input (-i) and --output (-o) options are valid.
  
  Argument options for QUEEN functions: Please use the appropriate argument options corresponding to the QUEEN function.
  --database {ncbi,addgene,benchling}, -db {ncbi,addgene,benchling}
                        For '--db ncbi', set the NCBI accession number. 
			For '--db addgene', set plasmid ID. Sometimes different full sequence maps are provided by the depositor 
			and adgene, respectively, for a single plasmid. In this case, please specify the plasmid ID followed by 
			'addgene' or 'depositor' (Ex. 50005:addgene or 50005:depositor) If you set only plasmid ID, the value will 
			be specified as 'plsmidID:addgene'. For 'benchling', set a benchling shaared link.
  --seqid SEQID, -si SEQID
                        Sequence ID for the corresponding database specified by '--database'.
  --input INPUT [INPUT ...], -i INPUT [INPUT ...] 
                        Input file with FASTA or GenBank format. The file type is estimated based on the file extension. 
			The value on stdin can also be used as a input.
  --output OUTPUT, -o OUTPUT
                        Output file. The file type is estimated based on the file extension.
  --positions POSITIONS [POSITIONS ...], -pos POSITIONS [POSITIONS ...]
                        List of cut positions. A cut position should be provided by `int`. 
			For generating sticy-ends, please use the QUEEN functions as python commands instead of this CLI.
  --start START, -s START 
                        Start position of the target range in the GenBank/Fasta input. 
  --end END, -e END     
  			End position of the target range in the GenBank/Fasta input.
  --attribute ATTRIBUTE, -a ATTRIBUTE
                        Attribute type to be searched (feature_id, feature_type, 'qualifier:*', or sequence). 
			If the value is not provided, all sequence features will be to subjected to the operation by the specified command.
  --query QUERY, -q QUERY
                        Sequence features with the attribute values that match to the query will be searched.
  --columns COLUMNS [COLUMNS ...], -c COLUMNS [COLUMNS ...]
                        List of feature attributes to be displayed in the output table. 
			If the value is 'all', it will generate a table for all the attributes held by the sequence features 
			in the GenBank input except for `sequence`.
  --sequence, -seq      
  			If True when --feature_description is specified, the sequence of each feature for its encoded direction will be 
			displayed in the output table. If True when --dnamap_visualization is specified, a color map representing 
			the QUEEN_object sequence will be displayed below the sequence map.
  --rcseq, -rs          
  			If True when --feature_description, the sequence of each feature for its encoded direction will be displayed in 
			the output table.
  --linebreak LINEBREAK, -lb LINEBREAK
                        Sequence length for line break.
  --map_view MAP_VIEW, -m MAP_VIEW
                        Visualization style. ('linear' or 'circular')

```

### Example commands
Please move to the demo/CLI directory and execute the following commands.  
1. Download a GenBank file with a queried sequence ID or URL of a specified database.   
   By executing the following command, the GenBank file of pUC19 can be downloaded.  
   `QUEEN --get_gbk --database addgene --seqid 50005`

2. Describe the 'Materials and Methods' of the DNA construct in a QUEEN-generated GenBank input.  
   `cat input/pCMV-Target-AID.gbk | QUEEN --protocol_description`  
   The below command return the same result.  
   `QUEEN --protocol_description --input input/pCMV-Target-AID.gbk` 
   
   **Output** 
   
   ``` 
   1. The N-terminus half of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair RS045/HM129.
   2. The C-terminus half of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair HM128/RS046.
   3. A backbone fragment was amplified from pCMV-ABE7.10 using the primer pair RS047/RS048.
   4. The three fragments were assembled by Gibson Assembly.
   ```

3. Describe the python script to simulate the DNA construction process of a QUEEN-generated GenBank input.  
   `QUEEN --script_description --input input/pCMV-Target-AID.gbk`  
   
   **Output** 

   ```
   project='pCMV_Target_AID'
   import sys
   sys.path.append("/usr/local/lib/python3.9/site-packages")
   from QUEEN.queen import *
   from QUEEN import cutsite as cs

   pCMV_ABE = QUEEN(record='https://benchling.com/s/seq-K4HkSd2E8WiTAulJUeBf', dbtype='benchling', product='pCMV_ABE', process_id='pCMV_Target_AID-8VZ56URU5AVQ9VYZDJJUEMGH', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   pcDNA31_Target_AID = QUEEN(record='https://benchling.com/s/seq-cfnGDU0Mq8cUwn185LPF', dbtype='benchling', product='pcDNA31_Target_AID', process_id='pCMV_Target_AID-86JH1MQ2J28S2J1AKDNMTO9P', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   RS045 = QUEEN(seq='GAGAGCCGCCACCATGGCACCGAAGAAGAAGCG', supfeature={'feature_type': 'primer_bind', 'qualifier:label': 'RS045'}, product='RS045', process_id='pCMV_Target_AID-7OW3OF1RLWC5448Q9NCN80R9', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   HM129 = QUEEN(seq='CTGGGGCACGATATGATCCACGTCGTAGTCGGAGA', supfeature={'feature_type': 'primer_bind', 'qualifier:label': 'HM129'}, product='HM129', process_id='pCMV_Target_AID-1ARL2MY8QCUM1HKG1XKM984W', original_ids=[], _sourcefile='pCMV_Target_AID_construction')

   process1={'name':'PCR', 'description':'1. The N-terminus half of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair RS045/HM129.'}
   FW1 = pcDNA31_Target_AID.searchsequence(query=RS045.seq[-18:], product='FW1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-4NFVHEFMWRN16QPABJLD7ISH', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   RV1 = pcDNA31_Target_AID.searchsequence(query=HM129.seq[-18:], product='RV1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-5WVF01KPHJRH13M9YT8JSG3P', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   extract1 = cropdna(pcDNA31_Target_AID, start=FW1[0].end, end=RV1[0].start, product='extract1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-A8UBN1SA1NGJ6MRDO04MJN5D', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   fragment1 = modifyends(extract1, left=RS045.seq, right=HM129.rcseq, supfeature={'feature_id': 'f1', 'qualifier:label': 'fragment-1'}, product='fragment1', process_name=process1['name'], process_description=process1['description'], process_id='pCMV_Target_AID-1XSH4HNL9P3W5XDJJ26N8V6B', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
  
   HM128 = QUEEN(seq='CTACGACGTGGATCATATCGTGCCCCAGTCTTTTC', supfeature={'feature_type': 'primer_bind', 'qualifier:label': 'HM128'}, product='HM128', process_id='pCMV_Target_AID-15MOJYKN6VG5B9AE07WDO4E9', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   RS046 = QUEEN(seq='TTTAAACTCATTATAGCATCTTGATCTTGTTCTCTC', supfeature={'feature_type': 'primer_bind', 'qualifier:label': 'RS046'}, product='RS046', process_id='pCMV_Target_AID-7ATW8VGEDIDK3TV6OLPAPNAO', original_ids=[], _sourcefile='pCMV_Target_AID_construction')

   process2={'name':'PCR', 'description':'2. The C-terminus half of Target-AID was amplified from pcDNA3.1_pCMV-nCas-PmCDA1-ugi pH1-gRNA(HPRT) using the primer pair HM128/RS046.'}
   FW2 = pcDNA31_Target_AID.searchsequence(query=HM128.seq[-18:], product='FW2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-7HC669JZUKZW7W9DQYU0EN74', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   RV2 = pcDNA31_Target_AID.searchsequence(query=RS046.seq[-18:], product='RV2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-AD3ESLVWRLTZB4HX3GY338EA', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   extract2 = cropdna(pcDNA31_Target_AID, start=FW2[0].end, end=RV2[0].start, product='extract2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-6K9AZGU8RT9QA2VDROZRK94F', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   fragment2 = modifyends(extract2, left=HM128.seq, right=RS046.rcseq, supfeature={'feature_id': 'f2', 'qualifier:label': 'fragment-2'}, product='fragment2', process_name=process2['name'], process_description=process2['description'], process_id='pCMV_Target_AID-BDKJOCPZRSFR6F268ZFB9X4E', original_ids=[], _sourcefile='pCMV_Target_AID_construction')

   RS047 = QUEEN(seq='ATCAAGATGCTATAATGAGTTTAAACCCGCTGATC', supfeature={'feature_type': 'primer_bind', 'qualifier:label': 'RS047'}, product='RS047', process_id='pCMV_Target_AID-ALZ0SU0JKMDZ42C7JGBVOPDV', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   RS048 = QUEEN(seq='CTTCGGTGCCATGGTGGCGGCTCTCCCTATAG', supfeature={'feature_type': 'primer_bind', 'qualifier:label': 'RS048'}, product='RS048', process_id='pCMV_Target_AID-2D0MSKUXK0M04Z2BVOWFLPQA', original_ids=[], _sourcefile='pCMV_Target_AID_construction')

   process3={'name':'PCR', 'description':'3. A backbone fragment was amplified from pCMV-ABE7.10 using the primer pair RS047/RS048.'}
   FW3 = pCMV_ABE.searchsequence(query=RS047.seq[-18:], product='FW3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-6Z7TJ8ZUETN4R8SGLIXCNHF4', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   RV3 = pCMV_ABE.searchsequence(query=RS048.seq[-18:], product='RV3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-WI28INSGUMCZ2NUOUFKKWIKE', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   extract3 = cropdna(pCMV_ABE, start=FW3[0].end, end=RV3[0].start, product='extract3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-36NO0BR31UOD3CM08BGT6WC7', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   fragment3 = modifyends(extract3, left=RS047.seq, right=RS048.rcseq, supfeature={'feature_id': 'f3', 'qualifier:label': 'fragment-3'}, product='fragment3', process_name=process3['name'], process_description=process3['description'], process_id='pCMV_Target_AID-98SJNGJWWQ3J8B8RAOU4YQF0', original_ids=[], _sourcefile='pCMV_Target_AID_construction')

   process4={'name':'Gibson Assembly', 'description':'4. The three fragments were assembled by Gibson Assembly.'}
   fragment1_mod = modifyends(fragment1, left='*{30}/-{30}', right='-{30}/*{30}', product='fragment1_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-6JNVDR69HUIVEYX7Z85NO11A', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   fragment2_mod = modifyends(fragment2, left='*{30}/-{30}', right='-{30}/*{30}', product='fragment2_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-AFMFFB673FWQBW3GA20LXBFH', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   fragment3_mod = modifyends(fragment3, left='*{30}/-{30}', right='-{30}/*{30}', product='fragment3_mod', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-6867T0Q4KU2M9WCPFDME0HA1', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   pCMV_Target_AID = joindna(*[fragment1_mod, fragment2_mod, fragment3_mod], topology='circular', product='pCMV_Target_AID', process_name=process4['name'], process_description=process4['description'], process_id='pCMV_Target_AID-8IXV92W0BW934FILCZ1TOH7R', original_ids=[], _sourcefile='pCMV_Target_AID_construction')
   if __name__ == '__main__':
       pCMV_Target_AID.outputgbk()
   ```

4. Describe the sequence features of the DNA construct in a QUEEN-generated GenBank input.  
   `QUEEN --feature_description --input input/pCMV-Target-AID.gbk`  
   
   **Output**
   
   ```
   feature_id  feature_type  qualifier:label     start  end   strand
   0           source        source              0      3308  +
   100         primer_bind   M13 Reverse         275    292   -
   200         primer_bind   M13/pUC Reverse     288    311   -
   300         protein_bind  lac operator        299    316   +
   400         promoter      lac promoter        323    354   -
   500         protein_bind  CAP binding site    368    390   +
   600         primer_bind   L4440               506    524   -
   700         rep_origin    ori                 677    1266  -
   800         primer_bind   pBR322ori-F         757    777   -
   900         CDS           AmpR                1436   2297  -
   1000        primer_bind   Amp-R               2059   2079  +
   1100        promoter      AmpR promoter       2297   2402  -
   1200        primer_bind   pRS-marker          2480   2500  -
   1300        enhancer      CMV enhancer        2671   3051  +
   1400        promoter      CMV promoter        3051   3255  +
   1500        primer_bind   CMV-F               3205   3226  +
   1600        promoter      T7 promoter         3296   3315  +
   1700        primer_bind   RS048               3308   3340  -
   1800        primer_bind   RS045               3315   3348  +
   1900        misc_feature  fragment-1          3315   5911  +
   2000        CDS           SV40 NLS            3334   3355  +
   2100        source        source              3348   8678  +
   2200        CDS           Cas9(D10A)          3379   7483  +
   2300        primer_bind   HM129               5876   5911  -
   2400        primer_bind   HM128               5883   5918  +
   2500        misc_feature  fragment-2          5883   8714  +
   2600        CDS           SV40 NLS            7495   7516  +
   2700        CDS           3xFLAG              7723   7789  +
   2800        CDS           PmCDA1              7795   8422  +
   2900        CDS           SV40 NLS            8422   8443  +
   3000        CDS           UGI                 8449   8701  +
   3100        primer_bind   RS046               8678   8714  -
   3200        primer_bind   RS047               8689   8724  +
   3300        misc_feature  fragment-3          8689   3340  +
   3400        source        source              8724   8752  +
   3500        primer_bind   BGH-rev             8726   8744  -
   3600        polyA_signal  bGH poly(A) signal  8732   205   +
   ```
	
   By using --query and --attribute options, you can describe the seauence features holding a queried value in a designated attribute.
   `QUEEN --feature_description --input input/pCMV-Target-AID.gbk --query \.\*NLS --attribute qualifier:label --sequence --separation ,`
   
   **Output** 
   
   ```
   feature_id,feature_type,qualifier:label,start,end,strand,sequence
   2000,CDS,SV40 NLS,3334,3355,+,CCGAAGAAGAAGCGTAAAGTC
   2600,CDS,SV40 NLS,7495,7516,+,CCCAAGAAGAAGAGGAAGGTG
   2900,CDS,SV40 NLS,8422,8443,+,CCCAAGAAGAAAAGAAAAGTC
   ```

5. Change the origin of a circular dna.  
   `QUEEN --cutdna --input input/pCMV-Target-AID.gbk --positions 3379 | QUEEN --joindna > output/pCMV-Target-AID_slided.gbk`  
   The output file is `demo/CLI/output/pCMV-Target-AID_slided.gbk`.

6. Extract the Cas9 fragment and generate its reverse complement.  
   `QUEEN --cropdna --input input/pCMV-Target-AID.gbk --start 3379 --end 7483 | QUEEN --flipdna > output/Cas9_rc.gbk`  
   The output file is [`demo/CLI/output/Cas9_rc.gbk`](https://github.com/yachielab/QUEEN/blob/master/demo/CLI/output/Cas9_rc.gbk).

7. Generate the annotated dna sequence map of the GenBank input.  
   `QUEEN --dnamap_visualization --input input/pCMV-Target-AID.gbk --map_view circular --output output/pCMV-Target-AID_map.pdf`  
   `QUEEN --dnamap_visualization --input output/pCMV-Target-AID_slided.gbk --map_view circular --output output/pCMV-Target-AID_slided_map.pdf`  
   `QUEEN --dnamap_visualization --input output/Cas9_rc.gbk --linebreak 200 --sequence --rcseq --output output/Cas9_rc_map.pdf`  
   The output files are [`demo/CLI/output/pCMV-Target-AID_map.pdf`](https://github.com/yachielab/QUEEN/blob/master/demo/CLI/output/pCMV-Target-AID_map.pdf), [`demo/CLI/output/pCMV-Target-AID_map_slided.pdf`](https://github.com/yachielab/QUEEN/blob/master/demo/CLI/output/pCMV-Target-AID_slided_map.pdf), and [`demo/output/Cas9_rc_map.pdf`](https://github.com/yachielab/QUEEN/blob/master/demo/CLI/output/Cas9_rc_map.pdf).  

8. Generate the plotocol flow to construct a QUEEN-generated GenBank input.  
   `QUEEN --protocolflow_visualization --input input/pCMV-Target-AID.gbk --output output/pCMV-Target-AID_flow.pdf"`  
   The output file is [`demo/CLI/output/pCMV-Target-AID_flow.pdf`](https://github.com/yachielab/QUEEN/blob/master/demo/CLI/output/pCMV-Target-AID_flow.pdf).  

