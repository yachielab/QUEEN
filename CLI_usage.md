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
  --protocol_description, -pd
                        Describe the 'Materials and Methods' of the DNA construct in a QUEEN-generated GenBank input. 
			When --protocol_description (-pd) is specified, --input (-i) and --output (-o) options are valid.
  --script_description, -sd
                        Describe the python script to simulate the DNA construction process of a QUEEN-generated GenBank input. 
			When --script_description (-sd) is specified, --input (-i) and --output (-o) options are valid.
  --feature_description, -fd
                        Print a table of the sequence features in a GenBank input. 
			When --feauture_description (-fd) is specified, --input (-i), --attribute (-a), --query(-q), --columns (-c) --sequence (-seq), and --output (-o) options are valid.
  --dnamap_visualization, -dv
                        Generate the annotated DNA sequene map in a GenBank input. 
			When --dnamap_visualization (-dv) is specified, --input (-i), --map_view (-m), --attribute (-a), --query(-q), --sequence(-seq), --rcseq(-rs), --linebreak (-lb), and --output (-o) options are valid.
  --protocolflow_visualization, -pv
                        Generate the flow chart representing the DNA construction processes of a QUEEN-generated GenBank input. When --protocolflow_visualization (-pv) is specified, --input (-i) and --output (-o).
  --cutdna,  -cu        
  			Cut the DNA construct given of a GenBank/Fasta input. 
  			When --cropdna (-c) is specified, --input (-i), --positions (-pos), and --output (-o) options are valid.
  --cropdna, -cr        
  			Extract a partial DNA fragment from a GenBank/Fasta input. 
			When --cropdna (-c) is specified, --input (-i), --start (-s), --end (-e), and --output (-o) options are valid.
  --flipdna, -fl        
  			Generate the revese complement of a GenBank/Fasta input. 
  			When --flipdna (-f) is specified, --input (-i) and --output (-o) options are valid.
  --joindna, -jo        
  			Join multiple GenBank inputs and generate the single assembled GenBank/Fasta output. 
			When --joindna (-j) is specified, --input (-i) and --output (-o) options are valid.
  
  Aragument options for QUEEN functions:
  --input INPUT [INPUT ...], -i INPUT [INPUT ...] 
                        Input file with FASTA or GenBank format. The file type is estimated based on the file extension. The value on stdin can also be used as a input.
  --output OUTPUT, -o OUTPUT
                        Output file. The file type is estimated based on the file extension.
  --positions POSITIONS [POSITIONS ...], -pos POSITIONS [POSITIONS ...]
                        List of cut positions. A cut position should be provided by `int`.For generating sticy-ends, please use the QUEEN functions as python commands instead of this CLI.
  --start START, -s START 
                        Start position of the target range in the GenBank/Fasta input. 
  --end END, -e END     
  			End position of the target range in the GenBank/Fasta input.
  --attribute ATTRIBUTE, -a ATTRIBUTE
                        Attribute type to be searched (feature_id, feature_type, 'qualifier:*', or sequence). If the value is not provided, all sequence features will be to subjected to the operation by the specified command.
  --query QUERY, -q QUERY
                        Sequence features with the the attribute values that match to the query will be searched.
  --columns COLUMNS [COLUMNS ...], -c COLUMNS [COLUMNS ...]
                        List of feature attributes to be displayed in the output table. If the value is 'all', it will generate a table for all the attributes held by the sequence features in the GenBank input except for `sequence`.
  --sequence, -seq      
  			If True when --feature_description is specified, the sequence of each feature for its encoded direction will be displayed in the output table. If True when --dnamap_visualization is specified, a color map representing the QUEEN_object sequence will be displayed below the sequence map.
  --rcseq, -rs          
  			If True when --feature_description, the sequence of each feature for its encoded direction will be displayed in the output table.
  --linebreak LINEBREAK, -lb LINEBREAK
                        Sequence length for line break.
  --map_view MAP_VIEW, -m MAP_VIEW
                        Visualization style. ('linear' or 'circular')

```

### Example commands
Please move to the demo/CLI directory and execute the following commands.

1. Describe the 'Materials and Methods' of the DNA construct in a QUEEN-generated GenBank input.  
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

2. Describe the sequence features of the DNA construct in a QUEEN-generated GenBank input
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

3. Change the origin of a circular dna. 
   `QUEEN --cutdna --input input/pCMV-Target-AID.gbk --positions 3379 | QUEEN --joindna > output/pCMV-Target-AID_slided.gbk`

4. Extract the Cas9 fragment and generate its reverse complement.
   `QUEEN --cropdna --input input/pCMV-Target-AID.gbk --start 3379 --end 7483 | QUEEN --flipdna > output/Cas9_rc.gbk`

5. Generate the annotated dna sequence map of the GenBank input.
   `QUEEN --dnamap_visualization --input input/pCMV-Target-AID.gbk --output output/pCMV-Target-AID_map.pdf`
   `QUEEN --dnamap_visualization --input output/pCMV-Target-AID_slided.gbk --output output/pCMV-Target-AID_map_slided.pdf`
