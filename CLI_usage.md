## CLI usage
```
SYNOPSIS
QUEEN [--help] 
      (--protocol_description | --script_description | --feature_description | --dnamap_visualization | --protocolflow_visualization | --cutdna | --cropdna | --flipdna | --joindna) 
      [--input INPUT [INPUT ...]] [--output OUTPUT] [--positions POSITIONS [POSITIONS ...]] [--start START] [--end END] [--attribute ATTRIBUTE] [--query QUERY]
      [--columns COLUMNS [COLUMNS ...]] [--sequence] [--rcseq] [--linebreak LINEBREAK] [--map_view MAP_VIEW]

OPTIONS:
  -h, --help            show this help message and exit

  QUEEN function options:
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
  
  Aragument options for QUEEN functions
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
