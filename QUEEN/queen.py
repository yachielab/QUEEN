#!/usr/bin/env python
import os 
import io
import functools
import sys
sys.path.append(os.path.dirname(__file__))
from qobj import QUEEN
import qfunction
import quine
import qgraph 
import argparse
from Bio import SeqIO

__version__ = "1.2.0"

_namespace       = None
QUEEN._namespace = {}
replaceattribute = qfunction.replaceattribute
createattribute  = qfunction.createattribute
removeattribute  = qfunction.removeattribute
cutdna           = qfunction.cutdna
cropdna          = qfunction.cropdna
modifyends       = qfunction.modifyends
joindna          = qfunction.joindna
flipdna          = qfunction.flipdna
editsequence     = qfunction.editsequence
editfeature      = qfunction.editfeature
visualizemap     = qfunction.visualizemap
visualizeflow    = qgraph.visualizeflow
quine            = quine.quine


commands = [
            "protocol_description", 
            "script_description", 
            "feature_description",
            "dnamap_visualization",
            "protcolflow_visualization",
            "cropdna",
            "flipdna",
            "joindna"
            ]

command_abbr_dict = {
                     "protcol_description"  :"pd",
                     "script_description"   :"sd",
                     "feature_description"  :"fd", 
                     "dnamap_visualization ":"dv", 
                     "protocolflow_visualization":"pv",
                     "cutdna":"c",
                     "cropdna":"e",
                     "flipdna":"f",
                     "joindna":"j"
                     }

abbr_command_dict = dict(zip(*[list(command_abbr_dict.values()), list(command_abbr_dict.keys())]))

def set_namespace(_globals=None):
    if _globals is None:
        _namespace = None
        QUEEN._namespace     = {}
        QUEEN._namespaceflag = 0
    else:
        _namespace = _globals
        QUEEN._namespace     = _globals
        QUEEN._namespaceflag = 1

def main(args):
    if args.joindna and len(args.input) > 1:
        qinput = [] 
        for ainput in args.input:
            content = ainput.read()
            if content[0] == ">":
                aqinput = QUEEN(record=SeqIO.read(io.StringIO(content), format="fasta"))
            else:
                aqinput = QUEEN(record=SeqIO.read(io.StringIO(content), format="genbank"))
            qinput.append(aqinput) 
    else:
        content = args.input[0].read()
        if content[0] == ">":
            qinput = QUEEN(record=SeqIO.read(io.StringIO(content), format="fasta"))
        else:
            qinput = QUEEN(record=SeqIO.read(io.StringIO(content), format="genbank"))

    if args.protocol_description:
        quine(qinput, output=args.output, process_description=True) 

    elif args.script_description:
        quine(qinput, output=args.output) 

    elif args.feature_description:
        output = args.output
        seq    = args.sequence   
        if output is None:
            separation = args.separation 
        else:
            extension = output.split(".")[-1] 
            if extension == "csv":
                if args.separation is None:
                    separation = ","
                else:
                    separation = args.separation
            
            elif extension == "tsv": 
                if args.separation is None:
                    separation = "\t"
                else:
                    separation = args.separation
            else:
                separation = args.separation

        if args.attribute == "all" and args.query == ".+":
            features = qinput.dnafeatures
        else:
            features = qinput.searchfeature(key_attribute=args.attribute, query=args.query)
        qinput.printfeature(features, attribute=args.columns, seq=seq, separation=separation, output=args.output)   
    
    elif args.dnamap_visualization:
        if args.attribute == "all" and args.query == ".+":
            features = None
        else:
            features = qinput.searchfeature(key_attribute=args.attribute, query=args.query)
        
        if args.map_view == "linear":
            fig = visualizemap(qinput, map_view=args.map_view, feature_list=features, linebreak=args.linebreak, seq=args.sequence, rcseq=args.rcseq)  
        else:
            fig = visualizemap(qinput, map_view=args.map_view, feature_list=features)  

        if args.output is None:
            fig.savefig(qinput.project + "_map.pdf") 
        else:
            fig.savefig(args.output) 

    elif args.protocolflow_visualization:
        graph     = visualizeflow(qinput)
        if args.output is None:
            output    = qinput.project + "_flow"
            extension = "pdf"
        elif "." in args.output:
            extension = args.output.split(".")[-1]
            output = ".".join(args.output.split(".")[:-1])
        else:
            output = args.output
            extension = "pdf"

        if extension == "pdf":
            pass 
        else: 
            pass 
        graph.render(output) 
    else:
        output = args.output
        if output is None:
            format="genbank"
        else:
            if "." in output:
                extension = output.split(".")[-1] 
                if extension in ("fasta", "fna"):
                    format = "fasta"
                else:
                    format = "genbank" 
            else:
                format    = "genbank"
                extension = ".gbk"
                output    = output + extension

        if args.cutdna:
            positions = [int(pos) for pos in args.positions] 
            qoutputs  = cutdna(qinput, *positions) 
            
            if len(qoutputs) > 1:
                for i, aqoutput in enumerate(qoutputs):
                    if output is None:
                        aqoutput.outputgbk()
                    else:
                        suboutput = ".".join(output.split(".")[:-1]) + "_" + str(i) + extension
                        aqoutput.outputgbk(suboutput)
            else:
                qoutputs[0].outputgbk(output, format=format)
            
        elif args.cropdna:
            qoutput = cropdna(qinput, start=args.start, end=args.end, product="hoge") 
            qoutput.outputgbk(output, format=format)

        elif args.flipdna:
            qoutput = flipdna(qinput) 
            qoutput.outputgbk(output, format=format)

        elif args.joindna:
            if type(qinput) == list: 
                qoutput = joindna(*qinput)
            else:
                qoutput = joindna(qinput) 
            qoutput.outputgbk(output, format=format)


if __name__ == "__main__": 
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)

    g.add_argument("--protocol_description", "-pd", action="store_true",
                    help="Describe the 'Materials and Methods' of the DNA construct in a QUEEN-generated GenBank input.\n" + 
                    "When --protocol_description (-pd) is specified, --input (-i) and --output (-o) options are valid.")
    
    g.add_argument("--script_description", "-sd", action="store_true",
                    help="Describe the python script to simulate the DNA construction process of a QUEEN-generated GenBank input.\n" +  
                    "When --script_description (-sd) is specified, --input (-i) and --output (-o) options are valid.")
    
    g.add_argument("--feature_description", "-fd", action="store_true",
                    help="Print a table of the sequence features in a GenBank input.\n" + 
                    "When --feauture_description (-fd) is specified, --input (-i), --attribute (-a), --query(-q), --columns (-c) " + 
                    "--sequence (-seq), and --output (-o) options are valid.")
                    
    g.add_argument("--dnamap_visualization",  "-dv", action="store_true",
                    help="Generate the annotated DNA sequene map in a GenBank input.\n" +
                    "When --dnamap_visualization (-dv) is specified, --input (-i), --map_view (-m), --attribute (-a), --query(-q), " + 
                    "--sequence(-seq), --rcseq(-rs), --linebreak (-lb), and --output (-o) options are valid.") 

    g.add_argument("--protocolflow_visualization", "-pv", action="store_true",
                    help="Generate the flow chart representing the DNA construction processes of a QUEEN-generated GenBank input.\n" +
                    "When --protocolflow_visualization (-pv) is specified, --input (-i) and --output (-o).")
    
    g.add_argument("--cutdna", "-cu", action="store_true",
                    help="Cut the DNA construct given of a GenBank/Fasta input.\n" +
                    "When --cropdna (-c) is specified, --input (-i), --positions (-pos), and --output (-o) options are valid.")

    g.add_argument("--cropdna", "-cr", action="store_true",
                    help="Extract a partial DNA fragment from a GenBank/Fasta input.\n" +
                    "When --cropdna (-c) is specified, --input (-i), --start (-s), --end (-e), and --output (-o) options are valid.")
    
    g.add_argument("--flipdna", "-fl", action="store_true",
                    help="Generate the revese complement of a GenBank/Fasta input.\n" + 
                    "When --flipdna (-f) is specified, --input (-i) and --output (-o) options are valid.") 

    g.add_argument("--joindna", "-jo", action="store_true",
                    help="Join multiple GenBank inputs and generate the single assembled GenBank/Fasta output.\n" + 
                    "When --joindna (-j) is specified, --input (-i) and --output (-o) options are valid.")

    p.add_argument("--input", "-i", type=argparse.FileType('r'), nargs="+", default=[sys.stdin,] , 
                   help="Input file with FASTA or GenBank format. The file type is estimated based on the file extension.\n" +
                   "The value on stdin can also be used as a input.")

    p.add_argument("--output", "-o", type=str, default=None, 
                   help="Output file. The file type is estimated based on the file extension.")
    
    p.add_argument("--separation", "-sep", type=str, default=None, 
                    help="String to separate values of each line. If the value is not given, space(s) to generate a well-formatted table is used as separators.")
    
    p.add_argument("--positions", "-pos", nargs = "+", 
                    help="List of cut positions. A cut position should be provided by `int`." +  
                    "For generating sticy-ends, please use the QUEEN functions as python commands instead of this CLI.")

    p.add_argument("--start", "-s", type=int, default=0, 
                   help="Start position of the target range in the GenBank/Fasta input.")

    p.add_argument("--end", "-e", type=int, default=None, 
                   help="End position of the target range in the GenBank/Fasta input.")

    p.add_argument("--attribute", "-a", type=str, default="all",
                   help="Attribute type to be searched (feature_id, feature_type, 'qualifier:*', or sequence).\n" +  
                   "If the value is not provided, all sequence features will be to subjected to the operation by the specified command.")
    
    p.add_argument("--query", "-q", type=str, default=".+",
                   help="Sequence features with the the attribute values that match to the query will be searched.")

    p.add_argument("--columns", "-c", nargs = "+", default=["feature_id", "feature_type", "qualifier:label", "start", "end", "strand"], 
                   help="List of feature attributes to be displayed in the output table.\nIf the value is 'all', it will generate " + 
                   "a table for all the attributes held by the sequence features in the GenBank input except for `sequence`.")
    
    p.add_argument("--sequence", "-seq", action="store_true", 
                   help="If True when --feature_description is specified, the sequence of each feature for its encoded direction will be displayed in the output table.\n" + 
                   "If True when --dnamap_visualization is specified, a color map representing the QUEEN_object sequence will be displayed below the sequence map." )
 
    p.add_argument("--rcseq", "-rs", action="store_true", 
                   help="If True when --feature_description, the sequence of each feature for its encoded direction will be displayed in the output table.")

    p.add_argument("--linebreak", "-lb", type=int, default=None,
                   help="Sequence length for line break.")

    p.add_argument("--map_view", "-m", default="linear", 
                   help="Visualization style. ('linear' or 'circular')")

    args = p.parse_args()

    main(args) 

