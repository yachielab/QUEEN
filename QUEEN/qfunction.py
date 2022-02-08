import os 
import io 
import sys
import copy
import math
import random
import string
import hashlib
import collections
import functools
import warnings
import sre_parse
import regex as re
import numpy as np 

from Bio import SeqIO, pairwise2, BiopythonParserWarning
try:
    from Bio import Alphabet #For Biopython < 1.78
except ImportError:
    Alphabet = False #For Biopyhton ≥ 1.78 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, FeatureLocation, ExactPosition
warnings.simplefilter('ignore', BiopythonParserWarning)

sys.path.append("/".join(__file__.split("/")[:-1]))
import cutsite as cs
import visualize_circular_dna as vc
import visualize_linear_dna as vl 
import qgraph as qg 

def _assigndnafeatures(dnafeatures):
    features = [] 
    for feat in dnafeatures:
        if feat.location.start.position == -1:
            pass 
        else:
            features.append(feat) 
    return features 

def _slide(feats,slide):
    new_feats = []
    for feat in feats:
        feat = copy.deepcopy(feat)
        strand = feat.location.strand
        for p in range(len(feat.location.parts)):
            feat.location.parts[p]._start = ExactPosition(feat.location.parts[p].start.position + slide)
            feat.location.parts[p]._end   = ExactPosition(feat.location.parts[p].end.position + slide)
        feat.location.strand = strand
        new_feats.append(feat.__class__(feat))
    return new_feats 

def _combine_history(dna, history_features):
    history_feature = SeqFeature(FeatureLocation(0, len(dna.seq), strand=0), type="source")
    history_feature = history_features[0].__class__(history_feature, subject=dna.seq)
    for feat in history_features:
        for key in feat.qualifiers:
            if "building_history" in key and feat.qualifiers[key] not in history_feature.qualifiers.values():
                history_feature.qualifiers[key] = feat.qualifiers[key]  
    return history_feature

def add_history(dna, histories, _sourcefile=None):  
    if _sourcefile is not None and sys.argv[0].split("/")[-1].rstrip(".py") != _sourcefile.rstrip(".py"):
        histories[1] = histories[1] + "; _source: {}".format(_sourcefile) 

    elif dna.__class__._source is not None and dna.__class__._source != "__main__":
        histories[1] = histories[1] + "; _source: {}".format(dna.__class__._source) 

    dna.__class__._num_history += 1 
    dna._history_feature.qualifiers["building_history_{}".format(dna.__class__._num_history)] = [] 
    for h, history in enumerate(histories):
        if h == 0:
            dna._history_feature.qualifiers["building_history_{}".format(dna.__class__._num_history)].append(history.replace(" ","–")) 
        else:
            dna._history_feature.qualifiers["building_history_{}".format(dna.__class__._num_history)].append(history)
    dna._features_dict  = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))
    dna.record.features = dna.dnafeatures
    if dna.__class__._keep == 1 and "_product_id" in dna.__dict__: 
        dna.__class__._products[dna._product_id] = dna
        dna._productids.append(dna._product_id) 

def archivehistory(dna):
    for feat in dna.dnafeatures:
        if feat.type == "source":
            old_keys = [] 
            qualifier_keys = list(feat.qualifiers.keys()) 
            for key in qualifier_keys:
                if "building_history" in key:
                    old_keys.append(key) 
                    history = feat.qualifiers[key]
                    new_key = "archived_" + key
                    feat.qualifiers[new_key] = history
            for key in old_keys:
                del feat.qualifiers[key] 

def deletehistory(dna):
    #new_features = [] 
    for feat in dna.dnafeatures:
        flag = 0 
        if feat.type == "source":
            old_keys = [] 
            for key in feat.qualifiers:
                if "building_history" in key or "broken_feature" in key:
                    flag = 1
                    old_keys.append(key)   
            for key in old_keys:
                del feat.qualifiers[key] 
        #if flag == 1:
        #    new_features.append(feat) 
    
def make_processid(dna, chars, process_id=None, original_ids=None):
    new_chars = chars
    for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", chars):
        key = key.group(0)
        new_chars = new_chars.replace(key, "queen")
    
    for key in re.finditer("QUEEN.queried_features_dict\['[^\[\]]+'\]", chars):
        key = key.group(0)
        new_chars = new_chars.replace(key, "dnafeature")
    
    id1 = int(hashlib.md5(chars.encode('UTF-8')).hexdigest(), 16)
    id2 = int(hashlib.md5(new_chars.encode('UTF-8')).hexdigest(), 16)
    id1 = str(np.base_repr(int(id1), base=36))[:12] 
    id2 = str(np.base_repr(int(id2), base=36))[:12]
    newprocess_id = id1 + id2 #str(dna.__class__._num_history)
    if process_id is not None: #and process_id not in dna.__class__._processes:
        if original_ids is None:
            original_ids = [] 
        
        if "-" in process_id:
            original_id = process_id.split("-")[1]
        else:
            original_id = process_id 

        #oid1, oid2, oid3 = original_id.split(":")
        if newprocess_id == original_id:
            newprocess_id = process_id
        else: 
            original_ids.append(process_id)
    else:
        pass
        #newprocess_id = id1 + ':' + id2 + ':' + str(dna.__class__._num_history) 

    dna.__class__._processes[newprocess_id] = new_chars 
    dna._processids.append(newprocess_id) 
    return newprocess_id, original_ids

def compile_cutsite(query):  
    cutsite = query.upper() 
    re_format1 = re.compile(r"([ATGCRYKMSWBDHVN]+)(\([\-0-9]+/[\-0-9]+\))")
    re_format2 = re.compile(r"(\([\-0-9]+/[\-0-9]+\))([ATGCRYKMSWBDHVN]+)(\([\-0-9]+/[\-0-9]+\))")
    re_format3 = re.compile(r"(\([\-0-9]+/[\-0-9]+\))([ATGCRYKMSWBDHVN]+)")    
    match1 = re_format1.fullmatch(query)
    match2 = re_format2.fullmatch(query)
    match3 = re_format3.fullmatch(query)
    if  match1:
        seq = match1.group(1)
        topr, bottomr = map(int, match1.group(2)[1:-1].split("/")) 
        topl, bottoml = "null", "null"
        
    elif match2:
        seq = match2.group(2) 
        topl, bottoml = map(int, match2.group(1)[1:-1].split("/")) 
        topr, bottomr = map(int, match2.group(3)[1:-1].split("/")) 

    elif match3: 
        seq = match3.group(2)
        topl, bottoml = map(int, match3.group(1)[1:-1].split("/")) 
        topr, bottomr = "null", "null"

    else:
        re_format = re.compile("(N*)([^N].+[^N])(N*)")
        match     = re_format.fullmatch(query) 
        nl, nr    = len(match.group(1)), len(match.group(3))
        seq       = match.group(2).replace("^","").replace("_","")  

        query_top    = query.replace("_","")
        query_bottom = query.replace("^","")
        topl = "null" 
        topr = "null" 
        bottoml = "null"
        bottomr = "null"

        count = 0 
        for t, char in enumerate(query_top):
            if count == 0 and char == "^":
                topl = -1 * t + nl
            if count == 1 and char == "^":
                topr = -1 * (len(query_top) - 1 - t) + nr
        
        count = 0  
        for t, char in enumerate(query_bottom):
            if count == 0 and char == "_":
                bottoml = -1 * t + nl 
            if count == 1 and char == "_":
                brromrr = -1 * (len(query_bottom) - 1 - t) + nr 
    return cutsite, seq, topl, topr, bottoml, bottomr

def cutdna(dna, *positions, crop=False, project=None, product=None, process_name=None, process_description=None, pn=None, pd=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1):   
    #Set process name, description and ID
    project = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description

    dna = copy.deepcopy(dna)
    #if process_description is None:
    #    process_description = dna.__class__.process_description
    #else:
    #    dna.__class__.process_description = process_description 
   
    def extract(dna, start, end, project=None): 
        start_top    = start[0] 
        start_bottom = start[1] 
        start = min(start)
        
        end_top    = end[0]
        end_bottom = end[1] 
        end = max(end)

        if start == 0 and end == len(dna.seq):
            new_dna = copy.copy(dna)
            new_dna._topology = "linear"
            return new_dna

        if dna.topology == "circular":
            start = len(dna.seq) + start if start < 0 else start
            start = start - len(dna.seq) if start > len(dna.seq) else start
            end   = end - len(dna.seq) if end > len(dna.seq) else end
        
        if (start >= end or (start_top == end_top and start_bottom == end_bottom)) and dna.topology == "circular":
            subdna1 = extract(dna, [start, start], [len(dna.seq), len(dna.seq)])
            if start == end and start == 0:
                subdna = subdna1
            else:
                subdna2 = extract(dna, [0,0], [end,end])
                subdna  = joindna(subdna1, subdna2, __direct=0)
        else:
            if start > end and dna.topology == "linear":
                raise ValueError("'end' position must be larger than 'start' position.")
            feats = []
            new_features = []
            
            #Linearize feature (Split feature covering zero position) 
            for feat in dna.dnafeatures:
                strand = feat.strand
                s = feat.start
                e = feat.end
                if "_original" not in feat.__dict__:
                    feat._original = dna.printsequence(s, e, feat.location.strand if feat.location.strand !=0 else 1) 
                if s > e:
                    if len(feat.location.parts) == 1:
                        length = len(dna.seq) - s + e
                        locations = [FeatureLocation(s,len(dna.seq)),FeatureLocation(0,e)]
                        if strand == -1:
                            locations.reverse()
                        feat.location = CompoundLocation(locations)
                        feat.location.strand = strand

                    strand = feat.strand
                    if len(feat.location.parts) == 2:
                        feat1 = copy.deepcopy(feat)
                        feat1.location = feat.location.parts[0]
                        feat1.location.strand = feat.location.strand
                        feat2 = copy.deepcopy(feat)
                        feat2.location = feat.location.parts[1]
                        feat2.location.strand = feat.location.strand

                    else:
                        feat1 = copy.deepcopy(feat)
                        new_locations = []
                        for part in feat1.location.parts:
                            if part.start.position > part.end.postion:
                                new_locations.append(FeatureLocation(part.start.position, len(dna.seq)))
                                break
                            else:
                                new_locations.append(part)
                        if strand == -1:
                            new_locations.reverse()
                        feat1.location = CompoundLocation(new_locations)
                        feat1.location.strand = strand
                        flag  = 0
                        feat2 = copy.deepcopy(feat)
                        new_locations = []
                        for part in feat1.location.parts:
                            if part.start.position > part.end.postion:
                                new_locations.append(FeatureLocation(0, part.end.position))
                                flag = 1

                            if flag == 1:
                                new_locations.append(part)

                        if strand == -1:
                            new_locations.reverse()
                        feat2.location = CompoundLocation(new_locations)
                        feat2.location.strnad = strand

                    if "broken_feature" not in feat1.qualifiers:
                        label = feat1._id
                        if feat1.feature_type == "source":
                            original_seq = "-"
                        else:
                            original_seq = feat1.original
                        
                        if feat1.feature_type == "CDS" and "translation" in feat1.qualifiers:
                            del feat1.qualifiers["translation"]

                        label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat1.original), original_seq, s, e))
                        if strand >= 0:
                            feat1.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, len(dna.seq)-s)]
                        else:
                            feat1.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(dna.seq)-s, 1)]

                    else:
                        note = feat.qualifiers["broken_feature"]
                        note = note[0] if type(note) is list else note 
                        if strand >= 0:
                            label  = ":".join(note.split(":")[:-1])  
                            pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s, pos_s + len(dna.seq)-s)
                        else:
                            label  = ":".join(note.split(":")[:-1])
                            pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s, pos_s - (len(dna.seq)-s))
                        feat1.qualifiers["broken_feature"] = [note]

                    if "broken_feature" not in feat2.qualifiers:
                        label = feat2._id
                        if feat2.feature_type == "source":
                            original_seq = "-"
                        else:
                            original_seq = feat2.original
                        
                        if feat2.feature_type == "CDS" and "translation" in feat2.qualifiers:
                            del feat2.qualifiers["translation"]

                        label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat2.original), original_seq, s, e))
                        if strand >= 0:
                            feat2.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(dna.seq)-s+1, len(dna.seq)-s+e)]
                        else:
                            feat2.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(dna.seq)-s+e, len(dna.seq)-s+1)]

                    else:
                        note   = feat.qualifiers["broken_feature"][0]
                        if strand >= 0:
                            label  = ":".join(note.split(":")[:-1])
                            length = int(note.split(":")[-4])
                            pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s + len(dna.seq)-s, pos_e)
                        else:
                            label  = ":".join(note.split(":")[:-1])
                            length = int(note.split(":")[-4]) 
                            pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s - (len(dna.seq)-s), pos_e)
                        feat2.qualifiers["broken_feature"] = [note]
                    new_features.append(feat.__class__(feature=feat1))
                    new_features.append(feat.__class__(feature=feat2))
                
                else:
                    #print(feat, start, end) 
                    new_features.append(feat.__class__(feature=feat))
            
            #Cropping
            for feat in new_features:
                strand = feat.strand
                s = feat.start
                e = feat.end
                if "_original" not in feat.__dict__:
                    feat._original = dna.printsequence(s, e, feat.location.strand if feat.location.strand !=0 else 1) 
                feat = copy.deepcopy(feat)
                if len(feat.location.parts) == 1 and s <= e:
                    if e > start and s < end:
                        if s - start < 0:
                            feat.location.parts[0]._start = ExactPosition(0)
                            if "broken_feature" not in feat.qualifiers:
                                label = feat._id
                                if feat.feature_type == "source" or len(feat.original) > 10000:
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), original_seq, s, e))
                                if strand >= 0:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, abs(s-start)+1, e-s)] 
                                else:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(feat.original) - abs(s-start), 1)] 
                            else:
                                note = feat.qualifiers["broken_feature"][0]
                                if strand >= 0:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ","")) + abs(s-start) 
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                else:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ","")) - abs(s-start) 
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                feat.qualifiers["broken_feature"] = [note]
                        else:
                            feat.location.parts[0]._start = ExactPosition(s - start) 
                    
                        feat.location.parts[-1]._end = ExactPosition(e - start)  
                        if feat.location.parts[-1]._end > end-start:
                            feat.location.parts[-1]._end = ExactPosition(end - start)
                            if "broken_feature" not in feat.qualifiers: 
                                label = feat._id
                                if feat.feature_type == "source" or len(feat.original) > 10000:
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), original_seq, s, e))
                                if strand >= 0: 
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, end-s)]
                                else:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(feat.original), len(feat.original)-(end-s)+1)]
                            else:
                                s = int(feat.location.parts[0].start.position)
                                note = feat.qualifiers["broken_feature"][0]
                                if strand >= 0:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s + (end-start-s)-1)
                                else:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s - (end-start-s)+1)
                                feat.qualifiers["broken_feature"] = [note]
                       
                        feat.location.strand = strand
                        feats.append(feat.__class__(feature=feat))
                
                else:
                    length = e-s
                    locations = []
                    sflag = 0 
                    eflag = 0
                    for apart in feat.location.parts:
                        s = apart.start.position 
                        e = apart.end.position
                        if e > start and s <= end:
                            _start = ExactPosition(s)
                            if s - start <= 0:
                                sflag = 1
                            _end = ExactPosition(e) 
                            if _end > end-start:
                                eflag = 1
                            locations.append([_start,_end,feat.location.strand])
                    
                    if len(locations) > 0:
                        s = int(locations[0][0])
                        e = int(locations[-1][1])
                        if s - start < 0 and sflag == 1:
                            locations[0][0] = ExactPosition(0)
                            if "broken_feature" not in feat.qualifiers:
                                label = feat._id
                                if feat.feature_type == "source" or len(feat.original) > 10000:
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), original_seq, s, e))
                                if strand >= 0:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, abs(s-start)+1, e-s)]
                                else:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, e-s,  abs(s-start)+1)] 
                            else:
                                note   = feat.qualifiers["broken_feature"][0]
                                if strand >= 0:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ","")) + abs(s-start) 
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                else:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ","")) - abs(s-start) 
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                feat.qualifiers["broken_feature"] = [note]
                        else:
                            locations[0][0] = ExactPosition(s - start)
                        
                        if e > end-start and eflag == 1:
                            locations[-1][1] = ExactPosition(end-start)
                            if "broken_feature" not in feat.qualifiers:
                                label = feat._id 
                                if feat.feature_type == "source" or len(feat.original) > 10000:
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), original_seq, s, e))
                                feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, end-s)]
                            else:
                                s      = int(locations[0][0])
                                note   = feat.qualifiers["broken_feature"][0]
                                if strand  >= 0:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s + (end-start-s)-1)
                                else:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s - (end-start-s)+1)
                                feat.qualifiers["broken_feature"] = [note]
                        else:
                            locations[-1][1] = ExactPosition(e - start)
                    
                        if len(locations) == 1:
                            feat.location = FeatureLocation(*locations[0])
                        else:
                            for l in range(len(locations)):
                                if l == 0:
                                    locations[l][1] = locations[l][1] - start
                                elif l == len(locations) - 1:
                                    locations[l][0] = locations[l][0] - start
                                else:
                                    locations[l][0] = locations[l][0] - start
                                    locations[l][1] = locations[l][1] - start
                            locations = [FeatureLocation(*loc) for loc in locations] 
                            if strand == -1:
                                locations.reverse()
                            feat.location = CompoundLocation(locations)
                        feats.append(feat.__class__(feature=feat))
        
            feats.sort(key=lambda x:(x.location.parts[0].start.position, x.location.parts[-1].end.position))
            subdna = dna.__class__(seq=str(dna.seq[start:end]), _direct=0)
            subdna._history_feature = copy.deepcopy(dna._history_feature) 
            subdna._dnafeatures = feats

            subdna._features_dict = dict(list(map(lambda x:(x._id, x), subdna.dnafeatures)))
            subdna._topology = "linear"
            
            if start < len(dna._left_end) and dna.topology == "linear":
                subdna._left_end          = dna._left_end
                subdna._left_end_top      = dna._left_end_top
                subdna._left_end_bottom   = dna._left_end_bottom
            else:
                subdna._left_end          = subdna.seq[0:20] 
                subdna._left_end_top      = 1
                subdna._left_end_bottom   = 1
            
            if len(dna.seq) - end < len(dna._right_end) and dna.topology == "linear":
                subdna._right_end          = dna._right_end
                subdna._right_end_top      = dna._right_end_top
                subdna._right_end_bottom   = dna._right_end_bottom
            else:
                subdna._right_end          = subdna.seq[-20:]
                subdna._right_end_top      = 1
                subdna._right_end_bottom   = 1

            subdna.record.annotations["topology"] = subdna.topology
            subdna.record.features = subdna.dnafeatures
        
        #Generate sticky end
        if start_top != start_bottom or end_top != end_bottom:
            start_dif = start_top - start_bottom
            if start_dif > 0:
                left = "-" * start_dif + "/" +  "*" * start_dif
            elif start_dif < 0:
                left = "*" * abs(start_dif) + "/" +  "-" * abs(start_dif) 
            else:
                left = ""
            
            end_dif = end_top - end_bottom
            if end_dif > 0:
                right = "*" * end_dif + "/" + "-" * end_dif  
            elif end_dif < 0:
                right = "-" * abs(end_dif) + "/" + "*" * abs(end_dif) 
            else:
                right = ""
            subdna = modifyends(subdna, left, right, __direct=0)
        else:
            pass  
        
        for dnafeature in subdna.dnafeatures:
            dnafeature.subject = subdna 
        
        if start >= end:
            subdna._positions = dna._positions[start:] + dna._positions[:end]
        else:
            subdna._positions = dna._positions[start:end] 

        return subdna 
    
    dnas = [] 
    new_positions = [] 
    for pos in positions:
        if type(pos) is str:
            pos = tuple(map(int,pos.split("/")))
            spos, epos = pos
            spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
            epos = epos + len(dna.seq) if epos < 0 else epos
            new_positions.append((spos,epos))  
        
        elif type(pos) is int or ("__dict__" in dir(pos) and "_qint" in pos.__dict__):
            pos = (pos, pos)  
            spos, epos = pos
            spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
            epos = epos + len(dna.seq) if epos < 0 else epos
            new_positions.append((spos,epos))  
       
        elif type(pos) is SeqFeature or ("__dict__" in dir(pos) and "_dnafeature" in pos.__dict__):
            strand = pos.location.strand
            if "cutsite" not in pos.qualifiers:
                raise ValueError("DNAfeature object should hold 'qualifiers:cutsite' attribute.")
            
            if pos._digestion_topl == "null":
                _, _, pos._digestion_topl, pos._digestion_topr, pos._digestion_bottoml, pos._digestion_bottomr = compile_cutsite(pos.qualifiers["cutsite"][0])
            
            if strand != -1:
                if pos._digestion_topl != "null":
                    spos = pos.start - pos._digestion_topl
                    epos = pos.start - pos._digestion_bottoml 
                    spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
                    epos = epos + len(dna.seq) if epos < 0 else epos
                    new_positions.append((spos,epos))
                
                elif pos._digestion_topr != "null": 
                    spos = pos.end + pos._digestion_topr
                    epos = pos.end + pos._digestion_bottomr
                    spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
                    epos = epos + len(dna.seq) if epos < 0 else epos
                    new_positions.append((spos,epos))
            else:
                if pos._digestion_topr != "null":
                    spos = pos.start - pos._digestion_bottomr
                    epos = pos.start - pos._digestion_topr
                    spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
                    epos = epos + len(dna.seq) if epos < 0 else epos
                    new_positions.append((spos,epos))
                
                elif pos._digestion_topl != "null": 
                    spos = pos.end + pos._digestion_bottoml
                    epos = pos.end + pos._digestion_topl
                    spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
                    epos = epos + len(dna.seq) if epos < 0 else epos
                    new_positions.append((spos,epos))
    
    tmp_positions    = new_positions[:]
    tmp_positions.sort() 
    top_positions    = list(list(zip(*tmp_positions))[0])
    bottom_positions = list(list(zip(*tmp_positions))[1])
    for b in range(len(bottom_positions)-1):
        if bottom_positions[b] <= bottom_positions[b+1]:
            pass
        else:
            raise ValueError("Invalid cut pattern.")

    new_positions_original = new_positions[:] 
    new_positions_original = ["{}/{}".format(*posset) for posset in new_positions_original]
    if crop == True:
        crop_positions = (new_positions[0], new_positions[1])
    
    if dna.topology == "linear":
        if (0,0) not in new_positions:
            new_positions = [(0,0)]  + new_positions 
        else:
            pass
            #new_positions = new_positions     
        if (len(dna.seq),len(dna.seq)) not in new_positions:
            new_positions = new_positions + [(len(dna.seq), len(dna.seq))] 
        new_positions = list(new_positions) 
        new_positions.sort() 
    
    elif dna.topology == "circular":
        new_positions = list(new_positions) 
        tmp_positions = new_positions[:]
        new_positions.sort() 
        for pindex, pos in enumerate(new_positions):
            if pos == tmp_positions[0]:
                new_positions = new_positions[pindex:] + new_positions[:pindex]
                break 

    if dna.topology == "linear":
        if crop == True:
            dnas.append(extract(dna, crop_positions[0], crop_positions[1], project=project))
        else:
            for i, pos in enumerate(new_positions[0:-1]):
                dnas.append(extract(dna, pos, new_positions[i+1], project=project))
    
    elif dna.topology == "circular":
        if crop == True: 
            dnas.append(extract(dna, crop_positions[0], crop_positions[1], project=project))
        else:
            for i, pos in enumerate(new_positions[0:-1]):
                dnas.append(extract(dna, pos, new_positions[i+1], project=project))
            dnas.append(extract(dna, new_positions[-1], new_positions[0], project=project)) 

    if project is None:
        for subdna in dnas:
            subdna._unique_id = dna._unique_id
    else:
        for subdna in dnas:
            subdna._unique_id = project
    
    if __direct == 1:
        products = []
        dna_keys = list(dnas[0].__class__.dna_dict.keys())
        for i in range(len(dnas)):
            dnas[i]._product_id = dnas[i]._unique_id if product is None else product 
            products.append("QUEEN.dna_dict['{}']".format(dnas[i]._product_id))

        args = [] 
        history_features = [dnas[0]._history_feature] 
        for pos in positions:
            if "__dict__" in dir(pos) and "_dnafeature" in pos.__dict__:
                qkey = pos._qkey
                for qindex, qfeat in enumerate(dnas[0].__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == pos._second_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}]".format(qkey, qindex))
                history_features.append(pos.subject._history_feature) 

            elif "__dict__" in dir(pos) and "_qint" in pos.__dict__:
                qkey = pos.qkey
                for qindex, qfeat in enumerate(dnas[0].__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == pos.parental_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}].{}".format(pos.qkey, qindex, pos.name))
                history_features.append(pos.parent.subject._history_feature) 

            else:
                if type(pos) is int:
                    args.append(str(pos))
                else:
                    args.append("'" + str(pos) + "'")
        
        if crop == True:
            fcrop = ", crop=True"
        else:
            fcrop = "" 
        
        project             = "" #if project is None else ", project='" + project + "'"
        fproduct            = "" if product is None else ", product='" + product + "'"
        process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
        process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
        
        if len(products) > 1:
            building_history = "{} = cutdna(QUEEN.dna_dict['{}'], {}{}{}{}{}{})".format(", ".join(products), dna._product_id, ", ".join(args), fcrop, project, fproduct, process_name, process_description) 
        else:
            building_history = "{}, = cutdna(QUEEN.dna_dict['{}'], {}{}{}{}{}{})".format(", ".join(products), dna._product_id, ", ".join(args), fcrop, project, fproduct, process_name, process_description)
        
        for subdna in dnas:
            history_feature = _combine_history(subdna, history_features)
            subdna._history_feature = history_feature
            process_id, original_ids = make_processid(subdna, building_history, process_id, original_ids)
            add_history(subdna, [building_history, "positions: {}".format(",".join(list(map(str, new_positions_original)))) + "; num_products: {}".format(len(dnas)), ",".join([process_id] + original_ids)], _sourcefile=_sourcefile)
    else:
        for subdna in dnas:
            subdna.__dict__["_product_id"] = dna._product_id if "_product_id" in dna.__dict__ else dna._unique_id
    
    if product is None:
        pass 
    else:
        product = product.replace(" ","") 
        if "," in product:
            for name, subdna in zip(product.split(","), dnas):
                dnas[0].__class__._namespace[name] = subdna
        else:   
            dnas[0].__class__._namespace[product] = dnas

    if crop == True:
        return dnas[0], crop_positions 
    else:
        return dnas

def cropdna(dna, start=0, end=None, project=None, product=None, process_description=None, process_name=None, pd=None, pn=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1):
    project = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description

    #if start >= end and dna.topology == "linear":
    #    start, end = end, start

    if end is None or end == 0:
        end = len(dna.seq) 
    
    subdna, crop_positions = cutdna(dna, start, end, project=project, crop=True, __direct=0)  
    if __direct == 1:   
        args = []
        new_positions = [] 
        history_features = [subdna._history_feature] 
        for new_pos, pos in zip(crop_positions, (start, end)):
            new_pos = "/".join(list(map(str,new_pos)))
            new_pos = "'" + new_pos + "'"
            if "__dict__" in dir(pos) and "_dnafeature" in pos.__dict__:
                qkey = pos._qkey
                for qindex, qfeat in enumerate(dna.__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == pos._second_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}]".format(qkey, qindex))
                history_features.append(pos.subject._history_feature) 
            
            elif "__dict__" in dir(pos) and "_qint" in pos.__dict__:
                qkey = pos.qkey
                for qindex, qfeat in enumerate(dna.__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == pos.parental_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}].{}".format(pos.qkey, qindex, pos.name))
                history_features.append(pos.parent.subject._history_feature) 
        
            else:
                args.append(new_pos)

            pos1, pos2 = new_pos.split("/")  
            if pos1[1:] == pos2[:-1]:
                new_positions.append(pos1[1:]) 
            else:
                new_positions.append(new_pos) 

        project             = "" #if project is None else ", project='" + project + "'"
        fproduct            = "" if product is None else ", product='" + product + "'"
        process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
        process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
        
        subdna._product_id = subdna._unique_id if product is None else product 
        building_history = "QUEEN.dna_dict['{}'] = cropdna(QUEEN.dna_dict['{}'], start={}, end={}{}{}{}{})".format(subdna._product_id, dna._product_id, args[0], args[1], project, fproduct, process_name, process_description)
        
        history_feature = _combine_history(subdna, history_features)
        subdna._history_feature = history_feature
        process_id, original_ids = make_processid(subdna, building_history, process_id, original_ids)
        add_history(subdna, [building_history, "start: {}; end: {}".format(*new_positions), ",".join([process_id] + original_ids)], _sourcefile=_sourcefile) 
    else:
        subdna.__dict__["_product_id"] = dna._product_id if "_product_id" in dna.__dict__ else dna._unique_id

    if product is None:
        pass 
    else:
        product.replace(" ","") 
        match = re.fullmatch("(.+)\[(.+)\]", product)
        if match:
            if match.group(2).isdecimal() == True:
                subdna.__class__._namespace[match.group(1)][int(match.group(2))] = subdna
            else:
                subdna.__class__._namespace[match.group(1)][match.group(2)] = subdna
        else:    
            subdna.__class__._namespace[product] = subdna
    return subdna

def joindna(*dnas, topology="linear", project=None, product=None, process_name=None, process_description=None, pn=None, pd=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1):
    project             = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description

    new_dnas = [] 
    for i, dna in enumerate(dnas):
        if dna.topology == "circular":
            if i == 0:
                order = "first"
            elif i == 1:
                order = "second"
            elif i == 2:
                order = "third"
            else:
                order = str(i) + "th" 
            raise ValueError("The {} QUEEN object topology is 'circular.' Circular QUEEN objects cannot be connected with others.".format(order)) 
        new_dnas.append(dna) 
    
    dnas = new_dnas
    #Extract history information
    history_features = [] 
    for dna in dnas:
        history_features.append(dna._history_feature)
    
    construct = copy.deepcopy(dnas[0])
    positions_list = [construct._positions] 
    if len(dnas) > 1:
        for dna in dnas[1:]:
            anneal = False
            feats  = dna.dnafeatures
            if (dna._left_end_top * construct._right_end_bottom == 1 or dna._left_end_bottom * construct._right_end_top == 1) and ((dna._left_end_top == -1 or dna._left_end_bottom == -1) or (construct._right_end_top == -1 or construct._right_end_bottom == -1)):
                if dna._ssdna == False:
                    if dna._left_end_top == 1:
                        sticky_end = dna._left_end 
                    else:
                        sticky_end = construct._right_end
                    ovhg    = dna._left_end
                    new_dna = cropdna(dna,len(ovhg),len(dna.seq),__direct=0) 
                    if (dna._left_end == construct._right_end):
                        pass
                        #print("The DNA objects were joined based on complementary sticky end of each fragment. The sticky end is '{}'".format(sticky_end)) 
                    else:
                        raise ValueError("The QUEEN_objects cannot be joined due to the end structure incompatibility.")
                        return False
                else:
                    anneal=True
                    ovelap_queen_dict = {} 
                    ovhg = dna._left_end 
                    for i in range(1, len(construct._right_end)):
                        if construct._right_end[-1*i:] == ovhg[:i]:
                            top = construct._right_end + (len(ovhg)-i) * "-"
                            bottom = (len(construct._right_end)-i) * "-" + ovhg.translate(str.maketrans("ATGC","TACG"))
                            ovelap_queen_dict[i] = construct.__class__(seq=top.upper() + "/" + bottom.upper(), _direct=0)
                        else:
                            pass 
                    items = list(ovelap_queen_dict.items())
                    items.sort() 
                    new_q = items[-1][1] 
            else:
                if (construct._right_end == "" and ((dna._left_end == "") or (dna._left_end == dna.seq))) or (construct._right_end_top == 1 and construct._right_end_bottom == 1 and dna._left_end_top == 1 and dna._left_end_bottom == 1):
                    new_dna = dna
                    ovhg = ""
                else:
                    raise ValueError("The QUEEN_objects cannot be joined due to the end structure incompatibility.")
                    return False
            
            feats  = _slide(feats, len(construct.seq) - len(ovhg))
            feats1 = [feat for feat in construct.dnafeatures if "broken_feature" in feat.qualifiers]
            feats2 = [feat for feat in feats if "broken_feature" in feat.qualifiers]
            if anneal == True:
                construct._seq = new_q._seq
                construct._right_end = new_q._right_end
                construct._right_end_top = new_q._right_end_top
                construct._right_end_bottom = new_q._right_end_bottom
                construct._left_end = new_q._left_end
                construct._left_end_top = new_q._left_end_top
                construct._left_end_bottom = new_q._left_end_bottom
                construct._topology  = "linear"
                construct._positions = new_q._positions 
                construct._ssdna = False
                ovhg = new_q._right_end
                positions_list.append(construct._positions)
            else:
                construct._seq = construct.seq + new_dna.seq 
                construct._right_end        = dna._right_end
                construct._right_end_top    = dna._right_end_top
                construct._right_end_bottom = dna._right_end_bottom
                construct._topology = "linear"
                ovhg = dna._right_end
                positions_list.append(new_dna._positions) 

            const_features = copy.copy(construct.dnafeatures) 
            #Restore a original feature from fragmented features
            if len(feats1) > 0 and len(feats2) > 0:
                for feat1 in feats1:
                    if feat1.location.strand == -1:
                        s1, e1 = feat1.location.parts[-1].start.position, feat1.location.parts[0].end.position
                    else:
                        s1, e1 = feat1.location.parts[0].start.position, feat1.location.parts[-1].end.position

                    for feat2 in feats2:
                        if feat2.location.strand == -1:
                            s2, e2 = feat2.location.parts[-1].start.position - (len(construct.seq) - len(ovhg)), feat2.location.parts[0].end.position - (len(construct.seq) - len(ovhg))
                        else:
                            s2, e2 = feat2.location.parts[0].start.position - (len(construct.seq) - len(ovhg)), feat2.location.parts[-1].end.position - (len(construct.seq) - len(ovhg))
                        if feat1.type == feat2.type:
                            flag = 0
                            for key in feat1.qualifiers:
                                if key == "broken_feature":
                                    pass 
                                elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                                    flag = 1
                                else:
                                    flag = 0
                                    break    
                            if flag == 1:
                                note1   = feat1.qualifiers["broken_feature"][0]
                                label1  = ":".join(note1.split(":")[:-1])
                                length1 = int(note1.split(":")[-4]) 
                                pos_s1  = int(note1.split(":")[-1].split("..")[0].replace(" ",""))
                                pos_e1  = int(note1.split(":")[-1].split("..")[1].replace(" ",""))

                                note2   = feat2.qualifiers["broken_feature"][0]
                                label2  = ":".join(note2.split(":")[:-1])
                                length2 = int(note2.split(":")[-4]) 
                                pos_s2  = int(note2.split(":")[-1].split("..")[0].replace(" ",""))
                                pos_e2  = int(note2.split(":")[-1].split("..")[1].replace(" ",""))
                                
                                #Join fragmented features
                                if length1 == length2 and "_original" in feat1.__dict__ and "_original" in feat2.__dict__ and feat1.original == feat2.original and feat1.location.strand == feat2.location.strand:
                                    note    = "{}:{}..{}".format(label1, pos_s1, pos_e2)
                                    new_seq = construct.seq[s1:e1] + dna.seq[s2:e2]
                                    feat1_index = const_features.index(feat1)
                                    new_feat    = copy.deepcopy(const_features[feat1_index]) 
                                    strand      = new_feat.location.strand
                                    if len(feat1.location.parts) == 1 and len(feat2.location.parts) == 1:
                                        new_feat.location = FeatureLocation(feat1.location.parts[0].start.position, feat2.location.parts[-1].end.position, feat1.strand)
                                        new_feat.location.strand = strand
                                    else:
                                        locations = feat1.location.parts[0:-1] + [FeatureLocation(feat1.location.parts[-1].start.position, feat2.location.parts[0].end.position, feat1.strand)] + feat2.location.parts[0:-1]
                                        if strand == -1:
                                            locations.reverse() 
                                        new_feat.location = CompoundLocation(locations) 
                                        new_feat.location.strand = strand 
                                   
                                    new_feat  = feat1.__class__(feature=new_feat, subject=construct)
                                    new_feat1 = feat1.__class__(feature=feat1, subject=construct)
                                    new_feat2 = feat1.__class__(feature=feat2, subject=construct) 
                                    s = new_feat.start 
                                    e = new_feat.end if new_feat.end <= len(construct.seq) else new_feat.end - len(construct.seq) 
                                    if construct.printsequence(s, e, new_feat.location.strand if new_feat.location.strand !=0 else 1) in new_feat.original:
                                        construct._dnafeatures[feat1_index] = feat1.__class__(feature=new_feat)
                                        construct._dnafeatures[feat1_index].qualifiers["broken_feature"] = [note]
                                        if feat2 in feats:
                                            del feats[feats.index(feat2)] 
            
            construct._dnafeatures = construct.dnafeatures + feats

        construct._dnafeatures.sort(key=lambda x:x.location.parts[0].start.position)
        for feat in construct.dnafeatures:
            if "broken_feature" in feat.qualifiers:
                note   = feat.qualifiers["broken_feature"][0]
                label  = ":".join(note.split(":")[:-1])
                length = int(note.split(":")[-4]) 
                pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["broken_feature"]
        if Alphabet:
            new_record = SeqRecord(Seq(str(construct.seq), Alphabet.DNAAlphabet()))
        else:
            new_record = SeqRecord(Seq(str(construct.seq)))

        new_record.features = construct.dnafeatures
        new_record.annotations["topology"] = topology
        construct.record = new_record     
        
        if topology == "circular":
            construct = _circularizedna(construct)
            #Adjust zero positions
            if __direct == 1: 
                zero_positions = [] 
                for d, positions in enumerate(positions_list):
                    if 0 in positions:
                        zero_positions.append((len(positions),d,positions.index(0)))
                if len(zero_positions) > 0:
                    zero_positions.sort() 
                    zero_positions.reverse() 
                    zero_position = 0
                    for dna in dnas[0:zero_positions[0][1]]:
                        zero_position += len(dna.seq)
                    zero_position += zero_positions[0][2]
                    construct = cutdna(construct, zero_position, __direct=0)[0]
                    construct = _circularizedna(construct) 
                construct._positions = tuple(range(len(construct.seq)))    
            else:
                construct._positions = tuple(range(len(construct.seq)))    

        else:
            zero_positions = [] 
            for d, positions in enumerate(positions_list):
                if 0 in positions:
                    zero_positions.append((len(positions),d,positions.index(0)))
            
            if len(zero_positions) > 0:
                zero_positions.sort() 
                zero_positions.reverse() 
                zero_origin = zero_positions[0][1]
                new_positions = [] 
                for d, positions in enumerate(positions_list): 
                    if d == zero_origin:
                        new_positions.extend(positions)
                    else:
                        new_positions.extend([-1] * len(positions))
                construct._positions = tuple(new_positions) 
            else:
                construct._positions = tuple(range(len(construct.seq)))
        
        construct._setfeatureid() #Update feature ID
    else:
        construct = _circularizedna(dnas[0])
        construct._positions = construct._positions[0:len(construct.seq)]

    if project is None:
        construct._unique_id = dnas[0]._unique_id
    else:
        construct._unique_id = project

    #Recover fragmented features if complete sequence is in the construct.
    new_features = [] 
    remove_features = [] 
    for feat in construct.dnafeatures:
        if "broken_feature" in feat.qualifiers:
            note       = feat.qualifiers["broken_feature"][0]
            label      = ":".join(note.split(":")[:-1])
            poss, pose = list(map(int,note.split(":")[-1].split("..")))
            length     = int(note.split(":")[-4])  
            #print(feat.start, poss, len(construct.seq),  ) 
            if feat.location.strand != -1:
                sfeat = feat.start-(poss-1) 
                sfeat = sfeat if sfeat > 0 else len(construct.seq) + sfeat
                efeat = feat.end+(length-pose)
            else:
                sfeat = feat.start-(length-pose) 
                sfeat = sfeat if sfeat > 0 else len(construct.seq) + sfeat
                efeat = feat.end+(poss-1)    
            
            #print(feat.qualifiers["label"], len(feat.original), len(note.split(":")[-3]), sfeat, efeat, len(construct.printsequence(sfeat, efeat, strand=feat.location.strand))) 
            if note.split(":")[-3] == construct.printsequence(sfeat, efeat, strand=feat.location.strand):
                if sfeat < efeat:
                    location = FeatureLocation(sfeat, efeat, feat.location.strand) 
                else:
                    location = CompoundLocation([FeatureLocation(sfeat, len(construct.seq)), FeatureLocation(0, efeat, feat.location.strand)])  
                newfeat = feat.__class__(location=location, subject=construct)
                newfeat.type = feat.type
                newfeat.qualifiers = feat.qualifiers
                del newfeat.qualifiers["broken_feature"]
                newfeat._id = feat.feature_id
                new_features.append(newfeat)
                remove_features.append(feat)

    for feat in remove_features:
        del construct._dnafeatures[construct.dnafeatures.index(feat)]
    
    for feat in new_features:
        construct._dnafeatures.append(feat) 

    construct._features_dict = dict(list(map(lambda x:(x._id, x), construct.dnafeatures)))
    construct.record.feartures = construct.dnafeatures
    
    if __direct == 1:
        project             = "" #if project is None else ", project='" + project + "'"
        fproduct            = "" if product is None else ", product='" + product + "'"
        process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
        process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
        
        construct._product_id = construct._unique_id if product is None else product 
        construct.record.id   = construct.project
        dna_elements = "[" + ", ".join(["QUEEN.dna_dict['{}']".format(dna._product_id) for dna in dnas]) + "]"
        building_history = "QUEEN.dna_dict['{}'] = joindna(*{}, topology='{}'{}{}{}{})".format(construct._product_id, dna_elements, topology, project, fproduct, process_name, process_description) 
        history_feature = _combine_history(construct, history_features)         
        construct._history_feature = history_feature 
        process_id, original_ids = make_processid(construct, building_history, process_id, original_ids)
        add_history(construct, [building_history, "topology: {}".format(topology), ",".join([process_id] + original_ids)], _sourcefile=_sourcefile) 
    else:
        construct.__dict__["_product_id"] = dnas[0]._product_id if "_product_id" in dnas[0].__dict__ else dnas[0]._unique_id

    for dnafeature in construct.dnafeatures:
        dnafeature.subject = construct
    
    if product is None:
        pass 
    else:
        product = product.replace(" ","") 
        match   = re.fullmatch("(.+)\[(.+)\]", product) 
        if match:
            if match.group(2).isdecimal() == True:
                construct.__class__._namespace[match.group(1)][int(match.group(2))] = construct
            else:
                construct.__class__._namespace[match.group(1)][match.group(2)] = construct
        else:    
            construct.__class__._namespace[product] = construct
    return construct

def modifyends(dna, left="", right="", add=0, add_right=0, add_left=0, project=None, product=None, process_name=None, process_description=None, pn=None, pd=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1):
    project             = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description

    if dna.topology == "circular":
        raise ValueError("End sequence structures cannot be modified. The topology of the given QUEEN_object is circular.") 
    else:
        pass
    #if process_description is None:
    #    process_description = dna.__class__.process_description
    #else:
    #    dna.__class__.process_description = process_description 
    
    def parse(seq,count=0):
        l_bracket = re.compile("\(")
        r_bracket = re.compile("\)")
        l_brace   = re.compile("\{")
        r_brace   = re.compile("\}")
        if set(str(seq)) <= set("0123456789ATGCRYKMSWBDHVNatgcnrykmswbdhv{}()/-*"):
            lbk_list = [l.start() for l in re.finditer(l_bracket,seq)]
            lbc_list = [l.start() for l in re.finditer(l_brace,seq)]
            rbc_list = [l.start() for l in re.finditer(r_brace,seq)]
            rbk_list = []
            for bk in lbk_list:
                n = 0
                p = bk
                for c in seq[bk:]:
                    if c == "(":
                        n += 1
                    elif c == ")":
                        n -= 1
                    if n == 0:
                        rbk_list.append(p)
                        break
                    p += 1
                    
            if len(lbk_list) == len(rbk_list):
                bk_set  = list(zip(lbk_list,rbk_list))
                bc_set  = list(zip(lbc_list,rbc_list))
                if len(bk_set) > 0:
                    for bks in bk_set:
                        new_seq = seq[bks[0]+1:bks[1]] 
                        if "(" not in new_seq:                   
                            num  = 0
                            flag = 0
                            sub_lbc_list = [l.start() for l in re.finditer(l_brace,new_seq)]
                            sub_rbc_list = [l.start() for l in re.finditer(r_brace,new_seq)]
                            sub_bc_set  = list(zip(sub_lbc_list, sub_rbc_list))
                            new_new_seq = new_seq
                            sub_bc_set.sort()
                            for bcs in sub_bc_set:
                                try:
                                    num = int(new_seq[bcs[0]+1:bcs[1]])
                                except:
                                    return False
                                if num > 0:
                                    new_new_seq = new_new_seq[:new_new_seq.find("{")-1] + new_seq[bcs[0]-1] * num + new_seq[bcs[1]+1:] 
                                else:
                                    pass
                            new_seq = new_new_seq
                            num = 0
                            for bcs in bc_set:
                                if bcs[0] == bks[1]+1:
                                    try:
                                        num = int(seq[bcs[0]+1:bcs[1]])
                                    except:
                                        return False
                                    break
                                else:
                                    pass
                            if num > 0:
                                new_seq = new_seq * num
                            else:
                                pass
                            break
                        else:
                            pass
                    new_seq = seq[:bks[0]] + new_seq + seq[bcs[1]+1:]
                
                else:
                    new_seq = seq
                    bc_set.sort()
                    for bcs in bc_set:
                        try:
                            num = int(seq[bcs[0]+1:bcs[1]])
                        except:
                            return False
                        if num > 0:
                            new_seq = new_seq[:new_seq.find("{")-1] + seq[bcs[0]-1] * num + seq[bcs[1]+1:] 
                        else:
                            pass
                
                if set(str(new_seq)) <= set("ATGCRYKMSWBDHVNatgcnrykmswbdhv/-*"):
                    pass
                else:
                    new_seq = parse(new_seq,count+1)
                return new_seq
            else:
                return False
        else:
            return False
    
    def check_endseq(top,bottom):
        #Check if the top strand seqeunce is complement with the bottom strand.  
        new_top    = ""
        new_bottom = ""
        for t,b in zip(top,bottom):
            if t != b.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) and (t != "-" and b != "-"):
                return False, False
            new_top    += t
            new_bottom += b
        return new_top, new_bottom
    
    if left == "":
        left = "*/*"
    
    if right == "":
        right = "*/*"
    
    left_origin, right_origin = left, right
    left, right = parse(left.upper()), parse(right.upper())
    left, rihgt = str(left), str(right) 
    pattern1, pattern2, patternl1, patternl2, patternr1, patternr2 = "[ATGCRYKMSWBDHVN*-]*/?[ATGCRYKMSWBDHVN*-]*", "[ATGCRYKMSWBDHVN*]+-+[ATGCRYKMSWBDHVN*]+", "^[ATGCRYKMSWBDHVN*]+-+/", "/[ATGCRYKMSWBDHVN*]+-+$", "^-+[ATGCRYKMSWBDHVN*]+/", "/-+[ATGCRYKMSWBDHVN*]+$" 
    pattern1  = re.compile(pattern1) 
    pattern2  = re.compile(pattern2) 
    patternl1 = re.compile(patternl1)
    patternl2 = re.compile(patternl2)
    patternr1 = re.compile(patternr1)
    patternr2 = re.compile(patternr2)
    left_end, right_end  = left, right
    
    if pattern1.fullmatch(left_end) != None and pattern2.search(left_end) is None and patternl1.search(left_end) is None and patternl2.search(left_end) is None:
        pass 
    else:
        raise TypeError("Please sepcify a proper sequence pattern for the 'left' argument") 
    
    if pattern1.fullmatch(right_end) != None and pattern2.search(right_end) is None and patternr1.search(right_end) is None and patternr2.search(right_end) is None:
        pass
    else:
        raise TypeError("Please sepcify a proper sequence pattern for the 'right' argument") 
    
    if "/" in left_end:
        left_end_top, left_end_bottom = left_end.split("/")
    else:
        add = 1
        add_left = 1
        left_end_top, left_end_bottom = left_end, left_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
        left_end = left_end_top + "/" + left_end_bottom

    if len(left_end_top) != len(left_end_bottom):
        raise ValueError("Please specify a proper end sequence structure for the 'left' argument.")
        
    elif "-" in left_end or "*" in left_end:
        left_end_top, left_end_bottom = check_endseq(left_end_top, left_end_bottom)
        if left_end_top != False:
            left_end_length = len(left_end_top)
            if "*" in left_end_top or "*" in left_end_bottom:
                if set(left_end_top) <= set(["*","-"]) and set(left_end_bottom) <= set(["*","-"]):
                    left_end_top    = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[:left_end_length], left_end_top)])
                    left_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[:left_end_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")), left_end_bottom)])
                else:
                    raise ValueError("'*' cannot be used wih IUPAC nucleotide codes")
            
                left_end_top, left_end_bottom = left_end_top.replace("-",""), left_end_bottom.replace("-","") 
                if len(left_end_top) < len(left_end_bottom):
                    left_length     = len(left_end_bottom)
                    left_end        = left_end_bottom[0:len(left_end_bottom)-1*len(left_end_top)].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) 
                    left_end_top    = -1
                    left_end_bottom = 1

                elif len(left_end_top) > len(left_end_bottom):
                    left_length     = len(left_end_top)
                    left_end        = left_end_top[0:len(left_end_top)-1*len(left_end_bottom)]
                    left_end_top    = 1
                    left_end_bottom = -1
                
                else:
                    left_length     = len(left_end_top)
                    left_end        = left_end_top[0:len(left_end_top)-1*len(left_end_bottom)]
                    left_end_top    = 1
                    left_end_bottom = 1
            else:
                add = 1
                add_left = 1
                left_end = left_end_top + "/" + left_end_bottom

        else:
            raise ValueError("Please specify a proper end sequence structure for the 'left' argument.")
    else:
        left_end_length = len(left_end_top) 
        left_length = len(left_end_top)
        left_end_top    = 1
        left_end_bottom = 1

    if "/" in right_end:
        right_end_top, right_end_bottom = right_end.split("/")
    else:
        add = 1
        add_right = 1 
        right_end_top, right_end_bottom = right_end, right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
        right_end = right_end_top + "/" + right_end_bottom
    
    if len(right_end_top) != len(right_end_bottom):
        raise ValueError("Please specify a proper end sequence structure for the 'right' argument.")
    
    elif "-" in right_end or "*" in right_end:
        right_end_top, right_end_bottom = check_endseq(right_end_top, right_end_bottom)
        if right_end_top != False:
            right_end_length = len(right_end_top)
            if "*" in right_end_top or "*" in right_end_bottom:
                if set(right_end_top) <= set(["*","-"]) and set(right_end_bottom) <= set(["*","-"]):
                    right_end_top    = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[-1*right_end_length:], right_end_top)])
                    right_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[-1*right_end_length:].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")), right_end_bottom)])
                else:
                    raise ValueError("'*' cannot be used wih characters for DNA sequence.")

                right_end_top, right_end_bottom = right_end_top.replace("-",""), right_end_bottom.replace("-","") 
                if len(right_end_top) < len(right_end_bottom):
                    right_length     = len(right_end_bottom)
                    right_end        = right_end_bottom[len(right_end_top):].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) 
                    right_end_top    = -1
                    right_end_bottom = 1

                elif len(right_end_top) > len(right_end_bottom):
                    right_length     = len(right_end_top)
                    right_end        = right_end_top[len(right_end_bottom):] 
                    right_end_top    = 1
                    right_end_bottom = -1
               
                else:
                    right_length     = len(right_end_top)
                    right_end        = right_end_top[len(right_end_bottom):] 
                    right_end_top    = 1
                    right_end_bottom = 1
            else:
                add = 1
                add_right = 1 
                right_end = right_end_top + "/" + right_end_bottom

        else:
            raise ValueError("Please specify a proper end sequence structure for the 'right' argument.")
    else:
        right_end_length = len(right_end_top) 
        right_length = len(right_end_top)  
        right_end_top    = 1
        right_end_bottom = 1
 
    if add == 1 or (left_end != dna.seq[left_end_length-left_length:left_end_length-left_length+len(left_end)] 
                or right_end != str(dna[len(dna.seq)-right_end_length + right_length - len(right_end):len(dna.seq)-right_end_length + right_length].seq)):
        
        if add_left == 1 and add_right == 1:
            new_dna = dna.__class__(seq=left_end.split("/")[0] + dna.seq + right_end.split("/")[0] + "/" 
                    + left_end.split("/")[1] + dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + right_end.split("/")[1], _direct=0) 
            new_dna._dnafeatures     = _slide(dna.dnafeatures, len(left_end.split("/")[0])) 
            new_dna.record           = copy.copy(dna.record) 
            new_dna.record.features  = new_dna.dnafeatures
            new_dna._positions =  (-1,) * len(left_end.split("/")[0]) + new_dna._positions + (-1,) * len(right_end.split("/")[0])
            
        
        elif add_right == 1:
            left_end  = dna.__class__(seq=left_end,  _direct=0) 
            right_end = dna.__class__(seq=right_end, _direct=0)    
            new_dna   = cropdna(dna, start=left_end_length-left_length, end=len(dna.seq), __direct=0) + right_end
            new_dna._right_end       = right_end._right_end
            new_dna._left_end        = left_end.seq
            new_dna._left_end_top    = left_end_top 
            new_dna._left_end_bottom = left_end_bottom
            new_dna._positions = new_dna._positions + (-1,) * len(right_end.seq)
        
        else:
            left_end  = dna.__class__(seq=left_end,  _direct=0) 
            right_end = dna.__class__(seq=right_end, _direct=0)    
            new_dna = left_end + cropdna(dna, start=0, end=len(dna.seq)-right_end_length+right_length, __direct=0)
            new_dna._left_end         = left_end._left_end
            new_dna._right_end        = right_end.seq
            new_dna._right_end_top    = right_end_top 
            new_dna._right_end_bottom = right_end_bottom
            new_dna._positions = (-1,) * len(left_end.seq) + new_dna._positions
        
        new_dna._history_feature = dna._history_feature 
        
        #Recover fragmented features if complete sequence is in the construct.
        new_features = [] 
        remove_features = [] 
        for feat in new_dna.dnafeatures:
            if "broken_feature" in feat.qualifiers:
                note       = feat.qualifiers["broken_feature"][0]
                label      = ":".join(note.split(":")[:-1])
                poss, pose = list(map(int,note.split(":")[-1].split("..")))
                length     = int(note.split(":")[-4])  
                #print(feat.start, poss, len(construct.seq),  ) 
                if feat.location.strand != -1:
                    sfeat = feat.start-(poss-1) 
                    sfeat = sfeat if sfeat > 0 else len(new_dna.seq) + sfeat
                    efeat = feat.end+(length-pose)
                else:
                    sfeat = feat.start-(length-pose) 
                    sfeat = sfeat if sfeat > 0 else len(new_dna.seq) + sfeat
                    efeat = feat.end+(poss-1)    
                
                if note.split(":")[-3] == new_dna.printsequence(sfeat, efeat, strand=feat.location.strand):
                    if sfeat < efeat:
                        location = FeatureLocation(sfeat, efeat, feat.location.strand) 
                    else:
                        location = CompoundLocation([FeatureLocation(sfeat, len(new_dna.seq)), FeatureLocation(0, efeat, feat.location.strand)])  
                    newfeat = feat.__class__(location=location, subject=new_dna)
                    newfeat.type = feat.type
                    newfeat.qualifiers = feat.qualifiers
                    del newfeat.qualifiers["broken_feature"]
                    newfeat._id = feat.feature_id
                    new_features.append(newfeat)
                    remove_features.append(feat)
        
        for feat in remove_features:
            del new_dna._dnafeatures[new_dna.dnafeatures.index(feat)] 
        
        for feat in new_features:
            new_dna._dnafeatures.append(feat) 
    else:
        new_dna   = copy.deepcopy(dna)  
        new_dna._left_end  = left_end
        new_dna._right_end = right_end
        new_dna._left_end_top     = left_end_top 
        new_dna._left_end_bottom  = left_end_bottom
        new_dna._right_end_top    = right_end_top
        new_dna._right_end_bottom = right_end_bottom
    
    if project is None:
        new_dna._unique_id = dna._unique_id 
    else:
        new_dna._unique_id = project

    if __direct == 1:
        args = [] 
        history_features = [new_dna._history_feature] 
        args.append("'{}'".format(new_dna._unique_id))
        args.append("'{}'".format(dna._unique_id)) 
        ends = [] 
        if type(left_origin) == new_dna.seq.__class__:
            if left_origin.parental_class == "DNAFeature":
                qkey = left_origin.qkey
                for qindex, qfeat in enumerate(new_dna.__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == left_origin.parental_id:
                        break
                if type(left_origin.item)   == int:
                    args.append("QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , left_origin.item))
                elif type(left_origin.item) == slice:
                    sl_start = left_origin.item.start
                    sl_stop  = left_origin.item.stop 
                    sl_step  = left_origin.item.step
                    sl_start = "" if sl_start is None else sl_start
                    sl_stop  = "" if sl_stop is None else sl_stop
                    if sl_step == 1 or sl_step == None:
                        args.append("QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop))
                    else:
                        args.append("QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step))
                else:
                    args.append("QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex))
                history_features.append(left_origin.parent.subject._history_feature) 

            elif left_origin.parental_class == "QUEEN": 
                parental_id = left_origin.parental_id

                if left_origin.name != None: 
                    if "printsequence" in left_origin.name:
                        if len(left_origin.name.split("_")) == 2: 
                            seqname = "QUEEN.dna_dict['{}'].printsequence(strand={})".format(parental_id, left_origin.name.split("_")[-1]) 
                        else:
                            seqname = "QUEEN.dna_dict['{}'].printsequence(start={}, end={}, strand={})".format(parental_id, *left_origin.name.split("_")[1:])
                    elif lefet_origin.name == "rcseq":
                        seqname = "QUEEN.dna_dict['{}'].rcseq".format(parental_id) 
                else:
                    seqname = "QUEEN.dna_dict['{}'].seq".format(parental_id)

                if type(left_origin.item)   == int:
                    args.append("{}.seq[{}]".format(seqname, left_origin.item))
                elif type(left_origin.item) == slice:
                    sl_start = left_origin.item.start
                    sl_stop  = left_origin.item.stop 
                    sl_step  = left_origin.item.step
                    sl_start = "" if sl_start is None else sl_start
                    sl_stop  = "" if sl_stop is None else sl_stop
                    if sl_step == 1 or sl_step == None:
                        args.append("{}[{}:{}]".format(seqname, sl_start, sl_stop))
                    else:
                        args.append("{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step))
                else:
                    args.append("{}".format(seqname))
                history_features.append(left_origin.parent._history_feature) 
            
            elif left_origin.parental_class == "Cutsite":
                if left_origin.parent.name not in cs.defaultkeys:
                    cs.new_cutsites.append((left_origin.parent.name, left_origin.parent.cutsite)) 
                args.append("cs.lib['{}'].{}".format(left_origin.parent.name, left_origin.name)) 
            else:
                args.append("'{}'".format(left_origin)) 
            ends.append(left) 
        else:
            ends.append(left_origin) 
            args.append("'{}'".format(left_origin)) 
        
        if type(right_origin) == new_dna.seq.__class__:
            if right_origin.parental_class == "DNAFeature":
                qkey = right_origin.qkey
                for qindex, qfeat in enumerate(new_dna.__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == right_origin.parental_id:
                        break
                if type(right_origin.item)   == int:
                    args.append("QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , right_origin.item))
                elif type(right_origin.item) == slice:
                    sl_start = right_origin.item.start 
                    sl_stop  = right_origin.item.stop  
                    sl_step  = right_origin.item.step
                    sl_start = "" if sl_start is None else sl_start
                    sl_stop  = "" if sl_stop is None else sl_stop
                    if sl_step == 1 or sl_step == None:
                        args.append("QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop))
                    else:
                        args.append("QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step))
                else:
                    args.append("QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex))
                history_features.append(right_origin.parent.subject._history_feature) 

            elif right_origin.parental_class == "QUEEN": 
                parental_id = right_origin.parental_id
                
                if right_origin.name != None: 
                    if "printsequence" in right_origin.name:
                        if len(right_origin.name.split("_")) == 2: 
                            seqname = "QUEEN.dna_dict['{}'].printsequence(strand={})".format(parental_id, right_origin.name.split("_")[-1]) 
                        else:
                            seqname = "QUEEN.dna_dict['{}'].printsequence(start={}, end={}, strand={})".format(parental_id, *right_origin.name.split("_")[1:])
                    elif right_origin.name == "rcseq":
                        seqname = "QUEEN.dna_dict['{}'].rcseq".format(parental_id) 
                else:
                    seqname = "QUEEN.dna_dict['{}'].seq".format(parental_id)
                
                if type(right_origin.item)   == int:
                    args.append("{}[{}]".format(seqname, right_origin.item))
                elif type(right_origin.item) == slice:
                    sl_start = right_origin.item.start
                    sl_stop  = right_origin.item.stop 
                    sl_step  = right_origin.item.step
                    sl_start = "" if sl_start is None else sl_start
                    sl_stop  = "" if sl_stop is None else sl_stop
                    if sl_step == 1 or sl_step == None:
                        args.append("{}[{}:{}]".format(seqname, sl_start, sl_stop))
                    else:
                        args.append("{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step))
                else:
                    args.append("{}".format(seqname))
                history_features.append(right_origin.parent._history_feature) 
            
            elif right_origin.parental_class == "Cutsite":
                if right_origin.parent.name not in cs.defaultkeys:
                    cs.new_cutsites.add((right_origin.parent.name, right_origin.parent.cutsite)) 
                args.append("cs.lib['{}'].{}".format(right_origin.parent.name, right_origin.name)) 

            else:
                args.append("'{}'".format(right_origin)) 
            ends.append(right)
        else:
            ends.append(right_origin) 
            args.append("'{}'".format(right_origin)) 

        project             = "" #if project is None else ", project='" + project + "'"
        fproduct            = "" if product is None else ", product='" + product + "'"
        process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
        process_description = "" if process_description is None else ", process_description='" + process_description + "'" 

        new_dna._product_id = new_dna._unique_id if product is None else product 
        building_history    = "QUEEN.dna_dict['{}'] = modifyends(QUEEN.dna_dict['{}'], left={}, right={}{}{}{}{})".format(new_dna._product_id, dna._product_id, args[2], args[3], project, fproduct, process_name, process_description)  
        history_feature     = _combine_history(new_dna, history_features) 
        new_dna._history_feature = history_feature
        process_id, original_ids = make_processid(new_dna, building_history, process_id, original_ids)
        add_history(new_dna, [building_history, "left: {}; right: {}; leftobj: {}; rightobj: {}".format(*ends, args[2], args[3]), ",".join([process_id] + original_ids)], _sourcefile=_sourcefile) 
    else:
        new_dna.__dict__["_product_id"] = dna._product_id if "_product_id" in dna.__dict__ else dna._unique_id

    for dnafeature in new_dna.dnafeatures:
        dnafeature.subject = new_dna
    
    if product is None:
        pass 
    else:
        product = product.replace(" ","")
        match   = re.fullmatch("(.+)\[(.+)\]", product)
        if match:
            if match.group(2).isdecimal() == True:
                new_dna.__class__._namespace[match.group(1)][int(match.group(2))] = new_dna
            else:
                new_dna.__class__._namespace[match.group(1)][match.group(2)] = new_dna
        else:    
            new_dna.__class__._namespace[product] = new_dna
    return new_dna

def flipdna(dna, project=None, product=None, process_name=None, process_description=None, pn=None, pd=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1):
    project             = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description

    #if process_description is None:
    #    process_description = dna.__class__.process_description
    #else:
    #    dna.__class__.process_description = process_description 

    if type(dna) is str:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        return seq

    elif type(dna) is Seq:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 
    
    else:
        dna = copy.deepcopy(dna)
        seq = dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 
        for feat in dna.dnafeatures:
            strand = feat.location.strand
            for p in range(len(feat.location.parts)):
                s, e = feat.location.parts[p].start.position, feat.location.parts[p].end.position
                feat.location.parts[p]._start = ExactPosition(len(dna.seq) - e) 
                feat.location.parts[p]._end   = ExactPosition(len(dna.seq) - s) 
            
            if "original" in feat.__dict__:
                feat._original = feat.original.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]

            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
                if "broken_feature" in feat.qualifiers:
                    note   = feat.qualifiers["broken_feature"][0]
                    label  = ":".join(note.split(":")[:-1])
                    length = int(note.split(":")[-4]) 
                    pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                    pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                    note = "{}:{}..{}".format(label, pos_e, pos_s)
                    feat.qualifiers["broken_feature"] = [note]
            else:
                feat.location.strand = strand
            feats.append(feat.__class__(feature=feat,subject=seq))
    
    comp = dna.__class__(seq=seq, _direct=0) 
    feats.sort(key=lambda x:x.location.parts[0].start.position) 
    comp._dnafeatures = feats
    comp._history_feature = dna._history_feature
    comp._setfeatureid()
    comp._features_dict = dict(list(map(lambda x:(x._id, x), comp.dnafeatures)))
    comp.record.features = comp.dnafeatures
    comp._right_end, comp._left_end = dna._left_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], dna._right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
    comp._right_end_top, comp._left_end_bottom = dna._left_end_bottom, dna._right_end_top
    comp._right_end_bottom, comp._left_end_top = dna._left_end_top, dna._right_end_bottom
    
    if project is None:
        comp._unqiue_id = dna._unique_id
    else:
        comp._unique_id = project 

    if __direct == 1: 
        project             = "" #if project is None else ", project='" + project + "'"
        fproduct            = "" if product is None else ", product='" + product + "'"
        process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
        process_description = "" if process_description is None else ", process_description='" + process_description + "'" 

        comp._product_id = comp._unique_id if product is None else product 
        building_history = "QUEEN.dna_dict['{}'] = flipdna(QUEEN.dna_dict['{}']{}{}{}{})".format(comp._product_id, dna._product_id, project, fproduct, process_name, process_description) 
        process_id, original_ids = make_processid(comp, building_history, process_id, original_ids)
        add_history(comp, [building_history,"", ",".join([process_id] + original_ids)], _sourcefile=_sourcefile) 
    else:
        comp.__dict__["_product_id"] = dna._product_id if "_product_id" in dna.__dict__ else dna._unique_id

    comp._positions = dna._positions[::-1] 
    for dnafeature in comp.dnafeatures:
        dnafeature.subject = comp
    
    if product is None:
        pass 
    else:
        product = product.replace(" ","")
        match   = re.fullmatch("(.+)\[(.+)\]", product)
        if match:
            if match.group(2).isdecimal() == True:
                comp.__class__._namespace[match.group(1)][int(match.group(2))] = comp
            else:
                #print(match.group(1),match.group(2))
                comp.__class__._namespace[match.group(1)][match.group(2)] = comp
        else:    
            comp.__class__._namespace[product] = comp
    return comp

def get_matchlist_regex(dna, query, value=None, subject=None, s=None, e=None, strand=None):
    segment = dna.__class__(seq="", _direct=0) 
    match_list = [] 
    
    if type(query) == str:
        query_pattern = re.compile(query) 
    else:
        query_pattern = query

    if value is None:
        mode = "search"
    else:
        mode = "edit"
    
    match_iter = re.finditer(query, subject) 
    if strand == -1:
        pre_s = s 
        match_iter = list(match_iter)
        match_iter.reverse()
        match_set = set([]) 
        for match in match_iter:
            result  = {"start":None, "end":None, "strand":None, "match":None} 
            span    = match.span() 
            span    = (span[0] - e if span[0] > e else span[0], span[1] - e if span[1] > e else span[1])
            if span not in match_set:
                match_set.add(span) 
                span    = (e-span[1], e-span[0])
                result["start"]  = span[0]
                result["end"]    = span[1]
                result["strand"] = strand
                result["match"]  = match
                match_list.append(result)     
                if mode == "edit":
                    groups, literals = sre_parse.parse_template(value, query_pattern)
                    groups   = dict(groups) 
                    literals = list(literals) 
                    literals.reverse() 
                    if len(groups) == 0:
                        destination = flipdna(dna.__class__(seq=match.expand(value)), _direct=0) 
                    else:
                        destination = dna.__class__(seq="", _direct=0)
                        for index, literal in enumerate(literals): #["A","A",None,"A"]  
                            index = len(literals) - index - 1
                            if literal is None:
                                sub_span = match.span(groups[index])
                                sub_span = (e-sub_span[1], e-sub_span[0])
                                destination = destination + dna[sub_span[0]:sub_span[1]]
                            else:
                                destination = destination + flipdna(dna.__class__(seq=literal, _direct=0))
                    
                    if span[0] == pre_s:
                        segment = segment + destination
                    else:
                        segment = segment + dna[pre_s:span[0]] + destination
                    pre_s   = span[1]
        
        if mode == "edit" and pre_s != e:
            segment = segment + dna[pre_s:e]
    else:
        pre_s = 0  
        match_set = set([]) 
        for match in match_iter:
            result  = {"start":None, "end":None, "strand":None, "match":None} 
            span    = match.span() 
            span    = (span[0] - e if span[0] > e else span[0], span[1] - e if span[1] > e else span[1])
            if span not in match_set:
                match_set.add(span) 
                result["start"]  = s + span[0]
                result["end"]    = s + span[1]
                result["strand"] = strand
                result["match"]  = match
                match_list.append(result) 
                if mode == "edit":
                    groups, literals = sre_parse.parse_template(value, query_pattern)
                    groups = dict(groups)
                    if len(groups) == 0:
                        destination = dna.__class__(seq=match.expand(value), _direct=0)
                    else:
                        destination = dna.__class__(seq="", _direct=0)
                        for index, literal in enumerate(literals):
                            if literal is None:
                                sub_span    = match.span(groups[index])
                                destination = destination + dna[s+sub_span[0]:s+sub_span[1]]
                            else:
                                destination = destination + dna.__class__(seq=literal, _direct=0)
                    if span[0] == pre_s:
                        segment = segment + destination
                    else:
                        segment = segment + dna[s+pre_s:s+span[0]] + destination
                    pre_s = span[1]
        
        if mode == "edit" and s+pre_s != e:
            segment = segment + dna[s+pre_s:e]
    
    if mode == "edit":
        return segment 
    else:
        return match_list

def editsequence(dna, source_sequence, destination_sequence=None, start=0, end=None, strand=1, project=None, product=None, process_name=None, process_description=None, pn=None, pd=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1):
    project             = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description

    dna    = copy.deepcopy(dna) 
    start  = 0 if start == len(dna.seq) else start
    end    = len(dna.seq) if end is None else end
    strand = 1 if strand is None else strand 
    if start == 0 and end == len(dna.seq):
        subject = dna.seq
    else:
        subject = dna.printsequence(start, end, strand)

    _mode = "edit"
    feat_list = [] 
    if source_sequence is None:
        segment = dna.__class__(seq=re.sub(subject, value, subject, _direct=0))   
    else:
        source = source_sequence.upper() 
        query  = source 
        if strand == 1 or strand == -1:
            if destination_sequence is None:
                _mode = "search"
                feature_list = get_matchlist_regex(dna, query, value=None, subject=subject, s=start, e=end, strand=strand) 
            else:
                segment = get_matchlist_regex(dna, query, value=destination_sequence, subject=subject, s=start, e=end, strand=strand) 
        else:
            ValueError("When edit the sequence, the sequence strand to be edit should be '-1' or '+1.'")
    
    if _mode == "edit":
        segment._history_feature = dna._history_feature
        if start == 0 and end == len(dna.seq):
            new_dna = segment
        elif start == 0:
            new_dna = joindna(segment, cropdna(dna,e,len(dna.seq)))
        elif end == len(dna.seq):
            new_dna = joindna(cropdna(dna,0,s), segment)
        else:
            new_dna = joindna(cropdna(dna, 0, s), segment, cropdna(dna, e, len(dna.seq))) 
        
        if dna.topology == "circular":
            new_dna._topology = "circular"
        else:
            pass 

        original_id = dna._product_id
        if project is None:
            new_dna._unique_id = dna._unique_id
        else:
            new_dna._unique_id = project
       
    history_features = [new_dna._history_feature] 
    if type(source_sequence) == new_dna.seq.__class__:
        if source_sequence.parental_class == "DNAFeature":
            qkey = source_sequence.qkey
            for qindex, qfeat in enumerate(new_dna.__class__.queried_features_dict[qkey]):
                if qfeat._second_id == source_sequence.parental_id:
                    break
            if type(source_sequence.item)   == int:
                fsource = "QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , source_sequence.item)
            elif type(source_sequence.item) == slice:
                sl_start = source_sequence.item.start 
                sl_stop  = source_sequence.item.stop  
                sl_step  = source_seqeunce.item.step
                sl_start = "" if sl_start is None else sl_start
                sl_stop  = "" if sl_stop is None else sl_stop
                if sl_step == 1 or sl_step == None:
                    fsource = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop)
                else:
                    fsource = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step)
            else:
                fsource = "QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex)
            history_features.append(source_sequence.parent.subject._history_feature) 

        elif source_sequence.parental_class == "QUEEN": 
            parental_id = source_sequence.parental_id
            if source_sequence.name != None: 
                if "printsequence" in source_sequence.name:
                    if len(source_sequence.name.split("_")) == 2: 
                        seqname = "QUEEN.dna_dict['{}'].printsequence(strand={})".format(parental_id, source_sequence.name.split("_")[-1]) 
                    else:
                        seqname = "QUEEN.dna_dict['{}'].printsequence(start={}, end={}, strand={})".format(parental_id, *source_sequence.name.split("_")[1:])
                elif source_sequence.name == "rcseq":
                    seqname = "QUEEN.dna_dict['{}'].rcseq".format(parental_id) 
            else:
                seqname = "QUEEN.dna_dict['{}'].seq".format(parental_id)
            
            if type(source_sequence.item)   == int:
                fsource = "{}[{}]".format(seqname, source_sequence.item)
            
            elif type(source_sequence.item) == slice:
                sl_start = source_sequence.item.start
                sl_stop  = source_sequence.item.stop 
                sl_step  = source_sequence.item.step
                sl_start = "" if sl_start is None else sl_start
                sl_stop  = "" if sl_stop is None else sl_stop
                if sl_step == 1 or sl_step == None:
                    fsource = "{}[{}:{}]".format(seqname, sl_start, sl_stop)
                else:
                    fsource = "{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step)
            else:
                fsource = "{}".format(seqname)
            history_features.append(source_sequence.parent._history_feature) 
        
        elif source_sequence.parental_class == "Cutsite":
            if source_sequence.parent.name not in cs.defaultkeys:
                cs.new_cutsites.add((source_sequence.parent.name, source_origin.parent.cutsite)) 
            fsource = "cs.lib['{}'].{}".format(source_sequence.parent.name, source_sequence.name)
        else:
            fsourcee = "'{}'".format(source_sequence) 
    else:
        fsource = "'{}'".format(source_sequence)
    
    if __direct == 1:    
        source_sequence      = repr(source_sequence) if source_sequence is not None else None
        destination_sequence = repr(destination_sequence) if destination_sequence is not None else None
        project             = "" #if project is None else ", project='" + project + "'"
        fproduct            = "" if product is None else ", product='" + product + "'"
        process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
        process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
        
        new_dna._product_id = new_dna._unique_id if product is None else product 
        if start == 0 and end == len(dna.seq):
            building_history     = "QUEEN.dna_dict['{}'] = editsequence(QUEEN.dna_dict['{}'], source_sequence={}, destination_sequence={}, strand={}{}{}{}{})".format(new_dna._product_id, original_id, fsource, destination_sequence, strand, project, fproduct, process_name, process_description)
        else:  
            building_history     = "QUEEN.dna_dict['{}'] = editsequence(QUEEN.dna_dict['{}'], source_sequence={}, destination_sequence={}, start={}, end={}, strand={}{}{}{}{})".format(new_dna._product_id, original_id, fsource, destination_sequence, start, end, strand, project, fproduct, process_name, process_description)
        if len(history_features) > 1:
            history_feature = _combine_history(new_dna, history_features) 
            new_dna._history_feature = history_feature
        process_id, original_ids = make_processid(new_dna, building_history, process_id, original_ids)
        add_history(new_dna, [building_history, "source: {}; destination: {}; start: {}; end: {}; strand: {}".format(source_sequence, destination_sequence, start, end, strand), ",".join([process_id] + original_ids)], _sourcefile=_sourcefile)  

    if product is None:
        pass  
    else:
        if _mode == "edit":
            product = product.replace(" ","")
            match   = re.fullmatch("(.+)\[(.+)\]", product) 
            if match:
                if match.group(2).isdecimal() == True:
                    new_dna.__class__._namespace[match.group(1)][int(match.group(2))] = new_dna
                else:
                    new_dna.__class__._namespace[match.group(1)][match.group(2)] = new_dna
            else:    
                new_dna.__class__._namespace[product] = new_dna
        else:
            pass 
    if _mode == "edit":
        return new_dna
    else:
        return feature_list


def _replaceattribute(dna=None, feat_list=None, target_attribute=None, query_re=None, value=None):    
    _exec = 0
    project = dna._product_id 
    if target_attribute[0:8] == "sequence":
        attribute_regex = re.compile("sequence:\![0-9]+\.\.[0-9]+\!") 
        for feat in feat_list:
            tmpid = feat._tmpid
            feat = dna._features_dict[feat._id]
            feat._tmpid = tmpid
        
        topology = dna.topology 
        for feat in feat_list:
            for tmpfeat in dna.dnafeatures:
                if "_tmpid" in tmpfeat.__dict__:
                    if tmpfeat._tmpid == feat._tmpid:
                        feat = tmpfeat
                        break 
                    else:
                        pass 
                else:
                    pass 
            tmpid = feat._tmpid
            fs, fe = feat.start, feat.end 
            if target_attribute == "sequence":
                strand = feat.location.strand
                if feat.start > feat.end:
                    ss, se = 0, len(dna.seq) - feat.start + feat.end 
                else:
                    ss, se = 0, feat.end - feat.start
            
            elif attribute_regex.fullmatch(target_attribute) != None:
                ss, se = tuple(map(int,target_attribute.split("!")[1].split("..")))
                strand = target_attribute.split("!")[2]
                if strand == "":
                    strand = feat.location.strand 
                elif strnad == "+":
                    strand = 1 
                elif strand == "-": 
                    strand = -1 
                            
            if strand == -1:
                s,e = fe - se, fe - ss 
                target_seq = dna.seq[s:e].translate(str.maketrans("ATGC","TACG"))[::-1]
            else:
                s,e = fs + ss, fs + se
                target_seq = dna.seq[s:e] 

            if query_re == "" or query_re is None:
                if value is None:
                    segment = dna.__class__(seq=re.sub(target_seq, "", target_seq), _direct=0)
                else:
                    segment = dna.__class__(seq=re.sub(target_seq, value, target_seq), _direct=0)  
                if strand == -1:
                    segment = flipdna(segment, __direct==0)
            else:
                segment = get_matchlist_regex(dna, query_re, value=value, subject=target_seq, s=s, e=e, strand=feat.strand)
           
            if s == 0 and e == len(dna.seq):
                dna = segment
            elif s == 0:
                dna = joindna(segment, dna[e:len(dna.seq)], topology=dna.topology, __direct=0)
            elif e == len(dna.seq):
                dna = joindna(dna[0:s], segment, topology=dna.topology, __direct=0) 
            else:
                dna = joindna(dna[0:s], segment, dna[e:len(dna.seq)], topology=dna.topology, __direct=0) 
            
            tmpnum = 0 
            tmpremoves = []
            for tmpfeat in dna.dnafeatures:
                if "_tmpid" in tmpfeat.__dict__:
                    if tmpid == tmpfeat._tmpid:
                        if tmpnum == 0:
                            tmpfeat.set_position([s, s+len(segment.seq)], "locaiton")
                            if "broken_feature" in tmpfeat.qualifiers:
                                del tmpfeat.qualifiers["broken_feature"] 
                        else:
                            tmpremoves.append(tmpfeat) 
                        #print(tmpid, s,s+len(segment.seq)) 
                        #print(tmpid, tmpfeat) 
        
        new_dnafeatures = [] 
        for tmpfeat in dna.dnafeatures:
            if tmpfeat in tmpremoves:
                pass
            else:
                if "_tmpid" in tmpfeat.__dict__:
                    del tmpfeat._tmpid
                new_dnafeatures.append(tmpfeat) 
        dna._dnafeatures = new_dnafeatures
        dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))

        _exec += 1
        if set(dna.seq) < set("ATGCRYKMSWBDHVNatgcrykmsbdhvn-"):
            pass 
        else: 
            raise ValueError("Invalid sequence pattern was detected") 

    elif target_attribute == "start":
        for feat in feat_list:
            feat = dna._features_dict[feat._id] 
            feat.subject = dna 
            feat.set_position([value, feat.end], "start") 
            _exec += 1 
    
    elif target_attribute == "end":
        for feat in feat_list:
            feat = dna._features_dict[feat._id] 
            feat.subject = dna 
            feat.set_position([feat.start, value], "end") 
            _exec += 1 

    
    elif target_attribute == "strand":
        for feat in feat_list:
            feat = dna._features_dict[feat._id]
            if value == 1 or value == -1 or value == 0:
                feat.location.strand = value 
                _exec += 1
            else:
                raise ValueError("Strand attribute can take only +1, -1 or 0")
    else:
        new_dnafeatures = []
        _id_all  = [feat._id for feat in dna.dnafeatures] 
        _id_list = [feat._id for feat in feat_list] 
        for feat in dna.dnafeatures:
            _del = 0
            if feat._id in _id_list:
                if target_attribute == "feature_id" or target_attribute == "feature id":
                    if value is None or value == "":
                        _exec += 1
                        _del = 1 
                        pass #It means remove the feature from 'dna.dnafeatures'.
                    else:
                        if value not in _id_all:
                            if len(_id_list) == 1:
                                feat._id = value
                            else:
                                warnings.warn("Multiple features were detected. To ensure the uniqueness of the feature ID, a unique number will be added after the feature ID of each, scuh as '{}_[0-9]'".format(value)) 
                                unique_num = 1
                                while 1:
                                    new_value = value + "_" + str(unique_num) 
                                    if new_value not in _id_all:
                                        feat._id = new_value
                                        _id_all.append(new_value) 
                                        break 
                                    else:
                                        unique_num += 1
                                        pass            
                            _exec += 1
                        else:
                            warnings.warn("The feature ID is already used. please specify another ID")
                        
                elif target_attribute == "feature type" or target_attribute == "feature_type":
                    if query_re is None:
                        new_type  = re.sub(feat.type, value, feat.type)
                    else:
                        new_type  = re.sub(query_re, value, feat.type)
                    feat.type = new_type
                    _exec     += 1

                elif target_attribute[0:len("qualifier:")] == "qualifier:":
                    key = target_attribute.split(":")[-1] 
                    if key in feat.qualifiers:
                        if type(feat.qualifiers[key]) is list:
                            pass 
                        else:
                            feat.qualifiers[key] = [feat.qualifiers[key]]

                        new_elements = [] 
                        for element in feat.qualifiers[key]:
                            if value is None or value == "":
                                pass 
                            else:
                                if query_re is None:
                                    new_element = re.sub(element, value, element)
                                else:
                                    new_element = re.sub(query_re, value, element)
                                new_elements.append(new_element)
                        if len(new_elements) > 0:
                            feat.qualifiers[key] = new_elements
                        else:
                            del feat.qualifiers[key] 
                        _exec += 1 
                    else:
                        pass
            
            if _del == 0:
                new_dnafeatures.append(feat) 
        dna._dnafeatures = new_dnafeatures
        dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))
    
    if _exec == 0:
        warnings.warn("Warning : No target was detected") 
    else:
        return dna 

def replaceattribute(query_re, value=None):
    if value is None:
        value    = query_re
        query_re = None
    return functools.partial(_replaceattribute, query_re=query_re, value=value)

def _removeattribute(dna=None, feat_list=None, target_attribute=None):
    dna = _replaceattribute(dna=dna, feat_list=feat_list, target_attribute=target_attribute, query_re=None, value="")
    return dna 

def removeattribute():
    return functools.partial(_removeattribute)

def _createattribute(dna=None, feat_list=None, target_attribute=None, value=None): 
    new_dnafeatures = copy.copy(dna.dnafeatures)
    _id_all  = [feat._id for feat in dna.dnafeatures]            
    if target_attribute[0:8] == "sequence" or target_attribute == "start" or target_attribute == "end" or target_attribute == "strand" or target_attribute == "feature_type" or target_attribute== "feature type":
        warnings.warn("Warning : 'sequence,' 'start,' 'end,' 'strand,' and 'feature_type' attributes cannote be used in 'createattribute' operateion.")
        return dna 

    elif target_attribute[0:len("qualifier:")] == "qualifier:":
        if value is None:
            pass  
        else:
            key = target_attribute.split(":")[-1]
            _id_list = [feat._id for feat in feat_list if "_id" in feat.__dict__]
            for feat in new_dnafeatures:
                if feat._id in _id_list:
                    if key not in feat.qualifiers:
                        feat.qualifiers[key] = [] 
                    else:
                        if type(feat.qualifiers[key]) is not list:
                            feat.qualifiers[key] = [feat.qualifiers[key]]
                    feat.qualifiers[key].append(value)
            
            for feat in feat_list:
                if feat.location.strand == -1: 
                    s1, e1 = feat.location.parts[-1].start.position, feat.location.parts[0].end.position
                else:
                    s1, e1 = feat.location.parts[0].start.position, feat.location.parts[-1].end.position
                
                if "_id" not in feat.__dict__:
                    flag = 0
                    for feat2 in dna.dnafeatures:
                        if feat2.location.strand == -1: 
                            s2, e2 = feat.location.parts[-1].start, feat.location.parts[0].end
                        else:
                            s2, e2 = feat.location.parts[0].start, feat.location.parts[-1].end

                        if s1 >= e2 and feat2._id.isdecimal() == True:
                            unique_num = 1
                            new_id = str(int(feat2._id) + unique_num) 
                            while new_id in _id_all:
                                unique_num += 1
                                new_id = str(int(feat2._id) + unique_num) 
                            flag = 1
                            break
                 
                    if flag == 0:
                        unique_num = 1
                        new_id = str(unique_num) 
                        while new_id in _id_all:
                            unique_num += 1
                            new_id = str(unique_num) 
                        flag = 1
                    
                    _id_all.append(new_id)
                    feat._id = new_id
                    feat.qualifiers[key] = [] 
                    feat.qualifiers[key].append(value)
                    new_dnafeatures.append(feat) 
    
    elif target_attribute == "feature id" or target_attribute == "feature_id":
        if len(feat_list) > 1:
            _id_all.append(value) 

        for feat in feat_list:
            feat_type  = "misc_feature"
            if len(feat_list) == 1 and value not in _id_all == 1:
                new_id = value
                break

            if value is None:
                value = "" 
                if feat.location.strand == -1: 
                    s1, e1 = feat.location.parts[-1].start.position, feat.location.parts[0].end.position
                else:
                    s1, e1 = feat.location.parts[0].start.position, feat.location.parts[-1].end.position
                flag = 0 
                for feat2 in dna.dnafeatures:
                    if feat2.location.strand == -1: 
                        s2, e2 = feat2.location.parts[-1].start, feat2.location.parts[0].end
                    else:
                        s2, e2 = feat2.location.parts[0].start, feat2.location.parts[-1].end

                    if s1 >= e2 and feat2._id.isdecimal() == True:
                        unique_num = 1
                        new_id = str(int(feat2._id) + unique_num) 
                        while new_id in _id_all:
                            unique_num += 1
                            new_id = str(int(feat2._id) + unique_num) 
                        flag = 1
                        break
                
                if flag == 0:
                    unique_num = 1
                    new_id = str(unique_num) 
                    while new_id in _id_all:
                        unique_num += 1
                        new_id = str(unique_num) 
            
            else:
                unique_num = 1
                new_id = value 
                while new_id in _id_all:
                    new_id = value + "-" + str(unique_num)
                    unique_num += 1
            
            _id_all.append(new_id)
        
            if "_id" not in feat.__dict__ or feat._id is None:
                feat._id  = new_id
                feat.type = feat_type 
                new_dnafeatures.append(feat)
            else:
                new_feat      = feat.__class__(feature=feat, subject=dna) 
                new_feat._id  = new_id 
                new_feat.type = feat_type
                new_dnafeatures.append(new_feat)

    dna._dnafeatures = list(dict([(feat._id, feat) for feat in new_dnafeatures]).values()) 
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna._dnafeatures))) 
    return dna 

def createattribute(value=None):
    return functools.partial(_createattribute, value=value)

def _search(dna, source, query, attribute=None, strand=None): 
    #attributes = "feature ID", "feature type", "start", "end", "sequence", "query:" 
    re_digestion    = re.compile("[ATGCRYKMSWBDHVN]+[\^]?\([0-9]+/[0-9]+\)") 
    attribute_regex = re.compile("sequence:\|[0-9]+\.\.[0-9]+\|[+-]{0,1}")
    if query == ".+" and attribute == "all":
        return source 

    if type(query) is list or type(query) is tuple:
        if type(attribute) is list or type(attribute) is tuple:
            if len(query) != len(attribute):
                raise ValueError("If attribute is given as {}, the list length should be same with query list.".format(type(attribute)))
        else:
            features = []  
            queries  = query
            attributes = [attribute] * len(query) 
            for query, attribute in zip(queries, attributes):
                features.extend(search(dna, query, attribute, min_match, max_mismatch, strand)) 
            return features 
    
    if attribute == "all":
        if type(query) is str and (set(query) <= set("ATGCatgc")):
            attribute = "sequence"
        else:
            pass 

    #if attribute == "start" and attribute == "end":
    #    raise TypeError("Positional value cannot be specified as 'key_attribute'. The 'start', 'end', and 'strand' argument are avaiable to define the sequence range for the search")         
    
    feat_list  = [] 
    if attribute == "sequence":
        for feat in source:
            matches = dna.searchsequence(query, start=feat.start, end=feat.end, strand=feat.strand, _direct=0)
            if len(matches) > 0:
                feat_list.append(feat) 

    if attribute == "all" or attribute == "feature type" or attribute == "feature_type":
        if query is None:
            query = ".+"
        cquery = re.compile(query)
        for feat in source:
            if cquery.match(feat.type) != None or query == feat.type:
                feat_list.append(feat)  

    if attribute == "all" or attribute == "feature id" or attribute == "feature_id":
        if query is None:
            query = ".+"
        cquery       = re.compile(str(query))
        id_feat_list = [(feat.feature_id, feat) for feat in source]  
        for key, feat in id_feat_list:
            if cquery.match(key) != None or query == key:
                feat_list.append(copy.deepcopy(feat)) 

    if attribute == "all" or attribute == "strand":
        for feat in source:
            strand = feat.location.strand
            if query is None or strand == query or (query.isdecimal() and strand == query):
                feat_list.append(feat) 
    
    if attribute == "all" or attribute[0:len("qualifier:")] == "qualifier:":
        if query is None:
            query = ".+"
        cquery = re.compile(query)
        if attribute == "all" or attribute == "qualifier:*":
            for feat in source:
                flag = 0 
                for key in feat.qualifiers:
                    if type(feat.qualifiers[key]) is list:
                        pass 
                    else:
                        feat.qualifiers[key] = [feat.qualifiers[key]]
                    
                    for element in feat.qualifiers[key]:
                        if cquery.search(element) != None or query == element:
                            feat_list.append(feat) 
                            flag = 1
                            break
                    
                    if flag == 1:
                        break
        else:
            key = attribute.split(":")[-1] 
            for feat in source:
                if key in feat.qualifiers:
                    if type(feat.qualifiers[key]) is list:
                        pass 
                    else:
                        feat.qualifiers[key] = [feat.qualifiers[key]]
                    
                    for element in feat.qualifiers[key]:
                        if cquery.search(element) != None or query == element:
                            feat_list.append(feat) 
                            break
    feat_set  = set([])  
    new_feat_list = []
    for feat in feat_list:
        element = (feat.location.start, feat.location.end, feat.sequence, feat.feature_id) 
        if element in feat_set:
            pass 
        else:
            new_feat_list.append(feat)
            feat_set.add(element)
    return new_feat_list

def editfeature(dna, key_attribute="all", query=".+", source=None, start=0, end=None, strand=2, target_attribute=None, operation=None, quine=None, new_copy=True, project=None, product=None, process_name=None, process_description=None, pn=None, pd=None, process_id=None, original_ids=[], _sourcefile=None, __direct=1): 
    project = project if product is None else product
    process_name        = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description
    if quine is None:
        if target_attribute == "sequence":
            quine = True
        else:
            quine = False

    if new_copy == False or quine == False:
        pass 
    else:
        original_id = dna._product_id
        dna = copy.deepcopy(dna) 

    end = len(dna.seq) if end is None else end
    feature_list = dna.dnafeatures if source is None else source
    
    new_source = [] 
    if start > end:
        for feat in feature_list:
            if start <= feat.start <= len(dna.seq) and 0 <= feat.end <= end and (feat.strand == strand or strand == 2): 
                new_source.append(feat)  
    else:
        for feat in feature_list:
            if start <= feat.start <= end and (feat.strand == strand or strand == 2): 
                new_source.append(feat) 

    feat_list = _search(dna, new_source, query, attribute=key_attribute) 
    if operation is None or target_attribute is None:
        return feat_list
    
    elif operation.func.__name__ in ("_createattribute", "_removeattribute", "_replaceattribute"):
        for i, feat in enumerate(feat_list):
            feat._tmpid = i 
        
        dna   = operation(dna=dna, feat_list=feat_list, target_attribute=target_attribute)
        largs = [] 
        for item in operation.keywords.items():
            item = list(item) 
            if type(item[1]) is str:
                item[1] = repr(item[1]) 
                #"'" + item[1] + "'"
            else:
                item[1] = str(item[1]) 
            largs.append("=".join(item)) 
        command = operation.func.__name__[1:] + "(" + ",".join(largs) + ")"
        
        if type(query) == dna.seq.__class__:
            if query.parental_class == "DNAFeature":
                qkey = left_origin.qkey
                for qindex, qfeat in enumerate(dna.__class__.queried_features_dict[qkey]):
                    if qfeat._second_id == query.parental_id:
                        break
                if type(query.item)   == int:
                    fquery = "QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , query.item)
                elif type(query.item) == slice:
                    sl_start = query.item.start
                    sl_stop  = query.item.stop 
                    sl_step  = query.item.step
                    sl_start = "" if sl_start is None else sl_start
                    sl_stop  = "" if sl_stop is None else sl_stop
                    if sl_step == 1 or sl_step == None:
                        fquery = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop)
                    else:
                        fquery = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step)
                else:
                    fquery = "QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex)
                history_features.append(query.parent.subject._history_feature) 
            
            elif query.parental_class == "QUEEN": 
                parental_id = query.parental_id 
                if query.name != None: 
                    if "printsequence" in query.name:
                        if len(query.name.split("_")) == 2: 
                            seqname = "QUEEN.dna_dict['{}'].printsequence(strand={})".format(parental_id, query.name.split("_")[-1]) 
                        else:
                            seqname = "QUEEN.dna_dict['{}'].printsequence(start={}, end={}, strand={})".format(parental_id, *query.name.split("_")[1:])
                    elif query.name == "rcseq":
                        seqname = "QUEEN.dna_dict['{}'].rcseq".format(parental_id) 
                else:
                    seqname = "QUEEN.dna_dict['{}'].seq".format(parental_id)
                if type(query.item)   == int:
                    args.append("QUEEN.dna_dict['{}'].seq[{}]".format(parental_id, query.item))
                elif type(query.item) == slice:
                    sl_start = query.item.start
                    sl_stop  = query.item.stop 
                    sl_step  = query.item.step
                    sl_start = "" if sl_start is None else sl_start
                    sl_stop  = "" if sl_stop is None else sl_stop
                    if sl_step == 1 or sl_step == None:
                        fquery = "{}[{}:{}]".format(seqname, sl_start, sl_stop)
                    else:
                        fquery = "{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step)
                else:
                    fquery = "{}".format(seqname)
                history_features.append(query.parent._history_feature) 
            
            elif query.parental_class == "Cutsite":
                if query.parent.name not in cs.defaultkeys:
                    cs.new_cutsites.append((query.parent.name, query.parent.cutsite)) 
                fquery = "cs.lib['{}'].{}".format(qorigin.parent.name, qorigin.name) 
            else:
                fquery = "{}".format(repr(query))  

        else:
            fquery = "{}".format(repr(query))  
        
        if source is not None:
            qkeys = set([]) 
            for feat in source:
                if "_qkey" in feat.__dict__:
                    qkeys.add(feat._qkey)
            
            if len(set(qkeys)) == 1:
                source = "QUEEN.queried_features_dict['{}']".format(list(qkeys)[0])
            else:
                pass 

        if project is None:
            pass
        else:
            dna._unique_id = project
                            
        if __direct == 1 and new_copy == True and quine == True:
            project             = "" #if project is None else ", project='" + project + "'"
            fproduct            = "" if product is None else ", product='" + product + "'"
            process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
            process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
            dna._product_id = dna._unique_id if product is None else product 
            if start == 0 and end == len(dna.seq):
                args = [key_attribute, fquery, source, strand, target_attribute, command, new_copy]
                for i in range(len(args)):
                    if type(args[i]) is str and i != 5 and i != 1:
                        args[i] = "'" + args[i] + "'" 
                building_history = "QUEEN.dna_dict['{}'] = editfeature(QUEEN.dna_dict['{}'], key_attribute={}, query={}, source={}, strand={}, target_attribute={}, operation={}, new_copy={}{}{}{}{})".format(dna._product_id, original_id, *args, project, fproduct, process_name, process_description) 
            
            else:
                args = [key_attribute, fquery, source, start, end, strand, target_attribute, command, new_copy]
                for i in range(len(args)):
                    if type(args[i]) is str and i != 7 and i != 1:
                        args[i] = "'" + args[i] + "'" 
                building_history = "QUEEN.dna_dict['{}'] = editfeature(QUEEN.dna_dict['{}'], key_attribute={}, query={}, source={}, start={}, end={}, strand={}, target_attribute={}, operation={}, new_copy={}{}{}{}{})".format(dna._product_id, original_id, *args, project, fproduct, project, process_name, process_description) 
            
            process_id, original_ids = make_processid(dna, building_history, process_id, original_ids)
            add_history(dna, [building_history, "key_attribute: {}; query: {}; start: {}; end: {}; strand: {}; target_attribute: {}; operation: {}".format(key_attribute, fquery, start, end, strand, target_attribute, command), process_id, ",".join([process_id] + original_ids)], _sourcefile=_sourcefile)
            if product is None:
                pass 
            else:
                dna._product_id = product 
                match   = re.fullmatch("(.+)\[(.+)\]", product)
                if match:
                    if match.group(2).isdecimal() == True:
                        dna.__class__._namespace[match.group(1)][int(match.group(2))] = dna
                    else:
                        dna.__class__._namespace[match.group(1)][match.group(2)] = dna
                else:    
                    dna.__class__._namespace[product] = dna

        elif new_copy == False or quine == False:
            project             = "" if project is None else ", project='" + project + "'"
            fproduct            = "" if product is None else ", product='" + product + "'"
            process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
            process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
            args = [key_attribute, fquery, source, start, end, strand, target_attribute, command, new_copy]
            for i in range(len(args)):
                if type(args[i]) is str and i != 7 and i != 1:
                    args[i] = "'" + args[i] + "'" 
            building_history = "editfeature(QUEEN.dna_dict['{}'], key_attribute={}, query={}, source={}, start={}, end={}, strand={}, target_attribute={}, operation={}, new_copy={}{}{}{}{})".format(dna._unique_id, *args, project, fproduct, process_name, process_description)
            #process_id, original_ids = make_processid(dna, building_history, process_id, original_ids)
            #add_history(dna, [building_history, "key_attribute:{}; query:{}; start:{}; end:{}; strand:{}; target_attribute:{}; operation:{}".format(key_attribute, fquery, start, end, strand, target_attribute, command), process_id])
    
    else:
        raise ValueError("'operation' can take only one of 'createattribute,' 'removeattribute,' and 'replaceattribute.'")

    if new_copy == True:
        return dna 

def _circularizedna(dna):
    dna = copy.deepcopy(dna)
    seq_origin = dna.seq
    feats_origin = dna.dnafeatures
    if dna.topology == "circular" and dna.record.annotations["topology"] == "circular":
        print("The QUEEN object topology is circular")

    if (dna._right_end_top * dna._left_end_bottom == 1 and dna._right_end_bottom * dna._left_end_top == 1) and len(dna._right_end) > 0 and (dna._left_end_top == -1 or dna._left_end_bottom == -1):
        if str(dna._right_end) == str(dna._left_end): 
            ovhg       = dna._right_end
            subdna     = cropdna(dna,0,len(dna.seq)-len(dna._right_end),__direct=0)
            dna._seq   = subdna.seq
            dna.record = subdna.record
        else:
            return False
    else:
        ovhg = ""
    dna._topology = "circular"
    
    remove_list = [] 
    feats1      = [feat for feat in dna.dnafeatures if "broken_feature" in feat.qualifiers]
    for feat1 in feats1:
        s1, e1 = feat1.start, feat1.end
        for feat2 in feats1:
            s2, e2 = feat2.start, feat2.end
            if feat1 == feat2 or feat1 in remove_list:
                pass 
            
            elif feat1.type == feat2.type:
                flag = 0
                for key in feat1.qualifiers:
                    if key == "broken_feature":
                        pass 
                    elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                        flag = 1
                    else:
                        flag = 0
                        break   
                
                if flag == 1 and "broken_feature" in feat1.qualifiers and "broken_feature" in feat2.qualifiers:
                    note1   = feat1.qualifiers["broken_feature"][0]
                    label   = ":".join(note1.split(":")[:-1])
                    length1 = int(note1.split(":")[-4])
                    pos_s1  = int(note1.split(":")[-1].split("..")[0].replace(" ",""))
                    pos_e1  = int(note1.split(":")[-1].split("..")[1].replace(" ","")) 
                    
                    note2   = feat2.qualifiers["broken_feature"][0] 
                    length2 = int(note2.split(":")[-4])
                    pos_s2  = int(note2.split(":")[-1].split("..")[0].replace(" ",""))
                    pos_e2  = int(note2.split(":")[-1].split("..")[1].replace(" ","")) 
                    if (s1 >= e2) and length1 == length2 and "_original" in feat1.__dict__ and "_original" in feat2.__dict__ and feat1.original == feat2.original and feat1.location.strand == feat2.location.strand:
                        note        = "{}:{}..{}".format(label, pos_s1, pos_e2)
                        new_seq     = seq_origin[s1:e1] + seq_origin[s2:e2]
                        feat1_index = dna.dnafeatures.index(feat1)
                        new_feat    = copy.deepcopy(dna.dnafeatures[feat1_index]) 
                        strand      = new_feat.location.strand
                        if len(feat1.location.parts) == 1 and len(feat2.location.parts) == 1:
                            new_feat.location = FeatureLocation(feat1.location.parts[0].start.position, len(dna.seq) + feat2.location.parts[-1].end.position, feat1.strand)
                            new_feat.location.strand = strand
                        else:
                            feat2_parts = [(p.start.position + len(dna.seq), p.end.position + len(dna.seq), feat2.strand) for p in feat2.location.parts]
                            locations   = feat1.location.parts[0:-1] + [FeatureLocation(feat1.location.parts[-1].start.position, len(dna.seq) + feat2.location.parts[0].end.position, feat1.strand)] + feat2_parts[0:-1]
                            if strand == -1:
                                locations.reverse() 
                            new_feat.location = CompoundLocation(locations) 
                            new_feat.location.strand = strand 
                        
                        new_feat  = feat1.__class__(feature=new_feat, subject=dna)
                        new_feat1 = feat1.__class__(feature=feat1, subject=dna)
                        new_feat2 = feat1.__class__(feature=feat2, subject=dna) 
                        s = new_feat.start 
                        e = new_feat.end if new_feat.end <= len(dna.seq) else new_feat.end - len(dna.seq)
                        if new_feat._original == dna.printsequence(new_feat1.start, new_feat2.end, new_feat.location.strand if new_feat.location.strand !=0 else 1):
                            dna._dnafeatures[feat1_index].qualifiers["broken_feature"] = [note]
                            if len(new_seq) - len(ovhg) == length1:
                                del dna._dnafeatures[dna.dnafeatures.index(feat1)].qualifiers["broken_feature"]
                            dna._dnafeatures[feat1_index].location = new_feat.location
                            if feat2 in dna._dnafeatures:
                                dna._dnafeatures.remove(feat2) 
                                remove_list.append(feat2) 
    
    for i in range(len(dna.dnafeatures)):    
        if dna.dnafeatures[i].location.parts[-1].end.position > len(dna.seq):
            if dna.dnafeatures[i].location.parts[0].start.position >= len(dna.seq):
                strand                      = dna.dnafeatures[i].location.strand
                dna._dnafeatures[i].location = FeatureLocation(dna.dnafeatures[i].location.parts[0].start.position-len(dna.seq),dna.dnafeatures[i].location.parts[-1].end.position-len(dna.seq))
                dna._dnafeatures[i].location.strand = strand
            else:
                strand    = dna.dnafeatures[i].location.strand
                locations = [FeatureLocation(dna.dnafeatures[i].location.parts[0].start.position,len(dna.seq)), FeatureLocation(0,dna.dnafeatures[i].location.parts[-1].end.position-len(dna.seq))]
                if strand == -1:
                    locations.reverse()   
                dna._dnafeatures[i].location = CompoundLocation(locations)
                dna._dnafeatures[i].location.strand = strand
            dna._dnafeatures[i] = dna.dnafeatures[i].__class__(feature=dna.dnafeatures[i], location=dna.dnafeatures[i].location)
    
    #dna.record.features = _assigndnafeatures(dna.dnafeatures)
    dna._left_end  = ""
    dna._left_end_top    = 0 
    dna._left_end_bottom = 0 

    dna._right_end = ""
    dna._right_end_top    = 0 
    dna._right_end_bottom = 0 

    if Alphabet:
        dna.record.seq = Seq(str(dna.seq),Alphabet.DNAAlphabet())
    else:
        dna.record.seq = Seq(str(dna.seq))
    dna.record.annotations["topology"] = dna.topology
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))
    return dna

def visualizemap(dna, map_view="linear", feature_list=None, start=0, end=None,label_location=None, display_label=2, display_title=True, display_axis=True, fontsize=None, tick_interval="auto", labelcolor="k", title=None, width_scale="auto", height_scale=1.0, linebreak=None, seq=False, diamater_scale=1.0, fig= None):
    if fontsize is None and map_view == "linear":
        fontsize = 12
    elif fontsize is None and map_view == "circular":
        fontsize = 10
    else:
        pass 

    if title is None or title == "":
        display_titlee = False

    if feature_list is None:
        feature_list = dna.dnafeatures
        if map_view == "circular":
            feature_list.sort(key=lambda x:len(dna.printsequence(x.start, x.end)))
    
    standard_scale = 4000
    if map_view == "circular":
        fig, axes = vc.visualize(dna, format=0, feature_list=feature_list, unvisible_types=["source"], visible_types=[],
                                 bottom=400 * diamater_scale, label_visible=display_label, fontsize=fontsize, 
                                 title_visible=display_title, axis_visible=display_axis, tick_space=tick_interval, labelcolor=labelcolor, titlename=title, fig=fig)
    
    else:
        fig, ax = vl.visualize(dna, start=start, end=end, feature_list=feature_list, wrap_width=linebreak, annotation_loc=label_location, unvisible_types=["source"], 
                               visible_types=[], enlarge_w=width_scale, enlarge_h=height_scale, fontsize=fontsize, with_seq=seq, nucl_char=None, nucl_color_dict=None, label_visible=display_label, 
                               scale="fix", title_visible=display_title, axis_visible=display_axis, tick_space=tick_interval, labelcolor=labelcolor, titlename=title, fig=fig)
    return fig

