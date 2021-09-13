import os 
import io 
import sys
import copy
import math
import argparse
import collections
import configparser
import Levenshtein
import functools
import warnings
import sre_parse
import regex as re

from Bio import SeqIO, pairwise2, BiopythonParserWarning
try:
    from Bio import Alphabet #For Biopython < 1.78
except ImportError:
    Alphabet = False #For Biopyhton ≥ 1.78 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, FeatureLocation, ExactPosition

sys.path.append("/".join(__file__.split("/")[:-1]))
from Qobj import *
import visualize_circular_dna as vc
import visualize_linear_dna as vl 
warnings.simplefilter('ignore', BiopythonParserWarning)

_namespace = globals()
_namespaceflag = 0

def _combine_history(dna, history_features):
    history_feature = SeqFeature(FeatureLocation(0, len(dna.seq), strand=0), type="source")
    hsitory_feature = DNAfeature(history_feature, subject=dna.seq)
    for feat in history_features:
        for key in feat.qualifiers:
            if "building_history" in key:
                history_feature.qualifiers[key] = feat.qualifiers[key]  
    return history_feature

def _add_history(dna, history=""): 
    flag = 0 
    if dna._history_feature is not None:
        QUEEN._num_history += 1 
        dna._history_feature.qualifiers["building_history_{}".format(QUEEN._num_history)] = [history.replace(" ","–")] 
        flag = 1

    if flag == 0:
        _id_all = [int(feat._id) for feat in dna.dnafeatures if str(feat._id).isdecimal()]
        _id_all.append(0)
        if len(_id_all) == 1:
            max_id = str(0)
        else:
            max_id  = str(100 + math.floor(max(_id_all) / 100) * 100)
        feat = SeqFeature(FeatureLocation(0, len(dna.seq), strand=1), type="source") 
        feat._id = max_id
        feat.qualifiers["label"]       = [dna.project] 
        feat.qualifiers["description"] = ["Record of building history"]
        QUEEN._num_history += 1 
        feat.qualifiers["building_history_{}".format(QUEEN._num_history)] = [history.replace(" ","–")] 
        feat = DNAfeature(feature=feat, subject=dna)         
        dna._history_feature = feat 
        #dna.dnafeatures.append(feat)
        #dna._features_dict[feat._id] = feat
        topfeat = feat 
    
    dna._features_dict  = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))
    dna.record.features = dna.dnafeatures

def set_namespace(_globals=None):
    global _namespace
    global _namespaceflag  
    if _globals is None:
        pass
    else:
        _namespace = _globals
        _namespaceflag = 1 

def _slide(feats,slide):
    new_feats = []
    for feat in feats:
        feat = copy.deepcopy(feat)
        strand = feat.location.strand
        for p in range(len(feat.location.parts)):
            feat.location.parts[p]._start = ExactPosition(feat.location.parts[p].start.position + slide)
            feat.location.parts[p]._end   = ExactPosition(feat.location.parts[p].end.position + slide)
        feat.location.strand = strand
        new_feats.append(DNAfeature(feat))
    return new_feats 

def cutdna(dna, *positions, crop=False, project=None, product=None, process_description=None, __direct=1):
    dna = copy.deepcopy(dna)
    if process_description is None:
        process_description = QUEEN.process_description
    else:
        QUEEN.process_description = process_description 
   
    for feat in dna.dnafeatures:
        if "_start" not in feat.__dict__:
            print("hogehoge")
            print(type(feat)) 
            print(feat.__dict__) 
            print()

    def extract(dna, start, end, project=None): 
        start_top    = start[0] 
        start_bottom = start[1] 
        start = min(start)
        
        end_top    = end[0]
        end_bottom = end[1] 
        end = max(end)

        if start == 0 and end == len(dna.seq):
            new_dna = copy.copy(dna)
            new_dna.topology = "linear"
            if project is None:
                project = dna.project
            new_dna.project = project
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
            if start > end and dna.topoloy == "linear":
                raise ValueError("Start value should be larger than or equal to end value.")
            feats = []
            new_features = []
            
            #Linearization
            for feat in dna.dnafeatures:
                strand = feat.location.strand
                s = feat.start
                e = feat.end
                if "_original" not in feat.__dict__:
                    if s > e:
                        feat._original = str(dna.seq)[s:len(dna.seq)] + str(dna.seq)[:e]
                    else:
                        feat._original = str(dna.seq)[s:e].upper()
                    #feat.original_location = feat.location
                if s > e:
                    if len(feat.location.parts) == 1:
                        length = len(dna.seq) - s + e
                        locations = [FeatureLocation(s,len(dna.seq)),FeatureLocation(0,e)]
                        if strand == -1:
                            locations.reverse()
                        feat.location = CompoundLocation(locations)
                        feat.location.strand = strand

                    strand = feat.location.strand
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
                    new_features.append(DNAfeature(feature=feat1))
                    new_features.append(DNAfeature(feature=feat2))
                
                else:
                    #print(feat, start, end) 
                    new_features.append(DNAfeature(feature=feat))
            #Cropping
            for feat in new_features:
                strand = feat.strand
                s = feat.start
                e = feat.end
                if "_original" not in feat.__dict__:
                    if s > e:
                        feat._original = str(dna.seq)[s:len(dna.seq)] + str(dna.seq)[:e]
                    else:
                        feat._original = str(dna.seq)[s:e].upper()    
                    #feat.original_location = feat.location
                feat = copy.deepcopy(feat)
                if len(feat.location.parts) == 1 and s <= e:
                    if e > start and s < end:
                        if s - start < 0:
                            feat.location.parts[0]._start = ExactPosition(0)
                            if "broken_feature" not in feat.qualifiers:
                                label = feat._id
                                if feat.feature_type == "source":
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
                                if feat.feature_type == "source":
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(original_seq), feat.original, s, e))
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
                        feats.append(DNAfeature(feature=feat))
                
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
                                if feat.feature_type == "source":
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "{}".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), feat.original, s, e))
                                #if strand is None:
                                #    print(feat) 
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
                                
                                if feat.feature_type == "source":
                                    original_seq = "-"
                                else:
                                    original_seq = feat.original
                                
                                if feat.feature_type == "CDS" and "translation" in feat.qualifiers:
                                    del feat.qualifiers["translation"]

                                label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), feat.original, s, e))
                                feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, end-s)]
                            else:
                                s      = int(locations[0][0])
                                note   = feat.qualifiers["broken_feature"][0]
                                if strand  >= 0:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s + (end-start-s)-1)
                                else:
                                    label  = ":".join(note.split(":")[:-1])
                                    length = int(note.split(":")[-4]) 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
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
                        feats.append(DNAfeature(feature=feat))
        
            feats.sort(key=lambda x:(x.location.parts[0].start.position, x.location.parts[-1].end.position))
           
            
            subdna = QUEEN(seq=str(dna.seq[start:end]), _direct=0)
            subdna._history_feature = copy.deepcopy(dna._history_feature) 
            subdna.dnafeatures = feats

            subdna._features_dict = dict(list(map(lambda x:(x._id, x), subdna.dnafeatures)))
            subdna.topology = "linear"
            
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
            
        if project is None:
            project = dna.project
        subdna.project = project
      
        if start_top != start_bottom or end_top != end_bottom:
            if (start_dif := start_top - start_bottom) > 0:
                left = "-" * start_dif + "/" +  "*" * start_dif
            elif start_dif < 0:
                left = "*" * abs(start_dif) + "/" +  "-" * abs(start_dif) 
            else:
                left = ""
            
            if (end_dif := end_top - end_bottom) > 0:
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
        
        elif type(pos) is SeqFeature or type(pos) is DNAfeature:
            strand = pos.location.strand
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
        
        elif type(pos) is int or type(pos) is Qint:
            pos = (pos, pos)  
            spos, epos = pos
            spos = spos - len(dna.seq) if spos > len(dna.seq) else spos 
            epos = epos + len(dna.seq) if epos < 0 else epos
            new_positions.append((spos,epos))  
        
    if crop == True:
        crop_positions = (new_positions[0], new_positions[1])
    
    if dna.topology == "linear":
        if (0,0) not in new_positions:
            positions = [(0,0)]  + new_positions 
        else: 
            positions = new_positions 
        
        if (len(dna.seq),len(dna.seq)) not in positions:
            positions = positions + [(len(dna.seq), len(dna.seq))] 
        positions = list(positions) 
        positions.sort() 

    elif dna.topology == "circular":
        positions = new_positions
        positions = list(positions) 
        tmp_positions = positions[:]
        positions.sort() 
        for pindex, pos in enumerate(positions):
            if pos == tmp_positions[0]:
                new_positions = positions[pindex:] + positions[:pindex]
                break 
        positions = new_positions

    if dna.topology == "linear":
        if crop == True:
            dnas.append(extract(dna, crop_positions[0], crop_positions[1], project=project))
        else:
            for i, pos in enumerate(positions[0:-1]):
                dnas.append(extract(dna,pos,positions[i+1],project=project))
    
    elif dna.topology == "circular":
        if crop == True: 
            dnas.append(extract(dna, crop_positions[0], crop_positions[1], project=project))
        else:
            for i, pos in enumerate(positions[0:-1]):
                dnas.append(extract(dna, pos, positions[i+1], project=project))
            dnas.append(extract(dna, positions[-1], positions[0], project=project)) 

    if project is None:
        project = dna.project
    
    for new_dna in dnas:
        new_dna.project = project

    if __direct == 1:
        products = []
        dna_keys = list(QUEEN.dna_dict.keys())
        for i in range(len(dnas)):
            if dnas[i].project in QUEEN.dna_dict:
                if project.split("_")[-1].isdecimal() == True:
                    project = "_".join(project.split("_")[:-1])
                keys   = list(QUEEN.dna_dict.keys())
                unique = 0
                while project + "_" + str(unique) in dna_keys:
                    unique += 1
                dnas[i]._unique_id = project + "_" + str(unique)
            else:         
                dnas[i]._unique_id = project
            dna_keys.append(dnas[i]._unique_id)
            if product is not None:
                QUEEN.dna_dict[dnas[i]._unique_id] = product + "_" + str(i)  
            else:
                QUEEN.dna_dict[dnas[i]._unique_id] = product 
            products.append("QUEEN.dna_dict['{}']".format(dnas[i]._unique_id))

        args = [] 
        history_features = [dnas[0]._history_feature] 
        for new_pos, pos in zip(new_positions, positions):
            if type(pos) == DNAfeature:
                qkey = pos._qkey
                for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                    if qfeat._second_id == pos._second_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}]".format(qkey, qindex))
                history_features.append(pos.subject._history_feature) 

            elif type(pos) == Qint:
                qkey = pos.qkey
                for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                    if qfeat._second_id == pos.parental_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}].{}".format(pos.qkey, qindex, pos.name))
                history_features.append(pos.parent.subject._history_feature) 

            else:
                new_pos = "/".join(list(map(str,new_pos)))
                new_pos = "'" + new_pos + "'"
                args.append(new_pos)
        
        if type(project) is str:
            project = "'" + project + "'"
        if type(process_description) is str:
            process_description = "'" + process_description + "'"
        
        if len(products) > 1:
            building_history = "{}  = cutdna(QUEEN.dna_dict['{}'], {}, crop={}, project={}, product={}, process_description={})".format(", ".join(products), dna._unique_id, ", ".join(args), str(crop), project, product if product is None else "'" + product + "'", process_description) 
        else:
            building_history = "{}, = cutdna(QUEEN.dna_dict['{}'], {}, crop={}, project={}, product={}, process_description={})".format(", ".join(products), dna._unique_id, ", ".join(args), str(crop), project, product if product is None else "'" + product + "'", process_description)
        
        for subdna in dnas:
            history_feature = _combine_history(subdna, history_features)
            subdna._history_feature = history_feature
            _add_history(subdna, building_history)
    
    if product is None:
        pass 
    else:
        _namespace[product] = dnas

    if crop == True:
        return dnas[0], crop_positions 
    else:
        return dnas

def cropdna(dna, start=0, end=None, project=None, product=None, process_description=None, __direct=1):
    """
    Extract subsequence from the specified region in the genome.
    All features within the subsequence including truncated features are carried over to extracrted 'dna' object. 
    """ 
    if process_description is None:
        process_description = QUEEN.process_description
    else:
        QUEEN.process_description = process_description 

    if project is None:
        project = dna.project
    
    if end is None:
        end = len(dna.seq) 

    subdna, crop_positions = cutdna(dna, start, end, project=project, crop=True, __direct=0)  
    if __direct == 1:
        historis = [subdna._history_feature] 
        if project in QUEEN.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(QUEEN.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in QUEEN.dna_dict:
                unique += 1
            subdna._unique_id = project + "_" + str(unique)
        else:         
            subdna._unique_id = project
        
        args = []
        history_features = [subdna._history_feature] 
        for new_pos, pos in zip(crop_positions, (start, end)):
            if type(pos) == DNAfeature:
                qkey = pos._qkey
                for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                    if qfeat._second_id == pos._second_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}]".format(qkey, qindex))
                history_features.append(pos.subject._history_feature) 
            
            elif type(pos) == Qint:
                qkey = pos.qkey
                for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                    if qfeat._second_id == pos.parental_id:
                        break
                args.append("QUEEN.queried_features_dict['{}'][{}].{}".format(pos.qkey, qindex, pos.name))
                history_features.append(pos.parent.subject._history_feature) 
        
            else:
                new_pos = "/".join(list(map(str,new_pos)))
                new_pos = "'" + new_pos + "'"
                args.append(new_pos)
        
        args.append(project)
        args.append(process_description) 
        if type(args[-2]) is str:
            args[-2] = "'" + args[-2] + "'" 
        
        if type(args[-1]) is str:
            args[-1] = "'" + args[-1] + "'" 
        
        QUEEN.dna_dict[subdna._unique_id] = product
        building_history = "QUEEN.dna_dict['{}'] = cropdna(QUEEN.dna_dict['{}'], start={}, end={}, project={}, product={}, process_description={})".format(subdna._unique_id, dna._unique_id, args[0], args[1], args[2], product if product is None else "'" + product + "'", args[3])
        
        history_feature = _combine_history(subdna, history_features)
        subdna._history_feature = history_feature
        _add_history(subdna, building_history)
    
    if product is None:
        pass 
    else:
        _namespace[product] = subdna
    return subdna

def joindna(*dnas, topology="linear", project=None, product=None, process_description=None, __direct=1):
    """
    Join dna objects. 
    When ovhg_check argument is True, adjacent QUEEN objects should have common overhang sequences.  
    """
    if process_description is None:
        process_description = QUEEN.process_description
    else:
        QUEEN.process_description = process_description 

    new_dnas = [] 
    for i, dna in enumerate(dnas):
        if type(dna) is not QUEEN:
            dna = QUEEN(dna)
        if dna.topology == "circular":
            if i == 0:
                order = "first"
            elif i == 1:
                order = "second"
            elif i == 2:
                order = "third"
            else:
                order = str(i) + "th" 
            raise ValueError("The {} QUEEN object topology is 'circular'. Cicular QUEEN objects cannot be joined with others.".format(order)) 
        new_dnas.append(dna) 
    
    dnas = new_dnas
    #Extract history information
    history_features = [] 
    for dna in dnas:
        history_features.append(dna._history_feature)

    construct = copy.deepcopy(dnas[0])
    if len(dnas) > 1:
        for dna in dnas[1:]:
            feats   = dna.dnafeatures 
            if (dna._left_end_top * construct._right_end_bottom == 1 or dna._left_end_bottom * construct._right_end_top == 1) and ((dna._left_end_top == -1 or dna._left_end_bottom == -1) or (construct._right_end_top == -1 or construct._right_end_bottom == -1)):
                if dna._left_end_top == 1:
                    sticky_end = dna._left_end 
                else:
                    sticky_end = construct._right_end
                
                ovhg    = dna._left_end
                new_dna = cropdna(dna,len(ovhg),len(dna.seq),__direct=0) 
                if dna._left_end == construct._right_end:
                    pass
                    #print("The DNA objects were joined based on complementary sticky end of each fragment. The sticky end is '{}'".format(sticky_end)) 
                else:
                    print("These QUEEN objects were not able to be joined. Please check sticky end seqeunces of the DNA objects")
                    return False

            else:
                new_dna = dna
                ovhg = ""
            
            feats  = _slide(feats, len(construct.seq) - len(ovhg))
            feats1 = [feat for feat in construct.dnafeatures if "broken_feature" in feat.qualifiers]
            feats2 = [feat for feat in feats if "broken_feature" in feat.qualifiers]

            construct.seq  = construct.seq + new_dna.seq 
            const_features = copy.copy(construct.dnafeatures) 
            #Recover a original feature from fragmented features
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
                                   
                                    new_feat  = DNAfeature(feature=new_feat, subject=construct)
                                    new_feat1 = DNAfeature(feature=feat1, subject=construct)
                                    new_feat2 = DNAfeature(feature=feat2, subject=construct) 
                                    s = new_feat.start
                                    e = new_feat.end        
                                    if construct.getdnaseq(new_feat1.start, new_feat2.end, new_feat.location.strand if new_feat.location.strand !=0 else 1) in new_feat.original:
                                        construct.dnafeatures[feat1_index] = DNAfeature(feature=new_feat)
                                        construct.dnafeatures[feat1_index].qualifiers["broken_feature"] = [note]
                                        del feats[feats.index(feat2)] 
            
            construct.dnafeatures = construct.dnafeatures + feats
            construct._right_end        = dna._right_end
            construct._right_end_top    = dna._right_end_top
            construct._right_end_bottom = dna._right_end_bottom
            construct.topology = "linear"
            ovhg = dna._right_end

        construct.dnafeatures.sort(key=lambda x:x.location.parts[0].start.position)
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
            zero_positions = [] 
            for d, dna in enumerate(dnas):
                if 0 in dna._positions:
                    zero_positions.append((len(dna.seq),d,dna._positions.index(0)))
            if len(zero_positions) > 0:
                zero_positions.sort() 
                zero_positions.reverse() 
                zero_position = 0
                for dna in dnas[0:zero_positions[0][1]]:
                    zero_position += len(dna.seq)
                zero_position += zero_positions[0][2]
                construct = cutdna(construct, zero_position)[0]
                construct = _circularizedna(construct) 
            construct._positions = tuple(range(len(construct.seq)))    
        else:
            construct._positions = tuple(range(len(construct.seq)))
        
        construct._setfeatureid() #Update feature ID
    else:
        construct = _circularizedna(dnas[0])
        construct._positions = construct._positions[0:len(construct.seq)]

    if project is None:
        project = dnas[0].project

    #Recover fragmented features if complete sequence is in the construct.
    new_features = [] 
    remove_features = [] 
    for feat in construct.dnafeatures:
        if "broken_feature" in feat.qualifiers:
            note       = feat.qualifiers["broken_feature"][0]
            label      = ":".join(note.split(":")[:-1])
            poss, pose = list(map(int,note.split(":")[-1].split("..")))
            length = int(note.split(":")[-4]) 
            if feat.location.strand != -1:
                sfeat = feat.start-(poss-1) 
                sfeat = sfeat if sfeat > 0 else len(construct.seq) - sfeat
                efeat = feat.end+(length-pose)
            else:
                sfeat = feat.start-(length-pose) 
                sfeat = sfeat if sfeat > 0 else len(construct.seq) - sfeat
                efeat = feat.end+(poss-1)    
            
            if feat.original == construct.getdnaseq(sfeat, efeat, strand=1):
                if sfeat < efeat:
                    location = FeatureLocation(sfeat, efeat, feat.location.strand) 
                else:
                    location = CompoundLocation([FeatureLocation(sfeat, len(construct.seq)), FeatureLocation(0, efeat, feat.location.strand)])  
                newfeat = DNAfeature(location=location, subject=construct)
                newfeat.type = feat.type
                newfeat.qualifiers = feat.qualifiers
                del newfeat.qualifiers["broken_feature"]
                newfeat._id = feat.feature_id
                new_features.append(newfeat)
                remove_features.append(feat)

    for feat in remove_features:   
        del construct.dnafeatures[construct.dnafeatures.index(feat)]
    
    for feat in new_features:
        construct.dnafeatures.append(feat) 

    construct.record.id = project
    construct.project   = project
    construct._features_dict = dict(list(map(lambda x:(x._id, x), construct.dnafeatures)))
    construct.record.feartures = construct.dnafeatures
    if __direct == 1:
        if project in QUEEN.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(QUEEN.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in QUEEN.dna_dict:
                unique += 1
            construct._unique_id = project + "_" + str(unique)
        else:         
            construct._unique_id = project
        
        QUEEN.dna_dict[construct._unique_id] = product
        
        if type(process_description) is str:
            process_description = "'" + process_description + "'"
        
        if type(project) is str:
            project =  "'" + project + "'" 
        
        dna_elements = "[" + ", ".join(["QUEEN.dna_dict['{}']".format(dna._unique_id) for dna in dnas]) + "]"
        building_history = "QUEEN.dna_dict['{}'] = joindna(*{}, topology='{}', project={}, product={}, process_description={})".format(construct._unique_id, dna_elements, topology, project, product if product is None else "'" + product + "'", process_description) 
        
        history_feature = _combine_history(construct, history_features)         
        construct._history_feature = history_feature 
        _add_history(construct, building_history) 

    for dnafeature in construct.dnafeatures:
        dnafeature.subject = construct
    
    if product is None:
        pass 
    else:
        _namespace[product] = construct
    return construct

def modifyends(dna, left="", right="", add=0, add_right=0, add_left=0, project=None, product=None, process_description=None, __direct=1):
    """
    Set end sequence structures. 
    """
    if process_description is None:
        process_description = QUEEN.process_description
    else:
        QUEEN.process_description = process_description 
    
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
        raise TypeError("The length of top and bottom strand sequence should be same.")
    
    elif "-" in left_end or "*" in left_end:
        left_end_top, left_end_bottom = check_endseq(left_end_top, left_end_bottom)
        if left_end_top != False:
            left_end_length = len(left_end_top)
            if "*" in left_end_top or "*" in left_end_bottom:
                if set(left_end_top) <= set(["*","-"]) and set(left_end_bottom) <= set(["*","-"]):
                    left_end_top    = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[:left_end_length], left_end_top)])
                    left_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[:left_end_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")), left_end_bottom)])
                else:
                    TypeError("'*' cannot be used wih characters for DNA sequence.")
            
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
            raise ValueError("The bottom strand sequence should be complent with the top strand sequence.")
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
        raise ValueError("The length of top and bottom sequence should be same.")
    
    elif "-" in right_end or "*" in right_end:
        right_end_top, right_end_bottom = check_endseq(right_end_top, right_end_bottom)
        if right_end_top != False:
            right_end_length = len(right_end_top)
            if "*" in right_end_top or "*" in right_end_bottom:
                if set(right_end_top) <= set(["*","-"]) and set(right_end_bottom) <= set(["*","-"]):
                    right_end_top    = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[-1*right_end_length:], right_end_top)])
                    right_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[-1*right_end_length:].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")), right_end_bottom)])
                else:
                    raise TypeError("'*' cannot be used wih characters for DNA sequence.")

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
            raise ValueError("The bottom strand sequence should be complent with the top strand sequence.")
    else:
        right_end_length = len(right_end_top) 
        right_length = len(right_end_top)  
        right_end_top    = 1
        right_end_bottom = 1
 
    if add == 1 or (left_end != dna.seq[left_end_length-left_length:left_end_length-left_length+len(left_end)] 
                or right_end != str(dna[len(dna.seq)-right_end_length + right_length - len(right_end):len(dna.seq)-right_end_length + right_length].seq)):
        
        if add_left == 1 and add_right == 1:
            new_dna = QUEEN(seq=left_end.split("/")[0] + dna.seq + right_end.split("/")[0] + "/" 
                    + left_end.split("/")[1] + dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + right_end.split("/")[1], _direct=0) 
            new_dna._history_feature = dna._history_feature 
            new_dna.dnafeatures      = _slide(dna.dnafeatures, len(left_end.split("/")[0])) 
            new_dna.record           = copy.copy(dna.record) 
            new_dna.record.features  = new_dna.dnafeatures
            new_dna._positions =  (-1,) * len(left_end.split("/")[0]) + new_dna._positions + (-1,) * len(right_end.split("/")[0])
            
        
        elif add_right == 1:
            left_end  = QUEEN(seq=left_end,  _direct=0) 
            right_end = QUEEN(seq=right_end, _direct=0)    
            new_dna   = cropdna(dna, start=left_end_length-left_length, end=len(dna.seq), __direct=0) + right_end
            new_dna._right_end       = right_end._right_end
            new_dna._left_end        = left_end.seq
            new_dna._left_end_top    = left_end_top 
            new_dna._left_end_bottom = left_end_bottom
            new_dna._positions = new_dna._positions + (-1,) * len(right_end.seq)
        
        else:
            left_end  = QUEEN(seq=left_end,  _direct=0) 
            right_end = QUEEN(seq=right_end, _direct=0)    
            new_dna = left_end + cropdna(dna, start=0, end=len(dna.seq)-right_end_length+right_length, __direct=0)
            new_dna._left_end         = left_end._left_end
            new_dna._right_end        = right_end.seq
            new_dna._right_end_top    = right_end_top 
            new_dna._right_end_bottom = right_end_bottom
            new_dna._positions = (-1,) * len(left_end.seq) + new_dna._positions

    else:
        new_dna   = copy.deepcopy(dna)  
        new_dna._left_end  = left_end
        new_dna._right_end = right_end
        new_dna._left_end_top     = left_end_top 
        new_dna._left_end_bottom  = left_end_bottom
        new_dna._right_end_top    = right_end_top
        new_dna._right_end_bottom = right_end_bottom
    
    if project is None:
        project = dna.project
    new_dna.project = project 

    if __direct == 1:
        if project in QUEEN.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(QUEEN.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in QUEEN.dna_dict:
                unique += 1
            new_dna._unique_id = project + "_" + str(unique)
        else:         
            new_dna._unique_id = project
        QUEEN.dna_dict[new_dna._unique_id] = product
        
        args = [] 
        history_features = [new_dna._history_feature] 
        args.append("'{}'".format(new_dna._unique_id))
        args.append("'{}'".format(dna._unique_id)) 
        if type(left_origin) == Qseq:
            if left_origin.parental_class == "DNAFeature":
                qkey = left_origin.qkey
                for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
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
                if left_origin.name != None and "getdnaseq" in left_origin.name:
                    seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, left_origin.name.split("_")[-1]) 
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
            
            elif left_origin.parental_class == "RE":
                args.append("RE_sites['{}'].{}".format(left_origin.parent.name, left_origin.name)) 

            else:
                args.append("'{}'".format(left_origin)) 
        else:
            args.append("'{}'".format(left_origin)) 
        
        if type(right_origin) == Qseq:
            if right_origin.parental_class == "DNAFeature":
                qkey = right_origin.qkey
                for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
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
                if right_origin.name != None and "getdnaseq" in right_origin.name:
                    seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, right_origin.name.split("_")[-1]) 
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
            
            elif right_origin.parental_class == "RE":
                args.append("RE_sites['{}'].{}".format(right_origin.parent.name, right_origin.name)) 

            else:
                args.append("'{}'".format(right_origin)) 
        else:
            args.append("'{}'".format(right_origin)) 

        if type(process_description) is str:
            process_description = "'" + process_description + "'"
        
        if type(project) is str:
            project =  "'" + project + "'" 
        
        args.append(project) 
        args.append(product if product is None else "'" + product + "'") 
        args.append(process_description) 
        building_history = "QUEEN.dna_dict[{}] = modifyends(QUEEN.dna_dict[{}], left={}, right={}, project={}, product={}, process_description={})".format(args[0], args[1], args[2], args[3], args[4], args[5], args[6])  
        history_feature = _combine_history(new_dna, history_features) 
        new_dna._history_feature = history_feature
        _add_history(new_dna, building_history) 
    
    for dnafeature in new_dna.dnafeatures:
        dnafeature.subject = new_dna
    
    if product is None:
        pass 
    else:
        _namespace[product] = new_dna
    return new_dna

def flipdna(dna, project=None, product=None, process_description=None, __direct=1):
    """
    Return reverse complement sequence. 
    All feature infomation is also reversed with the sequence.
    """
    if process_description is None:
        process_description = QUEEN.process_description
    else:
        QUEEN.process_description = process_description 

    if type(dna) is QUEEN:
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
            feats.append(DNAfeature(feature=feat,subject=seq))
    
    elif type(dna) is str:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        return seq

    elif type(dna) is Seq:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 

    comp = QUEEN(seq=seq) 
    feats.sort(key=lambda x:x.location.parts[0].start.position) 
    comp.dnafeatures = feats
    comp._setfeatureid()
    comp._features_dict = dict(list(map(lambda x:(x._id, x), comp.dnafeatures)))
    comp.record.features = comp.dnafeatures
    comp._right_end, comp._left_end = dna._left_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], dna._right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
    comp._right_end_top, comp._left_end_bottom = dna._left_end_bottom, dna._right_end_top
    comp._right_end_bottom, comp._left_end_top = dna._left_end_top, dna._right_end_bottom
    if project is None:
        project = dna.project
    comp.project = project

    if __direct == 1:
        if project in QUEEN.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(QUEEN.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in QUEEN.dna_dict:
                unique += 1
            comp._unique_id = project + "_" + str(unique)
        else:         
            comp._unique_id = project

        QUEEN.dna_dict[comp._unique_id] = product
        if type(project) is str:
            project =  "'" + project + "'" 
        if type(process_description) is str:
            process_description =  "'" + process_description + "'" 
        building_history = "QUEEN.dna_dict['{}'] = flipdna(QUEEN.dna_dict['{}'], project={}, product={}, process_description={})".format(comp._unique_id, dna._unique_id, project, product if product is None else "'" + product + "'", process_description) 
        _add_history(comp, building_history) 
   
    comp._positions = dna._positions[::-1] 
    for dnafeature in comp.dnafeatures:
        dnafeature.subject = comp
    
    if product is None:
        pass 
    else:
        _namespace[product] = comp
    return comp

def _get_matchlist_regex(dna, query, value=None, subject=None, s=None, e=None, strand=None):
    segment       = QUEEN(seq="") 
    match_list    = [] 
    
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
        for match in match_iter:
            result  = {"start":None, "end":None, "strand":None, "match":None} 
            span    = match.span() 
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
                    destination = flipdna(QUEEN(seq=match.expand(value))) 
                else:
                    destination = QUEEN(seq="")
                    for index, literal in enumerate(literals): #["A","A",None,"A"]  
                        index = len(literals) - index - 1
                        if literal is None:
                            sub_span = match.span(groups[index])
                            sub_span = (e-sub_span[1], e-sub_span[0])
                            destination = destination + dna[sub_span[0]:sub_span[1]]
                        else:
                            destination = destination + flipdna(QUEEN(seq=literal))
                
                if span[0] == pre_s:
                    segment = segment + destination
                else:
                    segment = segment + dna[pre_s:span[0]] + destination
                pre_s   = span[1]
        
        if mode == "edit" and pre_s != e:
            segment = segment + dna[pre_s:e]
    else:
        pre_s = 0  
        for match in match_iter:
            result  = {"start":None, "end":None, "strand":None, "match":None} 
            span    = match.span() 
            result["start"]  = s + span[0]
            result["end"]    = s + span[1]
            result["strand"] = strand
            result["match"]  = match
            match_list.append(result) 
            
            if mode == "edit":
                groups, literals = sre_parse.parse_template(value, query_pattern)
                groups = dict(groups)
                if len(groups) == 0:
                    destination = QUEEN(seq=match.expand(value))
                else:
                    destination = QUEEN(seq="")
                    for index, literal in enumerate(literals):
                        if literal is None:
                            sub_span    = match.span(groups[index])
                            destination = destination + dna[s+sub_span[0]:s+sub_span[1]]
                        else:
                            destination = destination + QUEEN(seq=literal)
                
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

def _compile_cutsite(query):  
    cutsite = query.upper() 
    re_format1 = re.compile(r"([ATGCRYKMSWBDHVN]+)(\([\-0-9]+/[\-0-9]+\))")
    re_format2 = re.compile(r"(\([\-0-9]+/[\-0-9]+\))([ATGCRYKMSWBDHVN]+)(\([\-0-9]+/[\-0-9]+\))")
    re_format3 = re.compile(r"(\([\-0-9]+/[\-0-9]+\))([ATGCRYKMSWBDHVN]+)")    
    if  match := re_format1.fullmatch(query):
        seq = match.group(1)
        topr, bottomr = map(int, match.group(2)[1:-1].split("/")) 
        topl, bottoml = "null", "null"
        
    elif match := re_format2.fullmatch(query):
        seq = match.group(2) 
        topl, bottoml = map(int, match.group(1)[1:-1].split("/")) 
        topr, bottomr = map(int, match.group(3)[1:-1].split("/")) 

    elif match := re_format3.fullmatch(query): 
        seq = match.group(2)
        topl, bottoml = map(int, match.group(1)[1:-1].split("/")) 
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

def editsequence(dna, source_sequence, destination_sequence, start=0, end=None, strand=1, project=None, product=None, process_description=None, __direct=1):
    dna    = copy.deepcopy(dna) 
    start  = 0 if start == len(dna.seq) else start
    end    = len(dna.seq) if end is None else end
    strand = 1 if strand is None else strand 
    if start == 0 and end == len(dna.seq):
        subject = dna.seq
    else:
        subject = dna.getdnaseq(start, end, strand)

    feat_list = [] 
    if source_sequence is None:
        segment = QUEEN(seq=re.sub(subject, value, subject))  
    
    else:
        if type(source_sequence) == QUEEN:
            source = source_sequence.seq.upper() 
        else:
            source = source_sequence.upper() 
        
        query = source 
        if strand == 1 or strand == -1:
            segment = _get_matchlist_regex(dna, query, value=destination_sequence, subject=subject, s=start, e=end, strand=strand) 
        else:
            ValueError("When edit the sequence, the sequence strand to be edit should be '-1' or '+1.'")
            
    if start == 0 and end == len(dna.seq):
        new_dna = segment
    elif start == 0:
        new_dna = segment + dna[e:len(dna.seq)] 
    elif end == len(dna.seq):
        new_dna = dna[0:s] + segment
    else:
        new_dna = dna[0:s] + segment + dna[e:len(dna.seq)] 
    
    #dna = new_dna
    dna.seq = new_dna.seq
    dna.record = new_dna.record
    dna.dnafeatures = new_dna.dnafeatures
   
    if __direct == 1:
        original_id = dna._unique_id
        if project is None:
            project = dna.project
        if project in QUEEN.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(QUEEN.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in QUEEN.dna_dict:
                unique += 1
            dna._unique_id = project + "_" + str(unique)
        else:         
            dna._unique_id = project
        QUEEN.dna_dict[dna._unique_id] = product 
        source_sequence      = "'" + source_sequence + "'" if source_sequence is not None else None
        destination_sequence = "'" + destination_sequence + "'" if destination_sequence is not None else None
        process_description  = "'" + process_description + "'" if process_description is not None else None
        building_history     = "QUEEN.dna_dict['{}'] = editsequence(QUEEN.dna_dict['{}'], query='{}', source_sequence={}, destination_sequence={}, start={}, end={}, strand={}, project={}, product={}, process_description={})".format(dna._unique_id, original_id, query, source_sequence, destination_sequence, start, end, strand, project, product if product is None else "'" + product + "'", process_description)
        _add_history(dna, history=building_history)  

    if product is None:
        pass 
    else:
        _namespace[product] = new_dna
    return new_dna


def _replaceattribute(dna=None, feat_list=None, target_attribute=None, query_re=None, value=None):    
    _exec =  0
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
                    segment = QUEEN(seq=re.sub(target_seq, "", target_seq))
                else:
                    segment = QUEEN(seq=re.sub(target_seq, value, target_seq))  
                if strand == -1:
                    segment = flipdna(segment, __direct==0)
            else:
                segment = _get_matchlist_regex(dna, query_re, value=value, subject=target_seq, s=s, e=e, strand=feat.strand)
            
            if s == 0 and e == len(dna.seq):
                dna = segment
            elif s == 0:
                dna = joindna(segment, dna[e:len(dna.seq)], topology="circular", __direct=0)
            elif e == len(dna.seq):
                dna = joindna(dna[0:s], segment, topology="circular", __direct=0) 
            else:
                dna = joindna(dna[0:s], segment, dna[e:len(dna.seq)], topology="circular", __direct=0) 
            
        for tmpfeat in dna.dnafeatures:
            if "_tmpid" in tmpfeat.__dict__:
                del tmpfeat._tmpid

        _exec += 1
        if "label" in feat.qualifiers:
            label = feat.qualifiers["label"][0]
        else:
            label = "N.A"
        label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, fs, fe))
        
        if set(dna.seq) < set("ATGCRYKMSWBDHVNatgcrykmsbdhvn-"):
            pass 
        else: 
            raise ValueError("Invalid sequence pattern was detected") 

    elif target_attribute == "start":
        for feat in feat_list:
            feat = dna._features_dict[feat._id] 
            feat.subject = dna 
            feat.set_position(value, "start") 
            _exec += 1 
    
    elif target_attribute == "end":
        for feat in feat_list:
            feat = dna._features_dict[feat._id] 
            feat.subject = dna 
            feat.set_position(value, "end") 
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
                                warnings.warn("Some features were detected. To ensure the uniqueness of the feature ID, a unique number will be added after the feature ID of each, scuh as '{}_[0-9]'".foramt(value)) 
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
        dna.dnafeatures = new_dnafeatures
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
            else:
                new_feat      = DNAfeature(feature=feat, subject=dna) 
                new_feat._id  = new_id 
                new_feat.type = feat_type
            new_dnafeatures.append(new_feat)

    dna.dnafeatures = list(dict([(feat._id, feat) for feat in new_dnafeatures]).values()) 
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeatures))) 
    return dna 

def createattribute(value=None):
    return functools.partial(_createattribute, value=value)

def _search(dna, source, query, attribute=None, strand=None): 
    #attributes = "feature ID", "feature type", "start", "end", "sequence", "query:" 
    re_digestion    = re.compile("[ATGCRYKMSWBDHVN]+[\^]?\([0-9]+/[0-9]+\)") 
    attribute_regex = re.compile("sequence:\|[0-9]+\.\.[0-9]+\|[+-]{0,1}")
    if query == ".+" and attribute == None:
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
    
    if attribute is None:
        if type(query) is str and (set(query) <= set("ATGCatgc")):
            attribute = "sequence"
        else:
            pass 

    if attribute == "start" and attribute == "end":
        raise TypeError("Positional value cannot be specified as 'key_attribute'. The 'start', 'end', and 'strand' argument are avaiable to define the sequence range for the search")         
    
    feat_list  = [] 
    if attribute is not None and attribute == "sequence":
        for feat in source:
            matches = dna.searchsequence(query, start=feat.start, end=feat.end, strand=feat.strand, _direct=0)
            if len(matches) > 0:
                feat_list.append(feat) 

    if attribute is None or attribute == "feature type" or attribute == "feature_type":
        if query is None:
            query = ".+"
        cquery = re.compile(query)
        for feat in source:
            if cquery.match(feat.type) != None or query == feat.type:
                feat_list.append(feat)  
                
    if attribute is None or attribute == "feature id" or attribute == "feature_id":
        if query is None:
            query = ".+"
        cquery       = re.compile(str(query))
        id_feat_list = [(feat.feature_id, feat) for feat in source]  
        for key, feat in id_feat_list:
            if cquery.match(key) != None or query == key:
                feat_list.append(copy.deepcopy(feat)) 

    if attribute is None or attribute == "strand":
        for feat in source:
            strand = feat.location.strand
            if query is None or strand == query or (query.isdecimal() and strand == query):
                feat_list.append(feat) 
    
    if attribute is None or attribute[0:len("qualifier:")] == "qualifier:":
        if query is None:
            query = ".+"
        cquery = re.compile(query)
        if attribute is None or attribute == "qualifier:*":
            for feat in source:
                flag = 0 
                for key in feat.qualifiers:
                    if type(feat.qualifiers[key]) is list:
                        pass 
                    else:
                        feat.qualifiers[key] = [feat.qualifiers[key]]
                    
                    for element in feat.qualifiers[key]:
                        if cquery.match(element) != None or query == element:
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
                        if cquery.match(element) != None or query == element:
                            feat_list.append(feat) 
                            break
    feat_set  = set([])  
    new_feat_list = []
    for feat in feat_list:
        element = (feat.location.start, feat.location.end, feat.sequence) 
        if element in feat_set:
            pass 
        else:
            new_feat_list.append(feat)
            feat_set.add(element)
    return new_feat_list

def editfeature(dna, key_attribute=None, query=".+", source=None, start=0, end=None, strand=2, target_attribute=None, operation=None, new_copy=True, project=None, product=None, process_description=None, __direct=1): 
    if process_description is None:
        process_description = QUEEN.process_description
    else:
        QUEEN.process_description = process_description

    if new_copy == False:
        pass 
    else:
        original_id = dna._unique_id
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
                item[1] = "'" + item[1] + "'"
            else:
                item[1] = str(item[1]) 
            largs.append("=".join(item)) 
        command = operation.func.__name__[1:] + "(" + ",".join(largs) + ")"
        
        if __direct == 1 and new_copy == True:
            if type(query) == Qseq:
                if query.parental_class == "DNAFeature":
                    qkey = left_origin.qkey
                    for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                        if qfeat._second_id == query.parental_id:
                            break
                    if type(query.item)   == int:
                        query = "QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , query.item)
                    elif type(query.item) == slice:
                        sl_start = query.item.start
                        sl_stop  = query.item.stop 
                        sl_step  = query.item.step
                        sl_start = "" if sl_start is None else sl_start
                        sl_stop  = "" if sl_stop is None else sl_stop
                        if sl_step == 1 or sl_step == None:
                            query = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop)
                        else:
                            query = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step)
                    else:
                        query = "QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex)

                elif query.parental_class == "QUEEN": 
                    parental_id = query.parental_id 
                    if query.name != None and "getdnaseq" in query.name:
                        seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, query.name.split("_")[-1]) 
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
                            query = "{}[{}:{}]".format(seqname, sl_start, sl_stop)
                        else:
                            query = "{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step)
                    else:
                        query = "{}".format(seqname)
                
                elif qorigin.parental_class == "RE":
                    qorigin = "RE_sites['{}'].{}".format(qorigin.parent.name, qorigin.name) 

                else:
                    query = "'{}'".format(query) 
    
            else:
                query = "'{}'".format(query) 
            
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
                project = dna.project

            if project in QUEEN.dna_dict:
                if project.split("_")[-1].isdecimal() == True:
                    project = "_".join(project.split("_")[:-1])
                keys   = list(QUEEN.dna_dict.keys())
                unique = 0
                while project + "_" + str(unique) in QUEEN.dna_dict:
                    unique += 1
                dna._unique_id = project + "_" + str(unique)
            else:         
                dna._unique_id = project
            QUEEN.dna_dict[dna._unique_id] = product
            
            if start == 0 and end == len(dna.seq):
                args = [key_attribute, query, source, strand, target_attribute, command, project, new_copy, process_description]
                for i in range(len(args)):
                    if type(args[i]) is str and i != 5:
                        args[i] = "'" + args[i] + "'" 
                building_history = "QUEEN.dna_dict['{}'] = editfeature(QUEEN.dna_dict['{}'], key_attribute={}, query={}, source={}, strand={}, target_attribute={}, operation={}, project={}, new_copy={}, product={}, process_description={})".format(dna._unique_id, original_id, *args[0:-1], product if product is None else "'" + product + "'", args[-1]) 
            
            else:
                args = [key_attribute, query, source, start, end, strand, target_attribute, command, project, new_copy, process_description]
                for i in range(len(args)):
                    if type(args[i]) is str and i != 7:
                        args[i] = "'" + args[i] + "'" 
                building_history = "QUEEN.dna_dict['{}'] = editfeature(QUEEN.dna_dict['{}'], key_attribute={}, query={}, source={}, start={}, end={}, strand={}, target_attribute={}, operation={}, project={}, new_copy={}, product={}, process_description={})".format(dna._unique_id, original_id, *args[0:-1], product if product is None else "'" + product + "'", args[-1]) 

            _add_history(dna, history=building_history) 
            if product is None:
                pass 
            else:
                _namespace[product] = dna

        elif new_copy == False:
            args = [key_attribute, query, source, start, end, strand, target_attribute, command, project, new_copy, process_description]
            for i in range(len(args)):
                if type(args[i]) is str and i != 7:
                    args[i] = "'" + a/rgs[i] + "'" 
            building_history = "editfeature(QUEEN.dna_dict['{}'], key_attribute={}, query={}, source={}, start={}, end={}, strand={}, target_attribute={}, operation={}, project={}, new_copy={}, product={}, process_description={})".format(dna._unique_id, *args[0:-1], product if product is None else "'" + product + "'", args[-1])
            _add_history(dna, history=building_history) 

    else:
        raise ValueError("The operational function can be selected from only 'createattribute', 'removeattribute', 'replaceattribute'.")

    if new_copy == True:
        return dna 

def _circularizedna(dna):
    dna = copy.deepcopy(dna)
    seq_origin = dna.seq
    feats_origin = dna.dnafeatures
    if dna.topology == "circular" and dna.record.annotations["topology"] == "circular":
        print("The DNA object is already circularized")

    if (dna._right_end_top * dna._left_end_bottom == 1 and dna._right_end_bottom * dna._left_end_top == 1) and len(dna._right_end) > 0 and (dna._left_end_top == -1 or dna._left_end_bottom == -1):
        if str(dna._right_end) == str(dna._left_end): 
            ovhg       = dna._right_end
            subdna     = cropdna(dna,0,len(dna.seq)-len(dna._right_end),__direct=0)
            dna.seq    = subdna.seq
            dna.record = subdna.record
        else:
            return False
    else:
        ovhg = ""

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
                    if s1 > e2 and length1 == length2 and "_original" in feat1.__dict__ and "_original" in feat2.__dict__ and feat1.original == feat2.original and feat1.location.strand == feat2.location.strand:
                        #print(feat1)
                        #print(feat2) 
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
                        
                        new_feat  = DNAfeature(feature=new_feat, subject=dna)
                        new_feat1 = DNAfeature(feature=feat1, subject=dna)
                        new_feat2 = DNAfeature(feature=feat2, subject=dna) 
                        s = new_feat.start
                        e = new_feat.end                                    
                        #if (len(new_seq) - len(ovhg) == e - s and len(new_seq) - len(ovhg) <= len(feat1.original)) or 
                        if new_feat._original == dna.getdnaseq(new_feat1.start, new_feat2.end, new_feat.location.strand if new_feat.location.strand !=0 else 1):
                            dna.dnafeatures[feat1_index].qualifiers["broken_feature"] = [note]
                            if len(new_seq) - len(ovhg) == length1:
                                del dna.dnafeatures[dna.dnafeatures.index(feat1)].qualifiers["broken_feature"]
                            dna.dnafeatures[feat1_index].location = new_feat.location
                            dna.dnafeatures.remove(feat2) 
                            remove_list.append(feat2) 
    
    for i in range(len(dna.dnafeatures)):    
        if dna.dnafeatures[i].location.parts[-1].end.position > len(dna.seq):
            if dna.dnafeatures[i].location.parts[0].start.position >= len(dna.seq):
                strand                      = dna.dnafeatures[i].location.strand
                dna.dnafeatures[i].location = FeatureLocation(dna.dnafeatures[i].location.parts[0].start.position-len(dna.seq),dna.dnafeatures[i].location.parts[-1].end.position-len(dna.seq))
                dna.dnafeatures[i].location.strand = strand
            else:
                strand    = dna.dnafeatures[i].location.strand
                locations = [FeatureLocation(dna.dnafeatures[i].location.parts[0].start.position,len(dna.seq)), FeatureLocation(0,dna.dnafeatures[i].location.parts[-1].end.position-len(dna.seq))]
                if strand == -1:
                    locations.reverse()   
                dna.dnafeatures[i].location = CompoundLocation(locations)
                dna.dnafeatures[i].location.strand = strand
            dna.dnafeatures[i] = DNAfeature(feature=dna.dnafeatures[i], location=dna.dnafeatures[i].location)
    
    #dna.record.features = _assigndnafeatures(dna.dnafeatures)
    dna._left_end  = ""
    dna._left_end_top    = 0 
    dna._left_end_bottom = 0 

    dna._right_end = ""
    dna._right_end_top    = 0 
    dna._right_end_bottom = 0 

    dna.topology = "circular"
    if Alphabet:
        dna.record.seq = Seq(str(dna.seq),Alphabet.DNAAlphabet())
    else:
        dna.record.seq = Seq(str(dna.seq))
    dna.record.annotations["topology"] = dna.topology
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))
    return dna

def visualize(dna, map_view="linear", feature_list=None, start=0, end=None, standard_scale="auto", width_scale=1.0, height_scale=1.0, label_location=None, linebreak=None, seq=False, label_box=True, fontsize=None, diamater_scale=1.0, view_title=True, view_axis=True, tick_space="auto"):
    if fontsize is None and map_view == "linear":
        fontsize = 12
    elif fontsize is None and map_view == "circular":
        fontsize = 10
    else:
        pass 

    if feature_list is None:
        feature_list = dna.dnafeatures
    if map_view == "circular":
        fig, axes = vc.visualize(dna, format=0, feature_list=feature_list, unvisible_types=["source"], visible_types=[], bottom=400 * diamater_scale, label_box=label_box, fontsize=fontsize, 
        view_title=view_title, view_axis=view_axis, tick_space=tick_space)
    else:
        fig, ax = vl.visualize(dna, start=start, end=end, feature_list=feature_list, wrap_width=linebreak, annotation_loc=label_location, unvisible_types=["source"], visible_types=[], enlarge_w=width_scale, enlarge_h=height_scale, 
        fontsize=fontsize, with_seq=seq, nucl_char=None, nucl_color_dict=None, label_box=label_box, scale=standard_scale, view_title=view_title, view_axis=view_axis, tick_space=tick_space)
    return fig

def _assigndnafeatures(dnafeatures):
    features = [] 
    for feat in dnafeatures:
        if feat.location.start.position == -1:
            pass 
        else:
            features.append(feat) 
    return features 

def quine(dna, output=None, description_only=False):
    def export(descriptions, histories, o=None, do=False):
        num = -1
        for h, (process_description, history) in enumerate(zip(descriptions, histories)):
            if do == True:
                if process_description is None or str(process_description) == str(pre_process_description):
                    pass
                else:
                    print(process_description[1:-1], file=o) 
            else:
                if process_description is None or str(process_description) == str(pre_process_description):
                    pass
                else:
                    if h > 0:
                        print("", file=o)
                    num += 1
                    print("description{} = ".format(num) + "'" + process_description[1:-1] + "'", file=o) 

                if "QUEEN.queried_feature_dict" in history[1][0:len("QUEEN.queried_feature_dict")]:
                    print(history[1], file=o)
                elif process_description == None or num == -1:
                    print(history[1][:-1] + ", process_description={})".format("None"), file=o)
                else:
                    print(history[1][:-1] + ", process_description=description{})".format(num), file=o)
            pre_process_description = process_description 
        return o 

    histories  = dna.history   
    hindex = 0
    for index, history in enumerate(histories):
        if "QUEEN.dna_dict" in history[1].replace(" ","").replace("–"," "):
            hindex = index
        else:
            pass 
    result     = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]",histories[hindex][1].replace(" ","").replace("–"," "))[0] 
    _unique_id = result.split("['")[1][:-2]
    
    pre_pd = None
    descriptions  = [] 
    new_histories = [] 
    for history in histories:
        history = list(history)
        if re.search("process_description=None",history[1].replace(" ","").replace("–"," ")) is None:
            process_description = re.search("process_description='.*'",history[1].replace(" ","").replace("–"," "))
            if process_description is None:
                pd = pre_pd                
            else: 
                process_description = history[1].replace(" ","").replace("–"," ")[process_description.start():process_description.end()] 
                pd = process_description.split("=")[1]    
                if pd == "''" or pd == '""':
                    pd = None
            descriptions.append(pd)
        else:
            process_description = "process_description=None"
            pd = None
            descriptions.append(pd)
        
        if process_description is not None:
            history[1] = history[1].replace(" ","").replace("–"," ").replace(process_description,"")
            history[1] = history[1].replace(", )",")") 
        else:
            history[1] = history[1].replace(" ","").replace("–"," ")
        new_histories.append(history) 
        pre_pd = pd
    histories = new_histories

    #Remove nonu-sed variable 
    outtext = export(descriptions, histories, io.StringIO())
    text    = outtext.getvalue().rstrip()
    var_num_dict = collections.defaultdict(int) 
    for row in text.split("\n"):
        row = row.rstrip()
        matches = re.findall("QUEEN.queried_feature_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_num_dict[match] += 1
        
        matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_num_dict[match] += 1
    
    new_descriptions = [] 
    new_histories    = []
    texts = [row for row in text.split("\n") if row != "" and row[0:len("description")] != "description"] 
    for row, description, history in zip(texts, descriptions, histories):
        row      = row.rstrip() 
        var_nums = [] 
        matches = re.findall("QUEEN.queried_feature_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_nums.append(var_num_dict[match]) 
        
        matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_nums.append(var_num_dict[match])
        
        if len(var_nums) > 0:
            if max(var_nums) == 1:
                pass 
            else:
                new_descriptions.append(description) 
                new_histories.append(history) 
        else:
            new_descriptions.append(description) 
            new_histories.append(history) 

    if type(output) is str:
        o = open(output, "w") 
    elif output is None:
        o = None
   
    pre_process_description = "''"
    if description_only == False:
        print("import sys", file=o)  
        print("sys.path.append(\"{}".format("/".join(__file__.split("/")[:-2])  + "\")"), file=o)
        print("from QUEEN.queen import *", file=o)  
        if _namespaceflag == 1:
            print("set_namespace(globals())", file=o)   
    
    outtext = export(new_descriptions, new_histories, io.StringIO(), do=description_only)
    outtext = outtext.getvalue().rstrip()
    texts   = outtext.split("\n") 
    if description_only == False:
        texts.append(result + ".writedna()".format(dna.project))

    name_dict = {}
    for row in texts:
        match1 = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\] = ",row)
        match2 = re.findall("product='.+',",row)
        if len(match1) > 0 and len(match2) > 0:
            name_dict[match1[0].split("=")[0][:-1]] = match2[0][:-1].split("=")[1][1:-1] 

        match1 = re.findall("QUEEN.dna_dict\['[^\[\]]+'\] = ",row)
        match2 = re.findall("product='.+',",row)
        if len(match1) > 0 and len(match2) > 0:
            name_dict[match1[0].split("=")[0][:-1]] = match2[0][:-1].split("=")[1][1:-1] 
    
    new_rows = []
    for row in texts:
        if _namespaceflag == 1:
            product = re.match("product=(.+),",row)
            if product is None:
                pass 
            elif product.group(1) == "None":
                pass 
            elif product is not None:
                matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\] = ",row)
                if len(matches) > 0:
                    for match in matches:
                        row  = row.replace(match, "")
            
                matches = re.findall("QUEEN.dna_dict\['[^\[\]]+'\] = ",row)
                if len(matches) > 0:
                    for match in matches:
                        row  = row.replace(match, "")
            else: 
                pass 

        matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
        if len(matches) > 0:
            for match in matches:
                name = match
                if name in name_dict:
                    row = row.replace(match, name_dict[name]) 
        
        matches = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]",row)
        if len(matches) > 0:
            for match in matches:
                name = match
                if name in name_dict:
                    row = row.replace(match, name_dict[name])
        new_rows.append(row) 
    
    for row in new_rows:
        print(row, file=o) 

    if output is not None:
        o.close() 
        
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
    for feat in dna.dnafeatures:
        if feat.type == "source":
            old_keys = [] 
            for key in feat.qualifiers:
                if "building_history" in key:
                    old_keys.append(key)   
            for key in old_keys:
                del feat.qualifiers[key] 

class QUEEN():    
    #Class variables that manage the execution histories of operational and search function
    dna_dict                   = {}
    queried_feature_dict       = {}
    queried_features_dict      = {}
    queried_features_name_dict = {} 
    process_description = None
    _num_history = 1  
    _qnum = 0
    
    def _check_seq(seq):
        top, bottom = seq.split("/") 
        new_seq     = ""
        new_top     = ""
        new_bottom  = ""
        if len(top) != len(bottom):
            return False, False, False

        for t,b in zip(top,bottom):
            if t != b.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) and (t != "-" and b != "-"):
                return False, False, False
            else:
                new_top    += t
                new_bottom += b
                if t == "-":
                    new_seq += b.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
                else:
                    new_seq += t
        return new_top, new_bottom, new_seq 
    
    def __deepcopy__(self, memo):
        obj = QUEEN(seq=self.seq, _direct=0)
        for key in self.__dict__:
            if key == "seq":
                pass 
            elif key == "_history_feature":
                obj._history_feature = DNAfeature(self._history_feature, subject=obj)     
            elif key == "dnafeatures": 
                feats = [] 
                for feat in self.dnafeatures:
                    feats.append(DNAfeature(feat, subject=obj))
                obj.dnafeatures = feats
            else:
                obj.__dict__[key] = copy.copy(self.__dict__[key])
        return obj

    def __repr__(self):
        if len(self.seq) > 50:
            out = "<queen.QUEEN object; project='{}', length='{} bp', topology='{}' >".format(self.project, len(self.seq), self.topology)
        else:
            out = "<queen.QUEEN object; project='{}', length='{} bp', sequence='{}', topology='{}'>".format(self.project, len(self.seq), self.seq, self.topology)
        return out 
    
    def __getattribute__(self, name):
        if name == "seq":
            qseq = Qseq(super().__getattribute__(name))
            if "_unique_id" in self.__dict__:
                qseq.parental_id = self._unique_id 
            qseq.parent = self 
            qseq.parental_class = "QUEEN"
            return qseq
        else:
            return super().__getattribute__(name)  

    def __getattr__(self, name): 
        if name == "history":
            histories = [] 
            for key in self._history_feature.qualifiers:
                if "building_history" in key[0:18]: 
                    history = self._history_feature.qualifiers[key][0] 
                    histories.append((int(key.split("_")[-1]), history)) 
            histories.sort()
            return histories

        if name == "sequence":
            self.sequence = self.seq
        else:
            raise AttributeError("Queen obejct has no attribute '{}'".format(name))

    def __init__(self, seq=None, record=None, project=None, topology="linear", format=None, product=None, process_description=None, import_history=True, _direct=1):
        if process_description is None:
            process_description = "" 
        else:
            QUEEN.process_description = process_description 

        recname = record
        self.seq                = None
        self.record             = None
        self.dnafeatures        = None
        self._right_end         = None 
        self._left_end          = None
        self._history_feature   = None
        if seq is None and record is None:
            if project is None: 
                project = "dna"
            self.topology           = topology
            self.project            = project
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1
        
        elif seq is None or "." in seq:
            if "." in str(seq):
                record = seq
            
            if type(record) == str:
                if format != None:
                    fmt = format
                else:
                    if record.split(".")[-1] in ["gb","gbk","genbank"]:
                        fmt = "genbank"
                    elif record.split(".")[-1] in ["fasta","fna","fa","faa"]:
                        fmt = "fasta"
                    else:
                        fmt = "genbank"
                record  = SeqIO.parse(record,fmt)
                record  = next(record) 

            elif type(record) == SeqRecord:
                record = record 
            else:
                record  = SeqIO.parse(record,None)

            self.seq    = str(record.seq).upper()
            self.record = record            
            if "topology" in record.annotations:
                self.topology = record.annotations["topology"]
            else:
                self.topology = topology
            
            if self.topology == "linear":
                self._right_end = self.seq[0:10] 
                self._left_end  = self.seq[-10:]
            else:
                self._right_end = ""
                self._left_end  = ""
            
            if project is None:
                if record.id == "" or record.id == ".":
                    project = recname.split("/")[-1].split(".")
                    project = project if len(project) == 1 else ".".join(project[:-1])
                else:
                    project = record.id 
            self.project = project
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1
            
            #import features
            self.dnafeatures = [] 
            if len(record.features) > 0:
                for feat in record.features:
                    self.dnafeatures.append(DNAfeature(feature=feat, subject=self))
               
                pairs = [] 
                history_feature = None
                history_nums = [QUEEN._num_history] 
                for feat in self.dnafeatures:
                    if feat.type == "source" and feat.start == 0 and feat.end == len(self.seq):
                        for key in feat.qualifiers:
                            if "building_history" in key[0:18] and import_history == True:
                                history = feat.qualifiers[key][0]
                                results = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]", history) 
                                for result in results:
                                    _unique_id = result.split("['")[1][:-2] 
                                    QUEEN.dna_dict[_unique_id] = None
                                history_num = int(key.split("_")[-1]) 
                                pairs.append((feat, history_num, history))
                                history_feature_id = "0" 
                                history_feature    = feat  
                    else:
                        for key in feat.qualifiers:
                            if key == "broken_feature":
                                note     = feat.qualifiers["broken_feature"][0]
                                original = note.split(":")[-3]
                                feat._original = original
               
                if len(pairs) == 0:
                    import_history = False

                if import_history == True:
                    for pair in pairs:
                        feat = pair[0]
                        new_history_num = pair[1] + QUEEN._num_history
                        feat.qualifiers["building_history_{}".format(new_history_num)] = [pair[2]] 
                        del feat.qualifiers["building_history_{}".format(pair[1])]
                        history_nums.append(new_history_num)      
                else:
                    QUEEN.process_description = process_description
                    archivehistory(self)    
                
                QUEEN._num_history = max(history_nums) 
                
                if history_feature is not None:
                    self._history_feature = history_feature
                    self.dnafeatures.remove(history_feature) 
            
            if len(self.dnafeatures) == 0:
                import_history = False
                self.dnafeatures = []

        elif record is None:
            import_history = False
            if project is None:
                project = "dna"

            sticky   = False
            pattern1 = "[ATGCRYKMSWBDHVNatgcrykmsbdhvn*-]+/?[ATGCRYKMSWBDHVNatgcrykmsbdhvn*-]*"
            pattern2 = "[ATGCRYKMSWBDHVNatgcrykmsbdhvn]+-+[ATGCRYKMSWBDHVNatgcrykmbdhvn]+"
            pattern1 = re.compile(pattern1)   
            pattern2 = re.compile(pattern2)   
            if pattern1.fullmatch(seq) != None or seq == "": 
                if "/" in seq:
                    top, bottom, seq = QUEEN._check_seq(seq)
                    if "-" in top or "-" in bottom:
                        if top !=False and pattern2.search(top) is None and pattern2.search(bottom) is None:
                            sticky = True
                        else:
                            raise TypeError("Invaild sequence pattern was detected.")
                    else: 
                        sticky = False

                self.seq        = str(seq).upper()
                self.project    = project
                if Alphabet:
                    self.record = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
                else:
                    self.record = SeqRecord(Seq(str(seq)))
                
                self.dnafeatures = []
                self.topology   = topology
                
                if sticky == True:
                    self.topology   = "linear"
                    self._left_end  = "" 
                    self._right_end = "" 
                    if top[0] == "-":
                        for c, char in enumerate(top):
                            if char != "-":
                                break 
                            else:
                                self._left_end += bottom[c].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
                        self._left_end_top    = -1
                        self._left_end_bottom = 1
                    
                    else:
                        if bottom[0] == "-":
                            for c, char in enumerate(bottom):
                                if char != "-":
                                    break 
                                else:
                                    self._left_end += top[c]
                            self._left_end_bottom = -1
                            self._left_end_top    = 1
                        
                        else:
                            self._left_end  = ""
                            self._left_end_bottom = 1
                            self._left_end_top    = 1

                    if top[-1] == "-":
                        for c, char in enumerate(top[::-1]):
                            if char != "-":
                                break 
                            else:
                                self._right_end += bottom[::-1][c].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
                        self._right_end = self._right_end[::-1]
                        self._right_end_top    = -1
                        self._right_end_bottom = 1

                    else:
                        if bottom[-1] ==  "-":
                            for c, char in enumerate(bottom[::-1]):
                                if char != "-":
                                    break 
                                else:
                                    self._right_end += top[::-1][c]
                            self._right_end = self._right_end[::-1]
                            self._right_end_bottom = -1
                            self._right_end_top    = 1
                        
                        else:
                            self._right_end  = ""
                            self._right_end_bottom = 1
                            self._right_end_top    = 1
                
                else:
                    self._right_end = ""
                    self._left_end  = ""
                    self._left_end_top      = 1
                    self._left_end_bottom   = 1
                    self._right_end_top     = 1 
                    self._right_end_bottom  = 1

            else:
                raise TypeError("Invaild sequence pattern was detected.")
       
        self._setfeatureid()
        self._features_dict = dict(list(map(lambda x:(x._id, x), self.dnafeatures)))
        
        if _direct == 1 and import_history == False:
            if project in QUEEN.dna_dict:
                keys   = list(QUEEN.dna_dict.keys())
                unique = 0
                while project + "_" + str(unique) in QUEEN.dna_dict:
                    unique += 1
                self._unique_id = project + "_" + str(unique)
            else: 
                self._unique_id = project 
        
            QUEEN.dna_dict[self._unique_id] = product
            args = [seq, recname, project, topology, format, process_description]
            for i in range(len(args)):
                if type(args[i]) is str:
                    args[i] = "'" + args[i] + "'" 

            building_history = "QUEEN.dna_dict['{}'] = QUEEN(seq={}, record={}, project={}, topology={}, format={}, product={}, process_description={})".format(self._unique_id, args[0], args[1], args[2], args[3], args[4], product if product is None else "'" + product + "'", args[5]) 
            _add_history(self, building_history)
        else:
            self._unique_id = project 
        

        self._positions       = tuple(range(len(self.seq))) 
        self.record.feartures = self.dnafeatures
        self.seq.parental_id  = self._unique_id
        if product is None:
            pass 
        else:
            _namespace[product] = self

    def searchsequence(self, query, start=0, end=None, strand=2, product=None, process_description=None, _direct=1):
        qorigin = query
        start   = 0 if start == len(self.seq) else start
        end     = len(self.seq) if end is None else end
        strand  = 2 if strand is None else strand 
        if start == 0 and end == len(self.seq):
            subject = self.seq
        else:
            subject = self.getdnaseq(start, end, strand=1)

        feat_list = [] 
        if query is None:
            if start > end:
                locations = [(start, len(self.seq), strand), (0, end, strand)] 
                if strand == -1:
                    locations.reverse()
                new_feat = SeqFeature(CompoundLocation(list(map(FeatureLocation, locations))), type="misc_feature")
            else:
                new_feat = SeqFeature(FeatureLocation(start, end, strand=1), type="misc_feature")
            new_feat = DNAfeature(feature=new_feat, subject=dna, query=subject)
            feat_list.append(new_feat) 
        
        else:
            cut = 0 
            if set(str(query)) <= set("ATGCRYKMSWBDHVN^_/()-0123456789") and len(set(str(query)) & set("^_/()-0123456789")) > 0 and query.count("^") <= 2 and query.count("_") <= 2 and query.count("^") == query.count("_"):
                cut = 1
                cutsite, query, topl, topr, bottoml, bottomr = _compile_cutsite(query)  
            else:
                pass
            
            if len(set(str(query)) & set("RYKMSWBDHVN()0123456789")) > 0:
                query = query.replace("R","[RGA]")
                query = query.replace("Y","[YTC]") 
                query = query.replace("K","[KGT]") 
                query = query.replace("M","[MAC]") 
                query = query.replace("S","[SGC]") 
                query = query.replace("W","[WAT]")
                query = query.replace("B","[BGTC]")
                query = query.replace("D","[DGAT]")
                query = query.replace("H","[HACT]")
                query = query.replace("V","[VGCA]")
                query = query.replace("N","[NAGCT]")
            
            if strand == 2:
                #_get_matchlist_regex(dna, query, value=None, subject, s, e, strand)
                match_list = _get_matchlist_regex(self, query, value=None, subject=subject, s=start, e=end, strand=1) 
                match_list.extend(_get_matchlist_regex(self, query, value=None, subject=subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], s=start, e=end, strand=-1))
            elif strand == 1:
                match_list = _get_matchlist_regex(self, query, value=None, subject=subject, s=start, e=end, strand=strand) 
            elif strand == -1:
                match_list = _get_matchlist_regex(self, query, value=None, subject=subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], s=start, e=end, strand=strand)
            else:
                ValueError("When edit the sequence, the sequence strand to be edit should be '-1' or '+1.'")
            
            match_positions = set() 
            for match in match_list:
                if (span := (match["start"], match["end"]))  not in match_positions:
                    match_positions.add(span) 
                    if match["start"] > match["end"] and self.topology == "circular":
                        locations = [[match["start"], len(self.seq), match["strand"]], [0, match["end"], match["strand"]]]
                        if match.strand == -1:
                            locations.reverse()
                            new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type="misc_feature") 
                        else:
                            new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type="misc_feature")
                    else:
                        new_feat = SeqFeature(FeatureLocation(match["start"], match["end"], strand=match["strand"]), type="misc_feature")

                    new_feat = DNAfeature(feature=new_feat, subject=self, query=qorigin.replace("^","").replace("_","")) 
                    if cut == 1:
                        new_feat._digestion_topl        = topl 
                        new_feat._digestion_topr        = topr
                        new_feat._digestion_bottoml     = bottoml
                        new_feat._digestion_bottomr     = bottomr
                        new_feat.qualifiers["cut_site"] = [Qseq(cutsite)]
                    new_feat.qualifiers["note_searchseqeunce"] = ["query:{}".format(qorigin)]
                    if type(qorigin) == Qseq and qorigin.parental_class == "RE":
                        new_feat.qualifiers["label"] = [qorigin.parent.name]
                    feat_list.append(new_feat)  
                else:
                    pass 
                    #print(match) 
                    
        qkey = self._unique_id + "_" + str(QUEEN._qnum)
        for i, feature in enumerate(feat_list):
            if "label" in feature.qualifiers:
                label = "_" + feature.qualifiers["label"][0]
            else:
                label = "" 
            
            if "_id" in feature.__dict__:
                key = qkey + "_" + feature._id + label
            else:
                key = qkey + "_" + "q" + str(i) + label
            
            feature._qkey              = qkey
            feature._second_id         = key
            feature._start.qkey        = qkey
            feature._start.parental_id = key
            feature._end.qkey          = qkey
            feature._end.parental_id   = key
        
        if "_unique_id" not in self.__dict__:
            self._unique_id = self.project

        QUEEN.queried_features_dict[qkey]      = feat_list
        QUEEN.queried_features_name_dict[qkey] = product
        
        if _direct == 1:
            if type(qorigin) == Qseq:
                if qorigin.parental_class == "DNAFeature":
                    qkey = qorigin.qkey
                    for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                        if qfeat._second_id == qorigin.parental_id:
                            break
                    if type(qorigin.item)   == int:
                        qorigin = "QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , qorigin.item)
                    elif type(qorigin.item) == slice:
                        sl_start = qorigin.item.start
                        sl_stop  = qorigin.item.stop 
                        sl_step  = qorigin.item.step
                        sl_start = "" if sl_start is None else sl_start
                        sl_stop  = "" if sl_stop is None else sl_stop
                        if sl_step == 1 or sl_step == None:
                            qorigin = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop)
                        else:
                            qorigin = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step)
                    else:
                        qorigin = "QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex)

                elif qorigin.parental_class == "QUEEN": 
                    parental_id = qorigin.parental_id 
                    if qorigin.name != None and "getdnaseq" in qorigin.name:
                        seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, qorigin.name.split("_")[-1]) 
                    else:
                        seqname = "QUEEN.dna_dict['{}'].seq".format(parental_id)
                    if type(qorigin.item)   == int:
                        args.append("{}[{}]".format(seqname, qorigin.item))
                    elif type(qorigin.item) == slice:
                        sl_start = qorigin.item.start
                        sl_stop  = qorigin.item.stop 
                        sl_step  = qorigin.item.step
                        sl_start = "" if sl_start is None else sl_start
                        sl_stop  = "" if sl_stop is None else sl_stop
                        if sl_step == 1 or sl_step == None:
                            qorigin = "{}[{}:{}]".format(seqname, sl_start, sl_stop)
                        else:
                            qorigin = "{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step)
                    else:
                        qorigin = "{}".format(seqname)
                
                elif qorigin.parental_class == "RE":
                    qorigin = "RE_sites['{}'].{}".format(qorigin.parent.name, qorigin.name) 
                else:
                    qorigin = "'{}'".format(qorigin) 
            else:
                qorigin = "'{}'".format(qorigin) 
            if start == 0 and end == len(self.seq):
                if strand == 2:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, product={}, process_description='{}')".format(qkey, self._unique_id, qorigin, product if product is None else "'" + product + "'", process_description)
                else:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, strand={}, product={}, process_description='{}')".format(qkey, self._unique_id, qorigin, strand, product if product is None else "'" + product + "'", process_description)
            else:
                building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, start={}, end={}, strand={}, product={}, process_description='{}')".format(qkey, self._unique_id, qorigin, start, end, strand, product if product is None else "'" + product + "'", process_description)   
            _add_history(self, building_history)
            QUEEN._qnum += 1 
        
        if product is None:
            pass 
        else:
            _namespace[product] = feat_list
  
        return feat_list

    def searchfeature(self, key_attribute=None, query=".+", source=None, start=0, end=None, strand=2, product=None, process_description=None, _direct=1):
        if process_description is None:
            process_description = QUEEN.process_description
        else:
            QUEEN.process_description = process_description  
        
        qkey = self._unique_id + "_" + str(QUEEN._qnum)
        features = editfeature(self, key_attribute=key_attribute, query=query, source=source, start=start, end=end, strand=strand, target_attribute=None, operation=None, __direct=0, process_description=process_description) 
        feature_names = [] 
        for i, feature in enumerate(features):
            if "label" in feature.qualifiers:
                label = "_" + feature.qualifiers["label"][0]
            else:
                label = "" 
            
            if "_id" in feature.__dict__:
                key = qkey + "_" + feature._id + label
            else:
                key = qkey + "_" + "q" + str(i) + label
            
            feature._qkey              = qkey
            feature._second_id         = key
            feature._start.qkey        = qkey
            feature._start.parental_id = key
            feature._end.qkey          = qkey
            feature._end.parental_id   = key
        
        if "_unique_id" not in self.__dict__:
            self._unique_id = self.project

        feature_names = ", ".join(feature_names)
        QUEEN.queried_features_dict[qkey]      = features
        QUEEN.queried_features_name_dict[qkey] = product
        
        if _direct == 1:
            if type(query) == Qseq:
                if query.parental_class == "DNAFeature":
                    qkey = left_origin.qkey
                    for qindex, qfeat in enumerate(QUEEN.queried_features_dict[qkey]):
                        if qfeat._second_id == query.parental_id:
                            break
                    if type(query.item)   == int:
                        query = "QUEEN.queried_features_dict['{}'][{}].{}[{}]".format(qkey, qindex, "seq" , query.item)
                    elif type(query.item) == slice:
                        sl_start = query.item.start
                        sl_stop  = query.item.stop 
                        sl_step  = query.item.step
                        sl_start = "" if sl_start is None else sl_start
                        sl_stop  = "" if sl_stop is None else sl_stop
                        if sl_step == 1 or sl_step == None:
                            query = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}]".format(qkey, qindex, sl_start, sl_stop)
                        else:
                            query = "QUEEN.queried_features_dict['{}'][{}].seq[{}:{}:{}]".format(qkey, qindex, sl_start, sl_stop, sl_step)
                    else:
                        query = "QUEEN.queried_features_dict['{}'][{}].seq".format(qkey, qindex)

                elif query.parental_class == "QUEEN": 
                    parental_id = query.parental_id 
                    if query.name != None and "getdnaseq" in query.name:
                        seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, query.name.split("_")[-1]) 
                    else:
                        seqname = "QUEEN.dna_dict['{}'].seq".format(parental_id)
                    
                    if type(query.item) == int:
                        args.append("{}[{}]".format(seqname, query.item))
                    elif type(query.item) == slice:
                        sl_start = query.item.start
                        sl_stop  = query.item.stop 
                        sl_step  = query.item.step
                        sl_start = "" if sl_start is None else sl_start
                        sl_stop  = "" if sl_stop is None else sl_stop
                        if sl_step == 1 or sl_step == None:
                            query = "{}[{}:{}]".format(seqname, sl_start, sl_stop)
                        else:
                            query = "{}[{}:{}:{}]".format(seqname, sl_start, sl_stop, sl_step)
                    else:
                        query = "{}".format(seqname)
                
                elif qorigin.parental_class == "RE":
                    qorigin = "RE_sites['{}'].{}".format(qorigin.parent.name, qorigin.name) 
                else:
                    query = "'{}'".format(query) 
            else:
                query = "'{}'".format(query) 
            
            if source is not None:
                qkeys = set([]) 
                for feat in source:
                    if "_qkey" in feat.__dict__:
                        qkeys.add(feat._qkey)
                
                if len(set(qkeys)) == 1:
                    source = "QUEEN.queried_features_dict['{}']".format(list(qkeys)[0])
                else:
                    pass 
            if start == 0 and end == len(self.seq):
                if strand == 2 and source is None:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query='{}', product={}, process_description='{}')".format(qkey, self._unique_id, key_attribute, query, product if product is None else "'" + product + "'", process_description)
                elif strand == 2:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query='{}', source={}, product={}, process_description='{}')".format(qkey, self._unique_id, key_attribute, query, source, product if product is None else "'" + product + "'", process_description)
                else:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query='{}', source={}, strand={}, product={}, process_description='{}')".format(qkey, self._unique_id, key_attribute, query, source, strand, product if product is None else "'" + product + "'", process_description)
            else:
                building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query='{}', source={}, start={}, end={}, strand={}, product={}, process_description='{}')".format(qkey, self._unique_id, key_attribute, query, source, start, end, strand, product if product is None else "'" + product + "'", process_description) 
            _add_history(self, building_history)
            QUEEN._qnum += 1 
        
        if product is None:
            pass 
        else:
            _namespace[product] = features

        return features

    def __getitem__(self, item):
        if type(item) == slice:
            if item.step is None:
                strand = 1
            elif type(item.step) != int:
                raise TypeError("slice indices must be integers or None or have an __index__ method.")
            else:
                strand = item.step

            if item.start is None:
                start = 0 
            else:
                if type(item.start) == int:
                    start = item.start 
                else:
                    raise TypeError("slice indices must be integers or None or have an __index__ method.")
           
            if item.stop is None:
                end = len(self.seq) 
            else:
                if type(item.stop) == int:
                    end = item.stop
                
                else:
                    raise TypeError("slice indices must be integers or None or have an __index__ method.") 
            
            if type(start) == int and type(end) == int:
                if start < 0:
                    start = len(self.seq) - abs(start) 
                if end < 0: 
                    end = len(self.seq) - abs(end) 
                subdna = cropdna(self,start,end,__direct=0)
                             
            else:
                raise TypeError("slice indices must be integers or None or have an __index__ method")
            
            if strand == -1 or strand < 0:
                return flipdna(subdna, __direct==0)
            else:
                return subdna 

        else:
            raise ValueError("Invalid index type was specified.") 
    
    def __add__(self, other):
        if (type(other) == str and set(other) <= set("ATGCRYKMSWBDHVNatgcrykmswbdhvn-")) or type(other) == Seq:
            other = QUEEN(seq=other) 

        elif type(other) == SeqRecord:
            other = QUEEN(record=other) 

        elif type(other) == QUEEN:
            pass 
        
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized QUEEN object annot be joined with others.") 
        else:
            return joindna(self, other, __direct=0)

    def __radd__(self, other):
        if (type(other) == str and set(other) <= set("ATGCRYKMSWBDHVNatgcrykmswbdhvn-")) or type(other) == Seq:
            other = QUEEN(seq=other) 

        elif type(other) == Seq.SeqRecord:
            other = QUEEN(record=other) 

        elif type(other) == QUEEN:
            pass 
           
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized QUEEN object cannot be joined with others.") 
        else:
            return joindna(other, self, __direct==0) 
    
    def getdnaseq(self, start=None, end=None, strand=2, hide_middle=None, linebreak=None, display=False):
        whole = False
        if linebreak is None:
            width = len(self.seq) + 1
        else:
            width = linebreak

        if hide_middle is None or hide_middle > 0.5 * len(self.seq):
            hide_middle = int(0.5 * len(self.seq)) 
            whole = True

        if start is None and end is None and strand == 2:
            rcseq  = self.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
            if len(self._left_end) > hide_middle:
                left_length = hide_middle
            else:
                left_length = len(self._left_end)

            if self._left_end_top == 1:
                left_end_top = self.seq[:hide_middle]
            else:
                left_end_top = "-" * left_length + self.seq[left_length:hide_middle] 
            
            if self._left_end_bottom == 1:
                left_end_bottom = self.seq[:hide_middle].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) 
            else:
                left_end_bottom = "-" * left_length + self.seq[left_length:hide_middle].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
            
            if len(self._right_end) > hide_middle:
                right_length = hide_middle
            else:
                right_length = len(self._right_end)

            if self._right_end_top == 1:
                right_end_top = self.seq[len(self.seq)-hide_middle:]
            else:
                right_end_top = self.seq[len(self.seq)-hide_middle:len(self.seq)-right_length] + "-" * right_length 
            
            if self._right_end_bottom == 1:
                right_end_bottom = self.seq[len(self.seq)-hide_middle:].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
            else:
                right_end_bottom = self.seq[len(self.seq)-hide_middle:len(self.seq)-right_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + "-" * right_length     
            
            top = left_end_top + self.seq[hide_middle:len(self.seq)-hide_middle] + right_end_top
            bottom = left_end_bottom + self.seq[hide_middle:len(self.seq)-hide_middle].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + right_end_bottom 

        else: 
            if type(start) is DNAfeature:
                feature = start
                start   = feature.start
                end     = feature.end
                strand  = feature.strand
            if start is None:
                start = 0
            
            if end is None:
                end = len(self.seq)

            if strand is None:
                strand = 2
            
            if start >= end:
                top = self.seq[start:] + self.seq[:end]
            else:    
                top = self.seq[start:end]
            
            bottom = top.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))

        out = ""
        if display == True:
            if whole == False:
                if display == True:
                    print("5' {}...{} 3'".format(left_end_top, right_end_top))
                    print("3' {}...{} 5'".format(left_end_bottom, right_end_bottom))
                    print() 
                out += "5' {}...{} 3'\n".format(left_end_top, right_end_top)
                out += "3' {}...{} 5'".format(left_end_bottom, right_end_bottom)
            else:    
                if len(top) < width:
                    if display == True:
                        print("5' {} 3'".format(top))
                        print("3' {} 5'".format(bottom))
                        print() 
                    out += "5' {} 3'\n".format(top)
                    out += "3' {} 5'\n".format(bottom)
                else:
                    for i in range(0, len(top), width):
                        if display == True:
                            print("5' {} 3'".format(top[i:i+width]))
                            print("3' {} 5'".format(bottom[i:i+width]))
                            print()
                        out += "5' {} 3'\n".format(top[i:i+width]) 
                        out += "3' {} 5'\n".format(bottom[i:i+width])
                        out += "\n"
                    out = out.rstrip()
        
        if strand == 1:
            return_seq = top 
        elif strand == -1:
            return_seq = bottom[::-1] 
        elif strand == 2:
            return_seq = top + "/" + bottom
        
        return_seq = Qseq(return_seq) 
        return_seq.__dict__ = self.seq.__dict__ 
        return_seq.name     = "getdnaseq_{}".format(strand) 
        return return_seq

    def _setfeatureid(self):
        for i in range(0, len(self.dnafeatures)):
            if i == 0:
                self.dnafeatures[i]._id = str(1)
            else:
                self.dnafeatures[i]._id = str(i*100)
    
    def getdnafeatures(self,feature_id):
        return self._features_dict[str(feature_id)] 

    def printfeature(self, feature_list=None, attribute=None, separation=None, detail=False, seq=False, output=None, x_based_index=0):
        #“feature ID,” “qualifier:label,” “feature type,” “start position,” “end position,” and “strand"] 
        _ids        = ["feature_id"] 
        types       = ["feature_type"] 
        labels      = ["qualifier:label"]
        starts      = ["start"] 
        ends        = ["end"] 
        strands     = ["strand"]
        sequences   = ["sequence"]
        sep         = separation
        seqflag     = seq
        others_dict = {}
        
        if attribute is None:
            attribute = ["feature_id", "feature_type", "qualifier:label", "start", "end", "strand"]

        new_attribute = [] 
        for att in attribute:
            if att == "$DEFAULT":
                new_attribute += ["feature_id", "feature_type", "qualifier:label", "start", "end", "strand"]
            else:
                new_attribute.append(att)
        
        attribute = new_attribute
        if feature_list is None:
            features = list(self.dnafeatures)
            features.sort(key=lambda x:x.location.parts[0].start.position)
        else:
            features = feature_list 

        for feat in features:
            flag  = 0  
            label_keys = []
            for key in feat.qualifiers:
                if key == "label":
                    label = feat.qualifiers["label"][0] 
                    flag  = 1
                elif ("qualifier:") + key not in others_dict and ((key in attribute or ("qualifier:" + key) in attribute) or detail==True):
                    others_dict["qualifier:" + key] = ["qualifier:" + key] + (len(labels)-1) * ["null"]
                                 
            if flag == 0:
                label = "null"

            strand = feat.location.strand
            start  = feat.start
            end    = feat.end
            seq    = feat.sequence 
            
            if x_based_index == 1:
                start += 1 

            if len(attribute) > 0:
                for key in others_dict:
                    if key == "label":
                        pass 
                    elif key in feat.qualifiers or key.replace("qualifier:","") in feat.qualifiers:
                        if type(feat.qualifiers[key.replace("qualifier:","")]) == list:
                            others_dict[key].append(":".join(feat.qualifiers[key.replace("qualifier:","")]))
                        else:
                            others_dict[key].append(feat.qualifiers[key.replace("qualifier:","")])
                    else:
                        others_dict[key].append("null")
            
            if "_id" not in feat.__dict__:
                feat._id = "null"
            _ids.append(str(feat._id)) 
            labels.append(str(label)) 
            types.append(str(feat.type)) 
            starts.append(str(start)) 
            ends.append(str(end))
            sequences.append(seq) 

            if strand == 1:
                strands.append("+")
            elif strand == 0:
                strands.append("+") 
            else:
                strands.append("-")
        
        _idmax       = max(list(map(len,_ids))) + 2
        labelmax     = max(list(map(len,labels))) + 2
        ftypemax     = max(list(map(len,types)))  + 2
        startmax     = max(list(map(len,starts))) + 2
        endmax       = max(list(map(len,ends)))   + 2
        strandmax    = max(list(map(len,strands))) + 2
        sequencemax  = max(list(map(len,sequences))) + 2
        other_maxes = [max(list(map(len,others_dict[key]))) + 2 for key in others_dict]
        
        dkeys   = ["feature_id", "feature_type", "qualifier:label", "start", "end", "strand"] + list(others_dict.keys())
        dvalues = [_ids, types, labels, starts, ends, strands] + list(others_dict.values())
        dmaxes  = [_idmax, ftypemax, labelmax, startmax, endmax, strandmax] + other_maxes
        hogera = list(zip(dkeys, dvalues, dmaxes))
        rows   = [] 
        maxes  = []
        
        if detail == False:
            for key, value, max_w in zip(dkeys, dvalues, dmaxes):
                if key.replace("qualifier:","") in attribute or key in attribute:
                    rows.append(value) 
                    maxes.append(max_w) 
        else:
            rows  = [_ids, labels, types, starts, ends, strands] + list(others_dict.values())
            maxes = [_idmax, labelmax, ftypemax, startmax, endmax, strandmax] + other_maxes
        
        if seqflag == True and "sequence" not in attribute:
            attribute.append("sequence")
        
        if "sequence" in attribute:
            rows.append(sequences) 
            maxes.append(sequencemax) 
        
        if type(output) is str:
            output = open(output,"w")
            if sep == ",":
                import csv
                output = csv.writer(output) 
        for n, row in enumerate(zip(*rows)): 
            if sep is None:
                text = ""
                for m, x in enumerate(row):
                    text += x + " " * (maxes[m]-len(x)) 
                print(text, file=output) 
            else:
                if sep == ",":
                    output.writerow(row) 
                else:
                    print(*row, sep=sep, file=output)
                    
        if output is None:
            print() 

    def writedna(self, handle=None, format="genbank", record_id=None, separate_history=False):
        stdIOflag = 0 
        if handle is None:
            stdIOflag = 1
            handle    = io.StringIO()
        
        features = copy.deepcopy(self.dnafeatures)
        features.append(self._history_feature) 
        for feat in features:
            if "broken_feature" in feat.qualifiers:
                note   = feat.qualifiers["broken_feature"][0]
                label  = ":".join(note.split(":")[:-1])
                length = int(note.split(":")[-4]) 
                pos_s  = int(note.split(":")[-1].split("..")[0].replace(" ",""))
                pos_e  = int(note.split(":")[-1].split("..")[1].replace(" ",""))
                if pos_e > length:
                    note = label + ":" + str(pos_s) + ".." + str(pos_e-length)
                    feat.qualifiers["broken_feature"] = [note]

                elif (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["broken_feature"]
                
            if separate_history is not False and type(separate_history) is str:
                for feat in self.dnafeatures:
                    if feat.type == "source":
                        for key in feat.qualifiers:
                            if "building_history" in key[0:18] and import_history == True:
                                history = feat.qualifiers[key][0]
                                results = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]", history) 
                                history_num = int(key.split("_")[-1]) 
                                pairs.append((feat, history_num, history))     
                
                with open(separate_history, "w") as o:
                    print("This file describes the buiding history of '{}/{}'.".format(os.getcwd(), handle), file=o)
                    for pair in pairs:
                        del feat.qualifires["building_history" + "_" + str(pair[1])]
                        print(pair[1], pair[2], sep=",", file=o)
                self.record.annotaions["source"] = os.getcwd() + "/" + separate_history

        if type(handle) is str:
            handle = open(handle, "w") 

        self.record.features = features 
        if record_id is None:
            self.record.id = self.project
        else:
            self.record.id = record_id
        if Alphabet:
            self.record.seq = Seq(str(self.seq),Alphabet.DNAAlphabet()) 
        else:
            self.record.seq = Seq(str(self.seq))
            self.record.annotations["molecule_type"] = "DNA"
    
        SeqIO.write(self.record, handle, format)
        self.record.features = self.dnafeatures
        if stdIOflag == 1:
            print(handle.getvalue(), end="") 
