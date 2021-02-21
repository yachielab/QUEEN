import os 
import regex as re
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
from Bio import SeqIO

try:
    from Bio import Alphabet 
except ImportError:
    Alphabet = False

from Bio import Restriction
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

import Bio
import warnings
warnings.simplefilter('ignore', Bio.BiopythonParserWarning)


def _slide(feats, slide, ref_position=0, limit=None):
    new_feats = []
    for feat in feats:
        feat      = copy.deepcopy(feat)
        strand    = feat.location.strand
        new_parts = [] 
        for p in range(len(feat.location.parts)):
            flag = 0
            if feat.location.parts[p].start.position > ref_position:
                feat.location.parts[p]._start = ExactPosition(feat.location.parts[p].start.position + slide)
                if limit is not None and feat.location.parts[p].start.position >= limit:
                    feat.location.parts[p]._start = ExactPosition(int(feat.location.parts[p].start.position) - limit)
                    flag = 1
                else:
                    pass
            else:
                pass 
            
            if feat.location.parts[p].end.position > ref_position:
                feat.location.parts[p]._end = ExactPosition(feat.location.parts[p].end.position + slide)
                if limit is not None and feat.location.parts[p].end.position >= limit:
                    feat.location.parts[p]._end = ExactPosition(int(feat.location.parts[p].end.position) - limit)
                    if flag == 1:
                        flag = 0 
                    else:
                        flag = 1
                else:
                    pass
            else:
                pass 
            
            if flag == 1:
                new_parts.append(FeatureLocation(feat.location.parts[p].start.position, limit-1, feat.location.strand))
                new_parts.append(FeatureLocation(0, feat.location.parts[p].end.position, feat.location.strand))
            else:
                new_parts.append(FeatureLocation(feat.location.parts[p].start.position, feat.location.parts[p].end.position, feat.location.strand))

        if strand == -1:
            new_parts.reverse()

        if len(new_parts) > 1:
            feat.location = CompoundLocation(new_parts)
        else:
            feat.location = new_parts[0]
        
        feat.location.strand = strand
        new_feats.append(feat)
    return new_feats 

def _add_history(dna, history=""): 
    flag = 0 
    for feat in dna.dnafeature:
        if feat.type == "source":
            if "label" in feat.qualifiers and feat.qualifiers["label"][0] in dna.project and "description" in feat.qualifiers and feat.qualifiers["description"][0] == "Record of editing history" and feat.location.start == 0 and feat.location.end == len(dna.seq):
                DNA._num_history += 1 
                feat.qualifiers["construction_history_{}".format(dna._num_history)] = [history] 
                flag = 1
                break

    if flag == 0:
        _id_all = [int(feat._id) for feat in dna.dnafeature if str(feat._id).isdecimal()]
        _id_all.append(0) 
        max_id  = str(100 + math.floor(max(_id_all) / 100) * 100)
        dna._num_history = 0
        feat = SeqFeature(FeatureLocation(0, len(dna.seq), strand=1), type="source") 
        feat._id = max_id
        feat.qualifiers["label"]       = [dna.project] 
        feat.qualifiers["description"] = ["Record of editing history"]
        feat.qualifiers["construction_history_{}".format(dna._num_history)] = [history] 
        DNA._num_history += 1 
        dna.dnafeature.append(feat)
        dna._features_dict[feat._id] = feat

def exporthistory(dna):
    histories = [] 
    for feat in dna.dnafeature:
        if feat.type == "source":
            old_keys = [] 
            for key in feat.qualifiers:
                if "construction_history" in key:
                    history = feat.qualifiers[key][0] 
                    histories.append((int(key.split("_")[-1]), history)) 
    histories.sort() 
    for history in histories:
        print(history[1]) 

def archivehistory(dna):
    for feat in dna.dnafeature:
        if feat.type == "source":
            old_keys = [] 
            for key in feat.qualifiers:
                if "construction_history" in key:
                    old_keys.append(key) 
                    history = feat.qualifiers[key]
                    new_key = "archived_" + key
                    feat.qualifiers[new_key] = history
            for key in old_keys:
                del feat.qualifiers[key] 

def deletehistory(dna):
    for feat in dna.dnafeature:
        if feat.type == "source":
            old_keys = [] 
            for key in feat.qualifiers:
                if "construction_history" in key:
                    old_keys.append(key)   
            for key in old_keys:
                del feat.qualifiers[key] 

def cropdna(dna, start=0, end=0, project=None, __direct=1):
    """
    Extract subsequence from the specified region in the genome.
    All features within the subsequence including truncated features are carried over to extracrted 'dna' object. 
    """   
    if type(start) is tuple:
        start_top = start[0] 
        start_bottom = start[1] 
        start = min(start)
    else:
        start_top = start 
        start_bottom = start
    
    if type(end) is tuple: 
        end_top = end[0]
        end_bottom = end[1] 
        end = max(end)
    else:
        end_top = end
        end_bottom = end

    if start == 0 and end == len(dna.seq):
        new_dna = copy.deepcopy(dna)
        if project is None:
            project = dna.project
        new_dna.project = project
        return new_dna

    if dna.topology == "circular":
        start = len(dna.seq) + start if start < 0 else start
        start = start - len(dna.seq) if start > len(dna.seq) else start
        end   = end - len(dna.seq) if end > len(dna.seq) else end
    
    if start > end:
        if dna.topology == "circular":
            source       = copy.deepcopy(dna) 
            new_features = _slide(source.dnafeature, len(source.seq) - start, 0, limit=len(source.seq))
            new_seq      = source.seq[start:] + source.seq[:start] 
            source.seq   = new_seq
            source.dnafeature      = new_features
            source.record.features = _assigndnafeature(source.dnafeature)
            subdna = cropdna(source, 0, end+(len(new_seq)-start),__direct=0) 
            subdna.subject = SourceDNA() 
            subdna.subject.start   = start 
            subdna.subject.end     = end
            subdna.subject.strand  = 1
            subdna.subject.project = dna.project 
            subdna.subject.dna     = dna  
        else:
            raise ValueError("Start value should be larger than or equal to end value.")
    else:
        feats        = []
        for feat in dna.dnafeature:
            strand = feat.location.strand
            if strand == -1:
                s = feat.location.parts[-1].start.position
                e = feat.location.parts[0].end.position
            else:
                s = feat.location.parts[0].start.position
                e = feat.location.parts[-1].end.position 
            
            if "original" not in feat.__dict__:
                if s > e:
                    feat.original = str(dna.seq)[s:len(dna.seq)] + str(dna.seq)[:e]
                else:
                    feat.original = str(dna.seq)[s:e].upper()    
            
            feat = copy.deepcopy(feat)
            feat.full = 1
            if len(feat.location.parts) == 1 and s <= e:
                if e > start and s < end:
                    if s - start < 0:
                        feat.full = 0 #The feature is not completely included in the region. 
                        feat.location.parts[0]._start = ExactPosition(0)
                        if "cropped_region" not in feat.qualifiers:
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            if strand >= 0:
                                feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, len(feat.original) - abs(s-start), 1,  len(feat.original)) 
                        else:
                            note   = feat.qualifiers["cropped_region"]
                            if strand >= 0:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0]) + abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            else:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0]) - abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            feat.qualifiers["cropped_region"] = note
                    else:
                        feat.location.parts[0]._start = ExactPosition(s - start) 
                
                    feat.location.parts[-1]._end = ExactPosition(e - start)  
                    if feat.location.parts[-1]._end > end-start:
                        feat.full = 0 #The feature is not completely included in the region. 
                        feat.location.parts[-1]._end = ExactPosition(end - start)
                        if "cropped_region" not in feat.qualifiers: 
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            if strand >= 0: 
                                feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, 1, end-s, len(feat.original)) 
                            else:
                                feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, len(feat.original), len(feat.original)-(end-s)+1,  len(feat.original)) 
                        else:
                            s = int(feat.location.parts[0].start.position)
                            note = feat.qualifiers["cropped_region"]
                            if strand >= 0:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s + (end-start-s)-1, length)
                            else:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s - (end-start-s)+1, length)
                            feat.qualifiers["cropped_region"] = note
                    
                    feat.location.strand = strand
                    feats.append(feat)
            
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
                        feat.full = 0 
                        locations[0][0] = ExactPosition(0)
                        if "cropped_region" not in feat.qualifiers:
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            if strand >= 0:
                                feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, e-s,  abs(s-start)+1, len(feat.original)) 
                        else:
                            note   = feat.qualifiers["cropped_region"]
                            if strand >= 0:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0]) + abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            else:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0]) - abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            feat.qualifiers["cropped_region"] = note
                    else:
                        locations[0][0] = ExactPosition(s - start)
                    
                    if e > end-start and eflag == 1:
                        feat.full = 0 
                        locations[-1][1] = ExactPosition(end-start)
                        if "cropped_region" not in feat.qualifiers:
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            feat.qualifiers["cropped_region"] = "{}:{}..{}:{}".format(label, 1, end-s, len(feat.original)) 
                        else:
                            s      = int(locations[0][0])
                            note   = feat.qualifiers["cropped_region"]
                            if strand  >= 0:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s + (end-start-s)-1, length)
                            else:
                                label  = note.split("]")[0] + "]"
                                note   = note.split("]")[1] 
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s - (end-start-s)+1, length)
                            feat.qualifiers["cropped_region"] = note
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
                    feats.append(feat)
    
        feats.sort(key=lambda x:(x.location.parts[0].start.position, x.location.parts[-1].end.position))
        subdna     = copy.deepcopy(dna)
        subdna.seq = dna.seq[start:end]
        subdna.dnafeature = feats
        
        subdna._features_dict = dict(list(map(lambda x:(x._id, x), subdna.dnafeature)))
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

        new_rec          = copy.deepcopy(dna.record)
        if Alphabet:
            new_rec.seq      = Seq(dna.seq[start:end],Alphabet.DNAAlphabet()) 
        else:
            new_rec.seq      = Seq(dna.seq[start:end]) 
        new_rec.annotations["topology"] = subdna.topology
        subdna.record = new_rec
        subdna.subject = SourceDNA()
        subdna.subject.strand  = 1
        subdna.subject.start   = start 
        subdna.subject.end     = end
        subdna.subject.project = dna.project 
        subdna.subject.dna = dna 
        new_rec.features = _assigndnafeature(subdna.dnafeature)
        
    if project is None:
        project = dna.project
    subdna.project = project

    if __direct == 1:
        if project in DNA.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(DNA.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in DNA.dna_dict:
                unique += 1
            subdna._unique_id = project + "_" + str(unique)
        else:         
            subdna._unique_id = project
        
        args = [start, end, project]
        for i in range(len(args)):
            if type(args[i]) is str:
                args[i] = "'" + args[i] + "'" 
        DNA.dna_dict[subdna._unique_id] = subdna
        edit_history = "DNA.dna_dict['{}'] = cropdna('DNA.dna_dict['{}'], start={}, end={}, project={})".format(subdna._unique_id, dna._unique_id, args[0], args[1], args[2]) 
        _add_history(subdna, edit_history)
    
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
    return subdna 

def joindna(*dnas, topology="linear", project=None, __direct=1):
    """
    Join dna objects. 
    When ovhg_check argument is True, adjacent DNA objects should have common overhang sequences.  
    """
    def slide(feats,slide):
        new_feats = []
        for feat in feats:
            feat = copy.deepcopy(feat)
            strand = feat.location.strand
            for p in range(len(feat.location.parts)):
                feat.location.parts[p]._start = ExactPosition(feat.location.parts[p].start.position + slide)
                feat.location.parts[p]._end   = ExactPosition(feat.location.parts[p].end.position + slide)
            feat.location.strand = strand
            new_feats.append(feat)
        return new_feats 
    
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
            raise ValueError("The {} DNA object topology is 'circular'. Cicularized DNA objects cannot be joined with others.".format(order)) 

    construct = copy.deepcopy(dnas[0])
    if len(dnas) > 1:
        for dna in dnas[1:]:
            feats   = dna.dnafeature 
            if (dna._left_end_top * construct._right_end_bottom == 1 or dna._left_end_bottom * construct._right_end_top == 1) and (dna._left_end_top == -1 or dna._left_end_bottom == -1):
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
                    print("These DNA objects were not able to be joined. Please check sticky end seqeunces of the DNA objects")
                    return False

            else:
                new_dna = dna
                ovhg = ""
            
            feats  = slide(feats, len(construct.seq) - len(ovhg))
            feats1 = [feat for feat in construct.dnafeature if "cropped_region" in feat.qualifiers]
            feats2 = [feat for feat in feats if "cropped_region" in feat.qualifiers]
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
                                if key == "cropped_region":
                                    pass 
                                elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                                    flag = 1
                                else:
                                    flag = 0
                                    break    
                            if flag == 1:
                                note1   = feat1.qualifiers["cropped_region"]
                                label   = note1.split("]")[0] + "]"
                                note1   = note1.split("]")[1] 
                                pos_s1  = int(note1.split(":")[1].split("..")[0])
                                pos_e1  = int(note1.split(":")[1].split("..")[1]) 
                                length1 = int(note1.split(":")[2])
                                
                                note2   = feat2.qualifiers["cropped_region"] 
                                label   = note2.split("]")[0] + "]"
                                note2   = note2.split("]")[1] 
                                pos_s2  = int(note2.split(":")[1].split("..")[0])
                                pos_e2  = int(note2.split(":")[1].split("..")[1]) 
                                length2 = int(note2.split(":")[2])
                                if length1 == length2 and "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original and s2 == 0 and feat1.location.strand == feat2.location.strand:
                                    note     = "{}:{}..{}:{}".format(label, pos_s1, pos_e2, length1)
                                    new_seq  = construct.seq[s1:e1] + dna.seq[s2:e2] 
                                    new_feat = copy.deepcopy(construct.dnafeature[construct.dnafeature.index(feat1)]) 
                                    strand = new_feat.location.strand
                                    if len(feat1.location.parts) == 1 and len(feat2.location.parts) == 1:
                                        new_feat.location = FeatureLocation(feat1.location.parts[0].start.position, feat2.location.parts[-1].end.position, feat1.strand)
                                        new_feat.location.strand = strand
                                    else:
                                        locations = feat1.location.parts[0:-1] + [FeatureLocation(feat1.location.parts[-1].start.position, feat2.location.parts[0].end.position, feat1.strand)] + feat2.location.parts[0:-1]
                                        if strand == -1:
                                            locations.reverse() 
                                        new_feat.location = CompoundLocation(locations) 
                                        new_feat.location.strand = strand 
                                    
                                    if strand == -1:
                                        s = new_feat.location.parts[-1].start.position
                                        e = new_feat.location.parts[0].end.position
                                    else:
                                        s = new_feat.location.parts[0].start.position
                                        e = new_feat.location.parts[-1].end.position
                                    
                                    if len(new_seq) - len(ovhg) == e - s and len(new_seq) - len(ovhg) <= len(feat1.original):
                                        construct.dnafeature[construct.dnafeature.index(feat1)].qualifiers["cropped_region"] = note
                                        construct.dnafeature[construct.dnafeature.index(feat1)].location = new_feat.location
                                        del feats[feats.index(feat2)] 
            
            construct.dnafeature = construct.dnafeature + feats
            construct.seq               = construct.seq + new_dna.seq 
            construct._right_end         = dna._right_end
            construct._right_end_top     = dna._right_end_top
            construct._right_end_bottom  = dna._right_end_bottom
            construct.topology      = "linear"
            ovhg = dna._right_end

        construct.dnafeature.sort(key=lambda x:x.location.parts[0].start.position)
        for feat in construct.dnafeature:
            if "cropped_region" in feat.qualifiers:
                 note   = feat.qualifiers["cropped_region"]
                 label  = note.split("]")[0] + "]"
                 note   = note.split("]")[1] 
                 pos_s  = int(note.split(":")[1].split("..")[0])
                 pos_e  = int(note.split(":")[1].split("..")[1])
                 length = int(note.split(":")[2])
                 if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["cropped_region"]
        
        if Alphabet:
            new_record = SeqRecord(Seq(str(construct.seq), Alphabet.DNAAlphabet()))
        else:
            new_record = SeqRecord(Seq(str(construct.seq)))

        new_record.features = _assigndnafeature(construct.dnafeature) 
        new_record.annotations["topology"] = topology
        construct.record = new_record     
        
        if project is None:
            project = dnas[0].project 

        new_record.id     = project
        construct.project = project
        
        if topology == "circular":
            construct = _circularizedna(construct)
           
        construct.start = 0 
        construct.end = len(construct.seq) 
        
        feat = SeqFeature(FeatureLocation(0, len(construct.seq), strand=1), type="source") 
        feat.qualifiers["label"] = [construct.project] 
        construct.dnafeature.append(feat) 
        construct.record.feartures = construct.dnafeature
        construct._setfeatureid() #Update feature ID
        construct._features_dict = dict(list(map(lambda x:(x._id, x), construct.dnafeature)))
        construct.subject = SourceDNA()
        construct.query   = SourceDNA() 
    else:
        construct = _circularizedna(dnas[0])
    
    if __direct == 1:
        if project in DNA.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(DNA.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in DNA.dna_dict:
                unique += 1
            construct._unique_id = project + "_" + str(unique)
        else:         
            construct._unique_id = project
        DNA.dna_dict[construct._unique_id] = construct
        dna_elements = "[" + ",".join(["DNA.dna_dict['{}']".format(dna._unique_id) for dna in dnas]) + "]"
        edit_history = "DNA.dna_dict['{}'] = joindna(*{}, topology='{}', project='{}')".format(construct._unique_id, dna_elements, topology, project) 
        print(edit_history)  
        _add_history(construct, edit_history) 
    
    return construct

def modifyends(dna, left="", right="", add=0, add_right=0, add_left=0, project=None, __direct=1):
    """
    Set end sequence structures. 
    """
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
        left_end = left_end_top
        left_end_top    = 1
        left_end_bottom = 1

    if "/" in right_end:
        right_end_top, right_end_bottom = right_end.split("/")
    else:
        add = 1
        add_right = 1 
        right_end_top, right_end_bottom = right_end, right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))

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
        right_end = right_end_top
        right_end_top    = 1
        right_end_bottom = 1
    
    if add == 1 or (left_end != dna.seq[left_end_length-left_length:left_end_length-left_length+len(left_end)] or right_end != str(dna[len(dna.seq)-right_end_length + right_length - len(right_end):len(dna.seq)-right_end_length + right_length].seq)):
        left_end  = DNA(seq=left_end, __direct==0) 
        right_end = DNA(seq=right_end, __direct==0) 
        if add_left == 1 and add_right == 1:
            new_dna = joindna(left_end, dna, right_end, __direct=0)
            new_dna._left_end  = left_end._left_end
            new_dna._right_end = right_end._right_end
        
        elif add_right == 1:
            new_dna = cropdna(dna, start=left_end_length-left_length, end=len(dna.seq), __direct=0) + right_end
            new_dna._right_end = right_end._right_end
            new_dna._left_end  = left_end.seq
            new_dna._left_end_top  = left_end_top 
            new_dna._left_end_bottom = left_end_bottom

        else:
            new_dna = left_end + cropdna(dna, start=0, end=len(dna.seq)-right_end_length+right_length, __direct=0)
            new_dna._left_end  = left_end._left_end
            new_dna._right_end = right_end.seq
            new_dna._right_end_top = right_end_top 
            new_dna._right_end_bottom = right_end_bottom

    else:
        left_end  = DNA(seq=left_end, __direct==0) 
        right_end = DNA(seq=right_end, __direct==0)
        new_dna = cropdna(dna, start=left_end_length-left_length, end=len(dna.seq)-right_end_length+right_length) 
        new_dna._left_end  = left_end.seq
        new_dna._right_end = right_end.seq
        new_dna._left_end_top     = left_end_top 
        new_dna._left_end_bottom  = left_end_bottom
        new_dna._right_end_top    = right_end_top
        new_dna._right_end_bottom = right_end_bottom
    
    if project is None:
        project = dna.project
    new_dna.project = project 

    if __direct == 1:
        if project in DNA.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(DNA.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in DNA.dna_dict:
                unique += 1
            new_dna._unique_id = project + "_" + str(unique)
        else:         
            new_dna._unique_id = project

        DNA.dna_dict[new_dna._unique_id] = new_dna
        edit_history = "DNA.dna_dict['{}'] = modifyends(DNA.dna_dict['{}'], left='{}', right='{}', project='{}')".format(new_dna._unique_id, dna._unique_id, left, right, project) 
        _add_history(new_dna, edit_history) 

    return new_dna

def cutdna(dna, *positions, cut_stracture="", project=None, __direct=1):
    dnas = [] 
    new_positions = [] 
    for pos in positions:
        if type(pos) is int:
            pos = (pos, pos) 
        spos, epos = pos
        if spos > len(dna.seq):
            spos = spos - len(dna.seq)
        elif epos < 0:
            epos =  epos + len(dna.seq) 
        new_positions.append((spos,epos))  

    if 0 not in positions:
        positions = [(0,0)]  + new_positions 
    else: 
        positions = new_positions 

    if len(dna.seq) not in positions:
        positions = positions + [(len(dna.seq), len(dna.seq))] 
    positions = list(positions)
    positions.sort() 

    for i, pos in enumerate(positions[0:-1]):
        dnas.append(cropdna(dna,pos,positions[i+1],__direct=1))
    
    if dna.topology == "circular":
        new_dna = joindna(dnas[-1],dnas[0],__direct=0) 
        dnas    = [new_dna] + dnas[1:-1] 
    
    if project is None:
        project = dna.project
    
    for new_dna in dnas:
        new_dna.project = project

    return dnas

def flipdna(dna, project=None, __direct=1):
    """
    Return reverse complement sequence. 
    All feature infomation is also reversed with the sequence.
    """
    if type(dna) is DNA:
        dna = copy.deepcopy(dna)
        seq  = dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 
        for feat in dna.record.features:
            strand = feat.location.strand
            for p in range(len(feat.location.parts)):
                s, e = feat.location.parts[p].start.position, feat.location.parts[p].end.position
                feat.location.parts[p]._start = ExactPosition(len(dna.seq) - e) 
                feat.location.parts[p]._end   = ExactPosition(len(dna.seq) - s) 
                if "original" in feat.__dict__:
                    feat.original = feat.original.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]

            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
                if "cropped_region" in feat.qualifiers:
                    note   = feat.qualifiers["cropped_region"]
                    label  = note.split("]")[0] + "]"
                    note   = note.split("]")[1] 
                    pos_s  = int(note.split(":")[1].split("..")[0])
                    pos_e  = int(note.split(":")[1].split("..")[1])
                    length = int(note.split(":")[2])
                    note = "{}:{}..{}:{}".format(label, pos_e, pos_s, length)
                    feat.qualifiers["cropped_region"] = note
            else:
                feat.location.strand = strand
            feats.append(feat)
        seq = dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
    
    elif type(dna) is str or type(dna) is Seq:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 

    comp = DNA(seq=seq) 
    feats.sort(key=lambda x:x.location.parts[0].start.position) 
    comp.dnafeature = feats
    comp._setfeatureid()
    comp._features_dict = dict(list(map(lambda x:(x._id, x), comp.dnafeature)))
    comp.record.features = _assigndnafeature(comp.dnafeature) 
    comp._right_end, comp._left_end = dna._left_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], dna._right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
    comp._right_end_top, comp._left_end_bottom = dna._left_end_bottom, dna._right_end_top
    comp._right_end_bottom, comp._left_end_top = dna._left_end_top, dna._right_end_bottom
    if dna.subject.dna is None:
        pass 
    else:
        comp.subject = dna.subject
        comp.subject.strand = comp.subject.strand * -1 
    if project is None:
        project = dna.project
    comp.project = project

    if __direct == 1:
        if project in DNA.dna_dict:
            if project.split("_")[-1].isdecimal() == True:
                project = "_".join(project.split("_")[:-1])
            keys   = list(DNA.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in DNA.dna_dict:
                unique += 1
            comp._unique_id = project + "_" + str(unique)
        else:         
            comp._unique_id = project

        DNA.dna_dict[comp._unique_id] = comp
        edit_history = "DNA.dna_dict['{}'] = flipdna(DNA.dna_dict['{}'], project='{}')".format(comp._unique_id, dna._unique_id, project) 
        _add_history(comp, edit_history) 
    return comp

def _replacedna(dna=None, feat_list=None, target_attribute=None, query_re=None, value=None):    
    _exec =  0
    if target_attribute[0:8] == "sequence":
        attribute_regex = re.compile("sequence:\![0-9]+\.\.[0-9]+\!")
        if target_attribute == "sequence":
            ss, se = 0, len(target.seq)    
        elif attribute_regex.fullmatch(target_attribute) != None:
            ss, se = tuple(map(int,target_attribute[10:-1].split("..")))

        for feat in feat_list:
            feat = dna._features_dict[feat._id]
            if feat.location.strand == -1:
                fs, fe = feat.location.parts[-1].start.position, feat.location.parts[0].end.position
            else:
                fs, fe = feat.location.parts[0].start.position, feat.location.parts[-1].end.position
        
            s,e = fs + ss, fs + se
            target_seq = dna.seq[s:e]
            if set(value) < set("ATGCRYKMSWBDHVNatgcrykmsbdhv"):
                pass 
            else: 
                raise ValueError("Invalid sequence pattern was detected") 
            
            if query_re == "" or query_re is None:
                if value is None:
                    new_seq = re.sub(target_seq, "", target_seq)
                else:
                    new_seq = re.sub(target_seq, value, target_seq)
            else:
                new_seq = re.sub(query_re, value, target_seq)

            dna.seq = dna.seq[0:s] + new_seq + dna.seq[e:len(dna.seq)] 
            if len(new_seq) == len(target_seq): 
                pass 
            else:
                new_dnafeature = _slide(dna.dnafeature, len(new_seq) - len(target_seq), s) 
                dna.dnafeature = new_dnafeature
            _exec += 1
            
            if "label" in feat.qualifiers:
                label = feat.qualifiers["label"][0]
            else:
                label = "N.A"
            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, fs, fe))
            if "editing_history" not in feat.qualifiers:
                feat.qualifiers["editing_history"] = "{}:{}..{}:{}2{}:{}".format(label, ss, se, target_seq, new_seq, len(feat.original)) 
            else:
                feat.qualifiers["editing_history"].append("{}:{}..{}:{}2{}:{}".format(label, ss, se, target_seq, new_seq, len(feat.original)))

    else:
        new_dnafeature = []
        _id_all  = [feat._id for feat in dna.dnafeature] 
        _id_list = [feat._id for feat in feat_list] 
        for feat in dna.dnafeature:
            _del = 0
            if feat._id in _id_list:
                if target_attribute == "feature_ID" or target_attribute == "feature ID":
                    if value is None or value == "":
                        _exec += 1
                        _del = 1 
                        pass #It means remove the feature from 'dna.dnafeature'.
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
                        
                #elif target_attribute == "position":
                #   new_dnafeature.append(feat) 
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
                new_dnafeature.append(feat) 
        dna.dnafeature = new_dnafeature
        dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeature)))
    
    if _exec == 0:
        warnings.warn("Warning : No target was detected") 

def replacedna(tquery, tvalue=None):
    if tvalue is None:
        tvalue = tquery 
        tquery = None
    return functools.partial(_replacedna, query_re=tquery, value=tvalue)

def _removedna(dna=None, feat_list=None, target_attribute=None):
    dna = _replacedna(dna=dna, feat_list=feat_list, target_attribute=target_attribute, query_re=None, value="")
    
def removedna():
    return functools.partial(_removedna)

def _createdna(dna=None, feat_list=None, target_attribute=None, value=None): 
    new_dnafeature = copy.copy(dna.dnafeature)
    if target_attribute[0:8] == "sequence" or target_attribute == "strand" or target_attribute == "position":
        warnings.warn("Warning : 'sequence', 'strand,' and 'position' attributes cannote be used in 'createdna' operateion.")
        return dna 

    elif target_attribute[0:len("qualifier:")] == "qualifier:":
        if value is None:
            pass  
        else:
            key = target_attribute.split(":")[-1]
            _id_list = [feat._id for feat in feat_list if "_id" in feat.__dict__]
            for feat in new_dnafeature:
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
                    for feat2 in dna.dnafeature:
                        if feat2.location.strand == -1: 
                            s2, e2 = feat.location.part[-1].start, feat.location.part[0].end
                        else:
                            s2, e2 = feat.location.part[0].start, feat.location.part[-1].end

                        if s1 >= e2 and feat2._id.isdecimal() == True:
                            unique_num = 1
                            new_id = str(int(feat._id) + unique_num) 
                            while new_id in _id_all:
                                unique_num += 1
                                new_id = str(int(feat._id) + unique_num) 
                            flag = 1
                            break
                    feat._id = new_id
                    feat.qualifiers[key] = [] 
                    feat.qualifiers[key] = [feat.qualifiers[key]]
                    feat.qualifiers[key].append(value)
                    feat.qualifiers[key] = new_elements
                    new_dnafeature.append(feat) 
    else:
        _id_all  = [feat._id for feat in dna.dnafeature]            
        if len(feat_list) > 1:
            #warnings.warn("Some target regions were detected. To ensure the uniqueness of each feature ID, a unique number will be added after each feature ID, scuh as '{}_[0-9]'".format(value)) 
            pass
        
        for feat in feat_list:
            if value is None:
                value = ""
            if feat.location.strand == -1: 
                s1, e1 = feat.location.parts[-1].start.position, feat.location.parts[0].end.position
            else:
                s1, e1 = feat.location.parts[0].start.position, feat.location.parts[-1].end.position

            if target_attribute == "feature ID" or target_attribute == "feature_ID":
                feat_type  = "misc_feature"
                unique_num = 1
                new_id  = value
                while new_id in _id_all:
                    new_id = value + "_" + str(unique_num)
                    unique_num += 1
                
                if value == "":
                    new_id = new_id[1:]
                elif unique_num == 1: 
                    new_id = value

            elif target_attribute == "feature type" or target_attribute == "feature_type":
                feat_type = value
                flag = 1
                for feat2 in dna.dnafeature:
                    if feat2.location.strand == -1: 
                        s2, e2 = feat.location.part[-1].start, feat.location.part[0].end
                    else:
                        s2, e2 = feat.location.part[0].start, feat.location.part[-1].end

                    if s1 >= e2 and feat2._id.isdecimal() == True:
                        unique_num = 1
                        new_id = str(int(feat._id) + unique_num) 
                        while new_id in _id_all:
                            unique_num += 1
                            new_id = str(int(feat._id) + unique_num) 
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
            if "_id" not in feat.__dict__:
                feat._id  = new_id
                feat.type = feat_type 
                new_dnafeature.append(feat)  
            else:
                strand = feat.location.strand
                if dna.topology == "circular" and s1 > e1:
                    locations = [FeatureLocation(s1, len(dna.seq), strand=strand), FeatureLocation(0, e1, strand=strand)]
                    if strand == -1:
                        locations.reverse() 
                    new_feat = SeqFeature(CompoundLocation(locations), type=feat_type)
                else:
                    new_feat = SeqFeature(FeatureLocation(s1, e1, strand=strand), type=feat_type)
                new_feat._id = new_id 
                new_dnafeature.append(new_feat)

    dna.dnafeature = new_dnafeature
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeature))) 
    return dna 

def createdna(tvalue=None):
    return functools.partial(_createdna, value=tvalue)

def editdna(dna, key_attribute=None, query=None, target_attribute=None, operation=None, comment="", destructive=True, project=None, __direct=1):
    if destructive == True:
        pass 
    else:
        original_id = dna._unique_id
        dna = copy.deepcopy(dna) 

    feat_list = dna.finddna(query, attribute=key_attribute)
    if operation.func.__name__ in ("_createdna", "_removedna", "_replacedna"):
        operation(dna=dna, feat_list=feat_list, target_attribute=target_attribute)
        largs = [] 
        for item in operation.keywords:
            largs.append("=".join(item)) 
        command = operation.func.__name__ + "(" + ",".join(largs) + ")"
        if __direct == 1 and destructive == False:
            if project is None:
                project = dna.project
            
            if project in DNA.dna_dict:
                if project.split("_")[-1].isdecimal() == True:
                    project = "_".join(project.split("_")[:-1])
                keys   = list(DNA.dna_dict.keys())
                unique = 0
                while project + "_" + str(unique) in DNA.dna_dict:
                    unique += 1
                dna._unique_id = project + "_" + str(unique)
            else:         
                dna._unique_id = project
            DNA.dna_dict[dna._unique_id] = dna
            
            args = [key_attribute, query, target_attribute, command, comment]
            for i in range(len(args)):
                if type(args[i]) is str:
                    args[i] = "'" + args[i] + "'" 
            edit_history = "DNA.dna_dict['{}'] = editdna(DNA.dna_dict['{}'], key_attribute={}, query={}, target_attribute={}, operation={}, comment={})".format(dna._unique_id, original_id, args[0], args[1], args[2], args[3], args[4]) 
            _add_history(dna, history=edit_history) 

        elif destructive == True:
            edit_history = "editdna(DNA.dna_dict['{}'], key_attribute={}, query={}, target_attribute={}, operation={}, comment={})".format(args[0], args[1], args[2], args[3], args[4])
            _add_history(dna, history=edit_history) 

    else:
        raise ValueError("The operational function can be selected from only 'createdna', 'removedna', 'replacedna'.")

    if destructive == False:
        return dna 

def _circularizedna(dna):
    dna = copy.deepcopy(dna)
    seq_origin = dna.seq
    feats_origin = dna.dnafeature
    if dna.topology == "circular" and dna.record.annotations["topology"] == "circular":
        print("The DNA object is already circularized")

    if (dna._right_end_top * dna._left_end_bottom == 1 and dna._right_end_bottom * dna._left_end_top == 1) and len(dna._right_end) > 0 and (dna._left_end_top == -1 or dna._left_end_bottom == -1):
        if str(dna._right_end) == str(dna._left_end): 
            #print("The DNA object was circularized based on complementary sticky end between 3' end and 5' end. The sticky end is '{}'".format(dna._left_end)) 
            ovhg = dna._right_end
            subdna = cropdna(dna,0,len(dna.seq)-len(dna._right_end),__direct=0)
            dna.seq    = subdna.seq
            dna.record = subdna.record
        else:
            return False
    else:
        ovhg = ""

    remove_list = [] 
    feats1      = [feat for feat in dna.dnafeature if "cropped_region" in feat.qualifiers]
    for feat1 in feats1:
        if feat1.location.strand == -1:
            s1, e1 = feat1.location.parts[-1].start.position, feat1.location.parts[0].end.position
        else:
            s1, e1 = feat1.location.parts[0].start.position, feat1.location.parts[-1].end.position
        
        for feat2 in feats1:
            if feat2.location.strand == -1:
                s2, e2 = feat2.location.parts[-1].start.position, feat2.location.parts[0].end.position
            else:
                s2, e2 = feat2.location.parts[0].start.position, feat2.location.parts[-1].end.position
            
            if feat1 == feat2 or feat1 in remove_list:
                pass 
            
            elif feat1.type == feat2.type:
                flag = 0
                for key in feat1.qualifiers:
                    if key == "cropped_region":
                        pass 
                    elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                        flag = 1
                    else:
                        flag = 0
                        break    
                if flag == 1:
                    note1   = feat1.qualifiers["cropped_region"]
                    label   = note1.split("]")[0] + "]"
                    note1   = note1.split("]")[1] 
                    pos_s1  = int(note1.split(":")[1].split("..")[0])
                    pos_e1  = int(note1.split(":")[1].split("..")[1]) 
                    length1 = int(note1.split(":")[2])
                    
                    note2   = feat2.qualifiers["cropped_region"] 
                    label   = note2.split("]")[0] + "]"
                    note2   = note2.split("]")[1] 
                    pos_s2  = int(note2.split(":")[1].split("..")[0])
                    pos_e2  = int(note2.split(":")[1].split("..")[1]) 
                    length2 = int(note2.split(":")[2])
                    
                    if length1 == length2 and ((s1 == 0 and e2 == len(seq_origin)) or (s2 == 0 and e1 == len(seq_origin))):
                        if "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original:
                            note     = "{}:{}..{}:{}".format(label, pos_s1, pos_e2, length1)
                            new_seq  = seq_origin[s1:e1] + seq_origin[s2:e2] 
                            new_feat = copy.deepcopy(dna.dnafeature[dna.dnafeature.index(feat1)]) 
                            if s1 == 0:
                                new_feat.location = FeatureLocation(s2, e2 + (e1-len(ovhg)), feat1.strand)
                            else:
                                new_feat.location = FeatureLocation(s1, e1 + (e2-len(ovhg)), feat1.strand)
                            
                            if len(new_seq) - len(ovhg) <= len(feat1.original):
                                new_feat.qualifiers["cropped_region"] = note
                                if len(new_seq) - len(ovhg) == length1:
                                    del dna.dnafeature[dna.dnafeature.index(feat1)].qualifiers["cropped_region"]
                                dna.dnafeature[dna.dnafeature.index(feat1)].location = new_feat.location
                                dna.dnafeature.remove(feat2)  
                                remove_list.append(feat2) 
        
    for i in range(len(dna.dnafeature)):    
        if dna.dnafeature[i].location.parts[-1].end.position > len(dna.seq):
            if dna.dnafeature[i].location.parts[0].start.position >= len(dna.seq):
                strand                    = dna.dnafeature[i].location.strand
                dna.dnafeature[i].location = FeatureLocation(dna.dnafeature[i].location.parts[0].start.position-len(dna.seq),dna.dnafeature[i].location.parts[-1].end.position-len(dna.seq))
                dna.dnafeature[i].location.strand = strand
            else:
                strand    = dna.dnafeature[i].location.strand
                locations = [FeatureLocation(dna.dnafeature[i].location.parts[0].start.position,len(dna.seq)), FeatureLocation(0,dna.dnafeature[i].location.parts[-1].end.position-len(dna.seq))]
                if strand == -1:
                    locations.reverse() 
                dna.dnafeature[i].location = CompoundLocation(locations)
                dna.dnafeature[i].location.strand = strand
    
    dna.record.features = _assigndnafeature(dna.dnafeature)
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
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeature)))
    return dna

def _assigndnafeature(dnafeature):
    features = [] 
    for feat in dnafeature:
        if feat.location.start.position == -1:
            pass 
        else:
            features.append(feat) 
    return features 

class MatchDNA():
    def __getitem__(self, key):
        return self.__dict__[key]
    
    def __init__(self):
        self.start  = None
        self.end    = None
        self.qseq   = None
        self.sseq   = None
        self.strand = None
        self.match_count  = None

class SourceDNA():
    def __init__(self):
        self.strand  = None
        self.start   = None
        self.end     = None
        self.project = None
        self.dna     = None
        
class DNA():
    dna_dict = {}
    _num_history = 1  
    def __repr__(self):
        if len(self.seq) > 50:
            out = "<dna.DNA object; project='{}', length='{} bp', topology='{}'>".format(self.project, len(self.seq), self.topology)
        else:
            out = "<dna.DNA object; project='{}', length='{} bp', sequence='{}', topology='{}'>".format(self.project, len(self.seq), self.seq, self.topology)
        #out += self.printdnaseq(whole=False, end_length=max([10, len(str(self._left_end)), len(str(self._right_end))]), linebreak=None, display=False)
        return out 

    def __init__(self, seq=None, record=None, project=None, topology="linear", format=None, __direct=1): 
        recname = record
        if seq is None and record is None:
            if project is None: 
                project = "dna"
            self.seq                = None
            self.record             = None
            self.dnafeature         = None
            self._right_end         = None 
            self._left_end          = None
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
            else:
                fmt = None
            self.seq    = str(record.seq).upper()
            self.record = record
            if len(record.features) > 0:
                self.dnafeature = record.features
                history_nums = [1] 
                for feat in self.dnafeature:
                    if feat.type == "source":
                        pairs = []
                        for key in feat.qualifiers:
                            if "construction_history" in key:
                                history = feat.qualifiers[key][0]
                                results = re.findall("DNA.dna_dict\['.+'\]") 
                                for result in results:
                                    _unique_id = result.split("['")[1][:-2] 
                                    DNA.dna_dict[_unique_id] = None
                                history_num = int(key.split("_")[-1]) 
                                pairs.append((history_num, history)) 
                        
                        for pair in pairs:
                            new_history_num = pairs[0] + DNA._num_history
                            feat.qualifiers["construction_history_{}".format(new_history_num)] = [pairs[1]] 
                            del feat.qualifiers["construction_history_{}".format(pairs[0])]
                            history_nums.append(new_history_num) 
                DNA._num_history = max(history_nums)  
            else:
                self.dnafeature = []
            
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
            
        elif record is None:
            if project is None:
                project = "dna"
            def check_seq(seq):
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
            
            sticky   = False
            pattern1 = "[ATGCRYKMSWBDHVNatgcrykmsbdhvn*-]+/?[ATGCRYKMSWBDHVNatgcrykmsbdhvn*-]*"
            pattern2 = "[ATGCRYKMSWBDHVNatgcrykmsbdhvn]+-+[ATGCRYKMSWBDHVNatgcrykmbdhvn]+"
            pattern1 = re.compile(pattern1)   
            pattern2 = re.compile(pattern2)   
            if pattern1.fullmatch(seq) != None or seq == "": 
                if "/" in seq:
                    top, bottom, seq = check_seq(seq)
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
                

                self.dnafeature = self.record.features
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
       

        self.subject = SourceDNA()
        self.query   = SourceDNA()         
        if project in DNA.dna_dict:
            keys   = list(DNA.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in DNA.dna_dict:
                unique += 1
            self._unique_id = project + "_" + str(unique)
        else: 
            self._unique_id = project 
        
        DNA.dna_dict[self._unique_id] = self
        args = [seq, recname, project, topology, format]
        for i in range(len(args)):
            if type(args[i]) is str:
                args[i] = "'" + args[i] + "'" 
        self._setfeatureid()
        self._features_dict = dict(list(map(lambda x:(x._id, x), self.dnafeature)))
        if __direct == 1:
            edit_history = "DNA.dna_dict['{}'] = DNA(seq={}, record={}, project={}, topology={}, format={})".format(self._unique_id, args[0], args[1], args[2], args[3], args[4]) 
            _add_history(self, edit_history)
        self.record.feartures = _assigndnafeature(self.dnafeature)

    def _get_matchlist_normal(query, subject, strand, topology, min_match, error, seq_mismatch, constant):
        match_list      = []
        subject_origin  = subject 
        if error < 2:
            sub_regex     = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, error, 0, 0, error)
            ins_regex     = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, 0, error, 0, error)
            del_regex     = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, 0, 0, error, error)
        else:
            sub_regex     = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, error, 1, 1, error)
            ins_regex     = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, 1, error, 1, error)
            del_regex     = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, 1, 1, error, error)
        sub_regex     = re.compile(sub_regex)
        ins_regex     = re.compile(ins_regex)
        del_regex     = re.compile(del_regex)
         
        i = 0 
        if strand == -1:
            constant = constant.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
            if topology == "circular":
                subject = subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
                subject = subject + subject[0:len(query)-1].upper()
            else: 
                subject = subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
        else:
            if topology == "circular": 
                subject = subject.upper() + subject_origin[0:len(query)-1].upper()
            else:
                subject = subject.upper()
        
        reg = re.compile("[ATGCRYKMSWBDHV]+([ATGCRYKMSWBDHV]*-{{0,{}}})[ATGCRYKMSWBDHV]+".format(seq_mismatch))
        while 1:
            results  = [re.search(sub_regex, subject[i:]),re.search(ins_regex, subject[i:]),re.search(del_regex, subject[i:])]
            starts   = [] 
            for result in results:
                if result is not None:
                    match = MatchDNA()
                    match.sstart = result.start() + i
                    match.send   = result.end() + i 
                                    
                    if strand == -1:
                        match.sseq_na = subject[match.sstart:match.send].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper() 
                        match.sstart, match.send = len(subject_origin)-match.send, len(subject_origin)-match.sstart
                        match.qseq_na = query.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
                        match.strand = -1   
                        if match.sstart < 0:
                            match.sstart = match.sstart + len(subject_origin)    

                    else:
                        match.sseq_na   = subject[match.sstart:match.send] 
                        match.qseq_na   = query
                        match.strand    = 1
                        if match.send >= len(subject_origin):
                            match.send = match.send - len(subject_origin)    
                    
                    starts.append(result.start() + 1) 
                    if constant in match.sseq_na:
                        match_count = max([len(match.sseq_na), len(match.qseq_na)]) - Levenshtein.distance(match.sseq_na, match.qseq_na)
                        if match_count >= min_match:
                            gla        = pairwise2.align.globalms(match.sseq_na, match.qseq_na, 2, -3, -5, -2, one_alignment_only=True)
                            match.score = gla[0][2]
                            match.sseq  = gla[0][0]  
                            match.qseq  = gla[0][1]
                            new_seq    = "".join([q if s == q else "-" for s,q in zip(gla[0][0], gla[0][1])])

                            max_length = 0 
                            qresult = None 
                            for result in re.finditer(reg, new_seq):
                                s, e = result.span()
                                qseq = new_seq[s:e] 
                                if len(qseq) > max_length:
                                    max_length = len(qseq)
                                    qresult    = result
                            
                            if qresult is not None:
                                match.sseq = match.sseq[qresult.start():qresult.end()]     
                                match.qseq = match.qseq[qresult.start():qresult.end()]    
                                if strand == -1:
                                    match.qstart, match.qend = len(match.qseq_na)-(qresult.start()-gla[0][1][:qresult.start()].count("-")), len(match.qseq_na)-(qresult.end()-gla[0][1][:qresult.end()].count("-"))
                                else:
                                    match.qstart, match.qend = qresult.start()-gla[0][1][:qresult.start()].count("-"), qresult.end()-gla[0][1][:qresult.end()].count("-") 
                                
                                match.sstart, match.send = match.sstart + qresult.start()-gla[0][0][:qresult.start()].count("-"), match.sstart + qresult.end()-gla[0][0][:qresult.end()].count("-")
                                match.match_count = sum([1 if s == q else 0 for s,q in zip(match.sseq, match.qseq)]) 
                                match.sspan = (match.sstart, match.send)
                                match.qspan = (match.qstart, match.qend)
                                if len(match_list) > 0:
                                    check = sum(list(map(lambda x : (x.sspan[0] <= match.sspan[0] and match.sspan[1] <= x.sspan[1]), match_list)))
                                else:
                                    check = 0 
                                if match.match_count >= min_match and check == 0:
                                    if match.send > len(subject_origin): 
                                        match.send = match.send - len(subject_origin) 
                                    match.sspan = (match.sstart, match.send)
                                    match.qspan = (match.qstart, match.qend)
                                    match_list.append(match)
            
            if i >= len(subject) or len(starts) == 0:
                break
            i = i + min(starts)

        if len(match_list) > 0:
            match_list.sort(key=lambda x: x.match_count*-1.0)
            top_match      = match_list[0] 
            new_match_list = [] 
            new_match_list.append(top_match)

        return match_list

    def _get_matchlist_regex(query, subject, strand, topology, constant):
        match_list = []
        subject_origin = subject
        regex = re.compile(query)
        
        i = 0 
        if strand == -1:
            constant = constant.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
            if topology == "circular":
                subject = subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
                subject = subject + subject
            else: 
                subject = subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
        else:
            if topology == "circular": 
                subject = subject.upper() + subject.upper()
            else:
                subject = subject.upper()

        while 1:
            result  = re.search(regex, subject[i:]) 
            if result is not None:
                match = MatchDNA()
                match.sstart = result.start() + i
                match.send   = result.end() + i 
                                
                if strand == -1:
                    match.sseq_na = subject[match.sstart:match.send].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper() 
                    match.sstart, match.send = len(subject)-match.send, len(subject)-match.sstart
                    match.qseq_na = query.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1].upper()
                    match.strand = -1   
                    if match.sstart < 0:
                        match.sstart = match.sstart + len(subject_origin)    

                else:
                    match.sseq_na   = subject[match.sstart:match.send] 
                    match.qseq_na   = query
                    match.strand    = 1
                    if match.send >= len(subject_origin):
                        match.send = match.send - len(subject_origin)    

                           
                i = i + result.end() + 1
                if constant in match.sseq_na:
                    match.sseq  = match.sseq_na 
                    match.qseq  = query 
                    if strand == -1:
                        match.qend, match.qstart = 0, len(match.sseq)
                    else:
                        match.qstart, match.qend = 0, len(match.sseq)
                    
                    if match.sstart >= len(subject_origin):
                        match.sstart = match.sstart - len(subject_origin)
                    else:
                        pass
                    
                    if match.send >= len(subject_origin):
                        match.send = match.send - len(subject_origin)
                    else:
                        pass
                    
                    match.sspan = (match.sstart, match.send)
                    match.qspan = (match.qstart, match.qend)
                    match.match_count = None
                    if len(match_list) > 0:
                        check = sum(list(map(lambda x : (x.sspan[0] == match.sspan[0] and match.sspan[1] == x.sspan[1]), match_list)))
                    else:
                        check = 0
                    
                    if check == 0:
                        match.regex_info = result
                        match_list.append(match)
            else:
                break
        return match_list

    def finddna(self, query, attribute=None, min_match=None, max_mismatch=0, strand=None):
        #attributes = "feature ID", "feature type", "start", "end", "sequence", "query:" 
        if attribute == "start" and attribute == "end":
            raise TypeError("Numerical value cannot be specified in attrribute on 'finddna.'") 
        
        feat_list  = [] 
        if attribute is not None and attribute[0:8] == "sequence":
            return_feature_list = False
            match_list = [] 
            attribute_regex = re.compile("sequence:\|[0-9]+\.\.[0-9]+\|[+-]{0,1}")
            if attribute == "sequence":
                s = 0
                if strand is None:
                    strand  = 0
                subject = self.seq

            elif attribute_regex.fullmatch(attribute) != None:
                posinfo = attribute[10:].split("|")
                if len(posinfo[1]) == 0 and strand is None: 
                    strand = 0
                elif posinfo[1] == "+" and strand is None:
                    strand = 1
                elif strand is None:
                    strand = -1 
                s,e = tuple(map(int,posinfo[0].split("..")))
                subject = cropdna(self,s,e,__direct=0).seq
            else: 
                raise ValueError("Invalid attribute was detected.")

            if query is None:
                query = subject 
               
            if type(query) == DNA:
                query = query.seq 
            
            if set(str(query)) <= set("ATGCRYKMSWBDHVNatgcnrykmswbdhv"):
                mode = "normal"
            else:
                mode = "regex" 

            if mode == "normal": 
                error = max_mismatch
                seq_mismatch = error 
                if min_match is None:
                    min_match = len(query) - error
                else:
                    pass 

                constant = ""                
                if strand == 0:
                    match_list = DNA._get_matchlist_normal(query, subject, 1, self.topology, min_match, error, seq_mismatch, constant) 
                    match_list.extend(DNA._get_matchlist_normal(query, subject, -1, self.topology, min_match, error, seq_mismatch, constant)) 
                elif strand == 1:
                    match_list = DNA._get_matchlist_normal(query, subject, 1, self.topology, min_match, error, seq_mismatch, constant)
                elif strand == -1:
                    match_list = DNA._get_matchlist_normal(query, subject, -1, self.topology, min_match, error, seq_mismatch, constant)
                
                match_list.sort(key=lambda x: x.score*-1.0)
                querydna = DNA(seq=query) 
            
            elif mode == "regex":
                constant = "" 
                match_list = DNA._get_matchlist_regex(query, subject, 1, self.topology, constant) 
                match_list.extend(DNA._get_matchlist_regex(query, subject, -1, self.topology, constant))
            
            for match in match_list:
                flag = 0
                for feat in self.dnafeature:
                    strand = feat.location.strand
                    if strand == -1:
                        start = feat.location.parts[-1].start.position
                        end   = feat.location.parts[0].end.position 
                    else:
                        start = feat.location.parts[0].start.position
                        end   = feat.location.parts[-1].end.position 
                    
                    if start == match.sstart + s and end == match.send + s:
                        flag = 1  
                        feat_list.append(feat) 
                    else:
                        pass 
                
                if flag == 0:
                    if match.sstart + s > match.send + s and self.topology == "circular":
                        locations = [[match.sstart+s, len(dna.seq), match.strand], [0, match.send+s, match.strand]]
                        if match.strand == -1:
                            new_feat = SeqFeature(CompoundLocation(locations), type="misc_feature") 
                        else:
                            locations.reverse()
                            new_feat = SeqFeature(CompoundLocation(locations), type="misc_feature")
                    else:
                        new_feat = SeqFeature(FeatureLocation(match.sstart+s, match.send+s, strand=match.strand), type="misc_feature")
                    feat_list.append(new_feat) 

            del match_list

        if attribute is None or attribute == "feature type" or attribute == "feature_type":
            if query is None:
                query = ".+"
            cquery = re.compile(query)
            for feat in self.dnafeature:
                if cquery.fullmatch(feat.type) != None:
                    feat_list.append(copy.deepcopy(feat))  
                    
        if attribute is None or attribute == "feature ID" or attribute == "feature_ID":
            if query is None:
                query = ".+"
            cquery = re.compile(str(query))
            for key, feat in self._features_dict.items():
                if cquery.fullmatch(key) != None:
                    feat_list.append(copy.deepcopy(feat)) 

        if attribute is None or attribute == "strand":
            for feat in self.dnafeature:
                strand = feat.location.strand
                if query is None or strand == query or (query.isdecimal() and strand == query):
                    feat_list.append(copy.deepcopy(feat)) 

        if attribute is None or attribute[0:len("qualifier:")] == "qualifier:":
            if query is None:
                query = ".+"
            
            cquery = re.compile(query)
            if attribute is None or attribute == "qualifier:*":
                for feat in self.dnafeature:
                    flag = 0 
                    for key in feat.qualifiers:
                        if type(feat.qualifiers[key]) is list:
                            pass 
                        else:
                            feat.qualifiers[key] = [feat.qualifiers[key]]
                        
                        for element in feat.qualifiers[key]:
                            if cquery.fullmatch(element) != None:
                                feat_list.append(copy.deepcopy(feat)) 
                                flag = 1
                                break
                    if flag == 1:
                        break
            else:
                key = attribute.split(":")[-1] 
                for feat in self.dnafeature:
                    if key in feat.qualifiers:
                        if type(feat.qualifiers[key]) is list:
                            pass 
                        else:
                            feat.qualifiers[key] = [feat.qualifiers[key]]
                        
                        for element in feat.qualifiers[key]:
                            if cquery.fullmatch(element) != None:
                                feat_list.append(copy.deepcopy(feat)) 
                                break
        #else:
        #    raise ValueError("Invalid attribute was detected.") 
        
        return feat_list

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
                return reverse_complement(subdna)
            else:
                return subdna 

        else:
            raise ValueError("Invalid index type was specified.") 
    
    def __add__(self, other):
        if (type(other) == str and set(other) <= set("ATGCRYKMSWBDHVNatgcrykmswbdhvn")) or type(other) == Seq:
            other = DNA(seq=other) 

        elif type(other) == SeqRecord:
            other = DNA(record=other) 

        elif type(other) == DNA:
            pass 
        
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized DNA object cannot be joined with others.") 
        else:
            return joindna(self, other)

    def __radd__(self, other):
        if (type(other) == str and set(other) <= set("ATGCRYKMSWBDHVNatgcrykmswbdhvn")) or type(other) == Seq:
            other = DNA(seq=other) 

        elif type(other) == Seq.SeqRecord:
            other = DNA(record=other) 

        elif type(other) == DNA:
            pass 
           
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized DNA object cannot be joined with others.") 
        else:
            return joindna(other, self) 
    
    def getdnaseq(self, region=None, whole=True, end_length=10, linebreak=None, display=False):
        if region == None or ((type(region) is list or type(region) is tuple) and region[0] == 0 and region[1] == 1 and len(region) == 2):
            strand = 0 
            if linebreak is None:
                width = len(self.seq) + 1
            else:
                width = linebreak
            if end_length is None:
                if len(self.seq) - len(self._right_end) - len(self._left_end) > 100:
                    end_length = 5 
                elif len(self.seq) - len(self._right_end) - len(self._left_end) > 10:
                    end_length = 2 
                else:
                    end_length = 0 
            
            if end_length > 0.5 * len(self.seq):
                end_length = int(0.5 * len(self.seq)) 
                whole = True

            seq_rc = self.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
            if len(self._left_end) > end_length:
                left_length = end_length
            else:
                left_length = len(self._left_end)

            if self._left_end_top == 1:
                left_end_top = self.seq[:end_length]
            else:
                left_end_top = "-" * left_length + self.seq[left_length:end_length] 
            
            if self._left_end_bottom == 1:
                left_end_bottom = self.seq[:end_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) 
            else:
                left_end_bottom = "-" * left_length + self.seq[left_length:end_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
            
            if len(self._right_end) > end_length:
                right_length = end_length
            else:
                right_length = len(self._right_end)

            if self._right_end_top == 1:
                right_end_top = self.seq[len(self.seq)-end_length:]
            else:
                right_end_top = self.seq[len(self.seq)-end_length:len(self.seq)-right_length] + "-" * right_length 
            
            if self._right_end_bottom == 1:
                right_end_bottom = self.seq[len(self.seq)-end_length:].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
            else:
                right_end_bottom = self.seq[len(self.seq)-end_length:len(self.seq)-right_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + "-" * right_length     
            
            top = left_end_top + self.seq[end_length:len(self.seq)-end_length] + right_end_top
            bottom = left_end_bottom + self.seq[end_length:len(self.seq)-end_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + right_end_bottom 

             
        else:
            if type(region) is list or type(region) is tuple:
                if len(region) == 2: 
                    strand = 0 
                else:
                    strand = region[2]
                start = region[0]
                end   = region[1]

            elif type(region) == SeqFeature:
                strand = region.location.strand 
                if strand == -1:
                    start = region.location.parts[-1].start.position
                    end   = region.location.parts[0].end.position
                else: 
                    start = region.location.parts[0].start.position
                    end   = region.location.parts[-1].end.position
            else:
                start = 0
                end   = len(self.seq)
            top, bottom = cropdna(self, start, end).getdnaseq()
            bottom = bottom[::-1]
        
        out = ""
        if display == True:
            if whole == False:
                if display == True:
                    print("5' {}...{} 3'".format(left_end_top, right_end_top))
                    print("3' {}...{} 5'".format(left_end_bottom, right_end_bottom))
                out += "5' {}...{} 3'\n".format(left_end_top, right_end_top)
                out += "3' {}...{} 5'".format(left_end_bottom, right_end_bottom)
            else:    
                if len(top) < width:
                    if display == True:
                        print("5' {} 3'".format(top))
                        print("3' {} 5'".format(bottom))
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

        elif display == False:
            if strand == 1:
                return top
            elif strand == -1:
                return bottom[::-1]
            else:
                return top, bottom[::-1]
    
    def _setfeatureid(self):
        for i in range(0, len(self.dnafeature)):
            self.dnafeature[i]._id = str(i*100)
    
    def getdnafeature(self,feature_id):
        return self._features_dict[str(feature_id)] 

    def printfeature(self, separation=None, output=None, feature_type=None, feature_id=None, attribute=["feature ID", "feature type", "qualifier:label", "start position", "end position", "strand"], detail=False, seq=False, zero_based_index=True):
        #“feature ID,” “qualifier:label,” “feature type,” “start position,” “end position,” and “strand"] 
        _ids        = ["feature ID"] 
        types       = ["feature type"] 
        labels      = ["qualifier:label"]
        starts      = ["start position"] 
        ends        = ["end position"] 
        strands     = ["strand"]
        sequences   = ["sequence"]
        sep         = separation
        seqflag     = seq
        others_dict = {}
        
        new_attribute = [] 
        for att in attribute:
            if att == "$DEFAULT":
                new_attribute += ["feature ID", "feature type", "qualifier:label", "start position", "end position", "strand"]
            else:
                new_attribute.append(att)
        
        attribute = new_attribute
        feature_id_list = feature_id

        if feature_id_list is None:
            feature_id_list = list(self._features_dict.keys())
        else:
            feature_id_list = list(map(int, feature_id_list)) 
        
        feature_id_list.sort()  
        features = [self._features_dict[str(_id)] for _id in feature_id_list]
        features.sort(key=lambda x:x.location.parts[0].start.position)
        for feat in features:
            if feature_type is None or feat.type in feature_type:
                flag  = 0  
                label_keys = []
                for key in feat.qualifiers:
                    if key == "label":
                        label = feat.qualifiers["label"][0] 
                        flag  = 1
                    elif ("qualifier:") + key not in others_dict and ((key in attribute or ("qualifier:" + key) in attribute) or detail==True):
                        others_dict["qualifier:" + key] = ["qualifier:" + key] + (len(labels)-1) * ["N.A."]
                                     
                if flag == 0:
                    label = "N.A."
                
                strand = feat.location.strand
                if strand == -1:
                    start = feat.location.parts[-1].start.position
                    end   = feat.location.parts[0].end.position 
                    seq   = cropdna(self,start,end,__direct=0).seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                else:
                    start = feat.location.parts[0].start.position
                    end   = feat.location.parts[-1].end.position 
                    seq   = cropdna(self,start,end,__direct=0).seq
                
                if zero_based_index == False:
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
                            others_dict[key].append("N.A.")

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
        
        dkeys   = ["feature ID", "feature type", "qualifier:label", "start position", "end position", "strand"] + list(others_dict.keys())
        dvalues = [_ids, labels, types, starts, ends, strands] + list(others_dict.values())
        dmaxes  = [_idmax, labelmax, ftypemax, startmax, endmax, strandmax] + other_maxes
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

    def writedna(self, handle, format="genbank", record_id=None):
        features = copy.deepcopy(self.dnafeature)
        for feat in features:
            if "cropped_region" in feat.qualifiers:
                note   = feat.qualifiers["cropped_region"]
                label  = note.split("]")[0] + "]"
                note   = note.split("]")[1] 
                pos_s  = int(note.split(":")[1].split("..")[0])
                pos_e  = int(note.split(":")[1].split("..")[1])
                length = int(note.split(":")[2])
                
                if pos_e > length:
                    note = label + ":" + str(pos_s) + ".." + str(pos_e-length) + ":" + str(length)
                    feat.qualifiers["cropped_region"] = note

                elif (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["cropped_region"]
        
        self.record.features = features 
        if record_id is None:
            self.record.id = self.project
        else:
            self.record.id = record_id
        if Alphabet:
            self.record.seq = Seq(str(self.seq),Alphabet.DNAAlphabet()) 
        else:
            self.record.seq = Seq(str(self.seq)) 
        SeqIO.write(self.record, handle, format)
        self.record.features = self.dnafeature
