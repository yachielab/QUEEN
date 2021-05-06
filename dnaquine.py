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
import regex as re
from Bio import SeqIO

try:
    from Bio import Alphabet 
except ImportError:
    Alphabet = False

from Bio import pairwise2
from Bio import BiopythonParserWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, FeatureLocation, ExactPosition
import visualize_circular_dna as vc
import visualize_linear_dna as vl 

warnings.simplefilter('ignore', BiopythonParserWarning)

def _add_history(dna, history="", combine=0): 
    flag = 0 
    topfeat = None
    for feat in dna.dnafeatures:
        if feat.type == "source":
            if "label" in feat.qualifiers and feat.qualifiers["label"][0] in dna.project and "description" in feat.qualifiers and feat.qualifiers["description"][0] == "Record of building history" and feat.location.start == 0 and feat.location.end == len(dna.seq) and len(feat.location.parts) == 1:
                DNA._num_history += 1 
                feat.qualifiers["building_history_{}".format(DNA._num_history)] = [history.replace(" ","–")] 
                topfeat = feat
                flag = 1
                break

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
        DNA._num_history += 1 
        feat.qualifiers["building_history_{}".format(DNA._num_history)] = [history.replace(" ","–")] 
        feat = DNAfeature(feature=feat, subject=dna)         
        dna.dnafeatures.append(feat)
        dna._features_dict[feat._id] = feat
        topfeat = feat 

    remove_list = [] 
    if combine == 1:
        for feat in dna.dnafeatures:
            if feat.type == "source" and "description" in feat.qualifiers and feat.qualifiers["description"][0] == "Record of building history" and feat is not topfeat:
                for key in feat.qualifiers:
                    if "building_history" in key:
                        topfeat.qualifiers[key] = feat.qualifiers[key] 
                    else:
                        pass 
                remove_list.append(feat)
    for feat in remove_list:
        dna.dnafeatures.remove(feat) 
    
    dna._features_dict  = dict(list(map(lambda x:(x._id, x), dna.dnafeatures)))
    dna.record.features = dna.dnafeatures

def cutdna(dna, *positions, crop=False, project=None, process_description=None, __direct=1):
    dna = copy.deepcopy(dna) 
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description 

    def extract(dna, start, end, project=None): 
        start_top    = start[0] 
        start_bottom = start[1] 
        start = min(start)
        
        end_top    = end[0]
        end_bottom = end[1] 
        end = max(end)

        if start == 0 and end == len(dna.seq):
            new_dna = copy.deepcopy(dna)
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
                source   = copy.deepcopy(dna) 
                subdna1 = extract(source, [start, start], [len(dna.seq), len(dna.seq)])
                subdna2 = extract(source, [0,0], [end,end])
                subdna  = joindna(subdna1, subdna2, __direct=0)
                subdna.subject = SourceDNA() 
                subdna.subject.start   = start 
                subdna.subject.end     = end
                subdna.subject.strand  = 1
                subdna.subject.project = dna.project 
                subdna.subject.dna     = dna  
        else:
            if start > end and dna.topoloy == "linear":
                raise ValueError("Start value should be larger than or equal to end value.")
            feats = []
            new_features = []
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
                        label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat1.original), feat1.original, s, e))
                        if strand >= 0:
                            feat1.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, len(dna.seq)-s)]
                        else:
                            feat1.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(dna.seq)-s, 1)]

                    else:
                        note = feat.qualifiers["broken_feature"]
                        note = note[0] if type(note) is list else note 
                        if strand >= 0:
                            label  = note.split("]")[0] + "]"
                            length = int(note.split("]")[0].split(":")[2])
                            note   = note.split("]")[1]
                            pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s, pos_s + len(dna.seq)-s)
                        else:
                            label  = note.split("]")[0] + "]"
                            length = int(note.split("]")[0].split(":")[2])
                            note   = note.split("]")[1]
                            pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s, pos_s - (len(dna.seq)-s))
                        feat1.qualifiers["broken_feature"] = [note]

                    if "broken_feature" not in feat2.qualifiers:
                        label = feat2._id
                        label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat2.original), feat2.original, s, e))
                        if strand >= 0:
                            feat2.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(dna.seq)-s+1, len(dna.seq)-s+e)]
                        else:
                            feat2.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(dna.seq)-s+e, len(dna.seq)-s+1)]

                    else:
                        note   = feat.qualifiers["broken_feature"][0]
                        if strand >= 0:
                            label  = note.split("]")[0] + "]"
                            length = int(note.split("]")[0].split(":")[2])
                            note   = note.split("]")[1]
                            pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s + len(dna.seq)-s, pos_e)
                        else:
                            label  = note.split("]")[0] + "]"
                            length = int(note.split("]")[0].split(":")[2])
                            note   = note.split("]")[1]
                            pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                            pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                            note   = "{}:{}..{}".format(label, pos_s - (len(dna.seq)-s), pos_e)
                        feat2.qualifiers["broken_feature"] = [note]
                    new_features.append(DNAfeature(feature=feat1))
                    new_features.append(DNAfeature(feature=feat2))
                else:
                    #print(feat, start, end) 
                    new_features.append(DNAfeature(feature=feat))

            for feat in new_features:
                strand = feat.location.strand
                s = feat.start
                e = feat.end
                if "_original" not in feat.__dict__:
                    if s > e:
                        feat._original = str(dna.seq)[s:len(dna.seq)] + str(dna.seq)[:e]
                    else:
                        feat._original = str(dna.seq)[s:e].upper()    
                    #feat.original_location = feat.location
                feat = copy.deepcopy(feat)
                feat.full = 1
                if len(feat.location.parts) == 1 and s <= e:
                    if e > start and s < end:
                        if s - start < 0:
                            feat.full = 0 #The feature is not completely included in the region. 
                            feat.location.parts[0]._start = ExactPosition(0)
                            if "broken_feature" not in feat.qualifiers:
                                label = feat._id
                                label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), feat.original, s, e))
                                if strand >= 0:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, abs(s-start)+1, e-s)] 
                                else:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(feat.original) - abs(s-start), 1)] 
                            else:
                                note = feat.qualifiers["broken_feature"][0]
                                if strand >= 0:
                                    label  = note.split("]")[0] + "]"
                                    length = int(note.split("]")[0].split(":")[2])
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ","")) + abs(s-start) 
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                else:
                                    label  = note.split("]")[0] + "]"
                                    length = int(note.split("]")[0].split(":")[2])
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ","")) - abs(s-start) 
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                feat.qualifiers["broken_feature"] = [note]
                        else:
                            feat.location.parts[0]._start = ExactPosition(s - start) 
                    
                        feat.location.parts[-1]._end = ExactPosition(e - start)  
                        if feat.location.parts[-1]._end > end-start:
                            feat.full = 0 #The feature is not completely included in the region. 
                            feat.location.parts[-1]._end = ExactPosition(end - start)
                            if "broken_feature" not in feat.qualifiers: 
                                label = feat._id
                                label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), feat.original, s, e))
                                if strand >= 0: 
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, end-s)]
                                else:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, len(feat.original), len(feat.original)-(end-s)+1)]
                            else:
                                s = int(feat.location.parts[0].start.position)
                                note = feat.qualifiers["broken_feature"][0]
                                if strand >= 0:
                                    label  = note.split("]")[0] + "]"
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s + (end-start-s)-1)
                                else:
                                    label  = note.split("]")[0] + "]"
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
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
                            feat.full = 0 
                            locations[0][0] = ExactPosition(0)
                            if "broken_feature" not in feat.qualifiers:
                                label = "feat._id"
                                label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), feat.original, s, e))
                                if strand >= 0:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, abs(s-start)+1, e-s)]
                                else:
                                    feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, e-s,  abs(s-start)+1)] 
                            else:
                                note   = feat.qualifiers["broken_feature"][0]
                                if strand >= 0:
                                    label  = note.split("]")[0] + "]"
                                    length = int(note.split("]")[0].split(":")[2])
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ","")) + abs(s-start) 
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                else:
                                    label  = note.split("]")[0] + "]"
                                    length = int(note.split("]")[0].split(":")[2])
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ","")) - abs(s-start) 
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_e)
                                feat.qualifiers["broken_feature"] = [note]
                        else:
                            locations[0][0] = ExactPosition(s - start)
                        
                        if e > end-start and eflag == 1:
                            feat.full = 0 
                            locations[-1][1] = ExactPosition(end-start)
                            if "broken_feature" not in feat.qualifiers:
                                label = feat._id
                                label = "[{}]".format("{}:{}:{}:{}:{}..{}".format(dna.project, label, len(feat.original), feat.original, s, e))
                                feat.qualifiers["broken_feature"] = ["{}:{}..{}".format(label, 1, end-s)]
                            else:
                                s      = int(locations[0][0])
                                note   = feat.qualifiers["broken_feature"][0]
                                if strand  >= 0:
                                    label  = note.split("]")[0] + "]"
                                    length = int(note.split("]")[0].split(":")[2]) 
                                    note   = note.split("]")[1] 
                                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ","")) 
                                    note   = "{}:{}..{}".format(label, pos_s, pos_s + (end-start-s)-1)
                                else:
                                    label  = note.split("]")[0] + "]"
                                    length = int(note.split("]")[0].split(":")[2])
                                    note   = note.split("]")[1] 
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
            subdna     = copy.deepcopy(dna)
            subdna.seq = dna.seq[start:end]
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
            new_rec.features = _assigndnafeatures(subdna.dnafeatures)
            
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
        return subdna 
    
    dnas = [] 
    new_positions = [] 
    for pos in positions:
        if type(pos) is str:
            pos = tuple(map(int,pos.split("/")))
            
        elif type(pos) is SeqFeature or type(pos) is DNAfeature:
            strand = pos.location.strand
            if "_digestion_top" not in pos.__dict__:
                if strand != -1:
                    start = int(pos.location.parts[0].end) 
                    end   = int(pos.location.parts[0].end)
                else:
                    start = int(pos.location.parts[0].start)
                    end   = int(pos.location.parts[0].start)
            else:
                if strand != -1:
                    start = int(pos.location.parts[0].start) + pos._digestion_top
                    end   = int(pos.location.parts[0].start) + pos._digestion_bottom
                else:
                    start = int(pos.location.parts[0].end) - pos._digestion_bottom
                    end   = int(pos.location.parts[0].end) - pos._digestion_top 
            pos = (start, end) 

        elif type(pos) is int or type(pos) is Qint:
            pos = (pos, pos) 
        
        spos, epos = pos
        if spos > len(dna.seq):
            spos = spos - len(dna.seq)
        elif epos < 0:
            epos =  epos + len(dna.seq) 
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
    else:
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
        dna_keys = list(DNA.dna_dict.keys())
        for i in range(len(dnas)):
            if dnas[i].project in DNA.dna_dict:
                if project.split("_")[-1].isdecimal() == True:
                    project = "_".join(project.split("_")[:-1])
                keys   = list(DNA.dna_dict.keys())
                unique = 0
                while project + "_" + str(unique) in dna_keys:
                    unique += 1
                dnas[i]._unique_id = project + "_" + str(unique)
            else:         
                dnas[i]._unique_id = project
            dna_keys.append(dnas[i]._unique_id) 
            products.append("DNA.dna_dict['{}']".format(dnas[i]._unique_id))

        args = [] 
        assign_histories = [] 
        for new_pos, pos in zip(new_positions, positions):
            if type(pos) == DNAfeature:
                if "_second_id" in pos.__dict__:
                    qkey = pos._qkey
                    for qindex, qfeat in enumerate(DNA.queried_features_dict[qkey]):
                        if qfeat._second_id == pos._second_id:
                            break
                    assign_histories.append("DNA.queried_feature_dict['{}'] = DNA.queried_features_dict['{}'][{}]".format(second_id, qkey, qindex))
                    args.append("DNA.queried_feature_dict['{}']".format(pos._second_id))
                else:
                    if "label" in pos.qualifiers:
                        label = "_" + pos.qualifiers["label"][0] 
                    else:
                        label = ""
                    second_id     = dna._unique_id + "_" + pos._second_id + label
                    pos._second_id = second_id
                    assign_histories.append("DNA.queried_feature_dict['{}'] = DNA.finddna(key_attribute='feature_id', query='{}')".format(second_id, pos.feature_id))
                    args.append("DNA.queried_feature_dict['{}']".format(second_id))
            
            elif type(pos) == Qint:
                qkey = pos.qkey
                for qindex, qfeat in enumerate(DNA.queried_features_dict[qkey]):
                    if qfeat._second_id == pos.parent_id:
                        break
                assign_histories.append("DNA.queried_feature_dict['{}'] = DNA.queried_features_dict['{}'][{}]".format(pos.parent_id, pos.qkey, qindex))
                args.append("DNA.queried_feature_dict['{}'].{}".format(pos.parent_id, pos.name))

            else:
                new_pos = "/".join(list(map(str,new_pos)))
                new_pos = "'" + new_pos + "'"
                args.append(new_pos)
        
        if type(project) is str:
            project = "'" + project + "'"
        if type(process_description) is str:
            process_description = "'" + process_description + "'"
        if len(products) > 1:
            edit_history = "{} = cutdna(DNA.dna_dict['{}'], {}, crop={}, project={}, process_description={})".format(", ".join(products), dna._unique_id, ", ".join(args), str(crop), project, process_description) 
        else:
            edit_history = "{}, = cutdna(DNA.dna_dict['{}'], {}, crop={}, project={}, process_description={})".format(", ".join(products), dna._unique_id, ", ".join(args), str(crop), project, process_description)
        for subdna in dnas:
            for history in assign_histories:
                _add_history(subdna, history) 
            _add_history(subdna, edit_history)
 
    if crop == True:
        return dnas[0], crop_positions 
    else:
        return dnas

def cropdna(dna, start=0, end=None, project=None, process_description=None, __direct=1):
    """
    Extract subsequence from the specified region in the genome.
    All features within the subsequence including truncated features are carried over to extracrted 'dna' object. 
    """ 
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description 

    if project is None:
        project = dna.project
    
    if end is None:
        end = len(dna.seq) 

    subdna, crop_positions = cutdna(dna, start, end, project=project, crop=True, __direct=0)  
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
        
        args = []
        assign_histories = [] 
        for new_pos, pos in zip(crop_positions, (start, end)):
            if type(pos) == DNAfeature:
                if "second_id" in pos.__dict__:
                    qkey = pos._qkey
                    for qindex, qfeat in enumerate(DNA.queried_features_dict[qkey]):
                        if qfeat._second_id == pos._second_id:
                            break
                    assign_histories.append("DNA.queried_feature_dict['{}'] = DNA.queried_features_dict['{}'][{}]".format(second_id, qkey, qindex))
                    args.append("DNA.queried_feature_dict['{}']".format(pos._second_id))
                else:
                    if "label" in pos.qualifiers:
                        label = "_" + pos.qualifiers["label"][0] 
                    else:
                        label = ""
                    second_id     = dna._unique_id + "_" + pos._second_id + label
                    pos._second_id = second_id
                    assign_histories.append("DNA.queried_feature_dict['{}'] = DNA.finddna(key_attribute='feature_id', query='{}')".format(second_id, pos.feature_id))
                    args.append("DNA.queried_feature_dict['{}']".format(second_id))
            
            elif type(pos) == Qint:
                qkey = pos.qkey
                for qindex, qfeat in enumerate(DNA.queried_features_dict[qkey]):
                    if qfeat._second_id == pos.parent_id:
                        break
                assign_histories.append("DNA.queried_feature_dict['{}'] = DNA.queried_features_dict['{}'][{}]".format(pos.parent_id, pos.qkey, qindex))
                args.append("DNA.queried_feature_dict['{}'].{}".format(pos.parent_id, pos.name))

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
        
        DNA.dna_dict[subdna._unique_id] = None
        edit_history = "DNA.dna_dict['{}'] = cropdna(DNA.dna_dict['{}'], start={}, end={}, project={}, process_description={})".format(subdna._unique_id, dna._unique_id, args[0], args[1], args[2], args[3]) 
        
        for history in assign_histories:
            _add_history(subdna, history) 
        _add_history(subdna, edit_history)
    
    return subdna

def joindna(*dnas, topology="linear", project=None, process_description=None, __direct=1):
    """
    Join dna objects. 
    When ovhg_check argument is True, adjacent DNA objects should have common overhang sequences.  
    """
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description 

    def slide(feats,slide):
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
    
    new_dnas = [] 
    for i, dna in enumerate(dnas):
        if type(dna) is not DNA:
            dna = DNA(dna)

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
        new_dnas.append(dna) 
    
    dnas = new_dnas
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
                    print("These DNA objects were not able to be joined. Please check sticky end seqeunces of the DNA objects")
                    return False

            else:
                new_dna = dna
                ovhg = ""
            
            feats  = slide(feats, len(construct.seq) - len(ovhg))
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
                                note1    = feat1.qualifiers["broken_feature"][0]
                                label    = note1.split("]")[0] + "]"
                                length1  = int(note1.split("]")[0].split(":")[2])
                                note1    = note1.split("]")[1] 
                                pos_s1   = int(note1.split(":")[1].split("..")[0].replace(" ",""))
                                pos_e1   = int(note1.split(":")[1].split("..")[1].replace(" ","")) 
                                
                                note2   = feat2.qualifiers["broken_feature"][0]
                                label   = note2.split("]")[0] + "]"
                                length2 = int(note2.split("]")[0].split(":")[2]) 
                                note2   = note2.split("]")[1] 
                                pos_s2  = int(note2.split(":")[1].split("..")[0].replace(" ",""))
                                pos_e2  = int(note2.split(":")[1].split("..")[1].replace(" ","")) 
                                
                                #Join fragmented features
                                if length1 == length2 and "_original" in feat1.__dict__ and "_original" in feat2.__dict__ and feat1.original == feat2.original and feat1.location.strand == feat2.location.strand:
                                    note    = "{}:{}..{}".format(label, pos_s1, pos_e2)
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
                                    #if (len(new_seq) - len(ovhg) == e - s and len(new_seq) - len(ovhg) <= len(feat1.original))
                                    #if new_feat.original == construct.getdnaseq(new_feat1.start, new_feat2.end, new_feat.location.strand):
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
                 label  = note.split("]")[0] + "]"
                 length = int(note.split("]")[0].split(":")[2])
                 note   = note.split("]")[1] 
                 pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                 pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                 if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["broken_feature"]
        
        if Alphabet:
            new_record = SeqRecord(Seq(str(construct.seq), Alphabet.DNAAlphabet()))
        else:
            new_record = SeqRecord(Seq(str(construct.seq)))

        new_record.features = _assigndnafeatures(construct.dnafeatures) 
        new_record.annotations["topology"] = topology
        construct.record = new_record     
        
                
        if topology == "circular":
            construct = _circularizedna(construct)
        construct.start = 0 
        construct.end = len(construct.seq) 
        construct._setfeatureid() #Update feature ID
        construct.subject = SourceDNA()
    else:
        topology  = "circular"
        construct = _circularizedna(dnas[0])
    
    if project is None:
        project = dnas[0].project

    construct.record.id = project
    construct.project   = project
    construct._features_dict = dict(list(map(lambda x:(x._id, x), construct.dnafeatures)))
    construct.record.feartures = construct.dnafeatures
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
        #DNA.dna_dict[construct._unique_id] = construct
        DNA.dna_dict[construct._unique_id] = None
        
        if type(process_description) is str:
            process_description = "'" + process_description + "'"
        if type(project) is str:
            project =  "'" + project + "'" 
        dna_elements = "[" + ", ".join(["DNA.dna_dict['{}']".format(dna._unique_id) for dna in dnas]) + "]"
        edit_history = "DNA.dna_dict['{}'] = joindna(*{}, topology='{}', project={}, process_description={})".format(construct._unique_id, dna_elements, topology, project, process_description) 
        _add_history(construct, edit_history, combine=1) 
    
    for dnafeature in construct.dnafeatures:
        dnafeature.subject = construct
    return construct

def modifyends(dna, left="", right="", add=0, add_right=0, add_left=0, project=None, process_description=None, __direct=1):
    """
    Set end sequence structures. 
    """
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description 
    
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
        left_end  = DNA(seq=left_end,  _direct=0) 
        right_end = DNA(seq=right_end, _direct=0) 
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
        left_end  = DNA(seq=left_end,  _direct=0) 
        right_end = DNA(seq=right_end, _direct=0)
        new_dna = cropdna(dna, start=left_end_length-left_length, end=len(dna.seq)-right_end_length+right_length, __direct=0) 
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

        #DNA.dna_dict[new_dna._unique_id] = new_dna
        DNA.dna_dict[new_dna._unique_id] = None
        if type(process_description) is str:
            process_description = "'" + process_description + "'"
        if type(project) is str:
            project =  "'" + project + "'" 
        edit_history = "DNA.dna_dict['{}'] = modifyends(DNA.dna_dict['{}'], left='{}', right='{}', project={}, process_description={})".format(new_dna._unique_id, dna._unique_id, left_origin, right_origin, project, process_description) 
        _add_history(new_dna, edit_history, combine=1) 
    
    for dnafeature in new_dna.dnafeatures:
        dnafeature.subject = new_dna
    return new_dna

def flipdna(dna, project=None, process_description=None, __direct=1):
    """
    Return reverse complement sequence. 
    All feature infomation is also reversed with the sequence.
    """
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description 

    if type(dna) is DNA:
        dna = copy.deepcopy(dna)
        seq  = dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 
        for feat in dna.dnafeatures:
            strand = feat.location.strand
            for p in range(len(feat.location.parts)):
                s, e = feat.location.parts[p].start.position, feat.location.parts[p].end.position
                feat.location.parts[p]._start = ExactPosition(len(dna.seq) - e) 
                feat.location.parts[p]._end   = ExactPosition(len(dna.seq) - s) 
            
            if "original" in feat.__dict__:
                feat._original = feat.original.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                #feat.original_location.parts  = list(reversed(feat.original_location.parts))
                #feat.original_location.strand = -1 * feat.original_location.strand 

            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
                if "broken_feature" in feat.qualifiers:
                    note   = feat.qualifiers["broken_feature"][0]
                    label  = note.split("]")[0] + "]"
                    length = int(note.split("]")[0].split(":")[2]) 
                    note   = note.split("]")[1] 
                    pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                    pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                    note = "{}:{}..{}".format(label, pos_e, pos_s)
                    feat.qualifiers["broken_feature"] = [note]
            else:
                feat.location.strand = strand
            feats.append(feat)
        seq = dna.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
    
    elif type(dna) is str:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        return seq

    elif type(dna) is Seq:
        seq   = str(dna) 
        seq   = seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        feats = [] 

    comp = DNA(seq=seq) 
    feats.sort(key=lambda x:x.location.parts[0].start.position) 
    comp.dnafeatures = feats
    comp._setfeatureid()
    comp._features_dict = dict(list(map(lambda x:(x._id, x), comp.dnafeatures)))
    comp.record.features = _assigndnafeatures(comp.dnafeatures) 
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

        #dna.dna_dict[comp._unique_id] = comp
        DNA.dna_dict[comp._unique_id] = None
        if type(project) is str:
            project =  "'" + project + "'" 
        if type(process_description) is str:
            process_description =  "'" + process_description + "'" 
        edit_history = "DNA.dna_dict['{}'] = flipdna(DNA.dna_dict['{}'], project={}, process_description={})".format(comp._unique_id, dna._unique_id, project, process_description) 
        _add_history(comp, edit_history) 
    
    for dnafeature in comp.dnafeatures:
        dnafeature.subject = comp
    return comp

def _replaceattribute(dna=None, feat_list=None, target_attribute=None, query_re=None, value=None):    
    _exec =  0
    if target_attribute[0:8] == "sequence":
        attribute_regex = re.compile("sequence:\![0-9]+\.\.[0-9]+\!")
        
        for feat in feat_list:
            feat   = dna._features_dict[feat._id]
            fs, fe = feat.start, feat.end 
            overlap_feats = []
            for feat2 in dna.dnafeatures:
                strand = feat2.location.strand
                s2 = feat2.start
                e2 = feat2.end 
                if e2 >= fs and s2 <= fe:
                    overlap_feats.append(feat2) 
            
            if target_attribute == "sequence":
                strand = feat.location.strand
                if feat.location.strand == -1: 
                    ss, se = 0, abs(feat.location.parts[-1].start.position - feat.location.parts[0].end.position)
                else:
                    ss, se = 0, abs(feat.location.parts[0].start.position - feat.location.parts[-1].end.position) 
            
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
           
            if strand == -1:
                if s == 0 and e == len(dna.seq):
                    new_dna = flipdna(DNA(seq=new_seq)) 
                elif s == 0:
                    new_dna = flipdna(DNA(seq=new_seq)) + dna[e:len(dna.seq)] 
                elif e == len(dna.seq):
                    new_dna = dna[0:s] + flipdna(DNA(seq=new_seq))
                else:
                    new_dna = dna[0:s] + flipdna(DNA(seq=new_seq)) + dna[e:len(dna.seq)] 
            else:
                if s == 0 and e == len(dna.seq):
                    new_dna = DNA(seq=new_seq)
                elif s == 0:
                    new_dna = DNA(seq=new_seq) + dna[e:len(dna.seq)] 
                elif e == len(dna.seq):
                    new_dna = dna[0:s] + DNA(seq=new_seq)
                else:
                    new_dna = dna[0:s] + DNA(seq=new_seq) + dna[e:len(dna.seq)] 
            
            #dna = new_dna
            dna.seq = new_dna.seq
            dna.record = new_dna.record
            dna.dnafeatures = new_dna.dnafeatures
            _exec += 1
            if "label" in feat.qualifiers:
                label = feat.qualifiers["label"][0]
            else:
                label = "N.A"
            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, fs, fe))

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
                        
                #elif target_attribute == "position":
                #   new_dnafeatures.append(feat) 
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

def replaceattribute(query_re, value=None):
    if value is None:
        value    = query_re
        query_re = None
    return functools.partial(_replaceattribute, query_re=query_re, value=value)

def _removeattribute(dna=None, feat_list=None, target_attribute=None):
    dna = _replaceattribute(dna=dna, feat_list=feat_list, target_attribute=target_attribute, query_re=None, value="")
    
def removeattribute():
    return functools.partial(_removeattribute)

def _createattribute(dna=None, feat_list=None, target_attribute=None, value=None): 
    new_dnafeatures = copy.copy(dna.dnafeatures)
    _id_all  = [feat._id for feat in dna.dnafeatures]            
    if target_attribute[0:8] == "sequence" or target_attribute == "strand" or target_attribute == "position":
        warnings.warn("Warning : 'sequence', 'strand,' and 'position' attributes cannote be used in 'createattribute' operateion.")
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
    else:
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
                flag = 1

            if target_attribute == "feature id" or target_attribute == "feature_id":
                feat_type  = "misc_feature"
                unique_num  = 1
                unique_flag = 0 
                
                if len(value) > 0 and value[-1] == "*":
                    unique_flag = 1
                    new_id = value[:-1]
                elif value == "":
                    pass
                else: 
                    new_id = value 

                while (new_id in _id_all and unique_flag == 1):
                    new_id = value[:-1] + str(unique_num)
                    unique_num += 1
                
                if value == "":
                    pass
                
                elif unique_num == 1:
                    if unique_flag == 0:
                        new_id = value
                    else:
                        new_id = value[:-1]
            elif target_attribute == "feature type" or target_attribute == "feature_type":
                feat_type = value
                            
            _id_all.append(new_id)
            if "_id" not in feat.__dict__:
                feat._id  = new_id
                feat.type = feat_type 
                new_dnafeatures.append(feat)  
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
                new_feat = DNAfeature(feature=new_feat, subject=dna) 
                new_dnafeatures.append(new_feat)
    
    dna.dnafeatures = list(dict([(feat._id, feat) for feat in new_dnafeatures]).values()) 
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeatures))) 
    return dna 

def createattribute(value=None):
    return functools.partial(_createattribute, value=value)

def editsequence(dna, query=None, key_attribute=None, min_match=None, max_mismatch=0, target_range=None, old_value=None, new_value=None, project=None, process_description=None, __direct=1):
    if target_range is None:
        target_ragne     = (0, len(dna.seq)) 
    target_attribute = "sequence:!{}..{}!".format(target_range[0], target_range[1])        
    new_dna = editfeature(dna, query=query, key_attribute=key_attribute, min_match=None, max_mismatch=max_mismatch, target_attribute=target_attribute, operation=replaceattribute(old_value, new_value), project=project, new_copy=True, process_description=None, __funcname="editsequence", __direct=0) 
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description  
    
    if __direct == 1:
        original_id = dna._unique_id
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
        DNA.dna_dict[dna._unique_id] = None
        process_description = "'" + process_description + "'" if process_description is not None else None
        old_value = "'" + old_value + "'" if old_value is not None else None
        new_value = "'" + new_value + "'" if new_value is not None else None
        edit_history = "DNA.dna_dict['{}'] = editsequence(DNA.dna_dict['{}'], query='{}', key_attribute='{}', min_match={}, max_mismatch={}, trarget_ragne=({},{}), operation=replaceattribute({}, {}), project={}, process_description={})".format(dna._unique_id, original_id, query, key_attribute, min_match , max_mismatch, target_range[0], target_range[1], old_value, new_value, project, process_description)
        _add_history(dna, history=edit_history)  
    return new_dna

def editfeature(dna, query=None, key_attribute=None, min_match=None, max_mismatch=0, target_attribute=None, operation=None, project=None, new_copy=True, process_description=None, __funcname="editfeature", __direct=1): 
    if __funcname == "editfeature" and "sequence" in str(target_attribute):
        AttributeError("'editfeature' cannot edit 'seqeunce' attribute. Plase use 'editsequence' fucntion")
    if process_description is None:
        process_description = DNA.process_description
    else:
        DNA.process_description = process_description 

    def _get_matchlist_normal(query, subject, strand, topology, min_match, error, seq_mismatch, constant):
        match_list      = []
        subject_origin  = subject
        if error > 0:
            seed_length = 7 if len(query) > 7 else len(query)
            seed_length = min_match if min_match > seed_length else seed_length 
        else:
            seed_length = len(query) 

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
                        
        seeds = [] 
        for i in range(0, len(query)-seed_length+1):
            seeds.append(query[i:i+seed_length]) 
        
        match_positions = []
        for seed in seeds:
            i = 0
            while 1:
                match = re.search(seed, subject[i:]) 
                if match is not None:
                    s = match.start() + i - (len(query)-seed_length)
                    s = 0 if s < 0 else s 
                    e = match.end() + i + len(query)-seed_length
                    match_positions.append((s,e,seed)) 
                    i = match.start() + i + 1 
                else:
                    i += 1
                if i >= len(subject):
                    break
        
        #print(match_positions) 
        sub_regex = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, error, 1, 1, error)
        ins_regex = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, 1, error, 1, error)
        del_regex = "(?e)({}){{s<={},i<={},d<={},e<={}}}".format(query, 1, 1, error, error)
        sub_regex = re.compile(sub_regex)
        ins_regex = re.compile(ins_regex)
        del_regex = re.compile(del_regex) 

        reg = re.compile("[ATGCRYKMSWBDHV]+([ATGCRYKMSWBDHV]*-{{0,{}}})[ATGCRYKMSWBDHV]+".format(seq_mismatch))
        for region in match_positions:
            i = region[0]
            target = subject[region[0]:region[1]] 
            results  = [re.search(sub_regex, target),re.search(ins_regex, target),re.search(del_regex, target)]
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
                        if match.send >= len(subject_origin) and topology=="circular":
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
                                if (match.match_count >= min_match or len(match.qseq_na)-match.match_count <= error) and check == 0:
                                    if match.send > len(subject_origin) and topology=="circular": 
                                        match.send = match.send - len(subject_origin) 
                                    match.sspan = (match.sstart, match.send)
                                    match.qspan = (match.qstart, match.qend)
                                    match_list.append(match)
            

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
                    if match.send >= len(subject_origin) and topology == "circular":
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
                    
                    if match.send >= len(subject_origin) and topology=="circular":
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

    def search(dna, query, attribute=None, min_match=None, max_mismatch=0, strand=None): 
        #attributes = "feature ID", "feature type", "start", "end", "sequence", "query:" 
        re_digestion    = re.compile("[ATGCRYKMSWBDHVN]+[\^]?\([0-9]+/[0-9]+\)") 
        attribute_regex = re.compile("sequence:\|[0-9]+\.\.[0-9]+\|[+-]{0,1}")
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
            if type(query) is str and (set(query) <= set("ATGCWSN") or re_digestion.fullmatch(query) != None):
                attribute = "sequence"
            else:
                pass 

        if attribute == "start" and attribute == "end":
            raise TypeError("Positional value can be specified by 'sequence:|int..int|'")         
        
        feat_list  = [] 
        if attribute is not None and attribute[0:8] == "sequence":
            return_feature_list = False
            match_list = [] 
            if attribute == "sequence":
                s = 0
                e = len(dna.seq) 
                if strand is None:
                    strand  = 0
                subject = dna.seq

            elif attribute_regex.fullmatch(attribute) != None:
                posinfo = attribute[10:].split("|")
                if len(posinfo[1]) == 0 and strand is None: 
                    strand = 0
                elif posinfo[1] == "+" and strand is None:
                    strand = 1
                elif strand is None:
                    strand = -1 
                s,e = tuple(map(int,posinfo[0].split("..")))
                subject = cropdna(dna,s,e,__direct=0).seq
            else: 
                raise ValueError("Invalid attribute was detected.")

            if query is None:
                query = subject 
                new_feat = SeqFeature(FeatureLocation(s, e, strand=1), type="misc_feature")
                new_feat = DNAfeature(feature=new_feat, subject=dna, query=subject)
                feat_list.append(new_feat) 

            else:
                if type(query) == DNA:
                    query = query.seq 
                
                qorigin = query
                hat = 0
                hat_pos_top    = 0
                hat_pos_bottom = 0
                query = str(query).upper() 
                if set(str(query)) <= set("ATGCRYKMSWBDHVN^_/()0123456789"):
                    if "/" in query:
                        re_digestion = re.compile("[ATGCRYKMSWBDHVN]+[\^]?\([0-9]+/[0-9]+\)") 
                        if re_digestion.fullmatch(query) != None:
                            query  = query.replace("^","") 
                            common = query.split("(")[0]
                            ds, de = tuple(map(int, query.split("(")[1][0:-1].split("/")))
                            query  = common + (ds * "N") + "^" + ((de-ds) * "N") + "/" + common.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + (de * "N") + "_"
                        query_top    = query.split("/")[0] 
                        query_bottom = query.split("/")[1]
                        if len(query_top) == len(query_bottom) and query_top.replace("^","") == query_bottom.replace("_","").translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")):
                            pass 
                        else:
                            raise ValueError("Top strand sequence should be compatible with bottom strand sequence.") 
                    else:
                        query_top    = query.replace("_","")
                        query_bottom = query.replace("^","").translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))
                    
                    mode = "normal"
                    if "^" in query_top:
                        hat = 1 
                        hat_pos_top = query_top.index("^")
                        query = query_top.replace("^","")                
                    
                    if "_" in query_bottom:
                        hat = 1 
                        hat_pos_bottom = query_bottom.index("_") 

                    elif "^" in query_top and query_top.replace("^","") == query_bottom[::-1]:
                        hat = 1
                        hat_pos_bottom = len(query) - hat_pos_top

                    if len(set(str(query)) & set("RYKMSWBDHVN()0123456789")) > 0:
                        query = query.replace("R","[GA]")
                        query = query.replace("Y","[TC]") 
                        query = query.replace("K","[GT]") 
                        query = query.replace("M","[AC]") 
                        query = query.replace("S","[GC]") 
                        query = query.replace("W","[AT]")
                        query = query.replace("B","[GTC]")
                        query = query.replace("D","[GAT]")
                        query = query.replace("H","[ACT]")
                        query = query.replace("V","[GCA]")
                        query = query.replace("N","[AGCT]")
                        mode = "regex"
                else:
                    mode = "regex" 
                
                if mode == "normal": 
                    error = max_mismatch
                    seq_mismatch = error 
                    if min_match is None:
                        min_match = len(query)
                    else:
                        #error = len(query) - min_match
                        seq_mismatch = error 

                    constant = ""                
                    if strand == 0:
                        match_list = _get_matchlist_normal(query, subject, 1, dna.topology, min_match, error, seq_mismatch, constant) 
                        match_list.extend(_get_matchlist_normal(query, subject, -1, dna.topology, min_match, error, seq_mismatch, constant)) 
                    elif strand == 1:
                        match_list = _get_matchlist_normal(query, subject, 1, dna.topology, min_match, error, seq_mismatch, constant)
                    elif strand == -1:
                        match_list = _get_matchlist_normal(query, subject, -1, dna.topology, min_match, error, seq_mismatch, constant)
                    
                    match_list.sort(key=lambda x: x.score*-1.0)
                
                elif mode == "regex":
                    constant = "" 
                    match_list = _get_matchlist_regex(query, subject, 1, dna.topology, constant) 
                    match_list.extend(_get_matchlist_regex(query, subject, -1, dna.topology, constant))
                
                for match in match_list:
                    flag = 0
                    for feat in dna.dnafeatures:
                        strand = feat.location.strand
                        start  = feat.start
                        end    = feat.end 
                        if start == match.sstart + s and end == match.send + s:
                            flag = 1  
                            if "note" not in feat.qualifiers:
                                feat.qualifiers["note"] = [] 
                            feat.query_seq = qorigin.replace("^","").replace("_","")
                            feat.qualifiers["note"].append("query_sequence:{}".format(feat.query_seq)) 
                            feat_list.append(feat) 
                        else:
                            pass 
                    
                    if flag == 0:
                        if match.sstart + s > match.send + s and dna.topology == "circular":
                            locations = [[match.sstart+s, len(dna.seq), match.strand], [0, match.send+s, match.strand]]
                            if match.strand == -1:
                                locations.reverse()
                                new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type="misc_feature") 
                            else:
                                new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type="misc_feature")
                        else:
                            new_feat = SeqFeature(FeatureLocation(match.sstart+s, match.send+s, strand=match.strand), type="misc_feature")
                        new_feat = DNAfeature(feature=new_feat, subject=dna, query=qorigin.replace("^","").replace("_","")) 
                        if hat == 1:
                            new_feat._digestion_top    = hat_pos_top 
                            new_feat._digestion_bottom = hat_pos_bottom
                        new_feat.qualifiers["note"] = ["query_sequence:{}".format(new_feat.query_seq)] 
                        feat_list.append(new_feat) 
                        
                del match_list

        if attribute is None or attribute == "feature type" or attribute == "feature_type":
            if query is None:
                query = ".+"
            cquery = re.compile(query)
            for feat in dna.dnafeatures:
                if cquery.fullmatch(feat.type) != None or query == feat.type:
                    feat_list.append(copy.deepcopy(feat))  
                    
        if attribute is None or attribute == "feature id" or attribute == "feature_id":
            if query is None:
                query = ".+"
            cquery = re.compile(str(query))
            for key, feat in dna._features_dict.items():
                if cquery.fullmatch(key) != None or query == key:
                    feat_list.append(copy.deepcopy(feat)) 

        if attribute is None or attribute == "strand":
            for feat in dna.dnafeatures:
                strand = feat.location.strand
                if query is None or strand == query or (query.isdecimal() and strand == query):
                    feat_list.append(copy.deepcopy(feat)) 

        if attribute is None or attribute[0:len("qualifier:")] == "qualifier:":
            if query is None:
                query = ".+"
            cquery = re.compile(query)
            if attribute is None or attribute == "qualifier:*":
                for feat in dna.dnafeatures:
                    flag = 0 
                    for key in feat.qualifiers:
                        if type(feat.qualifiers[key]) is list:
                            pass 
                        else:
                            feat.qualifiers[key] = [feat.qualifiers[key]]
                        
                        for element in feat.qualifiers[key]:
                            if cquery.fullmatch(element) != None or query == element:
                                feat_list.append(copy.deepcopy(feat)) 
                                flag = 1
                                break
                    if flag == 1:
                        break
            else:
                key = attribute.split(":")[-1] 
                for feat in dna.dnafeatures:
                    if key in feat.qualifiers:
                        if type(feat.qualifiers[key]) is list:
                            pass 
                        else:
                            feat.qualifiers[key] = [feat.qualifiers[key]]
                        
                        for element in feat.qualifiers[key]:
                            if cquery.fullmatch(element) != None or query == element:
                                feat_list.append(copy.deepcopy(feat)) 
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

    if new_copy == False:
        pass 
    else:
        original_id = dna._unique_id
        dna = copy.deepcopy(dna) 

    feat_list = search(dna, query, attribute=key_attribute, min_match=min_match, max_mismatch=max_mismatch) 
    if operation is None:
        return feat_list

    elif operation.func.__name__ in ("_createattribute", "_removeattribute", "_replaceattribute"):
        operation(dna=dna, feat_list=feat_list, target_attribute=target_attribute)
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
            #DNA.dna_dict[dna._unique_id] = dna 
            DNA.dna_dict[dna._unique_id] = None
            args = [query, key_attribute, min_match, max_mismatch, target_attribute, command, project, new_copy, process_description]
            for i in range(len(args)):
                if type(args[i]) is str and i != 5:
                    args[i] = "'" + args[i] + "'" 
            edit_history = "DNA.dna_dict['{}'] = editfeature(DNA.dna_dict['{}'], query={}, key_attribute={}, min_match={}, max_mismatch={}, target_attribute={}, operation={}, project={}, new_copy={}, process_description={})".format(dna._unique_id, original_id, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8]) 
            _add_history(dna, history=edit_history) 

        elif new_copy == False:
            args = [query, key_attribute, min_match, max_mismatch, target_attribute, command, project, new_copy, process_description]
            for i in range(len(args)):
                if type(args[i]) is str and i != 5:
                    args[i] = "'" + args[i] + "'" 
            edit_history = "editfeature(DNA.dna_dict['{}'], query={}, key_attribute={}, min_match={}, max_mismatch={}, target_attribute={}, operation={}, project={}, new_copy={}, process_description={})".format(dna._unique_id, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])
            _add_history(dna, history=edit_history) 

    else:
        raise ValueError("The operational function can be selected from only 'createattribute', 'removeattribute', 'replaceattribute'.")

    if new_copy == True:
        return dna 

#def leavecomment(#comment):
#    DNA.process_description = comment 

def _circularizedna(dna):
    dna = copy.deepcopy(dna)
    seq_origin = dna.seq
    feats_origin = dna.dnafeatures
    if dna.topology == "circular" and dna.record.annotations["topology"] == "circular":
        print("The DNA object is already circularized")

    if (dna._right_end_top * dna._left_end_bottom == 1 and dna._right_end_bottom * dna._left_end_top == 1) and len(dna._right_end) > 0 and (dna._left_end_top == -1 or dna._left_end_bottom == -1):
        if str(dna._right_end) == str(dna._left_end): 
            #print("The DNA object was circularized based on complementary sticky end between 3' end and 5' end. The sticky end is '{}'".format(dna._left_end)) 
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
                    label   = note1.split("]")[0] + "]"
                    length1 = int(note1.split("]")[0].split(":")[2])
                    note1   = note1.split("]")[1] 
                    pos_s1  = int(note1.split(":")[1].split("..")[0].replace(" ",""))
                    pos_e1  = int(note1.split(":")[1].split("..")[1].replace(" ","")) 
                    
                    note2   = feat2.qualifiers["broken_feature"][0] 
                    label   = note2.split("]")[0] + "]"
                    length2 = int(note2.split("]")[0].split(":")[2])
                    note2   = note2.split("]")[1] 
                    pos_s2  = int(note2.split(":")[1].split("..")[0].replace(" ",""))
                    pos_e2  = int(note2.split(":")[1].split("..")[1].replace(" ","")) 
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

def _assigndnafeatures(dnafeatures):
    features = [] 
    for feat in dnafeatures:
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

class Qint(int):
    def __init__(self, num):
        self.qkey = None
        self.parent_id = None
        self.name = None 

class DNAfeature(SeqFeature):
    def __getattr__(self, name):
        if name == "feature_id":
            return self._id
        elif name == "feature_type":
            return self.type
        elif name == "seq":
            return self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        elif name == "_original":
            return self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        elif name == "original":
            if "_original" in self.__dict__:
                return self._original 
            else:
                return self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        elif name == "sequence":
            return self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        elif name == "strand":
            return self.location.stran
        elif name == "start":
            return self._start
        elif name == "end":
            return self._end 
        elif name == "span":
            return (self.start, self.end)
        else:
            raise AttributeError

    def __setattr__(self, key, value):
        if key in ["feature_id", "feature_type", "seq", "sequenece", "original", "strand", "start", "end", "span"]:
            raise AttributeError("Cannot assign to '{}' attribute. To set or change feature attribute value, plase use 'editfeature' or 'editseqeunce' function.")
        else:
            object.__setattr__(self, key, value)

    def __init__(self, feature=None, location=None, type="misc_feature", subject=None, query=None):
        if feature is None:
            SeqFeature.__init__(self, location, type)
        else:
            self.__dict__ = feature.__dict__
        
        #Start->End direction should be 5' to 3' on the top strand.
        if self.location.strand == -1:
            self._start = self.location.parts[-1].start.position
            self._end   = self.location.parts[0].end.position
        else:
            self._start = self.location.parts[0].start.position
            self._end   = self.location.parts[-1].end.position
        
        self._start = Qint(self._start)
        self._end   = Qint(self._end) 
        self._start.name = "start"
        self._end.name   = "end"

        if subject is None:
            pass 
        else:
            self._digestion_top    = self._end - self._start if self._start < self._end else len(subject.seq) - self._end + self._start
            self._digestion_bottom = self._digestion_top
        self.subject   = subject
        self.query_seq = query

    def set_position(self, position, attribute):
        if self.subject.topology == "linear" and  position[0] > position[1] and attribute == "start":
            raise ValueError("'start' value should be smaller than 'end' value. First, please change 'end' value") 
        
        elif self.subject.topology == "linear" and  position[0] > position[1] and attribute =="end":
            raise ValueError("'end' value should be larger than 'start' value. First, please change 'start' value") 

        strand   = self.location.strand
        if type(position[0]) is int and type(position[1]) is int:
            s = position[0]
            e = position[1] 
            if s < e:
                location = FeatureLocation(s, e, strand) 
            else:
                location = CompoundLocation([FeatureLocation(s, len(self.subject.seq)), FeatureLocation(0, e, strand)]) 
            self.location = location 

        elif type(value[0]) is list and type(value[1]) is list:
            locations = [] 
            strand = self.location.strand
            for s,e in zip(value[0], value[1]):
                if s < e:
                    location = FeatureLocation(s, e, strand)
                    locations.append(location)
                else:
                    loc1 = FeatureLocation(s, len(self.subject.seq), strand)
                    loc2 = FeatureLocation(0, e, strand)
                    locations.append(loc1) 
                    locations.append(loc2) 
            if strand == -1:
                locations.reverse()
            self.location = CompoundLocation(locations, type=feat.type) 
        
        if self.location.strand == -1:
            self._start = self.location.parts[-1].start.position
            self._end   = self.location.parts[0].end.position
        else:
            self._start = self.location.parts[0].start.position
            self._end   = self.location.parts[-1].end.position
        
        self._start = Qint(self._start)
        self._end   = Qint(self._end) 
        self._start.name = "start"
        self._end.name   = "end"

def quine(dna, output=None, description_only=False):
    histories  = dna.history   
    hindex = 0
    for index, history in enumerate(histories):
        if "DNA.dna_dict" in history[1].replace(" ","").replace("–"," "):
            hindex = index
        else:
            pass 
    result     = re.findall("DNA.dna_dict\['[^\[\]]+'\]",histories[hindex][1].replace(" ","").replace("–"," "))[0] 
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
    
    if type(output) is str:
        o = open(output, "w") 
    elif output is None:
        o = None
    
    pre_process_description = "''"
    if description_only == False:
        print("import sys", file=o)  
        #print("#sys.path.append(\"{}".format(os.getcwd().replace(" ","\ ").replace("(","\(").replace(")","\)") + "/" + __file__[:-6] + "\")"), file=o)
        print("sys.path.append(\"{}".format(__file__[:-11] + "\")"), file=o)
        print("from dnaquine import *", file=o)  
    
    num = -1
    for h, (process_description, history) in enumerate(zip(descriptions, new_histories)):
        if description_only == True:
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

            if "DNA.queried_feature_dict" in history[1][0:len("DNA.queried_feature_dict")]:
                print(history[1], file=o)
            elif process_description == None or num == -1:
                print(history[1][:-1] + ", process_description={})".format("None"), file=o)
            else:
                print(history[1][:-1] + ", process_description=description{})".format(num), file=o)
        pre_process_description = process_description 

    if description_only == False:
        print(result + ".writedna('reconstructed_{}.gbk')".format(dna.project), file=o)

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

class DNA():
    dna_dict = {}
    queried_feature_dict  = {}
    queried_features_dict = {}
    process_description   = None
    _num_history = 1  
    _qnum = 0
    def __repr__(self):
        if len(self.seq) > 50:
            out = "<dna.DNA object; project='{}', length='{} bp', topology='{}'>".format(self.project, len(self.seq), self.topology)
        else:
            out = "<dna.DNA object; project='{}', length='{} bp', sequence='{}', topology='{}'>".format(self.project, len(self.seq), self.seq, self.topology)
        #out += self.printdnaseq(whole=False, end_length=max([10, len(str(self._left_end)), len(str(self._right_end))]), linebreak=None, display=False)
        return out 
    
    """
    def __setattr__(self, key, value):
        if key == "project" and "_unique_id" in self.__dict__ and self._unique_id in DNA.dna_dict:
            old_id = self._unique_id 
            if value in DNA.dna_dict:
                if value.split("_")[-1].isdecimal() == True:
                    value = "_".join(value.split("_")[:-1])
                keys   = list(DNA.dna_dict.keys())
                unique = 0
                while value + "_" + str(unique) in DNA.dna_dict:
                    unique += 1
                self._unique_id = value + "_" + str(unique)
            else:         
                self._unique_id = value
            
            history = self.history
            for i in range(len(history)):
                if "['{}']".format(old_id) in self.history:
                    self.hisotry[i] = self.history[i].replace("['{}']".format(old_id), "['{}']".format(self._unique_id))
            print("hogehogehoge", self._unique_id)
            DNA.dna_dict[self._unique_id] = None        
            #del DNA.dna_dict[old_id]
        object.__setattr__(self, key, value)
    """
    
    def __getattr__(self, name): 
        if name == "history":
            histories = [] 
            for feat in self.dnafeatures:
                if feat.type == "source":
                    old_keys = [] 
                for key in feat.qualifiers:
                    if "building_history" in key[0:18]: 
                        history = feat.qualifiers[key][0] 
                        histories.append((int(key.split("_")[-1]), history)) 
            histories.sort()
            return histories

        if name == "sequence":
            self.sequence = self.seq
        else:
            raise AttributeError

    def __init__(self, seq=None, record=None, project=None, topology="linear", format=None, process_description=None, import_history=True, _direct=1):
        if process_description is None:
            process_description = "" 
        else:
            DNA.process_description = process_description 
        recname = record
        if seq is None and record is None:
            if project is None: 
                project = "dna"
            self.seq                = None
            self.record             = None
            self.dnafeatures        = None
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
            if len(record.features) > 0:
                self.dnafeatures = [] 
                for feat in record.features:
                    self.dnafeatures.append(DNAfeature(feature=feat, subject=self))
                history_nums = [DNA._num_history] 
                
                pairs = [] 
                for feat in self.dnafeatures:
                    if feat.type == "source":
                        for key in feat.qualifiers:
                            if "building_history" in key[0:18] and import_history == True:
                                history = feat.qualifiers[key][0]
                                results = re.findall("DNA.dna_dict\['[^\[\]]+'\]", history) 
                                for result in results:
                                    _unique_id = result.split("['")[1][:-2] 
                                    DNA.dna_dict[_unique_id] = None
                                history_num = int(key.split("_")[-1]) 
                                pairs.append((feat, history_num, history))    
                        
                if len(pairs) == 0:
                    import_history = False

                if import_history == True:
                    for pair in pairs:
                        feat = pair[0]
                        new_history_num = pair[1] + DNA._num_history
                        feat.qualifiers["building_history_{}".format(new_history_num)] = [pair[2]] 
                        del feat.qualifiers["building_history_{}".format(pair[1])]
                        history_nums.append(new_history_num) 
                else:
                    #process_description = None
                    DNA.process_description = process_description
                    archivehistory(self) 
                
                DNA._num_history = max(history_nums)  
            
            else:
                import_history = False
                self.dnafeatures = []

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
                
                #self.dnafeatures = self.record.features
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
       
        self.subject = SourceDNA()
        self._setfeatureid()
        self._features_dict = dict(list(map(lambda x:(x._id, x), self.dnafeatures)))
        if _direct == 1 and import_history == False:
            if project in DNA.dna_dict:
                keys   = list(DNA.dna_dict.keys())
                unique = 0
                while project + "_" + str(unique) in DNA.dna_dict:
                    unique += 1
                self._unique_id = project + "_" + str(unique)
            else: 
                self._unique_id = project 
        
            DNA.dna_dict[self._unique_id] = None
            args = [seq, recname, project, topology, format, process_description]
            for i in range(len(args)):
                if type(args[i]) is str:
                    args[i] = "'" + args[i] + "'" 

            edit_history = "DNA.dna_dict['{}'] = DNA(seq={}, record={}, project={}, topology={}, format={}, process_description={})".format(self._unique_id, args[0], args[1], args[2], args[3], args[4], args[5]) 
            _add_history(self, edit_history)
        else:
            self._unique_id = project 
        self.record.feartures = _assigndnafeatures(self.dnafeatures)

    def finddna(self, query, key_attribute=None, min_match=None, max_mismatch=0, process_description=None):
        if process_description is None:
            process_description = DNA.process_description
        else:
            DNA.process_description = process_description  
        
        qkey = self._unique_id + "_" + str(DNA._qnum)
        features = editfeature(self, query=query, key_attribute=key_attribute, min_match=min_match, max_mismatch=max_mismatch, target_attribute=None, operation=None, __direct=0, process_description=process_description) 
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
            
            feature._qkey            = qkey
            feature._second_id       = key
            feature._start.qkey      = qkey
            feature._start.parent_id = key
            feature._end.qkey        = qkey
            feature._end.parent_id   = key
        
        if "_unique_id" not in self.__dict__:
            self._unique_id = self.project

        feature_names = ", ".join(feature_names)
        DNA.queried_features_dict[qkey] = features
        edit_history  = "DNA.queried_features_dict['{}'] = DNA.dna_dict['{}'].finddna(query='{}', key_attribute='{}', min_match={}, max_mismatch={}, process_description='{}')".format(qkey, self._unique_id, query, key_attribute, min_match, max_mismatch, process_description) 
        _add_history(self, edit_history)
        DNA._qnum += 1 
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
        if (type(other) == str and set(other) <= set("ATGCRYKMSWBDHVNatgcrykmswbdhvn")) or type(other) == Seq:
            other = DNA(seq=other) 

        elif type(other) == SeqRecord:
            other = DNA(record=other) 

        elif type(other) == DNA:
            pass 
        
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized DNA object cannot be joined with others.") 
        else:
            return joindna(self, other, __direct=0)

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
            strand = 0  
            seq_rc = self.seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
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

            if strand is None:
                strand = None
            
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
            return top
        elif strand == -1:
            return bottom[::-1]
        else:
            return top, bottom[::-1]
    
    def _setfeatureid(self):
        for i in range(0, len(self.dnafeatures)):
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

    def writedna(self, handle, format="genbank", record_id=None):
        features = copy.deepcopy(self.dnafeatures)
        for feat in features:
            if "broken_feature" in feat.qualifiers:
                note   = feat.qualifiers["broken_feature"][0]
                label  = note.split("]")[0] + "]"
                length = int(note.split("]")[0].split(":")[2]) 
                note   = note.split("]")[1] 
                pos_s  = int(note.split(":")[1].split("..")[0].replace(" ",""))
                pos_e  = int(note.split(":")[1].split("..")[1].replace(" ",""))
                
                if pos_e > length:
                    note = label + ":" + str(pos_s) + ".." + str(pos_e-length)
                    feat.qualifiers["broken_feature"] = [note]

                elif (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["broken_feature"]
        
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
