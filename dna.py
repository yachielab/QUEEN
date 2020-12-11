import os 
import regex as re
import io 
import sys
import copy
import argparse
import collections
import configparser
import Levenshtein
from Bio import SeqIO
from Bio import Alphabet
from Bio import Restriction
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import Bio
import warnings
warnings.simplefilter('ignore', Bio.BiopythonParserWarning)

amb     = IUPACAmbiguousDNA()
#re_dict = Restriction.__dict__
def cropdna(dna, start=0, end=0):
    """
    Extract subsequence from the specified region in the genome.
    All features within the subsequence including truncated features are carried over to extracrted 'dna' object. 
    """   
    start = len(dna.seq) + start if start < 0 else start
    end   = end - len(dna.seq) if end > len(dna.seq) else end
    if start > end:
        if dna.topology == "circular":
            sub_dna1      = cropdna(dna, start, len(dna.seq))
            sub_dna2      = cropdna(dna, 0, end)
            sub_dna       = joindna(sub_dna1, sub_dna2)
            sub_dna.subject = SourceDNA()
            sub_dna.subject = SourceDNA() 
            sub_dna.subject.start   = start 
            sub_dna.subject.end     = end
            sub_dna.subject.project = dna.project 
            sub_dna.subject.dna     = dna  
        else:
            raise ValueError("Start value should be larger than or equal to end value.")
    else:
        feats   = []
        new_features = [] 
        for feat in dna.record.features:
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

                if "crop_info_dna.py" not in feat1.qualifiers: 
                    if "label" in feat1.qualifiers:
                        label = feat1.qualifiers["label"][0]
                    else:
                        label = "N.A"
                    label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                    if strand >= 0:
                        feat1.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, 1, len(dna.seq)-s, len(feat1.original)) 
                    else:
                        feat1.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, len(dna.seq)-s, 1, len(feat1.original)) 
                
                else:
                    note   = feat.qualifiers["crop_info_dna.py"]
                    if strand >= 0:
                        label  = note.split("]")[0] + "]"
                        note   = note.split("]")[1] 
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s, pos_s + len(dna.seq)-s, length)
                    else:
                        label  = note.split("]")[0] + "]"
                        note   = note.split("]")[1] 
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s, pos_s - (len(dna.seq)-s), length)
                    feat1.qualifiers["crop_info_dna.py"] = note
                
                if "crop_info_dna.py" not in feat2.qualifiers:
                    if "label" in feat2.qualifiers:
                        label = feat1.qualifiers["label"][0]
                    else:
                        label = "N.A"
                    label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                    if strand >= 0:
                        feat2.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, len(dna.seq)-s+1, len(dna.seq)-s+e, len(feat2.original)) 
                    else:
                        feat2.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, len(dna.seq)-s+e, len(dna.seq)-s+1, len(feat2.original)) 
                
                else:
                    note   = feat.qualifiers["crop_info_dna.py"]
                    if strand >= 0:
                        label  = note.split("]")[0] + "]"
                        note   = note.split("]")[1] 
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s + len(dna.seq)-s, pos_e, length)
                    else:
                        label  = note.split("]")[0] + "]"
                        note   = note.split("]")[1] 
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s - (len(dna.seq)-s), pos_e, length)
                    feat2.qualifiers["crop_info_dna.py"] = note
                new_features.append(feat1)
                new_features.append(feat2)
            else:
                new_features.append(feat) 

        for feat in new_features:
            strand = feat.location.strand
            if strand == -1:
                s = feat.location.parts[-1].start.position
                e = feat.location.parts[0].end.position
            else:
                s = feat.location.parts[0].start.position
                e = feat.location.parts[-1].end.position 
            
            if len(feat.location.parts) == 1 and s <= e:
                if e > start and s < end:
                    feat = copy.deepcopy(feat) 
                    if s - start <= 0:
                        feat.location.parts[0]._start = ExactPosition(0)
                        if "crop_info_dna.py" not in feat.qualifiers:
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            if strand >= 0:
                                feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, len(feat.original) - abs(s-start), 1,  len(feat.original)) 
                        else:
                            note   = feat.qualifiers["crop_info_dna.py"]
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
                            feat.qualifiers["crop_info_dna.py"] = note
                    else:
                        feat.location.parts[0]._start = ExactPosition(s - start) 
                
                    feat.location.parts[-1]._end = ExactPosition(e - start)  
                    if feat.location.parts[-1]._end > end-start:
                        feat.location.parts[-1]._end = ExactPosition(end - start)
                        if "crop_info_dna.py" not in feat.qualifiers: 
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            if strand >= 0: 
                                feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, 1, end-s, len(feat.original)) 
                            else:
                                feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, len(feat.original), len(feat.original)-(end-s)+1,  len(feat.original)) 
                        else:
                            s = int(feat.location.parts[0].start.position)
                            note = feat.qualifiers["crop_info_dna.py"]
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
                            feat.qualifiers["crop_info_dna.py"] = note
                    
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
                        feat = copy.deepcopy(feat) 
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
                    if s - start <= 0 and sflag == 1:
                        locations[0][0] = ExactPosition(0)
                        if "crop_info_dna.py" not in feat.qualifiers:
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            if strand >= 0:
                                feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, e-s,  abs(s-start)+1, len(feat.original)) 
                        else:
                            note   = feat.qualifiers["crop_info_dna.py"]
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
                            feat.qualifiers["crop_info_dna.py"] = note
                    else:
                        locations[0][0] = ExactPosition(s - start)
                    
                    if e > end-start and eflag == 1:
                        locations[-1][1] = ExactPosition(end-start)
                        if "crop_info_dna.py" not in feat.qualifiers:
                            if "label" in feat.qualifiers:
                                label = feat.qualifiers["label"][0]
                            else:
                                label = "N.A"
                            label = "[{}]".format("{}:{}:{}..{}".format(dna.project, label, s, e))
                            feat.qualifiers["crop_info_dna.py"] = "{}:{}..{}:{}".format(label, 1, end-s, len(feat.original)) 
                        else:
                            s      = int(locations[0][0])
                            note   = feat.qualifiers["crop_info_dna.py"]
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
                            feat.qualifiers["crop_info_dna.py"] = note
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
        sub_dna     = copy.deepcopy(dna)
        sub_dna.seq = dna.seq[start:end]
        sub_dna.dnafeature       = feats
        sub_dna._features_dict = dict(list(map(lambda x:(x._id, x), sub_dna.dnafeature)))
        sub_dna.topology       = "linear"
        if start == 0:
            sub_dna._left_end          = dna._left_end
            sub_dna._left_end_top      = dna._left_end_top
            sub_dna._left_end_bottom   = dna._left_end_bottom
        else:
            sub_dna._left_end          = sub_dna.seq[0:20] 
            sub_dna._left_end_top      = 1
            sub_dna._left_end_bottom   = 1
        
        if end == len(dna.seq):
            sub_dna._right_end          = dna._right_end
            sub_dna._right_end_top      = dna._right_end_top
            sub_dna._right_end_bottom   = dna._right_end_bottom

        else:
            sub_dna._right_end          = sub_dna.seq[-20:]
            sub_dna._right_end_top      = 1
            sub_dna._right_end_bottom   = 1

        new_rec          = copy.deepcopy(dna.record)
        new_rec.seq      = Seq(dna.seq[start:end],Alphabet.DNAAlphabet()) 
        new_rec.features = sub_dna.dnafeature
        new_rec.annotations["topology"] = sub_dna.topology
        sub_dna.record = new_rec
        sub_dna.subject = SourceDNA()
        sub_dna.subject.strand  = 1
        sub_dna.subject.start   = start 
        sub_dna.subject.end     = end
        sub_dna.subject.project = dna.project 
        sub_dna.subject.dna = dna 
    return sub_dna 


def joindna(*dnas, topology="linear", project=None):
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
                new_dna = cropdna(dna, len(ovhg), len(dna.seq)) 
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
            feats1 = [feat for feat in construct.dnafeature if "crop_info_dna.py" in feat.qualifiers]
            feats2 = [feat for feat in feats if "crop_info_dna.py" in feat.qualifiers]
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
                                if key == "crop_info_dna.py":
                                    pass 
                                elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                                    flag = 1
                                else:
                                    flag = 0
                                    break    
                            
                            if flag == 1:
                                note1   = feat1.qualifiers["crop_info_dna.py"]
                                label   = note1.split("]")[0] + "]"
                                note1   = note1.split("]")[1] 
                                pos_s1  = int(note1.split(":")[1].split("..")[0])
                                pos_e1  = int(note1.split(":")[1].split("..")[1]) 
                                length1 = int(note1.split(":")[2])
                                
                                note2   = feat2.qualifiers["crop_info_dna.py"] 
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
                                        construct.dnafeature[construct.dnafeature.index(feat1)].qualifiers["crop_info_dna.py"] = note
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
            if "crop_info_dna.py" in feat.qualifiers:
                 note   = feat.qualifiers["crop_info_dna.py"]
                 label  = note.split("]")[0] + "]"
                 note   = note.split("]")[1] 
                 pos_s  = int(note.split(":")[1].split("..")[0])
                 pos_e  = int(note.split(":")[1].split("..")[1])
                 length = int(note.split(":")[2])
                 if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["crop_info_dna.py"]
        new_record          = SeqRecord(Seq(str(construct.seq),Alphabet.DNAAlphabet()))
        new_record.features = construct.dnafeature
        new_record.annotations["topology"] = topology
        construct.record = new_record     
        
        if project is None:
            pass
        else:
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
        construct._setfeatureid() #Update feature id
        construct._features_dict = dict(list(map(lambda x:(x._id, x), construct.dnafeature)))
        construct.subject = SourceDNA()
        construct.query   = SourceDNA() 
    else:
        construct = _circularizedna(dnas[0])
    
    return construct

def modifyends(dna, left="", right="", add=0, add_right=0, add_left=0):
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
        left_end  = DNA(seq=left_end) 
        right_end = DNA(seq=right_end) 
        if add_left == 1 and add_right == 1:
            dna = joindna(left_end, dna, right_end)
            dna._left_end  = left_end._left_end
            dna._right_end = right_end._right_end
        
        elif add_right == 1:
            dna = dna[left_end_length-left_length:] + right_end
            dna._right_end = right_end._right_end
            dna._left_end  = left_end.seq
            dna._left_end_top = left_end_top 
            dna._left_end_bottom = left_end_bottom

        else:
            dna = left_end + dna[:len(dna.seq)-right_end_length + right_length]
            dna._left_end = left_end._left_end
            dna._right_end = right_end.seq
            dna._right_end_top = right_end_top 
            dna._right_end_bottom = right_end_bottom

    else:
        left_end  = DNA(seq=left_end) 
        right_end = DNA(seq=right_end)
        dna = dna[left_end_length-left_length:len(dna.seq)-right_end_length + right_length]
        dna._left_end  = left_end.seq
        dna._right_end = right_end.seq
        dna._left_end_top     = left_end_top 
        dna._left_end_bottom  = left_end_bottom
        dna._right_end_top    = right_end_top
        dna._right_end_bottom = right_end_bottom
    
    return dna

def cutdna(dna, *positions):
    dnas = [] 
    new_positions = [] 
    for pos in positions:
        if pos > len(dna.seq):
            pos = pos - len(dna.seq)
        elif pos < 0:
            pos = pos + len(dna.seq) 
        new_positions.append(pos) 

    if 0 not in positions:
        positions = (0,) + positions 
    if len(dna.seq) not in positions:
        positions = positions + (len(dna.seq),) 
    positions = list(positions)
    positions.sort() 

    for i, pos in enumerate(positions[0:-1]):
        dnas.append(cropdna(dna,pos,positions[i+1]))
    
    if dna.topology == "circular":
        new_dna = joindna(dnas[-1],dnas[0]) 
        dnas    = [new_dna] + dnas[1:-1] 
    return dnas

def flipdna(dna):
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
            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
                if "crop_info_dna.py" in feat.qualifiers:
                    note   = feat.qualifiers["crop_info_dna.py"]
                    label  = note.split("]")[0] + "]"
                    note   = note.split("]")[1] 
                    pos_s  = int(note.split(":")[1].split("..")[0])
                    pos_e  = int(note.split(":")[1].split("..")[1])
                    length = int(note.split(":")[2])
                    note = "{}:{}..{}:{}".format(label, pos_e, pos_s, length)
                    feat.qualifiers["crop_info_dna.py"] = note
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
    comp.record.features = feats 
    comp._right_end, comp._left_end = dna._left_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], dna._right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
    comp._right_end_top, comp._left_end_bottom = dna._left_end_bottom, dna._right_end_top
    comp._right_end_bottom, comp._left_end_top = dna._left_end_top, dna._right_end_bottom
    if dna.subject.dna is None:
        pass 
    else:
        comp.subject = dna.subject
        comp.subject.strand = comp.subject.strand * -1 
    return comp
    
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
            subdna = cropdna(dna, 0, len(dna.seq)-len(dna._right_end))
            dna.seq      = subdna.seq
            dna.record   = subdna.record
        else:
            return False
    else:
        ovhg = ""

    remove_list = [] 
    feats1      = [feat for feat in dna.dnafeature if "crop_info_dna.py" in feat.qualifiers]
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
                    if key == "crop_info_dna.py":
                        pass 
                    elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                        flag = 1
                    else:
                        flag = 0
                        break    
                if flag == 1:
                    note1   = feat1.qualifiers["crop_info_dna.py"]
                    label   = note1.split("]")[0] + "]"
                    note1   = note1.split("]")[1] 
                    pos_s1  = int(note1.split(":")[1].split("..")[0])
                    pos_e1  = int(note1.split(":")[1].split("..")[1]) 
                    length1 = int(note1.split(":")[2])
                    
                    note2   = feat2.qualifiers["crop_info_dna.py"] 
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
                                new_feat.qualifiers["crop_info_dna.py"] = note
                                if len(new_seq) - len(ovhg) == length1:
                                    del dna.dnafeature[dna.dnafeature.index(feat1)].qualifiers["crop_info_dna.py"]
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
    
    dna.record.features = dna.dnafeature
    dna._left_end  = ""
    dna._left_end_top    = 0 
    dna._left_end_bottom = 0 

    dna._right_end = ""
    dna._right_end_top    = 0 
    dna._right_end_bottom = 0 

    dna.topology = "circular"
    dna.record.seq = Seq(str(dna.seq),Alphabet.DNAAlphabet())
    dna.record.annotations["topology"] = dna.topology
    dna._features_dict = dict(list(map(lambda x:(x._id, x), dna.dnafeature)))
    return dna

class MatchDNA():
    def __repr__(self):
        return "<dna.MatchDNA object; sseq='{}', sstart={}, send={}, qseq='{}', qstart={}, qend={}, strand={}, match_count={}>".format(self.sseq, self.sstart, self.send, self.qseq, self.qstart, self.qend, self.strand, self.match_count)

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
    def __repr__(self):
        out = "<dna.DNA object; project='{}', sequence length='{} bp', topology='{}'>".format(self.project, len(self.seq), self.topology)
        #out += self.printdnaseq(whole=False, end_length=max([10, len(str(self._left_end)), len(str(self._right_end))]), linebreak=None, display=False)
        return out 

    def __init__(self, seq=None, record=None, project="None", topology="linear", format=None):
        if seq is None and record is None:
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
                record = SeqIO.parse(record,fmt)
                record = next(record) 
            else:
                fmt = None
            self.seq      = str(record.seq).upper()
            self.record   = record
            if len(record.features) > 0:
                self.dnafeature = record.features 
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
                self.project = record.id  
            else:
                self.project = project
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1

        elif record is None:
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
                self.record     = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
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
        flag = 0 
        for feat in self.dnafeature: 
            if feat.type == "source" and feat.location.parts[0].start == 0 and feat.location.parts[-1].end == len(self.seq):
                feat.qualifiers["label"] = [self.project] 
                flag = 1 
                break 
        
        if flag == 0:
            feat = SeqFeature(FeatureLocation(0, len(self.seq), strand=1), type="source") 
            feat.qualifiers["label"] = [self.project] 
            self.dnafeature.append(feat) 
            self.record.feartures = self.dnafeature

        self._setfeatureid()
        self._features_dict = dict(list(map(lambda x:(x._id, x), self.dnafeature)))
        self.subject = SourceDNA()
        self.query   = SourceDNA() 

    def _get_matchlist_normal(query, subject, strand, topology, min_match, error, seq_mismatch, constant):
        match_list      = []
        subject_origin  = subject 
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
                    
                    i = i + result.start() + 1 
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
            else:
                break
        
        if len(match_list) > 0:
            match_list.sort(key=lambda x: x.match_count*-1.0)
            top_match      = match_list[0] 
            new_match_list = [] 
            new_match_list.append(top_match) 
            for match in match_list[1:]:
                check = sum(list(map(lambda x : (x.sspan[0] <= match.sspan[0] and match.sspan[1] <= x.sspan[1] and match.sspan != x.sspan), match_list)))
                if check > 0:
                    pass 
                else:
                    new_match_list.append(match) 
        else:
            new_match_list = match_list 

        return new_match_list
    
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

    def finddna(self, query, attribute=None, min_match=None, max_mismatch=0):
        #attributes = "feature ID", "feature type", "start", "end", "sequence", "query:"
        if type(attribute) is list or type(attribute) is tuple:
            if "seqeunce" in attribute:
                raise ValueError("'sequence' attribute cannot be used with others") 
            else:
                dna_list = [self]
                attributes = attribute
                for attribute in attributes:
                    new_dna_list = []
                    for dna in dna_list:
                        new_dna_list.extend(finddna(dna, query, attribute))
                    dna_list = new_dna_list
            return new_dna_list

        if attribute == "start" and attribute == "end" and attribute == "strand":
            raise TypeError("Numerical value cannot be specified in attrribute on 'finddna.'") 
        
        dna_list = [] 
        if attribute is None:
            attribute = "sequence"
        
        if attribute == "sequence":
            match_list = [] 
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
                match_list = DNA._get_matchlist_normal(query, self.seq, 1, self.topology, min_match, error, seq_mismatch, constant) 
                match_list.extend(DNA._get_matchlist_normal(query, self.seq, -1, self.topology, min_match, error, seq_mismatch, constant)) 
                match_list.sort(key=lambda x: x.score*-1.0)
                querydna = DNA(seq=query) 

            elif mode == "regex":
                constant = "" 
                match_list = DNA._get_matchlist_regex(query, self.seq, 1, self.topology, constant) 
                match_list.extend(DNA._get_matchlist_regex(query, self.seq, -1, self.topology, constant))
            
            for match in match_list: 
                dna = cropdna(self, match.sstart, match.send) 
                if match.strand == -1:
                    dna = flipdna(dna) 
                
                if mode == "normal":
                    dna.query.start  = match.qstart
                    dna.query.end    = match.qend
                    dna.query.strand = 1 
                    dna.query.dna    = querydna
                
                dna_list.append(dna)
            del match_list

        elif attribute == "feature type" or attribute == "feature_type":
            query = re.compile(query)
            for feat in self.dnafeature:
                if query.fullmatch(feat.type) != None:
                    strand = feat.location.strand
                    if strand == -1:
                        start = feat.location.parts[-1].start.position
                        end   = feat.location.parts[0].end.position 
                    else:
                        start = feat.location.parts[0].start.position
                        end   = feat.location.parts[-1].end.position 
                    
                    if start == 0 and end == len(self.seq):
                        dna_list.append(self) 
                    else:
                        if strand == -1:
                            dna = flipdna(cropdna(self,start,end))
                        else:
                            dna = cropdna(self,start,end)
                        dna_list.append(dna) 

        elif attribute == "feature ID" or attribute == "feature_ID":
            query = re.compile(str(query))
            for key, feat in self._features_dict.items():
                if query.fullmatch(key) != None:
                    strand = feat.location.strand
                    if strand == -1:
                        start = feat.location.parts[-1].start.position
                        end   = feat.location.parts[0].end.position 
                    else:
                        start = feat.location.parts[0].start.position
                        end   = feat.location.parts[-1].end.position 
                
                    if start == 0 and end == len(self.seq):
                        dna_list.append(self) 
                    else:
                        if strand == -1:
                            dna = flipdna(cropdna(self,start,end))
                        else:
                            dna = cropdna(self,start,end)
                        dna_list.append(dna) 
        
        #elif attribute == "strand":
        #    for feat in self.dnafeature:
        #        strand = feat.location.strand
        #        if strand == query:
        #            if strand == -1:
        #                start = feat.location.parts[-1].start.position
        #                end   = feat.location.parts[0].end.position 
        #                dna   = flipdna(cropdna(self,start,end))
        #            else:
        #                start = feat.location.parts[0].start.position
        #                end   = feat.location.parts[-1].end.position 
        #                dna   = cropdna(self,start,end) 
        #            dna_list.append(dna) 
        
        elif attribute[0:len("qualifier:")] == "qualifier:":
            key   = attribute.split(":")[-1] 
            query = re.compile(query)
            for feat in self.dnafeature:
                if key in feat.qualifiers:
                    if type(feat.qualifiers[key]) is list:
                        pass 
                    else:
                        feat.qualifiers[key] = [feat.qualifiers[key]]

                for element in feat.qualifiers[key]:
                    if query.fullmatch(element) != None:
                        strand = feat.location.strand
                        if strand == -1:
                            start = feat.location.parts[-1].start.position
                            end   = feat.location.parts[0].end.position 
                        else:
                            start = feat.location.parts[0].start.position
                            end   = feat.location.parts[-1].end.position 
                    
                        if start == 0 and end == len(self.seq):
                            dna_list.append(self) 
                        else:
                            if strand == -1:
                                dna = flipdna(cropdna(self,start,end))
                            else:
                                dna = cropdna(self,start,end)
                            dna_list.append(dna) 
                        break
        else:
            raise ValueError("Invalid attribute was specified") 
        
        return dna_list 

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
                subdna = cropdna(self,start,end)
                             
            else:
                raise TypeError("slice indices must be integers or None or have an __index__ method")
            
            if strand == -1 or strand < 0:
                return reverse_complement(subdna)
            else:
                return subdna 
        
        elif type(item) == int:
            subdna = cropdna(self,item,item+1) 
            return subdna
        
        elif type(item) == str and item.isdecimal() == True: 
            subdna = cropdna(self, 0, 0, target=item)
            return subdna
        
        elif type(item) == Seq or (type(item) == str and set(item) <= set("ATGCRYKMSWBDHVNatgcrykmswbdhvn")) or type(item) == DNA: 
            if type(item) == DNA: 
                item = str(item.seq) 
            else:
                item = str(item)
            subdna = cropdna(self, 0, 0, target=item)
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
    
    def printdnaseq(self, whole=True, end_length=10, linebreak=None, display=True):
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
        
        out = ""
        if whole == False:
            if display == True:
                print("5' {}...{} 3'".format(left_end_top, right_end_top))
                print("3' {}...{} 5'".format(left_end_bottom, right_end_bottom))
            out += "5' {}...{} 3'\n".format(left_end_top, right_end_top)
            out += "3' {}...{} 5'".format(left_end_bottom, right_end_bottom)
        else:    
            top = left_end_top + self.seq[end_length:len(self.seq)-end_length] + right_end_top
            bottom = left_end_bottom + self.seq[end_length:len(self.seq)-end_length].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")) + right_end_bottom 
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
        if display == True:
            print()
            return None
        else:
            return out 
        
    def editfeature(self, feature_id=None, start=None, end=None, strand=1, feature_type="misc_featurte", seq=None, qualifier=None, overwrite_qualifier=True):
        if str(feature_id) not in list(self._features_dict):
            if start is None or end is None:
                raise ValueError("The fearure_id was not found. If you add new featrure, please specify 'start' and 'end'")

            if len(self.dnafeature) > 0:
                if type(start) is list or type(start) is tuple:
                    s = start[0] 
                self.dnafeature.sort(key=lambda x:x.location.parts[0].start.position)
                for feat in self.dnafeature:
                    if feat.location.strand == 1:
                        if s >= feat.location.parts[0].end.position:
                            feature_id = str(int(feat._id) + 1)
                    else:
                        if s >= feat.location.parts[-1].end.position:
                            feature_id = str(int(feat._id) + 1)

                if feature_id is None:
                    feature_id = str(int(feat._id) + 1) 
            else:
                feature_id = str(0)

            if (type(start) is list or type(start) is tuple) or (type(end) is list or type(end) is tuple):
                if len(start) != len(end):
                    raise ValueError("If 'start' and 'end' are specified as list objects, both length should be same.")

                locations = [] 
                for s,e in zip(start, end):
                    if s > e:
                        if self.topology == "linear":
                            raise ValueError("When you add feature to linear DNA object, start value should be larger than or equal to end value.")
                        else:
                            locations.append(FeatureLocation(start, len(self.seq), strand=strand))
                            locations.append(FeatureLocation(0, end, strand=strand))
                    else:
                        locations.append(FeatureLocation(start, end, strand=strand))
                
                if strand == -1:
                    locations.reverse() 

                feat = SeqFeature(CompoundLocation(locations), type=feature_key)
                feat = FeatureLocation(start, end, strand=strand)
                
            else:
                feat = FeatureLocation(start, end, strand=strand)
            
            feat.id = feature_id
            self.dnafeature.append(feat) 
            self._features_dict  = dict(list(map(lambda x:(x._id, x), self.dnafeature))) 
            self.dnafeature.sort(key=lambda x:x.location.parts[0].start.position)
            self.record.features = self.dnafeature

        else:
            if start is not None or end is not None:
                warnings.warn("'start', 'end' and 'strand' values are ignored")
            
            feat   = self._features_dict[feature_id] 
            strand = feat.location.strand
            if strand == -1:
                start        = feat.location.parts[-1].start.position
                end          = feat.location.parts[0].end.position 
                original_seq = cropdna(self,start,end).seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
            else:
                start        = feat.location.parts[0].start.position
                end          = feat.location.parts[-1].end.position 
                original_seq = cropdna(self,start,end).seq

            if seq is not None:
                warnings.warn("When you edit the sequence of a feature, we recommend you add a record of the edit in the feature's qualifiers.")
                if len(seq) != original_seq:
                    raise ValueError("'editfeature' can execute only substitution for the sequnence on the feature. If you add insertion or deletion, please use 'cutdna' and 'joindna.'") 
                else:
                    new_seq  = cropdna(self,0,start).seq + seq + cropdna(self, end, len(self.seq)) 
                    self.seq = new_seq 
            else: 
                seq = original_seq
            
            if qualifier is None:
                qualifier = {} 

            for key, value in qualifier.items():
                if key not in feat.qualifiers() or overwrite_qualifier is True:
                    if type(qualifiers[key]) is list or type(qualifiers[key]) is tuple: 
                        quarifier[key] = list(qualifier[key]) 
                    else: 
                        quarifier[key] = [str(qualifier[key])] 
                    feat.qualifiers[key] = quarifier[key] 
                
                else:
                    if type(qualifier[key]) is list or type(qualifier[key]) is tuple: 
                        feat.qualifiers[key].extend(list(quarifier[key])) 
                    else: 
                        feat.qualifiers[key].append(str(qualifier[key]))
            
            self._features_dict[feature_id] = feat
            self.dnafeature = list(self._features_dict.values()) 
            self.dnafeature.sort(key=lambda x:x.location.parts[0].start.position)
            self.record.features = self.dnafeature

    def removednafeature(self, feature_id_list):
        tmp_features = copy.copy(self.dnafeature)
        if type(feature_id_list) is str:
            feature_id_list = [feature_id_list]
        for _id in feature_id_list:
            tmp_features.remove(self._features_dict[_id]) 
            #print("Removed {}".format(_id))
        self.dnafeature        = tmp_features 
        self._features_dict  = dict(list(map(lambda x:(x._id, x), self.dnafeature))) 
        self.record.features = tmp_features
        return feature_id_list
    
    def _setfeatureid(self):
        for i in range(0, len(self.dnafeature)):
            self.dnafeature[i]._id = str(i*100)
    
    def getdnafeature(feature_id):
        return self._features_dict[str(feature_id)] 

    def printfeature(self, separation=None, output=None, feature_type=None, feature_id=None, attribute=["feature ID", "feature type", "qualifier:label", "start position", "end position"], detail=False, zero_based_index=True):
        #“feature ID,” “qualifier:label,” “feature type,” “start position,” “end position,” and “strand"] 
        _ids        = ["feature ID"] 
        types       = ["feature type"] 
        labels      = ["qualifier:label"]
        starts      = ["start position"] 
        ends        = ["end position"] 
        strands     = ["strand"]
        sequences   = ["sequence"]
        sep         = separation
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
                    seq   = cropdna(self,start,end).seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                else:
                    start = feat.location.parts[0].start.position
                    end   = feat.location.parts[-1].end.position 
                    seq   = cropdna(self,start,end).seq
                
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
            if "crop_info_dna.py" in feat.qualifiers:
                note   = feat.qualifiers["crop_info_dna.py"]
                label  = note.split("]")[0] + "]"
                note   = note.split("]")[1] 
                pos_s  = int(note.split(":")[1].split("..")[0])
                pos_e  = int(note.split(":")[1].split("..")[1])
                length = int(note.split(":")[2])
                
                if pos_e > length:
                    note = label + ":" + str(pos_s) + ".." + str(pos_e-length) + ":" + str(length)
                    feat.qualifiers["crop_info_dna.py"] = note

                elif (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["crop_info_dna.py"]
        
        self.record.features = features 
        if record_id is None:
            self.record.id = self.project
        else:
            self.record.id = record_id
        self.record.seq = Seq(str(self.seq),Alphabet.DNAAlphabet()) 
        SeqIO.write(self.record, handle, format)
        self.record.features = self.dnafeature
