import os 
import regex as re
import io 
import sys
import copy
import argparse
import collections
import configparser
from Bio import SeqIO
from Bio import Alphabet
from Bio import Restriction  
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

amb     = IUPACAmbiguousDNA()
#re_dict = Restriction.__dict__
class matchdna():
    def __init__(self):
        start  = None
        end    = None
        qseq   = None
        sseq   = None
        qseq_origin = None

def finddna(subject, query, mode="normal", option={"s":0,"d":0,"i":0,"c":""}, blast_option={"task":"blastn_short"}):
    if mode == "normal":
        if dna.topology == "circular": 
            subject = dna.seq + dna.seq[0:len(query)-1]
        constant    = ""
        query_regex = "({}){{s<={},i<={},d<={}}}".format(option["s"], option["d"], option["i"])
        regex       = re.compile(query_regex) 
        re.search(query_regex, subject )　
        
    elif mode == "regex":
        

    elif mode == "blast":

    
def finddna(dna, query, strand="both", min_match=1):
    if dna.topology == "circular": 
        subject= = dna.seq + dna.seq[0:len(query)-1]
    regex  = re.compile(query)
    for i in range(len(subject)-len(query)):
        sub_subject = subject[i:]
        sub_subject.match = 

def cropdna(dna, start=0, end=0):
    """
    Extract subsequence from the specified region in the genome.
    All features within the subsequence including truncated features are carried over to extracrted 'dna' object. 
    """   
    start = len(dna.seq) + start if start < 0 else start
    end   = end - len(dna.seq) if end > len(dna.seq) else end
    if start >= end:
        if dna.topology == "circular":
            sub_dna1      = cropdna(dna, start, len(dna.seq))
            sub_dna2      = cropdna(dna, 0, end)
            sub_dna       = joindna(sub_dna1, sub_dna2)
            sub_dna.start = start 
            sub_dna.end   = end
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

                if "note_dna.py" not in feat1.qualifiers: 
                    if feat1.type == "source":
                        label = dna.name
                    else:
                        flag = 0
                        for key in ["label", "gene", "product"]:
                            if key in feat1.qualifiers:
                                label = feat1.qualifiers[key][0]
                                flag  = 1
                                break 
                        if flag == 0:
                            label = feat1.type
                    
                    if strand >= 0:
                        feat1.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, 1, len(dna.seq)-s, len(feat1.original)) 
                    else:
                        feat1.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, len(dna.seq)-s, 1, len(feat1.original)) 
                
                else:
                    note   = feat.qualifiers["note_dna.py"]
                    if strand >= 0:
                        label  = note.split(":")[0]
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s, pos_s + len(dna.seq)-s, length)
                    else:
                        label  = note.split(":")[0]
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s, pos_s - (len(dna.seq)-s), length)
                    feat1.qualifiers["note_dna.py"] = note

                
                if "note_dna.py" not in feat2.qualifiers:
                    if feat2.type == "source":
                        label = dna.name
                    else:
                        flag = 0
                        for key in ["label", "gene", "product"]:
                            if key in feat2.qualifiers:
                                label = feat2.qualifiers[key][0]
                                flag  = 1
                                break 
                        if flag == 0:
                            label = feat2.type

                    if strand >= 0:
                        feat2.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, len(dna.seq)-s+1, len(dna.seq)-s+e, len(feat2.original)) 
                    else:
                        feat2.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, len(dna.seq)-s+e, len(dna.seq)-s+1, len(feat2.original)) 
                else:
                    note   = feat.qualifiers["note_dna.py"]
                    if strand >= 0:
                        label  = note.split(":")[0]
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s + len(dna.seq)-s, pos_e, length)
                    else:
                        label  = note.split(":")[0]
                        pos_s  = int(note.split(":")[1].split("..")[0]) 
                        pos_e  = int(note.split(":")[1].split("..")[1]) 
                        length = int(note.split(":")[2])
                        note   = "{}:{}..{}:{}".format(label, pos_s - (len(dna.seq)-s), pos_e, length)
                    feat2.qualifiers["note_dna.py"] = note
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
                        if "note_dna.py" not in feat.qualifiers:
                            if feat.type == "source":
                                label = dna.name
                            else:
                                flag = 0
                                for key in ["label", "gene", "product"]:
                                    if key in feat.qualifiers:
                                        label = feat.qualifiers[key][0]
                                        flag  = 1
                                        break 
                                if flag == 0:
                                    label = feat.type
                            
                            if strand >= 0:
                                feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, e-s, abs(s-start)+1, len(feat.original)) 
                        else:
                            if strand >= 0:
                                note   = feat.qualifiers["note_dna.py"]
                                label  = note.split(":")[0]
                                pos_s  = int(note.split(":")[1].split("..")[0]) + abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            else:
                                note   = feat.qualifiers["note_dna.py"]
                                label  = note.split(":")[0]
                                pos_s  = int(note.split(":")[1].split("..")[0]) - abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            feat.qualifiers["note_dna.py"] = note
                    else:
                        feat.location.parts[0]._start = ExactPosition(s - start) 
                
                    feat.location.parts[-1]._end = ExactPosition(e - start)  
                    if feat.location.parts[-1]._end > end-start:
                        feat.location.parts[-1]._end = ExactPosition(end - start)
                        if "note_dna.py" not in feat.qualifiers:
                            if feat.type == "source":
                                label = dna.name
                            else:
                                flag = 0
                                for key in ["label", "gene", "product"]:
                                    if key in feat.qualifiers:
                                        label = feat.qualifiers[key][0]
                                        flag  = 1
                                        break 
                                if flag == 0:
                                    label = feat.type

                            if strand >= 0: 
                                feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, 1, end-s, len(feat.original)) 
                            else:
                                feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, end-s, 1, len(feat.original)) 
                        else:
                            s = int(feat.location.parts[0].start.position)
                            if strand >= 0:
                                note   = feat.qualifiers["note_dna.py"]
                                label  = note.split(":")[0]
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s + (end-start-s)-1, length)
                            else:
                                note   = feat.qualifiers["note_dna.py"]
                                label  = note.split(":")[0]
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s - (end-start-s)+1, length)
                            feat.qualifiers["note_dna.py"] = note
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
                        if "note_dna.py" not in feat.qualifiers:
                            if feat.type == "source":
                                label = dna.name
                            else:
                                flag = 0
                                for key in ["label", "gene", "product"]:
                                    if key in feat.qualifiers:
                                        label = feat.qualifiers[key][0]
                                        flag  = 1
                                        break 
                                if flag == 0:
                                    label = feat.type

                            if strand >= 0:
                                feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, e-s, abs(s-start)+1, len(feat.original)) 
                        else:
                            if strand >= 0:
                                note   = feat.qualifiers["note_dna.py"]
                                label  = note.split(":")[0]
                                pos_s  = int(note.split(":")[1].split("..")[0]) + abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            else:
                                note   = feat.qualifiers["note_dna.py"]
                                label  = note.split(":")[0]
                                pos_s  = int(note.split(":")[1].split("..")[0]) - abs(s-start) 
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_e, length)
                            feat.qualifiers["note_dna.py"] = note
                    else:
                        locations[0][0] = ExactPosition(s - start)
                    
                    
                    if e > end-start and eflag == 1:
                        locations[-1][1] = ExactPosition(end-start)
                        if "note_dna.py" not in feat.qualifiers:
                            if feat.type == "source":
                                label = dna.name
                            else:
                                flag = 0
                                for key in ["label", "gene", "product"]:
                                    if key in feat.qualifiers:
                                        label = feat.qualifiers[key][0]
                                        flag  = 1
                                        break 
                                if flag == 0:
                                    label = feat.type
                            feat.qualifiers["note_dna.py"] = "{}:{}..{}:{}".format(label, 1, end-s, len(feat.original)) 
                        else:
                            s      = int(locations[0][0])
                            if strand  >= 0:
                                label  = note.split(":")[0]
                                note   = feat.qualifiers["note_dna.py"]
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s + (end-start-s)-1, length)
                            else:
                                label  = note.split(":")[0]
                                note   = feat.qualifiers["note_dna.py"]
                                pos_s  = int(note.split(":")[1].split("..")[0])
                                pos_e  = int(note.split(":")[1].split("..")[1]) 
                                length = int(note.split(":")[2])
                                note   = "{}:{}..{}:{}".format(label, pos_s, pos_s - (end-start-s)+1, length)
                            feat.qualifiers["note_dna.py"] = note
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
        sub_dna.features       = feats
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
        new_rec.features = sub_dna.features
        new_rec.annotations["topology"] = sub_dna.topology
        sub_dna.record = new_rec
        sub_dna.source = dna
        sub_dna.start  = start 
        sub_dna.end    = end 
            

    return sub_dna 


def joindna(*dnas, topology="linear", product_name=None)
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
    
    construct = copy.deepcopy(dnas[0])
    if len(dnas) > 0:
        for dna in dnas[1:]:
            feats   = dna.features 
            if (dna._left_end_top * construct._right_end_bottom == 1 or dna._left_end_bottom * construct._right_end_top == 1):
                if dna._left_end_top == 1:
                    sticky_end = dna._left_end 
                else:
                    sticky_end = construct._right_end
                
                ovhg    = dna._left_end
                new_dna = substr(dna, len(ovhg), len(dna.seq)) 
                if dna._left_end == construct._right_end
                    print("These DNA objects were not able to be joined. Please check sticky end seqeunces of the DNA objects")
                    return False
                else:
                    print("The DNA objects were joined based on complementary sticky end of each fragment. The sticky end is '{}'".format(sticky_end)) 
            else:
                new_dna = dna
            
            ovhg = ""
            feats  = slide(feats, len(construct.seq) - len(ovhg))
            feats1 = [feat for feat in construct.features if "note_dna.py" in feat.qualifiers]
            feats2 = [feat for feat in feats if "note_dna.py" in feat.qualifiers]
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
                                if key == "note_dna.py":
                                    pass 
                                elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                                    flag = 1
                                else:
                                    flag = 0
                                    break    
                            
                            if flag == 1:
                                note1   = feat1.qualifiers["note_dna.py"]
                                label   = note1.split(":")[0] 
                                pos_s1  = int(note1.split(":")[1].split("..")[0])
                                pos_e1  = int(note1.split(":")[1].split("..")[1]) 
                                length1 = int(note1.split(":")[2])
                                
                                note2   = feat2.qualifiers["note_dna.py"] 
                                pos_s2  = int(note2.split(":")[1].split("..")[0])
                                pos_e2  = int(note2.split(":")[1].split("..")[1]) 
                                length2 = int(note2.split(":")[2])
                                if length1 == length2 and "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original and s2 == 0 and feat1.location.strand == feat2.location.strand:
                                    note     = "{}:{}..{}:{}".format(label, pos_s1, pos_e2, length1)
                                    new_seq  = construct.seq[s1:e1] + dna.seq[s2:e2] 
                                    new_feat = copy.deepcopy(construct.features[construct.features.index(feat1)]) 
                                    strand = new_feat.location.strand
                                    if len(feat1.location.parts) == 1 and len(feat2.location.parts) == 1:
                                        new_feat.location = FeatureLocation(feat1.location.parts[0].start.position, feat2.location.parts[-1].end.position, feat1.strand)
                                        new_feat.location.strand = strand
                                    else:
                                        locations = feat1.location.parts[0:-1] + [FeatureLocation(feat1.location.parts[-1].start.position, feat2.location.parts[0].end.position, feat1.strand)] + feat2.location.parts[0:-1]
                                        if strand == -1:
                                            locations.reverse() 
                                        new_feat.location = CompoundLocation(locations) 
                                        new_feat.locaiton.strand = strand 
                                    
                                    if strand == -1:
                                        s = new_feat.location.parts[-1].start.position
                                        e = new_feat.location.parts[0].end.position
                                    else:
                                        s = new_feat.location.parts[0].start.position
                                        e = new_feat.location.parts[-1].end.position
                                    
                                    if len(new_seq) - len(ovhg) == e - s and len(new_seq) - len(ovhg) <= len(feat1.original):
                                        construct.features[construct.features.index(feat1)].qualifiers["note_dna.py"] = note
                                        construct.features[construct.features.index(feat1)].location = new_feat.location
                                        del feats[feats.index(feat2)] 
            
            construct.features = construct.features + feats
            construct.seq               = construct.seq + new_dna.seq 
            construct._right_end         = dna._right_end
            construct._right_end_top     = dna._right_end_top
            construct._right_end_bottom  = dna._right_end_bottom
            construct.source            = dna.source
            construct.topology      = "linear"
            ovhg = dna._right_end
    
    construct.features.sort(key=lambda x:x.location.parts[0].start.position)
    for feat in construct.features:
        if "note_dna.py" in feat.qualifiers:
             note = feat.qualifiers["note_dna.py"]
             pos_s  = int(note.split(":")[1].split("..")[0])
             pos_e  = int(note.split(":")[1].split("..")[1])
             length = int(note.split(":")[2])
             if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                del feat.qualifiers["note_dna.py"]
    new_record          = SeqRecord(Seq(str(construct.seq),Alphabet.DNAAlphabet()))
    new_record.features = construct.features
    new_record.annotations["topology"] = topology
    construct.record = new_record     
    if name is None:
        pass
    else:
        new_record.id = name
        construct.name = name 
    
    if topology == "circular":
        construct = _circularizedna(construct)
       
    construct.start = 0 
    construct.end = len(construct.seq) 
    construct._setfeatureid() #Update feature id
    return construct

def modifyends(dna, left="", right=""):
    """
    Set end sequence structures. 
    """
    def check_endseq(top,bottom):
        new_top    = ""
        new_bottom = ""
        for t,b in zip(top,bottom):
            if t != b.translate(str.maketrans("ATGC","TACG")) and (t != "-" and b != "-"):
                return False, False
            new_top    += t
            new_bottom += b
        return new_top, new_bottom
    
    left, right = left.upper(), right.upper() 
    pattern1, pattern2, patternl1, patternl2, patternr1, patternr2 = "[ATGCN*-]*/?[ATGCN*-]*", "[ATGCN*]+-+[ATGCN*]+", "^[ATGCN*]+-+/", "/[ATGCN*]+-+$", "^-+[ATGCN*]+/", "/-+[ATGCN*]+$" 
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
        raise ValueError("Please sepcify a proper sequence pattern for the 'left' argument") 
    
    if pattern1.fullmatch(right_end) != None and pattern2.search(right_end) is None and patternr1.search(right_end) is None and patternr2.search(right_end) is None:
        pass
    else:
        raise ValueError("Please sepcify a proper sequence pattern for the 'right' argument") 
    
    if "/" in left_end:
        left_end_top, left_end_bottom = left_end.split("/")
    else:
        left_end_top, left_end_bottom = left_end, left_end.translate(str.maketrans("ATGC","TACG"))

    if len(left_end_top) != len(left_end_bottom):
        raise ValueError("The length of top and bottom sequence should be same.")
    
    elif "-" in left_end:
        left_end_top, left_end_bottom = check_endseq(left_end_top, left_end_bottom)
        if left_end_top != False:
            left_end_length = len(left_end_top)
            if "*" in left_end_top or "*" in left_end_bottom:
                if set(left_end_top) <= set(["*","-"]) and set(left_end_bottom) <= set(["*","-"]):
                    left_end_top    = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[:left_end_length], left_end_top)])
                    left_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[:left_end_length].translate(str.maketrans("ATGC","TACG")), left_end_bottom)])
                else:
                    print("'*' cannot be used wih 'ATGCN'.")
            
            left_end_top, left_end_bottom = left_end_top.replace("-",""), left_end_bottom.replace("-","") 
            if len(left_end_top) < len(left_end_bottom):
                left_length     = len(left_end_bottom)
                left_end        = left_end_bottom[0:len(left_end_bottom)-1*len(left_end_top)].translate(str.maketrans("ATGC","TACG")) 
                left_end_top    = -1
                left_end_bottom = 1

            else:
                left_length     = len(left_end_top)
                left_end        = left_end_top[0:len(left_end_top)-1*len(left_end_bottom)]
                left_end_top    = 1
                left_end_bottom = -1

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
        right_end_top, right_end_bottom = right_end, right_end.translate(str.maketrans("ATGC","TACG"))

    if len(right_end_top) != len(right_end_bottom):
        raise ValueError("The length of top and bottom sequence should be same.")
    
    elif "-" in right_end:
        right_end_top, right_end_bottom = check_endseq(right_end_top, right_end_bottom)
        if right_end_top != False:
            right_end_length = len(right_end_top)
            if "*" in right_end_top or "*" in right_end_bottom:
                if set(right_end_top) <= set(["*","-"]) and set(right_end_bottom) <= set(["*","-"]):
                    right_end_top    = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[-1*right_end_length:], right_end_top)])
                    right_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(dna.seq[-1*right_end_length:].translate(str.maketrans("ATGC","TACG")), right_end_bottom)])
                else:
                    print("'*' cannot be used wih 'ATGCN'.")

            right_end_top, right_end_bottom = right_end_top.replace("-",""), right_end_bottom.replace("-","") 
            if len(right_end_top) < len(right_end_bottom):
                right_length     = len(right_end_bottom)
                right_end        = right_end_bottom[len(right_end_top):].translate(str.maketrans("ATGC","TACG")) 
                right_end_top    = -1
                right_end_bottom = 1

            else:
                right_length     = len(right_end_top)
                right_end        = right_end_top[len(right_end_bottom):] 
                right_end_top    = 1
                right_end_bottom = -1

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
        dna     = joindna(left_end, dna, right_end, ovhg_check=False)
        dna._left_end  = left_end.seq
        dna._right_end = right_end.seq
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
    if dna.topology == "linear":
        pass
    else:
        positions.append(positions[-1])

    for i, pos in enumerate(positions[0:-1]):
        dnas.append(cropdna(pos,positions[i+1])) 
    
    return dnas

def flipdna(dna):
    """
    Return reverse complement sequence. 
    All feature infomation is also reversed with the sequence.
    """
    if type(dna) == DNA:
        dna = copy.deepcopy(dna)
        seq  = dna.seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        feats = [] 
        for feat in dna.record.features:
            strand = feat.location.strand
            for p in range(len(feat.location.parts)):
                s, e = feat.location.parts[p].start.position, feat.location.parts[p].end.position
                feat.location.parts[p]._start = ExactPosition(len(dna.seq) - e) 
                feat.location.parts[p]._end   = ExactPosition(len(dna.seq) - s) 
            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
                if "note_dna.py" in feat.qualifiers:
                    note   = feat.qualifiers["note_dna.py"]
                    label  = note.split(":")[0]
                    pos_s  = int(note.split(":")[1].split("..")[0])
                    pos_e  = int(note.split(":")[1].split("..")[1])
                    length = int(note.split(":")[2])
                    note = "{}:{}..{}:{}".format(label, pos_e, pos_s, length)
                    feat.qualifiers["note_dna.py"] = note
            else:
                feat.location.strand = strand
            feats.append(feat)
    
    elif type(dna) == str or type(dna) == Seq:
        seq   = str(dna) 
        dna  = DNA(seq=seq) 
        seq   = seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        feats = [] 

    comp = DNA(seq=seq) 
    comp.features = feats
    comp.record.features = feats 
    comp._right_end, comp._left_end = dna._left_end.translate(str.maketrans("ATGC","TACG"))[::-1], dna._right_end.translate(str.maketrans("ATGC","TACG"))[::-1]
    comp._right_end_top, comp._left_end_bottom = dna._left_end_bottom, dna._right_end_top
    comp._right_end_bottom, comp._left_end_top = dna._left_end_top, dna._right_end_bottom
    return comp
    
def _circularize(dna):
    dna = copy.deepcopy(dna)
    seq_origin = dna.seq
    feats_origin = dna.features
    if dna.topology == "circular" and dna.record.annotations["topology"] == "circular":
        print("The DNA object is already circularized")

    if (dna._right_end_top * dna._left_end_bottom == 1 and dna._right_end_bottom * dna._left_end_top == 1) and len(dna._right_end) > 0:
        if str(dna._right_end) == str(dna._left_end) 
            print("The DNA object was circularized based on complementary sticky end between 3' end and 5' end. The sticky end is '{}'".format(dna._left_end)) 
            ovhg = dna._right_end
            subdna = substr(dna, 0, len(dna.seq)-len(dna._right_end))
            dna.seq      = subdna.seq
            dna.record   = subdna.record
        else:
            pass 
    else:
        ovhg = ""

    remove_list = [] 
    feats1      = [feat for feat in dna.features if "note_dna.py" in feat.qualifiers]
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
                    if key == "note_dna.py":
                        pass 
                    elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                        flag = 1
                    else:
                        flag = 0
                        break    
                if flag == 1:
                    note1   = feat1.qualifiers["note_dna.py"]
                    label   = note1.split(":")[0]
                    pos_s1  = int(note1.split(":")[1].split("..")[0])
                    pos_e1  = int(note1.split(":")[1].split("..")[1]) 
                    length1 = int(note1.split(":")[2])
                    
                    note2   = feat2.qualifiers["note_dna.py"] 
                    pos_s2  = int(note2.split(":")[1].split("..")[0])
                    pos_e2  = int(note2.split(":")[1].split("..")[1]) 
                    length2 = int(note2.split(":")[2])
                    
                    if length1 == length2 and ((s1 == 0 and e2 == len(seq_origin)) or (s2 == 0 and e1 == len(seq_origin))):
                        if "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original:
                            note     = "{}:{}..{}:{}".format(label, pos_s1, pos_e2, length1)
                            new_seq  = seq_origin[s1:e1] + seq_origin[s2:e2] 
                            new_feat = copy.deepcopy(dna.features[dna.features.index(feat1)]) 
                            if s1 == 0:
                                new_feat.location = FeatureLocation(s2, e2 + (e1-len(ovhg)), feat1.strand)
                            else:
                                new_feat.location = FeatureLocation(s1, e1 + (e2-len(ovhg)), feat1.strand)
                            
                            if len(new_seq) - len(ovhg) <= len(feat1.original):
                                new_feat.qualifiers["note_dna.py"] = note
                                if len(new_seq) - len(ovhg) == length1:
                                    del dna.features[dna.features.index(feat1)].qualifiers["note_dna.py"]
                                dna.features[dna.features.index(feat1)].location = new_feat.location
                                dna.features.remove(feat2)  
                                remove_list.append(feat2) 
        
    for i in range(len(dna.features)):    
        if dna.features[i].location.parts[-1].end.position > len(dna.seq):
            if dna.features[i].location.parts[0].start.position >= len(dna.seq):
                strand                    = dna.features[i].location.strand
                dna.features[i].location = FeatureLocation(dna.features[i].location.parts[0].start.position-len(dna.seq),dna.features[i].location.parts[-1].end.position-len(dna.seq))
                dna.features[i].location.strand = strand
            else:
                strand    = dna.features[i].location.strand
                locations = [FeatureLocation(dna.features[i].location.parts[0].start.position,len(dna.seq)), FeatureLocation(0,dna.features[i].location.parts[-1].end.position-len(dna.seq))]
                if strand == -1:
                    locations.reverse() 
                dna.features[i].location = CompoundLocation(locations)
                dna.features[i].location.strand = strand
    
    dna.record.features = dna.features
    dna._left_end  = ""
    dna._left_end_top    = 0 
    dna._left_end_bottom = 0 

    dna._right_end = ""
    dna._right_end_top    = 0 
    dna._right_end_bottom = 0 

    dna.topology = "circular"
    dna.record.seq = Seq(str(dna.seq),Alphabet.DNAAlphabet())
    dna.record.annotations["topology"] = dna.topology
    return dna

class DNA():
    order_qualifiers = ["label","gene","product"]
    min_overlap      = 10 
    max_overlap      = 300
    min_match        = 10
    
    def __init__(self, record=None, seq=None, name="None", topology="linear", format=None):
        if seq is None and record is None:
            self.seq            = None
            self.record         = None
            self.features       = None
            self._right_end     = None 
            self._left_end      = None
            self.topology       = topology
            self.source         = None
            self.name           = name
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1

        elif seq is None:
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
                self.features = record.features 
            else:
                self.features = []
            
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
            
            self.source = None
            self.name = record.id  
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
                    if t != b.translate(str.maketrans("ATGC","TACG")) and (t != "-" and b != "-"):
                        return False, False, False
                    else:
                        new_top    += t
                        new_bottom += b
                        if t == "-":
                            new_seq += b.translate(str.maketrans("ATGC","TACG"))
                        else:
                            new_seq += t
                return new_top, new_bottom, new_seq 
            
            sticky   = False
            pattern1 = "[ATGCNatgcn*-]+/?[ATGCNatgcn*-]*"
            pattern2 = "[ATGCNatgcn]+-+[ATGCNatgcn]+"
            pattern1 = re.compile(pattern1)   
            pattern2 = re.compile(pattern2)   
            if pattern1.fullmatch(seq) != None: 
                if "/" in seq:
                    top, bottom, seq = check_seq(seq)
                    if "-" in top or "-" in bottom:
                        if top !=False and pattern2.search(top) is None and pattern2.search(bottom) is None:
                            sticky = True
                        else:
                            raise TypeError("Invaild sequence pattern was detected.")
                    else: 
                        sticky = False

                self.seq      = str(seq).upper()
                self.record   = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
                self.features = self.record.features
                self.topology = topology
                
                if sticky == True:
                    self.topology   = "linear"
                    self._left_end  = "" 
                    self._right_end = "" 
                    if top[0] == "-":
                        for c, char in enumerate(top):
                            if char != "-":
                                break 
                            else:
                                self._left_end += bottom[c].translate(str.maketrans("ATGC","TACG"))
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
                                self._right_end += bottom[::-1][c].translate(str.maketrans("ATGC","TACG"))
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
                    self.source = None
                    self.name = name
                    self._left_end_top      = 1
                    self._left_end_bottom   = 1
                    self._right_end_top     = 1 
                    self._right_end_bottom  = 1

            else:
                raise TypeError("Invaild sequence pattern was detected.")
        
        self._setfeatureid()
        self._features_dict = dict(list(map(lambda x:(x._id, x), self.features)))
        self.start = 0 
        self.end   = len(self.seq) 
    
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
                subdna = substr(self,start,end)
                             
            else:
                raise TypeError("slice indices must be integers or None or have an __index__ method")
            
            subdna.source = self
            if strand == -1 or strand < 0:
                return reverse_complement(subdna)
            else:
                return subdna 
        
        elif type(item) == int:
            subdna = substr(self,item,item+1) 
            return subdna
        
        elif type(item) == str and item.isdecimal() == True: 
            subdna = substr(self, 0, 0, target=item)
            return subdna
        
        elif type(item) == Seq or (type(item) == str and set(item) <= set("ATGCNatgcn")) or type(item) == DNA: 
            if type(item) == DNA: 
                item = str(item.seq) 
            else:
                item = str(item)
            subdna = substr(self, 0, 0, target=item)
            return subdna
        
        else:
            raise ValueError("Invalid index type was specified.") 
    
    def __add__(self, other):
        if (type(other) == str and set(other) <= set("ATGCNatgcn")) or type(other) == Seq:
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
        if (type(other) == str and set(other) <= set("ATGCNatgcn")) or type(other) == Seq:
            other = DNA(seq=other) 

        elif type(other) == Seq.SeqRecord:
            other = DNA(record=other) 

        elif type(other) == DNA:
            pass 
           
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized DNA object cannot be joined with others.") 
        else:
            return join_dna(other, self) 
    
    def printdnaseq(self, whole=False, end_length=10):
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

        seq_rc = self.seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        if len(self._left_end) > end_length:
            left_length = end_length
        else:
            left_length = len(self._left_end)

        if self._left_end_top == 1:
            left_end_top = self.seq[:end_length]
        else:
            left_end_top = " " * left_length + self.seq[left_length:end_length] 
        
        if self._left_end_bottom == 1:
            left_end_bottom = self.seq[:end_length].translate(str.maketrans("ATGC","TACG")) 
        else:
            left_end_bottom = " " * left_length + self.seq[left_length:end_length].translate(str.maketrans("ATGC","TACG"))
        
        if len(self._right_end) > end_length:
            right_length = end_length
        else:
            right_length = len(self._right_end)

        if self._right_end_top == 1:
            right_end_top = self.seq[len(self.seq)-end_length:]
        else:
            right_end_top = self.seq[len(self.seq)-end_length:len(self.seq)-right_length] + " " * right_length 
        
        if self._right_end_bottom == 1:
            right_end_bottom = self.seq[len(self.seq)-end_length:].translate(str.maketrans("ATGC","TACG"))
        else:
            right_end_bottom = self.seq[len(self.seq)-end_length:len(self.seq)-right_length].translate(str.maketrans("ATGC","TACG")) + "–" * right_length     
        
        if whole == False:
            print("5'-{}...{}-3'".format(left_end_top, right_end_top))
            print("3'-{}...{}-5'".format(left_end_bottom, right_end_bottom))
        else:
            print("5'-{}{}{}-3'".format(left_end_top, self.seq[end_length:len(self.seq)-end_length], right_end_top))
            print("3'-{}{}{}-5'".format(left_end_bottom, self.seq[end_length:len(self.seq)-end_length].translate(str.maketrans("ATGC","TACG")), right_end_bottom))

    def adddnafeature(self, start, end, feature_id=None strand=1, feature_type="misc_feature", qualifiers={}):
        if feature_id is None:
            raise ValueError("Please set unique feature id") 
        
        elif str(feature_id) in list(self._features):
            raise ValueError("The feature id is already used. Please set unique feature id") 
        
            raise TypeError("Invalid type for 'qualifiers'. 'qualifiers' should be specified as 'dict' object.") 

        if type(start) == int and type(end) == int:
            if start > end:
                if self.topology == "linear":
                    raise ValueError("When you add feature to linear DNA object, start value should be larger than end value.")
                else:
                    locations = [FeatureLocation(start, len(self.seq)), FeatureLocation(0, end)]
                    if strand == -1:
                        locations.reverse() 
                    feat = SeqFeature(CompoundLocation(locations), type=feature_type)
                    feat.location.strand = strand
            else:
                feat = SeqFeature(FeatureLocation(start, end, strand=strand), type=feature_type)
        
        elif type(start) == list and type(end) == list:
            locations = [] 
            for s,e in zip(start, end):
                if s > e:
                    if self.topology == "linear":
                        raise ValueError("When you add feature to linear DNA object, start value should be larger than end value.")
                    else:
                        locations.append(FeatureLocation(start, len(self.seq), strand=strand))
                        locations.append(FeatureLocation(0, end, strand=strand))
                else:
                    locations.append(FeatureLocation(start, end, strand=strand))
            if strand == -1:
                locations.reverse() 
            feat = SeqFeature(CompoundLocation(locations, type=feature_type))
            feat.location.strand = strand 

        else:
            raise TypeError("'adddnafeature' arguments must be integers or list of integers.")

        for key, value in qualifiers.items():
            if type(qualifiers[key]) == str or type(qualifiers[key]) == int:
                qualifiers[key] = [value] 
            else:
                pass
                
        feat.qualifiers = qualifiers
        feat._id        = feature_id 
        self.features.append(feat) 
        self._features_dict  = dict(list(map(lambda x:(x._id, x), self.features))) 
        self.features.sort(key=lambda x:x.location.parts[0].start.position)
        self.record.features = self.features

    def removednafeature(self, feature_id_list):
        tmp_features = copy.copy(self.features)
        if type(feature_id_list) is str:
            feature_id_list = [feature_id_list]
        for _id in feature_id_list:
            tmp_features.remove(self._features_dict[_id]) 
        self.features = tmp_features 
        self.record.features = tmp_features
        return feature_id_list
    
    def _setfeatureid(self):
        for i in range(0,self.features):
            feat._id = str(i*100)
    
    def getdnafeature(feature_id):  
        return self._features_dict[str(feature_id] 

    def printdnafeature(self, sep=None, output=None, feature_types=None, detail=False, with_seq=False, zero_index=True, feature_id_list=None):
        _ids        = ["Feature_ID"] 
        labels      = ["Label"]
        types       = ["Type"] 
        starts      = ["Start"] 
        ends        = ["End"] 
        strands     = ["Strand"]
        sequences   = ["Seq"
        others_dict = {}

        if feature_id_list is None:
            feature_id_list = list(self._features_dict.keys())
        else:
            feature_id_list = list(map(int, feature_id_list)) 
        
        feature_id_list.sort()  
        features = [self._features_dict[str(_id] for _id in feature_id_list]
        features.sort(key=lambda x:x.location.parts[0].start.position)
        
        for feat in features:
            if feature_types is None or feat.type in feature_types:
                flag  = 0  
                label_keys = []
                for key in feat.qualifiers:
                    if key == "label" or key == "note_dna.py":
                        pass 
                    elif key not in others_dict:
                        others_dict[key] = [key]

                    if key in  DNA.order_qualifiers:
                        label_keys.append((DNA.order_qualifiers.index(key),key))
                        flag = 1
                
                if flag == 0:
                    label = feat.type
                else:
                    label_keys.sort()
                    label_key = label_keys[0][1]
                    if type(feat.qualifiers[label_key]) == list:
                        label = feat.qualifiers[label_key][0]
                    else:
                        label = feat.qualifiers[label_key]
                
                if "note_dna.py" in feat.qualifiers:
                    note   = feat.qualifiers["note_dna.py"]
                    label  = note.split(":")[0]
                    pos_s  = int(note.split(":")[1].split("..")[0])
                    pos_e  = int(note.split(":")[1].split("..")[1])
                    length = int(note.split(":")[2])
                    if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                        pass 
                    else:
                        label = note

                strand = feat.location.strand
                if strand == -1:
                    start = feat.location.parts[-1].start.position
                    end   = feat.location.parts[0].end.position 
                    seq   = substr(self,start,end).seq.translate(str.maketrans("ATGC","TACG"))[::-1]
                else:
                    start = feat.location.parts[0].start.position
                    end   = feat.location.parts[-1].end.position 
                    seq   = substr(self,start,end).seq
                
                if zero_index == False:
                    start += 1 

                if detail == True:
                    for key in others_dict:
                        if key == "label" or key == "note_dna.py":
                            pass 
                        elif key in feat.qualifiers:
                            if type(feat.qualifiers[key]) == list:
                                others_dict[key].append(":".join(feat.qualifiers[key]))
                            else:
                                others_dict[key].append(feat.qualifiers[key])
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
        
        if detail == True:
            rows  = [_ids, labels, types, starts, ends, strands] + list(others_dict.values()) 
            maxes = [_idmax, labelmax, ftypemax, startmax, endmax, strandmax] + other_maxes
        else:
            rows  = [_ids, labels, types, starts, ends, strands] 
            maxes = [_idmax, labelmax, ftypemax, startmax, endmax, strandmax] 

        if with_seq == True:
            rows.append(sequences) 
            maxes.append(sequencemax) 

        for n, row in enumerate(zip(*rows)):
            if type(output) is str:
                output = open(output,"w") 
            
            if sep is None:
                text = ""
                for m, x in enumerate(row):
                    text += x + " " * (maxes[m]-len(x)) 
                print(text, file=output) 
            
            else:
                print(*row, sep=",", file=output)
            
    def reindex(self, zero_position):
        if self.topology == "linear":
            print("The method cannot be used for linear DNA object.")
        else:
            tmp_dna = self[zero_position:] + self[:zero_position] 
            tmp_dna = circularize(adna, ovhg_check=False) 
            for key in tmp_dna.__dict__.keys():
                self.__dict__[key] = tmp_dna.__dict__[key] 
            
    def writedna(self, handle, format="genbank", record_id=None):
        features = copy.deepcopy(self.features)
        for feat in features:
            if "note_dna.py" in feat.qualifiers:
                note = feat.qualifiers["note_dna.py"]
                pos_s  = int(note.split(":")[1].split("..")[0])
                pos_e  = int(note.split(":")[1].split("..")[1])
                length = int(note.split(":")[2])
                if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["note_dna.py"]
        self.record.features = features 
        if record_id is None:
            self.record.id = self.name
        else:
            self.record.id = record_id
        self.record.seq = Seq(str(self.seq),Alphabet.DNAAlphabet()) 
        self.record.features = self.features
        SeqIO.write(self.record, handle, format)
