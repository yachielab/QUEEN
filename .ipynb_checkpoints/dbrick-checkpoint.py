import os 
import re
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
re_dict = Restriction.__dict__

def substr(brick, start=0, end=0, info=None):
    """
    Extract subsequence from the specified region in the genome.
    All features within the subsequence including truncated features are also extracted.
    The information that which regions were extracted is appended to the qualifiers of truncated features.
    """
    if info != None:
        if type(info) == list or type(info) == tuple:
            if (type(info[0]) == str and info[0].isdecimal()) or type(info[0]) == int:
                feature_s = brick.features[int(info[0])]
                feature_e = brick.features[int(info[-1])]
                strand    = feature_s.location.strand
                if strand == -1:
                    start = feature_s.location.parts[-1].start.position + start
                else:
                    start = feature_s.location.parts[0].start.position + start

                strand = feature_e.location.strand
                if strand == -1:
                    end = feature_e.location.parts[0].end.position + end
                else:
                    end = feature_e.location.parts[-1].end.position + end
            
            elif type(info[0]) == Seq or (type(info[0]) == str and set(info[0]) <= set("ATGCNatgcn")) or type(info[0]) == Dbrick:
                info = list(info)
                if type(info[0]) == Dbrick:
                    info[0] = str(info[0].seq).upper()
                else:
                    info[0] = str(info[0]).upper() 
                
                if brick.topology == "circular":
                    seq = brick.seq + brick.seq[:len(brick.seq)-1]
                else: 
                    seq = brick.seq
                
                s = seq.find(info[0])
                if s == -1:
                    s = seq.find(info[0].translate(str.maketrans("ATGC","TACG"))[::-1])
                    if s == -1:
                        raise ValueError("The sequence was not found.")
                
                if s > len(brick.seq):
                    s = s - len(brick.seq)
                    
                if type(info[-1]) == Dbrick:
                    info[-1] = str(info[-1].seq).upper()
                else:
                    info[-1] = str(info[-1]).upper()                 

                e = seq.find(info[-1])
                if e == -1:
                    e = seq.find(info[-1].translate(str.maketrans("ATGC","TACG"))[::-1])
                    if e == -1:
                        raise ValueError("The sequence was not found.")
                
                e = e + len(info[-1]) 
                if e > len(brick.seq):
                    e = e - len(brick.seq)
                
                start = s + start
                end   = e + end 
        elif (type(info) == str and info.isdecimal()) or type(info) == int:
            feature = brick.features[int(info)] 
            strand  = feature.location.strand
            if strand == -1:
                start = feature.location.parts[-1].start.position + start
                end   = feature.location.parts[0].end.position + end
            else:
                start = feature.location.parts[0].start.position + start
                end   = feature.location.parts[-1].end.position + end
        
        elif type(info) == Seq or (type(info) == str and set(info) <= set("ATGCNatgcn")) or type(info) == Dbrick:
            if type(info) == Dbrick:
                info = str(info.seq).upper()
            else:
                info = str(info).upper() 
            
            if brick.topology == "circular":
                seq = brick.seq + brick.seq[:len(brick.seq)-1]
            else: 
                seq = brick.seq
            
            s   = seq.find(info)
            if s == -1:
                s = seq.find(info.translate(str.maketrans("ATGC","TACG"))[::-1])
                if s == -1:
                    raise ValueError("The sequence was not found.")
            
            if s > len(brick.seq):
                s = s - len(brick.seq)
            
            start = s + start
            end   = s + len(info) + end
        else:
            raise TypeError("Invalid value for 'info'. 'info' should be specified as 'str' or 'list' object.")
        
        if brick.topology == "linear":
            start = 0 if start < 0 else start
            end   = len(brick.seq) if end > len(brick.seq) else end
        else:
            start = start + len(brick.seq) if start < 0 else start
            end   = end - len(brick.seq) if end > len(brick.seq) else end
            

    if (type(start) == str or type(start) == Seq or type(start) == Dbrick) and (type(end) == str or type(end) == Seq or type(end) == Dbrick):
        sub_brick = pcr(brick, start, end)

    elif start > end:
        if brick.topology == "circular":
            sub_brick1 = substr(brick, start, len(brick.seq))
            sub_brick2 = substr(brick, 0, end)
            sub_brick  = join_dbricks(sub_brick1, sub_brick2, ovhg_check=False)
            sub_brick.start = start 
            sub_brick.end = end
        else:
            raise ValueError("Start value should be larger than or equal to end value.")
    else:
        feats   = []
        new_features = [] 
        for feat in brick.record.features:
            strand = feat.location.strand
            if strand == -1:
                s = feat.location.parts[-1].start.position
                e = feat.location.parts[0].end.position
            else:
                s = feat.location.parts[0].start.position
                e = feat.location.parts[-1].end.position
            if "original" not in feat.__dict__:
                if s > e:
                    feat.original = str(brick.seq)[s:len(brick.seq)] + str(brick.seq)[:e]
                else:
                    feat.original = str(brick.seq)[s:e].upper() 
            if s > e:
                if len(feat.location.parts) == 1:
                    length = len(brick.seq) - s + e
                    locations = [FeatureLocation(s,len(brick.seq)),FeatureLocation(0,e)] 
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
                            new_locations.append(FeatureLocation(part.start.position, len(brick.seq)))   
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

                if "note_dbrick" not in feat1.qualifiers:
                    if strand >= 0:
                        feat1.qualifiers["note_dbrick"] = "{}..{}:{}".format(1, len(brick.seq)-s, len(feat1.original)) 
                    else:
                        feat1.qualifiers["note_dbrick"] = "{}..{}:{}".format(len(brick.seq)-s, 1, len(feat1.original)) 
                else:
                    note   = feat.qualifiers["note_dbrick"]
                    if strand >= 0:
                        pos_s  = int(note.split(":")[0].split("..")[0]) 
                        pos_e  = int(note.split(":")[0].split("..")[1]) 
                        length = int(note.split(":")[1])
                        note   = "{}..{}:{}".format(pos_s, pos_s + len(brick.seq)-s, length)
                    else:
                        pos_s  = int(note.split(":")[0].split("..")[0]) 
                        pos_e  = int(note.split(":")[0].split("..")[1]) 
                        length = int(note.split(":")[1])
                        note   = "{}..{}:{}".format(pos_s, pos_s - (len(brick.seq)-s), length)
                    feat1.qualifiers["note_dbrick"] = note

                
                if "note_dbrick" not in feat2.qualifiers:
                    if strand >= 0:
                        feat2.qualifiers["note_dbrick"] = "{}..{}:{}".format(len(brick.seq)-s+1, len(brick.seq)-s+e, len(feat2.original)) 
                    else:
                        feat2.qualifiers["note_dbrick"] = "{}..{}:{}".format(len(brick.seq)-s+e, len(brick.seq)-s+1, len(feat2.original)) 
                else:
                    note   = feat.qualifiers["note_dbrick"]
                    if strand >= 0:
                        pos_s  = int(note.split(":")[0].split("..")[0]) 
                        pos_e  = int(note.split(":")[0].split("..")[1]) 
                        length = int(note.split(":")[1])
                        note   = "{}..{}:{}".format(pos_s + len(brick.seq)-s, pos_e, length)
                    else:
                        pos_s  = int(note.split(":")[0].split("..")[0]) 
                        pos_e  = int(note.split(":")[0].split("..")[1]) 
                        length = int(note.split(":")[1])
                        note   = "{}..{}:{}".format(pos_s - (len(brick.seq)-s), pos_e, length)
                    feat2.qualifiers["note_dbrick"] = note
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
                        if "note_dbrick" not in feat.qualifiers:
                            if strand >= 0:
                                feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(e-s, abs(s-start)+1, len(feat.original)) 
                        else:
                            if strand >= 0:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0]) + abs(s-start) 
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_e, length)
                            else:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0]) - abs(s-start) 
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_e, length)
                            feat.qualifiers["note_dbrick"] = note
                    else:
                        feat.location.parts[0]._start = ExactPosition(s - start) 
                
                    #s = int(feat.location.parts[0].start.position) 
                    feat.location.parts[-1]._end = ExactPosition(e - start)  
                    if feat.location.parts[-1]._end > end-start:
                        feat.location.parts[-1]._end = ExactPosition(end - start)
                        if "note_dbrick" not in feat.qualifiers:
                            if strand >= 0: 
                                feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(1, end-s, len(feat.original)) 
                            else:
                                feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(end-s, 1, len(feat.original)) 
                        else:
                            s = int(feat.location.parts[0].start.position)
                            if strand >= 0:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0])
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_s + (end-start-s)-1, length)
                            else:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0])
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_s - (end-start-s)+1, length)
                            feat.qualifiers["note_dbrick"] = note
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
                        if "note_dbrick" not in feat.qualifiers:
                            if strand >= 0:
                                feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(abs(s-start)+1, e-s, len(feat.original)) 
                            else:
                                feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(e-s, abs(s-start)+1, len(feat.original)) 
                        else:
                            if strand >= 0:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0]) + abs(s-start) 
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_e, length)
                            else:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0]) - abs(s-start) 
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_e, length)
                            feat.qualifiers["note_dbrick"] = note
                    else:
                        locations[0][0] = ExactPosition(s - start)
                    
                    if e > end-start and eflag == 1:
                        locations[-1][1] = ExactPosition(end-start)
                        if "note_dbrick" not in feat.qualifiers:
                            feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(1, end-s, len(feat.original)) 
                        else:
                            s      = int(locations[0][0])
                            if strand  >= 0:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0])
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_s + (end-start-s)-1, length)
                            else:
                                note   = feat.qualifiers["note_dbrick"]
                                pos_s  = int(note.split(":")[0].split("..")[0])
                                pos_e  = int(note.split(":")[0].split("..")[1]) 
                                length = int(note.split(":")[1])
                                note   = "{}..{}:{}".format(pos_s, pos_s - (end-start-s)+1, length)
                            feat.qualifiers["note_dbrick"] = note
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
        sub_brick     = copy.deepcopy(brick)
        sub_brick.seq = brick.seq[start:end]
        sub_brick.features       = feats
        sub_brick.topology       = "linear"
        if start == 0:
            sub_brick.left_end          = brick.left_end
            sub_brick.left_end_top      = brick.left_end_top
            sub_brick.left_end_bottom   = brick.left_end_bottom
        else:
            sub_brick.left_end          = sub_brick.seq[0:20] 
            sub_brick.left_end_top      = 1
            sub_brick.left_end_bottom   = 1
        
        if end == len(brick.seq):
            sub_brick.right_end          = brick.right_end
            sub_brick.right_end_top      = brick.right_end_top
            sub_brick.right_end_bottom   = brick.right_end_bottom

        else:
            sub_brick.right_end          = sub_brick.seq[-20:]
            sub_brick.right_end_top      = 1
            sub_brick.right_end_bottom   = 1

        new_rec          = copy.deepcopy(brick.record)
        new_rec.seq      = Seq(brick.seq[start:end],Alphabet.DNAAlphabet()) 
        new_rec.features = sub_brick.features
        new_rec.annotations["topology"] = sub_brick.topology
        sub_brick.record  = new_rec
        sub_brick.source  = brick
        sub_brick.start = start 
        sub_brick.end   = end 
    
    return sub_brick 

def join_dbricks(*bricks, topology="linear", name=None, ovhg_check=True, min_overlap=15, max_overlap=150):
    """    
    Join Dbrick objects. 
    When ovhg_check argument is True, adjacent Dbrick objects should have common overhang sequences.  
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
        
    LIG = 0
    GIB = 0 
    ovhg  = bricks[0].right_end
    ovhgs = [] 
    construct = copy.deepcopy(bricks[0])

    if len(bricks) > 0:
        for brick in bricks[1:]:
            feats   = brick.features 
            if ovhg_check == True:
                if brick.left_end == construct.right_end and (brick.left_end_top * construct.right_end_bottom == 1 or brick.left_end_bottom * construct.right_end_top == 1):
                    if LIG == 0:
                        if brick.left_end_top == 1:
                            sticky_end = brick.left_end 
                        else:
                            sticky_end = construct.right_end
                        print("Based on complementary sticky end of each fragment, the Dbrick objects were joined. The sticky end is '{}'".format(sticky_end)) 
                        LIG = 1
                    
                    ovhg = brick.left_end
                    new_brick = substr(brick, len(ovhg), len(brick.seq)) 
                    ovhgs.append(str(ovhg)) 

                elif brick.left_end_top == 1 and brick.left_end_bottom == 1 and brick.right_end_top == 1 and brick.right_end_bottom == 1:
                    if GIB == 0:
                        GIB = 1
                    
                    if brick.left_end == construct.right_end:
                        ovhg = brick.left_end
                        new_brick = substr(brick, len(ovhg), len(brick.seq)) 
                        ovhgs.append(str(ovhg))
                        print("Based on sequence homology at the end of each fragment, the Dbrick objects were joined. The overhang sequence is '{}'".format(ovhg)) 
                    
                    else:
                        flag = 0 
                        for i in range(min_overlap,max_overlap):
                            if brick.seq[:i] == construct.seq[-1*i:]:
                                flag = 1
                                break
                        if flag == 1:
                            brick.left_end = brick.seq[:i]
                            ovhg       = brick.left_end 
                            construct.right_end = brick.left_end 
                            new_brick = substr(brick, len(ovhg), len(brick.seq)) 
                            ovhgs.append(str(ovhg))
                            print("Based on sequence homology at the end of each fragment, the Dbrick objects were joined. The overhang sequence is '{}'".format(ovhg)) 

                        else:
                            print("The adjacent Dbrick objects were not able to be joined. Please set a common overhang sequence between the Dbrick objects.")
                            return False
                else:
                    print(brick.left_end, brick.left_end_top, brick.left_end_bottom) 
                    print(brick.right_end, brick.right_end_top, brick.right_end_bottom)
                    print("The adjacent Dbrick objects were not able to be joined. Please set a common overhang sequence between the Dbrick objects")
                    return False
            else:
                new_brick = brick
                ovhg = ""
            feats  = slide(feats, len(construct.seq) - len(ovhg))
            feats1 = [feat for feat in construct.features if "note_dbrick" in feat.qualifiers]
            feats2 = [feat for feat in feats if "note_dbrick" in feat.qualifiers]
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
                                if key == "note_dbrick":
                                    pass 
                                elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                                    flag = 1
                                else:
                                    flag = 0
                                    break    
                            if flag == 1:
                                note1   = feat1.qualifiers["note_dbrick"]
                                pos_s1  = int(note1.split(":")[0].split("..")[0])
                                pos_e1  = int(note1.split(":")[0].split("..")[1]) 
                                length1 = int(note1.split(":")[1])
                                
                                note2   = feat2.qualifiers["note_dbrick"] 
                                pos_s2  = int(note2.split(":")[0].split("..")[0])
                                pos_e2  = int(note2.split(":")[0].split("..")[1]) 
                                length2 = int(note2.split(":")[1])
                                if length1 == length2 and "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original and s2 == 0 and feat1.location.strand == feat2.location.strand:
                                    note     = "{}..{}:{}".format(pos_s1, pos_e2, length1)
                                    new_seq  = construct.seq[s1:e1] + brick.seq[s2:e2] 
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
                                        construct.features[construct.features.index(feat1)].qualifiers["note_dbrick"] = note
                                        construct.features[construct.features.index(feat1)].location = new_feat.location
                                        del feats[feats.index(feat2)] 
            construct.features = construct.features + feats
            construct.seq               = construct.seq + new_brick.seq 
            construct.right_end         = brick.right_end
            construct.right_end_top     = brick.right_end_top
            construct.right_end_bottom  = brick.right_end_bottom
            construct.source            = brick.source
            construct.topology      = "linear"
            ovhg = brick.right_end
    
    construct.features.sort(key=lambda x:x.location.parts[0].start.position)
    for feat in construct.features:
        if "note_dbrick" in feat.qualifiers:
             note = feat.qualifiers["note_dbrick"]
             pos_s  = int(note.split(":")[0].split("..")[0])
             pos_e  = int(note.split(":")[0].split("..")[1])
             length = int(note.split(":")[1])
             if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                del feat.qualifiers["note_dbrick"]
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
        construct.circularize(ovhg_check=ovhg_check)
    
    construct.start = 0 
    construct.end = len(construct.seq) 
    return construct

def shell(brick, left="", right="", add=0):
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
        if len(left_end_top) != len(left_end_bottom):
            raise ValueError("The length of top and bottom sequence should be same.")
        
        elif "-" in left_end:
            left_end_top, left_end_bottom = check_endseq(left_end_top, left_end_bottom)
            if left_end_top != False:
                left_end_length = len(left_end_top)
                if "*" in left_end_top or "*" in left_end_bottom:
                    if set(left_end_top) <= set(["*","-"]) and set(left_end_bottom) <= set(["*","-"]):
                        left_end_top    = "".join([s if q != "-" else "-" for s,q in zip(brick.seq[:left_end_length], left_end_top)])
                        left_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(brick.seq[:left_end_length].translate(str.maketrans("ATGC","TACG")), left_end_bottom)])
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

    else:
        left_end_length = len(left_end) 
        left_lengh = len(left_end)
        left_end_top    = 1 
        left_end_bottom = 1
    
    if "/" in right_end:
        right_end_top, right_end_bottom = right_end.split("/")
        if len(right_end_top) != len(right_end_bottom):
            raise ValueError("The length of top and bottom sequence should be same.")
        
        elif "-" in right_end:
            right_end_top, right_end_bottom = check_endseq(right_end_top, right_end_bottom)
            if right_end_top != False:
                right_end_length = len(right_end_top)
                if "*" in right_end_top or "*" in right_end_bottom:
                    if set(right_end_top) <= set(["*","-"]) and set(right_end_bottom) <= set(["*","-"]):
                        right_end_top    = "".join([s if q != "-" else "-" for s,q in zip(brick.seq[-1*right_end_length:], right_end_top)])
                        right_end_bottom = "".join([s if q != "-" else "-" for s,q in zip(brick.seq[-1*right_end_length:].translate(str.maketrans("ATGC","TACG")), right_end_bottom)])
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

    else:
        right_end_length = len(right_end)
        right_length = len(right_end)
        rihgt_end_top    = 1 
        right_end_bottom = 1
    
    if add == 1 or (left_end != brick.seq[left_end_length-left_length:left_end_length-left_length+len(left_end)] or right_end != str(brick[len(brick.seq)-right_end_length + right_length - len(right_end):len(brick.seq)-right_end_length + right_length].seq)):
        left_end  = Dbrick(seq=left_end) 
        right_end = Dbrick(seq=right_end) 
        brick     = join_dbricks(left_end, brick, right_end, ovhg_check=False)
        brick.left_end  = left_end.seq
        brick.right_end = right_end.seq
    else:
        left_end  = Dbrick(seq=left_end) 
        right_end = Dbrick(seq=right_end)
        brick = brick[left_end_length-left_length:len(brick.seq)-right_end_length + right_length]
        brick.left_end  = left_end.seq
        brick.right_end = right_end.seq
    
    brick.left_end_top     = left_end_top 
    brick.left_end_bottom  = left_end_bottom
    brick.right_end_top    = right_end_top
    brick.right_end_bottom = right_end_bottom
    return brick

def get_target(bricks, target_seq=None, target_annotation=None):
    """
    Extract a Dbrick object specified by the specific sequence or annotation. 
    """
    if target_seq != None:
        flag = 0 
        target_seq = target_seq.upper()
        for brick in bricks:
            if target_seq in brick.seq or target_seq.translate(str.maketrans("ATGC","TACG"))[::-1] in brick.seq:
                return brick
            else:
                pass 
        print("Target sequence was not found.")
        return False

    elif target_annotation != None:
        for brick in bricks:
            for feat in abrick.features:
                if "label" in feat.qualifiers:
                    if target_annotation in feat.qualifiers["label"][0]:
                        return brick
                else:
                    pass 
        print("Target annotation was not found.")
        return False
    else:
        print("Please specify the arguments 'target_seq' or 'target_annotation'.")
        return False

def reverse_complement(brick):
    """
    Return reverse complement sequence. 
    All feature information is also reversed with the sequence.
    """
    if type(brick) == Dbrick:
        brick = copy.deepcopy(brick)
        seq  = brick.seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        feats = [] 
        for feat in brick.record.features:
            strand = feat.location.strand
            for p in range(len(feat.location.parts)):
                s, e = feat.location.parts[p].start.position, feat.location.parts[p].end.position
                feat.location.parts[p]._start = ExactPosition(len(brick.seq) - e) 
                feat.location.parts[p]._end   = ExactPosition(len(brick.seq) - s) 
            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
                if "note_dbrick" in feat.qualifiers:
                    note = feat.qualifiers["note_dbrick"]
                    pos_s  = int(note.split(":")[0].split("..")[0])
                    pos_e  = int(note.split(":")[0].split("..")[1])
                    length = int(note.split(":")[1])
                    note = "{}..{}:{}".format(pos_e, pos_s, length)
                    feat.qualifiers["note_dbrick"] = note
            else:
                feat.location.strand = strand
            feats.append(feat)
    
    elif type(brick) == str or type(brick) == Seq:
        seq   = str(brick) 
        brick  = Dbrick(seq=seq) 
        seq   = seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        feats = [] 

    comp = Dbrick(seq=seq) 
    comp.features = feats
    comp.record.features = feats 
    comp.right_end, comp.left_end = brick.left_end.translate(str.maketrans("ATGC","TACG"))[::-1], brick.right_end.translate(str.maketrans("ATGC","TACG"))[::-1]
    comp.right_end_top, comp.left_end_bottom = brick.left_end_bottom, brick.right_end_top
    comp.right_end_bottom, comp.left_end_top = brick.left_end_top, brick.right_end_bottom
    return comp
    
    def circularize(brick, ovhg_check=True):
        brick = copy.deepcopy(brick)
        seq_origin = brick.seq
        feats_origin = brick.features
        if brick.topology == "circular" and brick.record.annotations["topology"] == "circular":
            print("The dbrick object is already circularized")

        if ovhg_check == True:
            if str(brick.right_end) == str(brick.left_end) and (brick.right_end_top * brick.left_end_bottom == 1 and brick.right_end_bottom * brick.left_end_top == 1) and len(brick.right_end) > 0:
                print("Based on complementary sticky end between 3' end and 5' end, the Dbrick object was circularized. The sticky end is '{}'".format(brick.left_end)) 
                ovhg = brick.right_end
                subbrick = substr(brick, 0, len(brick.seq)-len(brick.right_end))
                brick.seq      = subbrick.seq
                brick.record   = subbrick.record

            elif str(brick.right_end) == str(brick.left_end) and len(brick.right_end) > 0 and brick.right_end_top == 1 and brick.right_end_bottom == 1 and brick.left_end_top == 1 and brick.left_end_bottom == 1:
                print("Based on sequence homology between 3' end and 5' end, the Dbrick object was circularized. The overhang sequence is '{}'.".format(brick.left_end)) 
                ovhg = brick.right_end
                subbrick = substr(brick, 0, len(brick.seq)-len(brick.right_end))
                brick.seq      = subbrick.seq
                brick.record   = subbrick.record
            
            elif brick.right_end_top == 1 and brick.right_end_bottom == 1 and brick.left_end_top == 1 and brick.left_end_bottom == 1:
                #Gibson assembly
                flag = 0
                for i in range(Dbrick.min_overlap,Dbrick.max_overlap):
                    if brick.seq[:i] == brick.seq[-1*i:]:
                        flag = 1
                        break
                    else:
                        pass
                if flag == 1:
                    ovhg = brick.seq[:i] 
                    subbrick  = substr(brick, 0, len(brick.seq)-i)
                    brick.seq      = subbrick.seq
                    brick.record   = subbrick.record
                    print("Based on sequence homology between 3' end and 5' end, the Dbrick object was circularized. The overhang sequence is '{}'".format(ovhg))
                else:
                    return False
            
            else:
                return False
        else:
            ovhg = ""

        remove_list = [] 
        feats1      = [feat for feat in brick.features if "note_dbrick" in feat.qualifiers]
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
                        if key == "note_dbrick":
                            pass 
                        elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                            flag = 1
                        else:
                            flag = 0
                            break    
                    if flag == 1:
                        note1   = feat1.qualifiers["note_dbrick"]
                        pos_s1  = int(note1.split(":")[0].split("..")[0])
                        pos_e1  = int(note1.split(":")[0].split("..")[1]) 
                        length1 = int(note1.split(":")[1])
                        
                        note2   = feat2.qualifiers["note_dbrick"] 
                        pos_s2  = int(note2.split(":")[0].split("..")[0])
                        pos_e2  = int(note2.split(":")[0].split("..")[1]) 
                        length2 = int(note2.split(":")[1])
                        
                        if length1 == length2 and ((s1 == 0 and e2 == len(seq_origin)) or (s2 == 0 and e1 == len(seq_origin))):
                            if "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original:
                                note     = "{}..{}:{}".format(pos_s1, pos_e2, length1)
                                new_seq  = seq_origin[s1:e1] + seq_origin[s2:e2] 
                                new_feat = copy.deepcopy(brick.features[brick.features.index(feat1)]) 
                                if s1 == 0:
                                    new_feat.location = FeatureLocation(s2, e2 + (e1-len(ovhg)), feat1.strand)
                                else:
                                    new_feat.location = FeatureLocation(s1, e1 + (e2-len(ovhg)), feat1.strand)
                                
                                if len(new_seq) - len(ovhg) <= len(feat1.original):
                                    new_feat.qualifiers["note_dbrick"] = note
                                    if len(new_seq) - len(ovhg) == length1:
                                        del brick.features[brick.features.index(feat1)].qualifiers["note_dbrick"]
                                    brick.features[brick.features.index(feat1)].location = new_feat.location
                                    brick.features.remove(feat2)  
                                    remove_list.append(feat2) 
            
        for i in range(len(brick.features)):    
            if brick.features[i].location.parts[-1].end.position > len(brick.seq): #and feats_origin[i].type == "source":
                if brick.features[i].location.parts[0].start.position >= len(brick.seq):
                    strand                    = brick.features[i].location.strand
                    brick.features[i].location = FeatureLocation(brick.features[i].location.parts[0].start.position-len(brick.seq),brick.features[i].location.parts[-1].end.position-len(brick.seq))
                    brick.features[i].location.strand = strand
                else:
                    strand    = brick.features[i].location.strand
                    locations = [FeatureLocation(brick.features[i].location.parts[0].start.position,len(brick.seq)), FeatureLocation(0,brick.features[i].location.parts[-1].end.position-len(brick.seq))]
                    if strand == -1:
                        locations.reverse() 
                    brick.features[i].location = CompoundLocation(locations)
                    brick.features[i].location.strand = strand
        
        brick.record.features = brick.features
        brick.left_end  = ""
        brick.left_end_top    = 0 
        brick.left_end_bottom = 0 

        brick.right_end = ""
        brick.right_end_top    = 0 
        brick.right_end_bottom = 0 

        brick.topology = "circular"
        brick.record.seq = Seq(str(brick.seq),Alphabet.DNAAlphabet())
        brick.record.annotations["topology"] = brick.topology
        
def digestion(brick, *enzymes, target_seq=None, target_annotation=None):
    """
    Digest the DNA sequences by specified restriction enzymes.
    """
    fragments = [] 
    brick = copy.deepcopy(brick)
    seq   = Seq(str(brick.seq),amb) 
    enzymes       = list(enzymes)
    position_list = [(-10,0,None)] 
    for re in enzymes:
        re = re_dict[re] 
        if brick.topology != "linear":
            positions = re.search(seq, linear=False)   
        else:
            positions = re.search(seq)
        
        if re.ovhg < 0:
            positions = [(pos-1, pos-1-re.ovhg, re) for pos in positions]
        else:
            positions = [(pos-1-re.ovhg, pos-1, re) for pos in positions]
        position_list.extend(positions)

    if brick.topology == "linear":
        position_list.append((len(brick.seq),len(brick.seq), None))

    elif brick.topology == "circular":
        position_list[0] = position_list[-1]

    for p, pos in enumerate(position_list[:-1]):
        start = pos[1]
        end   = position_list[p+1][0]
        re5   = pos[-1] 
        re3   = position_list[p+1][-1]
       
        if re5 != None:
            if start < 0:
                left_end = str(seq[start-abs(re5.ovhg):] + seq[:start])
            else:
                left_end = str(seq[start-abs(re5.ovhg):start])
            if re5.ovhg < 0:
                left_end_top     = left_end
                left_end_bottom  = "-" * len(left_end) 
            else:
                left_end_top     = "-" * len(left_end)
                left_end_bottom  = left_end.translate(str.maketrans("ATGC","TACG"))
        else:
            left_end        = ""
            left_end_top    = ""
            left_end_bottom = "" 

        if re3 != None:
            if end > len(brick.seq):
                right_end = str(seq[end:len(brick.seq)] + seq[:end+abs(re3.ovhg)])
            else:
                right_end = str(seq[end:end+abs(re3.ovhg)])
            if re3.ovhg < 0:
                right_end_top    = "-" * len(right_end) 
                right_end_bottom = right_end.translate(str.maketrans("ATGC","TACG"))
            else:
                right_end_top    = right_end
                right_end_bottom = "-" * len(right_end)
        else:
            right_end        = ""
            right_end_top    = "" 
            right_end_bottom = "" 
        fragment = substr(brick, start-len(left_end), end+len(right_end))
        #print("{}/{}".format(left_end_top,left_end_bottom), "{}/{}".format(right_end_top, right_end_bottom).format()) 
        fragment = shell(fragment, "{}/{}".format(left_end_top,left_end_bottom), "{}/{}".format(right_end_top, right_end_bottom), add=0)
        #print(fragment.left_end, fragment.right_end) 
        fragments.append(fragment)

    if target_annotation is None and target_seq is None:
        if len(fragments) > 1 :    
            print("Return multiple fragments") 
            return fragments
        else:
            return fragments[0]
    else:
        return get_target(fragments, taget_seq=target_seq, target_annotation=target_annotation)

def pcr(brick, fw=None, rv=None, adaptor_fw=None, adaptor_rv=None, name=None):
    """
    Amplify specific region in the DNA brick.
    """
    brick = copy.deepcopy(brick) 
    if type(fw) == str or type(fw) == Seq:
        fw = str(fw) 
        fw_brick = Dbrick(seq=fw)
    
    elif type(fw) == Dbrick:
        fw_brick = fw
        fw = str(fw.seq)
    
    if type(rv) == str or type(rv) == Seq:
        rv = str(rv) 
        rv_brick = Dbrick(seq=rv) 
    
    elif type(rv) == Dbrick:
        rv_brick = rv
        rv = str(rv.seq) 
        
    rv_rc = rv.translate(str.maketrans("ATGC","TACG"))[::-1]  
    if brick.topology == "circular":
        seq = brick.seq + brick.seq 
    else:
        seq = brick.seq

    #Search the forward primer-binding region. 
    if fw in seq:
        adaptor_fw = ""
    
    elif adaptor_fw is None:
        for i in range(len(fw)):
            if fw[i:] in seq:
                adaptor_fw = fw[0:i] 
                fw = fw[i:]
                break
    else:
        pass 
   
    #Search the reverse primer-binding region.
    if rv_rc in seq:
        adaptor_rv = ""
    
    elif adaptor_rv is None:
        for i in range(len(rv)):
            if rv_rc[:len(rv)-i] in seq: 
                adaptor_rv = rv[:i]
                rv = rv[i:] 
                break
    else:
        pass
    rv_rc = rv.translate(str.maketrans("ATGC","TACG"))[::-1] 

    if name is None:
        name = "_" + brick.name
    else:
        name = "_" + name
    
    seq = brick.seq + brick.seq[:len(fw)-1]
    if seq.count(fw) > 1:
        raise ValueError("First sequence was not unique.")
    elif seq.count(fw) == 0:
        raise ValueError("First seqeucne was not found.")

    start = seq.find(fw) 
    
    seq = brick.seq + brick.seq[:len(rv)-1]
    if seq.count(rv_rc) > 1:
        raise ValueError("Second seqeunce was not unique")
    elif seq.count(rv_rc) == 0: 
        raise ValueError("Second sequence was not found")
    
    end = seq.find(rv_rc) + len(rv) 
    if end > len(brick.seq):
        end = end - len(brick.seq) 

    if len(adaptor_fw) > 0:
        adaptor_fw_brick = substr(fw_brick, start = 0,  end = len(adaptor_fw)) 
    
    if len(adaptor_rv) > 0:     
        adaptor_rv_brick = substr(rv_brick, start = 0,  end = len(adaptor_rv))
        adaptor_rv_brick = reverse_complement(adaptor_rv_brick) 
    
    amplicon   = substr(brick, start, end) 
    if len(adaptor_fw) > 0 and len(adaptor_rv) > 0:
        product = join_dbricks(adaptor_fw_brick, amplicon, adaptor_rv_brick, ovhg_check=False)
        print("The region from start {} to end {} was extracted. Adapter sequneces were detected at both ends. Right redundant sequence is {}, Left redundant sequence is {}.".format(start+1, end, adaptor_fw, adaptor_rv))

    elif len(adaptor_fw) > 0:
        product = join_dbricks(adaptor_fw_brick, amplicon, ovhg_check=False)
        print("The region from start {} to end {} was extracted. Adapter sequnece was detected at right end. The sequence is {}.".format(start+1, end, adaptor_fw))

    elif len(adaptor_rv) > 0:
        product = join_dbricks(amplicon, adaptor_rv_brick, ovhg_check=False)
        print("The region from start {} to end {} was extracted. Adapter sequnece was detected at left end. The nsequence is {}.".format(start+1, end, adaptor_rv))

    else:
        product = amplicon       
        print("The region from start {} to end {} was extracted.".format(start+1, end))
    
    if len(adaptor_fw) > 0 or len(adaptor_rv) > 0:
        product.start = 0 
        product.end   = len(product.seq)
    else:
        product.start = start
        product.end   = end
    product.name = amplicon.name
    product.record.id = amplicon.name
    return product

class Dbrick():
    order_qualifiers = ["label","gene","product"]
    min_overlap      = 10 
    max_overlap      = 500

    def __init__(self, record=None, seq=None, name="None", topology="linear", format=None):
        if seq is None and record is None:
            self.seq            = None
            self.record         = None
            self.features       = None
            self.right_end      = None 
            self.left_end       = None
            self.topology          = topology
            self.source            = None
            self.name              = name
        
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
                self.right_end = self.seq[0:10] 
                self.left_end  = self.seq[-10:]
            else:
                self.right_end = ""
                self.left_end  = ""
            
            self.source = None
            self.name = record.id  

        elif record is None:
            self.seq      = str(seq).upper() 
            self.record   = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
            self.features = self.record.features
            self.topology = topology
            if self.topology == "linear":
                self.right_end = self.seq[0:10] 
                self.left_end  = self.seq[-10:]
            else:
                self.right_end = ""
                self.left_end  = ""
            self.source = None
            self.name = name
        
        self.left_end_top      = 1
        self.left_end_bottom   = 1
        self.right_end_top     = 1 
        self.right_end_bottom  = 1
        self.start = 0 
        self.end   = len(self.seq) 
    
    def __getitem__(self, item):
        if type(item) == slice:
            if item.step is None:
                strand = 1
            elif type(item.step) != int:
                raise TypeError("Invalid type for index.")
            else:
                strand = item.step

            if item.start is None:
                start = 0 
            else:
                if type(item.start) == int:
                    start = item.start

                elif type(item.start) == str or type(item.start) == Seq:
                    start = str(item.start)
                        
                elif type(item.start) == Dbrick:
                    start = item.start 

                else:
                    raise TypeError("Invalid type was detected in index. Please specify proper index.")
           
            if item.stop is None:
                end = len(self.seq) 
            else:
                if type(item.stop) == int:
                    end = item.stop
                
                elif type(item.stop) == str or type(item.stop) == Seq:
                    end = str(item.stop)
                        
                elif type(item.stop) == Dbrick:
                    end = item.stop
                
                else:
                    raise TypeError("Invalid type was detected in index. Please specify proper index.") 
            
            if type(start) == int and type(end) == int:
                if start < 0:
                    start = len(self.seq) - abs(start) 
                if end < 0: 
                    end = len(self.seq) - abs(end) 
                subbrick = substr(self,start,end)
                            
            elif type(start) == str and type(end) == str and start.isdecimal() and end.isdecimal():
                subbrick = substr(self, -1*abs(strand), abs(strand), feature_id=[start,end])
            
            elif (type(start) == str or type(start) == Seq or type(start) == Dbrick) or (type(end) == str or type(end) == Seq or type(end) == Dbrick):
                subbrick = substr(self, start, end)

            else:
                raise TypeError("Invalid type was detected in index. Please specify proper index.")
            
            subbrick.source = self
            if strand == -1 or strand < 0:
                return reverse_complement(subbrick)
            else:
                return subbrick 
        
        elif type(item) == int:
            subbrick = substr(self,item,item+1) 
            return subbrick
        
        elif type(item) == str and item.isdecimal() == True: 
            subbrick = substr(self, 0, 0, info=item)
            return subbrick
        
        elif type(item) == Seq or (type(item) == str and set(item) <= set("ATGCNatgcn")) or type(item) == Dbrick: 
            if type(item) == Dbrick: 
                item = str(item.seq) 
            else:
                item = str(item)
            subbrick = substr(self, 0, 0, info=item)
            return subbrick
        
        else:
            raise ValueError("Invalid index type was specified.") 
    
    def __setitem__(self, key, value):
        if (type(value)== str and set(value) <= set("ATGCNatgcn")) or type(value) == Seq or type(value) == Dbrick:
            if (type(value)== str and set(value) <= set("ATGCNatgcn")) or type(value) == Seq:
                value = Dbrick(seq=value) 
            else:
                pass
        else:
            raise TypeError("Invalid type for 'value'. 'value' should be specified as 'str', 'Seq', or 'Dbrick' object.") 

        if type(key) == int:
            if len(value.seq) == 1:
                origin = str(self.seq)[key]
                self.seq = self.seq[:key] + value.seq + self.seq[key+1:]
                self.record.seq = Seq(str(self.seq), Alphabet.DNAAlphabet())
                rec = SeqFeature(FeatureLocation(key,key,strand=0), type="misc_feature")
                rec.qualifiers = {"note":["{}>{}".foramt(origin, str(value.seq))]}
                self.features.append(rec) 
                self.features.sort(key=lambda x:x.location.parts[0].start.position)
                self.record.features = self.features

            else:
                raise TypeError("The length of the sequence must be one in case of assignment to a single position.")

        elif type(key) == str and key.isdecimal():
            start = self.features[int(key)].location.parts[0].start  
            end   = self.features[int(key)].location.parts[-1].end 
        
        elif type(key) == slice:
            if key.step is None:
                strand = 1
            elif type(key.step) != int:
                raise TypeError("Invalid type for index.")
            else:
                strand = key.step

            if key.start is None:
                start = 0 
            else:
                if type(key.start) == int:
                    start = key.start

                elif type(key.start) == str or type(key.start) == Seq:
                    start = str(key.start)
                        
                elif type(key.start) == Dbrick:
                    start = str(key.start.seq) 

                else:
                    raise TypeError("Invalid type was detected in index.")
           
            if key.stop is None:
                end = len(key.seq) 
            else:
                if type(key.stop) == int:
                    end = key.stop
                
                elif type(key.stop) == str or type(key.stop) == Seq:
                    end = str(key.stop)
                        
                elif type(key.stop) == Dbrick:
                    end = str(key.stop.seq)

                else:
                    raise TypeError("Invalid type was detected in index.") 
            
            if type(start) == int and type(end) == int:
                pass                 
            
            elif type(start) == str and type(end) == str and set(start) <= set("ATGCNatgcn") and set(end) <= set("ATGCNatgcn"):
                _end  = end
                start = str(self.seq).find(start) 
                end   = str(self.seq).find(str(reverse_complement(_end).seq)) + len(_end) 
                if start < 0:
                    raise ValueError("First index was not found in the sequence of the 'Dbrick' object.") 
                
                if end < len(_end):
                    raise ValueError("Second index was not found in the sequence of the 'Dbrick' object.") 
                
            if start > end:
                if self.topology == "linear":
                    raise ValueError("When you assign 'Dbrick' obejct into the other linear 'dbrick' object, start value should be larger than end value.")
                else:
                    origin = str((self[start:] + self[:end]).seq)
                    brick  = self[end:start] + value 
                    brick.circularize(ovhg_check=False)
                    rec = SeqFeature(FeatureLocation(len(brick.seq)-len(value.seq), len(brick.seq), strand=0), type="misc_feature")
                    rec.qualifiers = {"note":["{}>{}".format(origin, str(value.seq))]}

            else:
                origin = str(self[start:end].seq) 
                brick  = self[0:start] + value + self[end:]
                rec = SeqFeature(FeatureLocation(start, start+len(value.seq), strand=0), type="misc_feature")
                rec.qualifiers = {"note":["{}>{}".format(origin, str(value.seq))]}
            
            brick.topology = self.topology
            brick.record.annotations["topology"] = self.topology
            for key in brick.__dict__:
                self.__dict__[key] = brick.__dict__[key]
            self.features.append(rec) 
            self.features.sort(key=lambda x:x.location.parts[0].start.position)
            self.record.features = self.features

    def __add__(self, other):
        if (type(other) == str and set(other) <= set("ATGCNatgcn")) or type(other) == Seq:
            other = Dbrick(seq=other) 

        elif type(other) == SeqRecord:
            other = Dbrick(record=other) 

        elif type(other) == Dbrick:
            pass 
        
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized Dbrick object cannot be joined with others.") 
        else:
            return join_dbricks(self, other, ovhg_check=False) 

    def __radd__(self, other):
        if (type(other) == str and set(other) <= set("ATGCNatgcn")) or type(other) == Seq:
            other = Dbrick(seq=other) 

        elif type(other) == Seq.SeqRecord:
            other = Dbrick(record=other) 

        elif type(other) == Dbrick:
            pass 
           
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized Dbrick object cannot be joined with others.") 
        else:
            return join_dbricks(other, self, ovhg_check=False) 
    
    def __xor__(self, other):
        if (type(other) == str and set(other) <= set("ATGCNatgcn")) or type(other) == Seq:
            other = Dbrick(seq=other) 

        elif type(other) == SeqRecord:
            other = Dbrick(record=other) 

        elif type(other) == Dbrick:
            pass 
        
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicularized Dbrick object cannot be joined with others.") 
        else:
            return join_dbricks(self, other, ovhg_check=True, min_overlap=Dbrick.min_overlap, max_overlap=Dbrick.max_overlap) 
    
    def __rxor__(self, other):
        if (type(other) == str and set(other) <= set("ATGCNatgcn")) or type(other) == Seq:
            other = Dbrick(seq=other) 

        elif type(other) == Seq.SeqRecord:
            other = Dbrick(record=other) 

        elif type(other) == Dbrick:
            pass 
           
        if other.topology == "circular" or self.topology == "circular":
            raise ValueError("Cicular Dbrick object cannot be joined with others.") 
        else:
            return join_dbricks(other, self, ovhg_check=True, min_overlap=Dbrick.min_overlap, max_overlap=Dbrick.max_overlap) 

    def view_seq(self, whole=False, flength=None):
        if flength is None:
            if len(self.seq) - len(self.right_end) - len(self.left_end) > 100:
                flength = 5 
            elif len(self.seq) - len(self.right_end) - len(self.left_end) > 10:
                flength = 2 
            else:
                flength = 0 
            
        seq_rc = self.seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        if self.left_end_top == 1:
            left_end_top = self.seq[:len(self.left_end)+flength]
        else:
            left_end_top = " " * len(self.left_end) + self.seq[len(self.left_end):len(self.left_end)+flength] 
        
        if self.left_end_bottom == 1:
            left_end_bottom = self.seq[:len(self.left_end)+flength].translate(str.maketrans("ATGC","TACG")) 
        else:
            left_end_bottom = " " * len(self.left_end) + self.seq[len(self.left_end):len(self.left_end)+flength].translate(str.maketrans("ATGC","TACG"))
        
        if self.right_end_top == 1:
            right_end_top = self.seq[len(self.seq)-len(self.right_end)-flength:]
        else:
            right_end_top = self.seq[len(self.seq)-len(self.right_end)-flength:len(self.seq)-len(self.right_end)] + " " * len(self.right_end) 
        
        if self.right_end_bottom == 1:
            right_end_bottom = self.seq[len(self.seq)-len(self.right_end)-flength:].translate(str.maketrans("ATGC","TACG"))
        else:
            right_end_bottom = self.seq[len(self.seq)-len(self.right_end)-flength:len(self.seq)-len(self.right_end)].translate(str.maketrans("ATGC","TACG")) + " " * len(right_end)        
        
        if whole == False:
            print("{}...{}".format(left_end_top, right_end_top))
            print("{}...{}".format(left_end_bottom, right_end_bottom))
        else:
            print("{}{}{}".format(left_end_top, self.seq[len(self.left_end)+flength:len(self.seq)-len(self.right_end)-flength], right_end_top))
            print("{}{}{}".format(left_end_bottom, self.seq[len(self.left_end)+flength:len(self.seq)-len(self.right_end)-flength].translate(str.maketrans("ATGC","TACG")), right_end_bottom))

    def add_feature(self, start, end=None, strand=1, feature_type="misc_feature", qualifiers={}):
        if type(qualifiers) != dict:
            raise TypeError("Invalid type for 'qualifiers'. 'qualifiers' should be specified as 'dict' object.") 

        if (type(start) == str and set(start) <= set("ATGCNatgcn")) or type(start) == Seq or type(start) == Dbrick:
            if type(start) == Seq:
                seq = str(start).upper() 
            elif type(start) == Dbrick:
                seq = str(start.seq).upper()
            else:
                seq = start.upper() 

            if strand == -1:
                seq = seq.translate(str.maketrans("ATGC","TACG"))[::-1]

            i = 0 
            subject = self.seq + self.seq[:len(seq)-1]
            starts, ends = [], [] 
            while seq in subject[i:]:
                tmp = subject[i:].find(seq) 
                starts.append(i + tmp)
                if i + tmp + len(seq) > len(self.seq):
                    ends.append(i + tmp + len(seq) - len(self.seq)) 
                else:
                    ends.append(i + tmp + len(seq))
                i = i + tmp + 1
            
            for start, end in zip(starts, ends):
                self.add_feature(start=start, end=end, strand=strand, feature_type=feature_type, qualifiers=qualifiers)
                print("New feature was added in the range of start {} to end {}.".format(start,end))     
            
            if len(starts) == 0:
                print("The sequence was not found.") 
            return None
        
        elif type(start) == int and type(end) == int:
            if start > end:
                if self.topology == "linear":
                    raise ValueError("When you add feature to linear dbrick object, start value should be larger than end value.")
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
            for s,e in zip(start,end):
                if s > e:
                    if self.topology == "linear":
                        raise ValueError("When you add feature to linear dbrick object, start value should be larger than end value.")
                    else:
                        locations.append(FeatureLocation(start, len(self.seq), strand=strand))
                        locations.append(FeatureLocation(0, end, strand=strand))
                else:
                    locations.append(FeatureLocation(start, end, strand=strand))
            if strand == -1:
                locations.reverse() 
            feat = SeqFeature(CompoundLocation(locations, type=feature_type))
            feat.location.strand = strand 

        for key, value in qualifiers.items():
            if type(qualifiers[key]) == str or type(qualifiers[key]) == int:
                qualifiers[key] = [value] 
            else:
                pass

        feat.qualifiers = qualifiers
        self.features.append(feat) 
        self.features.sort(key=lambda x:x.location.parts[0].start.position)
        self.record.features = self.features

    def remove_feature(self, feature_ids):
        tmp_features = copy.copy(self.features)
        for _id in feature_ids:
            tmp_features.remove(self.features[_id]) 
        self.features = tmp_features 
        self.record.features = tmp_features
   
    def get_feature_index(self, start=None, end=None,  qualifiers={}, feature_type="all"): 
        if type(qualifiers) != dict:
            raise TypeError("Invalid type for 'qualifiers'. 'qualifiers' should be specified as dict.") 
        if feature_type == "all":
            all = True
        else:
            all = False
        if start is None:
            start = 0 

        if end is None:
            end = len(self.seq) 
        
        feats = [] 
        for feat in self.features:
            strand = feat.location.strand
            if strand == -1:
                s = feat.location.parts[-1].start.position 
                e = feat.location.parts[0].end.position
            else:
                s = feat.location.parts[0].start.position 
                e = feat.location.parts[-1].end.position
            if e > start and s < end:
                feats.append(feat) 
        
        feature_ids = [] 
        for _id, feat in feats:
            if all == True or feat.type == feature_type:
                if qualifiers == {}:
                    feature_ids.append(str(_id))
                else:
                    for key in qualifiers:
                        if key in feat.qualifiers: 
                            if type(feat.qualifiers[key]) == list and qualifiers[key] in feat.qualifiers[key]:
                                features.append(str(_id)) 
                            elif type(feat.qualifiers[key]) == str and qualifiers[key] == feat.qualifiers[key]:
                                features.append(str(_id)) 
                            else:
                                pass
                        else:
                            pass
        if len(features) == 0:
            print("There was no features that corresponded to the search")
        return feature_ids
   
    def view_features(self, feature_ids=None, detail=False, sep=None, o=None):
        _ids    = ["Feature_ID"] 
        labels  = ["Label"]
        types   = ["Type"] 
        starts  = ["Start"] 
        ends    = ["End"] 
        strands = ["Strand"]
        if feature_ids is None:
            feature_ids = list(range(len(self.features)))
        else:
            feature_ids = list(map(int, feature_ids))
        for _id in feature_ids:
            feat = self.features[_id]
            flag = 0 
            
            for key in Dbrick.order_qualifiers:
                if key in feat.qualifiers:
                    flag = 1
                    break 
                else:
                    pass 
            if flag == 0:
                label = feat.type
            else:
                if type(feat.qualifiers[key]) == list:
                    label = feat.qualifiers[key][0]
                else:
                    label = feat.qualifiers[key]
            strand = feat.location.strand
            if strand == -1:
                start = feat.location.parts[-1].start.position
                end   = feat.location.parts[0].end.position 
            else:
                start = feat.location.parts[0].start.position
                end   = feat.location.parts[-1].end.position 
            
            _ids.append(str(_id)) 
            labels.append(str(label)) 
            types.append(str(feat.type)) 
            starts.append(str(start)) 
            ends.append(str(end))
            if strand == 1:
                strands.append("+")
            elif strand == 0:
                strands.append("+") 
            else:
                strands.append("-")
            
        _idmax    = max(list(map(len,_ids)))   + 2
        labelmax  = max(list(map(len,labels))) + 2
        ftypemax  = max(list(map(len,types)))  + 2
        startmax  = max(list(map(len,starts))) + 2
        endmax    = max(list(map(len,ends)))   + 2
        strandmax = max(list(map(len,strands)))   + 2
        for n, (_id, label, ftype, start, end, strand) in enumerate(zip(_ids, labels, types, starts, ends, strands)):
            if sep is None:
                if o == io.TextIOWrapper:
                    print(_id + " " * (_idmax-len(_id)) + label + " " * (labelmax-len(label)) + ftype + " " * (ftypemax-len(ftype)) + start + " " * (startmax-len(start)) + end + " " * (endmax-len(end)) + strand + " " * (strandmax-len(strand)), file=o)
                else:
                    print(_id + " " * (_idmax-len(_id)) + label + " " * (labelmax-len(label)) + ftype + " " * (ftypemax-len(ftype)) + start + " " * (startmax-len(start)) + end + " " * (endmax-len(end)) + strand + " " * (strandmax-len(strand)))
            else:
                if o == io.TextIOWrapper:
                    print(_id, label, ftype, start, end, strand, sep=",", file=o)
                else:
                    print(_id, label, ftype, start, end, strand, sep=",") 
            
            if detail == True and n > 0:
                feat = self.features[int(_id)]
                for key,value in feat.qualifiers.items():
                    if type(value) == list:
                        for v in value:
                            print(" " * (_idmax) + "  ", "{}:{}".format(key,v)) 
                    else:
                        print(" " * (_idmax) + "  ",  "{}:{}".format(key,value))
                print() 

    def reindex(self, zero_position):
        if self.topology == "linear":
            print("The method cannot be used for linear Dbrick object.")
        else:
            abrick = self[zero_position:] + self[:zero_position] 
            abrick.circularize(ovhg_check=False) 
            for key in abrick.__dict__.keys():
                self.__dict__[key] = abrick.__dict__[key] 
            
    def write(self, handle, format="genbank"):
        features = copy.deepcopy(self.features)
        for feat in features:
            if "note_dbrick" in feat.qualifiers:
                note = feat.qualifiers["note_dbrick"]
                pos_s  = int(note.split(":")[0].split("..")[0])
                pos_e  = int(note.split(":")[0].split("..")[1])
                length = int(note.split(":")[1])
                if (pos_s == 1 and pos_e == length) or (pos_s == length and pos_e == 1):
                    del feat.qualifiers["note_dbrick"]
        self.record.features = features 
        self.record.id = self.name
        self.record.seq = Seq(str(self.seq),Alphabet.DNAAlphabet()) 
        SeqIO.write(self.record, handle, format)
        self.record.features = self.features

    def linearize(self):
        abrick = substr(self, 0, len(self.seq))
        for key in abrick.__dict__.keys():
            self.__dict__[key] = abrick.__dict__[key] 
        self.topology = "linear"
        self.record.annotations["topology"] = self.topology

    def circularize(self, ovhg_check=True):
        seq_origin = copy.deepcopy(self.seq)
        feats_origin = copy.deepcopy(self.features) 
        if self.topology == "circular" and self.record.annotations["topology"] == "circular":
            print("The dbrick object is already circularized")

        if ovhg_check == True:
            if str(self.right_end) == str(self.left_end) and (self.right_end_top * self.left_end_bottom == 1 and self.right_end_bottom * self.left_end_top == 1) and len(self.right_end) > 0:
                print("Based on complementary sticky end between 3' end and 5' end, the Dbrick object was circularized. The sticky end is '{}'".format(self.left_end)) 
                ovhg = self.right_end
                brick = substr(self, 0, len(self.seq)-len(self.right_end))
                self.seq      = brick.seq
                self.record   = brick.record

            elif str(self.right_end) == str(self.left_end) and len(self.right_end) > 0 and self.right_end_top == 1 and self.right_end_bottom == 1 and self.left_end_top == 1 and self.left_end_bottom == 1:
                print("Based on sequence homology between 3' end and 5' end, the Dbrick object was circularized. The overhang sequence is '{}'.".format(self.left_end)) 
                ovhg = self.right_end
                brick = substr(self, 0, len(self.seq)-len(self.right_end))
                self.seq      = brick.seq
                self.record   = brick.record
            
            elif self.right_end_top == 1 and self.right_end_bottom == 1 and self.left_end_top == 1 and self.left_end_bottom == 1:
                #Gibson assembly
                flag = 0
                for i in range(Dbrick.min_overlap,Dbrick.max_overlap):
                    if self.seq[:i] == self.seq[-1*i:]:
                        flag = 1
                        break
                    else:
                        pass
                if flag == 1:
                    ovhg = self.seq[:i] 
                    brick  = substr(self, 0, len(self.seq)-i)
                    self.seq      = brick.seq
                    self.record   = brick.record
                    print("Based on sequence homology between 3' end and 5' end, the Dbrick object was circularized. The overhang sequence is '{}'".format(ovhg))
                else:
                    return False
            
            else:
                return False
        else:
            ovhg = ""

        remove_list = [] 
        feats1      = [feat for feat in self.features if "note_dbrick" in feat.qualifiers]
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
                        if key == "note_dbrick":
                            pass 
                        elif key in feat2.qualifiers and feat1.qualifiers[key] == feat2.qualifiers[key]:
                            flag = 1
                        else:
                            flag = 0
                            break    
                    if flag == 1:
                        note1   = feat1.qualifiers["note_dbrick"]
                        pos_s1  = int(note1.split(":")[0].split("..")[0])
                        pos_e1  = int(note1.split(":")[0].split("..")[1]) 
                        length1 = int(note1.split(":")[1])
                        
                        note2   = feat2.qualifiers["note_dbrick"] 
                        pos_s2  = int(note2.split(":")[0].split("..")[0])
                        pos_e2  = int(note2.split(":")[0].split("..")[1]) 
                        length2 = int(note2.split(":")[1])
                        
                        if length1 == length2 and ((s1 == 0 and e2 == len(seq_origin)) or (s2 == 0 and e1 == len(seq_origin))):
                            if "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original:
                                note     = "{}..{}:{}".format(pos_s1, pos_e2, length1)
                                new_seq  = seq_origin[s1:e1] + seq_origin[s2:e2] 
                                new_feat = copy.deepcopy(self.features[self.features.index(feat1)]) 
                                if s1 == 0:
                                    new_feat.location = FeatureLocation(s2, e2 + (e1-len(ovhg)), feat1.strand)
                                else:
                                    new_feat.location = FeatureLocation(s1, e1 + (e2-len(ovhg)), feat1.strand)
                                
                                if len(new_seq) - len(ovhg) <= len(feat1.original):
                                    new_feat.qualifiers["note_dbrick"] = note
                                    if len(new_seq) - len(ovhg) == length1:
                                        del self.features[self.features.index(feat1)].qualifiers["note_dbrick"]
                                    self.features[self.features.index(feat1)].location = new_feat.location
                                    self.features.remove(feat2)  
                                    remove_list.append(feat2) 
            
        for i in range(len(self.features)):    
            if self.features[i].location.parts[-1].end.position > len(self.seq): #and feats_origin[i].type == "source":
                if self.features[i].location.parts[0].start.position >= len(self.seq):
                    strand                    = self.features[i].location.strand
                    self.features[i].location = FeatureLocation(self.features[i].location.parts[0].start.position-len(self.seq),self.features[i].location.parts[-1].end.position-len(self.seq))
                    self.features[i].location.strand = strand
                else:
                    strand    = self.features[i].location.strand
                    locations = [FeatureLocation(self.features[i].location.parts[0].start.position,len(self.seq)), FeatureLocation(0,self.features[i].location.parts[-1].end.position-len(self.seq))]
                    if strand == -1:
                        locations.reverse() 
                    self.features[i].location = CompoundLocation(locations)
                    self.features[i].location.strand = strand
        
        self.record.features = self.features
        self.left_end  = ""
        self.left_end_top    = 0 
        self.left_end_bottom = 0 

        self.right_end = ""
        self.right_end_top    = 0 
        self.right_end_bottom = 0 

        self.topology = "circular"
        self.record.seq = Seq(str(self.seq),Alphabet.DNAAlphabet())
        self.record.annotations["topology"] = self.topology
