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

def substr(brick, start=0, end=0, feature_id=None):
    """
    Extract record information within the specified region.a
    """
    if feature_id != None:
        if type(feature_id) == list: 
            feature_s = brick.features[int(feature_id[0])]
            feature_e = brick.features[int(feature_id[-1])]
            start     = feature_s.location.parts[0].start.position + start
            end       = feature_e.location.parts[-1].end.position + end

        elif type(feature_id) == str and feature_id.isdecimal():
            feature = brick.features[int(feature_id)] 
            start   = feature.location.parts[0].start.position + start
            end     = feature.location.parts[-1].end.position + end

        else:
            raise TypeError("Invalid value for 'feture_id'.'feature_id' should be specified as 'str' or 'list' object.")
            
    if start > end:
        if brick.topology == "circular":
            sub_brick1 = substr(brick, start, len(brick.seq))
            sub_brick2 = substr(brick, 0, end)
            sub_brick  = join_dbricks(sub_brick1, sub_brick2, ovhg_check=False)
        else:
            raise ValueError("When you extract specific region from linear dbrick object, start value should be larger than end value.")
    else:
        feats   = []
        new_features = [] 
        for feat in brick.record.features:
            s = feat.location.parts[0].start.position
            e = feat.location.parts[-1].end.position
            strand = feat.location.strand
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

                length = len(brick.seq) - s + e
                if "note_dbrick" not in feat1.qualifiers:
                    feat1.qualifiers["note_dbrick"] = "{}..{}:{}".format(1, len(brick.seq)-s, length) 
                else:
                    note   = feat.qualifiers["note_dbrick"]
                    pos_s  = int(note.split(":")[0].split("..")[0]) 
                    pos_e  = int(note.split(":")[0].split("..")[1]) 
                    length = int(note.split(":")[1])
                    note   = "{}..{}:{}".format(pos_s, pos_s + len(brick.seq)-s, length)
                    feat1.qualifiers["note_dbrick"] = note
                
                if "note_dbrick" not in feat2.qualifiers:
                    feat2.qualifiers["note_dbrick"] = "{}..{}:{}".format(len(brick.seq)-s+1, len(brick.seq)-s+e, length) 
                else:
                    note   = feat.qualifiers["note_dbrick"]
                    pos_s  = int(note.split(":")[0].split("..")[0]) 
                    pos_e  = int(note.split(":")[0].split("..")[1]) 
                    length = int(note.split(":")[1])
                    note   = "{}..{}:{}".format(pos_s + len(brick.seq)-s, pos_e, length)
                    feat2.qualifiers["note_dbrick"] = note
                new_features.append(feat1)
                new_features.append(feat2)
            else:
                new_features.append(feat) 

        for feat in new_features:
            s = feat.location.parts[0].start.position 
            e = feat.location.parts[-1].end.position            
            if len(feat.location.parts) == 1 and s <= e:
                strand = feat.location.strand 
                if e > start and s < end:
                    feat = copy.deepcopy(feat) 
                    if s - start <= 0:
                        feat.location.parts[0]._start = ExactPosition(0)
                        if "note_dbrick" not in feat.qualifiers:
                            feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(abs(s-start)+1, e-s, e-s) 
                        else:
                            note   = feat.qualifiers["note_dbrick"]
                            pos_s  = int(note.split(":")[0].split("..")[0]) + abs(s-start) 
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
                            feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(1, end-s, e-s) 
                        else:
                            s      = int(feat.location.parts[0].start.position)
                            note   = feat.qualifiers["note_dbrick"]
                            pos_s  = int(note.split(":")[0].split("..")[0])
                            pos_e  = int(note.split(":")[0].split("..")[1]) 
                            length = int(note.split(":")[1])
                            note   = "{}..{}:{}".format(pos_s, pos_s + (end-start-s) - 1, length)
                            feat.qualifiers["note_dbrick"] = note
                    
                    feat.location.strand = strand
                    feats.append(feat)
            else:
                strand = feat.location.strand
                length = e-s
                locations = []
                sflag = 0 
                eflag = 0
                for lbrick in feat.location.parts:
                    s = lbrick.start.position 
                    e = lbrick.end.position
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
                            feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(abs(s-start)+1, e-s, length) 
                        else:
                            note   = feat.qualifiers["note_dbrick"]
                            pos_s  = int(note.split(":")[0].split("..")[0]) + abs(s-start) 
                            pos_e  = int(note.split(":")[0].split("..")[1]) 
                            length = int(note.split(":")[1])
                            note   = "{}..{}:{}".format(pos_s, pos_e, length)
                            feat.qualifiers["note_dbrick"] = note
                    else:
                        locations[0][0] = ExactPosition(s - start)
                    
                    if e > end-start and eflag == 1:
                        locations[-1][1] = ExactPosition(end-start)
                        if "note_dbrick" not in feat.qualifiers:
                            feat.qualifiers["note_dbrick"] = "{}..{}:{}".format(1, end-s, length) 
                        else:
                            s      = int(locations[0][0])
                            note   = feat.qualifiers["note_dbrick"]
                            pos_s  = int(note.split(":")[0].split("..")[0])
                            pos_e  = int(note.split(":")[0].split("..")[1]) 
                            length = int(note.split(":")[1])
                            note   = "{}..{}:{}".format(pos_s, pos_s + (end-start-s)-1, length)
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
            sub_brick.left_end            = brick.left_end
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
            sub_brick.right_end          = sub_brick.seq[-20:0]
            sub_brick.right_end_top      = 1
            sub_brick.right_end_bottom   = 1

        new_rec          = copy.deepcopy(brick.record)
        new_rec.seq      = Seq(brick.seq[start:end],Alphabet.DNAAlphabet()) 
        new_rec.features = sub_brick.features
        new_rec.annotations["topology"] = sub_brick.topology
        sub_brick.record  = new_rec
        sub_brick.source  = brick
    return sub_brick 

def join_dbricks(*bricks, topology="linear", name=None, ovhg_check=True, min_overlap=15, max_overlap=150):
    """    
    Join the multiple DNA bricks. Each DNA brick should have a common overhang sequence with the adjacent element.  
    """
    def slide(feats,slide):
        new_feats = []
        for feat in feats:
            feat = copy.deepcopy(feat) 
            s = feat.location.parts[0].start.position 
            e = feat.location.parts[-1].end.position
            strand = feat.location.strand
            feat.location.parts[0]._start = ExactPosition(s + slide) 
            feat.location.parts[-1]._end  = ExactPosition(e + slide) 
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
                        print("Joined dbrick objects based on complementary sticky end of each DNA fragmment. Sticky end is '{}'".format(sticky_end)) 
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
                        print("Joined dbrick objects based on sequence homology at the end of each fragment. Overhang sequence is '{}'".format(ovhg)) 
                    
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
                            print("Joined dbrick objects based on sequence homology at the end of each fragment. Overhang sequence is '{}'".format(ovhg)) 

                        else:
                            print("These DNA fragments cannot be joined. Please reconsider the design of the construction")
                            return False
                else:
                    print(brick.left_end, brick.left_end_top, brick.left_end_bottom) 
                    print(brick.right_end, brick.right_end_top, brick.right_end_bottom)
                    print("These DNA fragments cannot be joined. Please reconsider the design of the construction")
                    return False
            else:
                new_brick = brick
                ovhg = ""
            feats  = slide(feats, len(construct.seq) - len(ovhg))
            feats1 = [feat for feat in construct.features if "note_dbrick" in feat.qualifiers]
            feats2 = [feat for feat in feats if "note_dbrick" in feat.qualifiers]
            if len(feats1) > 0 and len(feats2) > 0:
                for feat1 in feats1:
                    s1, e1 = feat1.location.parts[0].start.position, feat1.location.parts[-1].end.position
                    for feat2 in feats2:
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
                                
                                if length1 == length2:
                                    #print(feat1,feat2,"original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original)
                                    if "original" in feat1.__dict__ and "original" in feat2.__dict__ and feat1.original == feat2.original and s2 == 0:
                                        note     = "{}..{}:{}".format(pos_s1, pos_e2, length1)
                                        new_seq  = construct.seq[s1:e1] + brick.seq[s2:e2] 
                                        new_feat = copy.deepcopy(construct.features[construct.features.index(feat1)]) 
                                        strand   = new_feat.location.strand
                                        
                                        if len(feat1.location.parts) == 1 and len(feat2.location.parts) == 1:
                                            new_feat.location = FeatureLocation(feat1.location.parts[0].start.position, feat2.location.parts[-1].end.position, feat1.strand)
                                            new_feat.location.strand = strand
                                        else:
                                            locations = feat1.location.parts[0:-1] + [FeatureLocation(feat1.location.parts[-1].start.position, feat2.location.parts[0].end.position, feat1.strand)] + feat2.location.parts[0:-1]
                                            if strand == -1:
                                                locations.reverse() 
                                            new_feat.location = CompoundLocation(locations) 
                                            new_feat.locaiton.strand = strand 
                                        #print(len(new_seq) - len(ovhg), new_feat.location.parts[-1].end.position, new_feat.location.parts[0].start.position, feat1.original)  
                                        if len(new_seq) - len(ovhg) == new_feat.location.parts[-1].end.position - new_feat.location.parts[0].start.position and len(new_seq) - len(ovhg) <= len(feat1.original):
                                            new_feat.qualifiers["note_dbrick"] = note
                                            if pos_s1 == 1 and pos_s2 == length1:
                                                del construct.features[construct.features.index(feat1)].qualifiers["note_dbrick"]
                                            construct.features[construct.features.index(feat1)].location = new_feat.location
                                            del feats[feats.index(feat2)] 
                                            #print("delete", feat2)  
            construct.features = construct.features + feats
            construct.seq               = construct.seq + new_brick.seq 
            construct.right_end         = brick.right_end
            construct.right_end_top     = brick.right_end_top
            construct.right_end_bottom  = brick.right_end_bottom
            construct.source            = brick.source
            construct.topology      = "linear"
            ovhg = brick.right_end
    
    if len(ovhgs) != len(set(ovhgs)): 
        print("The same overhangs sequences were found from multiple joints. The process will generate unexpected constructs.") 

    construct.features.sort(key=lambda x:x.location.parts[0].start.position)
    new_record          = SeqRecord(Seq(str(construct.seq),Alphabet.DNAAlphabet()))
    new_record.features = construct.features
    new_record.annotations["topology"] = topology
    construct.record = new_record 

    
    if name == None:
        pass
    else:
        new_record.id = name
        construct.name = name 
    
    if topology == "circular":
        construct.circularize(ovhg_check=ovhg_check)
    return construct

def shell(brick, left="", right="", add=0): #top5 = 1, bottom5 = 1, top3 = 1, bottom3 = 1, add=0):
    """
    Set overhang sequence.
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
    if pattern1.fullmatch(left_end) != None and pattern2.search(left_end) == None and patternl1.search(left_end) == None and patternl2.search(left_end) == None:
        pass 
    else:
        raise ValueError("Please sepcify a proper sequence pattern for 'left'") 
    
    if pattern1.fullmatch(right_end) != None and pattern2.search(right_end) == None and patternr1.search(right_end) == None and patternr2.search(right_end) == None:
        pass
    else:
        raise ValueError("Please sepcify a proper sequence pattern for 'right'") 

    if "/" in left_end:
        left_end_top, left_end_bottom = left_end.split("/")
        if len(left_end_top) != len(left_end_bottom):
            raise ValueError("Top strand length and bottom strand length should be same size")
        
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
                raise ValueError("Top strand and bottom strand should be complent each other.")
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
            raise ValueError("Top strand length and bottom strand length should be same size")
        
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
                raise ValueError("Top strand and bottom strand should be complent each other.")
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
    #print(brick.seq[left_end_length-left_length:left_end_length-left_length+len(left_end)], str(brick[len(brick.seq)-right_end_length + right_length - len(right_end):len(brick.seq)-right_end_length + right_length].seq))
    if add == 1 or (left_end != brick.seq[left_end_length-left_length:left_end_length-left_length+len(left_end)] or right_end != str(brick[len(brick.seq)-right_end_length + right_length - len(right_end):len(brick.seq)-right_end_length + right_length].seq)):
        #print(1, left_end_length, right_end_length, len(brick.seq))
        left_end  = Dbrick(seq=left_end) 
        right_end = Dbrick(seq=right_end) 
        brick     = join_dbricks(left_end, brick, right_end, ovhg_check=False)
        brick.left_end  = left_end.seq
        brick.right_end = right_end.seq
    else:
        #print(2, left_end_length, right_end_length)
        left_end  = Dbrick(seq=left_end) 
        right_end = Dbrick(seq=right_end)
        #brick = join_dbricks(left_end, brick[left_end_length-left_length+len(left_end.seq):len(brick.seq)-right_end_length + right_length - len(right_end.seq)], right_end, ovhg_check=False)
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
    Extract a DNA brick specified by the specific sequence or annotation from the multiple DNA elements. 
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

def complement(brick):
    if type(brick) == Dbrick:
        brick = copy.deepcopy(brick)
        seq  = brick.seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        feats = [] 
        for feat in brick.record.features:
            s = feat.location.parts[0].start.position
            e = feat.location.parts[-1].end.position
            strand = feat.location.strand
            feat.location.parts[0]._start = ExactPosition(len(brick.seq) - e) 
            feat.location.parts[-1]._end  = ExactPosition(len(brick.seq) - s) 
            if strand == 1 or strand == -1:
                feat.location.strand = -1 * feat.location.strand
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


def digestion(brick, *enzymes, target_seq=None, target_annotation=None):
    """
    Digest the DNA brick by specified restriction enzymes.
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

    if target_annotation == None and target_seq == None:
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
    if type(fw) == str:
        fw_brick = Dbrick(seq=fw) 
    else:
        fw_brick = fw
        fw = str(fw.seq) 
    
    if type(rv) == str:
        rv_brick = Dbrick(seq=rv) 
    else:
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
    
    elif adaptor_fw == None:
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
    
    elif adaptor_rv == None:
        for i in range(len(rv)):
            if rv_rc[:len(rv)-i] in seq: 
                adaptor_rv = rv[:i]
                rv = rv[i:] 
                break
    else:
        pass
    rv_rc = rv.translate(str.maketrans("ATGC","TACG"))[::-1] 

    if name == None:
        name = "_" + brick.name
    else:
        name = "_" + name
    
    seq = brick.seq + brick.seq[:len(fw)-1]
    if seq.count(fw) > 1:
        raise ValueError("Any unique sequence was not detected in the first index.")
    elif seq.count(fw) == 0:
        raise ValueError("Any seqeunce in the first index was not found in subject.")
    start = seq.find(fw) 
    
    seq = brick.seq + brick.seq[:len(rv)-1]
    if seq.count(rv_rc) > 1:
        raise ValueError("Any unique sequence was not detected in the second index.")
    elif seq.count(rv_rc) == 0: 
        raise ValueError("Any seqeunce in the second index was not found in subject.")
    
    end = seq.find(rv_rc) + len(rv) 
    if end > len(brick.seq):
        end = end - len(brick.seq) 

    if len(adaptor_fw) > 0:
        adaptor_fw_brick = substr(fw_brick, start = 0,  end = len(adaptor_fw)) 
    
    if len(adaptor_rv) > 0:     
        adaptor_rv_brick = substr(rv_brick, start = 0,  end = len(adaptor_rv))
        adaptor_rv_brick = complement(adaptor_rv_brick) 
    
    amplicon   = substr(brick, start, end) 
    if len(adaptor_fw) > 0 and len(adaptor_rv) > 0:
        product = join_dbricks(adaptor_fw_brick, amplicon, adaptor_rv_brick, ovhg_check=False)
        print("Start:{}, End:{}, Redundant sequneces were detected as both ends. Right redundant sequence is {}, Left redundant sequence is {}.".format(start+1, end, adaptor_fw, adaptor_rv))

    elif len(adaptor_fw) > 0:
        product = join_dbricks(adaptor_fw_brick, amplicon, ovhg_check=False)
        print("Start:{}, End:{}, Redundant sequnece was detected at right end. The sequence is {}.".format(start+1, end, adaptor_fw))

    elif len(adaptor_rv) > 0:
        product = join_dbricks(amplicon, adaptor_rv_brick, ovhg_check=False)
        print("Start:{}, End:{}, Redundant sequnece was detected at left end. The nsequence is {}.".format(start+1, end, adaptor_rv))

    else:
        product = amplicon       
        print("Start:{}, End:{}".format(start+1, end))

    return product

class Dbrick():
    order_qualifiers = ["label","gene","product"]
    min_overlap      = 10 
    max_overlap      = 500

    def __init__(self, record=None, seq=None, name="None", topology="linear", format=None):
        if seq == None and record == None:
            self.seq            = None
            self.record         = None
            self.features       = None
            self.right_end      = None 
            self.left_end       = None
            self.left_end_top      = 0 
            self.left_end_bottom   = 0
            self.right_end_top     = 0 
            self.right_end_bottom  = 0
            self.topology          = topology
            self.source            = None
            self.name              = name
        
        elif seq == None:
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
            self.seq      = str(record.seq).upper()
            self.record   = record
            if fmt == "genbank":
                self.features = record.features 
            else:
                self.features = []
            
            if "topology" in record.annotations:
                self.topology = record.annotations["topology"]
            else:
                self.topology = topology
            
            if self.topology == "linear":
                self.right_end         = self.seq[0:10] 
                self.left_end          = self.seq[-10:]
            else:
                self.right_end          = ""
                self.left_end          = ""
            
            self.left_end_top      = 1
            self.left_end_bottom   = 1
            self.right_end_top     = 1
            self.right_end_bottom  = 1
            self.source = None
            self.name = record.id  

        elif record == None:
            self.seq      = str(seq).upper() 
            self.record   = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
            self.features = self.record.features
            self.topology = topology
            if self.topology == "linear":
                self.right_end          = self.seq[0:10] 
                self.left_end          = self.seq[-10:]
            else:
                self.right_end       = ""
                self.left_end        = ""
            self.left_end_top      = 1
            self.left_end_bottom   = 1
            self.right_end_top     = 1
            self.right_end_bottom  = 1
            self.source = None
            self.name = name
        
    def __getitem__(self, item):
        if type(item) == slice:
            if item.step == None:
                strand = 1
            elif type(item.step) != int:
                raise TypeError("Invalid type for index.")
            else:
                strand = item.step

            if item.start == None:
                start = 0 
            else:
                if type(item.start) == int:
                    start = item.start

                elif type(item.start) == str or type(item.start) == Seq:
                    start = str(item.start)
                        
                elif type(item.start) == Dbrick:
                    start = str(item.start.seq) 

                else:
                    raise TypeError("Invalid type was detected in index. Please specify proper index.")
           
            if item.stop == None:
                end = len(self.seq) 
            else:
                if type(item.stop) == int:
                    end = item.stop
                
                elif type(item.stop) == str or type(item.stop) == Seq:
                    end = str(item.stop)
                        
                elif type(item.stop) == Dbrick:
                    end = item.stop.seq
                
                else:
                    raise TypeError("Invalid type was detected in index. Please specify proper index.") 
            
            if type(start) == int and type(end) == int:
                if start < 0:
                    start = len(self.seq) - abs(start) 
                if end < 0: 
                    end = len(self.seq) - abs(end) 
                subbrick = substr(self,start,end)
                
            elif type(start) == str and type(end) == str and set(start) <= set("ATGCNatgcn") and set(end) <= set("ATGCNatgcn"):
                subbrick = pcr(self,start.upper(),end.upper())
            
            elif type(start) == str and type(end) == str and start.isdecimal() and end.isdecimal():
                subbrick = substr(self, -1*abs(strand), abs(strand), feature_id=[start,end])

            else:
                raise TypeError("Invalid type was detected in index. Please specify proper index.")
            
            subbrick.source = self
            if strand == -1 or strand < 0:
                return complement(subbrick)
            else:
               return subbrick 
        
        elif type(item) == int:
            subbrick = substr(self,item,item+1) 

        elif type(item) == str:
            subbrick = substr(self, -1*abs(strand), abs(strand), feature_id=item)

        else:
            raise ValueError("Dbrick class only accept slice index.") 
    
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
            if key.step == None:
                strand = 1
            elif type(key.step) != int:
                raise TypeError("Invalid type for index.")
            else:
                strand = key.step

            if key.start == None:
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
           
            if key.stop == None:
                end = len(key.seq) 
            else:
                if type(key.stop) == int:
                    end = key.stop
                
                elif type(key.stop) == str or type(key.stop) == Seq:
                    end = str(key.stop)
                        
                elif type(key.stop) == Dbrick:
                    end = key.stop.seq

                else:
                    raise TypeError("Invalid type was detected in index.") 
            
            if type(start) == int and type(end) == int:
                pass                 
            
            elif type(start) == str and type(end) == str and set(start) <= set("ATGCNatgcn") and set(end) <= set("ATGCNatgcn"):
                _end  = end
                start = str(self.seq).find(start) 
                end   = str(self.seq).find(str(complement(_end).seq)) + len(_end) 
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

    def view_seq(self):
        pass 

    def add_feature(self, start=0, end=0, seq=None, strand=1, feature_type="misc_feature", qualifiers={}):
        if seq == None:
            pass 
        
        elif (type(other) == str and set(other) <= set("ATGCNatgcn")):
            i = 0 
            starts, ends = [], [] 
            while seq in self.seq[i:]:
                tmp = self.seq[i:].find(seq) 
                i = i + tmp + 1
                starts.append(i + tmp)
                ends.append(i + tmp + len(seq))
            
            for start, end in zip(starts, ends):
                add_feature(self, start=start, end=end, seq=None, strand=strand, feature_type=feature_type, qualifiers={})
                print("New feature was added. Start:{}, End:{}:".format(start,end)) 
        else:
            pass

        if type(qualifiers) != dict:
            raise TypeError("Invalid type for 'qualifiers'. 'qualifiers' should be specified as 'dict' object.") 

        if type(start) == int and type(end) == int:
            if start > end:
                if self.topology == "linear":
                    raise ValueError("When you add feature to linear dbrick object, start value should be larger than end value.")
                else:
                    locations = [FeatureLocation(start, len(self.length)), FeatureLocation(0, end)]
                    if strand == -1:
                        locations.reverse() 
                    feat = SeqFeature(CompoundLocation(locations), type=feature_type)
                    feat.location.strand = strand
            else:
                feat = SeqFeature(FeatureLocation(start, end, strand=strand), type=feature_type)
        
        if type(start) == list and type(end) == list:
            locations = [] 
            for s,e in zip(start,end):
                if s > e:
                    if self.topology == "linear":
                        raise ValueError("When you add feature to linear dbrick object, start value should be larger than end value.")
                    else:
                        locations.append(FeatureLocation(start, len(self.length), strand=strand))
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
        if start == None:
            start = 0 

        if end == None:
            end = len(self.seq) 
        
        feats = [] 
        for feat in self.features:
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
        _ids   = ["Feature_ID"] 
        labels = ["Label"]
        types  = ["Type"] 
        starts = ["Start"] 
        ends   = ["End"] 

        if feature_ids == None:
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
            
            start = feat.location.parts[0].start.position
            end   = feat.location.parts[-1].end.position 
            
            _ids.append(str(_id)) 
            labels.append(str(label)) 
            types.append(str(feat.type)) 
            starts.append(str(start)) 
            ends.append(str(end)) 
            
        _idmax   = max(list(map(len,_ids)))   + 2
        labelmax = max(list(map(len,labels))) + 2
        ftypemax = max(list(map(len,types)))  + 2
        startmax = max(list(map(len,starts))) + 2
        endmax   = max(list(map(len,ends)))   + 2
        for n, (_id, label, ftype, start, end) in enumerate(zip(_ids, labels, types, starts, ends)):
            if sep == None:
                if o == io.TextIOWrapper:
                    print(_id + " " * (_idmax-len(_id)) + label + " " * (labelmax-len(label)) + ftype + " " * (ftypemax-len(ftype)) + start + " " * (startmax-len(start)) + end + " " * (endmax-len(end)), file=o)
                else:
                    print(_id + " " * (_idmax-len(_id)) + label + " " * (labelmax-len(label)) + ftype + " " * (ftypemax-len(ftype)) + start + " " * (startmax-len(start)) + end + " " * (endmax-len(end)))
            else:
                if o == io.TextIOWrapper:
                    print(_id, label, ftype, start, end, sep=",", file=o)
                else:
                    print(_id, label, ftype, start, end, sep=",") 
            
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
            print("The method cannot work for linear Dbrick object.")
        else:
            feats = [] 
            for feat in self.features:
                s = feat.location.parts[0].start.position 
                e = feat.location.parts[-1].end.position
                strand = feat.location.strand
                
                feat = copy.deepcopy(feat) 
                new_start = ExactPosition(s - zero_position) 
                new_end   = ExactPosition(e - zero_position) 
                if new_start < 0:
                    if new_end < 0:
                        feat.location = FeatureLocation(len(self.seq) + new_start, len(self.seq) + new_end)
                        feat.location.strand = strand 
                    else:
                        locations = [FeatureLocation(len(self.seq) + new_start, len(self.seq)), FeatureLocation(0, new_end)] 
                        if strand == -1: 
                            locations.reverse() 
                        feat.location = CompoundLocation(locations)
                        feat.location.strand = strand 
                else:
                    feat.location = FeatureLocation(new_start, new_end)
                feats.append(feat)
            
            self.record.features = feats
            self.features        = self.record.features
            self.record.seq      = Seq(str(self.seq)[cor:] + str(self.seq)[:cor],Alphabet.DNAAlphabet())
            self.seq             = str(self.record.seq)  

    def write(self, handle, format="genbank", reindex=0):
        self.record.id = self.name
        if reindex == 0:
            SeqIO.write(self.record, handle, format)
        else:
            self.reindex(reindex)
            SeqIO.write(self.record, handle, format)
    
    def circularize(self, ovhg_check=True):
        seq_origin = copy.deepcopy(self.seq)
        feats_origin = copy.deepcopy(self.features) 
        if self.topology == "circular" and self.record.annotations["topology"] == "circular":
            print("The dbrick object is already circularized")

        if ovhg_check == True:
            if str(self.right_end) == str(self.left_end) and (self.right_end_top * self.left_end_bottom == 1 and self.right_end_bottom * self.left_end_top == 1) and len(self.right_end) > 0:
                print("Circularization based on complementary sticky end between 3' end and 5' end. Sticky end is '{}'".format(self.left_end)) 
                ovhg = self.right_end
                brick = substr(self, 0, len(self.seq)-len(self.right_end))
                self.seq      = brick.seq
                self.record   = brick.record

            elif str(self.right_end) == str(self.left_end) and len(self.right_end) > 0 and self.right_end_top == 1 and self.right_end_bottom == 1 and self.left_end_top == 1 and self.left_end_bottom == 1:
                print("Circularization based on sequence homology between 3' end and 5' end. Overhang sequence is '{}'.".format(self.left_end)) 
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
                    print("Circularization based on sequence homology between 3' end and 5' end. Overhang sequence is '{}'".format(ovhg))
                else:
                    return False
            
            else:
                return False
        else:
            ovhg = ""

        if len(ovhg) > 0:
            remove_list = [] 
            feats1      = [feat for feat in self.features if "note_dbrick" in feat.qualifiers]
            for feat1 in feats1:
                s1, e1 = feat1.location.parts[0].start.position, feat1.location.parts[-1].end.position
                for feat2 in feats1:
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
                                        #print("delete", feat2) 
            
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
                        #print(feats_origin[i].location) 
            self.record.features = self.features
        
        self.left_end  = ""
        self.left_end_top    = 0 
        self.left_end_bottom = 0 

        self.right_end = ""
        self.right_end_top    = 0 
        self.right_end_bottom = 0 

        self.topology = "circular"
        self.record.seq      = Seq(str(self.seq),Alphabet.DNAAlphabet())
        self.record.annotations["topology"] = self.topology
