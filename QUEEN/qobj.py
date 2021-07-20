import sys
import copy 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, FeatureLocation, ExactPosition
from functools import total_ordering

sys.path.append("/".join(__file__.split("/")[:-1]))
from qfunction import *

class Qint(int):
    def __init__(self, num):
        self.qkey        = None
        self.parent      = None
        self.parental_id = None
        self.name        = None
        self._qint       = True

class Qseq(str):
    def __init__(self, seq):
        self.qkey           = None
        self.parent         = None
        self.parental_id    = None
        self.parental_class = None
        self.name           = None 
        self.item           = None
        
    def __getitem__(self, item):
        value = super().__getitem__(item) 
        if self.parental_class == "QUEEN" and self.parent is not None and self.parent.topology == "circular":
            if type(item) == slice:
                if item.start is None:
                    start = 0 
                elif item.start < 0:
                    start = len(self) + item.start
                else:
                    start = item.start
                
                if item.stop is None:
                    stop = len(self)
                elif item.stop < 0:
                    stop = len(self) + item.stop
                else:
                    stop = item.stop
                
                if item.step is None:
                    step = 1

                if start > stop: 
                    value = str(self)[start:] + str(self)[:stop]
                    value = value[::step] 
                else:
                    value = super().__getitem__(item)
            else:
                pass 
            value = Qseq(value)
        else:
            value = Qseq(value)
        value.qkey           = self.qkey
        value.parent         = None
        value.parental_id    = self.parental_id
        value.parental_class = self.parental_class 
        value.name           = self.name
        value.item           = item
        
        return value

@total_ordering
class DNAfeature(SeqFeature):
    def __deepcopy__(self, memo):
        obj = DNAfeature(self, subject=self.subject)
        return obj 

    def __getattr__(self, name):
        if name == "feature_id":
            if "_id" in self.__dict__:
                return self._id
            else:
                return None 

        elif name == "feature_type":
            return self.type
        
        elif name == "_original":
            seq = self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        
        elif name == "original":
            if "_original" in self.__dict__:
                return self._original 
            else:
                return self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        
        elif name == "seq" or name == "sequence":
            seq           = self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
            seq           = Qseq(seq)
            seq.qkey      = self._start.qkey
            seq.parental_id = self._start.parental_id 
            seq.parent    = self 
            seq.parental_class = "DNAfeature"
            seq = self.subject.getdnaseq(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
            return seq
        
        elif name == "strand":
            return self.location.strand
        
        elif name == "start":
            return self._start
        
        elif name == "end":
            return self._end 
        
        elif name == "span":
            return (self.start, self.end)
        
        else:
            raise AttributeError("DNAfeature obejct has no attribute '{}'".format(name))
    
    def __eq__(self, other):
        if not isinstance(other, DNAfeature):
            return NotImplemented
        else:
            if self.qualifiers == other.qualifiers and (self.subject is not None and other.subject is not None and self.seq == other.seq):
                return None
            else: 
                return False

    def __lt__(self, other):
        return NotImplemented

    def __setattr__(self, key, value):
        if key in ["feature_id", "feature_type", "seq", "sequenece", "original", "strand", "start", "end", "span"]:
            raise ValueError("Cannot assign to '{}' attribute. To set or change feature attribute value, plase use 'editfeature' function.".format(key))
        else:
            super.__setattr__(self, key, value)

    def __init__(self, feature=None, location=None, type="misc_feature", subject=None, query=None):
        if feature is None:
            SeqFeature.__init__(self, location, type)
        else:
            for key in feature.__dict__:
                if key in ["_start", "_end", "__digestion_topl", "_digestion_topr", "_digestion_bottomr", "_digestion_bottoml", "subject", "query"]: 
                    pass
                else:
                    self.__dict__[key] = copy.deepcopy(feature.__dict__[key]) 
        
        #Start->End direction should be 5' to 3' on the top strand.
        if self.location.strand == -1:
            self._start = self.location.parts[-1].start.position
            self._end   = self.location.parts[0].end.position
        else:
            self._start = self.location.parts[0].start.position
            self._end   = self.location.parts[-1].end.position
        
        self._start = Qint(self._start)
        self._end   = Qint(self._end)
        self._start.parent = self
        self._end.parent   = self
        self._start.name   = "start"
        self._end.name     = "end"

        if subject is None:
            pass 
        else:
            self._digestion_topl = "null" 
            self._digestion_topr = "null"
            self._digestion_bottoml = "null"
            self._digestion_bottomr = "null" 
        
        self.subject   = subject
        self.query_seq = query
        self._dnafeature = 1

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

def _compile_cutsite(query):  
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


class QUEEN():    
    #Class variables that manage the execution histories of operational and search function
    dna_dict                   = {}
    queried_feature_dict       = {}
    queried_features_dict      = {}
    queried_features_name_dict = {} 
    process_description = None
    _namespace     = {} 
    _namespaceflag = 0 
    _num_history   = 1  
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
                if self._history_feature is None:
                    obj._history_feature = None
                else:
                    obj._history_feature = DNAfeature(self._history_feature, subject=obj)     
            
            elif key == "_dnafeatures": 
                feats = [] 
                for feat in self.dnafeatures:
                    feats.append(DNAfeature(feat, subject=obj))
                obj._dnafeatures = feats
            else:
                obj.__dict__[key] = copy.copy(self.__dict__[key])
        return obj

    def __repr__(self):
        if len(self.seq) > 50:
            out = "<queen.QUEEN object; project='{}', length='{} bp', topology='{}'>".format(self.project, len(self.seq), self.topology)
        else:
            out = "<queen.QUEEN object; project='{}', length='{} bp', sequence='{}', topology='{}'>".format(self.project, len(self.seq), self.seq, self.topology)
        return out 
    
    def __setattr__(self, key, value):
        if key == "_project":
            if "_unique_id" in self.__dict__:
                if value in self.__class__.dna_dict:
                    if value.split("_")[-1].isdecimal() == True:
                        value = "_".join(self.project.split("_")[:-1]) 
                    unique = 0
                    while value + "_" + str(unique) in self.__class__.dna_dict:
                        unique += 1    
                    _unique_id = value + "_" + str(unique)
                else:         
                    _unique_id = value
                QUEEN.dna_dict[_unique_id] = None
                self._unique_id = _unique_id
                super.__setattr__(self, key, value)
            else:
                QUEEN.dna_dict[value] = None
                super.__setattr__(self, key, value)
           
        elif key == "project":
            if "_unique_id" in self.__dict__:
                if value in self.__class__.dna_dict:
                    if value.split("_")[-1].isdecimal() == True:
                        value = "_".join(dna.project.split("_")[:-1]) 
                    unique = 0
                    while value + "_" + str(unique) in self.__class__.dna_dict:
                        unique += 1    
                    _unique_id = value + "_" + str(unique)
                else:         
                    _unique_id = value
                
                if "_history_feature" in self.__dict__ and self._history_feature is not None:
                    for qualifier in self._history_feature.qualifiers:
                        if "building_history" in qualifier[0:18]: 
                            history = self._history_feature.qualifiers[qualifier][0]
                            history = re.sub("QUEEN.dna_dict\['{}'\]".format(self._unique_id), "QUEEN.dna_dict['{}']".format(_unique_id), history)
                            self._history_feature.qualifiers[qualifier][0] = history
                            #re.sub("(QUEEN.queried_features_dict\['){}(_.*'\])".format(self._unique_id), "\1{}\2".format(_unique_id), history)
         
                QUEEN.dna_dict[_unique_id] = None
                self._unique_id = _unique_id
                super.__setattr__(self, "_project", value)
            else:
                QUEEN.dna_dict[value] = None
                super.__setattr__(self, "_project", value)

        elif key == "seq" or key == "dnafeatures":
            raise ValueError("Cannot assign to '{}' attribute.".format(key))   
        else:
            super.__setattr__(self, key, value) 
    
    def __getattribute__(self, name):
        if name == "_seq":
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
                    history = self._history_feature.qualifiers[key]
                    histories.append((int(key.split("_")[-1]), history[0], history[1])) 
            histories.sort()
            return histories
        elif name == "seq":
            return self._seq
        elif name == "rcseq":
            rcseq = self._seq.upper().translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
            rcseq = Qseq(rcseq)
            if "_unique_id" in self.__dict__:
                rcseq.parental_id = self._unique_id 
            rcseq.parent = self 
            rcseq.parental_class = "QUEEN"
            rcseq.name = "rcseq"
            return rcseq 
        elif name == "project":
            return self._project
        elif name == "sequence":
            return self._seq 
        elif name == "dnafeatures":
            return self._dnafeatures
        else:
            raise AttributeError("Queen obejct has no attribute '{}'".format(name))

    def __init__(self, seq=None, record=None, project=None, topology="linear", fileformat=None, product=None, process_description=None, import_history=True, _direct=1):
        if process_description is None:
            pass  
        else:
            QUEEN.process_description = process_description 
        
        fseq      = seq 
        frecord   = record
        fproject  = project 
        ftopology = topology
        fproduct  = product 

        self._seq               = None
        self.record             = None
        self._dnafeatures       = None
        self._right_end         = None 
        self._left_end          = None
        self._history_feature   = None
        if seq is None and record is None:
            if project is None: 
                project = "dna"
            self.topology           = topology
            if "_" in project:
                project = project.replace("_","-")
            self._project           = project
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1
        
        elif seq is None or "." in seq:
            if "." in str(seq):
                record = seq
            
            if type(record) == str:
                if fileformat != None:
                    fmt = fileformat
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
                record = SeqIO.parse(record,None)

            self._seq = str(record.seq).upper()
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
                    project = frecord.split("/")[-1].split(".")
                    project = project if len(project) == 1 else ".".join(project[:-1])
                else:
                    project = record.id 
            if "_" in project:
                project = project.replace("_","-")
            self._project = project
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1
            
            #import features
            self._dnafeatures = [] 
            if len(record.features) > 0:
                for feat in record.features:
                    self._dnafeatures.append(DNAfeature(feature=feat, subject=self))
               
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
                                pairs.append((feat, history_num, feat.qualifiers[key]))
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
                        feat.qualifiers["building_history_{}".format(new_history_num)] = pair[2]
                        del feat.qualifiers["building_history_{}".format(pair[1])]
                        history_nums.append(new_history_num)      
                else:
                    QUEEN.process_description = process_description
                    archivehistory(self)    
                
                QUEEN._num_history = max(history_nums) 
                
                if history_feature is not None:
                    self._history_feature = history_feature
                    self._dnafeatures.remove(history_feature) 
            
            if len(self.dnafeatures) == 0:
                import_history = False
                self._dnafeatures = []

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
                
                self._seq       = str(seq).upper()
                if "_" in project:
                    project = project.replace("_","-")
                self._project   = project
                if Alphabet:
                    self.record = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
                else:
                    self.record = SeqRecord(Seq(str(seq)))
                
                self._dnafeatures = []
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
        
        if project in QUEEN.dna_dict:
            keys   = list(QUEEN.dna_dict.keys())
            unique = 0
            while project + "_" + str(unique) in QUEEN.dna_dict:
                unique += 1
            self._unique_id = project + "_" + str(unique)
        else: 
            self._unique_id = project 
        QUEEN.dna_dict[self._unique_id] = product

        if _direct == 1 and import_history == False:                
            fseq       = "" if fseq is None else "seq='{}'".format(fseq)
            frecord    = "" if frecord is None else "record='{}'".format(frecord) if fseq == "" else ", record='{}'".format(frecord)
            fproject   = "" if fproject is None else ", project='{}'".format(fproject)
            ftopology  = "" if topology == "linear" else ", topology='{}'".format(topology)
            fproduct   = "" if fproduct is None else ", product='{}'".format(fproduct)
            fileformat = "" if fileformat is None else ", fileformat='{}'".format(fileformat) 
            process_description = "" if process_description is None else ", process_description='{}'".format(process_description)
            args = [fseq, frecord, fproject, ftopology, fileformat, fproduct, process_description]
            
            building_history = "QUEEN.dna_dict['{}'] = QUEEN({}{}{}{}{}{}{})".format(self._unique_id, *args) 
            
            QUEEN._num_history += 1 
            feat = SeqFeature(FeatureLocation(0, len(self.seq), strand=1), type="source") 
            feat._id = "0"
            feat.qualifiers["label"]       = [self.project] 
            feat.qualifiers["description"] = ["Record of building history"]
            feat.qualifiers["building_history_{}".format(QUEEN._num_history)] = [building_history.replace(" ","–"), ""] 
            feat = DNAfeature(feature=feat, subject=self)         
            self._history_feature = feat

        else:
            pass 

        self._positions       = tuple(range(len(self.seq))) 
        self.record.feartures = self.dnafeatures
        self.seq.parental_id  = self._unique_id
        if product is None:
            pass 
        else:
            QUEEN._namespace[product] = self

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
                match_list = get_matchlist_regex(self, query, value=None, subject=subject, s=start, e=end, strand=1) 
                match_list.extend(get_matchlist_regex(self, query, value=None, subject=subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], s=start, e=end, strand=-1))
            elif strand == 1:
                match_list = get_matchlist_regex(self, query, value=None, subject=subject, s=start, e=end, strand=strand) 
            elif strand == -1:
                match_list = get_matchlist_regex(self, query, value=None, subject=subject.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1], s=start, e=end, strand=strand)
            else:
                ValueError("When edit the sequence, the sequence strand to be edit should be '-1' or '+1.'")
            
            match_positions = set() 
            for match in match_list:
                span = (match["start"], match["end"])
                if span not in match_positions:
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
                    if qorigin.name != None: 
                        if "getdnaseq" in qorigin.name:
                            if len(qorigin.name.split("_")) == 2: 
                                seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, qorigin.name.split("_")[-1]) 
                            else:
                                seqname = "QUEEN.dna_dict['{}'].getdnaseq(start={}, end={}, strand={})".format(parental_id, *qorigin.name.split("_")[1:])
                        if qorigin.name == "rcseq":
                            seqname = "QUEEN.dna_dict['{}'].rcseq".format(parental_id)
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
            
            fproduct = "" if product is None else ", product='{}'".format(product)
            process_description = "" if process_description is None else ", process_description='{}'".format(process_description)
            if start == 0 and end == len(self.seq):
                if strand == 2:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}{}{})".format(qkey, self._unique_id, qorigin, fproduct, process_description)
                    add_history(self, [building_history, "query:{}".format(query)])
                else:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, strand={}{}{})".format(qkey, self._unique_id, qorigin, strand, fproduct, process_description)
                    add_history(self, [building_history, "query:{}; strand{}".format(query, strand)])
            else:
                building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, start={}, end={}, strand={}{}{})".format(qkey, self._unique_id, qorigin, start, end, strand, fproduct, process_description)   
                add_history(self, [building_history, "query:{}; start:{}; end:{}; strand:{}".format(query, start, end, strand)])
            QUEEN._qnum += 1 
        
        if product is None:
            pass 
        else:
            QUEEN._namespace[product] = feat_list
  
        return feat_list

    def searchfeature(self, key_attribute="all", query=".+", source=None, start=0, end=None, strand=2, product=None, process_description=None, _direct=1):
        if process_description is None:
            process_description = QUEEN.process_description
        else:
            QUEEN.process_description = process_description  
        
        start   = 0 if start == len(self.seq) else start
        end     = len(self.seq) if end is None else end
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
                    if query.name != None:
                        if "getdnaseq" in query.name:
                            if len(query.name.split("_")) == 2: 
                                seqname = "QUEEN.dna_dict['{}'].getdnaseq(strand={})".format(parental_id, query.name.split("_")[-1]) 
                            else:
                                seqname = "QUEEN.dna_dict['{}'].getdnaseq(start={}, end={}, strand={})".format(parental_id, *query.name.split("_")[1:])
                        if query.name == "rcseq":
                            seqname = "QUEEN.dna_dict['{}'].rcseq".format(parental_id)
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
            
            fproduct = "" if product is None else ", product='{}'".format(product)
            process_description = "" if process_description is None else ", process_description='{}'".format(process_description)
            if start == 0 and end == len(self.seq):
                if strand == 2 and source is None:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}{}{})".format(qkey, self._unique_id, key_attribute, query, fproduct, process_description)
                    add_history(self, [building_history, "key_attribute:{}; query:{}".format(key_attribute, query)])
                elif strand == 2:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}, source={}{}{})".format(qkey, self._unique_id, key_attribute, query, source, fproduct, process_description)
                    add_history(self, [building_history, "key_attribute:{}; query:{}; soruce:{}".format(key_attribute, query, source)])

                else:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}, source={}, strand={}{}{})".format(qkey, self._unique_id, key_attribute, query, source, strand, fproduct, process_description)
                    add_history(self, [building_history, "key_attribute:{}; query:{}; soruce:{}; strand:{}".format(key_attribute, query, source, strand)])

            else:
                building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}, source={}, start={}, end={}, strand={}{}{})".format(qkey, self._unique_id, key_attribute, query, source, start, end, strand, fproduct, process_description) 
                add_history(self, [building_history, "key_attribute:{}; query:{}; soruce:{}; start:{}; end:{}; strand:{}".format(key_attribute, query, source, start, end, strand)]) 
            QUEEN._qnum += 1 

        if product is None:
            pass 
        else:
            QUEEN._namespace[product] = features

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

        elif type(other) ==QUEEN:
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
        if start == 0 or None and end == len(self.seq) or None:
            return_seq.name = "getdnaseq_{}".format(strand)
        else:
            return_seq.name = "getdnaseq_{}_{}_{}".format(start, end, strand)
        return return_seq

    def _setfeatureid(self):
        for i in range(0, len(self.dnafeatures)):
            if i == 0:
                self._dnafeatures[i]._id = str(1)
            else:
                self._dnafeatures[i]._id = str(i*100)
    
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
