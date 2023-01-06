import sys
import copy
import urllib
import tempfile
import requests
from bs4 import BeautifulSoup 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, FeatureLocation, ExactPosition
from functools import total_ordering

sys.path.append("/".join(__file__.split("/")[:-1]))
from qfunction import *
from quine import * 
from qint import Qint
from qseq import Qseq

def _combine_history(dna, history_features):
    history_feature = SeqFeature(FeatureLocation(0, len(dna.seq), strand=0), type="source")
    history_feature = history_features[0].__class__(history_feature, subject=dna.seq)
    for feat in history_features:
        for key in feat.qualifiers:
            if "building_history" in key and feat.qualifiers[key] not in history_feature.qualifiers.values():
                history_feature.qualifiers[key] = feat.qualifiers[key]  
    return history_feature

@total_ordering
class DNAfeature(SeqFeature):
    """DNA feature object

    Each DNAfeature object with the following attributes provides an annotation for a 
    given range of DNA sequence in a QUEEN object.

    """
    
    def __deepcopy__(self, memo):
        obj = DNAfeature(self, subject=self.subject)
        return obj 
    
    def __getattribute__(self, name):
        if name == "strand":
            return self.location.strand if self.location.strand is not None else 0
        else:
            return super().__getattribute__(name)  

    def __getattr__(self, name):
        if name == "feature_id":
            if "_id" in self.__dict__:
                return self._id
            else:
                return None 
        elif name == "feature_type":
            return self.type
        
        elif name == "_original":
            seq = self.subject.printsequence(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
        
        elif name == "original":
            if "_original" in self.__dict__:
                return self._original 
            else:
                return self.subject.printsequence(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
         
        elif name == "seq" or name == "sequence":
            seq             = self.subject.printsequence(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
            seq             = Qseq(seq)
            seq.qkey        = self._start.qkey
            seq.parental_id = self._start.parental_id 
            seq.parent      = self 
            seq.parental_class = "DNAfeature"
            seq = self.subject.printsequence(self.start, self.end, self.location.strand if self.location.strand !=0 else 1) 
            self._seq = seq
            return seq
        
        elif name == "strand":
            return self.location.strand if self.location.strand is not None else 0
        
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
            if self.type == other.type and self.qualifiers == other.qualifiers and (self.subject is not None and other.subject is not None and self.seq == other.seq) and self.start == other.start and self.end == other.end:
                return True
            else: 
                return False

    def __lt__(self, other):
        return NotImplemented

    def __setattr__(self, key, value):
        if key in ["feature_id", "feature_type", "seq", "sequenece", "original", "strand", "start", "end", "span"]:
            raise AttributeError("'DNAfeature' object attribute '{}' is read-only".format(key))
        else:
            super.__setattr__(self, key, value)

    def __init__(self, feature=None, location=None, type="misc_feature", subject=None, query=None):
        if feature is None:
            SeqFeature.__init__(self, location, type)
        else:
            for key in feature.__dict__:
                if key in ["_start", "_end", "__digestion_topl", "_digestion_topr", "_digestion_bottomr", "_digestion_bottoml", "_seq", "subject", "query"]: 
                    pass
                else:
                    if key == "_original":
                        self._original = str(feature._original)
                    else:
                        self.__dict__[key] = copy.deepcopy(feature.__dict__[key]) 
        
        #start->end direction should be 5' to 3' on the top strand.
        if self.location.strand == -1:
            self._start = self.location.parts[-1].start.position
            self._end   = self.location.parts[0].end.position
        else:
            self._start = self.location.parts[0].start.position
            self._end   = self.location.parts[-1].end.position
        
        self._seq       = None
        self._qkey      = None #ID for features_dict
        self._second_id = None #ID for a single feature 
        self._start     = Qint(self._start)
        self._end       = Qint(self._end)
        self._start.parent = self
        self._end.parent   = self
        self._start.name   = "start"
        self._end.name     = "end"

        if subject is None:
            pass 
        else:
            self._digestion_topl    = "null" 
            self._digestion_topr    = "null"
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

        elif type(position[0]) is list and type(position[1]) is list:
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

class QUEEN():    
    """QUEEN object
    
    The QUEEN class can define a dsDNA object with sequence annotations. It can be 
    created by specifying a DNA sequence or importing a sequence data file in 
    GenBank or FASTA file format (single sequence entry). When a GenBank format file 
    is imported, its NCBI accession number, Addgene plasmid ID, or Benchling share 
    link can be provided instead of downloading the file to your local environment.
    
    """
    
    #Class variables that manage the execution histories of operational and search function
    dna_dict                   = {}
    queried_feature_dict       = {}
    queried_features_dict      = {}
    queried_features_name_dict = {} 
    process_description = None
    _namespace     = {}
    _products      = {} 
    _processes     = {}  
    _namespaceflag = 0 
    _num_history   = 1  
    _qnum    = 0
    _keep    = 1 
    _source  = None
    _project = None 
    
    def set_project(value=None):
        QUEEN._project = value 
                     
    def _get_genbank(_id, dbtype="ncbi"):
        """
        
        Dbtype can be selected from "ncbi", "addgene", "benchling".   
        For "ncbi", set NCBI accession number.  
        For "addgene", set plasmid ID. Sometimes different full sequence maps are 
        provided by the depositor and adgene, respectively, for a single plasmid.  
        In this case, please specify the plasmid ID followed by "addgene" or 
        "depositor" (Ex. 50005:addgene or 50005:depositor). If you set only plasmid ID, 
        the value will be specified as "plsmidID:addgene".  
        For "benchling", set a benchling shaared link
        
        """ 
        outb = io.BytesIO()
        outs = io.StringIO() 
        headers = {"User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"} 
        if dbtype == "ncbi":
            url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbwithparts&id={}&withparts=on".format(_id) 
        
        elif dbtype == "addgene": 
            if ":" not in _id:
                _id = _id + ":addgene"
            site = "https://www.addgene.org/{}/sequences/".format(_id.split(":")[0]) 
            html = requests.get(site)
            soup = BeautifulSoup(html.content, "html.parser")
            url  = soup.find(id="{}-full".format(_id.split(":")[1])).find(class_="genbank-file-download").get("href")

        elif dbtype == "benchling":
            if "https://benchling.com/" not in _id:
                raise ValueError("Please specify a proper benchling share link")
            url = _id + ".gb"
        
        elif dbtype == "googledrive":
            if "https://drive.google.com" not in _id:
                raise ValueError("Please specify a proper googledrive share link")
            match = re.search("https://drive\.google\.com/file/d/(.+)/view\?usp=sharing", url)
            fileid = match.group(1) 
            url = "https://drive.google.com/uc?export=download&id=" + fileid
        else: 
            raise ValueError("'datatype' can take only one of 'ncbi,' 'addgeen,' and 'benchling.'") 

        request = urllib.request.Request(url, headers=headers) 
        with urllib.request.urlopen(request) as u:
            outb.write(u.read())
        outs.write(outb.getvalue().decode())
        return outs

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
        obj = QUEEN(seq=self.seq, quinable=0)
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
    
    def __eq__(self, other):
        if not isinstance(other, QUEEN):
            return NotImplemented
        else:
            if self.seq == other.seq:
                pass
            else:
                return False, 1
                            
            if self.topology == "linear" and self.topology == other.topology:
                if self._left_end == other._left_end and self._left_end_top == other._left_end_top and self._left_end_bottom == other._left_end_bottom: 
                    pass
                else:
                    return False, 2

                if self._right_end == other._right_end and self._right_end_top == other._right_end_top and self._right_end_bottom == other._right_end_bottom:
                    pass
                else:
                    return False, 3

            if len(self.dnafeatures) == len(other.dnafeatures):
                flag = 0 
                for feat1, feat2 in zip(self.dnafeatures, other.dnafeatures):
                    if feat1.type == "source" or feat2.type == "source":
                        pass 
                    else:
                        if feat1 == feat2:
                            flag = 1
                        else:
                            flag = 0 
                    
                if flag == 1:
                    return True

                if flag == 0:
                    for feat1 in self.dnafeatures:
                        flag = 0 
                        for feat2 in other.dnafeatures:
                            if feat1 == feat2:
                                flag = 1
                                break
                            else:
                                pass 
                        if flag == 0:
                            return False, 4
                        else:
                            pass 
                    return True
            else:
                return False, 5

    def __repr__(self):
        if len(self.seq) > 50:
            out = "<queen.QUEEN object; project='{}', length='{} bp', topology='{}'>".format(self.project, len(self.seq), self.topology)
        else:
            out = "<queen.QUEEN object; project='{}', length='{} bp', sequence='{}', topology='{}'>".format(self.project, len(self.seq), self.seq, self.topology)
        return out 
    
    def __setattr__(self, key, value):
        if key == "_unique_id":
            if "_unique_id" in self.__dict__:
                if value in self.__class__.dna_dict:
                    if value.split("_")[-1].isdecimal() == True:
                        value = "_".join(value.split("_")[:-1]) 
                    unique = 0
                    while value + "_" + str(unique) in self.__class__.dna_dict:
                        unique += 1    
                    _unique_id = value + "_" + str(unique)
                else:         
                    _unique_id = value
                QUEEN.dna_dict[_unique_id] = None
                super.__setattr__(self, "_unique_id", _unique_id)
            else:
                QUEEN.dna_dict[value] = None
                super.__setattr__(self, "_unique_id", value)
        
        elif key == "_product_id":
            if value in self.__class__._products:
                if value.split("_")[-1].isdecimal() == True:
                    value = "_".join(value.split("_")[:-1]) 
                unique = 0
                while value + "_" + str(unique) in self.__class__._products:
                    unique += 1    
                _product_id = value + "_" + str(unique)
            else:         
                _product_id = value 
            self.__class__._products[_product_id] = None
            super.__setattr__(self, "_product_id", _product_id)
           
        elif key in ("seq", "rcseq", "dnafeatures", "productdict", "prcessdict", "project", "topology"): 
            raise AttributeError("'QUEEN' object attribute '{}' is read-only".format(key)) 
            #raise ValueError("Cannot assign to '{}' attribute.".format(key))   
        
        else:
            super.__setattr__(self, key, value) 
    
    def __getattribute__(self, name):
        if name == "features_dict":
            return dict(list(map(lambda x:(x._id, x), self.dnafeatures)))

        elif name == "_seq":
            qseq = Qseq(super().__getattribute__(name))
            if "_product_id" in self.__dict__:
                qseq.parental_id = self._product_id 
            elif "_unique_id" in self.__dict__:
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
                    histories.append((int(key.split("_")[-1]), history[0], history[1], history[2])) 
            histories.sort()
            return histories
        
        elif name == "seq":
            return self._seq
        
        elif name == "topology":
            return self._topology

        elif name == "rcseq":
            rcseq = self._seq.upper().translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
            rcseq = Qseq(rcseq)
            if "_product_id" in self.__dict__:
                rcseq.parental_id = self._product_id 
            elif "_unique_id" in self.__dict__:
                rcseq.parental_id = self._unique_id 
            rcseq.parent = self 
            rcseq.parental_class = "QUEEN"
            rcseq.name = "rcseq"
            return rcseq 
        
        elif name == "project":
            return self._product_id
        
        elif name == "sequence":
            return self._seq 
        
        elif name == "dnafeatures":
            return self._dnafeatures
        
        elif name == "productdict":
            if self._load_history == 1 or self._load_history == -1:
                pass 
            elif self._load_history == 0:
                try:
                    quine(self, execution=True, _io=False)
                except:
                    pass
            return dict(zip(self._productids, [QUEEN._products[key] for key in self._productids]))

        elif name == "processdict":
            return dict(zip(self._processids, [QUEEN._processes[key] for key in self._processids]))
        
        else:
            raise AttributeError("QUEEN obejct has no attribute '{}'".format(name))
    

    def __init__(self, seq=None, record=None, fileformat=None, dbtype="local", topology="linear", ssdna=False, import_history=True, supfeature=False, project=None, product=None, process_name=None, process_description=None, 
        pd=None, pn=None, process_id=None, original_ids=[], quinable=True, **kwargs):
        
        """

        Parameters
        ----------
        seq : str
            DNA seqeunce of `QUEEN_ object` to be generated. 
        record : str
            Source GenBank or FASTA file path of `QUEEN_object` to be generated.  
            When a GenBank format file is imported, its NCBI accession number, Addgene 
            plasmid ID, or Benchling share link can be provided instead of downloading 
            the file to your local environment.
        fileformat : str ("fasta" or "genbank"), default: None
            If the value is specified the file specified on `record` interpereted as 
            `fileformat`. Otherwise, the file format automatically detected based on 
            its file name.  
        dptype : str ("local", "NCBI", "addgene", or "benchling") 
            Online database location of the GenBank file to be imported.  
        topology : str ("linear" or "circular") 
            Sequence topology. When a `QUEEN_object` is created by loading from a GenBank 
            file, the topology is set according to the description in the GenBank file. 
        ssdna : bool, default: False
            If True, QUEEN object will handled as ssDNA. ssDNA QUEEN object cannot be 
            processed with modify ends and, be joined with dsDNA QUEEN object. By annealing 
            ssDNA QUEEN object trough `joindna` function, dsDNA QUEEN object can be generated.
        import_history : bool, default: True
            If False, it disable the inheritance of operational process histories of previously 
            generated `QUEEN_objects` to a newly producing `QUEEN_object`.   
        product : str 
            This parameter enables users to provide label names for producing `QUEEN_objects`.  
            The provided labels are stored in `QUEEN_objects.project`.  
            If tha value is not specified, and a `QUEEN_object` is created from a GenBank or 
            FASTA file, its sequence ID will be inherited here.  

        Attribuites
        -----------
        project : str
            Project name of `QUEEN_object` construction. In QUEEN, this property is also used 
            as a dictionary key to access the `.productdict` described below.   
            If a `QUEEN_object` is created from a GenBank or FASTA file, its sequence ID will 
            be inherited here. Otherwise, the project name is automatically generated or 
            defined based on the `product` value  to be unique amongst the existing 
            `.productdict` keys.
        seq : str  
            Top strand sequence (5′→3′). This property cannot be directly edited; only the 
            built-in operational functions of QUEEN described below can edit this property.
        rcseq : str  
            Bottom strand sequence (5′→3′). This property cannot be directly edited; only the 
            built-in operational functions of QUEEN described below can edit this property.
        topology : str ("linear" or "circular")   
            Sequence topology. When a `QUEEN_object` is created by loading from a GenBank file, 
            the topology is set according to the description in the GenBank file. 
            Only the built-in operational functions of QUEEN described below can edit this 
            property.
        dnafeatures : list of `DNAfeature_objects`  
            When a `QUEEN_object` is loaded from a GenBank file, `.dnafeatures` will 
            automatically be generated from the GenBank file's sequence features.  
            Otherwise,`.dnafeatures` will be an empty list. Each `DNAfeature_object` with the 
            following attributes provides an annotation for a given range of DNA sequence in 
            a `QUEEN_object`.
            
            - feature_id : str    
              Unique identifier. It is automatically determined to each feature when a `QUEEN_object` 
              is loaded from a GenBank file. 
            - feature_type : str    
              Biological nature. Any value is acceptable. The GenBank format requires registering a 
              biological nature for each feature.
            - start : int  
              Start position of `DNAfeature_object` in `QUEEN_object`. 
            - end : int
              Start position of `DNAfeature_object` in `QUEEN_object`.
            - strand : int (1 or -1) 
              Direction of `DNAfeature_object` in `QUEEN_object`. 
              Top strand (`1`) or bottom strand (`-1`).  
            - sequence : str  
              Sequence of the `DNAfeature_object` for its encoded direction.
            - qualifiers : dict 
              Qualifiers. When a GenBank file is imported, qualifiers of each feature will be 
              registered here. Qualifier names and values will serve as dictionary keys and 
              values, respectively. 
            
            `DNAfeature_object` can be edited only by the `editfeature()` function described 
            below. DNAfeature class is implemented as a subclass of the Biopython SeqFeature 
            class. Therefore, apart from the above attributes, DNAfeature class inherits all 
            the attributes and methods of SeqFeature class. For details about SeqFeature class,  
            see (https://biopython.org/docs/dev/api/Bio.SeqFeature.html) 
        productdict : dict  
            Dictionary for all of the inherited `QUEEN_objects` used to construct the present 
            `QUEEN_object`. The `.project` of each `QUEEN_object` serves as a key of this 
            dictionary.

        """
        
        if ssdna == True and topology == "circular":
            raise TypeError("The parameters of 'topology=circulr' and 'ssdna=True' cannot be specified simultaneously.")

        fseq      = seq 
        frecord   = record
        fdbtype   = dbtype  
        fproject  = project 
        ftopology = topology
        fproduct  = product 

        project = project if product is None else product
        process_name        = pn if process_name is None else process_name
        process_description = pd if process_description is None else process_description
       
        self._seq               = None
        self.record             = None
        self._dnafeatures       = None
        self._right_end         = None 
        self._left_end          = None
        self._history_feature   = None
        self._ssdna             = ssdna
        self._productids        = []
        self._processids        = []
        if seq is None and record is None:
            if project is None: 
                project = "dna"
            self._topology           = topology
            if "_" in project:
                project = project.replace("_","-")
            self._left_end_top      = 1
            self._left_end_bottom   = 1
            self._right_end_top     = 1 
            self._right_end_bottom  = 1
        
        elif seq is None or "." in seq:
            if "." in str(seq):
                record = seq
            
            if type(record) == str:
                if dbtype == "local":
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
                
                elif dbtype in ("ncbi", "addgene", "benchling", "googledrive"):
                    fmt    = "genbank"
                    record = QUEEN._get_genbank(record, dbtype)
                    with tempfile.TemporaryFile(mode="w+") as o:
                        content = record.getvalue() 
                        o.write(content)
                        o.seek(0)  
                        record = SeqIO.parse(o,fmt)
                        record  = next(record)

            elif type(record) == SeqRecord:
                record = record 
            else:
                record = SeqIO.parse(record,None)

            self._seq = str(record.seq).upper()
            self.record = record            
            if "topology" in record.annotations:
                self._topology = record.annotations["topology"]
            else:
                self._topology = topology
            
            if ssdna == False:
                self._right_end = ""
                self._left_end  = ""
                self._left_end_top      = 1
                self._left_end_bottom   = 1
                self._right_end_top     = 1 
                self._right_end_bottom  = 1
            else:
                self._right_end = self._seq
                self._left_end  = self._seq
                self._left_end_top      =  1
                self._left_end_bottom   = -1
                self._right_end_top     =  1 
                self._right_end_bottom  = -1

            if project is None:
                if record.id == "" or record.id == ".":
                    project = frecord.split("/")[-1].split(".")
                    project = project[0] if len(project) == 1 else ".".join(project[:-1])
                else:
                    project = record.id 
            self._unique_id         = project
            
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
                                feat.qualifiers[key][0] = feat.qualifiers[key][0].replace(" ","") 
                                feat.qualifiers[key][2] = feat.qualifiers[key][2].replace(" ","") 
                                results = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]", history) 
                                for result in results:
                                    _unique_id = result.split("['")[1][:-2] 
                                    QUEEN.dna_dict[_unique_id] = None
                                history_num = int(key.split("_")[-1]) 
                                pairs.append((feat, history_num, feat.qualifiers[key]))
                                history_feature_id = "0" 
                                history_feature = copy.deepcopy(feat) 
                    else:
                        for key in feat.qualifiers:
                            if key == "broken_feature":
                                feat.qualifiers["broken_feature"][0] = feat.qualifiers["broken_feature"][0].replace(" ","") 
                                note = feat.qualifiers["broken_feature"][0]
                                original = note.split(":")[-3]
                                feat._original = original
               
                if len(pairs) == 0:
                    import_history = False

                if import_history == True:
                    for pair in pairs:
                        del history_feature.qualifiers["building_history_{}".format(pair[1])]
                    for pair in pairs:
                        feat = pair[0]
                        new_history_num = pair[1] + QUEEN._num_history
                        history_feature.qualifiers["building_history_{}".format(new_history_num)] = pair[2]
                        history_nums.append(new_history_num) 
                     
                else:
                    #QUEEN.process_description = process_description
                    deletehistory(self)    
                
                QUEEN._num_history = max(history_nums)   
                if history_feature is not None:
                    self._history_feature = history_feature
                    self._dnafeatures.remove(pairs[0][0]) 
            
            if len(self.dnafeatures) == 0:
                import_history = False
                self._dnafeatures = []

        elif record is None:
            import_history = False
            if project is None:
                project = "dna"

            seq = seq.upper()
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
                            raise TypeError("An invalid nucleotide sequence pattern was found.")
                        
                        nucl_set_top    = list(set(list(top)))
                        nucl_set_bottom = list(set(list(bottom)))
                        if (len(nucl_set_top) == 1 and nucl_set_top[0] == "-"): 
                            self._ssdna = True
                            seq = bottom[::-1] 
                            top = bottom[::-1] 
                            bottom = "-" * len(bottom)  
                        
                        elif (len(nucl_set_bottom) == 1 and nucl_set_bottom[0] == "-"):
                            self._ssdna = True 
                        else:
                            pass 
                    else: 
                        sticky = False
                else:
                    if ssdna == True:
                        top    = seq 
                        bottom = "-" * len(seq) 
                        sticky = True 
                    else:
                        pass 

                self._seq       = str(seq).upper()
                if "_" in project:
                    project = project.replace("_","-")
                
                self._unique_id =  project
                if Alphabet:
                    self.record = SeqRecord(Seq(str(seq),Alphabet.DNAAlphabet()))
                else:
                    self.record = SeqRecord(Seq(str(seq)))
                
                self._dnafeatures = []
                self._topology   = topology
                
                if sticky == True:
                    self._topology   = "linear"
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
                raise TypeError("An invalid nucleotide sequence pattern was found.")
        self._supfeatureids()
        
        if type(supfeature) in (tuple, list) and type(supfeature[0]) == dict:
            for feature_dict in supfeature: 
                self.setfeature(feature_dict) 
        elif type(supfeature) == dict:
            self.setfeature(supfeature)
        
        if quinable == True:
            self._product_id = self._unique_id if product is None else product
            if import_history == False:                
                fseq         = "" if fseq is None else "seq='{}'".format(fseq)
                frecord      = "" if frecord is None else "record='{}'".format(frecord) if fseq == "" else ", record='{}'".format(frecord)
                fdbtype      = "" if fdbtype == "local" else ", dbtype='{}'".format(fdbtype)
                fproject     = "" if fproject is None else ", project='{}'".format(fproject)
                fssdna       = "" if ssdna == False else ", ssdna='{}'".format(ssdna)
                ftopology    = "" if topology == "linear" else ", topology='{}'".format(topology)
                fproduct     = "" if fproduct is None else ", product='{}'".format(fproduct)
                fileformat   = "" if fileformat is None else ", fileformat='{}'".format(fileformat) 
                fsupfeature  = "" if supfeature == False else ", supfeature={}".format(str(supfeature)) 
                process_name = "" if process_name is None else ", process_name='" + process_name + "'"
                process_description = "" if process_description is None else ", process_description='" + process_description + "'" 
                
                args = [fseq, frecord, fdbtype, fproject, fssdna, ftopology, fileformat, fsupfeature, fproduct, process_name, process_description]
                building_history = "QUEEN.dna_dict['{}'] = QUEEN({}{}{}{}{}{}{}{}{}{})".format(self._product_id, *args)  
                process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                QUEEN._num_history += 1 
                
                feat = SeqFeature(FeatureLocation(0, len(self.seq), strand=1), type="source") 
                feat._id = "0"
                feat.qualifiers["label"]       = [self.project] 
                feat.qualifiers["description"] = ["Record of building history"]
                feat.qualifiers["building_history_{}".format(QUEEN._num_history)] = [building_history.replace(" ","–"), "", ",".join([process_id] + original_ids)] 
                feat = DNAfeature(feature=feat, subject=self)         
                self._history_feature = feat
            else:
                nums = []
                for key in self._history_feature.qualifiers: 
                    if "building_history" in key[0:18]:
                        nums.append(int(key.split("_")[-1])) 
                nums.sort() 
                for key in self._history_feature.qualifiers:
                    if "building_history" in key[0:18]:
                        num = int(key.split("_")[-1]) 
                        if num != nums[-1]:
                            if self._history_feature.qualifiers[key][1] == "":
                                self._history_feature.qualifiers[key][1] = "_source: " + self.project + " construction"
                            else:
                                self._history_feature.qualifiers[key][1] += "; _source: " + self.project + " construction"
                        else:
                            if self._history_feature.qualifiers[key][1] == "":
                                self._history_feature.qualifiers[key][1] = "_source: " + self.project + " construction" + "; _load: " + self.project
                            else:
                                self._history_feature.qualifiers[key][1] += "; _source: " + self.project + " construction" + "; _load: " + self.project
        
        if QUEEN._keep == 1 and quinable == True: 
            QUEEN._products[self._product_id] = self
        
        if quinable == 1 and import_history == True:
            self._load_history = 0
        else:
            self._load_history = -1 
     
        self._positions       = tuple(range(len(self.seq))) 
        self.record.feartures = self.dnafeatures
        self.seq.parental_id  = self._unique_id
        if product is None:
            pass 
        else:
            QUEEN._namespace[product] = self

    def searchsequence(self, query, start=0, end=None, strand=2, unique=False, product=None, process_name=None, process_description=None, pn=None, pd=None, quinable=True, **kwargs):
        """Search for specific sequences from a QUEEN object.
        
        Search for specific sequences from a user-defined region of a `QUEEN_object` and 
        return a list of `DNAfeature_objects`. Start and end attributes of returned 
        `DNAfeature_objects` represent the sequence regions of the `QUEEN_object` that 
        matched the user's query. Note that the returned `DNAfeature_objects` will not be 
        generated with `.feature_id` and reflected to the parental `QUEEN_object`.   
        The returned `DNAfeature_objects` can be added to `QUEEN_object.dnafeatures` by 
        `editfeature()` with the `createattribute` option as explained below.

        Parameters
        ----------
        query : regex or str, default: ".+"  
            Search query sequence. If the value is not provided, the user-specified search 
            region of the `QUEEN_object` sequence with `start` and `end` explained below will 
            be returned. It allows fuzzy matching and regular expression. For details,   
            see https://pypi.org/project/regex/.  
            All IUPAC nucleotide symbols can be used. Restriction enzyme cut motif 
            representation can be used to define a query with `"^"` and `"_"` or 
            `"(int/int)"`.   For example, EcoRI cut motif can be provided by `"G^AATT_C"`, 
            where `"^"` and `"_"` represent the cut positions on the top and bottom strands, 
            respectively, or by `"GAATTC(-5/-1)"` or `"(-5/-1)GAATTC"`, where the left and 
            right integers between the parentheses represent the cut positions on the top and 
            bottom strands, respectively. Similarly, the cut motif of a Type-IIS restriction 
            enzyme BsaI can be given by `"GGTCTCN^NNN_N"`, `"N^NNN_NGAGACC"`, `"GGTCTC(1/5)"` 
            or `"(5/1)GAGACC"`. The returned `DNAfeature_objects` obtained for a query 
            restriction enzyme cut motif will hold the cutting rule in the 
            `"qualifier:cutsite"` attribute, which can be added to `QUEEN_object.dnafeatures` 
            by `editfeature()` with the `createattribute` option as explained below.  
            Regular expression is disabled for restriction enzyme cut motifs.  
        start : int (zero-based indexing), default: 0  
            Start position of the target range of the `QUEEN_object` sequence for the search. 
        end : int (zero-based indexing; default: the last sequence position of `QUEEN_object`)  
            End position of the target range of the `QUEEN_object` sequence for the search. 
        strand : int: 1 (top strand only), -1 (bottom strand only), or 2 (both strands), default: 2  
            Sequence strand to be searched.
        unique : bool (True or False), default: False
            If the value is `True` and multiple (more than a single) sequence region are detected 
            in the search, it will raise error. If False, multiple seaquence detections could be 
            acceptable.  
        
        Returns
        -------
        list of QUEEN.qobj.DNAfeature object
        
        """
        kwargs.setdefault("_sourcefile", None) 
        kwargs.setdefault("process_id", None)
        kwargs.setdefault("original_ids", []) 
        _sourcefile  = kwargs["_sourcefile"]  
        process_id   = kwargs["process_id"] 
        original_ids = kwargs["original_ids"]

        process_name        = pn if process_name is None else process_name
        process_description = pd if process_description is None else process_description
        
        history_features = [self._history_feature] 
        if "__dict__" in dir(query) and "cutsite" in query.__dict__:
            query = query.cutsite
        qorigin = query
        start   = 0 if start == len(self.seq) else start
        end     = len(self.seq) if end is None else end
        strand  = 2 if strand is None else strand 
        if start == 0 and end == len(self.seq):
            if self.topology == "circular":
                subject = self.seq + self.seq 
            else:
                subject = self.seq
        else:
            subject = self.printsequence(start, end, strand=1)
            

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
                cutsite, query, topl, topr, bottoml, bottomr = compile_cutsite(query)  
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
                        if match["strand"] == -1:
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
                        new_feat.qualifiers["cutsite"] = [Qseq(cutsite)]
                    
                    if type(qorigin) == Qseq and qorigin.parental_class == "Cutsite":
                        new_feat.qualifiers["label"] = [qorigin.parent.name]
                    feat_list.append(new_feat)  
                
                else:
                    pass 
                    
        qkey = self._unique_id + "_f" + str(QUEEN._qnum)
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
        
        if quinable == True:
            if type(qorigin) == Qseq:
                if qorigin.parental_class == "DNAfeature":
                    history_features.append(qorigin.parent.subject._history_feature)
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
                    history_features.append(qorigin.parent._history_feature)
                    parental_id = qorigin.parental_id 
                    if qorigin.name != None: 
                        if "printsequence" in qorigin.name:
                            if len(qorigin.name.split("_")) == 2: 
                                seqname = "QUEEN.dna_dict['{}'].printsequence(strand={})".format(parental_id, qorigin.name.split("_")[-1]) 
                            else:
                                seqname = "QUEEN.dna_dict['{}'].printsequence(start={}, end={}, strand={})".format(parental_id, *qorigin.name.split("_")[1:])
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

                elif qorigin.parental_class == "Cutsite":
                    if qorigin.parent.name not in cs.defaultkeys:
                        cs.new_cutsites.add((qorigin.parent.name, qorigin.parent.cutsite)) 
                    if qorigin.name == "cutsite":
                        qorigin = "cs.lib['{}'].{}".format(qorigin.parent.name, qorigin.name) 
                    else:
                        qorigin = "cs.lib['{}']".format(qorigin.parent.name) 
                else:
                    qorigin = "{}".format(repr(qorigin)) 
            else:
                qorigin = "{}".format(repr(qorigin)) 
            
            if len(history_features) > 1:
                history_feature = _combine_history(self, history_features) 
                self._history_feature = history_feature
            
            funique  = "" if unique == False else ", unique=True"
            fproduct = "" if product is None else ", product='{}'".format(product)
            process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
            process_description = "" if process_description is None else ", process_description='{}'".format(process_description)
            
            if start == 0 and end == len(self.seq):
                if strand == 2:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}{}{}{}{})".format(qkey, self._product_id, qorigin, funique, fproduct, process_name, process_description)
                    process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                    add_history(self, [building_history, "query: {}".format(qorigin), ",".join([process_id] + original_ids)], _sourcefile)
                else:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, strand={}{}{}{}{})".format(qkey, self._product_id, qorigin, strand, funique, fproduct, process_name, process_description)
                    process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                    add_history(self, [building_history, "query: {}; strand: {}".format(qorigin, strand), ",".join([process_id] + original_ids)], _sourcefile)
            else:
                building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchsequence(query={}, start={}, end={}, strand={}{}{}{}{})".format(qkey, self._product_id, qorigin, start, end, strand, funique, fproduct, process_name, process_description)   
                process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                add_history(self, [building_history, "query: {}; start: {}; end: {}; strand: {}".format(qorigin, start, end, strand), ",".join([process_id] + original_ids)], _sourcefile)
            QUEEN._qnum += 1 
        
        if product is None:
            pass 
        else:
            QUEEN._namespace[product] = feat_list
        
        if unique == True and len(feat_list) > 1:
            raise ValueError("Mutiple sequence regions were detected. If `unique` is True, the search results should be unique.") 

        return feat_list

    def searchfeature(self, key_attribute="all", query=".+", source=None, start=0, end=None, strand=2, product=None, process_name=None, process_description=None, pn=None, pd=None, quinable=True, **kwargs):
        
        """Search for `DNAfeature_objects` holding a queried value in a designated `key_attribute`.
        
        Parameters
        ----------
        key_attribute : str, default: "all" 
            Attribute type to be searched (`feature_id`, `feature_type`, `"qualifier:*"`, 
            or `sequence`). If the value is not provided, it will be applied to all of the 
            attributes in the `QUEEN_object`, excluding `sequence`. However, if the `query` 
            value is provided with only the four nucleotide letters (A, T, G, and C), this 
            value will be automatically set to `sequence`.
        query : regex or str, default: ".+"
            Query term. `DNAfeature_objects` that have a value matches to this query for 
            `key_attribute` designated above will be returned. It allows fuzzy matching and 
            regular expression. For details, see https://pypi.org/project/regex/. If the 
            `key_attribute` is `sequence`, all IUPAC nucleotide symbols can be used.
        source : list of `DNAfeature_objects`, default: `QUEEN_object.dnafeatures`
            Source `DNAfeature_objects` to be searched. `DNAfeature_objects` outside the 
            search range defined by `start`, `end`, and `strand` will be removed from the 
            source. Any `DNAfeature_objects` can be provided here. For example, a list of 
            `DNAfeature_objects` _returned from another `searchsequence()` or `searchfeature()` 
            operation can be used as the source to achieve an AND search with multiple queries.
        start : int (zero-based indexing), default: 0 
            Start position of the target range of the `QUEEN_object` sequence for the search. 
        end : int (zero-based indexing), default: the last sequence position of `QUEEN_object`
            End position of the target range of the `QUEEN_object` sequence for the search. 
        strand : int 1 (top strand), -1 (bottom strand), or 2 (both strands), default: 2  
            Sequence strand to be searched.  
        process_name : str
            For detailes, see `QUEEN.queen.qobj.__doc__`. 
        process_description : str
            For detailes, see `QUEEN.queen.qobj.__doc__`.
        process_id : str 
            For detailes, see `QUEEN.queen.qobj.__doc__`.
        quinable : bool
            For detailes, see `QUEEN.queen.qobj.__doc__`.

        Returns
        ------- 
        list of QUEEN.qobj.DNAfeature object

        """
        kwargs.setdefault("_sourcefile", None) 
        kwargs.setdefault("process_id", None)
        kwargs.setdefault("original_ids", []) 
        _sourcefile  = kwargs["_sourcefile"]  
        process_id   = kwargs["process_id"] 
        original_ids = kwargs["original_ids"]

        process_name        = pn if process_name is None else process_name
        process_description = pd if process_description is None else process_description
        
        history_features = [self._history_feature] 
        start   = 0 if start == len(self.seq) else start
        end     = len(self.seq) if end is None else end
        qkey = self._unique_id + "_f" + str(QUEEN._qnum)
        features = editfeature(self, key_attribute=key_attribute, query=query, source=source, start=start, end=end, strand=strand, target_attribute=None, operation=None, quinable=0, process_description=process_description) 
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
        
        if quinable == 1:
            if type(query) == Qseq:
                if query.parental_class == "DNAfeature":
                    history_features.append(query.parent.subject._history_feature)
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
                    history_features.append(query.parent._history_feature)
                    parental_id = query.parental_id 
                    if query.name != None:
                        if "printsequence" in query.name:
                            if len(query.name.split("_")) == 2: 
                                seqname = "QUEEN.dna_dict['{}'].printsequence(strand={})".format(parental_id, query.name.split("_")[-1]) 
                            else:
                                seqname = "QUEEN.dna_dict['{}'].printsequence(start={}, end={}, strand={})".format(parental_id, *query.name.split("_")[1:])
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

                elif query.parental_class == "Cutsite":
                    if query.parent.name not in cs.defaultkeys:
                        cs.new_cutsites.add((query.parent.name, query.parent.cutsite)) 
                    query = "cs.lib['{}'].{}".format(query.parent.name, query.name) 
                else:
                    query = "{}".format(repr(query))
            else:
                query = "{}".format(repr(query))  
            
            if source is not None:
                qkeys = set([]) 
                for feat in source:
                    if "_qkey" in feat.__dict__:
                        qkeys.add(feat._qkey)
                if len(set(qkeys)) == 1:
                    source = "QUEEN.queried_features_dict['{}']".format(list(qkeys)[0])
                else:
                    pass
            
            fproduct            = "" if product is None else ", product='" + product + "'"
            process_name        = "" if process_name is None else ", process_name='" + process_name + "'"
            process_description = "" if process_description is None else ", process_description='" + process_description + "'" 

            if len(history_features) > 1:
                history_feature = _combine_history(self, history_features) 
                self._history_feature = history_feature

            if start == 0 and end == len(self.seq):
                if strand == 2 and source is None:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}{}{}{})".format(qkey, self._product_id, key_attribute, query, fproduct, process_name, process_description)
                    process_id, oiginal_ids = make_processid(self, building_history, process_id, original_ids)
                    add_history(self, [building_history, "key_attribute: {}; query: {}".format(key_attribute, query), ",".join([process_id] + original_ids)], _sourcefile)
                elif strand == 2:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}, source={}{}{}{})".format(qkey, self._product_id, key_attribute, query, source, fproduct, process_name, process_description)
                    process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                    add_history(self, [building_history, "key_attribute: {}; query: {}; : {}".format(key_attribute, query, source), ",".join([process_id] + original_ids)], _sourcefile)

                else:
                    building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}, source={}, strand={}{}{}{})".format(qkey, self._product_id, key_attribute, query, source, strand, fproduct, process_name, process_description)
                    process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                    add_history(self, [building_history, "key_attribute: {}; query: {}; : {}; strand: {}".format(key_attribute, query, source, strand), ",".join([process_id] + original_ids)], _sourcefile)

            else:
                building_history  = "QUEEN.queried_features_dict['{}'] = QUEEN.dna_dict['{}'].searchfeature(key_attribute='{}', query={}, source={}, start={}, end={}, strand={}{}{}{})".format(qkey, self._product_id, key_attribute, query, source, start, end, strand, fproduct, process_name, process_description) 
                process_id, original_ids = make_processid(self, building_history, process_id, original_ids)
                add_history(self, [building_history, "key_attribute: {}; query: {}; : {}; start: {}; end: {}; strand: {}".format(key_attribute, query, source, start, end, strand), ",".join([process_id] + original_ids)], _sourcefile) 
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
                subdna = cropdna(self, start, end, quinable=0)
                             
            else:
                raise TypeError("slice indices must be integers or None or have an __index__ method")
            
            if strand == -1 or strand < 0:
                return flipdna(subdna, quinable==0)
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
            raise ValueError("Cicularized QUEEN object cannot be joined with others.") 
        else:
            return joindna(self, other, quinable=0)

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
            return joindna(other, self, quinable=0) 
    
    def printsequence(self, start=None, end=None, strand=2, hide_middle=None, linebreak=None, display=False):
        """Returns and displays partial or the entire dsDNA sequence and sequence end structures.
        
        Parameters
        ----------
        start : int (zero-based indexing), default: 0  
            Start position of the sequence.   
        end : int (zero-based indexing),  default: the last sequence position of `QUEEN_object`
            End position of the sequence.  
        strand : int: 1 (top strand only), -1 (bottom strand only), or 2 (both strands), default: 2  
            Sequence strand(s) to be returned.
        display : bool (True or False), default: True   
            If `True`, the output will be displayed in `STDOUT`.  
        hide_middle: int or None, default: None  
            Length of both end sequences to be displayed.  
        linebreak: int (default: length of the `QUEEN_object` sequence)  
            Length of sequence for linebreak.

        Returns
        -------
        If `strand` == `1` or `-1`,  
        sequence of the defined strand (5’→3’)    
        If `strand` == `2`,  
        "str/str" (top strand sequence (5’→3’)/bottom strand sequence (3’→5’))  

        """

        if self._ssdna == True: 
            if display == True:
                print("5' {} 3'".format(self.seq)) 
            return self.seq

        whole = False
        if linebreak is None:
            width = len(self.seq) + 1
        else:
            width = linebreak

        if hide_middle is None or hide_middle > 0.5 * len(self.seq):
            hide_middle = int(0.5 * len(self.seq)) 
            whole = True

        if start is None and end is None and strand == 2:
            tl = len(self._left_end)  if self._left_end_top  == -1 else 0
            tr = len(self._right_end) if self._right_end_top == -1 else 0
            truetop = "-" * tl + self.seq[tl:(-1*tr if tr != 0 else None)] + "-" * tr
            
            bl = len(self._left_end)  if self._left_end_bottom  == -1 else 0
            br = len(self._right_end) if self._right_end_bottom == -1 else 0
            truebottom = "-" * bl + self.rcseq[::-1][bl:(-1*br if br != 0 else None)] + "-" * br

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
                if self.topology == "circular":
                    top = self.seq[start:] + self.seq[:end]
                else:
                    return ""
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
        
        if start is None and end is None: 
            if strand == 1 or strand == 0:
                return_seq = truetop
            elif strand == -1:
                return_seq = truebottom[::-1] 
            elif strand == 2:
                return_seq = truetop + "/" + truebottom
        
        else:
            if strand == 1 or strand == 0:
                return_seq = top
            elif strand == -1:
                return_seq = bottom[::-1]
            elif strand == 2:
                return_seq = top + "/" + bottom
        
        return_seq = Qseq(return_seq) 
        return_seq.__dict__ = self.seq.__dict__ 
        if start == 0 or None and end == len(self.seq) or None:
            return_seq.name = "printsequence_{}".format(strand)
        else:
            return_seq.name = "printsequence_{}_{}_{}".format(start, end, strand)
        return return_seq

    def _supfeatureids(self):
        keys = [feat._id for feat in self.dnafeatures if "_id" in feat.__dict__] 
        keys = tuple(keys) 
        
        for i in range(0, len(self.dnafeatures)):
            if "_id" in self._dnafeatures[i].__dict__ and (self._dnafeatures[i]._id.isdecimal()==False and keys.count(self._dnafeatures[i]._id)>0):
                pass 
            else:
                self._dnafeatures[i]._id = str(i*100)

    def _check_uniqueness(self):
        keys_list = [feat._id for feat in self.dnafeatures if "_id" in feat.__dict__] 
        keys_list = tuple(keys_list) 
        keys_set  = set(keys_list) 
        if len(keys_list) == len(keys_set):
            pass 
        else:
            for i in range(0, len(self.dnafeatures)):
                if "_id" in self._dnafeatures[i].__dict__ and keys_list.count(self._dnafeatures[i]._id)==1:
                    pass 
                else:
                    j = 0
                    while str(i*100 + j) in keys_set:
                        j += 1
                    self._dnafeatures[i]._id = str(i*100 + j)
    
    def getdnafeatures(self,feature_id):
        return self.features_dict[str(feature_id)] 
    
    def setfeature(self, attribute_dict=None):
        """Create a new feature in the QUEEN object.
        
        Parameters  
        ----------
        attribute_dict : dict 
            Dictionaly with key-value pairs of the attributes of DNAfeature objects: 
            "feature_id", "feature_type", "start", "end", "strand", and  "qualifier:*".  
            The following attributes have default value, so if they are not specified in the 
            dictionary, the value would be set with the default values. 
            - `feature_id` : `str`, 
              The default value is a random unique ID which is not used in `.dnafeatures` of 
              the QUEEN object.  
	    - `feature_type` : `str` (default: `"misc_feature"`) 
	    - `start` : `int` (default: 0) 
	    - `end` : `int` (default: length of the `QUEEN_object` sequence)
	    - `strand` : `int` (-1, 0 or 1, default: 1) 

        Returns 
        -------
        self.dnafeatures: list of DNAfeature object
    
        """
        if attribute_dict is None:
            attribute_dict = {}
        elif type(attribute_dict) == dict:
            attribute_dict = copy.copy(attribute_dict) 

        attribute_dict.setdefault("feature_type", "misc_feature")
        attribute_dict.setdefault("start", 0)  
        attribute_dict.setdefault("end", len(self.seq))  
        attribute_dict.setdefault("strand", 1)
         
        if type(attribute_dict["start"]) == int:
            if attribute_dict["start"] > attribute_dict["end"] and self.topology == "circular":
                locations = [[attribute_dict["start"], len(self.seq), attribute_dict["strand"]], [0, attribute_dict["end"], attribute_dict["strand"]]]
                if attribute_dict["strand"] == -1:
                    locations.reverse()
                    new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type=attribute_dict["feature_type"]) 
                else:
                    new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type=attribute_dict["feature_type"])
            else:
                new_feat = SeqFeature(FeatureLocation(attribute_dict["start"], attribute_dict["end"], strand=attribute_dict["strand"]), type=attribute_dict["feature_type"])
        
        elif type(attribute_dict["start"]) in (tuple, list):
            locations = [] 
            for s, e in zip(attribute_dict["start"], attribute_dict["end"]):
                if s > e:
                    if self.topology == "circular":
                        locations.append((s, len(self.seq), attribute_dict["strand"])) 
                        locations.append((0, e, attribute_dict["strand"]))  
                    else:
                        raise ValueError("'end' position must be larger than 'start' position.") 
                else:
                    locations.append((s, e, attribute_dict["strand"])) 
            
            if attribute_dict["strand"] == -1:
                locations.reverse()
                new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type=attribute_dict["feature_type"]) 
            else:
                new_feat = SeqFeature(CompoundLocation(list(map(lambda x:FeatureLocation(*x), locations))), type=attribute_dict["feature_type"])
        
        for key in attribute_dict:
            if key[:len("qualifier:")] == "qualifier:":
                subkey = key[len("qualifier:"):]
                if subkey not in new_feat.qualifiers:
                    new_feat.qualifiers[subkey] = [attribute_dict[key]]
                else:
                    new_feat.qualifiers[subkey].append(attribute_dict[key])
            else:
                pass 

        features_dict = self.features_dict 
        new_feat = DNAfeature(feature=new_feat, subject=self) 
        
        if "feature_id" not in attribute_dict:
            feature_ids = [feat.feature_id for feat in self.dnafeatures if feat.feature_id.isdecimal()] 
            
            if len(feature_ids) == 0:
                new_id = str(1)  
            
            for _id in feature_ids:
                start  = features_dict[_id].start
                new_id = str(int(_id) + 1)
                if attribute_dict["start"] > start:
                    idnum = 1
                    while new_id in features_dict: 
                        new_id = str(int(_id) + 1 + idnum)
                        idnum += 1
                    break
        
            new_feat._id = new_id
        
        else:
            if attribute_dict["feature_id"] in set(features_dict.keys()):
                raise ValueError("feature_id value should be uniqule in .dnafeatures.") 
            else:
                new_feat._id = attribute_dict["feature_id"]

        self._dnafeatures.append(new_feat)

    def printfeature(self, feature_list=None, attribute=None, separation=None, seq=False, output=None, x_based_index=0):
        """ Print a tidy data table of annotation features/attributes of `QUEEN_object`. 
        
        Default output attributes are `"feature_id"`, `"feature_type"`, 
        `"qualifier:label"`, `"start"`, `"end"`, and `"strand"`.
    
        Parameters    
        ----------
        feature_list : list of QUEEN.qobj.DNAfeature objects (default: self.dnafeatures)  
            List of features to be displayed in the output table. If not given, all features 
            held by the QUEEN_object will be the subject.
        attribute : list of feature attributes  
            The default value is `["feature_id", "feature_type", "qualifier:label", "start", 
            "end", "strand"]`.   List of feature attributes to be displayed in the output 
            table. If the value is `"all"`, it will generate a table for all the attributes 
            held by the `QUEEN_object` except for `"sequence"`.
        seq : bool (True or False), default: `False`  
            If `True`, the sequence of each feature for its encoded direction will be 
            displayed in the output table.   
        separation : str, default: space(s) to generate a well-formatted table   
            String to separate values of each line.
        output : str, default: STDOUT   
            Output file name. 
        x_based_index : 0 or 1, default: 0 
            As a default, positions of all features are given in the zero-based indexing in 
            QUEEN (same as Python). If this parameter is set to `1`, they will be shown in 
            the 1-based indexing (as seen in the GenBank format).
        
        Returns
        -------
        None

        """
        
        _ids        = ["feature_id"] 
        types       = ["feature_type"] 
        labels      = ["qualifier:label"]
        starts      = ["start"] 
        ends        = ["end"] 
        strands     = ["strand"]
        sequences   = ["sequence"]
        sep         = separation
        seqflag     = seq
        detail      = False
        others_dict = {}
        

        if attribute is None:
            attribute = ["feature_id", "feature_type", "qualifier:label", "start", "end", "strand"]
        
        elif attribute == "all" or attribute == ["all"]:
            attribute = ["feature_id", "feature_type", "qualifier:label", "start", "end", "strand"]
            detail    = True

        else:
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
        
        else:
            if sep == ",":
                import csv
                output = csv.writer(sys.stdout) 
            else:
                pass 

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

    def outputgbk(self, output=None, format="genbank", record_id=None, annotation=None, export_history=True, _return=False):
        """Output `QUEEN_object` to a GenBank file. 

        In addition to all of the `DNAfeature_objects` in the input `QUEEN_object`, 
        a `DNAfeature_object` encoding the entire construction processes that generated 
        the `QUEEN_object` in `qualifier:building_history` will also be output to the 
        GenBank file.  
        
        Parameters
        ----------
        output : str, default: STDOUT
            Output file name. 
        format : str, default: "genbank"
            Output file format
        record_id : str, default: None
            Record ID of a output file, If the value is `None`, self.project value is used as `record_id`.
        annotation : dict, default: None
            Dictionary of annotations for the genbank.   
            For details, please see https://biopython.org/docs/latest/api/Bio.SeqRecord.html.
        export_history : bool, default: True
            If False, construnction history of the `QUEEN_object` will not be output.

        Returns
        -------
        None
        
        """

        handle = output
        export_input = False
        separate_history = False
        
        stdIOflag = 0 
        if handle is None:
            stdIOflag = 1
            handle    = io.StringIO()
       
        product_dict    = {}
        histories       = quine(self, _return=True) 
        history_nums    = [history[0] for history in histories]
        history_feature = copy.deepcopy(self._history_feature)
        
        remove_keys = [] 
        for key in history_feature.qualifiers:
            if "building_history" in key[0:18]:
                num = int(key.split("_")[-1]) 
                if num in history_nums:
                    process_id = history_feature.qualifiers[key][2].split(",")[0].replace(" ", "")
                    if "-" in process_id:
                        pass 
                    else:
                        history_feature.qualifiers[key][2] = self.project + "-" + self._history_feature.qualifiers[key][2]
                else:
                    remove_keys.append(key) 
        
        for key in remove_keys:
            del history_feature.qualifiers[key] 

        features = copy.deepcopy(self.dnafeatures)
        if export_history is True:
            features.append(history_feature) 
        
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
                
            if separate_history is not False and type(separate_history) is str and export_history == True:
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
                self.record.annotations["source"] = os.getcwd() + "/" + separate_history
            
            if export_input == True:
                for hisoty in histories:
                    if "QUEEN(record=" in history[1] and ("dbtype='local'" in history or "dbtype" not in history[1]):
                        match = re.search("QUEEN.dna_dict\['([^\[\]]+)'\]", history[1])
                        if match is not None:
                            key = match.group(1)
                            if dans[0].productdict[key] is not None:
                                product_dict[key] = dnas[0].prodcutdict[key]
                                product_dict[key].record.annotations["keyword"]    = "QUEEN input"
                                product_dict[key].record.annotations["accession"] = re.search("record='([^=]+)'[,\)]", history[1]).group(1) 
                        else:
                            pass 
                    else:
                        pass 

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
        
        if annotation is None:
            pass 
        else:
            self.record.annotations = annotation

        #Add DATE
        import datetime
        dt = datetime.datetime.now() 
        self.record.annotations["date"] = dt
        SeqIO.write(self.record, handle, format)
        self.record.features = self.dnafeatures
        if stdIOflag == 1:
            if _return == True:
                return(handle.getvalue()) 
            else:
                print(handle.getvalue(), end="") 

        if export_input == True and len(product_dict) > 0:
            for key in product_dict:
                value = outputgbk(product_dict[key], export_input=False, _return=True)
                if stdIOflag == 1:
                    #print("//") 
                    print(value, end="") 
                else:
                    #handle.write("//\n") 
                    handle.write(value) 

        if stdIOflag == 0:
            handle.close() 
