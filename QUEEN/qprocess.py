class ProcessFlow():
    def get_depth(adict, target, depth=0):
        depth_list = []
        if target in adict:
            return depth
        else:
            depth += 1
            if len(list(adict.keys())) > 0:
                for key in adict:
                    depth_list.append(ProcessFlow.get_depth(adict[key], target, depth)) 
            else:
                return -1
        return max(depth_list)
    
    def __init__(self, tree, adict):
        self.tree         = tree
        self.process_dict = adict
        self.code         = [] 
        self.exec_dict    = {} 
        for key in self.process_dict:
            self.exec_dict[key] = 0  

    def decode(self):
        processid_depth_list = [] 
        for processid in self.process_dict:
            depth = ProcessFlow.get_depth(self.tree, processid) 
            processid_depth_list.append((depth, processid))
        
        processid_depth_list.sort() 
        for depth, processid in processid_depth_list:
            self.process_dict[processid].decode(self.process_dict, self.exec_dict, self.code)  
        return self.code

class qprocess:
    def __init__(self, product, process_name, process_description, process_id):
        self.funclabel = "entry"
        self.input    = [] 
        self.index    = 0 
        
        #Common arguments
        self.product = product
        self.process_name = process_name
        self.process_description = process_description
        self.process_id = process_id 
        self.arguments = set(["input", "product", "process_name", "process_description", "process_id"]) 

    #def __setattr__(self, key, value):
        #if key in self.arguements:
        #    super.__setattr__(self, key, value) 
        #else:
        #    raise KeyError("{}".format(key))   
    
    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]
        
    def decode(self, process_dict, exec_dict, codes):
        if exec_dict[self.process_id] == 1:
            pass
        else:
            parameters = [] 
            entry = []
            for key in self.input:
                key         = key[0] 
                index       = key[1] 
                funclabel   = process_dict[key]["funclabel"] 
                productname = process_dict[key]["product"] 
                if funclabel == "cutdna":
                    entry.append("{}[{}]".format(productname, index))
                else:
                    entry.append("{}".format(productname))
            entry = ", ".join(entry) 

            for key in self.arguments:
                #print(self.funclabel, key, self.__dict__[key]) 
                if key  == "input":
                    pass
                
                elif self.__dict__[key] == None:
                    pass 
                
                elif type(self.__dict__[key]) == tuple:
                    processid   = self.__dict__[key][0] 
                    index       = self.__dict__[key][1]
                    funclabel   = process_dict[processid]["funclabel"] 
                    productname = process_dict[processid]["product"] 
                    process_dict[processid].decode(process_dict, exec_dict, codes) 
                    
                    if funclabel == "cutdna":
                        parameters.append("{}={}[{}]".format(key, productname, index))
                    elif (funclabel == "searchsequence" or funclabel == "searchfeature") and index.isdecimal() == True:
                        parameters.append("{}={}[{}]".format(key, productname, index))
                    elif index is None or index.isdecimal() == True:
                        parameters.append("{}={}".format(key, productname)) 
                    else:
                        parameters.append("{}={}{}".format(key, productname, index))
                
                else:
                    if self.__dict__[key].isdecimal() == True or self.__dict__[key][0:2] == "*[":
                        parameters.append("{}={}".format(key, self.__dict__[key]))
                    else:
                        parameters.append("{}='{}'".format(key, self.__dict__[key]))
            
            parameters = ", ".join(parameters) 
            if self.funclabel == "entry":
                code = self.product + " = " + "QUEEN" + "(" + parameters + ")" 
            elif self.funclabel == "searchsequence" or self.funclabel == "searchfeature":
                code = self.product + " = "+ entry + "."  + self.funclabel + "(" + parameters + ")"
            else:
                code = self.product + " = " + self.funclabel + "(" + entry + ", " + parameters + ")"

            exec_dict[self.process_id] = 1
            codes.append(code) 

class qentry(qprocess):
    def __init__(self, product, process_name, process_description, process_id): 
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "entry" 
        self.seq       = None
        self.record    = None
        self.dbtype    = None
        self.topology  = None
        self.arguments = self.arguments | set(["seq", "record", "dbtype", "topology"]) 
    
class qcutdna(qprocess): 
    def __init__(self, product, process_name, process_description, process_id): 
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "cutdna" 
        self.positions = None
        self.arguments = self.arguments | set(["positions"]) 

class qcropdna(qprocess): 
    def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "cropdna" 
        self.start = None
        self.end   = None
        self.arguments = self.arguments | set(["start", "end"]) 

class qjoindna(qprocess):
    def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "joindna" 
        self.topology  = None
        self.arguments = self.arguments | set(["topology"]) 

class qmodifyends(qprocess): 
    def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "modifyends" 
        self.left  = None
        self.right = None
        self.arguments = self.arguments | set(["left", "right"]) 

class qflipnda(qprocess):
    def __init__(self, product, process_name, process_description, process_id): 
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "flipdna"
        pass 
    
class qeditsequence(qprocess): 
     def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "editsequence"
        self.source_sequence = None
        self.destination_sequence = None
        self.start  = None
        self.end    = None
        self.strand = None
        self.arguments = self.arguments | set(["source_sequence", "destination_sequence", "start", "end", "strand"]) 

class qeditsequence(qprocess): 
    def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel        = "editfeature"
        self.key_attribute    = None
        self.query            = None
        self.source           = None
        self.start            = None
        self.end              = None
        self.strand           = None
        self.target_attribute = None 
        self.operation        = None
        self.arguments = self.arguments | set(["key_attribute", "query", "source", "start", "end", "strand", "target_attribute", "operation"]) 

class qsearchsequence(qprocess): 
     def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel = "searchsequence"
        self.query = None
        self.start = None
        self.end   = None
        self.arguments = self.arguments | set(["query", "start", "end"]) 

class qsearchfeature(qprocess): 
    def __init__(self, product, process_name, process_description, process_id):
        super().__init__(product, process_name, process_description, process_id) 
        self.funclabel        = "searchfeature"
        self.key_attribute    = None
        self.query            = None
        self.source           = None
        self.start            = None
        self.end              = None
        self.strand           = None
        self.arguments = self.arguments | set(["key_attribute", "query", "source", "start", "end", "strand"]) 

