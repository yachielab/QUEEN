class Qseq(str):
    def __init__(self, seq):
        self.qkey           = None
        self.parent         = None
        self.parental_id    = None
        self.parental_class = None
        self.name           = None 
        self.item           = None        
    
    """
    def __eq__(self, other):
        if type(other) == str:
            return super().__eq__(other) 

        elif not isinstance(other, Qseq):
            return NotImplemented
        
        else:
            flag1 = 0 
            if self.parental_class == "QUEEN" and self.parent is not None and self.parent.topology == "circular":
                flag1 = 1 
            
            flag2 = 0 
            if other.parental_class == "QUEEN" and other.parent is not None and other.parent.topology == "circular":
                flag2 = 1  
            
            if flag1 == 1 and flag2 == 1 and len(self) == len(other) and str(self) != str(other):
                flag = 0 
                for i in range(1 , len(self)):
                    if str(self[i:] + self[:i]) == str(other):
                        return True 
                    else:
                        pass 
                
                if flag == 0:
                    return False
            else:
                return super().__eq__(other) 
    """

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
        value.parent         = self.parent
        value.parental_id    = self.parental_id
        value.parental_class = self.parental_class 
        value.name           = self.name
        if value.item is None:
            value.item = item
        else:
            value.item = slice(value.item.start + item.start, value.item.start + item.stop)
        return value
    

