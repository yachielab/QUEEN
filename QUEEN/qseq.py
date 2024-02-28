class Qseq(str):
    def __init__(self, seq):
        self.qkey           = None
        self.parent         = None
        self.parental_id    = None
        self.parental_class = None
        self.name           = None 
        self.item           = None        
   
    def __contains__(self, key):
        if self.parental_class == "QUEEN" and self.parent is not None and self.parent.topology == "circular":
            return key in self + self
        else: 
            return key in str(self) 

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
                elif item.stop > len(self):
                    stop = item.stop - (len(self) * (item.stop // len(self)))
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
    
    def find(self, *args, **kwargs): 
        if self.parent.topology == "linear" or (self.parental_class == "QUEEN" and self.parent is not None and self.parent.topology == "circular"):
            value = (self + self).find(*args, **kwargs) 
            return value 
        else:
            return super().find(*args, **kwargs)
