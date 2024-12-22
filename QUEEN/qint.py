class Qint(int):
    def __init__(self, num):
        self.qkey        = None
        self.parent      = None
        self.parental_id = None
        self.name        = None
        self._qint       = True
