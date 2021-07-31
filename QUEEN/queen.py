import os 
import functools
import sys
sys.path.append("/".join(__file__.split("/")[:-1]))

from qobj import QUEEN
import qfunction 

def set_namespace(_globals=None):
    if _globals is None:
        _namespace = None
        QUEEN._namespace     = {}
        QUEEN._namespaceflag = 0
    else:
        _namespace = _globals
        QUEEN._namespace     = _globals
        QUEEN._namespaceflag = 1

_namespace   = None
QUEEN._namespace = {}
replaceattribute = qfunction.replaceattribute
createattribute  = qfunction.createattribute
removeattribute  = qfunction.removeattribute
quine        = qfunction.quine
cutdna       = qfunction.cutdna
cropdna      = qfunction.cropdna
modifyends   = qfunction.modifyends
joindna      = qfunction.joindna
flipdna      = qfunction.flipdna
editsequence = qfunction.editsequence
editfeature  = qfunction.editfeature
visualizedna = qfunction.visualizedna 
