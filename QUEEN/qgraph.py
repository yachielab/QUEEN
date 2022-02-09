import os 
import sys 
import collections 
import regex as re
from graphviz import Digraph
sys.path.append("/".join(__file__.split("/")[:-1]))
import quine
from qprocess import * 

def make_newhistories(histories, search_function=True):
    new_histories = [] 
    name_dict              = {} 
    unique_name_dict       = {}
    unique_name_name_dict  = {}
    process_notes          = set([])   
    for history in histories:
        process_id = history[3].split(",")[0]
        info       = history[2] 
        history    = history[1]

        match0  = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = (.*)",history)
        match1  = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = ", history)
        match2  = re.search("product='([^=]+)'", history)
    
        match_name = re.search("process_name='([^=]*)'", history) 
        match_description = re.search("process_description='([^=]*)'", history) 
        process_name = match_name.group(1) if match_name is not None else None
        process_description = match_description.group(1) if match_description is not None else None

        if match1 is not None: 
            match3 = list(re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", match1.group(1))) 
            if match2 is not None:
                name = match2.group(1)
                if "," not in name and len(match3) > 1:
                    for m, match in enumerate(match3):
                        name_dict[match.group(0)] = match2.group(1)  + "[{}]".format(m) 
                else:
                    for aname, match in zip(name.split(","), match3):
                        name_dict[match.group(0)] = name 
            else:
                for m, match in enumerate(match3): 
                    name_dict[match.group(0)] = match.group(1) 

            for m, match in enumerate(match3): 
                unique_name_dict[match.group(0)] = match.group(1) 
                unique_name_name_dict[match.group(1)] = name_dict[match.group(0)] 
            
            infos = [history, info] 
            infos.append(process_name)
            infos.append(process_description) 
            infos.append(process_id) 
            new_histories.append(infos) 
            process_notes.add((process_name, process_description)) 

        elif search_function==True and match0 is not None:
            match3 = re.search("QUEEN.queried_features_dict\['([^\[\]]+)'\]", match0.group(1))
            if match2 is not None:
                name_dict[match3.group(0)] = match2.group(1) 
            else:
                name_dict[match3.group(0)] = match3.group(1) 

            unique_name_dict[match3.group(0)] = match3.group(1) 
            unique_name_name_dict[match3.group(1)] = name_dict[match3.group(0)] 
            
            infos = [history, info] 
            infos.append(process_name)
            infos.append(process_description)  
            infos.append(process_id) 
            new_histories.append(infos) 
            process_notes.add((process_name, process_description)) 
        else:
            pass 
    return new_histories, process_notes, name_dict, unique_name_dict, unique_name_name_dict


def visualizeflow(*dnas, search_function=None, grouping=True, inherited_process=None, process_description=None, split_input=None, sf=None, ip=None, pd=None, si=None, alias_dict=None):
    inherited_process = ip if inherited_process is None else inherited_process
    inherited_process = False if inherited_process is None else inherited_process
    
    search_function = sf if search_function is None else search_function
    search_function = True if search_function is None else search_function
    
    if grouping == True:
        process_classification = True 
    else:
        process_classification = False 

    process_description = pd if process_description is None else process_description
    process_description = False if process_description is None else process_description
    
    split_input = si if split_input is None else split_input
    if split_input is None:
        if len(dnas) > 2:
            split_input = None
        else:
            split_input = True
    #split_input = True if split_input is None else split_input
    
    
    pd_visible = process_description
    if alias_dict is None:
        alias_dict = {}

    histories     = quine.quine(*dnas, _return=True)
    sdgs          = {} 
    sdgs_nodes    = collections.defaultdict(set) 
    clusters      = {} 
    
    new_histories, process_notes, name_dict, unique_name_dict, unique_name_name_dict = make_newhistories(histories, search_function=search_function) 
    for key in unique_name_dict:
        if unique_name_dict[key] in alias_dict:
            unique_name_dict[key] = alias_dict[unique_name_dict[key]] 
    
    for key in name_dict:
        if name_dict[key] in alias_dict:
            name_dict[key] = alias_dict[name_dict[key]] 


    process_notes = list(process_notes)
    inputs = [] 
    gbks   = set([]) 
    nodes  = set([])  
    dg    = Digraph(name="cluster_operation")
    dg.attr(rankdir='LR')
    dg.attr(fontname="arial") 
    dg.attr(nodesep="0.1")
    dg.attr(ranksep="0.2")
    sdgs["__main__"] = dg 
    
    sourcenames_all        = []
    productnames_all       = [] 
    removedproducts        = [] 
    product_funcname_dict  = {}
    unique_name_count_dict = collections.defaultdict(int)  
    for h, history in enumerate(new_histories):
        process_description = history[3] 
        process_name        = history[2]
        info        = history[1] 
        history     = history[0]
       
        matchi = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = QUEEN(\(record.*)", history) 
        matcho = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = (.*)", history) 
        matchs = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = (.*)",history)
        process_index = process_notes.index((process_name, process_description))

        if matcho is not None:
            source = matcho.group(2) 
        else:
            source = matchs.group(2)  
            
        for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", source):
            key = key.group(0)
            sourcename = unique_name_dict[key]
            unique_name_count_dict[sourcename] += 1
    
    for h, history in enumerate(new_histories):
        process_description = history[3] 
        process_name        = history[2]
        info        = history[1] 
        history     = history[0]
       
        if len(info) > 1:
            info = info.split("; ")
            info_dict = dict([item.split(": ") for item in info])
            for key1 in info_dict:
                if key1 ==  "_source":
                    flag = 0
                    for key2 in alias_dict:
                        if key2 == info_dict["_source"]:
                            info_dict["_source"] = info_dict["_source"].replace(key2, alias_dict[key2])
                            flag = 1
                            break

                    if flag == 0:
                        for key2 in alias_dict:
                            if key2 in info_dict["_source"]:
                                info_dict["_source"] = info_dict["_source"].replace(key2, alias_dict[key2])

        else:
            info_dict = None

        matchi = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = QUEEN(\(record.*)", history) 
        matcho = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = (.*)", history) 
        matchs = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = (.*)",history)
        process_index = process_notes.index((process_name, process_description))

        sourcenames = []
        if matchi is not None: 
            for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", matchi.group(1)):
                key = key.group(0) 
                sourcename = unique_name_dict[key]
                if sourcename not in gbks:
                    gbks.add(sourcename) 
        
        if matcho is not None:
            source = matcho.group(2) 
        else:
            source = matchs.group(2)  
            
        for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", source):
            key = key.group(0)
            sourcename = unique_name_dict[key]
            if sourcename in gbks:
                if split_input == True and unique_name_count_dict[sourcename] > 1:
                    sourcename = sourcename + "–" + str(process_index)
                if sourcename not in nodes and sourcename not in removedproducts:
                    dg.node(sourcename, label=name_dict[key], margin="0.05", shape="note", fontname="Arial") 
                    nodes.add(sourcename) 
                    if sourcename not in inputs:
                        inputs.append(sourcename) 
                else:
                    pass
            else:
                if sourcename not in nodes and sourcename not in removedproducts:
                    dg.node(sourcename, label=name_dict[key], margin="0.05", shape="oval", fontname="Arial") 
                    nodes.add(sourcename) 
                else:
                    pass
            sourcenames.append(sourcename) 
        
        sourcenames_all.extend(sourcenames) 
        if re.match("cutdna", source) is not None:
            funclabel = "cutdna"
            funcname  = "cutdna_{}".format(h)
        elif re.match("cropdna", source) is not None:
            funclabel = "cropdna"
            funcname  = "cropdna_{}".format(h)
        elif re.match("modifyends", source) is not None:
            funclabel = "modifyends"
            funcname  = "modifyends_{}".format(h)
        elif re.match("flipdna", source) is not None:
            funclabel = "flipdna"
            funcname  = "flipdna_{}".format(h)
        elif re.match("joindna", source) is not None:
            funclabel = "joindna"
            funcname  = "joindna_{}".format(h)
        elif re.match("editsequence", source) is not None:
            funclabel = "editsequence"
            funcname  = "editsequence_{}".format(h)
        elif re.search("searchsequence", source) is not None:
            funclabel = "searchsequence"
            funcname  = "searchsequence_{}".format(h)
        elif re.search("searchfeature", source) is not None:
            funclabel = "searchfeature"
            funcname  = "searchfeature_{}".format(h)
        else:
            funclabel = None
            funcname  = None
        
        if matcho is not None:
            productnames = [] 
            product = matcho.group(1)
            for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", product):
                key = key.group(0) 
                #if funclabel == "cutdna":
                #    print(key) 
                productname = unique_name_dict[key]
                productnames.append(productname)
                if info_dict is not None and "_load" in info_dict:
                    #print(info_dict, productname, history) 
                    if productname not in gbks:
                        gbks.add(productname) 
                        if productname not in inputs:
                            inputs.append(productname) 
                else:
                    pass

                if productname not in nodes:
                    if inherited_process == True:
                        if productname not in nodes and productname not in gbks:
                            dg.node(productname, label=productname, margin="0.05", shape="oval", fontname="Arial") 
                            nodes.add(productname) 
                    else:
                        if funclabel is not None:
                            removedproducts.append(productname) 
                else:
                    pass 
            
            if funclabel is not None:
                productnames_all.extend(productnames)  
        else:
            productnames = [] 
            product = matchs.group(1)
            for key in re.finditer("QUEEN.queried_features_dict\['[^\[\]]+'\]", product):
                key = key.group(0) 
                productname = unique_name_dict[key]
                productnames.append(productname) 
                if productname not in nodes:
                    if len(name_dict[key]) > 4:
                        margin = "0.05"
                    else:
                        margin = "0.0"
                    dg.node(productname, label=name_dict[key], margin=margin, shape="box", fontname="Arial") 
                    nodes.add(productname) 
                else:
                    pass
            productnames_all.extend(productnames)  
        if funcname is None:
            pass 
        
        else:
            if funclabel == "searchsequence" or funclabel == "searchfeature":
                temp = '<tr><td port="{}" border="1" align="left"><b> </b><i>{} </i> = {}</td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                      '<tr>',
                      #'<td port="func" border="1" bgcolor="#BDD0FC"><font color="black" point-size="16">{}<b> </b></font></td>',
                      '<td port="func" border="1" bgcolor="#66FFCC"><font color="black" point-size="16">{}<b> </b></font></td>'
                      '</tr>',
                      '{}',
                      '</table>>'])
                infotext = ""
                
                query_flag      = 0 
                query_match     = re.search("query=([^=]+),", source)
                queryname_match = re.search("QUEEN.dna_dict\['([^\[\]]+)'\]", query_match.group(1))
                if queryname_match is not None:
                    queryname = queryname_match.group(1)
                    info_dict["query"] = query_match.group(1).replace(queryname_match.group(0), name_dict[queryname_match.group(0)])
                    #queryname = sourcenames[1] 
                    query_flag = 1 
                    
                if info_dict is not None:
                    for key, value in info_dict.items():
                        if key != "_source" and key != "_load":
                            infotext += temp.format(key, key, value)

                
                if sourcenames[0] not in nodes:
                    dg.node(sourcenames[0], label=sourcenames[0], margin="0.05", shape="oval", fontname="Arial") 
                    nodes.add(sourcenames[0]) 
                
                dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
                dg.edge(sourcenames[0], funcname+":func", arrowhead="dot")

                if query_flag == 1:
                    if queryname not in productnames_all and unique_name_count_dict[queryname] > 1:
                        if queryname + "_search" not in nodes:
                            dg.node(queryname + "_search", label=queryname, margin="0.05", shape="oval", fontname="Arial") 
                            nodes.add(queryname + "_search") 
                        else:
                            pass
                        dg.edge(queryname + "_search", funcname+":query", arrowhead="odot")
                    else:
                        if queryname in product_funcname_dict:
                            if inherited_process == True:
                                if queryname in product_funcname_dict:
                                    if queryname not in nodes:
                                        dg.node(queryname, label=queryname, margin="0.05", shape="oval", fontname="Arial") 
                                        nodes.add(queryname)
                                    if "cutdna" in product_funcname_dict[queryname][0]:
                                        dg.edge(product_funcname_dict[queryname][0], queryname) 
                                    else:
                                        dg.edge(product_funcname_dict[queryname][0] + ":func", queryname) 
                                dg.edge(queryname, funcname+":query", arrowhead="dot")

                            else:
                                if "cutdna" in product_funcname_dict[queryname][0]:
                                    dg.edge(product_funcname_dict[queryname][0], funcname+":query", arrowhead="dot")
                                else:
                                    dg.edge(product_funcname_dict[queryname][0] + ":func", funcname+":query", arrowhead="dot")
                        else:
                            dg.edge(queryname, funcname+":query", arrowhead="odot")

                for product in productnames:
                    dg.edge(funcname+":func", product)

            elif funclabel == "joindna":
                temp = '<tr><td border="1" color="#FFFDC7" bgcolor="#FFFDC7" port="f{}"> </td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="0">',
                      '<tr>',
                      '<td>',
                      '<table cellpadding="0" cellspacing="0" border="0">',
                      '<tr>',
                      #'<td border="1" color="#FFCECB" bgcolor="#FFCECB" port="f0"> </td>',
                      '<td border="1" color="#FFFDC7" bgcolor="#FFFDC7" port="f0"> </td>',
                      #'<td port="func" rowspan="{}" border="1" color="#FFCECB" bgcolor="#FFCECB" align="left"><font color="black" point-size="16">joindna</font><b> </b></td>',
                      '<td port="func" rowspan="{}" border="1" color="#FFFDC7" bgcolor="#FFFDC7" align="left"><font color="black" point-size="16">joindna</font><b> </b></td>',
                      '</tr>'
                      '{}',
                      '</table>'
                      '</td>'
                      '</tr>'
                      '<tr>',
                      '<td border="1" align="left" height="20"><b> </b><i>topology </i> = \'{}\'<b> </b></td>'
                      '</tr>',
                      '</table>>'])
                
                sourcelabel = ""
                for s, source in enumerate(sourcenames[1:]):
                    sourcelabel += temp.format(s+1) 
            
                dg.node(funcname, label.format(len(sourcenames), sourcelabel, info_dict["topology"]), shape="plaintext", fontname="Arial", margin="0.05")
                for s, sourcename in enumerate(sourcenames):
                    if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                        if sourcename in product_funcname_dict:
                            if sourcename not in nodes:
                                dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                                nodes.add(sourcename)
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                            else:
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                        dg.edge(sourcename, funcname+":f"+str(s), arrowhead="dot")
                    else:
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], funcname+":f"+str(s), arrowhead="dot")
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":f"+str(s), arrowhead="dot")
                 
                product_funcname_dict[productname] = (funcname, process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 


            elif funclabel == "cropdna":
                temp = '<tr><td port="{}" border="1" align="left"><b> </b><i>{} </i> = {}</td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                      '<tr>',
                      '<td port="func" border="1" bgcolor="#FFFDC7"><font color="black" point-size="16">{}<b> </b></font></td>',
                      '</tr>',
                      '{}',
                      '</table>>'])
                
                if search_function == False:
                    infotext = ""
                    if info_dict is not None:
                        for key, value in info_dict.items():
                            if key != "_source" and key != "_load":
                                infotext += temp.format(key, key, value)
                    dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
                    
                    sourcename = sourcenames[0]
                    if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                        if sourcename in product_funcname_dict:
                            if sourcename not in nodes:
                                dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                                nodes.add(sourcename)
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                            else:
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                        dg.edge(sourcename, funcname+":func", arrowhead="dot")
                    else:  
                        dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")

                else:
                    startunique = None
                    start_match = re.search("start=(QUEEN.queried_features_dict\['[^\[\]]+'\])",source)
                    if start_match is not None:    
                        uniquename   = start_match.group(1) 
                        startunique  = unique_name_dict[uniquename]
                        startname    = name_dict[uniquename]
                        rsource      = source.replace(uniquename, startname) 
                        smatch       = re.search("start=([^=]+),", rsource) 
                        info_dict["start"] = smatch.group(1) 
                
                    endunique = None 
                    end_match = re.search("end=(QUEEN.queried_features_dict\['[^\[\]]+'\])",source) 
                    if end_match is not None:    
                        uniquename = end_match.group(1) 
                        endunique  = unique_name_dict[uniquename] 
                        endname    = name_dict[uniquename] 
                        rsource    = source.replace(uniquename, endname) 
                        ematch     = re.search("end=([^=]+)[,\)]", rsource) 
                        info_dict["end"] = ematch.group(1) 
                   
                    infotext = ""
                    if info_dict is not None:
                        for key, value in info_dict.items():
                            if key != "_source" and key != "_load":
                                infotext += temp.format(key, key, value)
                    
                    dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
                   
                    sourcename = sourcenames[0]
                    
                    #print(sourcenames) 
                    if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                        if sourcename in product_funcname_dict:
                            if sourcename not in nodes:
                                dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                                nodes.add(sourcename)
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                            else:
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                        dg.edge(sourcename, funcname+":func", arrowhead="dot")
                    else:
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")
                   
                    if startunique is not None:
                        dg.edge(startunique, funcname+":start", arrowhead="odot")
                    if endunique is not None:
                        dg.edge(endunique, funcname+":end", arrowhead="odot")
                 
                product_funcname_dict[productname] = (funcname, process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 

            elif funclabel == "cutdna":
                temp_pros = '<tr><td border="1" color="#00000" bgcolor="#FFFDC7" port="f{}">[{}]</td></tr>'
                temp = '<tr><td port="{}" border="1" align="left"><b> </b><i>{} = {}</td></tr>'     
                
                label="".join(['<<table border="0" cellborder="0" cellspacing="0" cellpadding="0">',
                      '<tr>',
                      '<td>',
                      '<table cellpadding="0" cellspacing="0" border="0">',
                      '<tr>',
                      '<td port="func" rowspan="{}" border="1" color="#000000" bgcolor="#FFFDC7" align="center"><font color="black" point-size="16">cutdna</font><b> </b></td>',
                      '<td border="1" color="#000000" bgcolor="#FFFDC7" port="f0">[0]</td>',
                      '</tr>'
                      '{}',
                      '</table>'
                      '</td>'
                      '</tr>'
                      '{}',
                      '</table>>'])
                
                productlabel = ""
                for s, productname in enumerate(productnames[1:]):
                    productlabel += temp_pros.format(s+1, s+1) 

                #label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                #      '<tr>',
                #      '<td port="func" border="1" bgcolor="#FFCECB"><font color="black" point-size="16">{}<b> </b></font></td>',
                #      '</tr>',
                #      '{}',
                #      '</table>>'])
                
                infotext = ""
                if search_function == False:
                    for key, value in info_dict.items():
                        if key  == "positions":
                            for v, pos in enumerate(value.split(",")):
                                infotext += temp.format(key+"_"+str(v), key+"</i>[{}]".format(v), pos)
                        else:
                            pass
                    dg.node(funcname, label.format(len(productnames), productlabel, infotext), shape="plaintext", fontname="Arial")

                    sourcename = sourcenames[0] 
                    if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                        if sourcename in product_funcname_dict:
                            if sourcename not in nodes:
                                dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                                nodes.add(sourcename)
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                            else:
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                        dg.edge(sourcename, funcname+":func", arrowhead="dot")
                    else:
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")
                else:   
                    cutsitenames = [] 
                    uniquenames  = [] 
                    rsource = source
                    for pos_match in re.finditer("(QUEEN.queried_features_dict\['[^\[\]]+'\])(\[[0-9]+\])",source):    
                        uniquename  = pos_match.group(1) 
                        cutsitename = name_dict[uniquename] 
                        rsource     = rsource.replace(uniquename, cutsitename)
                        uniquename  = unique_name_dict[uniquename]
                        uniquenames.append(uniquename) 
                        #uniquenames.append(uniquename+"_"+pos_match.group(2)[1:-1]) 
                        cutsitenames.append(cutsitename+pos_match.group(2))
                    rsource  = rsource.split(",") 
                    elements = [] 
                    for element in rsource[1:]:
                        if "=" in element:
                            break
                        else:
                            elements.append(element.replace(" ","")) 

                    info_dict["positions"] = ",".join(elements)  
                    for key, value in info_dict.items():
                        if key  == "positions":
                            for v, pos in enumerate(value.split(",")):
                                infotext += temp.format(key+"_"+str(v), key+"</i>[{}]".format(v), pos)
                        else:
                            pass
                    
                    dg.node(funcname, label.format(len(productnames), productlabel, infotext), shape="plaintext", fontname="Arial")
                    #dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
                    
                    sourcename = sourcenames[0] 
                    if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                        if sourcename in product_funcname_dict:
                            if sourcename not in nodes:
                                dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                                nodes.add(sourcename)
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                            else:
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                        dg.edge(sourcename, funcname+":func", arrowhead="dot")
                    else:
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")
                    
                    for key, value in info_dict.items():
                        if key  == "positions":
                            for v, pos in enumerate(value.split(",")):
                                for cutsitename, uniquename in zip(cutsitenames, uniquenames):
                                    if cutsitename in pos:
                                        dg.edge(uniquename, funcname+":" + key + "_" + str(v), arrowhead="odot")
                                else:
                                    pass 
                        else:
                            pass
                
                for p, productname in enumerate(productnames):
                    product_funcname_dict[productname] = (funcname + ":f" + str(p), process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 

            elif funclabel == "modifyends" and len(sourcenames) > 1:  
                temp = '<tr><td border="1" align="left" port="{}"><b> </b><i>{} </i> = {}</td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                      '<tr>',
                      '<td port="func" border="1" bgcolor="#FFFDC7"><font color="black" point-size="16">{}<b> </b></font></td>',
                      '</tr>',
                      '{}',
                      '</table>>'])
                
                infotext = ""
                matchl = re.search("left=\.*(QUEEN.dna_dict\['[^\[\]]+'\])", history) 
                if matchl is not None:
                    uniquename = re.search("QUEEN.dna_dict\['([^\[\]]+)'\]", info_dict["leftobj"]).group(0)
                    objname    = name_dict[uniquename]
                    infotext += temp.format("left", "left", info_dict["leftobj"].replace(uniquename, objname))
                else:
                    infotext += temp.format("left", "left", info_dict["left"])
                
                matchr = re.search("right=\.*(QUEEN\.dna[^=]*)", history)
                if matchr is not None:
                    uniquename = re.search("QUEEN.dna_dict\['([^\[\]]+)'\]", info_dict["rightobj"]).group(0)
                    objname    = name_dict[uniquename]
                    infotext += temp.format("right", "right", info_dict["rightobj"].replace(uniquename, objname))
                else:
                    infotext += temp.format("right", "right", info_dict["right"])

                dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
                
                for sourcename in sourcenames:
                    if matchl is not None and sourcename in matchl.group(1):
                        dg.edge(sourcename, funcname+":left", arrowhead="odot")
                    elif matchr is not None and sourcename in matchr.group(1):
                        dg.edge(sourcename, funcname+":right", arrowhead="odot")
                    else:
                        if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                            if sourcename in product_funcname_dict:
                                if sourcename not in nodes:
                                    dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                                    nodes.add(sourcename)
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                            else:
                                dg.edge(sourcename, funcname+":func", arrowhead="dot")
                        else:  
                            if "cutdna" in product_funcname_dict[sourcename][0]:
                                dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                            else:
                                dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")
 
                product_funcname_dict[productname] = (funcname, process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 

            elif funclabel == "modifyends":
                temp = '<tr><td border="1" align="left"><b> </b><i>{} </i> = {}</td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                      '<tr>',
                      '<td port="func" border="1" bgcolor="#FFFDC7"><font color="black" point-size="16">{}<b> </b></font></td>',
                      '</tr>',
                      '{}',
                      '</table>>'])
                infotext = ""
                if info_dict is not None:
                    for key in ["left", "right"]:
                        value = info_dict[key] 
                        infotext += temp.format(key, value)
                dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
               
                sourcename = sourcenames[0] 
                if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                    if sourcename in product_funcname_dict:
                        if sourcename not in nodes:
                            dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                            nodes.add(sourcename)
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                    dg.edge(sourcename, funcname+":func", arrowhead="dot")
                else:
                    if "cutdna" in product_funcname_dict[sourcename][0]:
                        dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                    else:
                        dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")

                product_funcname_dict[productname] = (funcname, process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 

            elif funclabel == "flipdna":
                temp = '<tr><td border="1" align="left"><b> </b><i>{} </i> = {}</td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                      '<tr>',
                      '<td port="func" border="1" bgcolor="#FFFDC7"><font color="black" point-size="16">{}</font></td>',
                      '</tr>',
                      '{}',
                      '</table>>'])
                
                infotext = ""
                if info_dict is not None:
                    for key, value in info_dict.items():
                        if key != "_source" and key != "load":
                            infotext += temp.format(key, value)
                dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
               
                sourcename = sourcenames[0] 
                if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                    if sourcename in product_funcname_dict:
                        if sourcename not in nodes:
                            dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                            nodes.add(sourcename)
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                    dg.edge(sourcename, funcname+":func", arrowhead="dot")
                else:
                    if "cutdna" in product_funcname_dict[sourcename][0]:
                        dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                    else:
                        dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")
                
                if info_dict is None:
                    product_funcname_dict[productname] = (funcname, process_name, process_description, None) 
                else:
                    product_funcname_dict[productname] = (funcname, process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 

            else:
                temp = '<tr><td border="1" align="left"><b> </b><i>{} </i> = {}</td></tr>'
                label="".join(['<<table border="0" cellborder="1" cellspacing="0" cellpadding="1">',
                      '<tr>',
                      '<td port="func" border="1" bgcolor="#FFFDC7"><font color="black" point-size="16"><b></b>{}<b> </b></font></td>',
                      '</tr>',
                      '{}',
                      '</table>>'])
                infotext = ""
                if info_dict is not None:
                    for key, value in info_dict.items():
                        if key != "_source" and key != "_laod":
                            infotext += temp.format(key, value)
                dg.node(funcname, label.format(funclabel, infotext), shape="plaintext", fontname="Arial")
                

                sourcename = sourcenames[0] 
                if inherited_process == True or (sourcename not in product_funcname_dict) or (process_name != product_funcname_dict[sourcename][1] or process_description != product_funcname_dict[sourcename][2]):
                    if sourcename in product_funcname_dict:
                        if sourcename not in nodes:
                            dg.node(sourcename, label=sourcename, margin="0.05", shape="oval", fontname="Arial") 
                            nodes.add(sourcename)
                        if "cutdna" in product_funcname_dict[sourcename][0]:
                            dg.edge(product_funcname_dict[sourcename][0], sourcename) 
                        else:
                            dg.edge(product_funcname_dict[sourcename][0] + ":func", sourcename) 
                    dg.edge(sourcename, funcname+":func", arrowhead="dot")
                else:
                    if "cutdna" in product_funcname_dict[sourcename][0]:
                        dg.edge(product_funcname_dict[sourcename][0], funcname+":func", arrowhead="dot")
                    else:
                        dg.edge(product_funcname_dict[sourcename][0] + ":func", funcname+":func", arrowhead="dot")
                
                product_funcname_dict[productname] = (funcname, process_name, process_description, info_dict["_source"] if "_source" in info_dict else None) 

            if info_dict is not None and "_source" in info_dict:
                name = info_dict["_source"]  
                if name not in sdgs:
                    sdgs[name] = dg.subgraph(name="cluster_source{}".format(len(list(sdgs.keys()))))
                sdg_obj = sdgs[name]
                sdg     = sdgs[name].graph
                sdg.attr(style='dashed')
                sdg.attr(rankdir='LR')
                sdg.attr(label=name) 
                sdgflag = 1
                for sourcename in sourcenames:
                    sdgs_nodes[name].add(sourcename) 
                
                for productname in productnames:
                    if "_load" in info_dict:
                        if productname not in nodes:
                            dg.node(productname, label=productname, margin="0.05", shape="oval", fontname="Arial") 
                            nodes.add(productname)
                        #dg.edge(funcname+":func", productname, arrowhead="dot")
                    sdgs_nodes[name].add(productname)
            else:
                sdg = dg 
                sdgflag = 0 
            
            
            if process_classification == True:
                if process_description is not None or process_name is not None: 
                    key = "cluster_{}".format(process_index)
                    if key not in clusters.keys():
                        clusters[key] = sdg.subgraph(name=key)
                        #with clusters[key] as subg:
                        subg = clusters[key].graph
                        if pd_visible == True:
                            subg.attr(style="solid")
                            subg.attr(labeljust="l")
                            text = process_name + ": " + process_description if process_name is not None else process_description
                            text = text.split(" ") 
                            charnum = 0
                            new_text = ""
                            for unit in text:
                                new_text += unit
                                if charnum > 60:
                                    new_text += "\l"
                                    charnum = 0 
                                else:
                                    new_text += " "
                                charnum += len(unit) 
                            subg.attr(label=new_text+"\l") 
                        else:
                            subg.attr(style="solid")
                            subg.attr(label=process_name if process_name is not None else process_description + "\l") 
                        subg.attr(rankdir='LR') 
                    
                    #with clusters[key] as subg:
                    subg = clusters[key].graph
                    subg.node(funcname)
                    for sourcename in sourcenames:
                        if "–" in sourcename:
                            gbkname      = "–".join(sourcename.split("–")[:-1]) 
                            source_index = sourcename.split("–")[-1] 
                            if gbkname in gbks and source_index == str(process_index):
                                subg.node(sourcename)
                        
                        elif sourcename in product_funcname_dict: 
                            if product_funcname_dict[sourcename][1] == process_name and product_funcname_dict[sourcename][2] == process_description and inherited_process == True:
                                subg.node(sourcename)

                    if funclabel == "modifyends" and len(sourcenames) > 1:
                        for sourcename in sourcenames:
                            if sourcename not in product_funcname_dict:
                                subg.node(sourcename)

                    if funclabel in ["searchsequence", "searchfeature"]: 
                        for productname in productnames:
                            subg.node(productname)
                    else: 
                        pass
                    
                    if funclabel in ["searchsequence", "searchfeature"] and query_flag == 1:
                        if queryname + "_search" in nodes:
                            subg.node(queryname + "_search") 
            
            #if sdgflag == 1:                
            #    sdg_obj.parent.subgraph(sdg_obj.graph)
            #else:
            #    pass 

    for name in clusters: 
        clusters[name].parent.subgraph(clusters[name].graph)

    for name in sdgs:
        if name != "__main__":
            for nodename in list(sdgs_nodes[name]):
                if nodename in nodes:
                    sdgs[name].graph.node(nodename) 
            sdgs[name].parent.subgraph(sdgs[name].graph)
    
    for productname in product_funcname_dict:
        if productname not in sourcenames_all:
            funcname            = product_funcname_dict[productname][0] 
            process_name        = product_funcname_dict[productname][1]
            process_description = product_funcname_dict[productname][2]
            importname          = product_funcname_dict[productname][3]
            if "cutdna" in funcname:
                sdg_flag = 0 
                if importname is not None and importname != "__main__":
                    sdg_flag = 1 
                    sdg_obj  = sdgs[importname]
                    sdg      = sdg_obj.graph
                else:
                    sdg = dg
                process_index = process_notes.index((process_name, process_description)) 
                key = "cluster_{}".format(process_index)
                if inherited_process == True: 
                    clusters[key] = sdg.subgraph(name=key)
                    with clusters[key] as subg:
                        subg.attr(style='dashed') 
                        subg.node(productname, label=productname, margin="0.05", shape="oval", fontname="Arial")
                        subg.edge(funcname, productname, arrowhead="dot") 
                    if sdgflag == 1:                
                        sdg_obj.parent.subgraph(sdg_obj.graph)
            else:
                dg.node(productname, label=productname, margin="0.05", shape="oval", fontname="Arial") 
                dg.edge(funcname + ":func", productname) 
    
    for node in nodes:
        if "–" in node:
            key = "–".join(node.split("–")[:-1])
            if key in nodes:
                 dg.edge(key, node)

    if process_classification == False:
        subdg = dg.subgraph() 
        with subdg as c:
            c.attr(rankdir='LR')
            c.attr(rank='same')
            if len(inputs) > 1:
                for i in range(len(inputs)):
                    c.node(inputs[i])
                #subdg.parent.subgraph(subdg.graph) 
    return dg

def add_key(adict, parent=None, child=None):
    if parent is None:
        adict[child] = {}
        return adict  
    elif parent not in adict:
        for key in adict:
            adict[key] = add_key(adict[key], parent, child) 
    else:
        adict[parent][child] = {} 
    return adict

def get_depth(adict, target, depth=0):
    depth_list = []
    if target in adict:
        return depth
    else:
        depth += 1
        if len(list(adict.keys())) > 0:
            for key in adict:
                depth_list.append(get_depth(adict[key], target, depth)) 
        else:
            return -1
    return max(depth_list)

def generate_processflow(*dnas):
    histories = quine.quine(*dnas, _return=True)
    new_histories, process_notes, name_dict, unique_name_dict, unique_name_name_dict = make_newhistories(histories) 
    
    process_tree = {}
    process_dict = {}  
    product_funcname_dict = {} 
    for h, history in enumerate(new_histories):
        process_id          = history[4] 
        process_description = history[3] 
        process_name        = history[2]
        info                = history[1] 
        history             = history[0]
        
        matchi = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = QUEEN(\(record.*)", history) 
        matcho = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = (.*)", history) 
        matchs = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = (.*)",history)
        matchp = re.search("product='([^=]+)'[,\)]", history) 
         
        if matcho is not None:
            source = matcho.group(2) 
        else:
            source = matchs.group(2)  
         
        if matchp is None:
            product = None
        else:
            product = matchp.group(1) 

        if re.match("cutdna", source) is not None:
            process_dict[process_id] = qcutdna(product, process_name, process_description, process_id) 
        elif re.match("cropdna", source) is not None:
            process_dict[process_id] = qcropdna(product, process_name, process_description, process_id) 
        elif re.match("modifyends", source) is not None:
            process_dict[process_id] = qmodifyends(product, process_name, process_description, process_id) 
        elif re.match("flipdna", source) is not None:
            process_dict[process_id] = qflipdna(product, process_name, process_description, process_id) 
        elif re.match("joindna", source) is not None:
            process_dict[process_id] = qjoindna(product, process_name, process_description, process_id) 
        elif re.match("editfeature", source) is not None:
            process_dict[process_id] = qeditfeature(product, process_name, process_description, process_id) 
        elif re.match("editsequence", source) is not None:
            process_dict[process_id] = qeditsequence(product, process_name, process_description, process_id) 
        elif re.search("searchsequence", source) is not None:
            process_dict[process_id] = qsearchsequence(product, process_name, process_description, process_id) 
        elif re.search("searchfeature", source) is not None:
            process_dict[process_id] = qsearchfeature(product, process_name, process_description, process_id) 
        else:
            process_dict[process_id] = qentry(product, process_name, process_description, process_id) 
        
        sourcenames = set([])
        for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", source):
            key = key.group(0)
            sourcename = unique_name_dict[key]
            sourcenames.add(sourcename) 
        
        connected_keys = [] 
        subsource_dict = {}
        for key in re.finditer("([^\(, =]*)=(QUEEN.dna_dict\['[^\[\]]+'\])([^,]*)", source):
            key1, key2, key3 = key.group(1), key.group(2), key.group(3) 
            sourcename = unique_name_dict[key2]
            subsource_dict[sourcename] = (key1, key3) if len(key3) > 0 else (key1, None) 
            connected_keys.append(key1) 
            sourcenames.add(sourcename)
        
        for key in re.finditer("([^\(, =]*)=(QUEEN.queried_features_dict\['[^\[\]]+'\])([^,]*)",source):
            key1, key2, key3 = key.group(1), key.group(2), key.group(3)  
            sourcename = unique_name_dict[key2]
            subsource_dict[sourcename] = (key1, key3) if len(key3) > 0 else (key1, None) 
            connected_keys.append(key1)
            sourcenames.add(sourcename)

        for key in re.finditer("([^\( ,=]*)=([^=]*)[,)]", source):
            key1, key2 = key.group(1), key.group(2)
            if key1 in connected_keys or key1 == "product":
                pass 
            else:
                process_dict[process_id][key1] = key2.replace("'", "") 

        if len(sourcenames) == 0:
            add_key(process_tree, None, process_id)

        for sourcename in list(sourcenames):
            if sourcename in subsource_dict:
                if subsource_dict[sourcename][1] is None:
                    process_dict[process_id][subsource_dict[sourcename][0]] = product_funcname_dict[sourcename] 
                else:
                    process_dict[process_id][subsource_dict[sourcename][0]] = (product_funcname_dict[sourcename][0], subsource_dict[sourcename][1])
            else:
                add_key(process_tree, product_funcname_dict[sourcename][0], process_id)
                process_dict[process_id].input.append(product_funcname_dict[sourcename]) 
        
        if product is None:
            process_dict[process_id]["product"] = productnames[0]
        else:   
            process_dict[process_id]["product"] = product
        
        if matcho is not None:
            n = 0 
            productnames = [] 
            product = matcho.group(1)
            for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", product):
                key = key.group(0) 
                productname = unique_name_dict[key]
                #if process_dict[process_id].funclabel != "entry":
                product_funcname_dict[productname] = (process_id, n)
                productnames.append(productname)  
                n += 1
        else:   
            n = 0 
            productnames = [] 
            product = matchs.group(1)
            for key in re.finditer("QUEEN.queried_features_dict\['[^\[\]]+'\]", product):
                key = key.group(0) 
                productname = unique_name_dict[key]
                productnames.append(productname) 
                product_funcname_dict[productname] = (process_id, n) 
                n += 1

    return process_tree, process_dict        

if __name__ == "__main__":
    adict = {} 
    add_key(adict, None, "A") 
    add_key(adict, None, "B") 
    add_key(adict, None, "C") 
    add_key(adict, "A", "D") 
    add_key(adict, "B", "E") 
    add_key(adict, "D", "E") 
    add_key(adict, "C", "E") 
    add_key(adict, "E", "F") 
    add_key(adict, "E", "G")
    add_key(adict, "E", "H")
    add_key(adict, "H", "I") 

    import pprint 
    pp = pprint.PrettyPrinter(width=10,compact=True)
    pp.pprint(adict)
    
    for target in ["A", "B", "C", "D", "E", "F", "G", "H", "I"]:
        depth = get_depth(adict, target)
        print(depth, target)
