import os
import io
import sys
import random
import datetime
import tempfile
import collections
import regex as re
sys.path.append("/".join(__file__.split("/")[:-1]))
import cutsite as cs

def check_processids(dna):
    for process_id in dna._processids:
        chars = dna.__class__.processes[process_id]
        
        num = 0
        for char in chars:
            num += ord(char) 
        
        random.seed(num) 
        new_process_id = ''.join(random.choices(string.ascii_letters + string.digits, k=10)) 
        while process_id in dna.__class__._processes:
            new_process_id = ''.join(random.choices(string.ascii_letters + string.digits, k=10)) 

        if new_process_id != process_id:
            return new_process_id
        else:
            False
    
def export(names, descriptions, histories, o=None, do=False, qexp_only=False):
    num = -1
    pre_process_name = None
    pre_process_description = None
    for h, (process_name, process_description, history) in enumerate(zip(names, descriptions, histories)):
        if do == True:
            if (process_description is None and process_name is None) or (str(process_name) == str(pre_process_name) and str(process_description) == str(pre_process_description)):
                pass
            else:
                print("{}:{}".format(process_name, process_description, file=o))
        
        else:
            #if (process_description is None and process_name is None) or (str(process_name) == str(pre_process_name) and str(process_description) == str(pre_process_description)):
            if str(process_name) == str(pre_process_name) and str(process_description) == str(pre_process_description):
                pass
            else:
                if h > 0:
                    print("", file=o)
                
                if process_name != None and process_description != None:
                    num += 1
                    if qexp_only == False:
                        print("process{}={{'name':{}, 'description':{}}}".format(num+1, str(process_name), str(process_description)), file=o) 

            if "QUEEN.queried_feature_dict" in history[1][0:len("QUEEN.queried_feature_dict")]:
                print(history[1], file=o)
            elif (process_name is None and process_description is None) or num == -1:
                print(history[1], file=o)
            else:
                if qexp_only == False:
                    print(history[1][:-1] + ", process_name=process{}['name'], process_description=process{}['description'])".format(num+1, num+1), file=o)
                else:
                    print(history[1][:-1] + ", process_name={}, process_description={})".format(process_name, process_description), file=o)
        pre_process_name = process_name 
        pre_process_description = process_description
    return o 

def quine(*dnas, output=None, author=None, project=None, process_description=False, qexperiment_only=True, execution=False, _return_histories=False, _return_script=False, _io=False): 
    """Generate "quine code" of `QUEEN_object` that produces the same `QUEEN_object`. A quine code can be executed as a Python script.
    
    Parameters
    ----------
    *dnas: QUEEN.qobj.QUEEN object
    output: str ,default: STDOUT   
        Output file name.
    process_description: bool, default: False)
        If True, this will output the process_descriptions registered to quinable operations along with the process flows. 
        The output can be used for the "Material and methods" of the `QUEEN_object` construction process.
    qexperiment_only: bool, default: True) 
        If True, this will output only the qexperiment commands and not output the internal quine commmands used in the qexperiment commands. 
    execution: bool, default: False  
        If True, this will reconstruct the `QUEEN_object` by generating and executing its quine code and confirm if the reconstructed 
        `QUEEN_object` is identical to the original one. If `execution` is `True` and `output` is `None`, the quine code will be output 
        into a temporary file instead of `STDOUT`; the temporary file will be removed after the operation. 
        The execution won't happen if `process_description` is `True`.

    Returns
    -------
    if `execution` is `False`, `None`.
    If `execution` is `True`, `True` if the reconstructed `QUEEN_object` is identical to the original one. Otherwise, `False`.

    """
    def extract_qexd(rows): 
        extracted_rows = [] 
        for row in rows:
            if "qexd =" not in row and "qexd=" not in row:
                extracted_rows.append(row) 
            else:
                pattern1 = r"qexd='(.*?)'"
                pattern2 = r"qexd = '(.*?)'"
                match1 = re.search(pattern1, row) 
                match2 = re.search(pattern2, row)
                if match1 is not None:
                    txt = match1.group(1)
                elif match2 is not None:
                    txt = match2.group(1)
                else:
                    txt = None
                if txt is None or "qexd=True" in txt:
                    pass 
                else:
                    outdna = row.split("=")[0] 
                    if "product=" in row:
                        index_start = row.find("product=") 
                    elif "process_name=" in row:
                        index_start = row.find("process_name=") 
                    elif "process_description" in row:
                        index_start = row.find("process_description=") 
                    elif "process_id=" in row:
                        index_start = row.find("process_id=")
                    txt = outdna.rstrip() + " = " +  txt[:-1].replace('"',"'") + ", " + row[index_start:]
                    extracted_rows.append(txt)
        return extracted_rows

    if execution == True and output is None:
        output  = tempfile.NamedTemporaryFile(mode="w+", delete=False) 
        outname = output.name 
    else:
        outname = output
    
    if project == None:
        project = dnas[0].project
    
    description_only = process_description

    commands     = []
    history_dict = collections.defaultdict(dict) 
    for dna in dnas:
        for key in dna.history:
            if "_script" in key:
                command = dna.history[key] 
                if command in commands:
                    pass 
                else:
                    commands.append(command) 
            key1 = int(key.split("_")[0]) 
            key2 = key.split("_")[1]
            history_dict[key1][key2] = dna.history[key] 
    
    histories = [] 
    for key in history_dict:
        histories.append([key, history_dict[key]["script"], history_dict[key]["args"], history_dict[key]["id"]])

    hindex = 0
    histories.sort() 
    for index, history in enumerate(histories):
        history1 = history[1].replace(" ","").replace("–"," ") if ",–" in history[1] else history[1]
        if "QUEEN.dna_dict" in history1:
            hindex = index
        else:
            pass 

    history1 = histories[hindex][1].replace(" ","").replace("–"," ") if ",–" in histories[hindex][1] else histories[hindex][1]
    result   = re.findall(r"QUEEN.dna_dict\['[^\[\]]+'\]", history1.replace(" ","").replace("–"," "))[0] 
    _unique_id = result.split("['")[1][:-2]
    
    pre_pd = None
    pre_pn = None
    edited_flag   = 0 
    names         = [] 
    descriptions  = []
    new_histories = [] 
    processid_originalids_dict = {}  
    for history in histories:
        history = list(history)
        history1 = history[1].replace(" ","").replace("–"," ") if ",–" in history[1] else history[1]
        if re.search(r"process_description=None",history1) is None:
            process_description = re.search(r"process_description='[^']*'",history1)
            if process_description is None:
                pd = None                
            else: 
                process_description = history1[process_description.start():process_description.end()] 
                pd = process_description.split("=")[1]    
                if pd == "''" or pd == '""':
                    pd = None   
            descriptions.append(pd)
        else:
            pd = None
            process_description = "process_description=None"
            descriptions.append(pd)
        
        if re.search(r"process_name=None",history1) is None:
            process_name = re.search(r"process_name='[^']*'",history1)
            if process_name is None:
                pn = pre_pn            
            else: 
                process_name = history1[process_name.start():process_name.end()] 
                pn = process_name.split("=")[1]    
                if pn == "''" or pn == '""':
                    pn = None    
            names.append(pn)
        else:
            pn = None
            process_name = "process_name=None"
            names.append(pn)
        
         

        pnflag = 0 
        if process_name is not None and _return_histories == False:
            pnflag = 1
            history1 = history1.replace(process_name+", ","").replace(process_name+",","").replace(process_name,"")
            history1 = history1.replace(", )",")").replace(",)",")")
        else:
            pass 
        original_id = history[3].split(",")[0]
        if "-" not in original_id:
            process_id = original_id 
            if len(history[3].split(",")) > 1:
                original_ids = "[" + ",".join(list(map(lambda x:"'{}'".format(x), history[3].split(",")[1:]))) + "]"
            else:
                original_ids = []
        else:
            process_id = original_id
            if len(history[3].split(",")) > 1:
                original_ids = "[" + ",".join(list(map(lambda x:"'{}'".format(x), history[3].split(",")[1:]))) + "]"
            else:
                original_ids = []
        
        processid_originalids_dict[process_id] = original_ids
        if process_description is not None and _return_histories == False:
            history1 = history1.replace(process_description+", ","").replace(process_description+",","").replace(process_description,"")
            history1 = history1.replace(", )",")").replace(",)",")")         
            history[1] = history1[:-1] + ", process_id='" + process_id + "')"
        else:
            history[1] = history1[:-1] + ", process_id='" + process_id + "')"
        new_histories.append(history) 
        pre_pd = pd
    histories = new_histories

    #Remove non-used variables 
    outtext = export(names, descriptions, histories, io.StringIO(), qexp_only=qexperiment_only)
    text    = outtext.getvalue().rstrip()
    var_num_dict = collections.defaultdict(int) 
    for row in text.split("\n"):
        row = row.rstrip()
        matches = re.findall(r"QUEEN.queried_feature_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_num_dict[match] += 1
        
        matches = re.findall(r"QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_num_dict[match] += 1
    
    new_names        = [] 
    new_descriptions = []
    new_histories    = []
    
    source_descriptions_dict = collections.defaultdict(list) 
    texts = [row for row in text.split("\n") if row[0:len("QUEEN.")] == "QUEEN."] 
    for row, name, description, history in zip(texts, names, descriptions, histories):
        row      = row.rstrip() 
        var_nums = [] 
        matches = re.findall(r"QUEEN.queried_feature_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_nums.append(var_num_dict[match]) 
        
        matches = re.findall(r"QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_nums.append(var_num_dict[match])
        
        if len(var_nums) > 0:
            if max(var_nums) == 1:
                pass 
            else:
                new_names.append(name) 
                new_descriptions.append(description) 
                new_histories.append(history)  
                info = history[2]
                if len(info) > 1:
                    info = info.split("; ")
                    info_dict = dict([item.split(": ") for item in info])
                    if "_source" in info_dict:
                        if description is not None:
                            source_descriptions_dict[info_dict["_source"]].append(description) 
                    else:
                        if description is not None:
                            source_descriptions_dict[dnas[0].project].append(description) 
                else:
                    if description is not None:
                        source_descriptions_dict[dnas[0].project].append(description) 
        else:
            new_names.append(name) 
            new_descriptions.append(description) 
            new_histories.append(history) 
            
            info = history[2]
            if len(info) > 1:
                info = info.split("; ")
                info_dict = dict([item.split(": ") for item in info])
                if "_source" in info_dict:
                    if description is not None:
                        source_descriptions_dict[info_dict["_source"]].append(description) 
                else: 
                    if description is not None:
                        source_descriptions_dict[dnas[0].project].append(description) 
    
    for key in source_descriptions_dict:
        source_descriptions_dict[key] = list(sorted(set(source_descriptions_dict[key]), key=source_descriptions_dict[key].index)) 
    
    if type(output) is str:
        o = open(output, "w") 
    elif output is None:
        o = None
    elif type(output) is io.TextIOWrapper:
        o = output
    elif type(output) is tempfile._TemporaryFileWrapper:
        o = output
    
    now = datetime.datetime.now()
    pre_process_description = "''"
    if description_only == False and _return_histories == False: 
        outtext = export(new_names, new_descriptions, new_histories, io.StringIO(), do=description_only, qexp_only=qexperiment_only)
    outtext = outtext.getvalue().rstrip()
    texts   = outtext.split("\n") 
    
    if _return_histories == True:
        return new_histories

    if description_only == False:
        texts.append("if __name__ == '__main__':") 
        texts.append("    " + result + ".outputgbk()")

    name_dict = {}
    for row in texts:
        match1 = re.search(r"(QUEEN.queried_features_dict\['[^\[\]]+'\]) = ",row)
        match2 = re.search(r"product='([^=]+)'[,\)]",row) 
        if match1 is not None and match2 is not None:
            name_dict[match1.group(1)] = match2.group(1) 

        match1 = re.search(r"(.*QUEEN.dna_dict\['[^\[\]]+'\]) = ", row)
        if match1 is not None and match2 is not None:
            match3 = re.findall(r"QUEEN.dna_dict\['[^\[\]]+'\]", match1.group(1))
            name_dict[match1.group(1)] = match2.group(1) 
            if len(match3) == 1:
                pass 
            else:
                if "," in match2.group(1):
                    for match, name in zip(match3, match2.group(1).split(",")):
                        name_dict[match] = name
                else:
                    for m, match in enumerate(match3):
                        name_dict[match] = match2.group(1) + "[{}]".format(m)
    
        
    new_rows = texts 
    #Check quine code is identical with original file.
    identical    = 1
    new_new_rows = [] 
    for row in new_rows: 
        match        = re.search(r"process_id='([^=]*)'", row)
        if match is not None:
            source1      = match.group(0) + ", "
            source2      = match.group(0) + ")"
            process_id   = match.group(1) 
            original_ids = processid_originalids_dict[process_id] 
            if source1 in row:
                new_new_rows.append(row.replace(source1, "")[:-1] + ", process_id='{}', original_ids={})".format(process_id, original_ids))
            elif source2 in row:
                new_new_rows.append(row.replace(source2, "") + "process_id='{}', original_ids={})".format(process_id, original_ids))
        else:
            new_new_rows.append(row) 
    
    new_rows = new_new_rows
    if qexperiment_only == True:
        new_rows = extract_qexd(new_rows) 

    project_names = []
    if description_only == False:
        now = datetime.datetime.now()
        if _return_script == False:
            print("project='{}'".format(project), file=o)
            print("import sys", file=o)  
            print("sys.path = [\"{}] + sys.path".format("/".join(__file__.split("/")[:-2])  + "\""), file=o)
            print("from QUEEN.queen import *", file=o) 
            print("from QUEEN import cutsite as cs", file=o) 
            for cutsite in list(cs.new_cutsites):
                print("cs.lib[{}] = '{}'".format(cutsite[0], cutsite[1]), file=o) 
            if dna.__class__._namespaceflag == 1 and execution == False:
                print("set_namespace(globals())", file=o)
            print("", file=o) 
        
        scripts = [] 
        for row in new_rows:
            match = re.search(r"process_id='([^=]*)'", row)
            if match is not None:
                if "-" in match.group(1):
                    if outname is not None:
                        row = row[:-1] + ", _sourcefile='{}')".format(outname.split("/")[-1].rstrip(".py")) 
                    else:
                        process_id  = match.group(1).split("-")[-1] 
                        productname = "-".join(match.group(1).split("-")[:-1]) 
                        process_id.split("-")[-1] 
                        row = row[:-1] + ", _sourcefile='{}')".format(productname + "_construction")  
                else:
                    if dnas[0].project in project_names:
                        row = re.sub(r"process_id='([^=]*)'", "process_id='{}_modified-\\1'".format(project), row)
                        if outname is not None:
                            row = row[:-1] + ", _sourcefile='{}')".format(outname.split("/")[-1].rstrip(".py"))
                    else:
                        row = re.sub(r"process_id='([^=]*)'", "process_id='{}-\\1'".format(project), row)
                        if outname is not None:
                            row = row[:-1] + ", _sourcefile='{}')".format(outname.split("/")[-1].rstrip(".py"))
            scripts.append(row) 
            if _return_script == False:
                print(row, file=o)
    
    else:
        if len(source_descriptions_dict) == 1:
            for key in source_descriptions_dict:
                for value in source_descriptions_dict[key]:
                    print(value[1:-1], file=o) 
        else:
            for key in source_descriptions_dict:
                print("{}{}".format(key.rstrip("_construction").rstrip(), " construction" if "construction" not in key.rstrip("_construction") else "")) 
                for value in source_descriptions_dict[key]:
                    print("    " + value[1:-1], file=o) 

    if output is not None:
        o.close()
    
    if _return_script == True:
        return scripts 

    if execution == True:
        intermediate_products = {} 
        sys.path.append("/".join(outname.split("/")[:-1]) if len(outname.split("/")) > 1 else ".")
         
        if outname[-3:] != ".py":
            os.rename(outname, outname + ".py") 
        
        flag = 0 
        if dnas[0].__class__._namespaceflag == 1:
            _globalspace = dnas[0].__class__._namespace
            dnas[0].__class__._namespace     = {} 
            dnas[0].__class__._namespaceflag = 0
            flag = 1

        fname = outname.split("/")[-1]
        if fname[-3:] == ".py":
            fname = fname[:-3] 
        vardict = {}
        dnas[0].__class__._source = fname
        
        exec("import {}".format(fname), locals(), vardict)
        exec("queen_objects = {}.QUEEN._products".format(fname), locals(), vardict) 
        if flag == 1:
            dnas[0].__class__._namespaceflag = 1
            dnas[0].__class__._namespace = _globalspace
        
        if type(output) is tempfile._TemporaryFileWrapper: 
            os.remove(outname + ".py") 
        
        dnas[0].__class__._source = None 
        keys = list(vardict["queen_objects"].keys())
        if dnas[0] == vardict["queen_objects"][keys[-1]]:
            if _io == True:
                print("QUEEN object reconstructed from the quine code and the original QUEEN object are identical.".format(dnas[0].project))
            dnas[0]._productids = keys[:-1]   
            return True, vardict["queen_objects"] 
        else:
            if _io == True:
                raise ValueError("The {} QUEEN object could not be reconstructed using its quine code. There may be bugs in QUEEN implementation. It would be helpful if you could tell us the details about your code on the Github issue (https://github.com/yachielab/QUEEN/issues).".format(dnas[0].project)) 
            return False

def printprotocol(dna, execution=False, output=None):
    """
    This function retrieve and return the history of only the qexperiment functions.
    If there are no qexperiment functions in the QUEEN script. It will return an empty str object.
    
    Parameters
    ----------
    dna: QUEEN.qobj.QUEEN object
    
    output: str ,default: STDOUT   
        Output file name.    
    """
    print("#{} Construction".format(dna.project), file=output)
    pattern_dict = {
                    "pcr":       r"pcr\((.*)\)",
                    "digestion": r"digestion\((.*)\)",
                    "ligation":  r"ligation\((.*)\)",
                    "hba":       r"homology_based_assembly\((.*)\)",
                    "anneal":    r"annealing\((.*)\)",
                    "gga":       r"golden_gate_assembly\((.*)\)",
                    "gateway":   r"gateway_reaction\((.*)\)",
                    "topo":      r"topo_cloning\((.*)\)"
                    }
    
    if execution == True: 
        qobjects = quine(dna, execution=True)[1]  
    
    for row in quine(dna, _return_script=True):
        row = row.replace(" = ","=")
        pro_names = [] 
        for match in re.finditer(r"QUEEN.dna_dict\['([^\[\]]+)'\]", row): 
            qpro_name = match.group(0)
            pro_name  = match.group(1)
            pro_names.append(pro_name) 
       
        if (match := re.search(pattern_dict["pcr"], row)):
            product  = pro_names[0] 
            template = pro_names[1]
            row    = "=".join(row.split("=")[1:]) 
            row_sp = row.split(",")
            if "QUEEN" in row_sp[1]:
                fw_name = re.search(r"QUEEN.dna_dict\['([^\[\]]+)'\]", row_sp[1]).group(1) 
            else:
                fw_name = row_sp[1]

            if "QUEEN" in row_sp[2]:
                rv_name = re.search(r"QUEEN.dna_dict\['([^\[\]]+)'\]", row_sp[2]).group(1) 
            else:
                rv_name = row_sp[2]
            
            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">PCR{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Template: {}, {} bp".format(template, len(qobjects[template].seq)), file=output) 
                print("- Forward Primer: {}, {}, {} bp".format(fw_name, qobjects[fw_name].seq, len(qobjects[fw_name].seq)), file=output)
                print("- Reverse Primer: {}, {}, {} bp".format(rv_name, qobjects[rv_name].seq, len(qobjects[rv_name].seq)), file=output)   
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}, {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}, {}".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            
            else:
                print("Parameters:", file=output) 
                print("- Template: {}".format(template), file=output)  
                print("- Forward Primer: {}".format(fw_name), file=output)
                print("- Reverse Primer: {}".format(rv_name), file=output)    
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output) 

        elif (match := re.search(pattern_dict["digestion"], row)):
            product  = pro_names[0] 
            sample   = pro_names[1]
            row      = "=".join(row.split("=")[1:]) 
            row_sp   = row[:row.find(", selection=")].split(",")
            cutsites = ",".join(list([cutsite.replace(" ","").replace("'","") for cutsite in row_sp[1:]])) 

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">Digestion{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Sample: {}, {} bp".format(sample, len(qobjects[sample].seq)), file=output) 
                print("- Restriction enzyme(s): {}".format(cutsites), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {}".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            
            else:
                print("Parameters:", file=output) 
                print("- Sample: {}".format(sample), file=output) 
                print("- Restriction enzyme(s): {}".format(cutsites), file=output)   
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output)
        
        elif (match := re.search(pattern_dict["ligation"], row)):
            product  = pro_names[0] 
            sample   = ", ".join(pro_names[1:])

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">Ligation{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Sample(s): {}; {}".format(sample, ", ".join([len(qobjects[asample].seq) for asample in sample])), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("follow_order", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {}".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            else:
                print("Parameters:", file=output) 
                print("- Sample(s): {}".format(sample), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("follow_order", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {}".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
 
        elif (match := re.search(pattern_dict["hba"], row)):
            product  = pro_names[0] 
            sample   = ", ".join(pro_names[1:])

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            mode = re.search(r"mode='([^']*)'",match.group(1)).group(1)
            print(">Homology based Assembly{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Sample(s): {}; {}".format(sample, ", ".join([str(len(qobjects[asample].seq)) + " bp" for asample in sample])), file=output) 
                print("- Assembly method: {}".format(mode), file=output)
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("mode", "follow_order", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {}".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            else:
                print("Parameters:", file=output) 
                print("- Sample(s): {}".format(sample), file=output) 
                print("- Assembly method: {}".format(mode), file=output)
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("mode", "follow_order", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output) 

        elif (match := re.search(pattern_dict["anneal"], row)):
            product  = pro_names[0] 
            ssdna1   = pro_names[1] 
            ssdna2   = pro_names[2]

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">Annealing{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Top strand DNA: {}, {}".format(ssdna1, qobjects[ssdna1].seq), file=output) 
                print("- Bottom strand DNA: {}, {}".format(ssdna2, qobjects[ssdna2].seq), file=output)
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("homology_length", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {} bp".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            else:
                print("Parameters:", file=output) 
                print("- Top strand DNA: {}".format(ssdna1), file=output) 
                print("- Bottom strand DNA: {}".format(ssdna2), file=output)
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("homology_length", "follow_order", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output) 

        elif (match := re.search(pattern_dict["gga"], row)):
            product     = pro_names[0] 
            destination = pro_names[1] 
            entry       = pro_names[2:]

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">Goden Gate Assembly{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            
            cutsite = re.search(r", cutsite='([^']*)'", match.group(1)).group(1) 
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Destination sample: {}, {} bp".format(destination, len(qobjects[destination].seq)), file=output)
                print("- Entry sample(s): {}; {}".format(", ".join(entry), ", ".join([str(len(qobjects[asample].seq)) + "bp" for asample in entry])), file=output) 
                print("- Restriction enzyme: {}".format(cutsite), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("cutsite", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {} bp".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            else:
                print("Parameters:", file=output) 
                print("- Destination sample: {}".format(destination), file=output)
                print("- Entry Sample(s): {}".format(", ".join(entry)), file=output)
                print("- Restriction enzyme: {}".format(cutsite), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("cutsite", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output)
        
        elif (match := re.search(pattern_dict["gateway"], row)):
            product     = pro_names[0] 
            destination = pro_names[1] 
            entry       = pro_names[2]

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">Gateway Cloning{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            
            mode = re.search(r"mode='([^']*)'",match.group(1)).group(1)
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Destination sample: {}, {} bp".format(destination, len(qobjects[destination].seq)), file=output)
                print("- Entry sample: {}; {} bp".format(entry, len(qobjects[entry].seq)), file=output) 
                print("- BP/LR: {}".format(mode), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("mode", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {} bp".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            else:
                print("Parameters:", file=output) 
                print("- Destination sample: {}".format(destination), file=output)
                print("- Entry sample(s): {}".format(entry), file=output)
                print("- BP/LR: {}".format(mode), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("mode", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output)
        
        elif (match := re.search(pattern_dict["topo"], row)):
            product     = pro_names[0] 
            destination = pro_names[1] 
            entry       = pro_names[2]

            if "process_name=" in row:
                pn = ": " + re.search(r"process_name='([^']*)'", match.group(1)).group(1) 
            else:
                pn = "" 
            
            if "process_description=" in row:
                pdmatch = re.search(r"process_description='([^']*)'", match.group(1))
                pd = pdmatch.group(1) if pdmatch is not None else None 
            else:
                pd = None 
            
            print(">TOPO Cloning{}".format(pn), file=output)
            if pd is not None:
                print("Description:\n{}".format(pd), file=output) 
            else:
                pass
            
            mode = re.search(r"mode='([^']*)'",match.group(1)).group(1)
            if execution  == True:
                print("Parameters:", file=output) 
                print("- Destination sample: {}, {} bp".format(destination, len(qobjects[destination].seq)), file=output)
                print("- Entry sample: {}; {} bp".format(entry, len(qobjects[entry].seq)), file=output) 
                print("- Cloning method: {}".format(mode), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("mode", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}: {} bp".format(product, len(qobjects[product].seq)), file=output)  
                print("", file=output) 
            else:
                print("Parameters:", file=output) 
                print("- Destination sample: {}".format(destination), file=output)
                print("- Entry Sample(s): {}".format(entry), file=output)
                print("- Cloning method: {}".format(mode), file=output) 
                for arg in re.finditer(r", ([^'=]*)='([^']*)'|, ([^'=]*)=(None)", match.group(1)):
                    if arg.group(1) is not None:
                        key   = arg.group(1)
                        value = arg.group(2)
                    else:
                        key   = arg.group(3)
                        value = arg.group(4)
                    if key not in ("mode", "product", "process_id", "process_name", "process_description") and "original_ids" not in key and "_sourcefile" not in key:
                        print("- {}: {}".format(key.capitalize(), value), file=output) 
                print("Output:\n{}".format(product), file=output)  
                print("", file=output)
        
        #print(row) 

