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
    
def export(names, descriptions, histories, o=None, do=False):
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
                    print("process{}={{'name':{}, 'description':{}}}".format(num+1, str(process_name), str(process_description)), file=o) 

            if "QUEEN.queried_feature_dict" in history[1][0:len("QUEEN.queried_feature_dict")]:
                print(history[1], file=o)
            elif (process_name is None and process_description is None) or num == -1:
                print(history[1], file=o)
            else:
                print(history[1][:-1] + ", process_name=process{}['name'], process_description=process{}['description'])".format(num+1, num+1), file=o)

        pre_process_name = process_name 
        pre_process_description = process_description
    return o 

def quine(*dnas, output=None, author=None, project=None, process_description=False, execution=False, _check=False, _return=False, _io=False): 
    if execution == True and output is None:
        output  = tempfile.NamedTemporaryFile(mode="w+", delete=False) 
        outname = output.name 
    else:
        outname = output
    
    if project == None:
        project = dnas[0].project
    
    description_only = process_description

    commands  = []
    histories = [] 
    for dna in dnas:
        for history in dna.history:
            command =  history[1].replace(" ","").replace("–"," ")
            if command in commands:
                pass 
            else:
                histories.append(history) 
                commands.append(command) 
    hindex = 0
    histories.sort() 
    for index, history in enumerate(histories):
        if "QUEEN.dna_dict" in history[1].replace(" ","").replace("–"," "):
            hindex = index
        else:
            pass 
    result     = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]",histories[hindex][1].replace(" ","").replace("–"," "))[0] 
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
        if re.search("process_description=None",history[1].replace(" ","").replace("–"," ")) is None:
            process_description = re.search("process_description='[^=]*'",history[1].replace(" ","").replace("–"," "))
            if process_description is None:
                pd = None                
            else: 
                process_description = history[1].replace(" ","").replace("–"," ")[process_description.start():process_description.end()] 
                pd = process_description.split("=")[1]    
                if pd == "''" or pd == '""':
                    pd = None   
            descriptions.append(pd)
        else:
            pd = None
            process_description = "process_description=None"
            descriptions.append(pd)
        
        if re.search("process_name=None",history[1].replace(" ","").replace("–"," ")) is None:
            process_name = re.search("process_name='[^=]*'",history[1].replace(" ","").replace("–"," "))
            #print(process_name) 
            if process_name is None:
                pn = pre_pn            
            else: 
                process_name = history[1].replace(" ","").replace("–"," ")[process_name.start():process_name.end()] 
                pn = process_name.split("=")[1]    
                if pn == "''" or pn == '""':
                    pn = None    
            names.append(pn)
        else:
            pn = None
            process_name = "process_name=None"
            names.append(pn)
        
        pnflag = 0 
        if process_name is not None and _return == False:
            pnflag = 1
            history[1] = history[1].replace(" ","").replace("–"," ").replace(process_name,"")
            history[1] = history[1].replace(", )",")")
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

        #else:
            #process_id = original_id.split("-")[1]   
            #original_project = original_id.split("-")[0]  
            #process_id = project + "-" + process_id 
            #if original_project != project:
            #original_ids = "[" + ",".join(list(map(lambda x:"'{}'".format(x), history[3].split(",")))) + "]"
            #else: 
            #    _ids = history[3].split(",") 
            #    if len(_ids) == 1:
            #        original_ids = "[]"
            #    else:
            #        original_ids = "[" + ",".join(list(map(lambda x:"'{}'".format(x), _ids[1:]))) + "]"
        
        processid_originalids_dict[process_id] = original_ids
        if process_description is not None and _return == False:
            if pnflag == 1: 
                history[1] = history[1].replace(process_description,"")
                history[1] = history[1].replace(", )",")")          
                history[1] = history[1][:-1] + "process_id='" + process_id + "')"
            else:
                history[1] = history[1].replace(" ","").replace("–"," ").replace(process_description,"")
                history[1] = history[1].replace(", )",")")
                history[1] = history[1][:-1] + ", process_id='" + process_id + "')"
        else:
            history[1] = history[1].replace(" ","").replace("–"," ")
            history[1] = history[1][:-1] + ", process_id='" + process_id + "')"

        new_histories.append(history) 
        pre_pd = pd
    histories = new_histories

    #Remove non-used variables 
    outtext = export(names, descriptions, histories, io.StringIO())
    text    = outtext.getvalue().rstrip()
    var_num_dict = collections.defaultdict(int) 
    for row in text.split("\n"):
        row = row.rstrip()
        matches = re.findall("QUEEN.queried_feature_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_num_dict[match] += 1
        
        matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
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
        matches = re.findall("QUEEN.queried_feature_dict\['[^\[\]]+'\]",row)
        if matches is not None:
            for match in matches:
                var_nums.append(var_num_dict[match]) 
        
        matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
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
    if description_only == False and _return == False: 
        outtext = export(new_names, new_descriptions, new_histories, io.StringIO(), do=description_only)
    outtext = outtext.getvalue().rstrip()
    texts   = outtext.split("\n") 
    
    if _return == True:
        return new_histories

    if description_only == False:
        texts.append("if __name__ == '__main__':") 
        texts.append("    " + result + ".outputgbk()")

    name_dict = {}
    for row in texts:
        match1 = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = ",row)
        match2 = re.search("product='([^=]+)'[,\)]",row) 
        if match1 is not None and match2 is not None:
            name_dict[match1.group(1)] = match2.group(1) 

        match1 = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = ", row)
        if match1 is not None and match2 is not None:
            match3 = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]", match1.group(1))
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
    
    new_rows = []
    for row in texts:
        if dna.__class__._namespaceflag == 1 and execution == False:
            product = re.search("product=([^=]+)[,\)]",row)
            if product is None:
                pass 
            elif product.group(1) == "None":
                pass 
            else:
                matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\] = ",row)
                if len(matches) > 0:
                    for match in matches:
                        row  = row.replace(match, "")
            
                match = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = ",row)
                if match is not None:
                    row  = row.replace(match.group(0), "")
        
        matches = re.findall("QUEEN.queried_features_dict\['[^\[\]]+'\]",row)
        if len(matches) > 0:
            for match in matches:
                name = match
                if name in name_dict:
                    row = row.replace(match, name_dict[name]) 
        
        match = re.search("QUEEN.dna_dict\['[^\[\]]+'\], QUEEN.dna_dict\['[^\[\]]+'\]", row)
        if match is not None:
            for name in name_dict:
                if name in row and match.group(0) in name:
                    row = row.replace(name, name_dict[name])

        matches = re.findall("QUEEN.dna_dict\['[^\[\]]+'\]",row)
        if len(matches) > 0:
            for match in matches:
                name = match
                if name in name_dict:
                    row = row.replace(match, name_dict[name])
        new_rows.append(row) 
   
    #Check quine code is identical with original file.
    identical    = 1
    new_new_rows = [] 
    for row in new_rows: 
        match        = re.search("process_id='([^=]*)'", row)
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
    
    new_rows      = new_new_rows
    project_names = [] 
    if _check == True and sys.argv[0] != "":
        check_flag = 0 
        with open(sys.argv[0]) as f:
            lnum = 0
            for line in f:
                line = line.rstrip()
                if check_flag == 1:
                    match = re.search("process_id='([^=]*)-([^=]*)'", line)
                    query = new_rows[lnum]
                    if match is not None:
                        project_names.append(match.group(1))
                        query = re.sub(", _sourcefile=.*\)", ")".format(project), query)
                        line  = re.sub(", _sourcefile=.*\)", ")".format(project), line)
                        if line == query:
                            pass 
                        else: 
                            identical = 0 
                    else:
                        if line == query:
                            pass 
                        else:
                            identical = 0 
                    lnum += 1 
                
                if line == "":
                    check_flag = 1 
        
        if identical == 1: 
            return True
                
    if _check == False or identical == 0:            
        if description_only == False:
            now = datetime.datetime.now()
            #print("#" * 80, file=o)
            #print("#This source code was auto-generated by 'quine' funtion of QUEEN 1.0.0.", file=o)
            #print("#Project Name    :{}".format(", ".join(list(map(lambda x: x.project if x.project not in project_names else x.project + "_modified", dnas)))), file=o)
            #print("#File Name       :{}".format(outname), file=o) 
            #print("#Creation Date   :{}".format("{0:%Y-%m-%d}".format(now)), file=o)      
            #if author is not None:
            #    print("Copyright © {%Y} ".format(now) + "{} All rights reserved.".format(author), file=o)
            #print("#" * 80, file=o)
            print("project='{}'".format(project), file=o)
            print("import sys", file=o)  
            print("sys.path.append(\"{}".format("/".join(__file__.split("/")[:-2])  + "\")"), file=o)
            print("from QUEEN.queen import *", file=o) 
            print("from QUEEN import cutsite as cs", file=o) 
            for cutsite in list(cs.new_cutsites):
                print("cs.lib[{}] = '{}'".format(cutsite[0], cutsite[1]), file=o) 
            if dna.__class__._namespaceflag == 1 and execution == False:
                print("set_namespace(globals())", file=o)
            print("", file=o) 
            for row in new_rows:
                match = re.search("process_id='([^=]*)'", row)
                if match is not None:
                    if "-" in match.group(1):
                        if outname is not None:
                            row = row[:-1] + ", _sourcefile='{}')".format(outname.split("/")[-1].rstrip(".py")) 
                        else:
                            row = row[:-1] + ", _sourcefile='{}')".format(dnas[0].project + "_construction")  
                    else:
                        if dnas[0].project in project_names:
                            row = re.sub("process_id='([^=]*)'", "process_id='{}_modified-\\1'".format(project), row)
                            if outname is not None:
                                row = row[:-1] + ", _sourcefile='{}')".format(outname.split("/")[-1].rstrip(".py"))
                        else:
                            row = re.sub("process_id='([^=]*)'", "process_id='{}-\\1'".format(project), row)
                            if outname is not None:
                                row = row[:-1] + ", _sourcefile='{}')".format(outname.split("/")[-1].rstrip(".py"))
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
        exec("varnames = dir({})".format(fname), locals(), vardict) 
        exec("queen_objects = {}.QUEEN._products".format(fname), locals(), vardict) 
        
        if flag == 1:
            dnas[0].__class__._namespaceflag = 1
            dnas[0].__class__._namespace = _globalspace
        
        if type(output) is tempfile._TemporaryFileWrapper: 
            os.remove(outname + ".py") 
        
        dnas[0].__class__._source = None 
        for name, obj in vardict["queen_objects"].items():
            for key in obj._history_feature.qualifiers:
                if "building_history" in key[0:18]:
                    obj._history_feature.qualifiers[key][1] = obj._history_feature.qualifiers[key][1].replace(fname, dnas[0].project + " construction")
        keys = list(vardict["queen_objects"].keys())
        
        if dnas[0] == vardict["queen_objects"][keys[-1]]:
            if _io == True:
                print("QUEEN object reconstructed from the quine code and the original QUEEN object are identical.".format(dnas[0].project))
            dnas[0]._productids = keys[:-1]   
            return True
        else:
            if _io == True:
                raise ValueError("The {} QUEEN object could not be reconstructed using its quine code. There may be bugs in QUEEN implementation. It would be helpful if you could tell us the details about your code on the Github issue (https://github.com/yachielab/QUEEN/issues).".format(dnas[0].project)) 
            return False
