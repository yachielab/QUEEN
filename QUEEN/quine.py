import os
import io
import sys
import random
import datetime
import tempfile
import subprocess
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


def _execute_quine_script_isolated(script_path, result_key, cwd=None, python_executable=None):
    cwd = os.getcwd() if cwd is None else cwd
    python_executable = sys.executable if python_executable is None else python_executable
    module_name = os.path.splitext(os.path.basename(script_path))[0]
    fd, gbk_path = tempfile.mkstemp(suffix=".gbk")
    os.close(fd)
    runner = "\n".join(
        [
            "import importlib.util",
            "import sys",
            "script_path, module_name, result_key, output_gbk = sys.argv[1:5]",
            "spec = importlib.util.spec_from_file_location(module_name, script_path)",
            "if spec is None or spec.loader is None:",
            "    raise RuntimeError(f'failed to load {script_path}')",
            "module = importlib.util.module_from_spec(spec)",
            "spec.loader.exec_module(module)",
            "module.QUEEN.dna_dict[result_key].outputgbk(output_gbk)",
        ]
    )
    proc = subprocess.run(
        [python_executable, "-c", runner, script_path, module_name, result_key, gbk_path],
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    return proc, gbk_path

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
    def _extract_qexd_text(row):
        if all(token not in row for token in ("qexd =", "qexd=", "qexparam =", "qexparam=")):
            return None
        pattern1 = r"qexd='(.*?)'"
        pattern2 = r"qexd = '(.*?)'"
        pattern3 = r"qexparam='(.*?)'"
        pattern4 = r"qexparam = '(.*?)'"
        match1 = re.search(pattern1, row) 
        match2 = re.search(pattern2, row)
        match3 = re.search(pattern3, row)
        match4 = re.search(pattern4, row)
        if match1 is not None:
            txt = match1.group(1)
        elif match2 is not None:
            txt = match2.group(1)
        elif match3 is not None:
            txt = match3.group(1)
        elif match4 is not None:
            txt = match4.group(1)
        else:
            txt = None
        if txt is None or "qexd=True" in txt:
            return None
        return txt.strip()

    def _extract_qexd_row(row):
        txt = _extract_qexd_text(row)
        if txt is None:
            return None

        outdna = row.split("=")[0].strip()
        if "product=" in row:
            index_start = row.find("product=") 
        elif "process_name=" in row:
            index_start = row.find("process_name=") 
        elif "process_description" in row:
            index_start = row.find("process_description=") 
        elif "process_id=" in row:
            index_start = row.find("process_id=")
        else:
            index_start = len(row) - 1

        extracted = outdna.rstrip() + " = " + txt[:-1].replace('"',"'") + ", " + row[index_start:]
        extracted = extracted.replace("follow_order='True'", "follow_order=True")
        extracted = extracted.replace("follow_order='False'", "follow_order=False")
        extracted = extracted.replace('follow_order=\"True\"', "follow_order=True")
        extracted = extracted.replace('follow_order=\"False\"', "follow_order=False")
        return extracted

    def _row_has_qex_metadata(row):
        return any(token in row for token in ("qexd =", "qexd=", "qexparam =", "qexparam="))

    def _is_qexperiment_row(row):
        return re.search(r"=\s*(pcr|digestion|ligation|homology_based_assembly|annealing|gateway_reaction|golden_gate_assembly|goldengate_assembly|topo_cloning|intra_site_specific_recombination|homologous_recombination)\(", row) is not None

    def _is_qexperiment_helper_row(row):
        return re.search(r"QUEEN\.dna_dict\['[^']+'\]\s*,?\s*=\s*(modifyends|cropdna|flipdna|joindna|cutdna)\(", row.strip()) is not None

    def _is_seed_row(row):
        stripped = row.strip()
        return stripped.startswith("QUEEN.dna_dict[") and " = QUEEN(" in stripped

    def _strip_qex_metadata(row):
        patterns = (
            r",\s*qexd\s*=\s*'[^']*'",
            r',\s*qexd\s*=\s*"[^"]*"',
            r",\s*qexd\s*=\s*[A-Za-z_][A-Za-z0-9_]*",
            r",\s*qexparam\s*=\s*'[^']*'",
            r',\s*qexparam\s*=\s*"[^"]*"',
            r",\s*qexparam\s*=\s*[A-Za-z_][A-Za-z0-9_]*",
        )
        for pattern in patterns:
            row = re.sub(pattern, "", row)
        row = re.sub(r",\s*,", ", ", row)
        row = re.sub(r"\(\s*,", "(", row)
        row = re.sub(r",\s*\)", ")", row)
        return row

    def _parse_args_info(args_text):
        info = {}
        if args_text is None or len(args_text) == 0:
            return info
        for item in args_text.split("; "):
            if ": " not in item:
                continue
            key, value = item.split(": ", 1)
            info[key] = value
        return info

    def _literal_arg(value):
        if value is None:
            return None
        value = value.strip()
        if len(value) == 0:
            return None
        if (value.startswith("'") and value.endswith("'")) or (value.startswith('"') and value.endswith('"')):
            return value
        if re.fullmatch(r"-?\d+", value):
            return value
        return repr(value)

    def _rewrite_lower_row_from_args(row, args_text, live_products):
        info = _parse_args_info(args_text)
        if len(info) == 0:
            return row

        def _replace_kwarg(text, key, value, next_keys):
            if value is None:
                return text
            next_pattern = "|".join(re.escape(next_key) for next_key in next_keys)
            pattern = r"{}\s*=.*?(?=,\s*(?:{})\s*=)".format(re.escape(key), next_pattern)
            return re.sub(pattern, "{}={}".format(key, value), text, count=1)

        if " = modifyends(" in row:
            left = _literal_arg(info.get("left"))
            right = _literal_arg(info.get("right"))
            row = _replace_kwarg(row, "left", left, ["right", "process_id", "original_ids", "product"])
            row = _replace_kwarg(row, "right", right, ["process_id", "original_ids", "product"])
            row = _normalize_identity_modifyends_row(row, live_products)
            return row

        if " = cropdna(" in row:
            start = _literal_arg(info.get("start"))
            end = _literal_arg(info.get("end"))
            row = _replace_kwarg(row, "start", start, ["end", "process_id", "original_ids", "product"])
            row = _replace_kwarg(row, "end", end, ["process_id", "original_ids", "product"])
            return row

        return row

    def _normalize_record_path_row(row):
        def repl(match):
            addgene_id = match.group(1)
            candidate = os.path.join(os.getcwd(), "gbks", f"new_{addgene_id}.gbk")
            if os.path.exists(candidate):
                return "record='gbks/new_{}.gbk'".format(addgene_id)
            return match.group(0)

        row = re.sub(r"record='[^']*addgene_(\d+)_addgene\.gbk'", repl, row)
        row = re.sub(r'record=\"[^\"]*addgene_(\d+)_addgene\.gbk\"', lambda m: repl(m).replace("'", '"'), row)
        return row

    def _normalize_none_join_row(row):
        match = re.search(r"^(.*?=\s*)joindna\(\*\[(.*)\](,\s*topology=.*)\)$", row)
        if match is None or "QUEEN.dna_dict['None']" not in row:
            return row

        lhs = match.group(1)
        args_blob = match.group(2)
        suffix = match.group(3)
        args = [part.strip() for part in args_blob.split(",") if part.strip()]
        args = [part for part in args if part != "QUEEN.dna_dict['None']"]
        if len(args) == 1:
            return lhs + args[0]
        if len(args) >= 2:
            return lhs + "joindna(*[" + ", ".join(args) + "]" + suffix + ")"
        return row

    def _normalize_identity_modifyends_row(row, live_products):
        match = re.search(
            r"^(QUEEN\.dna_dict\['[^']+'\]\s*=\s*)modifyends\(([^,]+),\s*left=(['\"])\*/\*\3,\s*right=(['\"])\*/\*\4(?:,\s*[^)]*)?\)$",
            row,
        )
        if match is None:
            return row
        src_expr = match.group(2).strip()
        src_match = re.fullmatch(r"QUEEN\.dna_dict\['([^']+)'\]", src_expr)
        if src_match is None:
            return row
        src_obj = live_products.get(src_match.group(1))
        if src_obj is None or getattr(src_obj, 'topology', None) != 'circular':
            return row
        return match.group(1) + src_expr

    def _normalize_ssdna_helper_row(row, live_products):
        match = re.search(r"^(QUEEN\.dna_dict\['([^']+)'\]\s*=\s*)(.+)$", row)
        if match is None:
            return row

        obj_name = match.group(2)
        obj = live_products.get(obj_name)
        if obj is None or hasattr(obj, "seq") is False or getattr(obj, "_ssdna", False) is False:
            return row

        seq = str(obj.seq)
        if len(seq) == 0:
            return row

        lhs = match.group(1)
        return "{}QUEEN(seq={}, ssdna={})".format(lhs, repr(seq), getattr(obj, "_ssdna", False))

    def _extract_load_alias(args_text):
        if args_text is None or "_load: " not in args_text:
            return None
        for item in args_text.split("; "):
            if item.startswith("_load: "):
                return item.split(": ", 1)[1]
        return None

    def _rewrite_load_alias_row(row, args_text):
        load_alias = _extract_load_alias(args_text)
        if load_alias is None or " = QUEEN(" not in row or "record=" not in row:
            return row

        match = re.search(r"^(QUEEN\.dna_dict\[')([^']+)('\]\s*=\s*QUEEN\()", row)
        if match is None:
            return row

        row = row[:match.start(2)] + load_alias + row[match.end(2):]
        if "product=" in row:
            row = re.sub(r"product='[^']*'", "product='{}'".format(load_alias), row, count=1)
            row = re.sub(r'product=\"[^\"]*\"', 'product="{}"'.format(load_alias), row, count=1)
        else:
            insert_pos = row.rfind(")")
            if insert_pos != -1:
                row = row[:insert_pos] + ", product='{}'".format(load_alias) + row[insert_pos:]
        return row

    def extract_qexd(rows): 
        extracted_rows = [] 
        for row in rows:
            if _row_has_qex_metadata(row) is False:
                if _is_seed_row(row) or _is_qexperiment_row(row) or _is_qexperiment_helper_row(row):
                    extracted_rows.append(_strip_qex_metadata(row))
                continue

            if "qexd = True" in row or "qexd=True" in row or "qexparam = True" in row or "qexparam=True" in row:
                continue

            extracted = _extract_qexd_row(row)
            if extracted is not None:
                extracted_rows.append(extracted)
        return extracted_rows

    def _prune_unused_rows(rows):
        dependency_counts = collections.defaultdict(int)
        for row in rows:
            if row.startswith("if __name__"):
                continue
            for match in re.findall(r"QUEEN\.dna_dict\[(?:'[^'\[\]]+'|\"[^\"\[\]]+\")\]", row):
                dependency_counts[match] += 1
            for match in re.findall(r"QUEEN\.queried_feature_dict\[(?:'[^'\[\]]+'|\"[^\"\[\]]+\")\]", row):
                dependency_counts[match] += 1
            for match in re.findall(r"QUEEN\.queried_features_dict\[(?:'[^'\[\]]+'|\"[^\"\[\]]+\")\]", row):
                dependency_counts[match] += 1

        pruned_rows = []
        for row in rows:
            lhs_match = re.match(r"^(QUEEN\.(?:dna_dict|queried_feature_dict|queried_features_dict)\['[^\[\]]+'\])\s*=", row.strip())
            if lhs_match is None:
                pruned_rows.append(row)
                continue
            lhs = lhs_match.group(1)
            lhs_alt = lhs.replace("['", '["').replace("']", '"]')
            if dependency_counts.get(lhs, 0) + dependency_counts.get(lhs_alt, 0) > 1:
                pruned_rows.append(row)
        return pruned_rows

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
        history[1] = _rewrite_load_alias_row(history[1], history[2])
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
    new_new_args = []
    row_args_iter = iter([history[2] for history in new_histories])
    for row in new_rows: 
        args_text = next(row_args_iter) if row.startswith("QUEEN.") else None
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
        new_new_args.append(args_text)
    
    live_products = getattr(dnas[0].__class__, "_products", {})
    if qexperiment_only == True:
        new_rows = extract_qexd(new_new_rows)
        new_rows = [_normalize_record_path_row(_normalize_identity_modifyends_row(_normalize_none_join_row(row), live_products)) for row in new_rows]
        #new_rows = _prune_unused_rows(new_rows)
    else:
        new_rows = []
        for row, args_text in zip(new_new_rows, new_new_args):
            if _is_qexperiment_row(row):
                continue
            row = _normalize_record_path_row(_normalize_identity_modifyends_row(_normalize_none_join_row(row), live_products))
            row = _rewrite_lower_row_from_args(row, args_text, live_products)
            new_rows.append(_strip_qex_metadata(row))

    project_names = []
    if description_only == False:
        now = datetime.datetime.now()
        if _return_script == False:
            print("project='{}'".format(project), file=o)
            print("import sys", file=o)  
            print("sys.path = [\"{}] + sys.path".format("/".join(__file__.split("/")[:-2])  + "\""), file=o)
            print("from QUEEN.queen import *", file=o) 
            print("import cutsite as cs", file=o) 
            custom_cutsites = {}
            for cutsite_name, _ in list(cs.new_cutsites):
                custom_cutsites[cutsite_name] = repr(cs.lib[cutsite_name].cutsite)
            for cutsite_name in sorted(custom_cutsites):
                print("cs.lib[{}] = {}".format(repr(cutsite_name), custom_cutsites[cutsite_name]), file=o) 
            if dna.__class__._namespaceflag == 1 and execution == False:
                print("set_namespace(globals())", file=o)
            print("", file=o) 
        
        scripts = [] 
        for n, row in enumerate(new_rows):
            row = re.sub(r",\s*,", ", ", row)
            row = re.sub(r"\(\s*,", "(", row)
            row = re.sub(r",\s*\)", ")", row)
            row = _normalize_record_path_row(row)
            row = re.sub(r"(?<=, )follow_order(?=,|\))", "follow_order=True", row)
            row = row.replace("follow_order='True'", "follow_order=True")
            row = row.replace("follow_order='False'", "follow_order=False")
            row = row.replace('follow_order=\"True\"', "follow_order=True")
            row = row.replace('follow_order=\"False\"', "follow_order=False")
            row = _normalize_ssdna_helper_row(row, live_products)
            primer_match = re.search(r"^QUEEN\.dna_dict\['([^']+)'\] = QUEEN\(seq=''[^)]*\)", row)
            if primer_match is not None:
                primer_name = primer_match.group(1)
                primer_obj = live_products.get(primer_name)
                if primer_obj is not None and hasattr(primer_obj, "seq") and len(primer_obj.seq) > 0:
                    primer_init = "QUEEN(seq={}, ssdna={})".format(repr(str(primer_obj.seq)), getattr(primer_obj, "_ssdna", False))
                    row = re.sub(r"QUEEN\(seq=''[^)]*\)", primer_init, row, count=1)
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
            if re.match(r"^QUEEN\.dna_dict\['[^']+'\]\s*=", row.strip()):
                last_line = row.strip()
            elif n == len(new_rows) - 1 and 'last_line' not in locals():
                last_line = row.strip() 
            
            if _return_script == False:
                print(row, file=o)
        
        result_match = re.search(r"^QUEEN\.dna_dict\['([^\[\]]+)'\]\s*=", last_line) if 'last_line' in locals() and last_line is not None else None
        result_expr = "QUEEN.dna_dict[{}]".format(repr(result_match.group(1))) if result_match is not None else result
        print("if __name__ == '__main__':", file=o) 
        print("    " + result_expr + ".outputgbk()", file=o)   
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
        script_path = outname
        if script_path[-3:] != ".py":
            os.rename(script_path, script_path + ".py")
            script_path = script_path + ".py"

        match = re.search(r"QUEEN.dna_dict\['([^\[\]]+)'\] = ", last_line) if last_line is not None else None
        key = match.group(1) if match is not None else dnas[0]._product_id

        proc, reconstructed_path = _execute_quine_script_isolated(script_path, key, cwd=os.getcwd())
        if type(output) is tempfile._TemporaryFileWrapper and os.path.exists(script_path):
            os.remove(script_path)

        dnas[0].__class__._source = None
        if proc.returncode != 0:
            if os.path.exists(reconstructed_path):
                os.remove(reconstructed_path)
            if _io == True:
                stderr = proc.stderr[-2000:] if proc.stderr is not None else ""
                raise ValueError(
                    "The {} QUEEN object could not be reconstructed using its quine code in an isolated subprocess. STDERR tail:\n{}".format(
                        dnas[0].project,
                        stderr,
                    )
                )
            return False

        reconstructed = dnas[0].__class__(record=reconstructed_path, dbtype="local", import_history=False)
        if os.path.exists(reconstructed_path):
            os.remove(reconstructed_path)

        original_seq = str(dnas[0].seq)
        reconstructed_seq = str(reconstructed.seq)
        same_seq = original_seq == reconstructed_seq
        if same_seq is False and getattr(dnas[0], "topology", None) == "circular" and getattr(reconstructed, "topology", None) == "circular" and len(original_seq) == len(reconstructed_seq):
            same_seq = original_seq in (reconstructed_seq + reconstructed_seq)

        if same_seq:
            if _io == True:
                print("QUEEN object reconstructed from the quine code and the original QUEEN object are identical.".format(dnas[0].project))
            return True, dnas[0]
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
                print("- Sample(s): {}; {}".format(sample, ", ".join([len(qobjects[asample].seq) for asample in pro_names[1:]])), file=output) 
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
                print("Output:\n{}".format(product), file=output)  
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
                print("- Sample(s): {}; {}".format(sample, ", ".join([str(len(qobjects[asample].seq)) + " bp" for asample in pro_names[1:]])), file=output) 
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
