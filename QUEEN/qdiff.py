import os 
import regex as re 
import sys 

def get_processsets(histories):
    #An element of process_set is composed of (sourcenames, func, productnames, process_id) 
    new_histories = []  
    for h, hisotry in enumerate(histories):
        info = history[2] 
        history = history[1]
        match0  = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = (.*)",history)
        match1  = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = ", history)
        match2  = re.search("product='([^=]+)'[,\)]", history)
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
            new_histories.append([history, info]) 
        
        elif match0 is not None:
            match3 = re.search("QUEEN.queried_features_dict\['([^\[\]]+)'\]", match0.group(1))
            if match2 is not None:
                name_dict[match3.group(0)] = match2.group(1) 
            else:
                name_dict[match3.group(0)] = match3.group(1) 

            unique_name_dict[match3.group(0)] = match3.group(1) 
            unique_name_name_dict[match3.group(1)] = name_dict[match3.group(0)] 
            new_histories.append([history, info]) 
        else:
            pass 
    
    process_sets = [] 
    for h, history in enumerate(new_histories):
        process_sets.append([[],[],[]])  
        info    = history[1] 
        history = history[0]
        matcho = re.search("(.*QUEEN.dna_dict\['[^\[\]]+'\]) = (.*)", history) 
        matchs = re.search("(QUEEN.queried_features_dict\['[^\[\]]+'\]) = (.*)", history)
        
        sourcenames = []  
        if matcho is not None:
            source = matcho.group(2) 
        else:
            source = matchs.group(2)  
        
        for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", source):
            key = key.group(0)
            sourcename = unique_name_dict[key]
            sourcenames.append(sourcename)
            process_sets[-1][0].append(("queen", sourcename)) 
            history =  history.replace(key, "queen") 
        
        #for key in re.finditer("QUEEN.queried_features_dict\['[^\[\]]+'\]", source):
        #    key = key.group(0)
        #    sourcename = unique_name_dict[key]
        #    sourcenames.append(sourcename)
        #    process_sets[-1][0].append(("dnafeature", sourcename)) 
        #    history = history.repalce(key, "dnafeature") 

        if matcho is not None:
            productnames = [] 
            product = matcho.group(1)
            for key in re.finditer("QUEEN.dna_dict\['([^\[\]]+)'\]", product):
                key = key.group(0) 
                productname = unique_name_dict[key]
                productnames.append(productname) 
                process_sets[-1][-1].append(("queen", productname)) 
                history = history.repalce(key, "queen") 
                
        else:
            productnames = [] 
            product = matchs.group(1)
            for key in re.finditer("QUEEN.queried_features_dict\['[^\[\]]+'\]", product):
                key = key.group(0) 
                productname = unique_name_dict[key]
                productnames.append(productname) 
                process_sets[-1][-1].append(("dnafeature", productname)) 
                history = history.repalce(key, "dnafeature")
        
        process_sets[-1][1].append(history) 
    return process_sets 

