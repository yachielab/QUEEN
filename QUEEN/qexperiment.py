import copy
import itertools as it
from qfunction import joindna, cropdna, cutdna, flipdna, modifyends
from qobj import QUEEN 
from qseq import Qseq 
import cutsite as cs
from cutsite import Cutsite
from Bio.SeqUtils import MeltingTemp as mt
import functools

def pcr(template, fw, rv, bindnum=15, mismatch=1, endlength=3, add_primerbind=False, tm_func=None, return_tm=False, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """
    Simulates a PCR (Polymerase Chain Reaction) process on a given DNA template using forward and reverse primers. This function does not provide the function to check cross dimer and homo dimer in primer design as default. If you wanna add such function, set the original `requirement` equiation. 
    Parameters
    ----------
    template : QUEEN
        The DNA template to be amplified.
    fw : QUEEN or str
        The forward primer. Can be a QUEEN object or a string representing the DNA sequence.
        If a dsDNA QUEEN object is given, its top strand sequence is used as a foward primer.
    rv : QUEEN or str
        The reverse primer. Can be a QUEEN object or a string representing the DNA sequence.
        If a dsDNA QUEEN object is given, its top strand sequence is used as reverse primer.
    bindnum : int, optional
        The minimum number of binding nucleotides for a primer, by default 15.
    mismatch : int, optional
        The maximum number of mismatches allowed in the primer binding, by default 1.
    endlength : int, optional
        The length of the end region of the primer to consider during binding, by default 3.
    add_primerbind : bool, optional
        If True, add DNAfeatures on the primer binding sites in the template DNA. 
    tm_func : str, function, optional
        Function to calculate the melting temperature of the primer pair 
        As `str` specfication, you can select `"Breslauer"` and `"SantaLucia"`. Default is `"SantaLucia"`.  
        Also, as built-in algorithms, `QUEEN.qexperiment.Tm_NN()`. This function is implemented 
        based on the `Bio.SeqUtils.MeltingTemp.Tm_NN()`, so the all parameters of 
        `Bio.SeqUtils.MeltingTemp.Tm_NN()`, excluding `seq` and `c_seq`, can be acceptable.
    return_tm : bool, optional
        If True, tm values of the primer pair are also returned.  
        The tm values will be calculated based on thier biding region excluding adaptor regions.
    product : str, optional 
        Product name of the PCR process.
    process_name : str, optional
        Brief label for the PCR process, by default "PCR"
    process_description : str, optional
        Additional description for the PCR process.
    pn : str, optional
        Alias for process_description.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    -------
    QUEEN (amplicon)
    Returns the PCR product (amplicon) as a QUEEN object.

    Examples
    --------
    >>> template = QUEEN("example_dna_sequence")
    >>> forward_primer = "xxxxxxxxx" #Specify a proper DNA sequence.
    >>> reverse_primer = "xxxxxxxxx" #Specify a proper DNA sequence.
    >>> pcr_product = pcr(template, forward_primer, reverse_primer)
    
    Notes
    -----
    The function internally uses a nested function `search_binding_site` to locate the binding sites of the primers.
    It performs DNA cropping and modification to simulate the PCR process. Errors are raised for invalid inputs or if 
    suitable binding sites are not found.

    See Also
    --------
    QUEEN, cropdna, modifyends

    """
    def search_binding_site(template, primer, strand=1, bindnum=15, endlength=3, flag=1, pn=None, pd=None, **kwargs): 
        site = [] 
        primer_end = primer.seq[-1*endlength:]
        for i in range(bindnum-endlength, len(primer.seq)-endlength+1):
            binding_site = primer.seq[-1*i + -1*endlength:-1*endlength]
            try:
                site = template.searchsequence(query="(?:{}){{s<={}}}{}".format(binding_site, mismatch, primer_end), quinable=False) 
            except Exception as e:
                site = [] 
            if len(primer_end) + len(binding_site) - mismatch >= bindnum:
                break
            else:
                pass 
        
        if len(site) == 1: #and site[0].strand == strand:
            site = template.searchsequence(query="(?:{}){{s<={}}}{}".format(binding_site, mismatch, primer_end), pn=pn, pd=pd)
            if flag == 1:
                return site[0]
            else:
                return site
        
        elif len(site) == 0:
            if flag == 1:
                raise ValueError("No primer binding sites were found.")
            else:
                return site 

        elif len(site) > 1: 
            raise ValueError("Multiple primer binding sites were detected.") 
        
        #else:
        #    raise ValueError("Primer binded to an unexpected strand.") 
    
    if type(template) != QUEEN: 
        raise TypeError("`template` object must be instance of QUEEN class.") 

    if type(fw) == Qseq or type(fw) == str:
        fw    = QUEEN(seq=fw, ssdna=True) 
        fwstr = fw.seq
    elif type(fw) == QUEEN:
        fwstr = "QUEEN.dna_dict['{}']".format(fw._product_id)
    else:
        raise TypeError("`fw` object must be instance of QUEEN or str class.") 

    if type(rv) == Qseq or type(rv) == str:
        rv    = QUEEN(seq=rv, ssdna=True) 
        rvstr = rv.seq
    elif type(rv) == QUEEN:
        rvstr = "QUEEN.dna_dict['{}']".format(rv._product_id)
    else:
        raise TypeError("`fw` object must be instance of QUEEN or str class.") 
    
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "PCR"
    
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    qexd = "pcr(QUEEN.dna_dict['{}'], {}, {}, bindnum={}, mismatch={}, endlegth={}, add_primerbind={}{})".format(template._product_id, fwstr, rvstr, bindnum, mismatch, endlength, add_primerbind, kwargs_str)

    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description
    if -1 in [template._left_end_top, template._left_end_bottom, template._right_end_top, template._right_end_bottom]:
        template = modifyends(template, pn=process_name, pd=process_description)
    
    if return_tm == True or add_primerbind == True:
        i = 0
        tmpsite1 = None
        while 1:
            tmpsite1 = search_binding_site(template, fw, 1, bindnum+i, endlength, flag=0, pn=process_name, pd=process_description)
            if len(tmpsite1) == 0:
                break
            site1 = tmpsite1[0]
            i += 1

        i = 0
        tmpsite2 = None
        while 1:
            tmpsite2 = search_binding_site(template, rv, -1, bindnum+i, endlength, flag=0, pn=process_name, pd=process_description)
            if len(tmpsite2) == 0:
                break
            site2 = tmpsite2[0] 
            i += 1
    else:
        site1 = search_binding_site(template, fw, 1, bindnum, endlength, pn=process_name, pd=process_description) 
        site2 = search_binding_site(template, rv, -1, bindnum, endlength, pn=process_name, pd=process_description) 
    
    if site1.strand == 1 and site2.strand == -1:
        fw_site = site1  
        rv_site = site2
    elif site1.strand == -1 and site2.strand == 1:
        fw, rv = rv, fw 
        fw_site = site2
        rv_site = site1
    else:
        raise ValueError("Both primers binded to the same strand.") 

    fw_bind_length = len(fw_site.sequence) 
    rv_bind_length = len(rv_site.sequence) 
    
    fw_feats = fw.searchfeature(key_attribute="feature_type", query="primer_bind") 
    rv_feats = rv.searchfeature(key_attribute="feature_type", query="primer_bind")
       
    if len(fw_feats) == 0:
        fw.setfeature({"qualifier:label":"{}".format(fw.project), "feature_type":"primer_bind"})  
    
    if len(rv_feats) == 0:
        rv.setfeature({"qualifier:label":"{}".format(rv.project), "feature_type":"primer_bind"})
    
    if fw_site.end > rv_site.start and fw_site.end <= rv_site.end:
        extract   = cropdna(template, rv_site.start,  fw_site.end, pn=process_name, pd=process_description)
        fw_index = len(fw.seq) - (fw_site.end - rv_site.start) 
        rv_index = fw_site.end - rv_site.start 
        amplicon = modifyends(extract, fw.seq[:fw_index], rv.rcseq[rv_index:], product=product, pn=process_name, pd=process_description, qexparam=qexd)
    else:
        extract  = cropdna(template, fw_site.end, rv_site.start, pn=process_name, pd=process_description)
        amplicon = modifyends(extract, fw.seq, rv.rcseq, product=product, pn=process_name, pd=process_description, qexparam=qexd)
    
    if add_primerbind == True:
        template.setfeature({"start": fw_site.start, "end": fw_site.end, "strand":1,  "feature_type":"primer_bind"}) 
        template.setfeature({"start": rv_site.start, "end": rv_site.end, "strand":-1, "feature_type":"primer_bind"})  
    
    if return_tm == True:
        if tm_func is None:
            tm_func = Tm_NN() 
        
        elif tm_func == "SantaLucia":
            tm_func = Tm_NN(nn_table=mt.DNA_NN3)

        elif tm_func == "Breslauer":
            tm_func = Tm_NN(nn_table=mt.DNA_NN1)

        fw_tm = tm_func(seq=template.seq[fw_site.start:fw_site.end]) 
        rv_tm = tm_func(seq=template.seq[rv_site.start:rv_site.end])
        return amplicon, (fw_tm, rv_tm) 
    else:
        return amplicon

def _select(fragments, selection=None): 
    if selection is None:
        return fragments  
    
    elif type(selection) == int:
        fragments.sort(key=lambda x:abs(len(x.seq)-selection))
        return fragments[0] 

    elif type(selection) == tuple:
        fragments = [fragment for fragment in fragments if min(size_range) <= len(fragment.seq) <= max(size_range)] 
        if len(fragments) > 1:
            raise ValueError("Multiple fragments holding the specified feature were detected") 
        return fragments[0] 

    elif selection in ("min", "max"):
        fragments.sort(key=lambda x:len(x.seq))
        if selection == "min":
            return fragments[0]
        else:
            return fragments[-1]
    
    elif selection.startswith("!") == False and ":" in selection: 
        query = ":".join(selection.split(":")[1:])
        fragments = [fragment for fragment in fragments if len(fragment.searchfeature(key_attribute="qualifier:{}".format(selection.split(":")[0]), query=query)) > 0]
        if len(fragments) > 1:
            raise ValueError("Multiple fragments holding the specified feature were detected") 
        return fragments[0]
    
    elif selection.startswith("!") and ":" in selection: 
        query = ":".join(selection.split(":")[1:])
        fragments = [fragment for fragment in fragments if len(fragment.searchfeature(key_attribute="qualifier:{}".format(selection.split(":")[0][1:]), query=query)) == 0]
        if len(fragments) > 1:
            raise ValueError("Multiple fragments holding the specified feature were detected") 
        return fragments[0] 

def digestion(dna, *cutsites, selection=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """
    Simulates a digestion of a DNA sequence using specified restriction enzymes (cutsites). 
    Optionally filters the resulting DNA fragments based on size.

    Parameters
    ----------
    dna : QUEEN
        The DNA sequence to be digested.
    *cutsites : Cutsite or str
        Variable number of Cutsite objects or names representing the restriction enzymes used for digestion.
    selection : "min", "max", "label:*", "!label:*", or tuple of int, optional
        The rule to select a specific digested fragment with the specified condition. 
        If "min" is provided, the minimum fragment of the digested fragments would be returned.
        If "max" is provided, the maximum fragment of the digested fragments would be returned.
        If "label:{feature_of_interest}" is provided, the unique fragment holding the DNAfeature with 
        `feature_of_interest` in "qualifer:label" would be returned. If multiple fragments holding the 
        specified feature are detected, a error will be raised.
        If "!label:{feature_of_interest}" is provided, the unique fragment not holding the DNAfeature with 
        `feature_of_interest` in "qualifer:label" would be returned. If multiple not fragments holding the 
        specified feature are detected, a error will be raised.
        If a `tuple` value is provided, The tuple (min_size, max_size) specifies the size range for filtering 
        the resulting fragments. If multiple fragments holding the are detected in the specified range, 
        a error will be raised.  
        If a `int` value is specified, the fragment with the nearest length to the specified value will be selected. 
        If None, no selection is done. Default is None.  
    product : str, optional 
        Product name of the digestion process.
    process_name : str, optional
        Brief label for the digestion process, by default "Digestion"
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the digestion process.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    -------
    list of QUEEN or QUEEN
        If `selection` is None, return list of QUEEN objects composed of the digested fragments.
        Otherwise, return a specific fragment filling the specified condition.
            
    Examples
    --------
    >>> dna_sequence = QUEEN("example_dna_sequence")
    >>> cutsite1 = Cutsite("restriction_enzyme_1")
    >>> cutsite2 = Cutsite("restriction_enzyme_2")
    >>> fragments = digestion(dna_sequence, cutsite1, cutsite2, selection=(100, 1000))

    Notes
    -----
    The function performs DNA digestion by searching for cut sites in the provided DNA sequence. 
    It then uses these sites to simulate cutting the DNA. Error handling is implemented to ensure 
    that the provided cut sites are valid Cutsite instances.

    See Also
    --------
    QUEEN, Cutsite, cutdna

    """
    if selection is not None:
        if (type(selection) not in (int, tuple)) and selection not in ("min", "max") and (selection.startswith("label:") == False and selection.startswith("!label:") == False):
            raise TypeError("`selection` should be `tuple` value, 'min', 'max', or `str` starting with 'label:' or '!label'.")
    
    cutsite_names = []
    cutsites = list(cutsites)
    for c in range(len(cutsites)):   
        if type(cutsites[c]) == str and cutsites[c] in cs.lib.keys():
            cutsites[c] = cs.lib[cutsites[c]]  
        elif type(cutsites[c]) == Cutsite or "cutsite" in cutsites[c].__dict__:
            pass 
        else:
            raise TypeError("Each element in `cutsites` must be instance of Cutsite class or its name must be included in `QUEEN.cutsite.lib`.")
        cutsite_names.append(cutsites[c].name) 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Digestion"
    
    cs_str = ", ".join(["'{}'".format(name) for name in cutsite_names])
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    if type(selection) == tuple:
        qexd = "digestion(QUEEN.dna_dict['{}'], {}, selection=[{}]{})".format(dna.project, cs_str, ",".join(map(str, selection)), kwargs_str) 
    else:
        qexd = "digestion(QUEEN.dna_dict['{}'], {}, selection='{}'{})".format(dna.project, cs_str, selection, kwargs_str) 
    
    if "qexparam" not in kwargs:
        pass
    else:
        qexd = kwargs["qexparam"]  
        del kwargs["qexparam"] 

    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description

    new_cutsites = [] 
    for cutsite in cutsites:
        sites = dna.searchsequence(query=cutsite, pn=process_name, pd=process_description)        
        new_cutsites.extend(sites) 
    
    fragments = cutdna(dna, *new_cutsites, product=product, pn=process_name, pd=process_description, qexparam=qexd)
    if len(fragments) == 1 and selection is None: 
        return _select(fragments, "max") 
    else:
        return _select(fragments, selection) 

def ligation(*fragments, unique=True, follow_order=False, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs): 
    """
    Simulates a ligation of DNA fragments, assembling them in various combinations and orientations.
    Can return either unique or multiple assembled DNA constructs.

    Parameters
    ----------
    *fragments : QUEEN
        Variable number of QUEEN objects representing DNA fragments to be ligated. 
    unique : bool, optional
        If True, ensures that only a unique assembled construct is returned. If multiple constructs 
        are possible, raises an error. Default is True.
    follow_order : bool, optional 
        If True, a ligation reaction will be simulated along with the given order of fragments.
        Default is False. 
    product : str, optional 
        Product name of the ligation process
    process_name : str, optional
        Brief label for the ligation process, by default "Ligation".
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the ligation process.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    -------
    QUEEN or list of QUEEN
        If `unique` is True and only one construct is possible, returns that construct as a QUEEN object.
        If `unique` is False, returns a list of all possible assembled constructs as QUEEN objects.
        If no constructs are possible, returns an empty list or None.

    Examples
    --------
    >>> fragment1 = QUEEN("dna_fragment_1")
    >>> fragment2 = QUEEN("dna_fragment_2")
    >>> assembled_product = ligation(fragment1, fragment2, unique=True)

    Notes
    -----
    The function attempts all permutations and orientations of the given fragments for ligation.
    If `unique` is set to True, it validates the uniqueness of the assembled product. The function
    handles situations where the assembly is not possible or results in multiple products.

    See Also
    --------
    QUEEN, flipdna, joindna

    """
    
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Ligation"  
    
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    fragments_str = ", ".join(["QUEEN.dna_dict['{}']".format(fragment._product_id) for fragment in fragments])
    
    if "qexparam" not in kwargs:
        qexd = "ligation({}, unique={}, follow_order={}{})".format(fragments_str, unique, follow_order, kwargs_str)
    else:
        qexd = kwargs["qexparam"]  
        del kwargs["qexparam"] 

    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description #"\n" + pd_suffix
    
    if follow_order == True:
        pass
    else:
        indexes_list = []  
        fragment1 = fragments[0] 
        f = 0
        orders = [] 
        flips  = []
        indexes = list(range(len(fragments))) 
        indexes.remove(f)
        orders.append(f) 
        flips.append(1) 
        while len(indexes) > 0:
            flag = 0 
            for g in indexes:
                fragment2 = fragments[g]
                rl = fragment1._right_end_top * fragment2._left_end_bottom 
                rr = fragment1._right_end_top * fragment2._right_end_top
                if rl == 1 and fragment1._right_end == fragment2._left_end: 
                    flag = 1
                    orders.append(g) 
                    flips.append(1) 
                    fragment1 = fragment2
                    indexes.remove(g) 
                    break
                elif rr == 1 and fragment1._right_end == fragment2._right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]:
                    flag = 1
                    orders.append(g) 
                    flips.append(-1)
                    fragment1 = flipdna(fragment2, quinable=0) 
                    indexes.remove(g) 
                    break  
            if flag == 0:
                break
            else:
                pass 
        
        if len(indexes) == 0:
            indexes_list.append((orders, flips))
    
    if follow_order == True:
        outobj = joindna(*fragments, topology="circular", autoflip=False, compatibility="complete", product=product, pn=process_name, pd=process_description, qexparam=qexd)
        return outobj

    elif unique == True:
        if len(indexes_list) == 1:
            orders, flips = indexes_list[-1]
            fragment_set  = [flipdna(fragments[ind], product=fragments[ind].project, pn=process_name, pd=process_description) if fl == -1 else fragments[ind] for ind, fl in zip(orders, flips)]
            #print([fragment.project for fragment in fragment_set])
            outobj = joindna(*fragment_set, topology="circular", autoflip=False, compatibility="complete", product=product, pn=process_name, pd=process_description, qexparam=qexd)
        elif len(indexes_list) == 0:
            raise ValueError("The QUEEN_objects cannot be joined due to the end structure incompatibility.") 
        else:
            raise ValueError("Multiple different constructs will be assembled. You should review your assembly design.")
        return outobj
    
    else:
        products = [] 
        for order, flips in indexes_list:
            fragment_set  = [flipdna(fragments[ind], product=fragments[ind].project, pn=process_name, pd=process_description) if fl == -1 else fragments[ind] for ind, fl in zip(orders, flips)]
            outobj = joindna(*fragment_set, topology="circular", autoflip=False, compatibility="complete", product=product, pn=process_name, pd=process_description, qexparam=qexd)
            products.append(outobj)
        return products 

def homology_based_assembly(*fragments, mode="gibson", homology_length=15, unique=True, follow_order=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs): #homology_based_assembly
    """
    Simulates a homology-based DNA assembly, supporting various modes like Gibson, Infusion, or Overlap PCR.

    Parameters
    ----------
    *fragments : QUEEN
        Variable number of QUEEN objects representing DNA fragments to be assembled.
    mode : str, optional
        The assembly mode to be used. Valid options are "gibson", "infusion", and "overlap pcr". Default is "gibson".
    homology_length : int, optional
        The minimum length of homology required for the assembly. Default is 20.
    unique : bool, optional
        If True, ensures that only a unique assembled construct is returned. If multiple constructs 
        are possible, raises an error. Default is True.
    follow_order : bool, optional 
        If True, a ligation reaction will be simulated along with the given order of fragments.  
        If the number of given fragments is larger than 4, default is True. Otherwise, False. 
    product : str, optional 
        Product name of the homology_based_assembly (hba) process. 
    process_name : str, optional
        Brief label for the `homology_based_assembly` process.
        If `mode` is `"gibson"`, default is "Gibson Assembly". 
        If `mode` is `"infusion"`, default is "In-Fusion Assembly". 
        If `mode` is `"overlappcr"`, default is "Overlap PCR".
    process_description : str, optional
        Additional description for the assembly process.
    pn : str, optional
        Alias for process_description.
    pd : str, optional
        Alias for process_name.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    -------
    QUEEN or list of QUEEN
        If `unique` is True and only one construct is possible, returns that construct as a QUEEN object.
        If `unique` is False, returns a list of all possible assembled constructs as QUEEN objects.
        If no constructs are possible, returns an empty list or None.

    Examples
    --------
    >>> fragment1 = QUEEN("dna_fragment_1")
    >>> fragment2 = QUEEN("dna_fragment_2")
    >>> assembled_product = homology_based_assembly(fragment1, fragment2, mode="gibson", unique=True)

    Notes
    -----
    The function considers different assembly modes, each with its specific requirements for fragment
    orientation and homology lengths. It handles permutations and orientations of the given fragments.
    Error handling is implemented for invalid inputs or assembly modes.

    See Also
    --------
    QUEEN, flipdna, joindna, modifyends

    """
    max_homology_length = 500 #max_homology_length
    
    process_name = pn if process_name is None else process_name
    if process_name is None:
        if mode == "gibson": 
            process_name = "Gibson Assembly"
        elif mode == "infusion":
            process_name = "In-Fusion Assembly"
        elif mode == "overlappcr":
            process_name = "Overlap PCR"
        else:
            process_name = "Homology based Assembly" 
    
    kwargs_str    = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str    = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    fragments_str = ", ".join(["QUEEN.dna_dict['{}']".format(fragment._product_id) for fragment in fragments])
    qexd = "homology_based_assembly({}, mode='{}', homology_length={}, unique={}, follow_order={}{})".format(fragments_str, mode, homology_length, unique, follow_order, kwargs_str)
    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description
    
    for fragment in fragments:
        if fragment.topology != "linear": 
            raise ValueError("A 'circular' fragment was detected. All fragments to be assembled should be 'linear' topology.")

    if mode not in ("gibson", "infusion", "overlappcr"):
        raise ValueError("Invalid mode value. The 'mode' variable can only take 'gibson', 'infusion' or 'overlappcr' as values.")
    
    if mode == "overlappcr":
        for fragment in fragments:
            if fragment._left_end_top * fragment._left_end_bottom == -1 or fragment._right_end_top * fragment._right_end_bottom == -1:
                raise ValueError("dsDNA ojbects with sticy ends cannot be handled in `overlappcr`. please use `modifyends` function to refrom sticy ends to blunt ends.")
            else:
                pass 

        nums = list(range(len(fragments)))
        if (len(fragments) < 5 and follow_order is None) or follow_order == False:
            nums_orders     = list(it.permutations(nums))
            new_nums_orders = [] 
            for nums_order in nums_orders:
                if tuple(reversed(nums_order)) in new_nums_orders:
                    pass
                else:
                    new_nums_orders.append(nums_order) 
            nums_orders = new_nums_orders
            flip_status_list = list(it.product(*[[1,-1] for i in range(len(fragments))]))   
        else:
            nums_orders      = [nums]  
            flip_status_list = [[1 for i in range(len(framgents))]]    
        errors = [] 
        products = [] 
        for numset in nums_orders:
            execed = [] 
            for flipset in flip_status_list:
                if tuple([state * -1 for state in flipset]) in execed: 
                    pass 
                else:
                    fragment_set = [fragments[num] if flip == 1 else flipdna(fragments[num], pn=process_name, pd=process_description) for num, flip in zip(numset, flipset)] 
                    try:
                        fragment1 = fragment_set[0]
                        for f, fragment2 in enumerate(fragment_set[1:]): 
                            fragment1 = modifyends(fragment1, "-{{{}}}/*{{{}}}".format(len(fragment1.seq),len(fragment1.seq)), "*{{{}}}/-{{{}}}".format(len(fragment1.seq),len(fragment1.seq)), pn=process_name, pd=process_description)
                            fragment2 = modifyends(fragment2, "-{{{}}}/*{{{}}}".format(len(fragment2.seq),len(fragment2.seq)), "*{{{}}}/-{{{}}}".format(len(fragment2.seq),len(fragment2.seq)), pn=process_name, pd=process_description)
                            fragment1 = joindna(fragment1, fragment2, topology="linear", homology_length=homology_length, pn=process_name, pd=process_description) 
                            if f == len(fragment_set) - 2:
                                fragment1 = modifyends(fragment1, "*/*", "*/*", product=product, pn=process_name, pd=process_description, qexparam=qexd) 
                            else:
                                fragment1 = modifyends(fragment1, "*/*", "*/*", pn=process_name, pd=process_description) 
                        products.append(fragment1)
                    except Exception as e:
                        errors.append(e) 
                execed.append(flipset) 

    else:
        nums = list(range(len(fragments)))
        if (len(fragments) < 6 and follow_order is None) or follow_order == False:
            nums_orders = list(map(list,it.permutations(nums[:-1])))
            nums_orders = [numlist + [nums[-1]] for numlist in nums_orders]
            flip_status_list = list(it.product(*[[1,-1] for i in range(len(fragments))]))   
        else:
            nums_orders      = [nums]  
            flip_status_list = [[1 for i in range(len(framgents))]]    
        errors = [] 
        products = [] 
        for numset in nums_orders:
            execed = [] 
            for flipset in flip_status_list:
                if tuple([state * -1 for state in flipset]) in execed: 
                    pass 
                else:
                    fragment_set = [fragments[num] if flip == 1 else flipdna(fragments[num]) for num, flip in zip(numset, flipset)]  
                    for f in range(len(fragment_set)):
                        fragment = fragment_set[f]
                        if len(fragment.seq) <= max_homology_length: 
                            mhl = int(len(fragment.seq)) - len(fragment._left_end) - len(fragment._right_end) - 1
                        else:
                            mhl = max_homology_length
                        if mode == "gibson":
                            fragment_set[f] = modifyends(fragment, "-{{{}}}/*{{{}}}".format(mhl,mhl), "*{{{}}}/-{{{}}}".format(mhl,mhl), pn=process_name, pd=process_description)
                        elif mode == "infusion":
                            fragment_set[f] = modifyends(fragment, "*{{{}}}/-{{{}}}".format(mhl,mhl), "-{{{}}}/*{{{}}}".format(mhl,mhl), pn=process_name, pd=process_description)
                    try:
                        outobj = joindna(*fragment_set, autoflip=False, homology_length=homology_length, topology="circular", product=product, pn=process_name, pd=process_description, qexparam=qexd) 
                        products.append(outobj) 
                    except Exception as e:
                        errors.append(e) 
                        #fragment_set[0].printsequence(display=True, hide_middle=30) 
                        #fragment_set[1].printsequence(display=True, hide_middle=30)
                        pass 
                execed.append(flipset) 

    if unique == True:
        if len(products) > 1: 
            raise ValueError("Multiple assembled constructs were detected. You should review your assembly design.")
        else:
            try:
                return products[0]
            except Exception as e:
                print(e, errors) 
                raise ValueError("Error") 
    else: 
        return products 

def annealing(ssdna1, ssdna2, homology_length=4, product=None, pn=None, pd=None, process_name=None, process_description=None, **kwargs):
    """
    Simulates an annealing of two single-stranded DNA (ssDNA) molecules based on homology length.
    If dsDNA objects are given, their top strand will be used for the annealing.
    Parameters
    ----------
    ssdna1 : QUEEN
        The first single-stranded DNA molecule to be annealed. If a dsDNA QUEEN object is given, it will be regarded as ssDNA QUEEN object with its top strand sequence. 
    ssdna2 : QUEEN
        The second single-stranded DNA molecule to be annealed. If a dsDNA QUEEN object is given, it will be regarded as ssDNA QUEEN object with its top strand sequence.
    homology_length : int, optional
        The length of the homologous region required for annealing. Default is 4.
    product : str, optional 
        Product name of the ligation process
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Gateway Reaction".
    process_description : str, optional
        Additional description for the gateway reaction process.
    pn : str, optional
        Alias for process_name.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    -------
    QUEEN
        The resulting double-stranded DNA molecule after annealing.

    Examples
    --------
    >>> ssdna1 = QUEEN("ATCG")
    >>> ssdna2 = QUEEN("CGTA")
    >>> dsdna = annealing(ssdna1, ssdna2, homology_length=4)

    See Also
    --------
    QUEEN, joindna

    """
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Annealing" 
    
    kwargs_str    = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str    = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    qexd = "annealing(QUEEN.dna_dict['{}'], QUEEN.dna_dict['{}'], homology_length={}{})".format(ssdna1._product_id, ssdna2._product_id, homology_length, kwargs_str)
    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix

    if type(ssdna1) == str:
        ssdna1 = QUEEN(seq=ssdna1, ssdna=True)
    
    if type(ssdna2) == str:
        ssdna2 = QUEEN(seq=ssdna2, ssdna=True)
    
    flag1 = 0
    if ssdna1._ssdna == False: 
        flag1 = 1
        ssdna1._ssdna = True

    flag2 = 0
    if ssdna2._ssdna == False:
        flag2 = 1
        ssdna2._ssdna = True

    if type(ssdna1) != QUEEN:
        raise TypeError("`ssdna_top` object must be a QUEEN or str object") 
    
    if type(ssdna2) != QUEEN:
        raise TypeError("`ssdna_down` object must be a QUEEN or str object") 
    
    annealed_dna  = joindna(ssdna1, ssdna2, homology_length=homology_length, product=product, pn=process_name, pd=process_description, qexparam=qexd)
    
    ssdna1._ssdna       = False if flag1 == 1 else True 
    ssdna2._ssdna       = False if flag2 == 1 else True

    return annealed_dna 

def gateway_reaction(destination, entry, mode="BP", product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """
    Simulates a gateway reaction of two DNA molecules.
    Basic `BP` and `LR` reactions are available.

    Parameters
    ----------
    destination : QUEE object N
        The destination QUEEN object holding the backbone DNA molecule.
    entry : QUEEN object
        The entry QUEEN object holding the insert DNA molecule.
    mode: str, tuple, or list
        The mode of the reaction, can be "BP" or "LR". Default is "BP".
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Gateway Reaction".
    process_description : str, optional
        Additional description for the gateway reaction process.
    pn : str, optional
        Alias for process_name.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    ----------
    QUEEN
        The QUEEN object representing the result of the gateway reaction process.
    """
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Gateway Reaction" 
    
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    qexd = "gateway_reaction(QUEEN.dna_dict['{}'], QUEEN.dna_dict['{}'], mode='{}'{})".format(destination._product_id, entry._product_id, mode, kwargs_str)
    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix

    if mode == "BP":
        cs.lib["attX1"] = "ACAAGTTT^GTACAAA_AAAGCAGGCT" #attB1
        cs.lib["attX2"] = "ACCCAGCTTT^CTTGTAC_AAAGTGGT" #attB2
        cs.lib["attY1"] = "CCAACTTT^GTACAAA_AAAGCTGAAC" #attP1
        cs.lib["attY2"] = "GTTCAGCTTT^CTTGTAC_AAAGTTGG" #attP2 
    
    elif mode == "LR":
        cs.lib["attX1"] = "CCAACTTT^GTACAAA_AAAGCAGGCT" #attL1
        cs.lib["attX2"] = "ACCCAGCTTT^CTTGTAC_AAAGTTGG" #attL2
        cs.lib["attY1"] = "ACAAGTTT^GTACAAA_AAAGCTGAAC" #attR1
        cs.lib["attY2"] = "GTTCAGCTTT^CTTGTAC_AAAGTGGT" #attR2

    elif type(mode) in (tuple, list) and len(type) == 4:
        cs.lib["attX1"] = mode[0]  
        cs.lib["attX2"] = mode[1] 
        cs.lib["attY1"] = mode[2] 
        cs.lib["attY2"] = mode[3] 
        mode = "XY"

    else:
        ValueError("Basically,'mode' value can take only 'BP' or 'LR' reaction at present. For executing a custom BP or LR reaction, please speicy [B1 or L1_sequence, B1 or L2_sequence, P1 or R1 sequnce, P2 or R2 sequnece] along with the QUEEN's cutsite format.") 

    attx1 = entry.searchsequence(cs.lib["attX1"], product="att{}1_site".format(mode[0]), pn=process_name, pd=process_description) 
    attx2 = entry.searchsequence(cs.lib["attX2"], product="att{}2_site".format(mode[0]), pn=process_name, pd=process_description)
    atty1 = destination.searchsequence(cs.lib["attY1"], product="att{}1_site".format(mode[1]), pn=process_name, pd=process_description)
    atty2 = destination.searchsequence(cs.lib["attY2"], product="att{}1_site".format(mode[1]), pn=process_name, pd=process_description) 
    if len(attx1) > 1:
        raise ValueError("Multiple att{}1 sites were detected.".format(mode[0]))
    else:
        attx1 = attx1[0] 

    if len(attx2) > 1:
        raise ValueError("Multiple att{}2 sites were detected.".format(mode[0]))
    else:
        attx2 = attx2[0] 
            
    if len(atty1) > 1:
        raise ValueError("Multiple att{}1 sites were detected.".format(mode[1]))
    else:
        atty1 = atty1[0] 

    if len(atty2) > 1:
        raise ValueError("Multiple att{}2 sites were detected.".format(mode[1]))
    else:
        atty2 = atty2[0] 
    
    #entry_fragments = cutdna(entry, attx1, attx2, pn=process_name, pd=process_description) 
    #entry_fragments.sort(key=lambda x:len(x.seq))
    #insert = _select(entry_fragments, entry_selection)
    if attx1.strand == 1 and attx2.strand == 1:
        insert = cropdna(entry, attx1, attx2, pn=process_name, pd=process_description) 
    
    elif attx1.strand == -1 and attx2.strand == -1:
        insert = cropdna(entry, attx2, attx1, pn=process_name, pd=process_description)

    #destination_fragments = cutdna(destination, atty1, atty2, pn=process_name, pd=process_description, qexparam=qexd)      
    #destination_fragments._fragments.sort(key=lambda x:len(x.seq)) 
    #_select(destination_fragments, destination_selection)
    if atty1.strand == 1 and atty2.strand == 1:
        destination = cropdna(destination, atty2, atty1, pn=process_name, pd=process_description) 
    elif atty1.strand == -1 and atty2.strand == -1:
        destination = cropdna(destination, atty1, atty2, pn=process_name, pd=process_description) 

    return ligation(insert, destination, product=product, pn=process_name, pd=process_description, qexparam=qexd) 

def goldengate_assembly(destination, entry, cutsite=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """
    Simulates a Golden Gate Assembly

    Parameters
    ----------
    destination : QUEEN object
        The destination QUEEN object holding the backbone DNA molecule.
    entry : list of QUEEN objects
        The entry QUEEN objects holding the insert DNA molecules.
    cutsite : Cutsite or str
        The restriction enzyme site used for this reaction.    
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Golden Gate Assembly".
    process_description : str, optional
        Additional description for the gateway reaction process.
    pn : str, optional
        Alias for process_name.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.
    """
    if type(cutsite) == str and cutsite in cs.lib.keys():
        cutsite = cs.lib[cutsite]  
    elif type(cutsite) == Cutsite or "cutsite" in cutsite.__dict__:
        pass 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Golden Gate Assembly" 
    
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    entry_str  = ", ".join(["QUEEN.dna_dict['{}']".format(aentry._product_id) for aentry in entry])
    entry_str  = "[{}]".format(entry_str)
    qexd = "golden_gate_assembly(QUEEN.dna_dict['{}'], {}, cutsite='{}'{})".format(destination._product_id, entry_str, cutsite.name, kwargs_str)
    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix

    if type(entry) == QUEEN:
        entry = [entry]
    
    fragments = [] 
    for aentry in entry:
        if aentry.topology == "linear" and len(aentry.searchsequence(query=cutsite, quinable=0)) == 0:
            insert = aentry
        else: 
            inserts = digestion(aentry, cutsite, product=aentry.project, pn=process_name, pd=process_description, qexparam="")
            for insert in inserts:
                if cutsite.seq in insert.seq or cutsite.rcseq in insert.seq:
                    pass
                else:
                    break
        
        fragments.append(insert) 

    backbones =  digestion(destination, cutsite, product=destination.project, pn=process_name, pd=process_description, qexparam="") 
    for backbone in backbones:
        if cutsite.seq in backbone.seq or cutsite.rcseq in backbone.seq:
            pass
        else:
            break
    
    fragments.append(backbone)
    outobj = ligation(*fragments, product=product, pn=process_name, pd=process_description, qexparam=qexd) 
    return outobj

def topo_cloning(destination, entry, mode="TA", product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """
    Simulate TOPO cloning 
    
    Parameters
    ----------
    destination : QUEEN object
        The destination QUEEN object holding the backbone DNA molecule. Basically, the sequence topology should be linear and proper end structure along with the `mode` value. However, if the sequence topology is circular, the sequence would be processed pre-defined restriction enzymes based on `mode` value.  
    entry : QUEEN object with blunt-linear topology
        The entry QUEEN objects holding the insert DNA molecules. If `mode` is `"directional"`, its 5' end sequnence should be `"CACC"`. 
    mode : "TA", "blunt", or "directional"
        Brief label for the `topo_cloning` process.
    process_description : str, optional
        Additional description for the gateway reaction process.
    pn : str, optional
        Alias for process_description.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.
    
    Returns
    -------
    QUEEN (construct)
    Returns the topo cloning construct.
    """
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    qexd = "topo_cloning(QUEEN.dna_dict['{}'], QUEEN.dna_dict['{}'], mode='{}'{})".format(destination._product_id, entry._product_id, mode, kwargs_str)

    if mode == "TA":
        if destination.topology == "circular":
            destination = digestion(destination, "AflII", selection="max", product=destination.project, pn=process_name, pd=process_description, qexparam="")
            destination = modifyends(destination, pn=process_name, pd=process_description) 
            destination = modifyends(destination, "-/*", "*/-", pn=process_name, pd=process_description) 
        
        entry  = modifyends(entry, "T", "A", pn=process_name, pd=process_description) 
        entry  = modifyends(entry, "-/*", "*/-", pn=process_name, pd=process_description)
        outobj = joindna(destination, entry, autoflip=False, compatibility="complete", homology_length=1, topology="circular", product=product, pn=process_name, pd=process_description, qexparam=qexd)

    elif mode == "blunt": 
        if destination.topology == "circular":
            destination = digestion(destination, "AflII", selection="max", product=destination.project, pn=process_name, pd=process_description, qexparam="")
            destination = modifyends(destination, pn=process_name, pd=process_description) 
        outobj = joindna(destination, entry, autoflip=False, compatibility="complete", topology="circular", product=product, pn=process_name, pd=process_description, qexparam=qexd)

    elif mode == "directional":
        if destination.topology == "circular":
            destination = digestion(destination, "StyI", selection="max", product=destination.project, pn=process_name, pd=process_description, qexparam="")
            destination = modifyends(destination, pn=process_name, pd=process_description)
            destination = cropdna(destination, 1, len(destination.seq)-3, pn=process_name, pd=process_description)
            destination = modifyends(destination, "", "----/****", pn=process_name, pd=process_description)
        entry = modifyends(entry, "****/----", "", pn=process_name, pd=process_description)
        outobj = joindna(destination, entry, autoflip=False, compatibility="complete", topology="circular", product=product, pn=process_name, pd=process_description, qexparam=qexd) 
    
    return outobj

def intra_site_specific_recombination(dna, site="loxP", product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """
    Simulates a intra molecule site-specific recombination.

    Parameters
    ----------
    dna : QUEEN object
        The target QUEEN object.
    site : str ("loxP", "lox2272", "FRT") 
        The target site for the recombination reaction. At least two identical   
        recombination sites in the DNA sequence to simulate the recombination process.
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Golden Gate Assembly".
    process_description : str, optional
        Additional description for the gateway reaction process.
    pn : str, optional
        Alias for process_description.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.
    
    Returns
    -------
    QUEEN or list of QUEEN objects.
        If the number of sites in the given DNA molecule is two, returns a single QUEEN object.
        Otherwise (the number of sites > 2), returns multiple QUEEN objects as all possible  
        recombination results.

    Notes
    -----
    If there are three or more sites in the DNA molecule, the function simulates   
    all possible combinations of site-specific recombination.  
    However, only the first stage of recombination is considered in the simulation.   
    It means that this function does not simulate multiple stages of recombination.  
    """ 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Intra site-specific recombination" 
    
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    qexd = "intra_site_specific_recombination({}, site={}{})".format(dna.project, site, kwargs_str)
    process_description = pd if process_description is None else process_description
    #process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix
    
    if site == "loxP":
        cs.lib["recsite"] = "ATAACTTCGTATAA^TGTATG_CTATACGAAGTTAT"
    elif site == "lox2272":
        cs.lib["recsite"] = "ATAACTTCGTATAA^AGTATC_CTATACGAAGTTAT"
    elif site == "loxN":
        cs.lib["recsite"] = "ATAACTTCGTATAG^TATACC_TTATACGAAGTTAT"
    elif site == "FRT":
        cs.lib["recsite"] = "GAAGTTCCTATTC^TCTAGAAA_GTATAGGAACTTC"
    else:
        cs.lib["recsite"] = site 
        site = "custom_site"

    recsites    = dna.searchsequence(cs.lib["recsite"], product=site, pn=process_name, pd=process_description)
    outobj_list = [] 
    for combi in it.combinations(recsites, 2):
        recsite1 = combi[0]
        recsite2 = combi[1]
        fragments = cutdna(dna, recsite1, recsite2, pn=process_name, pd=proces_description)
        
        if recsite1.strand == recsite2.strand:
            if dna.topology == "circular":
                fragment = _select(fragments, selection="max") 
                outobj = joindna(fragment, autoflip=False, compatibility="complete", topology="circular", product=product, pn=process_name, pd=process_description, qexparam=qexd)
            else:
                outobj = joindna(fragment[0], fragment[2], autoflip=False, compatibility="complete", product=product, pn=process_name, pd=process_description, qexparam=qexd) 

        else:
            if dna.topology == "circular":
                fragment1 = _select(fragments, selection="max") 
                fragment2 = _select(fragments, selection="min") 
                fragment2 = flipdna(fragment2, pn=process_name, pd=process_description)
                outobj = joindna(fragment1, fragment2, autoflip=False, compatibility="complete", topology="circular", product=product, pn=process_name, pd=process_description, qexparam=qexd)
            else:
                reversed_fragment = flipdna(fragment[1], pn=process_name, pd=process_description) 
                outobj = joindna(fragment[0], reversed_fragment, fragment[2], autoflip=False, compatibility="complete", product=product, pn=process_name, pd=process_description, qexparam=qexd) 
        outobj_list.append(outobj)
    
    if len(outobj_list) == 1:
        return outobj_list[0]
    else:
        return outobj_list 

def homologous_recombination(donor, entry, left_homology=None, right_homology=None, homology_legnth=100, product=None, process_name=None, pn=None, process_description=None, pd=None):
    """
    Simulates a homology recombination of two DNA molecules.
    
    Parameters
    ----------
    destination : QUEEN
        The destination QUEEN object holding the backbone DNA molecule.
    entry : QUEEN
        The entry QUEEN object holding the insert DNA molecule.
    left_homology : str or QUEEN, optional
        The homology sequence at 5' side on the top strand.
        The sequence will be used as the homology arm in the HR reaction.  
        If the value is not given, a proper homology sequence will be automatically detected. 
    right_homology : str or QUEEN, optional 
        The homology sequence at 3' side on the top strand.
        The sequence will be used as the homology arm in the HR reaction.  
        If the value is not given, a proper homology sequence will be automatically detected. 
    homology_length : int, optional
        The minimum length of homology required for the assembly. Default is 100.
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Gateway Reaction".
    process_description : str, optional
        Additional description for the gateway reaction process.
    pn : str, optional
        Alias for process_name.
    pd : str, optional
        Alias for process_description.
    **kwargs
        Additional keyword arguments for advanced configurations.

    Returns
    ----------
    QUEEN
        The QUEEN object representing the result of the gateway reaction process.
    """ 
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Homologous Recombination" 
    
    kwargs_str = ", ".join(["{}".format(str(key)) + "=" + "'{}'".format(str(kwargs[key])) for key in kwargs])
    kwargs_str = kwargs_str if kwargs_str == "" else ", " +  kwargs_str
    qexd = "homologous_recombination(QUEEN.dna_dict['{}'], QUEEN.dna_dict['{}']{})".format(destination._product_id, entry._product_id, mode, destination_selection, entry_selection, kwargs_str)
    process_description = pd if process_description is None else process_description

    if left_homology is None or right_homology is None:
        dstrand = 1
        if destination.topology == "circular":
            region = len(destination.seq) 
        else:
            region = len(destination.seq) - 2*homology_length

        flag = 0 
        for i in range(region):
            left_homology  = destination.seq[i:i+homology_length] 
            right_homology = destination.seq[i+homology_length:i+2*homology_length]
            results = entry.searchsequence(query=left_homology + ".+" + right_homology, quinable=False)
            if len(results) == 1:
                dtarget = i + homolgy_length
                flag = 1
                break
            else:
                pass 
    else:
        dresultl = destination.searchseqeunce(query=left_homology, unique=True, pn=process_name, pd=process_description)[0] 
        dresultr = destination.searchseqeunce(query=right_homology, unique=True, pn=process_name, pd=process_description)[0] 
        if dresultl.strand == -1 or dresultr.strand == -1:
            raise ValueError("`left_homology` and `right_homology` sequeneces should be on the top strand.") 
        dstrand = dresultl.strand 
        dtarget = dresultl.end

        results  = entry.searchsequence(query=left_homology + ".+" + right_homology)
        if len(results) == 1:
            flag = 1

    if flag == 0:
        raise ValueError("Any proper homology arms were not detected.")
    
    eresult = results[0]
    estrand = eresult.strand 
        
    left_dest  = cropdna(destination, 0, dtarget, pn=process_name, pd=process_description)
    right_dest = cropdna(destination, dtarget, len(destination.seq), pn=process_name, pd=process_description) 
    if estrand == dstrand:
        estart  = eresult.start + len(left_homology)
        eend    = eresult.end - len(right_homology) 
        insert = cropdna(entry, estart, eend, pn=process_name, pd=process_description) 
    else:
        estart  = eresult.start + len(right_homology)
        eend    = eresult.end - len(left_homology) 
        insert = flipdna(cropdna(entry, estart, eend, pn=process_name, pd=process_description), pn=process_name, pd=process_description)

    if destination.topology == "circular":
        outobj = joindna(left_dest, insert, right_dest, topology="circular", autoflip=False, product=product, pn=process_name, pd=process_description)
    else:
        outobj = joindna(left_dest, insert, right_dest, autoflio=False, pn=process_name, product=prodcut, pd=process_description) 
    return outobj

def Tm_NN(check=True, strict=True, nn_table=None, tmm_table=None, imm_table=None, de_table=None, dnac1=25, dnac2=25, selfcomp=False, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5): 
    return functools.partial(mt.Tm_NN, check=check, strict=strict, nn_table=nn_table, tmm_table=tmm_table, imm_table=imm_table, de_table=de_table, dnac1=dnac1, dnac2=dnac2, selfcomp=selfcomp, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)  

def primerdesign(template, target, fw_primer=None, rv_primer=None, fw_margin=0, rv_margin=0,
                 target_tm=60.0, tm_func=None, primer_length=(16, 25), 
                 design_num=1, fw_adapter=None, rv_adapter=None, cloning_mode="gibson", 
                 homology_length=20, nonspecific_limit=3, auto_adjust=1, 
                 requirement=None, fw_name="fw_primer", rv_name="rv_primer"):
    """
    Design forward and reverse primers for PCR amplification of a target region,
    allowing for introduction of specific mutations, checking primer specificity,
    and meeting additional user-defined requirements.

    Parameters
    ----------
    template : QUEEN object
        The QUEEN object to serve as the PCR template.
    target : QUEEN object
        The region of the template QUEEN object that needs to be included in the amplicon.
    fw_primer : ssDNA QUEEN object, optional
        If provided, this sequence will be used as the forward primer.
    rv_primer : ssDNA QUEEN object, optional
        If provided, this sequence will be used as the reverse primer.
    fw_margin : int, optional
        Additional base pairs to add to the 5' end of the target region when designing the forward primer. Default is 0.
    rv_margin : int, optiona
        Additional base pairs to add to the 3' end of the target region when designing the reverse primer. Default is 0.
    target_tm : float, optional
        Desired melting temperature (Tm) for the primers in degrees Celsius. Default is 60.0.
    tm_func : str, function, optional
        Function to calculate the melting temperature of primer candidates. 
        As `str` specfication, you can select `"Breslauer"` and `"SantaLucia"`. Default is `"SantaLucia"`.  
        Also, as built-in algorithms, `QUEEN.qexperiment.Tm_NN()`. This function is implemented 
        based on the `Bio.SeqUtils.MeltingTemp.Tm_NN()`, so the all parameters of 
        `Bio.SeqUtils.MeltingTemp.Tm_NN()`, excluding `seq` and `c_seq`, can be acceptable.
    primer_legnth : tuple of int, optional
        A tuple (min_size, max_size) specifying the primer length. Default is (16, 25).
    design_num : int, optional
        Number of primer pairs to design. Defaults to 1.
    fw_adapter : QUEEN object, str or Cutsite, optional
        If it's a string or a single-stranded DNA (ssDNA) QUEEN object, the sequence will be added at   
        the beginning of any forward primers designed.  
        If it's a double-stranded DNA (dsDNA) QUEEN object with linear topology, the homology sequence   
        to the 3' end of the QUEEN object will be added at the beginning of any forward primers.  
        Alternatively, you can specify the name of a restriction enzyme or 'attB'. In these cases,   
        the adapter sequence including the specified site will be added at the beginning of the primers.  
        For now, "attB" sequence as fw adapter is "GGGGACAAGTTTGTACAAAAAAGCAGGCT".
    rv_adapter : QUEEN object, str or Cutsite, optional
        If it's a string or a single-stranded DNA (ssDNA) QUEEN object, the sequence will be added at   
        the beginning of any reverse primers designed.  
        If it's a double-stranded DNA (dsDNA) QUEEN object with linear topology, the homology sequence   
        to the 3' end of the QUEEN object will be added at the beginning of any reverse primers.  
        Alternatively, you can specify the name of a restriction enzyme or 'attB'. In these cases,   
        the adapter sequence including the specified site will be added at the beginning of the primers.  
        For now, "attB" sequence as rv adapter is "GGGGACCACTTTGTACAAGAAAGCTGGGT".
    cloning_mode : str("gibson", "infusion", "overlappcr"), optional
        If the adapter is a dsDNA QUEEN object, the adapter sequence, which can be joined with the target amplicon   
        via the specified `cloning_mode`, will be automatically generated.
    homology_length : int, optional
        This parameter is active if either `fw_adapter` or `rv_adapter` is specified as a dsDNA QUEEN object 
        and `cloning_mode` is `"gibson"`, `"infusion"`, or `"overlappcr"`.  
        If an int value is provided, an adapter sequence including an overlapping end with the specified  
        dsDNA QUEEN object will be designed, such that the overlap is greater than or equal to the provided value.  
        Default value is 15.
    nonspecific_limit : int, optional
        The maximum number of mismatches allowed for primer binding outside of the designated primer design region  
        within the template sequence. Primer pairs that bind to any region of the template with a number of mismatches  
        equal to or less than this limit will be excluded from the design, to increase the specificity of the PCR reaction  
        and decrease the likelihood of nonspecific amplification.  
        Defaults to 3.
    auto_adjust : bool, optional
        If True and the adapter is dsDNA QUEEN object, the adapter sequence will be automatically adjusted to ensure  
        the reading frame of the genes in the target amplicon. 
    requirement : lambda function, optional
        Function that takes a dictionary representing a primer pair and returns True if the pair meets the specified conditions.
        Default requirement is as follows.
            - `x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T")`
            - `"AAAA" not in x["fw"] and "TTTT" not in x["fw"] and "GGGG" not in x["fw"] and "CCCC" not in x["fw"]`
            - `"AAAA" not in x["rv"] and "TTTT" not in x["rv"] and "GGGG" not in x["rv"] and "CCCC" not in x["rv"]`
    fw_name : str, optional
        The forward primer name designed in this function. If the value is not specified, the foward primer is named as `fw_primer`. 
    rv_name : str, optional 
        The reverse primer name designed in this function. If the value is not specified, the foward primer is named as `rv_primer`.
    
    Raises
    ------
    ValueError
        If the target sequence is not found within the template sequence.

    Returns
    -------
    list of dict
        A list of dictionaries where each dictionary represents a primer pair.
        Each dictionary contains two keys, "fw" and "rv", with the corresponding primer sequences formed by ssDNA QUEEN objects.
        The list is sorted by the closeness of the primer's Tm to the target_tm, with the closest pair first.

        Example return value:
        [
            {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
        ]
    
    Notes
    -----
    It is assumed that template and target are provided as QUEEN objects with appropriate annotations
    for target regions. The requirement for the target sequence to be within the template sequence
    ensures specificity of the primers to the region of interest. The function will not proceed if
    the target sequence is not a subset of the template.
    
    Example
    -------
    >>> from QUEEN.queen import *
    >>> template_Q = QUEEN('ATGC...')
    >>> target_Q = QUEEN('ATGC...')
    >>> # Assuming target_Q sequence is within template_Q
    >>> primers = primerdesign(template_Q, target_Q, target_tm=65.0, design_num=5)
    >>> primers
    [
        {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
        {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
        ...
    ]
    """
    def append_adapter(amplicon_region, filtered_primer_pairs, adapter, mode, adapter_form, strand, name):
        if adapter is None:
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq=filtered_primer_pairs[i][strand][0], ssdna=True, product=name)

        elif adapter == "attB":
            for i in range(len(filtered_primer_pairs)):
                if strand == "fw":
                    filtered_primer_pairs[i][strand][0] = QUEEN(seq="GGGGACAAGTTTGTACAAAAAAGCAGGCT" + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)
                else:
                    filtered_primer_pairs[i][strand][0] = QUEEN(seq="GGGGACCACTTTGTACAAGAAAGCTGGGT" + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)

        elif (type(adapter) == str and adapter in cs.lib.keys()) or "Cutsite" in type(adapter).__name__:
            if type(adapter) == str: 
                adapter = cs.lib[adapter]
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq="ATGC" + adapter.seq + filtered_primer_pairs[i][strand][0], ssdna=True, product=name) 
        
        elif type(adapter) == str:
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq=adapter + filtered_primer_pairs[i][strand][0], ssdna=True, product=name) 

        elif type(adapter) == QUEEN and adapter._ssdna == True:
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq=adapter.seq + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)

        elif type(adapter) == QUEEN and adapter._ssdna == False:
            adapter_features = [feat for feat in adapter.dnafeatures if feat.feature_type != "source"]
            for i in range(len(filtered_primer_pairs)):
                s = filtered_primer_pairs[i]["fw"][1]
                e = len(amplicon_region.seq) - filtered_primer_pairs[i]["rv"][1] 
                pcr_amplicon = amplicon_region[s:e] 
                amplicon_features = [feat for feat in pcr_amplicon.dnafeatures if feat.feature_type != "source"]
                amplicon_features.sort(key=lambda x: x.start) 
                
                if strand == "fw":
                    feat1 = adapter_features[-1] 
                    feat2 = amplicon_features[0] 
                    if feat1.feature_type == "promoter" and feat2.feature_type == "CDS":
                        if feat1.strand == 1 and feat2.strand == -1:
                            print("**Attention**: The directions of the gene in the adapter and the gene in the application are inconsistent. Could you confirm whether the direction of the target amplicon is as intended?")
                        else:
                            pass
                    if feat1.feature_type == "CDS" and feat2.feature_type == "promoter":
                        if feat1.strand == -1 and feat2.strand == 1:
                            print("**Attention**: The directions of the gene in the adapter and the gene in the application are inconsistent. Could you confirm whether the direction of the target amplicon is as intended?")
                        else:
                            pass 

                if strand == "rv":
                    feat1 = amplicon_features[-1] 
                    feat2 = adapter_features[0] 
                    if feat1.feature_type == "CDS" and feat2.feature_type == "promoter":
                        if feat1.strand == 1 and feat2.strand == -1:
                            print("**Attention**: The directions of the gene in the adapter and the gene in the application are inconsistent. Could you confirm whether the direction of the target amplicon is as intended?")
                        else:
                            pass
                    if feat1.feature_type == "promoter" and feat2.feature_type == "CDS":
                        if feat1.strand == -1 and feat2.strand == 1:
                            print("**Attention**: The directions of the gene in the adapter and the gene in the application are inconsistent. Could you confirm whether the direction of the target amplicon is as intended?")
                        else:
                            pass 

                if feat1.feature_type == "CDS" and feat2.feature_type == "CDS":
                    if ("broken_feature" in feat1.qualifiers or "broken_feature" in feat2.qualifiers) and (len(feat1.sequence)%3 != 0 or len(feat2.sequence)%3 != 0):
                        req = False         
                    else:
                        if feat1.strand == feat2.strand:
                            req = True
                        else:
                            req = True
                            print("**Attention**: The directions of the gene in the adapter and the gene in the application are inconsistent. Could you confirm whether the direction of the target amplicon is as intended?")
                else:
                    req = False 

                if strand == "fw": 
                    if mode == "gibson":
                        if adapter._right_end_bottom == 1 and adapter._right_end_top == -1: 
                            adapter = adapter[:len(adapter.seq) - len(adapter._right_end)] 
                        else:
                            pass 
                    elif mode == "infusion":
                        if adapter._right_end_bottom == -1 and adapter._right_end_top == 1:
                            adapter = adapter[:len(adapter.seq) - len(adapter._right_end)] 
                        else:
                            pass 
                if strand == "rv":
                    if mode == "gibson":
                        if adapter._left_end_bottom == -1 and adapter._left_end_top == 1: 
                            adapter = adapter[len(adapter._left_end):] 
                        else:
                            pass 
                    elif mode == "infusion":
                        if adapter._left_end_bottom == 1 and adapter._left_end_top == -1:
                            adapter = adapter[len(adapter._left_end):] 
                        else:
                            pass
                
                if req == True:
                    if type(adapter_form) == int:
                        if strand == "fw":
                            feat1 = adapter_features[-1] 
                            feat2 = amplicon_features[0] 
                            fragment1 = adapter[feat1.start:].seq  
                            fragment2 = pcr_amplicon[:feat2.end].seq 
                            gapseq = "".join([random.choice("ATGC") for _ in range((len(fragment1) + len(fragment2)) % 3)])
                            filtered_primer_pairs[i][strand][0] = QUEEN(seq=adapter.seq[-1*adapter_form:] + gapseq + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)
                        else:
                            feat1 = amplicon_features[-1] 
                            feat2 = adapter_features[0] 
                            fragment1 = pcr_amplicon[feat1.start:].seq  
                            fragment2 = adapter[:feat2.end].seq 
                            gapseq = "".join([random.choice("ATGC") for _ in range((len(fragment1) + len(fragment2)) % 3)])
                            filtered_primer_pairs[i][strand][0] = QUEEN(seq=adapter.rcseq[-1*adapter_form:] + gapseq + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)

                    elif type(adapter_form) == Cutsite:
                        pass

                else:
                    if type(adapter_form) == int:
                        if strand == "fw":
                            filtered_primer_pairs[i][strand][0] = QUEEN(seq=adapter.seq[-1*adapter_form:] + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)
                        else:
                            filtered_primer_pairs[i][strand][0] = QUEEN(seq=adapter.rcseq[-1*adapter_form:] + filtered_primer_pairs[i][strand][0], ssdna=True, product=name)
                    elif type(adapter_form) == Cutsite:
                        pass 
        else:
            raise TypeError("adapter object must be instance of QUEEN, Cutsite or str class.")
     
        return filtered_primer_pairs
    
    if requirement is None:
        def requirement(x):
            req1 =  x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T") 
            req4 = "AAAA" not in x["fw"] and "TTTT" not in x["fw"] and "GGGG" not in x["fw"] and "CCCC" not in x["fw"]
            req5 = "AAAA" not in x["rv"] and "TTTT" not in x["rv"] and "GGGG" not in x["rv"] and "CCCC" not in x["rv"]
            return req1 * req4 *  req5
    
    if type(homology_length) == int or homology_length is None:
        mode = "gibson"
        adapter_form = homology_length
    
    elif type(homology_length) == tuple:
        mode = homology_length[0] 
        adapter_form = homology_length[1] 

    if type(template) != QUEEN:
        raise TypeError("`template` object must be instance of QUEEN class.") 
    
    if type(target) != QUEEN:
        raise TypeError("`template` object must be instance of QUEEN class.") 

    if target.seq not in template.seq:
        raise ValueError("target sequence to be amplified is not included in template sequence.") 
    
    if fw_primer is not None:
        if type(fw_primer) == QUEEN:
            pass 
        elif type(fw_primer) == str:
            fw_primer = QUEEN(seq=fw_primer, ssdna=True) 
        else:
            raise TypeError("`fw_primer` object must be instance of ssDNA QUEEN class.")
    
    if rv_primer is not None:
        if type(rv_primer) == QUEEN:
            pass 
        elif type(rv_primer) == str:
            rv_primer = QUEEN(seq=rv_primer, ssdna=True) 
        else:
            raise TypeError("`rv_primer` object must be instance of ssDNA QUEEN class.")

    if tm_func is None:
        tm_func = Tm_NN() 
    elif tm_func == "SantaLucia":
        tm_func = Tm_NN(nn_table=mt.DNA_NN3)
    elif tm_func == "Breslauer":
        tm_func = Tm_NN(nn_table=mt.DNA_NN1)

    start = template.seq.find(target.seq) - fw_margin
    if start < 0:
        if template.topology == "circular":
            start = len(template.seq) + start
        else:
            start = 0
    
    end = template.seq.find(target.seq) + len(target.seq) + rv_margin
    if end > len(template.seq):
        if template.topology == "circular": 
            end = end - len(template.seq) * (end // len(template.seq))  
        else:
            pass
    amplicon_region = template[start:end]
    start = amplicon_region.seq.find(target.seq) 
    end   = start + len(target.seq)
    fw_candidates = [] 
    if fw_primer is None:
        for pos in range(start+1):
            for plen in range(primer_length[0], primer_length[1] + 1): 
                fw_candidate = amplicon_region.seq[pos:pos+plen]
                fw_candidates.append([str(fw_candidate), pos]) 
    else:
        site = amplicon_region.searchsequence(query=fw_primer, quinable=False)
        fw_candidates.append([fw_primer.seq, site.start])

    rv_candidates = [] 
    if rv_primer is None:
        for pos in range(end-len(target.seq)+1):
            for plen in range(primer_length[0], primer_length[1] + 1): 
                rv_candidate = amplicon_region.rcseq[pos:pos+plen]
                rv_candidates.append([str(rv_candidate), pos]) 
    else:
        site = amplicon_region.searchsequence(query=rv_primer, quinable=False)
        rv_candidates.append([rv_primer.seq, len(amplicon_region.seq) - site.end])

    checked_fw_candidates = [] 
    for candidate in fw_candidates:
        sites = template.searchsequence(query="(?:{}){{s<={}}}".format(candidate[0], nonspecific_limit), quinable=False) 
        if len(sites) > 1:
            pass 
        else:
            checked_fw_candidates.append(candidate) 
            
    checked_rv_candidates = [] 
    for candidate in rv_candidates:
        sites = template.searchsequence(query="(?:{}){{s<={}}}".format(candidate[0], nonspecific_limit), quinable=False) 
        if len(sites) > 1:
            pass
        else:
            checked_rv_candidates.append(candidate) 

    fw_tm_set = []
    for candidate in fw_candidates:
        tm = tm_func(seq=candidate[0])
        fw_tm_set.append([candidate, tm]) 
    
    rv_tm_set = [] 
    for candidate in rv_candidates: 
        tm = tm_func(seq=candidate[0])
        rv_tm_set.append([candidate, tm]) 
    
    primer_pairs = [] 
    for fw, rv in it.product(fw_tm_set, rv_tm_set):
        primer_pairs.append({"fw":copy.deepcopy(fw[0]), "rv":copy.deepcopy(rv[0]), "fw_tm":fw[1], "rv_tm":rv[1]}) 
    primer_pairs.sort(key=lambda x: abs(x["fw_tm"]-target_tm) + abs(x["rv_tm"]-target_tm))

    filtered_primer_pairs = [] 
    for primer_pair in primer_pairs:
        if requirement(primer_pair): 
            filtered_primer_pairs.append(primer_pair)  
        else:
            pass
    filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, fw_adapter, mode, adapter_form, "fw", fw_name)
    filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, rv_adapter, mode, adapter_form, "rv", rv_name)

    for i in range(len(filtered_primer_pairs)):
        filtered_primer_pairs[i]["fw"] = filtered_primer_pairs[i]["fw"][0] 
        filtered_primer_pairs[i]["fw"].setfeature({"feature_type":"primer_bind", "qualifier:label":fw_name})
        filtered_primer_pairs[i]["rv"] = filtered_primer_pairs[i]["rv"][0]
        filtered_primer_pairs[i]["rv"].setfeature({"feature_type":"primer_bind", "qualifier:label":rv_name})
    return filtered_primer_pairs[:design_num]
