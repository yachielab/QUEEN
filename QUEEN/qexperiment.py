import copy
import itertools as it
from qfunction import joindna, cropdna, cutdna, flipdna, modifyends
from qobj import QUEEN 
from qseq import Qseq 
import cutsite as cs
from cutsite import Cutsite
from Bio.SeqUtils import MeltingTemp as mt
import functools

def pcr(template, fw, rv, bindnum=16, mismatch=1, endlength=3, return_template=False, product=None, process_name=None, process_description=None, pn=None, pd=None, **args):
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
        The minimum number of binding nucleotides for a primer, by default 16.
    mismatch : int, optional
        The maximum number of mismatches allowed in the primer binding, by default 1.
    endlength : int, optional
        The length of the end region of the primer to consider during binding, by default 3.
    return_template : bool, optional
        If True, returns the modified template, by default False.
    product : str, optional 
        Product name of the PCR process.
    process_name : str, optional
        Brief label for the PCR process, by default "PCR"
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the PCR process.
    pd : str, optional
        Alias for process_description.
    **args
        Additional keyword arguments for advanced configurations.

    Returns
    -------
    QUEEN (amplicon) or QUEEN (template) and QUEEN (amplicon)
        If return_template is True, returns the modified DNA template with primer binding features. 
        Otherwise, returns the PCR product (amplicon) as a QUEEN object.

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
    def search_binding_site(template, primer, strand=1, endlength=3, pn=None, pd=None): 
        primer_end = primer.seq[-1*endlength:]
        for i in range(bindnum-endlength, len(primer.seq)-endlength):
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
            return site[0]
        
        elif len(site) == 0:
            raise ValueError("No primer binding sites were found.")
        
        elif len(site) > 1: 
            raise ValueError("Multiple primer binding sites were detected.") 
        
        #else:
        #    raise ValueError("Primer binded to an unexpected strand.") 
    
    if type(template) != QUEEN: 
        raise TypeError("`template` object must be instance of QUEEN class.") 

    if type(fw) == Qseq or type(fw) == str:
        fw = QUEEN(seq=fw, ssdna=True) 
    elif type(fw) == QUEEN:
        pass 
    else:
        raise TypeError("`fw` object must be instance of QUEEN or str class.") 

    if type(rv) == Qseq or type(rv) == str:
        rv = QUEEN(seq=rv, ssdna=True) 
    elif type(rv) == QUEEN:
        pass 
    else:
        raise TypeError("`fw` object must be instance of QUEEN or str class.") 
    
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "PCR"
    
    pd_suffix = "pcr({}, {}, {}, bindnum={}, mismatch={}, endlegth={}, {})".format(template.project, fw.project, rv.project, bindnum, mismatch, endlength, ", ".join([str(key) + "=" + str(args[key]) for key in args]))
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description
    if -1 in [template._left_end_top, template._left_end_bottom, template._right_end_top, template._right_end_bottom]:
        template = modifyends(template, pn=process_name, pd=process_description)
    
    site1 = search_binding_site(template, fw, 1, endlength, pn=process_name, pd=process_description) 
    site2 = search_binding_site(template, rv, -1, endlength, pn=process_name, pd=process_description) 
    if site1.strand == 1 and site2.strand == -1:
        fw_site = site1  
        rv_site = site2
    elif site1.strand == -1 and site2.strand == 1:
        fw_site = site2
        rv_site = site1
    else:
        raise ValueError("Both primers binded to the same strand.") 

    fw_bind_length = len(fw_site.sequence) 
    rv_bind_length = len(rv_site.sequence) 
    
    #fw.setfeature({"start":len(fw.seq)-fw_bind_length, "end":len(fw.seq), "qualifier:note":"primer_bind"})  
    #rv.setfeature({"start":len(rv.seq)-rv_bind_length, "end":len(rv.seq), "qualifier:note":"primer_bind"})
    
    extract  = cropdna(template, fw_site.end, rv_site.start, pn=process_name, pd=process_description)
    amplicon = modifyends(extract, fw.seq, rv.rcseq, pn=process_name, pd=process_description, product=product)
    
    if return_template == True:
        template.setfeature({"start": fw_site.start, "end": fw_site.end, "strand":1,  "feature_type":"primer_bind"}) 
        template.setfeature({"start": rv_site.start, "end": rv_site.end, "strand":-1, "feature_type":"primer_bind"})  
        return template, amplicon
    else:
        return amplicon

def _select(fragments, selection=None): 
    if selection is None:
        return fragments  
    
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

def digestion(dna, *cutsites, selection=None, requirement=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **args):
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
    **args
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
        if (type(selection) != tuple) and selection not in ("min", "max") and (selection.startswith("label:") == False and selection.startswith("!label:") == False):
            raise TypeError("`selection` should be `tuple` value, 'min', 'max', or `str` starting with 'label:' or '!label'.")
    
    cutsite_names = [] 
    for c in range(len(cutsites)):
        if type(cutsites[c]) == Cutsite or "cutsite" in cutsites[c].__dict__:
            pass 
        elif type(cutsites[c]) == str and cutsites[c] in cs.lib.keys():
            cutsites[c] = cs.lib[cutsites[c]]  
        else:
            raise TypeError("Each element in `cutsites` must be instance of Cutsite class or its name must be included in `QUEEN.cutsite.lib`.")
        cutsite_names.append(cutsites[c].name) 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Digestion"
    
    if type(selection) == tuple:
        pd_suffix = "#digestion({}, cutsites=*[{}], selection=[{}], {})".format(dna.project, ",".join(cutsite_names), ",".join(map(str, selection)), ", ".join([str(key) + "=" + str(args[key]) for key in args])) 
    else:
        pd_suffix = "#digestion({}, cutsites=*[{}], selection={}, {})".format(dna.project, ",".join(cutsite_names), selection, ", ".join([str(key) + "=" + str(args[key]) for key in args])) 
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description

    new_cutsites = [] 
    for cutsite in cutsites:
        sites = dna.searchsequence(query=cutsite, pn=process_name, pd=process_description)        
        new_cutsites.extend(sites) 
    
    fragments = cutdna(dna, *new_cutsites, product=product, pn=process_name, pd=process_description)
    
    if requirement is not None:
        fragments = [fragment for fragment in fragments if requirement] 
    else:
        pass 
    
    return _select(fragments, selection) 

def ligation(*fragments, unique=True, product=None, process_name=None, process_description=None, pn=None, pd=None, **args): 
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
    **args
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

    pd_suffix = "#ligation({}, unique={}, {})".format(", ".join([fragment.project for fragment in fragments]), unique, ", ".join([str(key) + "=" + str(args[key]) for key in args]))
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description #"\n" + pd_suffix
    nums = list(range(len(fragments)))
    nums_orders = list(map(list,it.permutations(nums[:-1])))
    nums_orders = [numlist + [nums[-1]] for numlist in nums_orders]
    flip_status_list = list(it.product(*[[1,-1] for i in range(len(fragments))]))   
    products = [] 
    
    for numset in nums_orders:
        execed = [] 
        for flipset in flip_status_list:
            if tuple([state * -1 for state in flipset]) in execed: 
                pass 
            else:
                fragment_set = [fragments[num] if flip == 1 else flipdna(fragments[num], product=fragments[num].project, pn=process_name, pd=process_description) for num, flip in zip(numset, flipset)] 
                try:
                    outobj = joindna(*fragment_set, topology="circular", autoflip=False, compatibility="complete", pn=process_name, pd=process_description, product=product) 
                    products.append(outobj) 
                except Exception as e:
                    pass 
            execed.append(flipset) 

    if unique == True:
        if len(products) == 0:
            raise ValueError("The QUEEN_objects cannot be joined due to the end structure incompatibility.") 
        elif len(products) > 1: 
            raise ValueError("Multiple assembled constructs were detected. You should review your assembly design.")
        else:
            return products[0] 
    else: 
        return products 

def homology_based_assembly(*fragments, mode="gibson", homology_length=20, unique=True, product=None, process_name=None, process_description=None, pn=None, pd=None, **args): #homology_based_assembly
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
    product : str, optional 
        Product name of the homology_based_assembly (hba) process. 
    process_name : str, optional
        Brief label for the `homology_based_assembly` process .
        If `mode` is `"gibson"`, default is "Gibson Assembly". 
        If `mode` is `"infusion"`, default is "In-Fusion Assembly". 
        If `mode` is `"overlappcr"`, default is "Overlap PCR".
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the assembly process.
    pd : str, optional
        Alias for process_description.
    **args
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
    max_homology_length = 200 #max_homology_length
    
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

    pd_suffix = "#homology_based_assembly({}, mode={}, homology_length={}, unique={}, {})".format(", ".join([fragment.project for fragment in fragments]), mode, homology_length, unique, ", ".join([str(key) + "=" + str(args[key]) for key in args]))
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description
    
    for fragment in fragments:
        if fragment.topology != "linear": 
            raise ValueError("A 'circular' fragment was detected. All fragments to be assembled should be 'linear' topology.")

    if mode not in ("gibson", "infusion", "overlappcr"):
        raise ValueError("Invalid mode value. The 'mode' variable can only take 'gibson', 'infusion' or 'overlappcr' as values.")
    
    if mode == "overlappcr":
        nums = list(range(len(fragments)))
        nums_orders = list(it.permutations(nums)) 
        flip_status_list = list(it.product(*[[1,-1] for i in range(len(fragments))]))   
        products = [] 
        for numset in nums_orders:
            execed = [] 
            for flipset in flip_status_list:
                if tuple([state * -1 for state in flipset]) in execed: 
                    pass 
                else:
                    fragment_set = [fragments[num] if flip == 1 else flipdna(fragments[num]) for num, flip in zip(numset, flipset)] 
                    try:
                        fragment1 = fragment_set[0]
                        for f, fragment2 in enumerate(fragment_set[1:]): 
                            fragment1 = modifyends(fragment1, "*{}/-{}".format(len(fragment1.seq), len(fragment1.seq)), "", pn=process_name, pd=process_description)
                            fragment2 = modifyends(fragment2, "-{}/*{}".format(len(fragment2.seq), len(fragment2.seq)), "", pn=process_name, pd=process_description)
                            fragment1 = joindna(fragment1, fragment2, topology="linear", homology_length=homology_length, pn=process_name, pd=process_description) 
                            if f == len(fragment_set) - 2:
                                fragment1 = modifyends(fragment1, "*/*", "*/*", pn=process_name, pd=process_description, product=product) 
                            else:
                                fragment1 = modifyends(fragment1, "*/*", "*/*", pn=process_name, pd=process_description) 
                        products.append(fragment1)
                    except:
                        pass 
                execed.append(flipset) 

    else: 
        nums = list(range(len(fragments)))
        nums_orders = list(map(list,it.permutations(nums[:-1])))
        nums_orders = [numlist + [nums[-1]] for numlist in nums_orders]
        flip_status_list = list(it.product(*[[1,-1] for i in range(len(fragments))]))   
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
                        if len(fragment.seq) <= 2 * max_homology_length: 
                            mhl = int(len(fragment.seq)/2) 
                        else:
                            mhl = max_homology_length
                        if mode == "gibson":
                            fragment_set[f] = modifyends(fragment, "-{{{}}}/*{{{}}}".format(mhl,mhl), "*{{{}}}/-{{{}}}".format(mhl,mhl), pn=process_name, pd=process_description)
                        elif mode == "infusion":
                            fragment_set[f] = modifyends(fragment, "*{{{}}}/-{{{}}}".format(mhl,mhl), "-{{{}}}/*{{{}}}".format(mhl,mhl), pn=process_name, pd=process_description)
                    try:
                        outobj = joindna(*fragment_set, autoflip=False, homology_length=homology_length, topology="circular", pn=process_name, pd=process_description, product=product) 
                        products.append(outobj) 
                    except Exception as e:
                        #fragment_set[0].printsequence(display=True, hide_middle=30) 
                        #fragment_set[1].printsequence(display=True, hide_middle=30)
                        pass 
                execed.append(flipset) 

    if unique == True:
        if len(products) > 1: 
            raise ValueError("Multiple assembled constructs were detected. You should review your assembly design.")
        else:
            return products[0] 
    else: 
        return products 

def annealing(ssdna1, ssdna2, homology_length=4, product=None, pn=None, pd=None, process_name=None, process_description=None, **args):
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
        Brief label for the ligation process, by default "Annealing".
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the annealing process.
    pd : str, optional
        Alias for process_description.
    **args
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
    pd_suffix = "annealing({}, {}, homology_length={})".format(ssdna1.project, ssdna2.project, homology_length)
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix

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
    
    annealed_dna  = joindna(ssdna1, ssdna2, homology_length=homology_length, pn=process_name, pd=process_description, product=product)
    
    ssdna1._ssdna = False if flag1 == 1 else True 
    ssdna2._ssdna = False if flag2 == 1 else True

    return annealed_dna 

def gateway_reaction(detination, entry, mode="BP", destination_selection="max", entry_selection="min", process_name=None, process_description=None, pn=None, pd=None, **args):
    """
    Simulates a gateway reaction of two DNA molecules.
    Basic `BP` and `LR` reactions are available.

    Parameters
    ----------
    destination : QUEEN
        The destination QUEEN object holding the backbone DNA molecule.
    entry : QUEEN
        The entry QUEEN object holding the insert DNA molecule.
    mode: str, tuple, or list
        The mode of the reaction, can be "BP" or "LR". Default is "BP".
    destination_selection : "min", "max", "label:*", "!label:*", or tuple of int, optional
        A rule to select a specific fragment from the destination DNA with the specified condition. Default is "max".
        The specification format is exactly the same as the `selection` parameter in the `digestion` function.
        For more details, refer to its description.
    entry_selection : "min", "max", "label:*", "!label:*", or tuple of int, optional
        A rule to select a specific fragment from the entry DNA with the specified condition. 
        If the entry's topology is circualr, default is "max". Otherwise "min".
        The specification format is exactly the same as the `selection` parameter in the `digestion` function.
        For more details, refer to its description.
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Gateway Reaction".
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the gateway reaction process.
    pd : str, optional
        Alias for process_description.
    **args
        Additional keyword arguments for advanced configurations.

    Returns
    ----------
    QUEEN
        The QUEEN object representing the result of the gateway reaction process.
    """
    if destination_selection is None:
        destination_selection = "max"

    if entry_selection is None:
        if entry.topology == "linear":
            entry_selection = "max"
        else:
            entry_selection = "min" 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Gateway Reaction" 
    pd_suffix = "gateway_reaction({}, {}, mode={}, destination_selection={}, entry_selection={})".format(destination.project, entry.project, mode, destination_selection, entry_selection)
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix

    if mode == "BP":
        cs.lib["attX1"] = "CAAGTTT^GTACAAA_AAAGCAGGCT"  #attB1
        cs.lib["attX2"] = "ACCCAGCTTT^CTTGTAC_AAAGTGG"  #attB2
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
    
    entry_fragments = cutdna(entry, attx1, attx2, pn=process_name, pd=process_description) 
    entry_fragments.sort(key=lambda x:len(x.seq))
    insert = _select(entry_fragments, entry_selection) 

    destination_fragments = cutdna(destination, atty1, atty2, pn=process_name, pd=process_description)      
    destination_fragments._fragments.sort(key=lambda x:len(x.seq)) 
    destination = _select(destination_fragments, destination_selection)

    return ligation(insert, destination, product=product, pn=process_name, pd=process_description) 

def golden_gate_assembly(destination, entry, enzyme=None, destination_selection="max", entry_selection="min", process_name=None, process_description=None, pn=None, pd=None, **args):
    """
    Simulates a Golden Gate Assembly.

    Parameters
    ----------
    destination : QUEEN
        The destination QUEEN object holding the backbone DNA molecule.
    entry : QUEEN object or list of QUEEN objects
        The entry QUEEN object(s) holding the insert DNA molecules.
    enzyme : Cutsite or str
        The restriction enzyme used for this reaction.   
    destination_selection : "min", "max", "label:*", "!label:*", or tuple of int, optional  
        A rule to select a specific fragment from the destination DNA with the specified condition. Default is "max".  
        The specification format is exactly the same as the `selection` parameter in the `digestion` function.  
        For more details, refer to its description.  
    entry_selection : "min", "max", "label:*", "!label:*", tuple of int, or list of these value, optional
        A rule to select a specific fragment from each entry DNA with the specified condition. 
        If the entry's topology is circualr, default is "max". Otherwise "min".
        For applying the independent rule to each fragment, the value should be given by `list`. 
        The specification format is exactly the same as the `selection` parameter in the `digestion` function.  
        For more details, refer to its description.  
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Golden Gate Assembly".
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the gateway reaction process.
    pd : str, optional
        Alias for process_description.
    **args
        Additional keyword arguments for advanced configurations.
    """ 
    if destination_selection is None:
        destination_selection = "max"

    if entry_selection is None:
        entry_selection_list = [] 
        for aentry in entry:
            if aentry.topology == "linear":
                selection = "max"
            else:
                selection = "min" 
            entry_selection_list.append(selection) 
    else:
        if type(entry_selection) == list:
            pass 
        else:
            entry_selection_list = [entry_selection] * len(entry)
    
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Golden Gate Assembly" 
    pd_suffix = "golden_gate_assembly({}, {}, enzyme={}, destination_selection={}, entry_detection={})".format(destination.project, "[" + ",".join([aentry.project for aentry in entry]) + "]", enzyme, destination_selection, entry_selection)
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix

    if type(entry) == QUEEN:
        entry = [entry]
    
    fragments = [] 
    for aentry, selection in zip(entry, entry_selection_list):
        insert = digestion(aentry, enzyme, selection=selection, pn=process_name, pd=process_description) 
        fragments.append(insert) 

    backbone =  digestion(destination, enzyme, selection=destination_selection, pn=process_name, pd=process_description) 
    fragments.append(backbone) 

    return ligation(*fragments, pn=process_name, pd=process_description) 

def intra_site_specific_recombination(dna, site="loxP", process_name=None, process_description=None, pn=None, pd=None, **args):
    """
    Simulates a intra molecule site-specific recombination.

    Parameters
    ----------
    dna : QUEEN
        The target QUEEN object.
    site : str ("loxP", "lox2272", "FRT") 
        The target site for the recombination reaction. At least two identical   
        recombination sites in the DNA sequence to simulate the recombination process.
    process_name : str, optional
        Brief label for the gateway reaction process. Default is "Golden Gate Assembly".
    pn : str, optional
        Alias for process_description.
    process_description : str, optional
        Additional description for the gateway reaction process.
    pd : str, optional
        Alias for process_description.
    **args
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
    pd_suffix = "intra_site_specific_recombination({}, site={})".format(dna.project, site)
    process_description = pd if process_description is None else process_description
    process_description = pd_suffix if process_description is None else process_description #+ "\n" + pd_suffix
    
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
                outobj = joindna(fragment, autoflip=False, compatibility="complete", topology="circular", pn=process_name, pd=process_description)
            else:
                outobj = joindna(fragment[0], fragment[2], autoflip=False, compatibility="complete", pn=process_name, pd=process_description) 

        else:
            if dna.topology == "circular":
                fragment1 = _select(fragments, selection="max") 
                fragment2 = _select(fragments, selection="min") 
                fragment2 = flipdna(fragment2, pn=process_name, pd=process_description)
                outobj = joindna(fragment1, fragment2, autoflip=False, compatibility="complete", topology="circular", pn=process_name, pd=process_description)
            else:
                reversed_fragment = flipdna(fragment[1], pn=process_name, pd=process_description) 
                outobj = joindna(fragment[0], reversed_fragment, fragment[2], autoflip=False, compatibility="complete", pn=process_name, pd=process_description) 
        outobj_list.append(outobj)
    
    if len(outobj_list) == 1:
        return outobj_list[0]
    else:
        return outobj_list 


def Tm_NN(check=True, strict=True, nn_table=None, tmm_table=None, imm_table=None, de_table=None, dnac1=25, dnac2=25, selfcomp=False, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5): 
    return functools.partial(mt.Tm_NN, check=check, strict=strict, nn_table=nn_table, tmm_table=tmm_table, imm_table=imm_table, de_table=de_table, dnac1=dnac1, dnac2=dnac2, selfcomp=selfcomp, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)  

def primerdesign(template, target, fw_primer=None, rv_primer=None, fw_margin=0, rv_margin=0,
                 target_tm=60.0, tm_func=None, primer_length=(16, 25), 
                 design_num=1, fw_adapter=None, rv_adapter=None, homology_length=20, nonspecific_limit=3, 
                 requirement=lambda x: x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T"),
                 fw_name="fw_primer", rv_name="rv_primer"):
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
    tm_func : function, optional
        Function to calculate the melting temperature of primer candidates. Default is xxxx. 
        As built-in algorithms, `QUEEN.qexperiment.Tm_NN()`. This function is implemented 
        based on the `Bio.SeqUtils.MeltingTemp.Tm_NN()`, so the all parameters of 
        `Bio.SeqUtils.MeltingTemp.Tm_NN()`, excluding `seq` and `c_seq`, can be acceptable.
    primer_legnth : tuple of int, optional
        A tuple (min_size, max_size) specifying the primer length. Default is (16, 25).
    design_num : int, optional
        Number of primer pairs to design. Defaults to 1.
    fw_adapter : ssDNA QUEEN object or str, optional
        Adapter sequence to prepend to any designed forward primer.
    rv_adapter : ssDNA QUEEN object or str, optional
        Adapter sequence to append to any designed reverse primer.
    homology_length : int, optional
        This parameter is active if either fw_adapter or rv_adapter is specified as a dsDNA QUEEN object.
        If an int value is provided, an adapter sequence including an overlapping end with the specified dsDNA QUEEN object
        will be designed, such that the overlap is greater than or equal to the provided value. Default value is 15.
    nonspecific_limit : int, optional
        The maximum number of mismatches allowed for primer binding outside of the designated primer design region within the template sequence.
        Primer pairs that bind to any region of the template with a number of mismatches equal to or less than this limit will be excluded from the design,
        to increase the specificity of the PCR reaction and decrease the likelihood of nonspecific amplification. Defaults to 3.
    requirement : lambda function, optional
        Function that takes a dictionary representing a primer pair and returns True if the pair meets the specified conditions.
        Ensures that the 3' end nucleotide of both primers is not A or T by default.
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
    >>> primers = primerdesign(template_Q, target_Q, target_tm=65.0, num_design=5)
    >>> primers
    [
        {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
        {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
        ...
    ]
    """
    def append_adapter(amplicon_region, filtered_primer_pairs, adapter, adapter_form, strand, name):
        if type(adapter) == str:
            adapter = QUEEN(seq=adapter, ssdna=True)
        
        if adapter is None:
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq=filtered_primer_pairs[i][strand][0], ssdna=True, product=name)
        
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
                
                fw_req = (adapter_features[-1].feature_type == "CDS" and amplicon_features[0].feature_type == "CDS")                
                rv_req = (adapter_features[0].feature_type == "CDS" and amplicon_features[-1].feature_type == "CDS") 
                req = {"fw":fw_req, "rv":rv_req} 
                
                if req[strand] == True: 
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
            raise TypeError("adapter object must be instance of QUEEN or str class.")
     
        return filtered_primer_pairs
    
    adapter_form = homology_length
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
    
    if fw_adapter is not None: 
        if type(fw_adapter) == QUEEN:
            pass 
        elif type(fw_adapter) == str:
            fw_adapter = QUEEN(seq=fw_adapter, ssdna=True)
        else:
            pass 

    if rv_adapter is not None: 
        if type(rv_adapter) == QUEEN:
            pass 
        elif type(fw_adapter) == str:
            rv_adapter = QUEEN(seq=rv_adapter, ssdna=True)
        else:
            pass

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
            for plen in range(*primer_length): 
                fw_candidate = amplicon_region.seq[pos:pos+plen]
                fw_candidates.append([str(fw_candidate), pos]) 
    else:
        site = amplicon_region.searchsequence(query=fw_primer, quinable=False)
        fw_candidates.append([fw_primer.seq, site.start])

    rv_candidates = [] 
    if rv_primer is None:
        for pos in range(end-len(target.seq)+1):
            for plen in range(*primer_length): 
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
    filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, fw_adapter, adapter_form, "fw", fw_name)
    filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, rv_adapter, adapter_form, "rv", rv_name)

    for i in range(len(filtered_primer_pairs)):
        filtered_primer_pairs[i]["fw"] = filtered_primer_pairs[i]["fw"][0] 
        filtered_primer_pairs[i]["fw"].setfeature({"feature_type":"misc_feature", "qualifier:label":fw_name})
        filtered_primer_pairs[i]["rv"] = filtered_primer_pairs[i]["rv"][0]
        filtered_primer_pairs[i]["rv"].setfeature({"feature_type":"misc_feature", "qualifier:label":rv_name})
    return filtered_primer_pairs[:design_num]
