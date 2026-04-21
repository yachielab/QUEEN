import copy
import regex as re
import itertools as it
from qfunction import joindna, cropdna, cutdna, flipdna, modifyends, editfeature, removeattribute, compile_cutsite 
from qobj import QUEEN 
from qseq import Qseq 
from quine import quine
import cutsite as cs
from cutsite import Cutsite
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Align import PairwiseAligner
import functools
import collections

HF_enzymes = ["AgeI", "ApoI", "BamHI", "BbsI", "BclI", "BmtI", "BsaI", "BsiWI", "BsrGI", "BstEII", "DraIII", 
              "EagI", "EcoRI", "EcoRV", "HindIII", "KpnI", "MfeI", "MluI", "NcoI", "NheI", "NotI", "NruI", 
              "NsiI", "PstI", "PvuI", "PvuII", "SacI", "SalI", "ScaI", "SpeI", "SphI", "SspI", "StyI"]

def _combine_history(dna, histories):
    combined_history = collections.defaultdict(dict) 
    combined_history["building_history"] = {} 
    for history in histories:
        for key in history["building_history"]: 
            combined_history["building_history"][key] = history["building_history"][key] 
    return combined_history 

def _convert_kwargs(arguments):
    out = [] 
    defaults = ["_sourcefile", "process_id", "original_ids"] 
    for key in arguments:
        if key in defaults:
            pass 
        else:
            if type(arguments[key]) == bool: 
                out.append('{}={}'.format(key, arguments[key]))
            else:
                arguments[key] = str(arguments[key]).replace('"', '\"')    
                out.append('{}="{}"'.format(key, arguments[key]))
    if len(out) == 0:
        out = ""
    else:
        out = ", " + ", ".join(out) 
    return out


def _deterministic_gap_seq(fragment1, fragment2, remseq, gap_len):
    if gap_len <= 0:
        return ""
    fragment1 = str(fragment1).upper()
    fragment2 = str(fragment2).upper()
    remseq = str(remseq).upper()
    stop_codons = {"TAA", "TAG", "TGA"}
    best_gap = None
    best_key = None
    for gap_tuple in it.product("ACGT", repeat=gap_len):
        gap = "".join(gap_tuple)
        junction = fragment1 + gap + remseq + fragment2
        stop_count = 0
        for idx in range(0, len(junction) - 2, 3):
            if junction[idx:idx+3] in stop_codons:
                stop_count += 1
        gc_count = gap.count("G") + gap.count("C")
        gc_distance = abs((gap_len / 2.0) - gc_count)
        key = (stop_count, gc_distance, gap)
        if best_key is None or key < best_key:
            best_key = key
            best_gap = gap
    return best_gap or ("A" * gap_len)


def _primer_pair_sort_key(pair, target_tm):
    if target_tm is None:
        tm_delta = 0.0
    else:
        tm_delta = abs(pair["fw_tm"] - target_tm) + abs(pair["rv_tm"] - target_tm)
    fw = pair["fw"]
    rv = pair["rv"]
    fw_seq = str(fw[0]) if type(fw) in (tuple, list) and len(fw) > 0 else str(fw)
    rv_seq = str(rv[0]) if type(rv) in (tuple, list) and len(rv) > 0 else str(rv)
    fw_pos = int(fw[1]) if type(fw) in (tuple, list) and len(fw) > 1 else -1
    rv_pos = int(rv[1]) if type(rv) in (tuple, list) and len(rv) > 1 else -1
    return (tm_delta, fw_pos, rv_pos, fw_seq, rv_seq)


def _oriented_feature_seq(dna, feat):
    strand = feat.strand if feat.strand not in (None, 0) else 1
    return str(dna.printsequence(feat.start, feat.end, strand=strand, display=False))


def _feature_parts_sorted(feat):
    parts = [(int(part.start), int(part.end)) for part in feat.location.parts]
    return sorted(parts, key=lambda x: (x[0], x[1]))


def _feature_contains(container_feat, inner_feat):
    container_parts = _feature_parts_sorted(container_feat)
    inner_parts = _feature_parts_sorted(inner_feat)
    for is_, ie in inner_parts:
        ok = False
        for cs, ce in container_parts:
            if cs <= is_ and ie <= ce:
                ok = True
                break
        if ok is False:
            return False
    return True


def _drop_redundant_subfeatures(product, rescued_feat):
    if rescued_feat.feature_type in ("source", "primer", "primer_bind"):
        return
    rescued_label = str((rescued_feat.qualifiers.get("label") or [""])[0])
    rescued_type = rescued_feat.feature_type
    rescued_strand = rescued_feat.strand if rescued_feat.strand not in (None, 0) else 1
    rescued_seq = _oriented_feature_seq(product, rescued_feat)
    kept = []
    for feat in product.dnafeatures:
        if feat is rescued_feat:
            kept.append(feat)
            continue
        label = str((feat.qualifiers.get("label") or [""])[0])
        feat_type = feat.feature_type
        strand = feat.strand if feat.strand not in (None, 0) else 1
        if label != rescued_label or feat_type != rescued_type or strand != rescued_strand:
            kept.append(feat)
            continue
        if _feature_contains(rescued_feat, feat) is False:
            kept.append(feat)
            continue
        feat_seq = _oriented_feature_seq(product, feat)
        if feat_seq == rescued_seq:
            continue
        if feat_seq != "" and feat_seq in rescued_seq:
            continue
        kept.append(feat)
    product._dnafeatures = kept


def _feature_equivalence_key(dna, feat):
    label = str((feat.qualifiers.get("label") or [""])[0])
    feat_type = feat.feature_type
    strand = feat.strand if feat.strand not in (None, 0) else 1
    source_seq = _oriented_feature_seq(dna, feat)
    return (label, feat_type, strand, source_seq)


def _count_equivalent_features(product, feature_key):
    count = 0
    for prod_feat in product.dnafeatures:
        if _feature_equivalence_key(product, prod_feat) == feature_key:
            count += 1
    return count


def _source_feature_requirements(source_dnas):
    requirements = collections.Counter()
    for source_dna in source_dnas:
        seen_keys = set()
        for feat in source_dna.dnafeatures:
            if feat.feature_type in ("source", "primer", "primer_bind"):
                continue
            feature_key = _feature_equivalence_key(source_dna, feat)
            # crop/modifyends can leave duplicate feature records that are
            # identical in label/type/strand/sequence. Count each unique
            # feature key only once per source molecule so the rescue loop
            # does not keep "re-adding" an already satisfied feature forever.
            if feature_key in seen_keys:
                continue
            seen_keys.add(feature_key)
            requirements[feature_key] += 1
    return requirements


def _nontrivial_features(dna):
    feats = [feat for feat in dna.dnafeatures if feat.feature_type not in ("source", "primer", "primer_bind")]
    feats.sort(key=lambda feat: (int(feat.start), int(feat.end), str((feat.qualifiers.get("label") or [""])[0])))
    return feats


def _neighbor_context(source_dna, feat):
    feats = _nontrivial_features(source_dna)
    for idx, candidate in enumerate(feats):
        if candidate is feat:
            prev_feat = feats[idx - 1] if idx > 0 else None
            next_feat = feats[idx + 1] if idx + 1 < len(feats) else None
            return prev_feat, next_feat
    return None, None


def _matching_product_features(product, feat, source_dna):
    feature_key = _feature_equivalence_key(source_dna, feat)
    matches = []
    for prod_feat in product.dnafeatures:
        if _feature_equivalence_key(product, prod_feat) == feature_key:
            matches.append(prod_feat)
    return matches


def _has_contextually_equivalent_feature(product, feat, source_dna):
    prev_feat, next_feat = _neighbor_context(source_dna, feat)
    if prev_feat is None or next_feat is None:
        return False

    prev_matches = _matching_product_features(product, prev_feat, source_dna)
    next_matches = _matching_product_features(product, next_feat, source_dna)
    if len(prev_matches) != 1 or len(next_matches) != 1:
        return False

    prod_prev = prev_matches[0]
    prod_next = next_matches[0]
    left_gap = int(feat.start) - int(prev_feat.end)
    right_gap = int(next_feat.start) - int(feat.end)
    expected_start = int(prod_prev.end) + left_gap
    expected_end = int(prod_next.start) - right_gap
    feature_key = _feature_equivalence_key(source_dna, feat)

    for prod_feat in product.dnafeatures:
        if _feature_equivalence_key(product, prod_feat) != feature_key:
            continue
        if int(prod_feat.start) == expected_start and int(prod_feat.end) == expected_end:
            return True
    return False


def _rescue_missing_feature_by_context(product, feat, source_dna):
    prev_feat, next_feat = _neighbor_context(source_dna, feat)
    if prev_feat is None or next_feat is None:
        return False

    prev_matches = _matching_product_features(product, prev_feat, source_dna)
    next_matches = _matching_product_features(product, next_feat, source_dna)
    if len(prev_matches) != 1 or len(next_matches) != 1:
        return False

    prod_prev = prev_matches[0]
    prod_next = next_matches[0]
    feat_seq = _oriented_feature_seq(source_dna, feat)
    if feat_seq == "":
        return False

    left_gap = int(feat.start) - int(prev_feat.end)
    right_gap = int(next_feat.start) - int(feat.end)
    new_start = int(prod_prev.end) + left_gap
    new_end = int(prod_next.start) - right_gap
    if new_end <= new_start:
        return False

    strand = feat.strand if feat.strand not in (None, 0) else 1
    product_seq = str(product.printsequence(new_start, new_end, strand=strand, display=False))
    if product_seq != feat_seq:
        return False

    feature_dict = {
        "feature_type": feat.feature_type,
        "start": int(new_start),
        "end": int(new_end),
        "strand": int(strand),
    }
    for key, value in feat.qualifiers.items():
        if key == "broken_feature":
            continue
        if isinstance(value, list):
            if len(value) == 0:
                continue
            feature_dict["qualifier:{}".format(key)] = value[0]
        else:
            feature_dict["qualifier:{}".format(key)] = value
    product.setfeature(feature_dict)
    return True


def _has_equivalent_feature(product, feat, source_dna):
    label = str((feat.qualifiers.get("label") or [""])[0])
    feat_type = feat.feature_type
    strand = feat.strand if feat.strand not in (None, 0) else 1
    source_seq = _oriented_feature_seq(source_dna, feat)
    for prod_feat in product.dnafeatures:
        prod_label = str((prod_feat.qualifiers.get("label") or [""])[0])
        prod_type = prod_feat.feature_type
        prod_strand = prod_feat.strand if prod_feat.strand not in (None, 0) else 1
        if prod_label != label or prod_type != feat_type or prod_strand != strand:
            continue
        if _oriented_feature_seq(product, prod_feat) == source_seq:
            return True
    return False


def _rescue_missing_features_by_exact_sequence(product, source_dnas):
    rescued = 0
    requirements = _source_feature_requirements(source_dnas)
    stalled_exact = set()
    changed = True
    while changed:
        changed = False
        product_counts = collections.Counter()
        for prod_feat in product.dnafeatures:
            if prod_feat.feature_type in ("source", "primer", "primer_bind"):
                continue
            product_counts[_feature_equivalence_key(product, prod_feat)] += 1

        for source_dna in source_dnas:
            for feat in source_dna.dnafeatures:
                if feat.feature_type in ("source", "primer", "primer_bind"):
                    continue
                feature_key = _feature_equivalence_key(source_dna, feat)
                if requirements[feature_key] == 1:
                    if _has_equivalent_feature(product, feat, source_dna):
                        continue
                else:
                    if _has_contextually_equivalent_feature(product, feat, source_dna):
                        continue

                feature_seq = _oriented_feature_seq(source_dna, feat)
                if feature_seq == "":
                    continue

                strand = feat.strand if feat.strand not in (None, 0) else 1
                hits = product.searchsequence(query=feature_seq, quinable=False)
                hits = [hit for hit in hits if (hit.strand if hit.strand not in (None, 0) else 1) == strand]
                if len(hits) == 1:
                    hit = hits[0]
                    exact_hit_key = (feature_key, int(hit.start), int(hit.end), int(strand))
                    if exact_hit_key in stalled_exact:
                        continue
                    feature_dict = {
                        "feature_type": feat.feature_type,
                        "start": int(hit.start),
                        "end": int(hit.end),
                        "strand": int(strand),
                    }
                    for key, value in feat.qualifiers.items():
                        if key == "broken_feature":
                            continue
                        if isinstance(value, list):
                            if len(value) == 0:
                                continue
                            feature_dict["qualifier:{}".format(key)] = value[0]
                        else:
                            feature_dict["qualifier:{}".format(key)] = value
                    product.setfeature(feature_dict)
                    _drop_redundant_subfeatures(product, product.dnafeatures[-1])
                    if requirements[feature_key] == 1:
                        satisfied = _has_equivalent_feature(product, feat, source_dna)
                    else:
                        satisfied = _has_contextually_equivalent_feature(product, feat, source_dna)
                    if satisfied:
                        rescued += 1
                        changed = True
                    else:
                        # Overlapping source fragments can legitimately share the
                        # same rescued sub-feature. If re-adding the exact same
                        # feature at the exact same product coordinates does not
                        # increase satisfaction, do not retry forever.
                        stalled_exact.add(exact_hit_key)
                    continue

                if _rescue_missing_feature_by_context(product, feat, source_dna):
                    rescued += 1
                    changed = True
    return product, rescued


def _primer_site_failure_message(template, target, amplicon_region, fw_candidates_total, rv_candidates_total,
                                 fw_candidates_kept, rv_candidates_kept, primer_length, fw_margin, rv_margin,
                                 nonspecific_limit):
    issues = []
    if fw_candidates_kept == 0:
        issues.append(
            "forward side has no unique primer candidates after specificity filtering "
            f"({fw_candidates_total} candidates tested, nonspecific_limit={nonspecific_limit})"
        )
    if rv_candidates_kept == 0:
        issues.append(
            "reverse side has no unique primer candidates after specificity filtering "
            f"({rv_candidates_total} candidates tested, nonspecific_limit={nonspecific_limit})"
        )

    advice = []
    if fw_candidates_kept == 0:
        advice.append(f"raise fw_margin from {fw_margin} to move the forward primer site away from repetitive terminal sequence")
    if rv_candidates_kept == 0:
        advice.append(f"raise rv_margin from {rv_margin} to move the reverse primer site away from repetitive terminal sequence")
    advice.append(f"increase primer_length beyond {primer_length}")

    try:
        target_seq = str(target.seq)
        template_seq = str(template.seq)
        target_rcseq = str(target.rcseq)
        if target_seq not in template_seq and target_rcseq in template_seq:
            advice.append("if this block is intended in the opposite orientation, try the reverse-complement donor/target view (for example [::-1])")
    except Exception:
        pass

    detail = "; ".join(issues) if len(issues) > 0 else "no unique primer candidates remained after specificity filtering"
    suggestions = "; ".join(advice)
    return (
        "No proper primer binding sites were found. "
        + detail
        + ". "
        + "Current parameters: "
        + f"primer_length={primer_length}, fw_margin={fw_margin}, rv_margin={rv_margin}. "
        + "Suggested next steps: "
        + suggestions
        + "."
    )

def sanger(template, primer, length=1000):
    """
    Return the template for sanger sequencing. Default length is 1000.
    This is not a quinable function.
    """
    if type(fw) == Qseq or type(fw) == str:
        fw    = QUEEN(seq=fw, ssdna=True) 
        fwstr = fw.seq
    elif type(fw) == QUEEN:
        fwstr = 'QUEEN.dna_dict["{}"]'.format(fw._product_id)
    site = template.searchsequence(query=primer.seq, quinable=False)
    
    if len(site) > 1:
        raise ValueError("Multiple primer binding sites were detected. You should re-design the primer sequneces.") 
    return template[site[0].start:site[0].start + length]


def pcr(template, fw, rv, bindnum=15, mismatch=0, endlength=3, add_primerbind=False, tm_func=None, return_tm=False, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """Simulate an in silico PCR reaction on one or more templates.

    This function performs a PCR (Polymerase Chain Reaction) on one or more
    `QUEEN` templates using a forward and reverse primer. Primer binding
    sites are searched on the template(s) according to the specified
    binding length, mismatch tolerance, and 3′-end constraints.

    If ``template`` is a list of `QUEEN` objects, the elements are treated
    as an ordered series of fragments for overlap‑extension (fusion) PCR:
    adjacent templates must share sufficient homology for the primers and
    overlaps to join them into a single amplicon.

    Parameters
    ----------
    template : QUEEN or list of QUEEN
        Template DNA molecule(s) to be amplified. When a list is given,
        it is interpreted as an ordered set of fragments for overlap‑
        extension PCR; the 3′ end of each fragment must be compatible
        with the 5′ end of the next fragment.
    fw : QUEEN or str
        Forward primer sequence. May be provided as:

        * a `QUEEN` object (ssDNA or dsDNA), in which case the top strand
          is used as the primer sequence (5′→3′), or
        * a plain DNA string representing the primer sequence (5′→3′).

    rv : QUEEN or str
        Reverse primer sequence. Same conventions as for ``fw``. If a
        dsDNA `QUEEN` is given, its top strand sequence is used and
        reverse‑complemented when searching for the binding site.
    bindnum : int, optional
        Minimum number of contiguous matching bases required for a
        candidate primer binding site. Default is ``15``.
    mismatch : int, optional
        Maximum number of mismatches allowed over the entire binding
        region. Default is ``0`` (perfect match only).
    endlength : int, optional
        Minimum number of perfectly matched bases required at the primer
        3′ end (the last bases toward the 3′ direction). Default is ``3``.
    add_primerbind : bool, optional
        If ``True``, primer binding sites are annotated as `DNAfeature`
        objects and added to both the template and product `QUEEN`
        objects. Default is ``False``.
    tm_func : callable, optional
        Function used to calculate primer melting temperatures. By
        default, this is equivalent to
        :func:`Bio.SeqUtils.MeltingTemp.Tm_NN` (SantaLucia nearest‑neighbor
        model) or :func:`QUEEN.qexperiment.Tm_NN`. The function must accept
        a ``seq`` keyword argument and may accept additional keyword
        arguments supplied via ``**kwargs``.
    return_tm : bool, optional
        If ``True``, also return the forward and reverse primer melting
        temperatures computed for the bound regions (adapter sequences
        are excluded). Default is ``False``.
    product : str, optional
        Human‑readable name for the PCR product. When provided, it is
        stored in the resulting `QUEEN` object (for example in the
        ``.product`` or ``.project`` attributes) and recorded in the
        construction history.
    process_name : str, optional
        Short label for this PCR step in the construction history. If
        ``None`` and ``pn`` is also ``None``, a default such as
        ``"PCR"`` is used.
    process_description : str, optional
        Free‑text description of the PCR step. Stored in the construction
        history for later export via :meth:`QUEEN.printprotocol` or
        :meth:`QUEEN.outputgbk`.
    pn : str, optional
        Alias for ``process_name``. Used only when ``process_name`` is
        ``None``.
    pd : str, optional
        Alias for ``process_description``. Used only when
        ``process_description`` is ``None``.
    **kwargs
        Additional keyword arguments forwarded directly to ``tm_func``.

    Returns
    -------
    amplicon : QUEEN
        The assembled PCR product spanning between the selected forward
        and reverse primer binding sites on the (possibly fused) template.
    tm : tuple of float, optional
        Only returned when ``return_tm`` is ``True``. A pair
        ``(fw_tm, rv_tm)`` giving the forward and reverse primer melting
        temperatures in degrees Celsius.

    Raises
    ------
    TypeError
        If ``template`` is neither a `QUEEN` object nor a list of `QUEEN`
        objects, or if ``fw``/``rv`` are neither `QUEEN` objects nor
        strings, or if ``tm_func`` is not callable.
    ValueError
        If no valid binding site is found for either primer; if multiple
        candidate binding sites remain and cannot be resolved
        unambiguously; if both primers bind to the same strand; or if the
        overlap‑extension templates are incompatible with the specified
        homology and primer constraints.

    Examples
    --------
    Amplify a region from a single circular template::

        amplicon = pcr(
            template=plasmid,
            fw="ACGTTGACT...",
            rv="TCAGCTTGA...",
            product="example_pcr"
        )

    Obtain primer melting temperatures in addition to the amplicon::

        amplicon, (fw_tm, rv_tm) = pcr(
            template=plasmid,
            fw=fw_primer,
            rv=rv_primer,
            return_tm=True
        )
    """
    def search_binding_site(template, primer, strand=1, bindnum=15, endlength=3, mismatch=1, flag=1, pn=None, pd=None, **kwargs): 
        site = [] 
        primer_end = primer.seq[-1*endlength:]
        for i in range(bindnum-endlength, len(primer.seq)-endlength+1):
            binding_site = primer.seq[-1*i + -1*endlength:-1*endlength]
            try:
                site = template.searchsequence(query="{}(?:{}){{s<={}}}{}".format(binding_site[0], binding_site[1:], mismatch, primer_end), quinable=False) 
            except Exception as e:
                site = [] 
            if len(primer_end) + len(binding_site) - mismatch >= bindnum:
                break
            else:
                pass 
        
        if len(site) == 1: #and site[0].strand == strand:
            site = template.searchsequence(query="{}(?:{}){{s<={}}}{}".format(binding_site[0], binding_site[1:], mismatch, primer_end), qexd=True, pn=pn, pd=pd)
            if flag == 1:
                return site[0]
            else:
                return site
        
        elif len(site) == 0:
            if flag == 1:
                raise ValueError("No primer binding sites were found. You should re-confirm the template-primer pair.")
            else:
                return site 

        elif len(site) > 1:
            if mismatch > 0:
                premismatch = mismatch
                site = search_binding_site(template, primer, strand, bindnum, endlength, 0, flag, pn, pd, **kwargs) 
                if premismatch > 1:
                    print("**Attention**: Multiple potential primer binding sites with 1-to-{} mismatches were found. It is recommended to redesign the primer sequences.".format(premismatch))
                else:
                    print("**Attention**: Multiple potential primer binding sites with a single mismatche were found. It is recommended to redesign the primer sequences.")
                return site 
            else:
                raise ValueError("Multiple primer binding sites were detected. You should re-design the primer sequneces.") 
         
    if type(template) != QUEEN: 
        if type(template) == list and False not in [type(element) == QUEEN for element in template]:
            pass
        else:
            raise TypeError("`template` object must be instance of QUEEN class or list of QUEEN objects") 

    if type(fw) == Qseq or type(fw) == str:
        fw    = QUEEN(seq=fw, ssdna=True) 
        fwstr = fw.seq
    elif type(fw) == QUEEN:
        fwstr = 'QUEEN.dna_dict["{}"]'.format(fw._product_id)
    else:
        raise TypeError("`fw` object must be instance of QUEEN or str class.") 

    if type(rv) == Qseq or type(rv) == str:
        rv    = QUEEN(seq=rv, ssdna=True) 
        rvstr = rv.seq
    elif type(rv) == QUEEN:
        rvstr = 'QUEEN.dna_dict["{}"]'.format(rv._product_id)
    else:
        raise TypeError("`rv` object must be instance of QUEEN or str class.") 
     
    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "PCR"
    
    #product = product.replace(" ","") if product is not None else None
    kwargs_str = _convert_kwargs(kwargs)
    
    if bindnum != 15:
        bindnumtxt = ", bindnum={}".format(bindnum) 
    else:
        bindnumtxt = ""
    
    if mismatch != 0:
        mismatchtxt = ", mismatch={}".format(mismatch) 
    else:
        mismatchtxt = ""
    
    if endlength != 3:
        endlengthtxt = ", endlength={}".format(endlegnth) 
    else:
        endlengthtxt = ""
    
    if add_primerbind == True:
        aptxt = ", add_primerbind={}".format(add_primerbind) 
    else:
        aptxt = ""
    
    if type(template) == list and False not in [type(element) == QUEEN for element in template]:
        temps = [] 
        for atemp in template:
            temps.append('QUEEN.dna_dict["{}"]'.format(atemp._product_id))
        temps = "[{}]".format(", ".join(temps))
        qexd  = 'pcr({}, {}, {}{}{}{}{}{})'.format(temps, fwstr, rvstr, bindnumtxt, mismatchtxt, endlengthtxt, aptxt, kwargs_str)
        template = homology_based_assembly(*template, mode="overlappcr") 
    else: 
        qexd = 'pcr(QUEEN.dna_dict["{}"], {}, {}{}{}{}{}{})'.format(template._product_id, fwstr, rvstr, bindnumtxt, mismatchtxt, endlengthtxt, aptxt, kwargs_str)
    process_description = pd if process_description is None else process_description
    
    if template._ssdna == True: 
        template = joindna(template, QUEEN(seq=template.rcseq, ssdna=True), qexd=True, pn=process_name, pd=process_description)

    if -1 in [template._left_end_top, template._left_end_bottom, template._right_end_top, template._right_end_bottom] and template._ssdna == False: 
        template = modifyends(template, qexd=True, pn=process_name, pd=process_description)
    
    i = 0
    site1      = None 
    tmpsite1   = None
    bindlength = len(fw.seq) 
    while bindlength-i >= bindnum:
        tmpsite1 = search_binding_site(template, fw, 1, bindlength-i, endlength, mismatch=mismatch, flag=0, qexd=True, pn=process_name, pd=process_description)
        if len(tmpsite1) == 1:
            site1 = tmpsite1[0]
            break
        i += 1

    i = 0
    site2      = None
    tmpsite2   = None
    bindlength = len(rv.seq) 
    while bindlength-i >= bindnum:
        tmpsite2 = search_binding_site(template, rv, -1, bindlength-i, endlength, mismatch=mismatch, flag=0, qexd=True, pn=process_name, pd=process_description)
        if len(tmpsite2) == 1:
            site2 = tmpsite2[0] 
            break
        i += 1
    
    if site1 is None:
        raise ValueError("No forward primer binding sites were found. You should re-confirm the template-primer pair.")
    
    if site2 is None:
        raise ValueError("No reverse primer binding sites were found. You should re-confirm the template-primer pair.")
    
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

    fw_feats = [feat for feat in fw.searchfeature(key_attribute="feature_type", query="primer_bind", qexd=True, pn=process_name, pd=process_description) if feat.end == len(fw.seq) and feat.start == 0] 
    rv_feats = [feat for feat in rv.searchfeature(key_attribute="feature_type", query="primer_bind", qexd=True, pn=process_name, pd=process_description) if feat.end == len(rv.seq) and feat.start == 0]

    rescue_sources = []
    if fw_site.end >= rv_site.start and fw_site.start < rv_site.start:
        if len(fw_feats) == 0:
            fw.setfeature({"qualifier:label":"{}".format(fw.project), "feature_type":"primer_bind"})  
        if len(rv_feats) == 0:
            rv.setfeature({"qualifier:label":"{}".format(rv.project), "feature_type":"primer_bind"})
        
        start = rv_site.start if rv_site.start < len(template.seq) else rv_site.start - len(template.seq)
        end   = fw_site.end if fw_site.end < len(template.seq) else fw_site.end - len(template.seq)
        extract  = cropdna(template, start, end, qexd=True, pn=process_name, pd=process_description)
        rescue_sources = [extract]
        fw_index = len(fw.seq) - (fw_site.end - rv_site.start) 
        rv_index = fw_site.end - rv_site.start 
        amplicon = modifyends(extract, fw.seq[:fw_index], rv.rcseq[rv_index:], qexd=qexd, product=product, pn=process_name, pd=process_description)
    else:
        req1 = (mismatch == 0 and ((fw_site.start < rv_site.end and fw_site.start >= rv_site.start) == False))
        req2 = (mismatch > 0 and ((fw_site.start < rv_site.end and fw_site.start >= rv_site.start) == False) and fw_bind_length == len(fw_site.sequence) and rv_bind_length == len(rv_site.sequence))
        
        if len(fw_feats) == 0:
            fw_label = fw.project
        else:
            fw_label = fw_feats[0].qualifiers["label"][0]
            editfeature(fw, query=fw_feats[0].feature_id, key_attribute="feature_id", target_attribute="feature_id", operation=removeattribute(), new_copy=False, quinable=False) 
        
        if len(rv_feats) == 0:
            rv_label = rv.project
        else:
            rv_label = rv_feats[0].qualifiers["label"][0] 
            editfeature(rv, query=rv_feats[0].feature_id, key_attribute="feature_id", target_attribute="feature_id", operation=removeattribute(), new_copy=False, quinable=False) 

        if req1 == True:
            fw_bind  = template.seq[fw_site.start:fw_site.end]
            rv_bind  = template.seq[rv_site.start:rv_site.end]
            start    = fw_site.start if fw_site.start < len(template.seq) else fw_site.start - len(template.seq)
            end      = rv_site.end if rv_site.end < len(template.seq) else rv_site.end - len(template.seq)
            extract  = cropdna(template, start, end, qexd=True, pn=process_name, pd=process_description)
            rescue_sources = [extract]
            
            amplicon = modifyends(extract, left=fw[0:len(fw.seq)-len(fw_bind)].seq, right=rv[0:len(rv.seq)-len(rv_bind)].rcseq, qexd=qexd, product=product, pn=process_name, pd=process_description)   
            amplicon.setfeature({"start":0, "end":len(fw.seq), "qualifier:label":"{}".format(fw_label), "feature_type":"primer_bind"})
            amplicon.setfeature({"start":len(amplicon.seq)-len(rv.seq), "end":len(amplicon.seq), "strand":-1, "qualifier:label":"{}".format(rv_label), "feature_type":"primer_bind"})  
        
        elif req2 == True:
            fw_mut, rv_mut = 0, 0
            start = fw_site.start if fw_site.start < len(template.seq) else fw_site.start - len(template.seq)
            end   = rv_site.end if rv_site.end <= len(template.seq) else rv_site.end - len(template.seq)

            if fw.seq != template.seq[start:start+len(fw.seq)]:
                extract_fw = cropdna(template, start, start+len(fw.seq), qexd=True, pn=process_name, pd=process_description)
                extract_fw._seq = fw.seq
                extract_fw._seq.qkey        = template.seq.qkey
                extract_fw._seq.parent      = template.seq.parent
                extract_fw._seq.parental_id = template.seq.parental_id
                extract_fw._seq.name        = template.seq.name
                extract_fw._seq.item        = template.seq.item
                fw_mut = 1
            
            if  rv.rcseq != template.seq[end-len(rv.seq):end]:
                extract_rv = cropdna(template, end-len(rv.seq), end, qexd=True, pn=process_name, pd=process_description)
                extract_rv._seq = rv.rcseq
                extract_rv._seq.qkey        = template.seq.qkey
                extract_rv._seq.parent      = template.seq.parent
                extract_rv._seq.parental_id = template.seq.parental_id
                extract_rv._seq.name        = template.seq.name
                extract_rv._seq.item        = template.seq.item
                rv_mut = 1
            
            if fw_mut == True and rv_mut == True:
                extract  = cropdna(template, start+len(fw.seq), end-len(rv.seq), qexd=True, pn=process_name, pd=process_description)
                amplicon = joindna(extract_fw, extract, extract_rv, qexd=qexd, product=product, pn=process_name, pd=process_description)
                rescue_sources = [extract_fw, extract, extract_rv]
            elif fw_mut == True:
                extract  = cropdna(template, start+len(fw.seq), end, qexd=True, pn=process_name, pd=process_description)
                amplicon = joindna(extract_fw, extract, qexd=qexd, product=product, pn=process_name, pd=process_description)
                rescue_sources = [extract_fw, extract]
            elif rv_mut == True:
                extract  = cropdna(template, start, end-len(rv.seq), qexd=True, pn=process_name, pd=process_description)
                amplicon = joindna(extract, extract_rv, qexd=qexd, product=product, pn=process_name, pd=process_description)
                rescue_sources = [extract, extract_rv]
            else:
                amplicon = cropdna(template, start, end, qexd=qexd, pn=process_name, pd=process_description)
                rescue_sources = [amplicon]

            amplicon.setfeature({"start":0, "end":len(fw.seq), "qualifier:label":"{}".format(fw_label), "feature_type":"primer_bind"})
            amplicon.setfeature({"start":len(amplicon.seq)-len(rv.seq), "end":len(amplicon.seq), "strand":-1, "qualifier:label":"{}".format(rv_label), "feature_type":"primer_bind"})  
        
        else:
            if len(fw_feats) == 0:
                fw.setfeature({"qualifier:label":"{}".format(fw.project), "feature_type":"primer_bind"})  
            if len(rv_feats) == 0:
                rv.setfeature({"qualifier:label":"{}".format(rv.project), "feature_type":"primer_bind"})
            
            start = fw_site.end if fw_site.end < len(template.seq) else fw_site.end - len(template.seq)
            end   = rv_site.start if rv_site.start < len(template.seq) else rv_site.start - len(template.seq)
            extract  = cropdna(template, start, end, qexd=True, pn=process_name, pd=process_description)
            rescue_sources = [extract]
            amplicon = modifyends(extract, fw.seq, rv.rcseq, qexd=qexd, product=product, pn=process_name, pd=process_description)

    histories = [amplicon._history, fw._history, rv._history]
    combined_history  = _combine_history(amplicon, histories)
    amplicon._history = combined_history
    amplicon, rescued_count = _rescue_missing_features_by_exact_sequence(amplicon, rescue_sources if len(rescue_sources) > 0 else [template])
    amplicon._rescued_feature_count = rescued_count

    if add_primerbind == True:
        template.setfeature({"start": fw_site.start, "end": fw_site.end, "strand":1,  "feature_type":"primer_bind"}) 
        template.setfeature({"start": rv_site.start, "end": rv_site.end, "strand":-1, "feature_type":"primer_bind"})  
    
    if return_tm == True:
        if tm_func is None:
            tm_func = Tm_NN() 
        
        elif tm_func == "SantaLucia" or "sa":
            tm_func = Tm_NN(nn_table=mt.DNA_NN3)

        elif tm_func == "Breslauer" or "br":
            tm_func = Tm_NN(nn_table=mt.DNA_NN1)

        fw_tm = tm_func(seq=template.seq[fw_site.start:fw_site.end]) 
        rv_tm = tm_func(seq=template.seq[rv_site.start:rv_site.end])
        return amplicon, (fw_tm, rv_tm) 
    else:
        return amplicon

def _select(fragments, selection=None, process_name=None, process_description=None): 
    if selection is None:
        return fragments  
    
    elif type(selection) == int:
        fragments.sort(key=lambda x:abs(len(x.seq)-selection))
        return fragments[0] 

    elif type(selection) == tuple:
        fragments = [fragment for fragment in fragments if min(size_range) <= len(fragment.seq) <= max(size_range)] 
        if len(fragments) > 1:
            raise ValueError("Multiple fragments were detected within the specified range.") 
        elif len(fragments) == 0:
            raise ValueError("No fragment was detected within the specified range.") 
        return fragments[0] 

    elif selection in ("min", "max"):
        fragments.sort(key=lambda x:len(x.seq))
        if selection == "min":
            return fragments[0]
        else:
            return fragments[-1]
    
    elif selection.startswith("!") == False and ":" in selection: 
        query = ":".join(selection.split(":")[1:])
        fragments = [fragment for fragment in fragments if len(fragment.searchfeature(key_attribute="qualifier:{}".format(selection.split(":")[0]), query=query, qexd=True, pn=process_name, pd=process_description)) > 0]
        if len(fragments) > 1:
            raise ValueError("Multiple fragments holding the specified feature were detected") 
        elif len(fragments) == 0:
            raise ValueError("No fragment holding the specified feature was detected.") 
        return fragments[0]
    
    elif selection.startswith("!") and ":" in selection: 
        query = ":".join(selection.split(":")[1:])
        fragments = [fragment for fragment in fragments if len(fragment.searchfeature(key_attribute="qualifier:{}".format(selection.split(":")[0][1:]), query=query, qexd=True, pn=process_name, pd=process_description)) == 0]
        if len(fragments) > 1:
            raise ValueError("Multiple fragments holding the specified feature were detected") 
        elif len(fragments) == 0:
            raise ValueError("No fragment holding the specified feature was detected.")
        return fragments[0] 


def _site_cut_interval(dna, site_feature):
    if "cutsite" not in site_feature.qualifiers:
        raise ValueError("DNAfeature object should hold 'qualifiers:cutsite' attribute.")

    if site_feature._digestion_topl == "null":
        _, _, site_feature._digestion_topl, site_feature._digestion_topr, site_feature._digestion_bottoml, site_feature._digestion_bottomr = compile_cutsite(site_feature.qualifiers["cutsite"][0])

    strand = site_feature.location.strand
    if strand != -1:
        if site_feature._digestion_topl != "null":
            pos1 = int(site_feature.start) - int(site_feature._digestion_topl)
            pos2 = int(site_feature.start) - int(site_feature._digestion_bottoml)
        else:
            pos1 = int(site_feature.end) + int(site_feature._digestion_topr)
            pos2 = int(site_feature.end) + int(site_feature._digestion_bottomr)
    else:
        if site_feature._digestion_topr != "null":
            pos1 = int(site_feature.start) - int(site_feature._digestion_bottomr)
            pos2 = int(site_feature.start) - int(site_feature._digestion_topr)
        else:
            pos1 = int(site_feature.end) + int(site_feature._digestion_bottoml)
            pos2 = int(site_feature.end) + int(site_feature._digestion_topl)

    length = len(dna.seq)
    pos1 %= length
    pos2 %= length
    return tuple(sorted((pos1, pos2)))


def _target_interval_in_source(source, target):
    positions = getattr(target, "_positions", None)
    if type(positions) in (tuple, list) and len(positions) == len(target.seq) and len(positions) > 0:
        source_len = len(source.seq)
        positions = [int(pos) % source_len for pos in positions]
        if len(positions) == 1:
            return positions[0], positions[0] + 1, 1

        diffs = [((positions[i + 1] - positions[i]) % source_len) for i in range(len(positions) - 1)]
        if all(diff == 1 for diff in diffs):
            start = positions[0]
            return start, start + len(positions), 1
        if all(diff == (source_len - 1) for diff in diffs):
            start = positions[-1]
            return start, start + len(positions), -1

    queries = [str(target.seq)]
    rcseq = str(target.rcseq)
    if rcseq != queries[0]:
        queries.append(rcseq)

    matches = []
    for query in queries:
        sites = source.searchsequence(query=query)
        for site in sites:
            start = int(site.start)
            end = int(site.end)
            strand = 1 if site.strand in (None, 0) else int(site.strand)
            if end <= start:
                end += len(source.seq)
            matches.append((start, end, strand))

    matches = sorted(set(matches))
    if len(matches) == 0:
        raise ValueError("`target` sequence was not found in `source`.")
    if len(matches) > 1:
        raise ValueError("`target` sequence mapped to multiple locations in `source`; provide a unique target fragment.")
    return matches[0]


def _canonical_pair_key(left_enzyme, right_enzyme):
    return "|".join(sorted((str(left_enzyme), str(right_enzyme))))


def _normalize_enzyme_constraints(enzyme_set):
    if enzyme_set is None:
        return None, None

    allowed_names = set()
    allowed_pair_keys = set()

    if hasattr(enzyme_set, "columns"):
        columns = set(str(col) for col in enzyme_set.columns)
        if "pair_key" in columns:
            for pair_key in enzyme_set["pair_key"].dropna().astype(str).tolist():
                parts = pair_key.split("|")
                if len(parts) != 2:
                    raise ValueError("`enzyme_set` DataFrame contains an invalid `pair_key` value.")
                left_name, right_name = parts
                if left_name not in cs.lib.keys() or right_name not in cs.lib.keys():
                    raise ValueError("`enzyme_set` DataFrame contains an unknown restriction enzyme name.")
                allowed_pair_keys.add(_canonical_pair_key(left_name, right_name))
                allowed_names.add(left_name)
                allowed_names.add(right_name)
            return allowed_names, allowed_pair_keys
        if "pair" in columns:
            for pair in enzyme_set["pair"].dropna().astype(str).tolist():
                parts = pair.split("|")
                if len(parts) != 2:
                    raise ValueError("`enzyme_set` DataFrame contains an invalid `pair` value.")
                left_name, right_name = parts
                if left_name not in cs.lib.keys() or right_name not in cs.lib.keys():
                    raise ValueError("`enzyme_set` DataFrame contains an unknown restriction enzyme name.")
                allowed_pair_keys.add(_canonical_pair_key(left_name, right_name))
                allowed_names.add(left_name)
                allowed_names.add(right_name)
            return allowed_names, allowed_pair_keys
        if "left_enzyme" in columns and "right_enzyme" in columns:
            left_series = enzyme_set["left_enzyme"].dropna().astype(str).tolist()
            right_series = enzyme_set["right_enzyme"].dropna().astype(str).tolist()
            for left_name, right_name in zip(left_series, right_series):
                if left_name not in cs.lib.keys() or right_name not in cs.lib.keys():
                    raise ValueError("`enzyme_set` DataFrame contains an unknown restriction enzyme name.")
                allowed_pair_keys.add(_canonical_pair_key(left_name, right_name))
                allowed_names.add(left_name)
                allowed_names.add(right_name)
            return allowed_names, allowed_pair_keys
        raise ValueError("`enzyme_set` DataFrame must contain `pair_key`, `pair`, or both `left_enzyme` and `right_enzyme` columns.")

    if type(enzyme_set) == str:
        enzyme_set = [enzyme_set]

    for enzyme in enzyme_set:
        if type(enzyme) in (tuple, list) and len(enzyme) == 2:
            left_name = enzyme[0].name if (type(enzyme[0]) == Cutsite or "cutsite" in getattr(enzyme[0], "__dict__", {})) else str(enzyme[0])
            right_name = enzyme[1].name if (type(enzyme[1]) == Cutsite or "cutsite" in getattr(enzyme[1], "__dict__", {})) else str(enzyme[1])
            if left_name not in cs.lib.keys() or right_name not in cs.lib.keys():
                raise ValueError("`enzyme_set` contains an unknown restriction enzyme name.")
            allowed_pair_keys.add(_canonical_pair_key(left_name, right_name))
            allowed_names.add(left_name)
            allowed_names.add(right_name)
            continue

        if type(enzyme) == Cutsite or "cutsite" in getattr(enzyme, "__dict__", {}):
            allowed_names.add(enzyme.name)
            continue

        enzyme_name = str(enzyme)
        if "|" in enzyme_name:
            parts = enzyme_name.split("|")
            if len(parts) != 2:
                raise ValueError("`enzyme_set` contains an invalid pair string.")
            left_name, right_name = parts
            if left_name not in cs.lib.keys() or right_name not in cs.lib.keys():
                raise ValueError("`enzyme_set` contains an unknown restriction enzyme name.")
            allowed_pair_keys.add(_canonical_pair_key(left_name, right_name))
            allowed_names.add(left_name)
            allowed_names.add(right_name)
            continue
        if enzyme_name not in cs.lib.keys():
                raise ValueError("`enzyme_set` contains an unknown restriction enzyme name.")
        allowed_names.add(enzyme_name)
    return allowed_names, (allowed_pair_keys if len(allowed_pair_keys) > 0 else None)


def _normalize_max_distance(max_distance):
    if max_distance is None:
        return None, None
    if type(max_distance) == int:
        if max_distance < 0:
            raise ValueError("`max_distance` must be >= 0.")
        return max_distance, max_distance
    if type(max_distance) in (tuple, list) and len(max_distance) == 2:
        left_max, right_max = max_distance
        if type(left_max) != int or type(right_max) != int:
            raise TypeError("`max_distance` tuple values must be integers.")
        if left_max < 0 or right_max < 0:
            raise ValueError("`max_distance` tuple values must be >= 0.")
        return left_max, right_max
    raise TypeError("`max_distance` must be None, an integer, or a tuple/list of two integers.")


def _apply_preferred_max_distance(df, preferred_max_distance):
    left_pref_distance, right_pref_distance = _normalize_max_distance(preferred_max_distance)
    mask = (
        (df["left_site_boundary_distance_bp"] <= left_pref_distance)
        & (df["right_site_boundary_distance_bp"] <= right_pref_distance)
    )
    if bool(mask.any()) is True:
        return df.loc[mask].reset_index(drop=True)
    return df


def _feature_primary_label(feat):
    for key in ("label", "ApEinfo_label", "standard_name", "gene"):
        if key in feat.qualifiers and len(feat.qualifiers[key]) > 0:
            return str(feat.qualifiers[key][0])
    return ""


def _feature_location_intervals(feat, source_len):
    parts = getattr(feat.location, "parts", None)
    if parts is None or len(parts) == 0:
        parts = [feat.location]

    intervals = []
    for part in parts:
        start = int(part.start)
        end = int(part.end)
        if end < start:
            intervals.append((start, source_len))
            intervals.append((0, end))
        elif end > start:
            intervals.append((start, end))
    return intervals


def _interval_overlaps(start1, end1, start2, end2):
    return start1 < end2 and start2 < end1


def _feature_guard_conflicts(
    source,
    left_cut,
    right_cut,
    target_start,
    target_end,
    *,
    ignored_feature_types=("primer_bind", "source", "primer"),
    ignored_feature_labels=("MCS",),
):
    source_len = len(source.seq)
    left_extra = (left_cut, target_start)
    right_extra = (target_end, right_cut)
    ignored_types = {str(x) for x in ignored_feature_types}
    ignored_labels = {str(x) for x in ignored_feature_labels}
    conflicts = set()

    for feat in source.dnafeatures:
        if "broken_feature" in feat.qualifiers:
            continue
        if str(feat.type) in ignored_types:
            continue

        label = _feature_primary_label(feat)
        if label in ignored_labels:
            continue

        for base_start, base_end in _feature_location_intervals(feat, source_len):
            for shift in (-source_len, 0, source_len):
                feat_start = int(base_start + shift)
                feat_end = int(base_end + shift)
                if feat_end <= left_cut or feat_start >= right_cut:
                    continue

                if feat_start >= target_start and feat_end <= target_end:
                    continue

                reason = None
                if _interval_overlaps(feat_start, feat_end, target_start, target_end):
                    reason = "target_partial_feature_overlap"
                elif feat_start < left_cut < feat_end or feat_start < right_cut < feat_end:
                    reason = "cut_inside_feature"
                elif _interval_overlaps(feat_start, feat_end, *left_extra) or _interval_overlaps(feat_start, feat_end, *right_extra):
                    reason = "extra_feature_overlap"

                if reason is not None:
                    conflicts.add((str(feat.type), label, int(base_start), int(base_end), reason))
                    break

    return conflicts


def _infer_cutsite_candidates(
    source,
    target,
    cuttype="single",
    enzyme_set=None,
    max_distance=None,
    preferred_max_distance=500,
    protect_features=True,
    ignored_feature_types=("primer_bind", "source", "primer"),
    ignored_feature_labels=("MCS",),
):
    """Return a ranked DataFrame of flanking cutsite-pair candidates around a target core."""

    if type(source) != QUEEN:
        raise TypeError("`source` must be a QUEEN object.")
    if type(target) != QUEEN:
        raise TypeError("`target` must be a QUEEN object.")
    if cuttype not in ("single", "dual", "all", "typeIIS"):
        raise ValueError("`cuttype` should be one of 'single', 'dual', 'all', or 'typeIIS'.")

    import pandas as pd

    allowed_enzyme_names, allowed_pair_keys = _normalize_enzyme_constraints(enzyme_set)
    left_max_distance, right_max_distance = _normalize_max_distance(max_distance)

    source_len = len(source.seq)
    target_start, target_end, target_strand = _target_interval_in_source(source, target)

    site_rows = []
    if allowed_enzyme_names is None:
        enzyme_items = list(cs.lib.items())
    else:
        enzyme_items = [(enzyme_name, cs.lib[enzyme_name]) for enzyme_name in sorted(allowed_enzyme_names)]

    for enzyme_name, enzyme in enzyme_items:
        if cuttype == "typeIIS" and enzyme.IIS is not True:
            continue

        sites = source.searchsequence(query=enzyme.cutsite)
        if cuttype == "single" and len(sites) != 1:
            continue
        if cuttype == "dual" and len(sites) != 2:
            continue

        for site in sites:
            cut_lo, cut_hi = _site_cut_interval(source, site)
            for shift in (-source_len, 0, source_len):
                shifted_lo = cut_lo + shift
                shifted_hi = cut_hi + shift
                shifted_site_start = int(site.start) + shift
                shifted_site_end = int(site.end) + shift
                if shifted_hi <= target_start:
                    side = "left"
                    cut_boundary_distance = target_start - shifted_hi
                    site_boundary_distance = target_start - shifted_site_end
                elif shifted_lo >= target_end:
                    side = "right"
                    cut_boundary_distance = shifted_lo - target_end
                    site_boundary_distance = shifted_site_start - target_end
                else:
                    continue

                site_rows.append({
                    "enzyme": enzyme_name,
                    "site_start": int(site.start),
                    "site_end": int(site.end),
                    "site_strand": 1 if site.strand in (None, 0) else int(site.strand),
                    "cut_lo": int(cut_lo),
                    "cut_hi": int(cut_hi),
                    "shifted_cut_lo": int(shifted_lo),
                    "shifted_cut_hi": int(shifted_hi),
                    "shifted_site_start": int(shifted_site_start),
                    "shifted_site_end": int(shifted_site_end),
                    "side": side,
                    "cut_boundary_distance_bp": int(cut_boundary_distance),
                    "site_boundary_distance_bp": int(site_boundary_distance),
                })

    left_sites = [row for row in site_rows if row["side"] == "left"]
    right_sites = [row for row in site_rows if row["side"] == "right"]
    if len(left_sites) == 0 or len(right_sites) == 0:
        raise ValueError("No flanking cutsite pair was found around the target region.")

    rows = []
    for left in left_sites:
        for right in right_sites:
            if left["shifted_cut_hi"] > right["shifted_cut_lo"]:
                continue
            extra_left = int(target_start - left["shifted_cut_hi"])
            extra_right = int(right["shifted_cut_lo"] - target_end)
            left_site_distance = int(left["site_boundary_distance_bp"])
            right_site_distance = int(right["site_boundary_distance_bp"])
            if left_max_distance is not None and left_site_distance > left_max_distance:
                continue
            if right_max_distance is not None and right_site_distance > right_max_distance:
                continue
            conflicts = set()
            if protect_features is True:
                conflicts = _feature_guard_conflicts(
                    source,
                    int(left["shifted_cut_hi"]),
                    int(right["shifted_cut_lo"]),
                    int(target_start),
                    int(target_end),
                    ignored_feature_types=ignored_feature_types,
                    ignored_feature_labels=ignored_feature_labels,
                )
            rows.append({
                "pair": "{}|{}".format(left["enzyme"], right["enzyme"]),
                "pair_key": _canonical_pair_key(left["enzyme"], right["enzyme"]),
                "left_enzyme": left["enzyme"],
                "right_enzyme": right["enzyme"],
                "same_enzyme": left["enzyme"] == right["enzyme"],
                "target_start": int(target_start),
                "target_end": int(target_end),
                "target_strand": int(target_strand),
                "left_site_start": int(left["site_start"]),
                "left_site_end": int(left["site_end"]),
                "right_site_start": int(right["site_start"]),
                "right_site_end": int(right["site_end"]),
                "left_cut": int(left["shifted_cut_hi"]),
                "right_cut": int(right["shifted_cut_lo"]),
                "left_boundary_distance_bp": int(extra_left),
                "right_boundary_distance_bp": int(extra_right),
                "left_site_boundary_distance_bp": int(left_site_distance),
                "right_site_boundary_distance_bp": int(right_site_distance),
                "extra_span_bp": int(extra_left + extra_right),
                "fragment_span_bp": int(right["shifted_cut_lo"] - left["shifted_cut_hi"]),
                "cuttype": cuttype,
                "protected_feature_conflict_count": int(len(conflicts)),
                "protected_feature_conflict_labels": ";".join(sorted({label for _, label, _, _, _ in conflicts if label != ""})),
                "protected_feature_conflict_types": ";".join(sorted({feature_type for feature_type, _, _, _, _ in conflicts})),
                "protected_feature_conflict_reasons": ";".join(sorted({reason for _, _, _, _, reason in conflicts})),
            })

    if len(rows) == 0:
        if max_distance is None:
            raise ValueError("No ordered flanking cutsite pair was found around the target region.")
        raise ValueError("No flanking cutsite pair was found within the requested target-boundary distance.")

    df = pd.DataFrame(rows)
    if allowed_pair_keys is not None:
        df = df[df["pair_key"].isin(allowed_pair_keys)].reset_index(drop=True)
        if len(df) == 0:
            raise ValueError("No flanking cutsite pair matched the requested enzyme-pair constraint.")
    if protect_features is True:
        df = df[df["protected_feature_conflict_count"] == 0].reset_index(drop=True)
        if len(df) == 0:
            raise ValueError("No flanking cutsite pair survived protected-feature filtering around the target region.")
    df = df.sort_values(
        by=[
            "extra_span_bp",
            "left_site_boundary_distance_bp",
            "right_site_boundary_distance_bp",
            "left_boundary_distance_bp",
            "right_boundary_distance_bp",
            "fragment_span_bp",
            "left_enzyme",
            "right_enzyme",
        ],
        ascending=[True, True, True, True, True, True, True, True],
    ).reset_index(drop=True)
    if preferred_max_distance is not None:
        df = _apply_preferred_max_distance(df, preferred_max_distance)
    return df


def infer_cutsites(
    source,
    target=None,
    cuttype="single",
    enzyme_set=None,
    max_distance=None,
    preferred_max_distance=500,
    display=False,
    return_df=False,
    protect_features=True,
    ignored_feature_types=("primer_bind", "source", "primer"),
    ignored_feature_labels=("MCS",),
    **kwargs
):
    """Infer a restriction pair flanking a target core region.

    Parameters
    ----------
    source : QUEEN
        Source DNA containing the target core region to be excised. In many
        cloning workflows this will be the donor plasmid, but the helper is
        generic and can also be applied to any other source molecule that
        carries the region of interest.
    target : QUEEN
        Core region that must be retained inside the excised fragment.
        This does not need to be the exact final released fragment; it is the
        must-include region that the inferred cutsite pair should flank from
        the outside. In typical use this is a `QUEEN` slice or feature-derived
        subobject obtained from ``source``.
    cuttype : {"single", "dual", "all", "typeIIS"}, optional
        Restrict candidate enzymes to the same categories used by
        :meth:`QUEEN.printcutsite`. Default is ``"single"``.
    enzyme_set : sequence of Cutsite or str, optional
        Restrict the search space to a user-specified enzyme set. This is
        useful when another molecule, such as a backbone, has already
        constrained the permissible enzymes.
    max_distance : int or tuple(int, int), optional
        Maximum allowed distance from the target-core boundary to the left and
        right recognition sites. Distances are measured from the target core
        to the nearest edge of each enzyme recognition site, not to the exact
        cleavage positions. If an integer is provided, it is applied to both
        sides. If a two-element tuple/list is provided, it is interpreted as
        ``(left_max_bp, right_max_bp)``. If no candidate pair satisfies this
        threshold, a ``ValueError`` is raised.
    preferred_max_distance : int or tuple(int, int), optional
        Soft distance preference used after any hard ``max_distance`` filter is
        applied. If one or more candidate pairs fall within this distance on
        both sides, only those candidates are retained in the ranking. If no
        candidate pair satisfies this preference, the full candidate ranking is
        returned unchanged. Default is ``500``.
    display : bool, optional
        If ``True``, print the ranked candidate table to standard output.
        Default is ``False`` so helper use does not consume unnecessary
        context unless explicitly requested.
    return_df : bool, optional
        If ``True``, return the full ranked candidate DataFrame instead of the
        best cutsite list.
    protect_features : bool, optional
        If ``True`` (default), reject candidate pairs whose extra span or cut
        boundary would intersect protected features outside the target core.
        By default, ``primer_bind``, ``source``, and ``primer`` feature types
        are ignored, and features labeled ``MCS`` are also ignored.
    ignored_feature_types : tuple/list of str, optional
        Feature types to ignore during protected-feature filtering. Default is
        ``("primer_bind", "source", "primer")``.
    ignored_feature_labels : tuple/list of str, optional
        Feature labels to ignore during protected-feature filtering regardless
        of feature type. Default is ``("MCS",)``.

    Returns
    -------
    list of Cutsite or pandas.DataFrame
        By default, return the best-ranked restriction enzyme list ready to
        pass directly to :func:`digestion`. For dual-cutter same-enzyme cases
        this may be a one-element list, because ``digestion(source, enzyme)``
        already cuts all occurrences of that enzyme in ``source``. If
        ``return_df=True``, return the full ranked candidate DataFrame.

    Notes
    -----
    This helper is intended for a must-keep target core on a source molecule. It can be applied to donor-side excision or to backbone site selection, as long as the target region is already represented as a single contiguous `QUEEN` object on ``source``.

    Examples
    --------
    Get the best-ranked cutsite pair as `Cutsite` objects ready for
    :func:`digestion`::

        cutsites = infer_cutsites(
            donor,
            target=payload,
            cuttype="single",
        )

    Inspect the full ranked candidate table without printing it::

        donor_df = infer_cutsites(
            donor,
            target=payload,
            cuttype="single",
            return_df=True,
        )

    Reuse the donor-side ranked pairs as a pair constraint on the backbone
    side. This answers the question "which donor-valid pairs are also valid on
    the backbone?" by passing the donor ranking DataFrame directly into
    ``enzyme_set``::

        backbone_df = infer_cutsites(
            backbone,
            target=removee_window,
            cuttype="single",
            enzyme_set=donor_df,
            return_df=True,
        )

    If you want the ranked table in stdout for interactive reasoning, use
    ``display=True``::

        infer_cutsites(
            donor,
            target=payload,
            cuttype="single",
            display=True,
        )
    """

    payload = kwargs.pop("payload", None)
    if kwargs:
        raise TypeError("Unexpected keyword arguments: {}".format(", ".join(sorted(kwargs.keys()))))
    if target is None:
        target = payload
    elif payload is not None:
        raise TypeError("Pass either `target` or deprecated `payload`, not both.")
    if target is None:
        raise TypeError("`target` must be provided.")

    df = _infer_cutsite_candidates(
        source=source,
        target=target,
        cuttype=cuttype,
        enzyme_set=enzyme_set,
        max_distance=max_distance,
        preferred_max_distance=preferred_max_distance,
        protect_features=protect_features,
        ignored_feature_types=ignored_feature_types,
        ignored_feature_labels=ignored_feature_labels,
    )

    if display is True:
        print(df.to_string(index=False))

    if return_df is True:
        return df

    top_row = df.iloc[0]
    if bool(top_row["same_enzyme"]) is True:
        return [cs.lib[str(top_row["left_enzyme"])]]
    return [cs.lib[str(top_row["left_enzyme"])], cs.lib[str(top_row["right_enzyme"])]]

def digestion(dna, *cutsites, selection=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """Simulate restriction digestion of a `QUEEN` object.

    The function digests an input `QUEEN` object using one or more
    restriction endonuclease specificities (``Cutsite`` objects or
    enzyme names). The resulting fragments can optionally be filtered
    or reduced to a single fragment based on size or feature content.

    Parameters
    ----------
    dna : QUEEN
        Template DNA molecule to digest.
    *cutsites : Cutsite or str
        One or more restriction sites. Each element may be a
        :class:`Cutsite` instance or a string key present in
        ``QUEEN.cutsite.lib``.
    selection : {"min", "max"} or int or tuple of int or str, optional
        Rule for selecting a subset or a single digested fragment.

        * ``None`` (default)  
          Return a list of all digested fragments.
        * ``"min"``  
          Return only the shortest fragment.
        * ``"max"``  
          Return only the longest fragment.
        * integer ``n``  
          Return the fragment whose length is closest to ``n`` base pairs.
        * tuple of two integers ``(min_len, max_len)``  
          Return fragment(s) whose lengths fall within the given range.
        * string of the form ``"{qualifier}:{label}"``  
          Return the unique fragment containing a `DNAfeature` whose
          qualifier ``qualifier`` contains the substring ``label``.
        * string of the form ``"!{qualifier}:{label}"``  
          Return the unique fragment that **does not** contain such a
          feature.

        If the selection rule matches multiple fragments when only one is
        expected, a :class:`ValueError` is raised.
    product : str, optional
        Human‑readable name for the selected fragment(s). Stored in the
        resulting `QUEEN` object(s) and recorded in the construction
        history.
    process_name : str, optional
        Short label for this digestion step in the construction history.
        If ``None`` and ``pn`` is also ``None``, a default such as
        ``"Digestion"`` is used.
    process_description : str, optional
        Free‑text description of the digestion step.
    pn : str, optional
        Alias for ``process_name``. Used only when ``process_name`` is
        ``None``.
    pd : str, optional
        Alias for ``process_description``. Used only when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions. Ignored by the current
        implementation but included in the recorded history string.

    Returns
    -------
    QUEEN or list of QUEEN
        If ``selection`` is ``None``, a list of digested `QUEEN` fragments
        is returned. Otherwise, a single `QUEEN` fragment is returned.

    Raises
    ------
    TypeError
        If ``selection`` is a tuple but not of length 2; or if any element
        in ``cutsites`` is neither a :class:`Cutsite` instance nor a
        string key present in ``QUEEN.cutsite.lib``.
    ValueError
        If a label‑based selection string matches zero or multiple
        fragments where a unique fragment is required.

    Examples
    --------
    Basic double digestion and retrieval of all fragments::

        fragments = digestion(
            dna=plasmid,
            *[cutsite_A, cutsite_B]
        )

    Select the fragment closest to 2000 bp::

        fragment = digestion(plasmid, cutsite_A, cutsite_B, selection=2000)

    Select the fragment containing a specific feature label::

        fragment = digestion(plasmid, cutsite_A, cutsite_B,
                             selection="qualifier:promoter")
    """

    if selection is not None:
        if type(selection) == str and selection not in ("min", "max") and bool(re.match(r"^.+:.+$", selection)) == False and bool(re.match(r"^!.+:.+$", selection)) == False:
            if selection.startswith("!"):
                selection = "!label:" + selection[1:]
            else:
                selection = "label:" + selection 
        if (type(selection) not in (int, tuple)) and selection not in ("min", "max") and (bool(re.match(r"^.+:.+$", selection)) == False and bool(re.match(r"^!.+:.+$", selection)) == False):
            raise TypeError("`selection` should be `tuple` value, 'min', 'max', or `str` starting with '!{qualifier_key}:{feature_of_interest}' or '!{qualifier_key}:{feature_of_interest}'.")
    
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
   
    cs_str  = ", ".join(['"{}"'.format(name) for name in cutsite_names])
    kwargs_str = _convert_kwargs(kwargs) 
    
    if type(selection) == tuple:
        qexd = 'digestion(QUEEN.dna_dict["{}"], {}, selection=[{}]{})'.format(dna.project, cs_str, ",".join(map(str, selection)), kwargs_str) 
    elif type(selection) == str:
        qexd = 'digestion(QUEEN.dna_dict["{}"], {}, selection="{}"{})'.format(dna.project, cs_str, selection, kwargs_str) 
    elif type(selection) == int:
        qexd = 'digestion(QUEEN.dna_dict["{}"], {}, selection={}{})'.format(dna.project, cs_str, selection, kwargs_str) 
    else:
        qexd = 'digestion(QUEEN.dna_dict["{}"], {}, selection="max"{})'.format(dna.project, cs_str, kwargs_str)
    process_description = pd if process_description is None else process_description

    new_cutsites = [] 
    for cutsite in cutsites:
        sites = dna.searchsequence(query=cutsite, qexd=True, pn=process_name, pd=process_description)        
        new_cutsites.extend(sites) 
    
    fragments = cutdna(dna, *new_cutsites, qexd=True, product=None, pn=process_name, pd=process_description)

    if len(fragments) == 1 and selection is None: 
        dfragment = _select(fragments, "max") 
    else:
        dfragment = _select(fragments, selection)
     
    if type(dfragment) == list:
        fragments = []
        for d, afragment in enumerate(dfragment):
            if d == len(dfragment) - 1:
                afragment = modifyends(afragment, left="", right="", qexd=qexd, product=product, pn=process_name, pd=process_description) 
            else:
                afragment = modifyends(afragment, left="", right="", qexd=True, product=product, pn=process_name, pd=process_description) 
            fragments.append(afragment) 
        return fragments
    else:
        #if product is not None:
        dfragment = modifyends(dfragment, left="", right="", qexd=qexd, product=product, pn=process_name, pd=process_description) 
        return dfragment 

def ligation(*fragments, unique=True, follow_order=False, auto_select=True, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs): 
    """Simulate ligation of one or more `QUEEN` fragments.

    This function joins multiple `QUEEN` fragments by matching their end
    structures (compatible cohesive or blunt ends). Different permutations
    and orientations are explored unless restricted by the arguments.

    Parameters
    ----------
    *fragments : QUEEN or list of QUEEN
        `QUEEN` fragment(s) to be ligated. Fragments must have compatible
        end structures (generated, for example, by :func:`digestion`,
        :func:`pcr`, or :func:`modifyends`).  Passing a single list or
        tuple of `QUEEN` objects is not allowed; fragments should be
        supplied as positional arguments, e.g. ``ligation(a, b, c)``.
    unique : bool, optional
        If ``True`` (default), require that exactly one construct can be
        assembled from the provided fragments. If multiple distinct
        constructs are possible, a :class:`ValueError` is raised. If
        ``False``, all valid constructs are returned as a list.
    follow_order : bool, optional
        If ``True``, ligation is restricted to the given fragment order
        (no permutations). This is useful when the intended order is known
        a priori. Default is ``False``.
    auto_select : bool, optional
        If ``unique`` is ``False`` and multiple constructs are generated,
        setting ``auto_select=True`` allows the function to pick a single
        construct using internal heuristics (for example, favoring
        constructs without truncated features). When ``auto_select=False``
        and multiple products exist, a list of all products is returned.
    product : str, optional
        Human‑readable name for the ligation product(s). Recorded in the
        returned `QUEEN` object(s) and in the construction history.
    process_name : str, optional
        Short label for this ligation step in the construction history.
        If ``None`` and ``pn`` is also ``None``, a default such as
        ``"Ligation"`` is used.
    process_description : str, optional
        Free‑text description of the ligation step.
    pn : str, optional
        Alias for ``process_name``. Used only when ``process_name`` is
        ``None``.
    pd : str, optional
        Alias for ``process_description``. Used only when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions. Ignored by the current
        implementation but included in the recorded history string.

    Returns
    -------
    QUEEN or list of QUEEN
        If a single unique construct is identified (either because
        ``unique=True`` or because ligation yields only one possible
        product), a single `QUEEN` is returned. Otherwise, a list of
        `QUEEN` objects representing all valid constructs is returned.

    Raises
    ------
    ValueError
        If fragments are passed as a single list or tuple instead of
        individual positional arguments; if no valid ligation products can
        be formed from the provided fragments; or if ``unique=True`` and
        multiple distinct products are possible.
    TypeError
        If the keyword argument ``fragments=...`` is used instead of
        positional arguments; or if any element of ``fragments`` is not a
        `QUEEN` object (for example, if a list of fragments was passed
        without selecting a single fragment from a previous digestion).

    Examples
    --------
    Ligate two compatible fragments into a unique construct::

        construct = ligation(fragment_a, fragment_b, product="joined")

    Allow multiple possible assemblies to be returned::

        constructs = ligation(fragment_a, fragment_b, fragment_c,
                              unique=False, follow_order=False)
    """

    if len(fragments) == 1 and isinstance(fragments[0], (list, tuple)):
        raise ValueError("Fragments must be given as positional arguments, not as a single list. Use ligation(a, b) or ligation(*[a,b]).")
    
    if "fragments" in kwargs:
        raise TypeError('"fragments" is not a valid keyword argument. Pass fragments as positional arguments (a, b). Use ligation(a, b) or ligation(*[a, b])."')

    def add_fragment(fragments, orders, remains, results, flip=1):
        flag = 0 
        fragment1 = fragments[orders[-1][0]] 
        if orders[-1][1] == -1:
            fragment1 = flipdna(fragment1, quinable=False) 
        
        for target in remains:
            fragment2 = fragments[target] 
            
            if flip == 1:
                rl = fragment1._right_end_top * fragment2._left_end_bottom 
                if rl == 1 and fragment1._right_end == fragment2._left_end:
                    flag = 1
                    orders.append((target, 1)) 
                    break
                elif fragment1._right_end_top == 1 and fragment1._right_end_bottom == 1 and fragment2._left_end_top == 1 and fragment2._left_end_bottom == 1:
                    flag = 1 
                    orders.append((target, 1))
                    break
                else:
                    pass 
            else:
                rr = fragment1._right_end_top * fragment2._right_end_top
                if rr == 1 and fragment1._right_end == fragment2._right_end.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]:
                    flag = 1
                    orders.append((target, -1))
                    break
                elif fragment1._right_end_top == 1 and fragment1._right_end_bottom == 1 and fragment2._right_end_top == 1 and fragment2._right_end_bottom == 1:
                    flag = 1 
                    orders.append((target, -1))
                    break
                else:
                    pass 
        
        if flag == 0:
            pass 
        else:
            remains.remove(target) 
            if len(remains) == 0 and len(orders) == len(fragments):
                results.append(orders) 
            else:
                add_fragment(fragments, orders[:], remains[:], results, flip=1)
                add_fragment(fragments, orders[:], remains[:], results, flip=-1)
        return results 

    for fragment in fragments:  
        if type(fragment) == list and type(fragment[0]) == QUEEN:
            raise TypeError("Each QUEEN object should be specified individually, not as a list. Perhaps you forgot to select a single fragment from the digestion results?") 
        else:
            pass 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Ligation"  
    
    kwargs_str = _convert_kwargs(kwargs)
    fragments_str = ", ".join(['QUEEN.dna_dict["{}"]'.format(fragment._product_id) for fragment in fragments])
    
    if unique == True:
        uniquetxt = ""
    else:
        uniquetxt = ", unique={},".format(unique) 

    if follow_order == False or follow_order is None: 
        fotxt = ""
    else:
        fotxt = ", follow_order={}".format(follow_order) 
    
    if follow_order == 'True':
        follow_order = True
    
    if len(fragments) == 1: 
        follow_order = True

    qexd = 'ligation({}, {}{}{})'.format(fragments_str, uniquetxt, fotxt, kwargs_str)
    process_description = pd if process_description is None else process_description
   
    if follow_order == True:
        outobj = joindna(*fragments, topology="circular", autoflip=False, compatibility="complete", qexd=qexd, product=product, pn=process_name, pd=process_description)
        outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragments)
        outobj.printfeature() 
        if len(fragments) == 1:
            if 0 in outobj._positions:
                zero_pos = outobj._positions.index(0)
                outobj   = cutdna(outobj, zero_pos, qexd=True)[0]
                outobj   = joindna(outobj, topology="circular", qexd=True)
                outobj._positions = tuple(range(len(outobj.seq)))
            else:
                pass
        return outobj
    
    else:
        orders   = [(0,1)] 
        results  = [] 
        remains  = list(range(1, len(fragments)))
        results1 = add_fragment(fragments, orders[:], remains[:], results[:], flip=1)
        results2 = add_fragment(fragments, orders[:], remains[:], results[:], flip=-1)
        results  = results1 + results2 
            
    if unique == True:     
        if len(results) == 1:
            orders, flips = list(zip(*results[-1])) 
            fragment_set  = [flipdna(fragments[ind], product=fragments[ind].project, qexd=True, pn=process_name, pd=process_description) if fl == -1 else fragments[ind] for ind, fl in zip(orders, flips)]
            outobj = joindna(*fragment_set, topology="circular", autoflip=False, compatibility="complete", qexd=qexd, product=product, pn=process_name, pd=process_description)
            outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragment_set)
        elif len(results) == 0:
            raise ValueError("The QUEEN_objects cannot be joined due to the end structure incompatibility. Please double-check that you haven't forgotten to perform the restriction enzyme digestion on the input fragments, that the fragments are digested with the appropriate restriction enzymes, and that you are using the correct primers for previous PCRs.") 
        else:
            tf_set = [] 
            new_results = []
            
            for result in results:
                tf_set.append([]) 
                index1, direction1 = result[0] 
                others = result[1:] + [result[0]]  
                if direction1 == 1:
                    fragment1 = fragments[index1]
                else:
                    fragment1 = fragments[index1][::-1]

                for (index2, direction2) in others:
                    if direction2 == 1:
                        fragment2 = fragments[index2]
                    else:
                        fragment2 = fragments[index2][::-1]
                    tf_set[-1].append(check_arrangement(fragment1, fragment2)) 
                    fragment1 = fragment2 
            
            for i, tf in enumerate(tf_set):
                if False in tf:
                    pass 
                else:
                    new_results.append(results[i]) 
            
            if len(new_results) == 1:
                orders, flips = list(zip(*new_results[-1])) 
                fragment_set  = [flipdna(fragments[ind], product=fragments[ind].project, qexd=True, pn=process_name, pd=process_description) if fl == -1 else fragments[ind] for ind, fl in zip(orders, flips)]
                outobj = joindna(*fragment_set, topology="circular", autoflip=False, compatibility="complete", qexd=qexd, product=product, pn=process_name, pd=process_description)
                outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragment_set)
            else:
                raise ValueError("Multiple different constructs will be assembled. You should review your assembly design.")
        
        if len(fragments) == 1:
            if 0 in outobj._positions:
                zero_pos = outobj._positions.index(0)
                outobj = cutdna(outobj, zero_pos, qexd=True)[0]
                outobj = joindna(outobj, topology="circular", qexd=True)
                outobj._positions = tuple(range(len(outobj.seq)))
            else:
                pass
        return outobj
    
    else:
        products = [] 
        for order, flips in indexes_list:
            fragment_set  = [flipdna(fragments[ind], qexd=True, product=fragments[ind].project, pn=process_name, pd=process_description) if fl == -1 else fragments[ind] for ind, fl in zip(orders, flips)]
            outobj = joindna(*fragment_set, topology="circular", autoflip=False, compatibility="complete", qexd=qexd, product=product, pn=process_name, pd=process_description)
            outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragment_set)
            products.append(outobj)
        return products 

def homology_based_assembly(*fragments, mode="gibson", homology_length=15, unique=True, follow_order=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs): #homology_based_assembly
    """Simulate homology‑based DNA assembly (e.g., Gibson or In‑Fusion).

    This function assembles multiple `QUEEN` fragments using overlapping
    homology at their ends. It supports common homology‑based cloning
    schemes such as Gibson Assembly, In‑Fusion cloning, and overlap‑PCR
    style fusion of fragments.

    Parameters
    ----------
    *fragments : QUEEN or list of QUEEN
        `QUEEN` fragment(s) to be assembled. Passing a single list or
        tuple is not allowed; fragments should be supplied as positional
        arguments, e.g. ``homology_based_assembly(a, b, c)``.
    mode : {"gibson", "infusion", "overlappcr"}, optional
        Assembly mode used to interpret and, if necessary, modify fragment
        ends.

        * ``"gibson"`` (default) – double‑stranded overlaps on both ends.
        * ``"infusion"`` – cohesive ends treated according to In‑Fusion
          style overlaps.
        * ``"overlappcr"`` – overlap‑PCR style assembly, typically with
          fragments generated by primers containing homology tails.

    homology_length : int, optional
        Minimum homology length required between adjacent fragments for a
        valid assembly. Default is ``20``.
    unique : bool, optional
        If ``True`` (default), require that exactly one product be formed.
        If multiple distinct constructs can be assembled, a
        :class:`ValueError` is raised. If ``False``, all valid constructs
        are returned.
    follow_order : bool, optional
        If ``True``, assembly is restricted to the given order of
        ``fragments``. If ``False``, permutations may be explored,
        depending on the mode. Default is ``None``, which lets the
        implementation choose a sensible behavior based on ``mode``.
    product : str, optional
        Human‑readable name for the assembled construct(s). Recorded in
        the returned `QUEEN` object(s) and in the construction history.
    process_name : str, optional
        Short label for this assembly step in the construction history.
        If ``None`` and ``pn`` is also ``None``, a default such as
        ``"Homology-based Assembly"`` is used.
    process_description : str, optional
        Free‑text description of the assembly step.
    pn : str, optional
        Alias for ``process_name``. Used when ``process_name`` is
        ``None``.
    pd : str, optional
        Alias for ``process_description``. Used when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions. Ignored by the current
        implementation but included in the recorded history string.

    Returns
    -------
    QUEEN or list of QUEEN
        If a single valid assembly product is obtained (or ``unique=True``),
        a single `QUEEN` object is returned. Otherwise, a list of `QUEEN`
        objects representing all valid assemblies is returned.

    Raises
    ------
    ValueError
        If fragments are passed as a single list or tuple instead of
        positional arguments; if ``mode`` is not one of ``"gibson"``,
        ``"infusion"``, or ``"overlappcr"``; if no valid product can be
        assembled; or if ``unique=True`` and multiple distinct constructs
        are possible.
    TypeError
        If the keyword argument ``fragments=...`` is used instead of
        positional arguments.
    ValueError
        If incompatible end structures are detected among the fragments
        (for example, missing or insufficient homology).

    Examples
    --------
    Assemble two fragments using Gibson Assembly::

        product = homology_based_assembly(
            fragment_a,
            fragment_b,
            mode="gibson",
            homology_length=30,
            product="gibson_product"
        )
    """

    if len(fragments) == 1 and isinstance(fragments[0], (list, tuple)):
        raise ValueError("Fragments must be given as positional arguments, not as a single list. Use homology_based_assembly(a, b, c) or homology_based_assembly(*[a,b,c]).")
    
    if "fragments" in kwargs:
        raise TypeError('"fragments" is not a valid keyword argument. Pass fragments as positional arguments (a, b, c). Use homology_based_assembly(a, b, c) or homology_based_assembly(*[a,b,c])."')

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
    
    #product = product.replace(" ","") if product is not None else None
    kwargs_str = _convert_kwargs(kwargs)
    fragments_str = ", ".join(['QUEEN.dna_dict["{}"]'.format(fragment._product_id) for fragment in fragments])
    
    if homology_length == 15:
        hltxt = ""
    else:
        hltxt = ", homology_length={},".format(homology_length) 

    if unique == True:
        uniquetxt = ""
    else:
        uniquetxt = ", unique={},".format(unique) 

    if follow_order == False or follow_order is None: 
        fotxt = ""
    else:
        fotxt = ", follow_order={}".format(follow_order) 

    if follow_order == 'True':
        follow_order = True

    if len(fragments) == 1:
        follow_order = True

    qexd = 'homology_based_assembly({}, mode="{}"{}{}{}{})'.format(fragments_str, mode, hltxt, uniquetxt, fotxt, kwargs_str)
    process_description = pd if process_description is None else process_description
    
    for fragment in fragments:
        if fragment.topology != "linear": 
            raise ValueError("A 'circular' fragment was detected. All fragments to be assembled should be 'linear' topology.")

    if mode not in ("gibson", "infusion", "overlappcr"):
        raise ValueError("Invalid mode value. The 'mode' variable can only take 'gibson' and 'infusion' values.")

    if mode in ("gibson", "infusion") and follow_order == True and len(fragments) > 1:
        try:
            for f in range(len(fragments)):
                fragment = fragments[f]
                if len(fragment.seq) <= max_homology_length: 
                    mhl = int(len(fragment.seq)) - len(fragment._left_end) - len(fragment._right_end) - 1
                else:
                    mhl = max_homology_length
                fragments[f] = modifyends(fragment, "-{{{}}}/*{{{}}}".format(mhl,mhl), "*{{{}}}/-{{{}}}".format(mhl,mhl), qexd=True, pn=process_name, pd=process_description)

            outobj = joindna(*fragments, autoflip=False, homology_length=homology_length, topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description)
            outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragments)
            if unique == True:
                return outobj
            else:
                return [outobj]
        except Exception:
            pass
    
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
            flip_status_list = [[1 for i in range(len(fragments))]]    
        errors = [] 
        products = [] 
        product_sources = []
        for numset in nums_orders:
            execed = [] 
            for flipset in flip_status_list:
                if tuple([state * -1 for state in flipset]) in execed: 
                    pass 
                else:
                    fragment_set = [fragments[num] if flip == 1 else flipdna(fragments[num], qexd=True, pn=process_name, pd=process_description) for num, flip in zip(numset, flipset)] 
                    for f in range(len(fragment_set)):
                        fragment = fragment_set[f]
                        if len(fragment.seq) <= max_homology_length: 
                            mhl = int(len(fragment.seq)) - len(fragment._left_end) - len(fragment._right_end) - 1
                        else:
                            mhl = max_homology_length
                        fragment_set[f] = modifyends(fragment, "-{{{}}}/*{{{}}}".format(mhl,mhl), "*{{{}}}/-{{{}}}".format(mhl,mhl), qexd=True, pn=process_name, pd=process_description)
                    try:
                        outobj = joindna(*fragment_set, autoflip=False, homology_length=homology_length, topology="linear", qexd=True, product=product, pn=process_name, pd=process_description) 
                        outobj = modifyends(outobj, qexd=True, product=product, pn=process_name, pd=process_description) 
                        outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragment_set)
                        products.append(outobj) 
                        product_sources.append(fragment_set)
                    except Exception as e:
                        errors.append(e) 
                execed.append(flipset) 

    else:
        nums = list(range(len(fragments)))
        if (len(fragments) < 5 and follow_order is None) or follow_order == False:
            nums_orders = list(map(list,it.permutations(nums[:-1])))
            nums_orders = [numlist + [nums[-1]] for numlist in nums_orders]
            flip_status_list = list(it.product(*[[1,-1] for i in range(len(fragments))]))   
        else:
            nums_orders      = [nums]  
            flip_status_list = [[1 for i in range(len(fragments))]]    
        errors = [] 
        products = [] 
        product_sources = []
        for numset in nums_orders:
            execed = [] 
            for flipset in flip_status_list:
                if tuple([state * -1 for state in flipset]) in execed: 
                    pass 
                else:
                    fragment_set = [fragments[num] if flip == 1 else flipdna(fragments[num], qexd=True, pn=process_name, pd=process_description) for num, flip in zip(numset, flipset)]  
                    for f in range(len(fragment_set)):
                        fragment = fragment_set[f]
                        if len(fragment.seq) <= max_homology_length: 
                            mhl = int(len(fragment.seq)) - len(fragment._left_end) - len(fragment._right_end) - 1
                        else:
                            mhl = max_homology_length
                        if mode == "gibson":
                            fragment_set[f] = modifyends(fragment, "-{{{}}}/*{{{}}}".format(mhl,mhl), "*{{{}}}/-{{{}}}".format(mhl,mhl), qexd=True, pn=process_name, pd=process_description)
                        elif mode == "infusion":
                            fragment_set[f] = modifyends(fragment, "*{{{}}}/-{{{}}}".format(mhl,mhl), "-{{{}}}/*{{{}}}".format(mhl,mhl), qexd=True, pn=process_name, pd=process_description)
                    try:
                        outobj = joindna(*fragment_set, autoflip=False, homology_length=homology_length, topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description) 
                        outobj, _ = _rescue_missing_features_by_exact_sequence(outobj, fragment_set)
                        products.append(outobj) 
                        product_sources.append(fragment_set)
                    except Exception as e:
                        errors.append(e) 
                        pass 
                execed.append(flipset) 

    if unique == True:
        if len(products) == 0:
            print(errors) 
            raise ValueError("Error, Incompatible ends were detected. Maybe you need to reflect the PCR primers or restriction enzymes used to generate the fragments.") 
        
        elif len(products) > 1: 
            raise ValueError("Multiple assembled constructs were detected. You should review your assembly design.")
        
        else:
            try:
                product = products[0] 
                if len(fragments) == 1:
                    if 0 in product._positions:
                        zero_pos = product._positions.index(0)
                        product  = cutdna(product, zero_pos, qexd=True)[0]
                        product  = joindna(product, topology="circular", qexd=True)
                        product._positions = tuple(range(len(product.seq)))
                    else:
                        pass
                
                product, _ = _rescue_missing_features_by_exact_sequence(product, fragments)
                return product
            
            except Exception as e:
                print(e, errors) 
                
    else: 
        return products 

def annealing(ssdna1, ssdna2, homology_length=4, product=None, pn=None, pd=None, process_name=None, process_description=None, **kwargs):
    """Simulate annealing of two complementary single‑stranded DNAs.

    The function joins two single‑stranded `QUEEN` objects (or DNA strings)
    into a double‑stranded DNA molecule based on a homologous overlap.

    Parameters
    ----------
    ssdna1 : QUEEN or str
        First single‑stranded DNA. If a dsDNA `QUEEN` is provided, its top
        strand is used. If a string is provided, it is internally
        converted to a ssDNA `QUEEN`.
    ssdna2 : QUEEN or str
        Second single‑stranded DNA. Same conventions as for ``ssdna1``.
    homology_length : int, optional
        Minimum length of the homologous region required for annealing.
        Default is ``4``.
    product : str, optional
        Human‑readable name for the annealed product. Recorded in the
        returned `QUEEN` object and in the construction history.
    process_name : str, optional
        Short label for this annealing step in the construction history.
        If ``None`` and ``pn`` is also ``None``, a default such as
        ``"Annealing"`` is used.
    process_description : str, optional
        Free‑text description of the annealing step.
    pn : str, optional
        Alias for ``process_name`` when ``process_name`` is ``None``.
    pd : str, optional
        Alias for ``process_description`` when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions.

    Returns
    -------
    QUEEN
        Double‑stranded `QUEEN` object representing the annealed duplex.

    Raises
    ------
    TypeError
        If ``ssdna1`` or ``ssdna2`` is neither a `QUEEN` object nor a
        string.
    ValueError
        If no suitable homologous region of length at least
        ``homology_length`` is found.

    Examples
    --------
    Anneal two complementary oligonucleotides into a short duplex::

        duplex = annealing(
            ssdna1=QUEEN(seq="ACGTACGT", ssdna=True),
            ssdna2=QUEEN(seq="ACGTACGT", ssdna=True),
            homology_length=8
        )
    """
    if type(ssdna1) == str:
        ssdna1 = QUEEN(seq=ssdna1, ssdna=True)
    
    if type(ssdna2) == str:
        ssdna2 = QUEEN(seq=ssdna2, ssdna=True)
        
    if type(ssdna1) != QUEEN:
        TypeError("`ssdna1` must be a QUEEN object or a string.")

    if type(ssdna2) != QUEEN:
        TypeError("`ssdna2` must be a QUEEN object or a string.")

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Annealing" 
    
    #product = product.replace(" ","") if product is not None else None
    kwargs_str = _convert_kwargs(kwargs)
    if homology_length == 4:
        hltxt = ""
    else:
        hltxt = ", homology_length={}".format(homology_length) 

    def _annealing_operand_text(ssdna):
        product_id = getattr(ssdna, "_product_id", None)
        if product_id:
            return 'QUEEN.dna_dict["{}"]'.format(product_id)
        return "QUEEN(seq={}, ssdna=True)".format(repr(str(ssdna.seq)))

    qexd = 'annealing({}, {}{}{})'.format(
        _annealing_operand_text(ssdna1),
        _annealing_operand_text(ssdna2),
        hltxt,
        kwargs_str,
    )
    process_description = pd if process_description is None else process_description
 
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
    
    annealed_dna = joindna(ssdna1, ssdna2, homology_length=homology_length, qexd=qexd, product=product, pn=process_name, pd=process_description)
    ssdna1._ssdna = False if flag1 == 1 else True 
    ssdna2._ssdna = False if flag2 == 1 else True

    return annealed_dna 

def gateway_reaction(destination, entry, mode="BP", product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """Simulate a Gateway-style recombination reaction.

    This function performs a Gateway‑like reaction between a destination
    backbone and an entry vector containing compatible recombination
    sites (att sites). Basic BP and LR‑style reactions are supported.

    Parameters
    ----------
    destination : QUEEN
        Circular `QUEEN` object representing the destination backbone.
    entry : QUEEN
        `QUEEN` object representing the entry construct containing the
        insert and att sites.
    mode : str or tuple of str, optional
        Reaction mode. Common values include:

        * ``"BP"`` – entry attL / attR to destination attP / attB
          style reaction.
        * ``"LR"`` – entry attB / attP to destination attL / attR
          style reaction.

        Internally, the mode may be treated as a tuple indicating the
        specific att site flavors; see implementation for details.
    product : str, optional
        Human‑readable name for the recombination product. Recorded in the
        returned `QUEEN` object and in the construction history.
    process_name : str, optional
        Short label for this Gateway step. If ``None`` and ``pn`` is also
        ``None``, the default ``"Gateway Reaction"`` is used.
    process_description : str, optional
        Free‑text description of the reaction.
    pn : str, optional
        Alias for ``process_name`` when ``process_name`` is ``None``.
    pd : str, optional
        Alias for ``process_description`` when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions.

    Returns
    -------
    QUEEN
        `QUEEN` object representing the recombination product.

    Raises
    ------
    TypeError
        If ``destination`` is not a circular `QUEEN` object.
    ValueError
        If the required att sites cannot be found or are found more than
        once in either ``destination`` or ``entry``; or if the specified
        ``mode`` is not supported by the implementation.

    Examples
    --------
    Perform a BP‑style reaction between a destination and an entry
    construct::

        product = gateway_reaction(
            destination=dest_plasmid,
            entry=entry_plasmid,
            mode="BP",
            product="bp_product"
        )
    """ 
    if type(destination) == QUEEN:
        if destination.topology == "circular":
            pass 
        else:
            raise TypeError("`destination` must be a QUEEN object with circular sequence topology.") 
    else:
        raise TypeError("`destination` must be a QUEEN object with circular sequence topology.") 

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Gateway Reaction" 
    
    #product = product.replace(" ","") if product is not None else None
    kwargs_str = _convert_kwargs(kwargs)
    qexd = 'gateway_reaction(QUEEN.dna_dict["{}"], QUEEN.dna_dict["{}"], mode="{}"{})'.format(destination._product_id, entry._product_id, mode, kwargs_str)
    process_description = pd if process_description is None else process_description

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

    attx1 = entry.searchsequence(cs.lib["attX1"], product="att{}1_site".format(mode[0]), qexd=True, pn=process_name, pd=process_description) 
    attx2 = entry.searchsequence(cs.lib["attX2"], product="att{}2_site".format(mode[0]), qexd=True, pn=process_name, pd=process_description)
    atty1 = destination.searchsequence(cs.lib["attY1"], product="att{}1_site".format(mode[1]), qexd=True, pn=process_name, pd=process_description)
    atty2 = destination.searchsequence(cs.lib["attY2"], product="att{}2_site".format(mode[1]), qexd=True, pn=process_name, pd=process_description) 
    if len(attx1) > 1:
        raise ValueError("Multiple att{}1 sites were detected.".format(mode[0]))
    elif len(attx1) == 1:
        attx1 = attx1[0] 
    else:
        raise ValueError("No att{}1 site was detected.".format(mode[0]))


    if len(attx2) > 1:
        raise ValueError("Multiple att{}2 sites were detected.".format(mode[0]))
    elif len(attx2) == 1:
        attx2 = attx2[0] 
    else:
        raise ValueError("No att{}2 site was detected.".format(mode[0]))
       
    if len(atty1) > 1:
        raise ValueError("Multiple att{}1 sites were detected.".format(mode[1]))
    elif len(atty1) == 1:
        atty1 = atty1[0] 
    else:
        raise ValueError("No att{}1 site was detected.".format(mode[1]))

    if len(atty2) > 1:
        raise ValueError("Multiple att{}2 sites were detected.".format(mode[1]))
    elif len(atty2) == 1:
        atty2 = atty2[0]
    else:
        raise ValueError("No att{}2 site was detected.".format(mode[1]))

    def _rc(seq):
        table = str.maketrans("ATGCRYKMSWBDHV", "TACGYRMKWSVHDB")
        return seq.translate(table)[::-1]

    insert = None
    if attx1.strand == 1:
        insert = cropdna(entry, attx1, attx2, qexd=True, pn=process_name, pd=process_description)
    elif attx1.strand == -1:
        insert = cropdna(entry, attx2, attx1, qexd=True, pn=process_name, pd=process_description)

    if insert is None:
        raise ValueError("Failed to crop the Gateway entry insert from the attX sites. Check the att-site orientation on the entry construct.")

    destination_candidates = []
    for dest_obj in (destination, flipdna(destination, quinable=0)):
        cand_y1 = dest_obj.searchsequence(cs.lib["attY1"], quinable=False)
        cand_y2 = dest_obj.searchsequence(cs.lib["attY2"], quinable=False)
        if len(cand_y1) != 1 or len(cand_y2) != 1:
            continue
        cand_y1 = cand_y1[0]
        cand_y2 = cand_y2[0]
        if cand_y1.strand == 1:
            destination_candidates.append(cropdna(dest_obj, cand_y2, cand_y1, qexd=True, pn=process_name, pd=process_description))
        elif cand_y1.strand == -1:
            destination_candidates.append(cropdna(dest_obj, cand_y1, cand_y2, qexd=True, pn=process_name, pd=process_description))

    if len(destination_candidates) == 0:
        raise ValueError("Failed to crop the Gateway destination backbone from the attY sites. Check the att-site orientation on the destination construct.")

    if mode == "BP":
        attl1_seq = "CCAACTTTGTACAAAAAAGCAGGCT"
        attl2_seq = "ACCCAGCTTTCTTGTACAAAGTTGG"
        core_left = len(insert._left_end) if len(insert._left_end) > 0 else 0
        core_right = len(insert.seq) - len(insert._right_end) if len(insert._right_end) > 0 else len(insert.seq)
        if core_right > core_left and (core_left > 0 or core_right < len(insert.seq)):
            entry_core = cropdna(insert, core_left, core_right, qexd=True, pn=process_name, pd=process_description)
        else:
            entry_core = insert
        gateway_insert = modifyends(entry_core, left=attl1_seq, right=attl2_seq, qexd=True, pn=process_name, pd=process_description)
        bp_products = []
        cs.lib["attL1"] = "CCAACTTT^GTACAAA_AAAGCAGGCT"
        cs.lib["attL2"] = "ACCCAGCTTT^CTTGTAC_AAAGTTGG"
        for destination_crop in destination_candidates:
            left_ovhg = len(destination_crop._left_end) if len(destination_crop._left_end) > 0 else 0
            right_ovhg = len(destination_crop._right_end) if len(destination_crop._right_end) > 0 else 0
            crop_end = len(destination_crop.seq) - right_ovhg
            if crop_end <= left_ovhg:
                continue
            backbone_internal = cropdna(destination_crop, left_ovhg, crop_end, qexd=True, pn=process_name, pd=process_description)
            try:
                product_obj = joindna(gateway_insert, backbone_internal, topology="circular", qexparam=qexd, product=product, pn=process_name, pd=process_description)
            except Exception:
                continue
            hits_l1 = product_obj.searchsequence(query=cs.lib["attL1"].cutsite, quinable=False)
            hits_l2 = product_obj.searchsequence(query=cs.lib["attL2"].cutsite, quinable=False)
            if len(hits_l1) == 1 and len(hits_l2) == 1:
                bp_products.append(product_obj)
        if len(bp_products) > 0:
            bp_products.sort(key=lambda obj: len(obj.seq), reverse=True)
            return bp_products[0]

    insert_candidates = [insert]
    if mode == "BP":
        normalized_insert = copy.deepcopy(insert)
        if len(normalized_insert._left_end) > 0 and attx1.strand == -1:
            normalized_insert._left_end = _rc(normalized_insert._left_end)
        if len(normalized_insert._right_end) > 0 and attx2.strand == -1:
            normalized_insert._right_end = _rc(normalized_insert._right_end)
        insert_candidates.append(normalized_insert)

    for insert_cand in insert_candidates:
        for destination_cand in destination_candidates:
            try:
                return joindna(insert_cand, destination_cand, topology="circular", compatibility="complete", autoflip=False, qexparam=qexd, product=product, pn=process_name, pd=process_description)
            except Exception:
                pass

    destination_crop = destination_candidates[0]
    outobj = ligation(insert_candidates[-1], destination_crop, qexd=True, pn=process_name, pd=process_description) 
    outobj = modifyends(outobj, left="", right="", qexd=qexd, product=product, pn=process_name, pd=process_description)
    return outobj

def goldengate_assembly(destination, entry, cutsite=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """Simulate Golden Gate Assembly.

    This function performs a Golden Gate Assembly between a circular
    destination backbone and one or more entry fragments using a type IIS
    restriction site.

    Parameters
    ----------
    destination : QUEEN
        Circular `QUEEN` object representing the destination backbone.
        The sequence topology must be ``"circular"``.
    entry : list of QUEEN
        List of `QUEEN` objects representing the insert fragment(s).
    cutsite : Cutsite or str, optional
        Type IIS restriction site used in the Golden Gate reaction.
        May be a :class:`Cutsite` object or a string key present in
        ``QUEEN.cutsite.lib``. If ``None``, a suitable default may be
        chosen by the implementation.
    product : str, optional
        Human‑readable name for the assembled construct. Recorded in the
        returned `QUEEN` object and in the construction history.
    process_name : str, optional
        Short label for this Golden Gate step in the construction history.
        If ``None`` and ``pn`` is also ``None``, the default
        ``"Golden Gate Assembly"`` is used.
    process_description : str, optional
        Free‑text description of the Golden Gate step.
    pn : str, optional
        Alias for ``process_name`` when ``process_name`` is ``None``.
    pd : str, optional
        Alias for ``process_description`` when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions.

    Returns
    -------
    QUEEN
        `QUEEN` object representing the assembled Golden Gate construct.

    Raises
    ------
    TypeError
        If ``destination`` is not a `QUEEN` object with circular topology;
        or if ``entry`` is not a list or tuple of `QUEEN` objects.
    ValueError
        If ``cutsite`` is neither a :class:`Cutsite` instance nor a string
        key in ``QUEEN.cutsite.lib``, or if no valid assembly can be
        generated using the given cut site and fragment set.

    Examples
    --------
    Assemble one insert into a backbone using a type IIS site::

        construct = goldengate_assembly(
            destination=backbone,
            entry=[insert],
            cutsite="BsaI",
            product="gg_product"
        )
    """
    
    if type(destination) == QUEEN:
        if destination.topology == "circular":
            pass 
        else:
            raise TypeError("`destination` must be a QUEEN object with circular sequence topology. You do not need to process the `digestion`.") 
    else:
        raise TypeError("`destination` must be a QUEEN object with circular sequence topology.") 
    
    if type(entry) in (tuple, list):
        pass 
    else:
        raise TypeError("`entry` must be a list composed of QUEEN object(s)")  

    if type(cutsite) == str and cutsite in cs.lib.keys():
        cutsite = cs.lib[cutsite]  
    elif type(cutsite) == Cutsite or "cutsite" in cutsite.__dict__:
        pass 
    else:
        raise ValueError("No valide restriction site was specified.")

    process_name = pn if process_name is None else process_name
    if process_name is None:
        process_name = "Golden Gate Assembly" 
    
    #product = product.replace(" ","") if product is not None else None
    kwargs_str = _convert_kwargs(kwargs)
    entry_str  = ", ".join(['QUEEN.dna_dict["{}"]'.format(aentry._product_id) for aentry in entry])
    entry_str  = "[{}]".format(entry_str)
    qexd = 'goldengate_assembly(QUEEN.dna_dict["{}"], {}, cutsite="{}"{})'.format(destination._product_id, entry_str, cutsite.name, kwargs_str)
    process_description = pd if process_description is None else process_description

    if type(entry) == QUEEN:
        entry = [entry]
    
    fragments = [] 
    for aentry in entry:
        # Pre-digested linear entry fragments may still carry terminal cutsite
        # annotations. For Golden Gate routing, what matters here is whether the
        # fragment sequence still contains an internal recognition site, not
        # whether a cutsite feature annotation is present.
        has_internal_cutsite = (cutsite.seq in aentry.seq) or (cutsite.rcseq in aentry.seq)
        if aentry.topology == "linear" and has_internal_cutsite is False:
            insert = aentry
        else: 
            inserts = digestion(aentry, cutsite, qexd=True, product=aentry.project, pn=process_name, pd=process_description)
            if type(inserts) != list:
                inserts = [inserts]
            for insert in inserts:
                if cutsite.seq in insert.seq or cutsite.rcseq in insert.seq:
                    pass
                else:
                    break
        
        fragments.append(insert) 

    backbones =  digestion(destination, cutsite, qexd=True, product=destination.project, pn=process_name, pd=process_description) 
    if type(backbones) != list:
        backbones = [backbones]
    for backbone in backbones:
        if cutsite.seq in backbone.seq or cutsite.rcseq in backbone.seq:
            pass
        else:
            break
    
    fragments.append(backbone)
    outobj = ligation(*fragments, qexd=True, pn=process_name, pd=process_description) 
    outobj = modifyends(outobj, left="", right="", qexd=qexd, product=product, pn=process_name, pd=process_description)
    return outobj

def topo_cloning(destination, entry, mode="TA", product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs):
    """Simulate TOPO cloning.

    This function models TOPO‑style cloning of an insert fragment into a
    destination backbone. Different modes correspond to TA cloning, blunt
    cloning, or directional cloning using a specific 5′ overhang.

    Parameters
    ----------
    destination : QUEEN
        `QUEEN` object representing the backbone molecule. Depending on
        ``mode``, its ends may be automatically processed to carry the
        required TOPO overhangs.
    entry : QUEEN
        Insert fragment to be cloned, typically a linear amplicon with
        appropriate end structures:

        * ``mode="TA"`` – insert should carry a 3′ A overhang on each
          end (typical of some polymerases).
        * ``mode="blunt"`` – insert should be blunt‑ended.
        * ``mode="directional"`` – insert should carry a specific 5′ end
          sequence (for example, a short directional tag).

    mode : {"TA", "blunt", "directional"}, optional
        Cloning mode. Default is ``"TA"``.
    product : str, optional
        Human‑readable name for the TOPO cloning product. Recorded in the
        returned `QUEEN` object and in the construction history.
    process_name : str, optional
        Short label for this TOPO step. If ``None`` and ``pn`` is also
        ``None``, the default such as ``"TOPO Cloning"`` is used.
    process_description : str, optional
        Free‑text description of the TOPO cloning step.
    pn : str, optional
        Alias for ``process_name`` when ``process_name`` is ``None``.
    pd : str, optional
        Alias for ``process_description`` when
        ``process_description`` is ``None``.
    **kwargs
        Reserved for future extensions.

    Returns
    -------
    QUEEN
        `QUEEN` object representing the TOPO cloning product.

    Raises
    ------
    TypeError
        If ``destination`` or ``entry`` is not a `QUEEN` object.
    ValueError
        If ``mode`` is not one of ``"TA"``, ``"blunt"``, or
        ``"directional"``; or if the end structures of ``destination`` and
        ``entry`` are incompatible with the requested mode (for example,
        missing the required 3′ A overhangs for ``"TA"`` cloning).

    Examples
    --------
    Clone a blunt amplicon into a TOPO backbone::

        product = topo_cloning(
            destination=backbone,
            entry=insert,
            mode="blunt",
            product="topo_product"
        )
    """
    if type(destination) != QUEEN: 
        raise TypeError("`destination` object must be instance of QUEEN object.")  
    
    if type(entry) != QUEEN: 
        raise TypeError("`entry` object must be instance of QUEEN object.")  

    kwargs_str = _convert_kwargs(kwargs)
    qexd = 'topo_cloning(QUEEN.dna_dict["{}"], QUEEN.dna_dict["{}"], mode="{}"{})'.format(destination._product_id, entry._product_id, mode, kwargs_str)

    if mode == "TA":
        if destination.topology == "circular":
            destination = digestion(destination, "AflII", selection="max", qexd=True, product=destination.project, pn=process_name, pd=process_description)
            destination = modifyends(destination, qexd=True, pn=process_name, pd=process_description) 
            destination = modifyends(destination, "-/*", "*/-", qexd=True, pn=process_name, pd=process_description) 
        else: 
            if destination._left_end_top == 1 and destination._left_end_bottom == 1 and destination._right_end_top == 1 and destination._right_end_bottom == 1:
                if destination.seq[0] == "A" and  destination.seq[-1] == "T":
                    destination = flipdna(destination, qexd=True, pn=process_name, pd=process_description)
                elif destination.seq[0] == "T" and  destination.seq[-1] == "A":
                    pass
                else:
                    raise ValueError("Incompatible end structures for 'TA' cloning: the linear destination must be blunt-ended and begin with T and end with A after opening the TOPO backbone.")
                destination = modifyends(destination, "-/*", "*/-", qexd=True, pn=process_name, pd=process_description)
 
        entry  = modifyends(entry, "A", "T", qexd=True, pn=process_name, pd=process_description) 
        entry  = modifyends(entry, "-/*", "*/-", qexd=True, pn=process_name, pd=process_description)
        outobj = joindna(destination, entry, autoflip=False, compatibility="complete", homology_length=1, topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description)

    elif mode == "blunt": 
        if destination.topology == "circular":
            destination = digestion(destination, "AflII", selection="max", qexd=True, product=destination.project, pn=process_name, pd=process_description)
            destination = modifyends(destination, qexd=True, pn=process_name, pd=process_description) 
        outobj = joindna(destination, entry, autoflip=False, compatibility="complete", topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description)

    elif mode == "directional":
        if destination.topology == "circular":
            destination = digestion(destination, "StyI", selection="max", qexd=True, product=destination.project, pn=process_name, pd=process_description)
            destination = modifyends(destination, qexd=True, pn=process_name, pd=process_description)
            destination = cropdna(destination, 1, len(destination.seq)-3, qexd=True, pn=process_name, pd=process_description)
            destination = modifyends(destination, "*/*", "----/****", qexd=True, pn=process_name, pd=process_description)
        else: 
            if destination._left_end_top == 1 and destination._left_end_bottom == 1 and destination._right_end_top == 1 and destination._right_end_bottom == 1:
                if destination.seq[-4:0] == "CACC":
                    destination = modifyends(destination, "*/*", "----/****", qexd=True) 
                elif destination.seq[0:4] == "GGTG":
                    destination = modifyends(destination, "****/----", "*/*", qexd=True) 
                else:
                    pass 

        entry = modifyends(entry, "****/----", "*/*", qexd=True, pn=process_name, pd=process_description)
        outobj = joindna(destination, entry, autoflip=False, compatibility="complete", topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description) 
    
    else:
        allowed = ("TA", "blunt", "directional")
        raise ValueError(f"Invalid cloning mode: {mode!r}. Allowed modes: {', '.join(map(repr, allowed))}.")
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
    
    kwargs_str = _convert_kwargs(kwargs)
    qexd = "intra_site_specific_recombination({}, site={}{})".format(dna.project, site, kwargs_str)
    process_description = pd if process_description is None else process_description
    
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

    recsites    = dna.searchsequence(cs.lib["recsite"], qexd=True, product=site, pn=process_name, pd=process_description)
    outobj_list = [] 
    for combi in it.combinations(recsites, 2):
        recsite1 = combi[0]
        recsite2 = combi[1]
        fragments = cutdna(dna, recsite1, recsite2, qexd=True, pn=process_name, pd=proces_description)
        
        if recsite1.strand == recsite2.strand:
            if dna.topology == "circular":
                fragment = _select(fragments, selection="max") 
                outobj = joindna(fragment, autoflip=False, compatibility="complete", topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description)
            else:
                outobj = joindna(fragment[0], fragment[2], autoflip=False, compatibility="complete", qexd=qexd, product=product, pn=process_name, pd=process_description) 

        else:
            if dna.topology == "circular":
                fragment1 = _select(fragments, selection="max") 
                fragment2 = _select(fragments, selection="min") 
                fragment2 = flipdna(fragment2, pn=process_name, qexd=True, pd=process_description)
                outobj = joindna(fragment1, fragment2, autoflip=False, compatibility="complete", topology="circular", qexd=qexd, product=product, pn=process_name, pd=process_description)
            else:
                reversed_fragment = flipdna(fragment[1], qexd=True, pn=process_name, pd=process_description) 
                outobj = joindna(fragment[0], reversed_fragment, fragment[2], autoflip=False, compatibility="complete", qexd=qexd, product=product, pn=process_name, pd=process_description) 
        outobj_list.append(outobj)
    
    if len(outobj_list) == 1:
        return outobj_list[0]
    else:
        return outobj_list 

def homologous_recombination(donor, entry, left_homology=None, right_homology=None, homology_length=100, product=None, process_name=None, pn=None, process_description=None, pd=None):
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
    
    #product = product.replace(" ","") if product is not None else None
    kwargs_str = _convert_kwargs(kwargs)
    qexd = 'homologous_recombination(QUEEN.dna_dict["{}"], QUEEN.dna_dict["{}"]{})'.format(destination._product_id, entry._product_id, mode, destination_selection, entry_selection, kwargs_str)
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
        dresultl = destination.searchseqeunce(query=left_homology, unique=True, qexd=True, pn=process_name, pd=process_description)[0] 
        dresultr = destination.searchseqeunce(query=right_homology, unique=True, qexd=True, pn=process_name, pd=process_description)[0] 
        if dresultl.strand == -1 or dresultr.strand == -1:
            raise ValueError("`left_homology` and `right_homology` sequeneces should be on the top strand.") 
        dstrand = dresultl.strand 
        dtarget = dresultl.end

        results  = entry.searchsequence(query=left_homology + ".+" + right_homology, qexd=True, pn=process_name, pd=process_description)
        if len(results) == 1:
            flag = 1

    if flag == 0:
        raise ValueError("Any proper homology arms were not detected.")
    
    eresult = results[0]
    estrand = eresult.strand 
        
    left_dest  = cropdna(destination, 0, dtarget, qexd=True, pn=process_name, pd=process_description)
    right_dest = cropdna(destination, dtarget, len(destination.seq), qexd=True, pn=process_name, pd=process_description) 
    if estrand == dstrand:
        estart  = eresult.start + len(left_homology)
        eend    = eresult.end - len(right_homology) 
        insert = cropdna(entry, estart, eend, qexd=True, pn=process_name, pd=process_description) 
    else:
        estart  = eresult.start + len(right_homology)
        eend    = eresult.end - len(left_homology) 
        insert = flipdna(cropdna(entry, estart, eend, qexd=True, pn=process_name, pd=process_description), qexd=True, pn=process_name, pd=process_description)

    if destination.topology == "circular":
        outobj = joindna(left_dest, insert, right_dest, topology="circular", autoflip=False, qexd=qexd, product=product, pn=process_name, pd=process_description)
    else:
        outobj = joindna(left_dest, insert, right_dest, autoflio=False, qexd=qexd, product=product, pn=process_name, pd=process_description) 
    return outobj

def check_arrangement(fragment1, fragment2):
    features1 = [feat for feat in fragment1.dnafeatures if feat.feature_type not in ("source", "primer", "primer_bind")] 
    features2 = [feat for feat in fragment2.dnafeatures if feat.feature_type not in ("source", "primer", "primer_bind")]
    feat1 = features1[-1]
    feat2 = features2[0] 
    if feat1.feature_type == "promoter" and feat2.feature_type == "CDS":
        if feat1.strand == 1 and feat2.strand == -1:
            return False
        else:
            pass
    
    if feat1.feature_type == "CDS" and feat2.feature_type == "promoter":
        if feat1.strand == -1 and feat2.strand == 1:
            return False
        else:
            pass 
    
    if feat1.feature_type == "CDS" and feat2.feature_type == "CDS":
        if feat1.strand == feat2.strand:
            pass 
        else:
            return False
    else:
        pass 
    return True 

def Tm_NN(check=True, strict=True, nn_table=None, tmm_table=None, imm_table=None, de_table=None, dnac1=25, dnac2=25, selfcomp=False, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5): 
    return functools.partial(mt.Tm_NN, check=check, strict=strict, nn_table=nn_table, tmm_table=tmm_table, imm_table=imm_table, de_table=de_table, dnac1=dnac1, dnac2=dnac2, selfcomp=selfcomp, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)  

def primerdesign(template, target, fw_primer=None, rv_primer=None, fw_margin=0, rv_margin=0,
                 adapter_mode="standard", fw_adapter=None, rv_adapter=None, fw_partner=None, 
                 rv_partner=None, requirement=None, fw_name="fw_primer", rv_name="rv_primer",
                 mut_pattern=None, target_tm=60.0, nonspecific_limit=3, auto_adjust=1, 
                 homology_length=30, tm_func=None, primer_length=(16, 25), design_num=1,
                 gap=None, batch_process=False, product=None, process_name=None,
                 process_description=None, pn=None, pd=None):
    
    """
    Design PCR primers for a specified target region.

    This function designs forward and reverse primers to amplify a ``target``
    region from a ``template`` `QUEEN` object. Primer design can be constrained
    by target melting temperature (Tm), primer length range, and user-defined
    filters, and can optionally encode site-directed mutagenesis. It also
    supports batch design across multiple template/target pairs and provides
    utilities for homology-based cloning workflows (e.g., Gibson/In-Fusion/
    overlap-PCR style assemblies).

    Parameters
    ----------
    template : QUEEN or sequence of QUEEN
        PCR template DNA as a `QUEEN` object. If a list/tuple of templates is
        provided, primer design is performed in batch mode (one design per
        template). In batch mode, list-valued parameters must either match the
        length of ``template`` or be provided as scalars (scalars are broadcast).
    target : QUEEN or sequence of QUEEN
        Sub-region that must be included in the PCR amplicon. In batch mode,
        ``target`` must be a list/tuple whose length matches ``template``, and
        each element corresponds to the template at the same index.

        The function requires that the target sequence is contained in the
        template sequence. If ``target.seq`` is not found in ``template.seq`` but
        ``target.rcseq`` is found, the target may be flipped internally to match
        the template orientation. If neither is found, a ``ValueError`` is raised.
    fw_primer : QUEEN (ssDNA recommended) or str or sequence, optional
        Forward primer to use instead of designing one. If a string is provided,
        it is interpreted as a DNA sequence and may be converted to an ssDNA
        `QUEEN`. In batch mode, provide a sequence aligned to ``template``/``target``
        or a scalar to broadcast.
    rv_primer : QUEEN (ssDNA recommended) or str or sequence, optional
        Reverse primer to use instead of designing one. Same conventions as
        ``fw_primer``.
    fw_margin : int or "auto" or sequence of int/"auto", optional
        Additional bases to include upstream (5′ side) of the target region when
        choosing forward primer binding sites. Default is ``0``.
        If set to ``"auto"``, the function starts from ``0`` and increases the
        forward margin in small steps until unique primer candidates are found
        or the internal auto-margin limit is reached.
    rv_margin : int or "auto" or sequence of int/"auto", optional
        Additional bases to include downstream (3′ side) of the target region when
        choosing reverse primer binding sites. Default is ``0``.
        If set to ``"auto"``, the function starts from ``0`` and increases the
        reverse margin in small steps until unique primer candidates are found
        or the internal auto-margin limit is reached.
    adapter_mode : {"standard", "gibson", "infusion", "overlappcr", "RE"}, optional
        Specifies how ``fw_adapter`` / ``rv_adapter`` are interpreted and how
        partner-derived overlaps are constructed.

        - ``"standard"`` (default):
          Adapters are simply prepended to designed primers.
        - ``"gibson"``, ``"infusion"``, ``"overlappcr"``:
          Adapters and/or partner-derived sequences are treated as homology tails
          for homology-based assembly. In many workflows, these modes are mainly
          used in batch designs where multiple fragments are intended to be assembled
          together, and the mode is expected to be common across the batch.
        - ``"RE"``:
          Adapters represent restriction site logic. When partner-derived ends are
          needed, partners are expected to be digested `QUEEN` objects so the end
          structure can be inferred from digestion history.

        Note
            Some versions/uses treat ``"gibson"``, ``"infusion"``, and ``"overlappcr"``
            as intended primarily for batch workflows. If you use them in single-target
            mode, explicitly provide partners and verify the generated tails.
    fw_adapter : QUEEN or str or Cutsite or sequence, optional
        Adapter to prepend to the 5′ end of designed forward primers.

        - If a `QUEEN` or DNA string is provided, it is prepended as sequence.
          IUPAC bases may be accepted depending on implementation.
        - If a `Cutsite` object (or a cutsite-name key) is provided, the corresponding
          restriction site sequence is prepended using the internal cutsite wrapper.
        - If ``adapter_mode`` is one of ``{"gibson","infusion","overlappcr","RE"}``
          and ``fw_partner`` is provided, the partner segment immediately upstream of
          the target junction may be automatically joined with ``fw_adapter`` to form
          the final tail that is prepended to the forward primer.
    rv_adapter : QUEEN or str or Cutsite or sequence, optional
        Adapter to prepend to the 5′ end of designed reverse primers (the reverse-
        complemented adapter is prepended to the reverse primer sequence).

        Behavior parallels ``fw_adapter``. If ``adapter_mode`` is one of
        ``{"gibson","infusion","overlappcr","RE"}`` and ``rv_partner`` is provided,
        partner context may be used to derive the final tail.
    fw_partner : QUEEN or str or sequence, optional
        DNA sequence (or `QUEEN`) that will join to the 5′ end of the PCR product
        via ligation or homology-based assembly. In homology-based modes, this partner
        can be used to derive overlap/homology tails. For some modes (e.g., ``"RE"``),
        partners are expected to be dsDNA `QUEEN` objects with appropriate end structures.
    rv_partner : QUEEN or str or sequence, optional
        DNA sequence (or `QUEEN`) that will join to the 3′ end of the PCR product.
        Same conventions as ``fw_partner``.
    requirement : callable or sequence of callable, optional
        Filter(s) applied to candidate primer pairs. Each callable should accept a
        dictionary describing a primer pair and return ``True`` if the pair is acceptable.
        The dictionary typically contains keys like ``"fw"``, ``"rv"``, ``"fw_tm"``,
        and ``"rv_tm"``.

        If ``None`` (default), a built-in default filter is applied equivalent to:

        - the last base of each primer is not ``"A"`` or ``"T"``
        - neither primer contains any 4-base homopolymer runs:
          no ``"AAAA"``, ``"TTTT"``, ``"GGGG"``, or ``"CCCC"``
    fw_name : str or sequence of str, optional
        Forward primer name(s)/label(s). In batch mode, provide a list aligned to
        ``target`` or a scalar to broadcast. If not provided (or if explicitly set
        to ``None`` in some workflows), implementations may name primers as
        ``fw_primer{num}``, where ``num`` is the batch index.
    rv_name : str or sequence of str, optional
        Reverse primer name(s)/label(s). Same conventions as ``fw_name``. If not
        provided (or set to ``None``), implementations may name primers as
        ``rv_primer{num}``.
    mut_pattern : dict (MutSpec) or list of dict (MutSpec), optional
        Site-directed mutagenesis specification(s) within ``target``. Coordinates are
        relative to ``target`` (5′→3′, 0-based, half-open [start, end)) unless
        ``relative`` is specified. Multiple edit sites can be represented by supplying
        list values (vectorized fields) inside one MutSpec dict.

        Defaults
            If ``mut_pattern`` is ``None`` (default), no mutagenesis is applied and
            primers are designed for the unmodified target region.

        MutSpec keys (dict)
            operation : {"Q5","QuickChange","gibson","infusion","overlappcr"}, optional
                Mutagenesis strategy. Currently implemented mutagenesis support is
                limited to a single scalar edit with ``template == target``. In that
                supported case, if ``operation`` is omitted, the default strategy is
                ``"Q5"``.

                Current implementation limits
                    - vectorized ``MutSpec`` fields are not implemented
                    - batch ``primerdesign(template=[...], target=[...], mut_pattern=...)``
                      is not implemented
                    - ``template != target`` mutagenesis (for example junction /
                      overlap-PCR mutagenesis) is not implemented
                    - ``operation="overlappcr"`` is reserved but not yet implemented

            relative : {"target"} or QUEEN or str, optional
                Reference subsequence for interpreting ``loc``/``find``. If omitted or
                "target", coordinates are relative to ``target``.
                If a `QUEEN` object or DNA string is provided, it must map uniquely within
                ``target``. Coordinates are interpreted in the 5′→3′ direction of the
                ``relative`` sequence; if a match is on the reverse strand, mapping to
                ``target`` is handled automatically (including reverse-complementing
                ``find`` as needed). Non-unique matches should raise an error.

            loc : (start, end) or list of (start, end), optional
                Interval(s) to replace, relative to ``relative`` if provided, otherwise
                relative to ``target``. Use ``(i, i)`` to denote an insertion at index ``i``.

            find : str or list of str, optional
                Alternative to ``loc``. A unique subsequence within the ``relative`` sequence
                (or within ``target`` when ``relative`` is omitted) that specifies the edit site.

            to : str or list of str, optional
                Replacement/insert sequence(s) (IUPAC may be accepted). Use ``""`` or
                ``None`` for deletion.

        Vectorization rules (MutSpec)
            - ``relative``, ``loc``, ``find``, and ``to`` accept scalars or lists.
            - If multiple fields are lists, their (non-None) lengths must match; scalars
              broadcast.
            - For each site k, specify either ``loc[k]`` OR ``find[k]`` (not both).
            - Deletion: ``to[k]`` in {``""``, ``None``}. Insertion: ``loc[k] == (i, i)``.
            - Note: vectorized mutagenesis is documented for future extension but is
              not implemented in the current runtime.

        Batch mode (MutSpec)
            - Batch mutagenesis is not implemented in the current runtime.

    target_tm : float or sequence of float or None, optional
        Desired melting temperature (Tm) for primers in degrees Celsius. Default is ``60.0``.
        If ``None``, primer candidates are not ranked by Tm proximity and are returned in
        their generated order after specificity/requirement filtering.
    nonspecific_limit : int or sequence of int, optional
        Specificity filter threshold. Candidate primers that bind to any region of the template
        with mismatches <= this value (outside the intended binding) are excluded to reduce
        nonspecific amplification. Default is ``3``.
    auto_adjust : bool or int or sequence, optional
        If enabled (default is ``1`` which behaves like ``True``), and partner-derived tails
        are used, the function may adjust the junction to preserve reading frame when joining
        coding sequences (CDS–CDS contexts). Depending on the implementation, this may include
        inserting a short gap sequence (typically 0–2 bases) between partner-derived tail and
        primer binding region so the join length becomes a multiple of 3.

        Note
            If ``gap`` is not supplied, some implementations may generate the gap sequence
            automatically (and in some cases pseudo-randomly). Always verify junction
            sequences with ``printsequence()``.
    homology_length : int or sequence of int, optional
        Active when ``adapter_mode`` is one of ``{"gibson","infusion","overlappcr","RE"}``.
        Controls the minimum homology length used when constructing homology tails with
        partner context. Default is ``30`` (function default). (Some older docs referenced
        20; current default is 30.)
    tm_func : str or callable or sequence, optional
        Function/selector used to calculate primer Tm.

        - As a string selector, supported values include:
          ``"SantaLucia"`` / ``"sa"`` and ``"Breslauer"`` / ``"br"``.
        - As a callable, it should be compatible with Biopython’s
          ``Bio.SeqUtils.MeltingTemp.Tm_NN``-like interface (accepting ``seq`` and optional
          parameters; excluding ``seq`` and ``c_seq`` in some wrapper implementations).

        If ``None``, a SantaLucia-like nearest-neighbor model is used by default.
    primer_length : tuple of int or sequence of tuple, optional
        Primer length bounds as ``(min_len, max_len)``. Default is ``(16, 25)``.
    design_num : int or sequence of int, optional
        Number of primer pairs to return per template/target. Default is ``1``.
    gap : object or sequence, optional
        Optional gap specification used for frame adjustment logic in multi-fragment or
        partner-aware workflows. This is mainly for internal/batch workflows.
        If provided, it is typically a pair ``(gap_fw, gap_rv)`` where each element is a
        short DNA string (often 0–2 bases) or ``None``. Default is ``None``.
    batch_process : bool, optional
        Internal flag used in batch homology workflows. If ``True``, the function returns
        an intermediate amplicon region (with adapters applied if provided) and gap
        information instead of primer pairs. Default is ``False``.
    product : str, optional
        Reserved metadata field accepted for API symmetry with other QUEEN design/
        construction helpers. It does not change primer selection.
    process_name : str, optional
        Reserved metadata field accepted for API symmetry with other QUEEN design/
        construction helpers. It does not change primer selection.
    process_description : str, optional
        Reserved metadata field accepted for API symmetry with other QUEEN design/
        construction helpers. It does not change primer selection.
    pn : str, optional
        Alias for ``process_name``.
    pd : str, optional
        Alias for ``process_description``.

    Returns
    -------
    list of dict
        In single-template mode, returns a list of primer-pair records sorted by how close
        the primer Tm values are to ``target_tm`` (closest first). When ``target_tm=None``,
        records are returned in generated order after filtering. Each record contains
        at least:

        - ``"fw"`` : `QUEEN` (ssDNA) forward primer
        - ``"rv"`` : `QUEEN` (ssDNA) reverse primer

        Implementations may also include fields such as ``"fw_tm"`` and ``"rv_tm"``.
    list of list of dict
        In batch mode (``template`` is a list/tuple), returns a list where each element is
        the corresponding primer-pair list for that template/target pair (same order).
    tuple
        If ``batch_process=True``, returns ``(amplicon_region, gap_fw, gap_rv)``, where
        ``amplicon_region`` is a `QUEEN` representing the design region (possibly with
        adapters applied) and ``gap_fw`` / ``gap_rv`` are optional gap sequences inferred
        by the partner-aware logic.

    Raises
    ------
    TypeError
        If ``template``/``target`` are not `QUEEN` objects (or lists of `QUEEN` objects),
        or if list-valued arguments are of invalid types.
    ValueError
        If the target sequence (or its reverse complement) is not found within the template
        sequence.
    ValueError
        If mutagenesis specification is malformed (e.g., missing both ``loc`` and ``find``,
        non-unique mapping for ``relative``/``find``, inconsistent vectorized lengths), or
        if adapter/partner values are inconsistent with ``adapter_mode`` (e.g., partner is
        required but missing or of incompatible type).

    Notes
    -----
    It is assumed that template and target are provided as `QUEEN` objects representing
    the sequence context for primer design. The requirement that the target sequence is
    within the template sequence ensures specificity of primer binding. The function will
    not proceed when the target is not a subset (or reverse-complement subset) of the template.

    Examples
    --------
    Design primers for a target subregion within a template::

        template = QUEEN(seq="ATGC" * 200)
        target   = template[100:200]
        primers  = primerdesign(template, target, target_tm=65.0, design_num=2)
        fw = primers[0]["fw"]
        rv = primers[0]["rv"]
    
    Design primer pairs for a homology_based_assembly of two PCR amplicons::

        template1  = QUEEN(record="123456", dbtype="addgene")
        insert     = template["Gene1"]
        template2  = QUEEN(record="654231", dbtype="addgene")
        backbone   = template["!Gene2"]
        primer_pairs = primerdesign([template1, template2], [insert, backbone], target_tm=65.0)
        primer_pair1 = primer_pairs[0][0]
        primer_pair2 = primer_pairs[1][0] 
        
        #pcr using the desgined primer pairs
        fw1 = primer_pair1["fw"] 
        rv1 = primer_pair1["rv"]
        insert_amplicon = pcr(template1, fw1, rv1) 
        
        fw2 = primer_pair2["fw"] 
        rv2 = primer_pair2["rv"]
        backbone_amplicon = pcr(template2, fw2, rv2) 

        #homology_based_assembly of insert and backbone amplicons. 
        construct = homology_based_assembly(insert_amplicon, backbone_amplicon)
    
    Design primers for introducing a site-specific mutation::
        
        template     = QUEEN(record="111111", dbtype="addgene")
        goi_region   = template["GeneX"]
        mut_pattern  = {"relative": goi_region, "loc": (102,105), "to": "GCC", "operation": "gibson"}
        primers = primerdesign(template=template, target=template, mut_pattern=mut_pattern)
        
        #pcr using the designed primers. 
        fw = primers[0]["fw"]
        rv = primers[0]["rv"]
        amplicon = pcr(template, fw, rv) 

        #homology_based_assembly for joining a linear PCR amplicon holding a site-specific mutation.   
        construct = homoogy_based_assembly(amplicon) 
    
    """ 
    
    process_name = pn if process_name is None else process_name
    process_description = pd if process_description is None else process_description
    _ = product
    _ = process_name
    _ = process_description

    def search_qexps(dna):
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
        qexps = [] 
        for line in quine(dna, _return_script=True):
            for key in pattern_dict:
                amatch = re.search(pattern_dict[key], line) 
                if amatch is None:
                    pass 
                else:
                    qexps.append((key, amatch.group(1).split())) 
        return qexps 
    
    def append_adapter(amplicon_region, filtered_primer_pairs, adapter, partner, mode, homology_length, strand, name, gapseq, auto_adjust): 
        iupac_map = {
            "A": {"A"},
            "C": {"C"},
            "G": {"G"},
            "T": {"T"},
            "R": {"A", "G"},
            "Y": {"C", "T"},
            "K": {"G", "T"},
            "M": {"A", "C"},
            "S": {"C", "G"},
            "W": {"A", "T"},
            "B": {"C", "G", "T"},
            "D": {"A", "G", "T"},
            "H": {"A", "C", "T"},
            "V": {"A", "C", "G"},
            "N": {"A", "C", "G", "T"},
        }

        def _endseq_matches(site_endseq, observed_endseq):
            site_endseq = str(site_endseq).upper()
            observed_endseq = str(observed_endseq).upper()
            if site_endseq == observed_endseq:
                return True
            if len(site_endseq) != len(observed_endseq):
                return False
            for schar, ochar in zip(site_endseq, observed_endseq):
                if ochar not in iupac_map.get(schar, {schar}):
                    return False
            return True

        def _extract_digestion_cutsites(dna):
            qexps = search_qexps(dna)
            if len(qexps) == 0 or qexps[-1][0] != "digestion":
                raise ValueError("When 'adapter_mode' is 'RE', partner value must be a digested QUEEN object.")
            cutsites = []
            for arg in qexps[-1][1][1:]:
                if "selection" in arg:
                    break
                token = arg.rstrip(",")
                if len(token) >= 2 and token[0] == token[-1] and token[0] in ("'", '"'):
                    token = token[1:-1]
                cutsites.append(token)
            return cutsites

        def _resolve_re_partner_context(dna, strand):
            cutsites = _extract_digestion_cutsites(dna)
            if strand == "fw":
                endseq = dna._right_end
                observed_pair = (dna._right_end_top, dna._right_end_bottom)
            else:
                endseq = dna._left_end
                observed_pair = (dna._left_end_bottom, dna._left_end_top)

            exact_matches = []
            blunt_matches = []
            for cutsite in cutsites:
                site = cs.lib[cutsite]
                if _endseq_matches(site.endseq, endseq) == False:
                    continue
                if (site.top, site.bottom) == observed_pair:
                    exact_matches.append(cutsite)
                elif endseq == "":
                    blunt_matches.append(cutsite)

            matches = exact_matches
            if len(matches) == 0 and len(blunt_matches) == 1:
                matches = blunt_matches

            if len(matches) != 1:
                raise ValueError("Failed to infer the restriction site corresponding to the digested partner end in RE primer design.")

            cutsite = matches[0]
            if strand == "fw":
                partner_seq = QUEEN(seq="ATGC" + cs.lib[cutsite].seq, quinable=False)
                remseq = cutdna(partner_seq, *partner_seq.searchsequence(cs.lib[cutsite], quinable=False), quinable=0)[-1]
            else:
                partner_seq = QUEEN(seq="ATGC" + cs.lib[cutsite].rcseq, quinable=False)
                remseq = cutdna(partner_seq, *partner_seq.searchsequence(cs.lib[cutsite], quinable=False), quinable=0)[0]
            return cutsite, partner_seq.seq, remseq

        if (type(adapter) == str and adapter == "") or (adapter is None):
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq=filtered_primer_pairs[i][strand][0])
            
        else:
            if mode == "BP" and type(adapter) == str and adapter in ("attB1", "attB2"):
                attb_map = {
                    "attB1": "GGGGACAAGTTTGTACAAAAAAGCAGGCT",
                    "attB2": "GGGGACCACTTTGTACAAGAAAGCTGGGT",
                }
                adapter = QUEEN(seq=attb_map[adapter])
                adapter = flipdna(adapter, quinable=0) if strand == "rv" else adapter
                for i in range(len(filtered_primer_pairs)):
                    filtered_primer_pairs[i][strand][0] = joindna(
                        adapter,
                        QUEEN(seq=filtered_primer_pairs[i][strand][0], quinable=0),
                        homology_length=0,
                        quinable=0,
                    )

            elif type(adapter) == QUEEN or (type(adapter) in (str, Qseq) and set(adapter.upper()) <= set("ATGCRYKMSWBDHVN")):
                if type(adapter) == QUEEN:
                    pass 
                else:
                    adapter = QUEEN(seq=adapter)
                adapter = flipdna(adapter, quinable=0) if strand == "rv" else adapter
                for i in range(len(filtered_primer_pairs)):
                    filtered_primer_pairs[i][strand][0] = joindna(adapter, QUEEN(seq=filtered_primer_pairs[i][strand][0], quinable=0), homology_length=0, quinable=0) 
            
            elif (type(adapter) == str and adapter in cs.lib.keys()) or "Cutsite" in type(adapter).__name__:
                if type(adapter) == str: 
                    adapter = cs.lib[adapter]
                for i in range(len(filtered_primer_pairs)):
                    adseq = adapter.rcseq if strand == "rv" else adapter.seq 
                    filtered_primer_pairs[i][strand][0] = QUEEN(seq="ATGC" + adseq + filtered_primer_pairs[i][strand][0]) 
            else:
                raise ValueError("When 'adapter_mode' is 'standard' or 'BP', adapter value must be a QUEEN, DNA string, Cutsite object, or attB1/attB2 token.")
        
        if mode in ("gibson", "infusion", "overlappcr", "RE"):
            if partner is None:
                pass 
            else:
                if type(partner) == QUEEN and partner._ssdna == False:
                    pass
                elif type(partner) in (str, Qseq) and set(str(partner).upper()) <= set("ATGCRYKMSWBDHVN"):
                    partner = QUEEN(seq=str(partner))
                else:
                    raise ValueError("When 'adapter_mode' is 'gibson', 'infusion', or 'overlappcr', partner value must be a dsDNA QUEEN object or str object.")

                partner_features  = [feat for feat in partner.dnafeatures if feat.feature_type not in ("source", "primer", "primer_bind")]
                amplicon_features = [feat for feat in amplicon_region.dnafeatures if feat.feature_type not in ("source", "primer", "primer_bind")]
                amplicon_features.sort(key=lambda x: x.start) 
                if len(partner_features) == 0 or len(amplicon_features) == 0:
                    req = False
                else:
                    if strand == "fw":
                        feat1 = partner_features[-1] 
                        feat2 = amplicon_features[0] 

                    if strand == "rv":
                        feat1 = amplicon_features[-1] 
                        feat2 = partner_features[0] 

                    # Promoter/CDS strand combinations across a partner junction are
                    # not sufficient evidence that the requested overlap direction is
                    # wrong. In modular plasmid assemblies a valid boundary can be
                    # `CDS(-) -> promoter(+)` or `promoter(+) -> CDS(-)` depending on
                    # circular origin choice and which exact donor block is being
                    # preserved. Keep the stricter CDS/CDS check below, but do not
                    # reject promoter/CDS boundaries here.
                                
                    if feat1.feature_type == "CDS" and feat2.feature_type == "CDS":
                        if feat1.strand == feat2.strand:
                            if ("broken_feature" in feat1.qualifiers or "broken_feature" in feat2.qualifiers) and (len(feat1.sequence)%3 != 0 or len(feat2.sequence)%3 != 0):
                                req = False
                            else:
                                req = True
                        else:
                            # Opposite-strand CDS/CDS junctions can be valid in modular
                            # plasmid assemblies, but they do not provide a safe basis for
                            # codon-frame auto adjustment. Keep primer design going and
                            # simply skip the frame-aware gap logic for this boundary.
                            req = False
                    else:
                        req = False 

                gapflag = 0 
                for i in range(len(filtered_primer_pairs)):
                    remseq = ""
                    if strand == "fw": 
                        if mode == "gibson":
                            if partner._right_end_bottom == 1 and partner._right_end_top == -1: 
                                mod_partner = partner[:len(partner.seq) - len(partner._right_end)]
                            else:
                                mod_partner = partner
                            partner_seq = mod_partner.seq[-1*homology_length:]
                        
                        elif mode == "infusion":
                            if partner._right_end_bottom == -1 and partner._right_end_top == 1:
                                mod_partner = partner[:len(partner.seq) - len(partner._right_end)] 
                            else:
                                mod_partner = partner
                            partner_seq = mod_partner.seq[-1*homology_length:]
                        
                        elif mode == "RE":
                            cutsite, partner_seq, remseq = _resolve_re_partner_context(partner, strand)

                    if strand == "rv":
                        if mode == "gibson":
                            if partner._left_end_bottom == -1 and partner._left_end_top == 1: 
                                mod_partner = partner[len(partner._left_end):] 
                            else:
                                mod_partner = partner
                            partner_seq = mod_partner.rcseq[-1*homology_length:]
                        
                        elif mode == "infusion":
                            if partner._left_end_bottom == 1 and partner._left_end_top == -1:
                                mod_partner = partner[len(partner._left_end):] 
                                partner_seq = mod_partner.rcseq[-1*homology_length:]
                            else:
                                mod_partner = partner
                            partner_seq = mod_partner.rcseq[-1*homology_length:]
     
                        elif mode == "RE":
                            cutsite, partner_seq, remseq = _resolve_re_partner_context(partner, strand)
                    
                    if mode == "RE":
                        if strand == "fw":
                            mod_partner = partner[:-1*len(cs.lib[cutsite].endseq)] 
                        else:
                            mod_partner = partner[len(cs.lib[cutsite].endseq):] 
                   
                    if req == True and auto_adjust == True:
                        try:
                            if strand == "fw":
                                feat1 = partner_features[-1] 
                                feat2 = amplicon_features[0] 
                                fragment1 = mod_partner[feat1.start:].seq  
                                fragment2 = amplicon_region[:feat2.end].seq 
                                rem = (len(fragment1) + len(fragment2) + len(remseq)) % 3
                                if i == 0 and gapseq is None:
                                    if rem > 0:
                                        gapflag = 1
                                        gapseq = _deterministic_gap_seq(fragment1, fragment2, remseq, 3-rem)
                                    else:
                                        gapseq = ""
                                filtered_primer_pairs[i][strand][0] = QUEEN(seq="", product=name) + partner_seq + gapseq + filtered_primer_pairs[i][strand][0]
                            else:
                                feat1 = amplicon_features[-1] 
                                feat2 = partner_features[0] 
                                fragment1 = amplicon_region[feat1.start:].seq  
                                fragment2 = mod_partner[:feat2.end].seq 
                                rem = (len(fragment1) + len(fragment2) + len(remseq)) % 3 
                                if i == 0 and gapseq is None:
                                    if rem > 0:
                                        gapflag = 1
                                        gapseq = _deterministic_gap_seq(fragment1, fragment2, remseq, 3-rem)
                                    else:
                                        gapseq = ""
                                filtered_primer_pairs[i][strand][0] = QUEEN(seq="", product=name) + partner_seq + gapseq + filtered_primer_pairs[i][strand][0]
                        except ValueError:
                            filtered_primer_pairs[i][strand][0] = QUEEN(seq="", product=name) + partner_seq + filtered_primer_pairs[i][strand][0]
                    else:
                        filtered_primer_pairs[i][strand][0] = QUEEN(seq="", product=name) + partner_seq + filtered_primer_pairs[i][strand][0]
                
                if gapflag == 1: 
                    if strand == "fw":
                        return "fw", gapseq
                    else:
                        return "rv", gapseq 
        else:
            for i in range(len(filtered_primer_pairs)):
                filtered_primer_pairs[i][strand][0] = QUEEN(seq="", product=name) + filtered_primer_pairs[i][strand][0]  
        
        for i in range(len(filtered_primer_pairs)):
            filtered_primer_pairs[i][strand][0]._ssdna = True
        
        return filtered_primer_pairs
     
    if requirement is None:
        def requirement(x):
            req1 =  x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T") 
            req4 = "AAAA" not in x["fw"] and "TTTT" not in x["fw"] and "GGGG" not in x["fw"] and "CCCC" not in x["fw"]
            req5 = "AAAA" not in x["rv"] and "TTTT" not in x["rv"] and "GGGG" not in x["rv"] and "CCCC" not in x["rv"]
            return req1 and req4 and req5
        
    if type(template) != QUEEN: 
        if type(template) in (tuple, list) and list(set(map(type, template)))[0] == QUEEN:
            pass 
        else: 
            raise TypeError("`template` object must be instance of QUEEN class or a list of QUEEN objects")  
    else:
        pass 

    if type(target) != QUEEN:
        if type(target) in (tuple, list) and list(set(map(type, target)))[0] == QUEEN:
            if len(template) == len(target): 
                new_target = [] 
                for te, ta in zip(template, target):
                    if -1 in (ta._left_end_top, ta._left_end_bottom, ta._right_end_top, ta._right_end_bottom):
                        ta = modifyends(ta, quinable=False)
                    else:    
                        pass 
                    if ta.seq in te.seq: 
                        pass
                    elif ta.rcseq in te.seq:
                        ta = flipdna(ta, quinable=0)
                    else:
                        raise ValueError("target sequence to be amplified is not included in template sequence.") 
                    new_target.append(ta)
                target = new_target
            else:
                raise ValueError("The length of target should be same with the template.")
        else: 
            raise TypeError("`target` object must be instance of QUEEN class or a list of QUEEN objects.") 
    else:
        if -1 in (target._left_end_top, target._left_end_bottom, target._right_end_top, target._right_end_bottom):
            target = modifyends(target, quinable=False) 
        
        if target.seq in template.seq: 
            pass
        elif target.rcseq in template.seq:
            target = flipdna(target, quinable=0)
        else:
            raise ValueError("target sequence to be amplified is not included in template sequence.") 

    if type(template) in (tuple, list):
        if mut_pattern is not None:
            raise NotImplementedError(
                "Batch primerdesign with mut_pattern is not implemented yet. "
                "For now, MutSpec is supported only in single-template, single-target calls."
            )
        fw_primers     = [fw_primer] * len(template) if type(fw_primer) not in (tuple, list) else fw_primer
        rv_primers     = [rv_primer] * len(template) if type(rv_primer) not in (tuple, list) else rv_primer
        fw_margins     = [fw_margin] * len(template) if type(fw_margin) not in (tuple, list) else fw_margin
        rv_margins     = [rv_margin] * len(template) if type(rv_margin) not in (tuple, list) else rv_margin
        target_tms     = [target_tm] * len(template) if type(target_tm) not in (tuple, list) else target_tm
        tm_funcs       = [tm_func] * len(template) if type(tm_func) not in (tuple, list) else tm_func
        primer_lengths = [primer_length] * len(template) if type(primer_length[0]) not in (tuple, list) else primer_length 
        design_nums    = [design_num] * len(template) if type(design_num) not in (tuple, list) else design_num
        gaps           = [gap] * len(template) if type(design_num) not in (tuple, list) else design_num
        fw_adapters    = [fw_adapter] * len(template) if type(fw_adapter) not in (tuple, list) else fw_adapter
        rv_adapters    = [rv_adapter] * len(template) if type(rv_adapter) not in (tuple, list) else rv_adapter
        adapter_modes  = [adapter_mode] * len(template) if type(adapter_mode) not in (tuple, list) else adapter_mode
        
        if adapter_mode in ("gibson", "infusion", "RE", "overlappcr"):
            if fw_partner is None and rv_partner is None:
                fw_partners = []  
                rv_partners = [] 
                if adapter_mode in ("gibson", "infusion", "RE"): 
                    for i in range(0, len(target)):
                        if i < len(target) - 1:
                            fw_partners.append(target[i-1])
                            rv_partners.append(target[i+1]) 
                        else:
                            fw_partners.append(target[i-1])
                            rv_partners.append(target[0])

                if adapter_mode == "overlappcr":
                    for i in range(0, len(target)):
                        if i == 0:
                            fw_partners.append(None)
                            rv_partners.append(target[i+1]) 
                        elif i == len(target) - 1:
                            fw_partners.append(target[i-1])
                            rv_partners.append(None) 
                        else:
                            fw_partners.append(target[i-1])
                            rv_partners.append(target[i+1]) 
            else:
                fw_partners = fw_partner
                rv_partners = rv_partner
        else:
            fw_partners = [fw_partner] * len(template) if type(fw_partner) not in (tuple, list) else fw_partner
            rv_partners = [rv_partner] * len(template) if type(rv_partner) not in (tuple, list) else rv_partner 
 
        homology_lengths   = [homology_length] * len(template) if type(homology_length) not in (tuple, list) else homology_length
        nonspecific_limits = [nonspecific_limit] * len(template) if type(nonspecific_limit) not in (tuple, list) else nonspecific_limit
        auto_adjusts       = [auto_adjust] * len(template) if type(auto_adjust) not in (tuple, list) else auto_adjust
        requirements       = [requirement] * len(template) if type(requirement) not in (tuple, list) else requirement
        fw_names           = [fw_name] * len(template) if type(fw_name) not in (tuple, list) else fw_name
        rv_names           = [rv_name] * len(template) if type(rv_name) not in (tuple, list) else rv_name 
        mut_patterns       = [mut_pattern] * len(template) if type(mut_pattern) not in (tuple, list) else mut_pattern 
        arguments = list(zip(*[template, target, fw_primers, rv_primers, fw_margins, rv_margins, adapter_modes, fw_adapters, rv_adapters, fw_partners, rv_partners, requirements, fw_names, rv_names, mut_patterns, target_tms, nonspecific_limits, auto_adjusts, homology_lengths, tm_funcs, primer_lengths, design_nums, gaps]))
        arguments = list(map(list, arguments))

        primer_pair_set = []
        if adapter_mode in ("gibson", "infusion", "RE", "overlappcr"):
            new_targets = []  
            for i, argument in enumerate(arguments):
                argument.append(True) 
                new_target, gap_fw, gap_rv = primerdesign(*argument)
                new_targets.append(new_target)
                argument[-2] = [gap_fw, gap_rv] 
                argument[-1] = False
            
            # For Gibson/Infusion batch primer design, each junction must share
            # one common overlap on both adjacent fragments. The previous
            # partner-based implementation derived the left and right overlaps
            # independently from opposite partner ends, which can yield
            # incompatible junctions when fragment boundaries are not identical.
            # When no frame-adjustment gap is required, derive explicit shared
            # overlaps from the ordered batch targets instead.
            #
            # Important: preserve any user-specified payload adapters already
            # assigned to fw_adapter/rv_adapter. For example, a Gibson insert may
            # need both a homology overlap and a coding payload such as a tag or
            # linker. In that case the final amplicon must encode:
            #   fw end = shared overlap + user fw payload
            #   rv end = user rv payload + shared overlap
            # The earlier second-pass patch for Q003 replaced the user adapter
            # with the shared overlap and silently truncated such payloads.
            if adapter_mode in ("gibson", "infusion") and False not in [arguments[i][-2] == [None, None] for i in range(len(arguments))]:
                def _adapter_seq(adapter):
                    if adapter is None:
                        return ""
                    if type(adapter) == QUEEN:
                        return str(adapter.seq)
                    return str(adapter)

                for i, target in enumerate(new_targets):
                    hlen = int(homology_lengths[i] / 2)
                    # In ordered batch Gibson/Infusion design, the forward-primer
                    # 5' overlap for fragment i must match the end of the previous
                    # fragment, while the reverse-primer 5' overlap must match the
                    # start of the next fragment. Using the current fragment prefix
                    # on the forward side duplicates the fragment's own 5' sequence
                    # and diverges from the explicit fw_partner/rv_partner route.
                    if i > 0:
                        shared_fw = str(new_targets[i - 1].seq[-hlen:])
                    else:
                        shared_fw = str(new_targets[-1].seq[-hlen:])
                    arguments[i][-1] = False
                    arguments[i][7] = shared_fw + _adapter_seq(arguments[i][7])
                    if i < len(new_targets) - 1:
                        shared_rv = str(new_targets[i + 1].seq[:hlen])
                    else:
                        shared_rv = str(new_targets[0].seq[:hlen])
                    arguments[i][8] = _adapter_seq(arguments[i][8]) + shared_rv
                    arguments[i][9] = None
                    arguments[i][10] = None
                    primer_pair = primerdesign(*arguments[i])
                    primer_pair_set.append(primer_pair)
            else:
                for i, argument in enumerate(arguments):
                    if type(arguments[i][-2][0]) == str:
                        arguments[i-1][-2][1] = arguments[i][-2][0].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                    else:
                        pass 
                
                fw_partners = []  
                rv_partners = [] 
                for i, target in enumerate(new_targets):
                    if adapter_mode in ("gibson", "infusion", "RE"): 
                        if i < len(new_targets) - 1:
                            fw_partners.append(new_targets[i-1])
                            rv_partners.append(new_targets[i+1]) 
                        else:
                            fw_partners.append(new_targets[i-1])
                            rv_partners.append(new_targets[0])

                    elif adapter_mode == "overlappcr":
                        if i == 0:
                            fw_partners.append(None)
                            rv_partners.append(new_targets[i+1]) 
                        elif i == len(new_targets) - 1:
                            fw_partners.append(new_targets[i-1])
                            rv_partners.append(None) 
                        else:
                            fw_partners.append(new_targets[i-1])
                            rv_partners.append(new_targets[i+1]) 

                for i, (fw_partner, rv_partner) in enumerate(zip(fw_partners, rv_partners)):
                    arguments[i][-1] = False
                    arguments[i][9]  = fw_partner
                    arguments[i][10] = rv_partner 
                    primer_pair = primerdesign(*arguments[i]) 
                    primer_pair_set.append(primer_pair) 
                 
        else:
            for i, argument in enumerate(arguments):
                argument.append(False) 
                primer_pair = primerdesign(*argument) 
                primer_pair_set.append(primer_pair) 

        return primer_pair_set 
        
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
    elif tm_func == "SantaLucia" or "sa":
        tm_func = Tm_NN(nn_table=mt.DNA_NN3)
    elif tm_func == "Breslauer" or "br":
        tm_func = Tm_NN(nn_table=mt.DNA_NN1)

    fw_margin_auto = isinstance(fw_margin, str) and fw_margin.lower() == "auto"
    rv_margin_auto = isinstance(rv_margin, str) and rv_margin.lower() == "auto"
    current_fw_margin = 0 if fw_margin_auto else fw_margin
    current_rv_margin = 0 if rv_margin_auto else rv_margin
    auto_margin_step = 5
    auto_margin_max = 100

    target_start = template.seq.find(target.seq)
    if target_start < 0 and template.topology == "circular":
        doubled_template = str(template.seq) + str(template.seq)
        target_start = doubled_template.find(str(target.seq))
        if target_start >= len(template.seq):
            target_start = -1

    start = target_start - current_fw_margin
    if start < 0:
        if template.topology == "circular":
            start = len(template.seq) + start
        else:
            start = 0
    
    end = target_start + len(target.seq) + current_rv_margin
    if end > len(template.seq):
        if template.topology == "circular": 
            end = end - len(template.seq) * (end // len(template.seq))  
        else:
            pass
    
    amplicon_region = template[start:end]
    if batch_process == True:
        if type(fw_adapter) == QUEEN or (type(fw_adapter) == str and set(fw_adapter.upper()) <= set("ATGCRYKMSWBDHVN")):
            if type(fw_adapter) == QUEEN:
                pass 
            else:
                fw_adapter = QUEEN(seq=fw_adapter, quinable=0)

        elif (type(fw_adapter) == str and fw_adapter in cs.lib.keys()) or "Cutsite" in type(fw_adapter).__name__:
            if type(fw_adapter) == str: 
                fw_adapter = cs.lib[fw_adapter].seq
            else:
                fw_adatepr = fw_adapter.seq
        
        if type(rv_adapter) == QUEEN or (type(rv_adapter) == str and set(rv_adapter.upper()) <= set("ATGCRYKMSWBDHVN")):
            if type(rv_adapter) == QUEEN:
                pass 
            else:
                rv_adapter = QUEEN(seq=rv_adapter, quinable=0)

        elif (type(rv_adapter) == str and rv_adapter in cs.lib.keys()) or "Cutsite" in type(rv_adapter).__name__:
            if type(rv_adapter) == str: 
                rv_adapter = cs.lib[rv_adapter].seq
            else:
                rv_adapter = rv_adapter.seq
        
        if fw_adapter is None and rv_adapter is None:
            amplicon_region = amplicon_region 
        elif fw_adapter is None:
            amplicon_region = amplicon_region + rv_adapter 
        elif rv_adapter is None:
            amplicon_region = fw_adapter + amplicon_region
        else:
            amplicon_region = fw_adapter + amplicon_region + rv_adapter  
        
        dammy = [{"fw":["ATGC", 4], "rv":["ATGC", 4], "fw_tm":50, "rv_tm":50}]
        gapinfo_fw = append_adapter(amplicon_region, dammy, fw_adapter, fw_partner, adapter_mode, int(homology_length/2), "fw", fw_name, None, auto_adjust)
        gapinfo_rv = append_adapter(amplicon_region, dammy, rv_adapter, rv_partner, adapter_mode, int(homology_length/2), "rv", rv_name, None, auto_adjust)
        if gapinfo_fw[0] == "fw":
            gap_fw = gapinfo_fw[1] 
        else:
            gap_fw = None
        
        if gapinfo_rv[0] == "rv":
            gap_rv = gapinfo_rv[1]
        else:
            gap_rv = None
        
        return amplicon_region, gap_fw, gap_rv 
    
    else: 
        if mut_pattern is not None:
            flip = 0 
            if "relative" not in mut_pattern or mut_pattern["relative"] is None:
                rs = amplicon_region.seq.find(target.seq)
            else:
                rs = amplicon_region.seq.find(mut_pattern["relative"].seq)
                if rs == -1:
                    amplicon_region = flipdna(amplicon_region, quinable=0)
                    template = flipdna(template, quinable=0)
                    rs = amplicon_region.seq.find(mut_pattern["relative"].seq)
                    flip = 1

            aligner = PairwiseAligner(mode="global")
            aligner.mismatch_score = 1.0
            aligner.target_end_gap_score = -0.1
            aligner.query_end_gap_score  = -0.1

            if type(mut_pattern.get("to")) in (tuple, list):
                raise NotImplementedError(
                    "Vectorized MutSpec is not implemented yet. Pass one scalar edit per primerdesign() call."
                )

            if template.seq != amplicon_region.seq:
                raise NotImplementedError(
                    "MutSpec with template != target is not implemented yet. "
                    "Junction/overlap-PCR mutagenesis should be modeled explicitly rather than through mut_pattern."
                )

            def _slice_context(seq, start, length, circular):
                if length <= 0:
                    return ""
                if circular:
                    seq_len = len(seq)
                    if seq_len == 0:
                        return ""
                    start = start % seq_len
                    end = start + length
                    if end <= seq_len:
                        return seq[start:end]
                    return seq[start:] + seq[:end-seq_len]
                if start < 0 or start + length > len(seq):
                    raise ValueError("Deletion primer design exceeds the template boundary. Increase margins or use a circular template.")
                return seq[start:start+length]
          
            if template.seq == amplicon_region.seq and type(mut_pattern["to"]) not in (tuple, list):
                if "operation" not in mut_pattern:
                    operation = "Q5"
                else:
                    operation = mut_pattern["operation"]  

                if operation == "overlappcr":
                    raise NotImplementedError(
                        "MutSpec operation='overlappcr' is documented but not implemented yet."
                    )
                
                if "loc" in mut_pattern:
                    loc = mut_pattern["loc"]
                    loc = (rs + loc[0], rs + loc[1])
                elif "find" in mut_pattern:
                    s   = amplicon_region.seq.find(mut_pattern["find"]) 
                    e   = s + len(mut_pattern["find"]) 
                    loc = (rs + s, rs + e) 
                else:
                    raise ValueError("Mutation location is not specified") 
                
                origins   = amplicon_region.seq[loc[0]:loc[1]] 
                mutations = mut_pattern["to"]
                deletion_mode = mutations in ("", None)

                if deletion_mode:
                    target_seq = amplicon_region.seq
                    mutate_seq = amplicon_region.seq[:loc[0]] + amplicon_region.seq[loc[1]:]
                    mut_pos = loc[0]
                else:
                    aln  = aligner.align(origins, mutations)[0] 
                    for i, (o, m) in enumerate(zip(aln[0], aln[1])):
                        if o != m:
                            mut_pos = loc[0]+i
                            break 
                
                    target_seq = QUEEN(seq=amplicon_region.seq[:loc[0]] + aln[0] + amplicon_region.seq[loc[1]:], topology="circular", quinable=0).seq
                    mutate_seq = QUEEN(seq=amplicon_region.seq[:loc[0]] + aln[1] + amplicon_region.seq[loc[1]:], topology="circular", quinable=0).seq
                
                if operation in ("Q5", "QuickChange"):
                    if deletion_mode:
                        raise ValueError("Deletion mutagenesis with 'Q5'/'QuickChange' is not supported. Use operation='gibson' or 'infusion'.")
                    fw_tm_set = [] 
                    rv_tm_set = [] 
                    for plen in range(primer_length[0], primer_length[1]):
                        flen = int(plen/2) + plen%2
                        rlen = int(plen/2)
                        fw_candidate = mutate_seq[mut_pos-flen:mut_pos+rlen] 
                        rv_candidate = target_seq[mut_pos-flen:mut_pos+rlen] 
                        fw_candidate = fw_candidate.replace("-","") 
                        rv_candidate = rv_candidate.replace("-","")
                        aln = aligner.align(fw_candidate, rv_candidate)[0]
                        try:
                            fwseq = ""
                            rvseq = "" 
                            for f,r in zip(aln[0], aln[1]):
                                if f == "-" or r == "-":
                                    pass 
                                else:
                                    fwseq += f
                                    rvseq += r
                            tm = tm_func(seq=fwseq, c_seq=rvseq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")))
                        except:
                            fwseq = ""
                            rvseq = "" 
                            for f,r in zip(aln[0], aln[1]):
                                if f != "-":
                                    fwseq += f
                                    rvseq += r
                            tm = tm_func(seq=fwseq, c_seq=rvseq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")))
                        fw_tm_set.append([[fw_candidate, 0], tm]) 
                        
                        if operation == "QuickChange": 
                            rv_candidate = rv_candidate.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                        else:
                            rv_candidate = target_seq[mut_pos-flen-plen:mut_pos-flen].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                            tm = tm_func(seq=rv_candidate) 
                        rv_tm_set.append([[rv_candidate, 0], tm])
                    
                    if flip == 1:
                        fw_tm_set, rv_tm_set       = rv_tm_set, fw_tm_set
                        fw_candidate, rv_candidate = rv_candidate, fw_candidate

                elif operation in ("gibson", "infusion"):
                    fw_tm_set = [] 
                    rv_tm_set = [] 
                    if deletion_mode:
                        for plen in range(primer_length[0], primer_length[1]):
                            fw_candidate = _slice_context(target_seq, loc[1], plen, template.topology == "circular")
                            rv_candidate = _slice_context(target_seq, loc[0]-plen, plen, template.topology == "circular").translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                            fw_tm = tm_func(seq=fw_candidate)
                            rv_tm = tm_func(seq=rv_candidate)
                            fw_tm_set.append([[fw_candidate, 0], fw_tm])
                            rv_tm_set.append([[rv_candidate, 0], rv_tm])
                        fw_adapter = _slice_context(target_seq, loc[0]-homology_length, homology_length, template.topology == "circular")
                        rv_adapter = _slice_context(target_seq, loc[1], homology_length, template.topology == "circular")
                    else:
                        for plen in range(primer_length[0], primer_length[1]):
                            flen = int(homology_length/2) 
                            rlen = int(homology_length/2)
                            fw_candidate = target_seq[mut_pos+flen:mut_pos+flen+plen] 
                            rv_candidate = target_seq[mut_pos-rlen-plen:mut_pos-rlen].translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
                            fw_tm = tm_func(seq=fw_candidate)
                            rv_tm = tm_func(seq=rv_candidate)
                            fw_tm_set.append([[fw_candidate, 0], fw_tm])
                            rv_tm_set.append([[rv_candidate, 0], rv_tm])
                        fw_adapter = mutate_seq[mut_pos-rlen:mut_pos+flen]
                        rv_adapter = mutate_seq[mut_pos-rlen:mut_pos+flen]
                    
                    if flip == 1:
                        fw_tm_set, rv_tm_set       = rv_tm_set, fw_tm_set
                        fw_candidate, rv_candidate = rv_candidate, fw_candidate
                        fw_adapter, rv_adapter     = rv_adapter, fw_adapter
                        #adapter_mode = operation
                else:
                    raise ValueError("If locations for mutations are single, operation should be 'gibson', 'infusion', 'QuickChange', 'Q5'.")

            else:
                pass
        else:
            amp_start = amplicon_region.seq.find(target.seq) 
            amp_end   = amp_start + len(target.seq)
            while True:
                fw_candidates = [] 
                if fw_primer is None:
                    for pos in range(amp_start+1):
                        for plen in range(primer_length[0], primer_length[1] + 1): 
                            fw_candidate = amplicon_region.seq[pos:pos+plen]
                            fw_candidates.append([str(fw_candidate), pos]) 
                else:
                    site = amplicon_region.searchsequence(query=fw_primer, quinable=False)
                    fw_candidates.append([fw_primer.seq, site.start])

                rv_candidates = [] 
                if rv_primer is None:
                    rv_window = len(amplicon_region.seq) - amp_end
                    for pos in range(rv_window + 1):
                        for plen in range(primer_length[0], primer_length[1] + 1): 
                            rv_candidate = amplicon_region.rcseq[pos:pos+plen]
                            rv_candidates.append([str(rv_candidate), pos]) 
                else:
                    site = amplicon_region.searchsequence(query=rv_primer, quinable=False)
                    rv_candidates.append([rv_primer.seq, len(amplicon_region.seq) - site.end])

                fw_candidates_total = len(fw_candidates)
                checked_fw_candidates = [] 
                for candidate in fw_candidates:
                    sites = template.searchsequence(query="(?:{}){{s<={}}}".format(candidate[0], nonspecific_limit), quinable=False) 
                    if len(sites) > 1:
                        pass 
                    else:
                        checked_fw_candidates.append(candidate) 
                        
                rv_candidates_total = len(rv_candidates)
                checked_rv_candidates = [] 
                for candidate in rv_candidates:
                    sites = template.searchsequence(query="(?:{}){{s<={}}}".format(candidate[0], nonspecific_limit), quinable=False) 
                    if len(sites) > 1:
                        pass
                    else:
                        checked_rv_candidates.append(candidate) 
                
                fw_candidates = checked_fw_candidates
                rv_candidates = checked_rv_candidates
                if len(fw_candidates) == 0 or len(rv_candidates) == 0:
                    grew = False
                    if len(fw_candidates) == 0 and fw_margin_auto and current_fw_margin < auto_margin_max:
                        current_fw_margin = min(current_fw_margin + auto_margin_step, auto_margin_max)
                        grew = True
                    if len(rv_candidates) == 0 and rv_margin_auto and current_rv_margin < auto_margin_max:
                        current_rv_margin = min(current_rv_margin + auto_margin_step, auto_margin_max)
                        grew = True

                    if grew:
                        start = target_start - current_fw_margin
                        if start < 0:
                            if template.topology == "circular":
                                start = len(template.seq) + start
                            else:
                                start = 0

                        end = target_start + len(target.seq) + current_rv_margin
                        if end > len(template.seq):
                            if template.topology == "circular":
                                end = end - len(template.seq) * (end // len(template.seq))
                        amplicon_region = template[start:end]
                        amp_start = amplicon_region.seq.find(target.seq)
                        amp_end = amp_start + len(target.seq)
                        continue

                    raise ValueError(
                        _primer_site_failure_message(
                            template=template,
                            target=target,
                            amplicon_region=amplicon_region,
                            fw_candidates_total=fw_candidates_total,
                            rv_candidates_total=rv_candidates_total,
                            fw_candidates_kept=len(fw_candidates),
                            rv_candidates_kept=len(rv_candidates),
                            primer_length=primer_length,
                            fw_margin=(f"auto->{current_fw_margin}" if fw_margin_auto else current_fw_margin),
                            rv_margin=(f"auto->{current_rv_margin}" if rv_margin_auto else current_rv_margin),
                            nonspecific_limit=nonspecific_limit,
                        )
                    )
                break
            
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
        primer_pairs.sort(key=lambda x: _primer_pair_sort_key(x, target_tm))
        filtered_primer_pairs = [] 
        for primer_pair in primer_pairs:
            if requirement(primer_pair): 
                filtered_primer_pairs.append(primer_pair)  
            else:
                pass
        
        filtered_primer_pairs = filtered_primer_pairs[:design_num]
        gap_is_empty = gap is None or (type(gap) in (tuple, list) and all(g is None for g in gap))
        if gap_is_empty:
            dammy = [{"fw":["ATGC", 4], "rv":["ATGC", 4], "fw_tm":50, "rv_tm":50}]
            gapinfo_fw = append_adapter(amplicon_region, dammy, fw_adapter, fw_partner, adapter_mode, int(homology_length/2), "fw", fw_name, None, auto_adjust)
            gapinfo_rv = append_adapter(amplicon_region, dammy, rv_adapter, rv_partner, adapter_mode, int(homology_length/2), "rv", rv_name, None, auto_adjust)
            gap_fw = gapinfo_fw[1] if type(gapinfo_fw) in (tuple, list) and len(gapinfo_fw) == 2 and gapinfo_fw[0] == "fw" else None
            gap_rv = gapinfo_rv[1] if type(gapinfo_rv) in (tuple, list) and len(gapinfo_rv) == 2 and gapinfo_rv[0] == "rv" else None
            filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, fw_adapter, fw_partner, adapter_mode, int(homology_length/2), "fw", fw_name, gap_fw, auto_adjust)
            filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, rv_adapter, rv_partner, adapter_mode, int(homology_length/2), "rv", rv_name, gap_rv, auto_adjust)
        else:
            filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, fw_adapter, fw_partner, adapter_mode, int(homology_length/2), "fw", fw_name, gap[0], auto_adjust)
            filtered_primer_pairs = append_adapter(amplicon_region, filtered_primer_pairs, rv_adapter, rv_partner, adapter_mode, int(homology_length/2), "rv", rv_name, gap[1], auto_adjust)

        for i in range(len(filtered_primer_pairs)):
            filtered_primer_pairs[i]["fw"] = filtered_primer_pairs[i]["fw"][0] 
            filtered_primer_pairs[i]["fw"].setfeature({"feature_type":"primer_bind", "qualifier:label":fw_name})
            filtered_primer_pairs[i]["rv"] = filtered_primer_pairs[i]["rv"][0]
            filtered_primer_pairs[i]["rv"].setfeature({"feature_type":"primer_bind", "qualifier:label":rv_name})
        
        return filtered_primer_pairs
