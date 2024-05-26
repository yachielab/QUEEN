# Usage of `qexperiment` module

---

## PCR simulation

Simulate PCR (Polymerase Chain Reaction) on a given DNA template using specific primers.

#### `pcr(template, fw, rv, bindnum=16, mismatch=1, endlength=3, return_template=False, product=None, product=None, process_description=None, pn=None, pd=None)`

#### Parameters

- **template**: `QUEEN object`  
  DNA template to be amplified.

- **fw**: `QUEEN object` or `str`  
  Forward primer. Can be a `QUEEN` object or a string.

- **rv**: `QUEEN object` or `str`  
  Reverse primer. Can be a `QUEEN` object or a string.

- **bindnum**: `int`, optional  
  Minimum binding nucleotides for a primer. Default is 16.

- **mismatch**: `int`, optional  
  Maximum mismatches allowed in primer binding. Default is 1.

- **endlength**: `int`, optional  
  Length of primer's end region for binding consideration. Default is 3.

- **add_primerbind** `bool`, optional  
  If True, add DNAfeature on the primer binding regions in the template DNA. 

- **tm_func** `function`, optional
  Function to calculate the melting temperature of the primer pair. 
  As `str` specfication, you can select `"Breslauer"` and `"SantaLucia"`. Default is `"SantaLucia"`.  
  Also, as built-in algorithms, `QUEEN.qexperiment.Tm_NN()`. This function is implemented 
  based on the `Bio.SeqUtils.MeltingTemp.Tm_NN()`, so the all parameters of 
  `Bio.SeqUtils.MeltingTemp.Tm_NN()`, excluding `seq` and `c_seq`, can be acceptable.

- **return_tm** : `bool`, optional
  If True, tm values of the primer pair are also returned.    
  The tm values will be calculated based on thier biding region excluding adaptor regions.

- **product** `str`, optional  
  Product name.

- **process_name** or **pn**:  `str`, optional  
  Brief label for the `pcr` process. Default is `"PCR"`.

- **process_description** or **pd**:  `str`, optional  
  Additional process description.

- **pn** `str`, optional
  Alias for `process_name`.

- **pd** : `str`, optional
  Alias for `process_description`.


#### Returns

- `QUEEN object` if `return_tm` is False, Otherwise `QUEEN object, (tm_fw, tm_rv)`
  Returns the PCR product (amplicon) as a QUEEN object. 
  If `return_tm` is True, the Tm values of the given primer pair is also returned.  

#### Example Usage

```python
>>> template = QUEEN("example_dna_sequence")
>>> forward_primer = "xxxxxxxxx"
>>> reverse_primer = "xxxxxxxxx"
>>> pcr_product = pcr(template, forward_primer, reverse_primer)
```

---

## Digestion simulation

Simulate DNA digestion using restriction enzymes and optionally perform size selection.

#### `digestion(dna, *cutsites, size_selection=None, requirement=None, product=None, process_name=None, process_description=None, pn=None, pd=None)`

#### Parameters

- **dna**: `QUEEN object`  
  DNA sequence to be digested. 

- **cutsites**: `Cutsite` or `str`  
  Restriction enzymes for digestion.

- **size_selection**: `"min"`, `"max"`, `"label:*"`, `"!label:*"`, or `tuple` of `int`, optional
  Criteria for fragment selection.
  If `"min"` is provided, the minimum fragment of the digested fragments would be returned.  
  If `"max"` is provided, the maximum fragment of the digested fragments would be returned.  
  If `"label:{feature_of_interest}"` is provided, the unique fragment holding the DNAfeature with `feature_of_interest` in "qualifer:label" would be returned. If multiple fragments holding the specified feature are detected, a error will be raised.  
  If `"!label:{feature_of_interest}"` is provided, the unique fragment not holding the DNAfeature with `feature_of_interest` in "qualifer:label" would be returned. If multiple not fragments holding the specified feature are detected, a error will be raised.
  If a `tuple` value is provided, The tuple (`min_size`, `max_size`) specifies the size range for filtering the resulting fragments. If multiple fragments holding the are detected in the specified range, a error will be raised.  
  If `None`, no filtering is done. Default is `None`.

- **product**: `str`, optional  
  Product name.

- **process_name** or **pn**:  `str`, optional  
  Brief label for the `digestion` process. Default is "Digestion".

- **process_description** or **pd** (`str`, optional)  
  Additional process description.

- **pn** `str`, optional
  Alias for `process_name`.

- **pd** : `str`, optional
  Alias for `process_description`.

#### Returns

- `list` of `QUEEN objects` or `QUEEN object` 
  If `selection` is None, return list of QUEEN objects composed of the digested fragments.  
  Otherwise, return a specific fragment filling the specified condition.

#### Example Usage

```python
>>> dna_sequence = QUEEN("example_dna_sequence")
>>> fragments = digestion(dna_sequence, cs.lib["BamHI"], cutsite["AgeI"], size_selection=(100, 1000))
```

---

## DNA Ligation Simulation

Simulate the ligation of DNA fragments into unique or multiple constructs.

#### `ligation(*fragments, unique=True, follow_order=False, product=None, process_name=None, process_description=None, pn=None, pd=None):`

#### Parameters

- **fragments**: `list` of `QUEEN object`  
  DNA fragments to be ligated.

- **unique**: `bool`, optional  
  Whether to return only a unique construct. Default is True.

- **follow_order** : `bool`, optional 
  If True, a ligation reaction will be simulated along with the given order of fragments.  
  Default is False. 

- **product**: `str`, optional  
  Product name.

- **process_name** or **pn**:  `str`, optional  
  Brief label for the ligation process. Default is "Ligation".

- **process_description** or **pd** :`str`, optional  
  Additional process description.

- **pn** `str`, optional
  Alias for `process_name`.

- **pd** : `str`, optional
  Alias for `process_description`.

#### Returns

- ｀QUEEN object｀ or `list` of `QUEEN objects`  
  Construct(s) resulted from ligation.  
  If `unique` is True and only one construct is possible, returns that construct as a QUEEN object.   
  If `unique` is False, returns a list of all possible assembled constructs as QUEEN objects.  
  If no constructs are possible, returns an empty list or None.  

#### Example Usage

```python
>>> fragment1 = QUEEN("dna_fragment_1")
>>> fragment2 = QUEEN("dna_fragment_2")
>>> assembled_product = ligation(fragment1, fragment2, unique=True)
```

#### Notes

The function attempts all permutations and orientations of the given fragments for ligation.  
If `unique` is set to True, it validates the uniqueness of the assembled product.   
The function handles situations where the assembly is not possible or results in multiple products.

---

## Homology-Based DNA Assembly Simulation

Simulate DNA assembly using homology for various assembly modes.

#### `homology_based_assembly(*fragments, mode="gibson", homology_length=20, unique=True, follow_order=None, product=None, process_name=None, process_description=None, pn=None, pd=None)`

#### Parameters

- **fragments**: `QUEEN object`  
  DNA fragments for assembly.

- **mode**: `str`, optional  
  Assembly mode. Valid options are `"gibson"`, `"infusion"`, anb `"overlappcr"`. Default is `"gibson"`.

- **homology_length**: `int`, optional  
  Required homology length. Default is 20.

- **unique**: `bool`, optional  
  Return only a unique construct. Default is True.

- **follow_order** : `bool`, optional 
  If True, a ligation reaction will be simulated along with the given order of fragments. 
  If the number of given fragments is larger than 4, default is True. Otherwise, False.   

- **product**: `str`, optional  
  Product name.

- **process_name** or **pn**:  `str`, optional  
  Brief label for the `homology_based_assembly` process. 
  If `mode` is `"gibson"`, default is `"Gibson Assembly"`. If `mode` is `"infusion"`, default is `"In-Fusion Assembly"`. 

- **process_description** or **pd**: `str`, optional  
  Additional description.

- **pn** `str`, optional
  Alias for `process_name`.

- **pd** : `str`, optional
  Alias for `process_description`.


#### Returns

- `QUEEN object` or list of `QUEEN object`  
  Construct(s) resulted from assembly.   
  If `unique` is True and only one construct is possible, returns that construct as a QUEEN object.  
  If `unique` is False, returns a list of all possible assembled constructs as QUEEN objects.   
  If no constructs are possible, returns an empty list or None.

#### Example Usage

```python
>>> fragment1 = QUEEN("dna_fragment_1")
>>> fragment2 = QUEEN("dna_fragment_2")
>>> assembled_product = homology_based_assembly(fragment1, fragment2, mode="gibson", unique=True)
```

#### Notes

The function considers different assembly modes, each with its specific requirements for fragment or ientation and homology lengths. It handles permutations and orientations of the given fragments.

---

## DNA Annealing Simulation

Simulate the annealing of two single-stranded DNA molecules based on homology.

#### `annealing(ssdna1, ssdna2, homology_length=4, product=None, process_name=None, process_description=None, pn=None, pd=None)`

#### Parameters

- **ssdna1** & **ssdna2**: `QUEEN object`  
  Single-stranded DNAs for annealing.

- **homology_length**: `int`, optional  
  Required homology length. Default is 4.

- **product**: `str`, optional  
  Product name.

- **process_name** or **pn**:  `str`, optional  
  Brief label for the `annealing` process. Default is `"Annealing"`. 

- **process_description** or **pd**: `str`, optional  
  Additional description.

- **pn** `str`, optional
  Alias for `process_name`.

- **pd** : `str`, optional
  Alias for `process_description`.

#### Returns

- `QUEEN object`  
  Double-stranded DNA molecule after annealing.

#### Example Usage

```python
>>> ssdna1 = QUEEN("ATCG")
>>> ssdna2 = QUEEN("CGAT")
>>> dsdna = annealing(ssdna1, ssdna2)
```

---

## Gateway Reaction Simulation 

Simulates a gateway reaction of two DNA molecules. For now, basic `BP` and `LR` reactions are available.  

#### `gateway_reaction(destination, entry, mode="BP", product=None, process_name=None, process_description=None, pn=None, pd=None)`

#### Parameters

- **destination** : `QUEEN object`
  The destination QUEEN object holding the backbone DNA molecule.

- **entry** : `QUEEN object`
  The entry QUEEN object holding the insert DNA molecule.

- **mode**: `str`, `tuple`, or `list
  The mode of the reaction, can be "BP" or "LR". Default is "BP".  
  For executing a custom BP or LR reaction, please speicy `[B1 or L1 sequence, B1 or L2 sequence, P1 or R1 sequnce, P2 or R2 sequnece]` along with the QUEEN's cutsite format."
- **process_name** or **pn**:  `str`, optional  
  Brief label for the `annealing` process. Default is `"Gateway Reaction"`. 

- **process_description** or **pd**: `str`, optional  
  Additional description.

- **pn** `str`, optional
  Alias for `process_name`.

- **pd** : `str`, optional
  Alias for `process_description`.

#### Returns

- `QUEEN object`
  The QUEEN object representing the result of the gateway reaction process.

## Golden Gate Assembly Simulation 

Simulates a Golden Gate Assembly.

#### `goldengate_assembly(destination, entry, cutsite=None, product=None, process_name=None, process_description=None, pn=None, pd=None, **kwargs)`

#### Parameters

- **destination** : `QUEEN object`
  The destination QUEEN object holding the backbone DNA molecule.

- **entry** : `list` of `QUEEN objects`
  The entry QUEEN object(s) holding the insert DNA molecules.

- **cutsite** : `Cutsite` or `str`
  The restriction enzyme site used for this reaction.    

- **process_name** : str, optional
  Brief label for the gateway reaction process. Default is "Golden Gate Assembly".

- **process_description** : str, optional
  Additional description for the gateway reaction process.

- **pn** : str, optional
  Alias for `process_name`.

- **pd** : str, optional
  Alias for `process_description`.

#### Returns

- `QUEEN object`
  The QUEEN object representing the result of the Golden Gate Assembly.

---

## Primer Design for PCR Amplification

Design forward and reverse primers for PCR amplification of a target region, allowing for introduction of specific mutations, checking primer specificity, and meeting additional user-defined requirements.

#### `primerdesign(template, target, fw_primer=None, rv_primer=None, fw_margin=0, rv_margin=0, target_tm=60.0, tm_func=None, primer_length=(16, 25), design_num=1, fw_adapter=None, rv_adapter=None, homology_length=20, nonspecific_limit=3, requirement=None, fw_name="fw_primer", rv_name="rv_primer")`

#### Parameters

- **template** & **target**: `QUEEN object`  
  Template and target regions for PCR.

- **fw_primer**:`ssDNA QUEEN object`, optional  
  A specific forward primer sequence.   
  If provided, this sequence will be used as the forward primer.

- **rv_primer** :`ssDNA QUEEN`, optional  
  A specific reverse primer sequence.   
  If provided, this sequence will be used as the reverse primer.

- **fw_margin** & **rv_margin**: `int`, optional  
  Margins for primer design around the target region.  

- **target_tm** `float`, optional   
  Desired primer Tm. Default is 60.0 °C.

- **tm_func**  `function`, optional  
  Function for calculating primer Tm.  As built-in algorithms, `QUEEN.qexperiment.Tm_NN()`.   
  This function is implemented based on the `Bio.SeqUtils.MeltingTemp.Tm_NN()`, so the all parameters of  `Bio.SeqUtils.MeltingTemp.Tm_NN()`, excluding `seq` and `c_seq`, can be acceptable.

- **primer_length**: `tuple` of `int`, optional  
  A tuple (min_size, max_size) specifying the primer length. Default is (16, 25).

- **design_num**: `int`, optional   
  Number of primer pairs to design. Default is 1.

- **fw_adapter** `ssDNA QUEEN object` or `str`, optional  
  If it's a string or a single-stranded DNA (ssDNA) QUEEN object, the sequence will be added at the beginning of any forward primers designed.  
  If it's a double-stranded DNA (dsDNA) QUEEN object with linear topology, the homology sequence to the 3' end of the QUEEN object will be added at the beginning of any forward primers.  
  Alternatively, you can specify the name of a restriction enzyme or 'attB'. In these cases, the adapter sequence including the specified site will be added at the beginning of the primers.  
  For now, "attB" sequence as fw adapter is "GGGGACAAGTTTGTACAAAAAAGCAGGCT".

- **rv_adapter**: `ssDNA QUEEN object` or `str`, optional  
  If it's a string or a single-stranded DNA (ssDNA) QUEEN object, the sequence will be added at the beginning of any reverse primers designed.  
  If it's a double-stranded DNA (dsDNA) QUEEN object with linear topology, the homology sequence to the 3' end of the QUEEN object will be added at the beginning of any reverse primers.  
  Alternatively, you can specify the name of a restriction enzyme or 'attB'. In these cases, the adapter sequence including the specified site will be added at the beginning of the primers.  
  For now, "attB" sequence as rv adapter is "GGGGACCACTTTGTACAAGAAAGCTGGGT".  

- **homology_length** `int`, optional  
  Required homology for adapters.  

- **nonspecific_limit** (`int`, optional): 
  The maximum number of mismatches allowed for primer binding outside of the designated primer design region within the template sequence.  
  Primer pairs that bind to any region of the template with a number of mismatches equal to or less than this limit will be excluded from the design, to increase the specificity of the PCR reaction and decrease the likelihood of nonspecific amplification. Default is  3.

- **requirement** (lambda function, optional):   
  Function that takes a dictionary representing a primer pair and returns True if the pair meets the specified conditions.  
  Ensures that the 3' end nucleotide of both primers is not A or T by default.   
  The detailed default function is `lambda x: x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T")`
  Default requirement is as follows.
  	- `x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T")`
	- `"AAAA" not in x["fw"] and "TTTT" not in x["fw"] and "GGGG" not in x["fw"] and "CCCC" not in x["fw"]`
	- `"AAAA" not in x["rv"] and "TTTT" not in x["rv"] and "GGGG" not in x["rv"] and "CCCC" not in x["rv"]`

- **fw_name**: `str`, optional  
  Foward primer name, Default is "fw_primer". 

- **rv_name**: `str`, optional  
  Reverse primer name, Default is "rv_primer". 

#### Returns

- `list` of `dict`
  A list of dictionaries where each dictionary represents a primer pair.  
  Each dictionary contains two keys, "fw" and "rv", with the corresponding primer sequences  formed by ssDNA QUEEN objects.   
  The list is sorted by the closeness of the primer's Tm to the target_tm, with the closest pair first.

#### Example Usage

```python
>>> from QUEEN.queen import *
>>> template_Q = QUEEN('ATGC...')
>>> target_Q = QUEEN('ATGC...')
>>> #Assuming target_Q sequence is within template_Q
>>> primers = primerdesign(template_Q, target_Q, target_tm=65.0, num_design=5)
>>> primers
[
    {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
    {"fw": QUEEN(seq="ATGCGT...", ssdna=True), "rv": QUEEN(seq="TACGCA...", ssdna=True)},
    ...
]
```

#### Notes

The requirement for the target sequence to be within the template sequence ensures specificity of the primers to the region of interest. The function will not proceed if the target sequence is not a subset of the template.

---

## Under implementation
- intra_site_specific_recombination
- homologous_recombination
