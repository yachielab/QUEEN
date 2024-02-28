# Usage of `qexperiment` module


## PCR simulation
Simulate PCR (Polymerase Chain Reaction) on a given DNA template using specific primers.

#### `pcr(template, fw, rv, bindnum=16, mismatch=1, endlength=3, return_template=False, product=None, process_description=None, pd=None)`

#### Parameters

- **template**: `QUEEN`  
  DNA template to be amplified.

- **fw**: `QUEEN` or `str`  
  Forward primer. Can be a `QUEEN` object or a string.

- **rv**: `QUEEN` or `str`  
  Reverse primer. Can be a `QUEEN` object or a string.

- **bindnum**: `int`, optional  
  Minimum binding nucleotides for a primer. Default is 16.

- **mismatch**: `int`, optional  
  Maximum mismatches allowed in primer binding. Default is 1.

- **endlength**: `int`, optional  
  Length of primer's end region for binding consideration. Default is 3.

- **return_template** `bool`, optional  
  Whether to return the modified template. Default is False.

- **product** `str`, optional  
  Product name.

- **process_description** or **pd**:  `str`, optional  
  Additional process description.

#### Returns
- `QUEEN object(amplicon)` or `QUEEN object(template)` and `QUEEN(amplicon)`  
  PCR product or both modified template and product.

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

#### `digestion(dna, *cutsites, size_selection=None, requirement=None, product=None, process_description=None, pd=None)`

#### Parameters
- **dna**: `QUEEN`  
  DNA sequence to be digested. 

- **cutsites**: `Cutsite` or `str`  
  Restriction enzymes for digestion.

- **size_selection**: `str` or `tuple` of `int`, optional  
  Criteria for fragment size selection.

- **product**: `str`, optional  
  Product name.

- **process_description** or **pd** (`str`, optional)  
  Additional process description.

#### Returns
- List of `QUEEN object` objects or a single `QUEEN object`  
  Depending on size selection, digested fragments.

#### Example Usage
```python
>>> dna_sequence = QUEEN("example_dna_sequence")
>>> cutsite1 = Cutsite("restriction_enzyme_1")
>>> cutsite2 = Cutsite("restriction_enzyme_2")
>>> fragments = digestion(dna_sequence, cutsite1, cutsite2, size_selection=(100, 1000))
```

---
## DNA Ligation Simulation
Simulate the ligation of DNA fragments into unique or multiple constructs.

#### `ligation(*fragments, unique=True, product=None, process_description=None, pd=None):`

#### Parameters
- **fragments**: `QUEEN object`  
  DNA fragments to be ligated.

- **unique**: `bool`, optional  
  Whether to return only a unique construct. Default is True.

- **product**: `str`, optional  
  Product name.

- **process_description** or **pd** :`str`, optional  
  Additional process description.

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

#### `homology_based_assembly(*fragments, mode="gibson", homology_length=20, unique=True, product=None, process_description=None, pd=None)`

#### Parameters
- **fragments**: `QUEEN object`  
  DNA fragments for assembly.

- **mode**: `str`, optional  
  Assembly mode. Default is "gibson".  
  Valid options are "gibson", "infusion", and "overlappcr". Default is "gibson".

- **homology_length**: `int`, optional  
  Required homology length. Default is 20.

- **unique**: `bool`, optional  
  Return only a unique construct. Default is True.

- **product**: `str`, optional  
  Process name.

- **process_description** or **pd**: `str`, optional  
  Additional description.

#### Returns
- `QUEEN object` or list of `QUEEN object`  
  Construct(s) resulted from assembly.

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

#### `annealing(ssdna1, ssdna2, homology_length=4, product=None, pd=None, process_description=None)`

#### Parameters
- **ssdna1** & **ssdna2**: `QUEEN object`  
  Single-stranded DNAs for annealing.

- **homology_length**: `int`, optional  
  Required homology length. Default is 4.

- **product**: `str`, optional  
  Product name.

- **process_description** or **pd**: `str`, optional  
  Additional description.

#### Returns
- `QUEEN object`  
  Double-stranded DNA molecule after annealing.

#### Example Usage

```python
>>> ssdna1 = QUEEN("ATCG")
>>> ssdna2 = QUEEN("CGTA")
>>> dsdna = annealing(ssdna1, ssdna2, homology_length=4)
```

---

## Primer Design for PCR Amplification

Design forward and reverse primers for PCR amplification of a target region, allowing for introduction of specific mutations, checking primer specificity, and meeting additional user-defined requirements.

#### `primerdesign(template, target, fw_primer=None, rv_primer=None, fw_margin=0, rv_margin=0, target_tm=60.0, tm_func=None, primer_length=(16, 25), design_num=1, fw_adapter=None, rv_adapter=None, homology_length=20, nonspecific_limit=3, requirement=lambda x: x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T"), fw_name="fw_primer", rv_name="rv_primer")`

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
  Adapter sequence to prepend to any designed forward primer.

- **rv_adapter**: `ssDNA QUEEN object` or `str`, optional  
  Adapter sequence to prepend to any designed reverse primer.

- **homology_length** `int`, optional  
  Required homology for adapters.  

- **nonspecific_limit** (`int`, optional): 
  The maximum number of mismatches allowed for primer binding outside of the designated primer design region within the template sequence.  
  Primer pairs that bind to any region of the template with a number of mismatches equal to or less than this limit will be excluded from the design, to increase the specificity of the PCR reaction and decrease the likelihood of nonspecific amplification. Default is  3.

- **requirement** (lambda function, optional):   
  Function that takes a dictionary representing a primer pair and returns True if the pair meets the specified conditions.  
  Ensures that the 3' end nucleotide of both primers is not A or T by default.   
  The detailed default function is `lambda x: x["fw"][-1] not in ("A", "T") and x["rv"][-1] not in ("A", "T")`

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
>>> # Assuming target_Q sequence is within template_Q
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
