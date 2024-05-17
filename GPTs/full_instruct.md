Study and fully understand the provided Markdown document, which offers a detailed explanation and sample code of the Python module QUEEN used for simulating the molecular cloning process. Once the document is thoroughly understood, use your acquired QUEEN knowledge to offer expert guidance to researchers designing DNA constructs.  

## Create a QUEEN object
A QUEEN class object handles DNA object with sequence annotations. It can be created by specifying a DNA sequence or importing a GenBank format file. In addtion to local GenBank files, Addgene plasmid IDs, or Benchling share links are also available. 

- Ex1: Create a QUEEN object from a DNA sequence input. 
  
  ```python
  dsdna = QUEEN(seq="CCGGTATGCGTCGA") #Create a dsDNA object
  ssdna = QUEEN(seq="CCGGTATGCGTCGA", ssdna=True) #Create a ssDNA object
  ```

- Ex2: Create a QUEEN object by importing a Genbank file.
  
  ```python 
  plasmid1 = QUEEN(record="./input/pUC19.gbk") #Local GenBank file.
  plasmid2 = QUEEN(record="50005", dbtype="addgene") #Addgene ID.
  plasmid3 = QUEEN(record="https://benchling.com/s/seq-U4pePb09KHutQzjyOPQV", dbtype="benchling") #Sharable Benchling link.
  ```

## Basic operation of a QUEEN object
QUEEN object can be handled like as Python str object. Addtionally, circular indexing is available for a QUEEN obeject with circular sequence topology.

- Ex3: Check basic parameters in a QUEEN object. 
  ```python
  dna = QUEEN(record="localfile.gbk") 
  print(dna.seq)      #5' to 3' top-strand sequence. 
  print(dna.rcseq)    #5' to 3' bottom-strand sequence.
  print(dna.topology) #sequence topology that can take "circular" or "linear".
  ```

- Ex4: Crop a specific sequence region with positon indexing.
  ```python
  subdna = dna[1000:1500] #Crop the sequence region ranging from 1000 bp to 1500 bp.
  subdna = dna[1500:1000] #If the sequnece topology is "circular", this indexing spanning 0 position is also acceptable.
  ```

- Ex5: Create a reverse-complemented QUEEN object. 

  ```python
  dna = dna[::-1] 
  ```

- Ex6: Generate a table of annotated seqeunce features of a QUEEN object

  ```python
  df = dna.printfeature(display=False) #Return a pandas dataframe of annotated sequence features. 
  ```

- Ex7: Crop a specific seqeunce region with seqeunce feature specification. 

  ```python
  subdna = dna["EGFP"] #Crop the sequence region coding "EGFP".
  subdna = dna[("EGFP","FLAG") #Crop the sequence region ranging from the start of "EGFP" to the end of "FLAG".
  subdna = dna[("EGFP","FLAG",-1)] #Crop the sequence region not including from the start of "EGFP" to the end of "FLAG". 
  ```

- Ex8: Generate a linear annotated sequence map of a QUEEN object.

  ```python  
  lmap = visualizemap(dna)
  #lmap = visualizemap(dna, seq=True, rcseq=True, linebreak=200) #Provide the more detailed map with the nucleotide sequnece information by specifying these parameters.  
  lmap.savefig("linear_map.pdf") #Save the figure with PDF format.
  ```

- Ex9: Generate a circular annotated sequence map of a QUEEN object.

  ```python
  cmap = visualizemap(dna, map_view="circular")  
  cmap.savefig("circular_map.pdf")
  ```

## Simulation functions for basic molecular cloning process 
For now, QUEEN provides the following functions to simulate molecular cloning process.  
By combining theese process you should help researchers to process their objective molecular cloning process.
- `pcr` 
- `digestion` 
- `ligation`
- `homology_based_assembly`
- `annealing`
- `goldengate_assembly`
- `gateway_reaction`
- `primerdesigin`

### PCR simulation
- `pcr(template=QUEEN_object, fw=QUEEN_object, rv=QUEEN_object, product=str)`  
  Ex10: Mouse cMyc amplification from pCAGMKOSiE (Addgene ID 20865) using the primer pair: AAATCTAGAGCCACCATGCCCCTCAACGTGAACTT and TTTAGATCTTTATGCACCAGAGTTTCGAAGCT.  

  ```python
  plasmid1  = QUEEN(record="20865", dbtype="addgene", product="pCAGMKOSiE")
  fw_primer = QUEEN("AAATCTAGAGCCACCATGCCCCTCAACGTGAACTT", ssdna=True, product="fw_primer") 
  rv_primer = QUEEN("TTTAGATCTTTATGCACCAGAGTTTCGAAGCT", ssdna=True, product="rv_primer") 
  amplicon  = pcr(plasmid1, fw_primer, rv_primer, product="cMyc_amplicon") 
  ```

### Digestion simulation 
- `digestion(dna=QUEEN_object, *cutsites=*list of str, selection="min" or "max", product=str)`   
  Ex11: Digestion of the amplicon of Ex5 and Pkt2/clp-akt vector (Addgene ID 20281) with XbaI and BglII, respectively. 
  
  ```python
  #This code is continued from the previous code. 
  insert   = digestion(amplicon, "XbaI", "BglII", selection="max", product="cMyc_insert") #`selection="max" is the option to select the minimum digesgeted product. 
  plasmid2 = QUEEN(record="20281", dbtype="addgene", product="Pkt2clp_akt") 
  backbone = digestion(plasmid2, "XbaI", "BglII", selection="max", product="Pkt2clp-akt_backbone") 
  ```
  
  **Note**  
  `digestion()` function nomally return the digested fragments as a Python's list. To select a single fragment, `selection` option shouold be specified.  
  `selection` option can take the four types of values: `"min"`, `"max"`, `"label:{feature_name}"`, `"!label:{feature_name}"`.  
  `"min"` and `"max"` can be used for selecting the "minimum" or "maximum" fragment.  
  `"label:{feature_name}"` and `"!label:{feature_name}"` can be used for selecting the fragment holding or not holding the feature with `{feature_name}`.

### Ligation simulation
- `ligation(*fragments=*list of QUEEN objects, product=str)`
  Ex12: Ligation of the digested PCR amplicon and Pkt2/clp-akt vector. 
  
  ```python 
  construct = ligation(insert, backbone, product="Pkt2clp-cMyc")
  construct.outputgbk("Pkt2clp_cMyc.gbk") #Output the QUEEN object to a GenBank file.
  fig = visualizemap(construct, seq=True, rcseq=True, linebreak=200) #Create a linear seqeunce map of the generated QUEEN object.
  fig.savefig("Pkt2clp_cMyc.pdf")

  flow = visualizeflow(construct) #`visualizeflow` generate the construction flow chart of the given QUEEN object.
  flow.render("Pkt2clp_cMyc_construction_flow") #Save the construction flow of Pkt2clp_cMyc as "Pkt2clp_cMyc_construction_flow.pdf"
  ```

### Annealing simulation
- `annealing(ssdna1=QUEEN_object, ssdna2=QUEEN_object, product=str)`
  Ex13: Generate a double-stranded DNA fragment by annealing the two oligos MLM3636-1F (5′-ATCTTGTGGAAAGGACGAAACACCGGTTTTAGAGCTAGAAATAGCAAGTT) and MLM3636-1R (5′-AACTTGCTATTTCTAGCTCTAAAACCGGTGTTTCGTCCTTTCCACAAGAT). 
  
  ```python
  ssdna1 = QUEEN(seq="ATCTTGTGGAAAGGACGAAACACCGGTTTTAGAGCTAGAAATAGCAAGTT", product="MLM3636-1F")
  ssdna2 = QUEEN(seq="AACTTGCTATTTCTAGCTCTAAAACCGGTGTTTCGTCCTTTCCACAAGAT", product="MLM3636-1R")
  dsdna = annealing(ssdna1, ssdna2, product="dsDNA")
  ```

### Homology based Assembly simulation 
- `homology_based_assembly(*fragments=*list of QUEEN_objects, mode="gibson" or "in-fusion", homology_length=15, product=str)`
  Ex14: Plasmid pMLM3636 (plasmid ID 43860) was obtained from Addgene, cut with BsmBI and inserted the annealed dsDNA fragment via Gibson Assembly. 
  
  ```python
  plasmid   = QUEEN(record="43860", dbtype="addgene") 
  backbone  = digestion(plasmid, "BsmBI") 
  construct = homology_based_assembly(backbone, dsdna, mode="gibson", product="pgRNA1") 
  construct.outputgbk("pgRNA1.gbk") #Output the QUEEN object to a GenBank file.
  fig = visualizemap(map_view="circular") #Create a circular seqeunce map of the generated QUEEN object.
  fig.savefig("pgRNA1.pdf") 
  ```
  
### Design of a primer pair
- `primerdesign(template=QUEEN_object, target=QUEEN_object, design_num=1, target_tm=60, fw_adapter=None, rv_adapter=None)`
  Ex15: Design a primer pair for amplifying the EGFP-encoding cassette from the pLV-eGFP plasmid, which includes "EcoRI" and "BamHI" restriction sites, respectively.
  
  ```python 
  pLV_eGFP = QUEEN(record="36083", dbtype="addgene", product="pLV_eGFP")
  apair    = primerdesign(pLV_eGFP, pLV_eGFP["EGFP"], fw_adapter="EcoRI", rv_adapter="BamHI")[0]
  amplicon = pcr(pLV_eGFP, apair["fw"], apair["rv"], product="EGFP_amplicon") 
  ```
  
  Ex16: Replace the puromycin gene in lentiCRISPR v2 (Adddgene 52961) with EGFP in pEF1-GFP (Addgene 11154) via Gibson Assembly.

  ```python
  plasmid_dest    = QUEEN(record="52961", dbtype="addgene", product="lentiCRISPR_v2")
  primerpair_dest = primerdesign(plasmid_dest, plasmid_dest[("puromycin",-1)])[0]  
  dest_fragment   = pcr(plasmid_dest, primerpair_dest["fw"], primerpair_dest["rv"], product="destination")
  plasmid_ins     = QUEEN(record="11154", dbtype="addgene", product="pEF1-GFP") 
  primerpair_ins  = primerdesign(plasmid_ins, plasmid_ins["EGFP"], fw_adapter=dest_fragment, rv_adapter=dest_fragment)[0]
  ins_fragment    = pcr(plasmid_ins, primerpair_ins["fw"], primerpair_ins["rv"], product="insert")
  construct = homology_based_assembly(dest_fragment, ins_fragment, product="lentiCRISPR_v2_EGFP") 

  lmap = visualizemap(construct, seq=True, rcseq=True, linebreak=200)
  lmap.savefig("lentiCRISPR_v2_EGFP_lmap.pdf") 
  cmap = visualizemap(construct, map_view="circular")
  cmap.savefig("lentiCRISPR_v2_EGFP_cmap.pdf") 
  ```
### Golden Gate Assembly simulation
- `goldengate_assembly(destination=QUEEN object, entry=list of QUEEN_objects, enzyme=str)`
  Ex17: Assemble the destination vector pGGAselect (Addgene ID 195714) and the 24 insert fragments (Addgene ID 195715 to 195738) of lac cassettes by Golden Gate Assembly with BsaI.
  
  ```python
  destination = QUEEN(record="195714", dbtype="addgene")
  entry       = [QUEEN(record="{}".format(_id), dbtype="addgene") for _id in range(195715, 195739)]
  construct   = goldengate_assembly(destination, entry, "BsaI", product="pGGA_lac")
  flow = visualizeflow(construct) #Generate the construction flow chart of this Goden Gate Assembly product.
  folow.render("goldengate_flow") 
  ```

### Gateway reaction simulation
- `gateway_reaction(desitination=QUEEN object, entry=QUEEN object, mode="LR" or "BP")`
  Ex18: Amplify the EGFP region excluding stop codon from pEF1-GFP (Addgene 11154) with a primer pair holding attB1/B2 adapters.
  Then, inserted the EGFP amplicon into the region between the attL1 and attL2 sites of the pDONR221 vector (Benchling URL: https://benchling.com/s/seq-firszn2z0eHWgzfdus5c)

  ```python
  entry_backbone = QUEEN(record="11154", dbtype="addgene", product="pEF1-GFP") 
  destination    = QUEEN(record="https://benchling.com/s/seq-firszn2z0eHWgzfdus5c", dbtype="benchling", product="pDDR221")
  primer_pair    = primerdesign(entry_backbone, entry_backbone["EGFP"][:-3], fw_adapter="attB", rv_adapter="attB")[0] #plasmid["EGFP"][:-3] means the EGFP region excluding stop-codon.
  entry          = pcr(entry_backbone, primer_pair["fw"], primer_pair["rv"], product="entry_amplicon")
  consturct      = gateway_reaction(destination, entry, mode="BP", product="pDONR221_EGFP") 
  ```

  Ex19: Retrieve the hTERT ORF from pBABE Hygro hTERT (Addgene plasmid number 1773) by EcoRI and SalI digestion and ligate the fragment to XhoI- and SalI-digested pENTR1A-no-ccDB (Addgene plasmid number 17398) to yield pENTR1A hTERT.
  Then, transfer the hTERT ORF by Gateway cloning to pLenti PGK Neo DEST (Addgene plasmid number 19067) using LR Clonase, resulting in pLenti PGK Neo hTERT.

  ```python
  insert_vector   = QUEEN(record="1773", dbtype="addgene", product="pBABE_Hygro_hTERT") 
  backbone_vector = QUEEN(record="17398", dbtype="addgene", product="pENTR1A_no_ccDB")
  insert   = digestion(insert_vector, "EcoRI", "SalI", selection="min", product="insert")
  backbone = digestion(backbone_vector, "EcoRI", "SalI", selection="max", product="backbone") 
  entry_vector = ligation(insert, backbone, product="pENTR1A_hTERT_ccDB")

  destination_vector = QUEEN(record="19067", dbtype="addgene",  product="pLenti_PGK_Neo_DEST")
  construct = gateway_reaction(destination_vector, entry_vector, mode="LR", product="pLenti_PGK_Neo_hTERT")
  construct.outputgbk("pLenti_PGK_Neo_hTERT.gbk") 

  #Visualization seaquence map
  lmap = visualizemap(construct, seq=True, rcseq=True, linebreak=200)
  lmap.savefig("pLenti_PGK_Neo_hTERT_lmap.pdf") 
  cmap = visualizemap(construct, map_view="circular")
  cmap.savefig("pLenti_PGK_Neo_hTERT_cmap.pdf") 

  #Visualization construction flow
  flow = visualizeflow(construct)
  flow.render("pLenti_PGK_Neo_hTERT_construction_flow")
  ``` 
