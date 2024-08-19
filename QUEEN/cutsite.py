import os 
import sys 
import regex as re 
sys.path.append("/".join(__file__.split("/")[:-1]))
from qseq import Qseq

def compilecutsite(site):
    ref1 = r"([ATGCRYKMSWBDHVN]+)\([\-0-9]+/[\-0-9]+\)"
    ref2 = r"\([\-0-9]+/[\-0-9]+\)([ATGCRYKMSWBDHVN]+)\([\-0-9]+/[\-0-9]+\)"
    ref3 = r"\([\-0-9]+/[\-0-9]+\)([ATGCRYKMSWBDHVN]+)"
    match1 = re.fullmatch(ref1, site)
    match2 = re.fullmatch(ref2, site)
    match3 = re.fullmatch(ref3, site)
    if match1 is not None:
        seq = match1.group(1)
    elif match2 is not None:
        seq = match2.group(1)
    elif match3 is not None:
        seq = match3.group(1)
    elif "_" in site and "^" in site and site.count("_") == site.count("^"):
        seq = site.replace("_","").replace("^","") 
    else:
        raise ValueError("The sequence does not match the format pattern for representing a cutting site.")  
    return seq 

class _CUTSITES:
    def __setitem__(self, key, item):
        self.__dict__[key] = Cutsite(compilecutsite(item), item, key)

    def __getitem__(self, key):
        return self.__dict__[key]
    
    def items(self):
        return self.__dict__.items()

    def keys(self): 
        return self.__dict__.keys() 
    
    def values(self): 
        return self.__dict__.values() 

class Cutsite:
    def __repr__(self):
        return self.cutsite 

    def __init__(self, seq, site, name):
        self.seq       = Qseq(seq)
        self.rcseq     = Qseq(seq.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB")))[::-1]
        self.cutsite   = Qseq(site) 
        self.name      = name
        self.seq.parent             = self
        self.seq.name               = "seq"
        self.seq.parental_class     = "Cutsite"
        self.rcseq.parent           = self
        self.rcseq.name             = "rcseq"
        self.rcseq.parental_class   = "Cutsite"
        self.cutsite.parent         = self
        self.cutsite.name           = "cutsite"
        self.cutsite.parental_class = "Cutsite"

lib = _CUTSITES()
lib["AatII"]     = "G_ACGT^C"
lib["AbaSI"]     = "CNNNNNNNNN^NN_NNNNNNNNNG"
lib["Acc65I"]    = "G^GTAC_C"
lib["AccI"]      = "GT^MK_AC"
lib["AciI"]      = "CCGC(-3/-1)"
lib["AclI"]      = "AA^CG_TT"
lib["AcuI"]      = "CTGAAG(16/14)"
lib["AfeI"]      = "AGC^_GCT"
lib["AflII"]     = "C^TTAA_G"
lib["AflIII"]    = "A^CRYG_T"
lib["AgeI"]      = "A^CCGG_T"
lib["AhdI"]      = "GACNN_N^NNGTC"
lib["AleI"]      = "CACNN^_NNGTG"
lib["AluI"]      = "AG^_CT"
lib["AlwI"]      = "GGATC(4/5)"
lib["AlwNI"]     = "CAG_NNN^CTG"
lib["Aor51HI"]   = "AGC^_GCT"
lib["ApaI"]      = "G_GGCC^C"
lib["ApaLI"]     = "G^TGCA_C"
lib["ApoI"]      = "R^AATT_Y"
lib["AscI"]      = "GG^CGCG_CC"
lib["AseI"]      = "AT^TA_AT"
lib["AsiSI"]     = "GCG_AT^CGC"
lib["AvaI"]      = "C^YCGR_G"
lib["AvaII"]     = "G^GWC_C"
lib["AvrII"]     = "C^CTAG_G"
lib["BaeGI"]     = "G_KGCM^C"
lib["BaeI"]      = "(10/15)ACNNNNGTAYC(12/7)"
lib["BamHI"]     = "G^GATC_C"
lib["BanI"]      = "G^GYRC_C"
lib["BanII"]     = "G_RGCY^C"
lib["BbsI"]      = "GAAGAC(2/6)"
lib["BbvCI"]     = "CCTCAGC(-5/-2)"
lib["BbvI"]      = "GCAGC(8/12)"
lib["BccI"]      = "CCATC(4/5)"
lib["BceAI"]     = "ACGGC(12/14)"
lib["BcgI"]      = "(10/12)CGANNNNNNTGC(12/10)"
lib["BciVI"]     = "GTATCC(6/5)"
lib["BclI"]      = "T^GATC_A"
lib["BfaI"]      = "C^TA_G"
lib["BglI"]      = "GCCN_NNN^NGGC"
lib["BglII"]     = "A^GATC_T"
lib["BlpI"]      = "GC^TNA_GC"
lib["BmgBI"]     = "CACGTC(-3/-3)"
lib["BmrI"]      = "ACTGGG(5/4)"
lib["BmtI"]      = "G_CTAG^C"
lib["BpmI"]      = "CTGGAG(16/14)"
lib["BpuEI"]     = "CTTGAG(16/14)"
lib["Bpu10I"]    = "CCTNAGC(-5/-2)"
lib["BsaAI"]     = "YAC^_GTR"
lib["BsaBI"]     = "GATNN^_NNATC"
lib["BsaHI"]     = "GR^CG_YC"
lib["BsaI"]      = "GGTCTC(1/5)"
lib["BsaJI"]     = "C^CNNG_G"
lib["BsaWI"]     = "W^CCGG_W"
lib["BsaXI"]     = "(9/12)ACNNNNNCTCC(10/7)"
lib["BseRI"]     = "GAGGAG(10/8)"
lib["BseYI"]     = "CCCAGC(-5/-1)"
lib["BsgI"]      = "GTGCAG(16/14)"
lib["BsiEI"]     = "CG_RY^CG"
lib["BsiHKAI"]   = "G_WGCW^C"
lib["BsiWI"]     = "C^GTAC_G"
lib["BslI"]      = "CCNN^NNN_NNGG"
lib["BsmAI"]     = "GTCTC(1/5)"
lib["BsmBI"]     = "CGTCTC(1/5)"
lib["BsmFI"]     = "GGGAC(10/14)"
lib["BsmI"]      = "GAATGC(1/-1)"
lib["BspCNI"]    = "CTCAG(9/7)"
lib["BspEI"]     = "T^CCGG_A"
lib["BspHI"]     = "T^CATG_A"
lib["Bsp1286I"]  = "G_DGCH^C"
lib["BspMI"]     = "ACCTGC(4/8)"
lib["BsrBI"]     = "CCGCTC(-3/-3)"
lib["BsrDI"]     = "GCAATG(2/0)"
lib["BsrFI"]     = "R^CCGG_Y"
lib["BsrGI"]     = "T^GTAC_A"
lib["BsrI"]      = "ACTGG(1/-1)"
lib["BssHII"]    = "G^CGCG_C"
lib["BssSI"]     = "CACGAG(-5/-1)"
lib["BstAPI"]    = "GCAN^NNN_NTGC"
lib["BstBI"]     = "TT^CG_AA"
lib["BstEII"]    = "G^GTNAC_C"
lib["BstNI"]     = "CC^W_GG"
lib["BstUI"]     = "CG^_CG"
lib["BstXI"]     = "CCAN_NNNN^NTGG"
lib["BstYI"]     = "R^GATC_Y"
lib["BstZ17I"]   = "GTA^_TAC"
lib["Bsu36I"]    = "CC^TNA_GG"
lib["BtgI"]      = "C^CRYG_G"
lib["BtgZI"]     = "GCGATG(10/14)"
lib["BtsCI"]     = "GGATG(2/0)"
lib["BtsIMutI"]  = "CAGTG(2/0)"
lib["BtsI"]      = "GCAGTG(2/0)"
lib["Cac8I"]     = "GCN^_NGC"
lib["ClaI"]      = "AT^CG_AT"
lib["CspCI"]     = "(11/13)CAANNNNNGTGG(12/10)"
lib["CviAII"]    = "C^AT_G"
lib["CviKI-1"]   = "RG^_CY"
lib["CviQI"]     = "G^TA_C"
lib["DdeI"]      = "C^TNA_G"
lib["DpnI"]      = "GA^_TC"
lib["DraI"]      = "TTT^_AAA"
lib["DraIII"]    = "CAC_NNN^GTG"
lib["DrdI"]      = "GACNN_NN^NNGTC"
lib["EaeI"]      = "Y^GGCC_R"
lib["EagI"]      = "C^GGCC_G"
lib["EarI"]      = "CTCTTC(1/4)"
lib["EciI"]      = "GGCGGA(11/9)"
lib["Eco47III"]  = "AGC^_GCT"
lib["Eco53kI"]   = "GAG^_CTC"
lib["EcoNI"]     = "CCTNN^N_NNAGG"
lib["EcoO109I"]  = "RG^GNC_CY"
lib["EcoP15I"]   = "CAGCAG(25/27)"
lib["EcoRI"]     = "G^AATT_C"
lib["EcoRV"]     = "GAT^_ATC"
lib["Esp3I"]     = "CGTCTC(1/5)"
lib["FatI"]      = "^_CATG"
lib["FauI"]      = "CCCGC(4/6)"
lib["Fnu4HI"]    = "GC^N_GC"
lib["FokI"]      = "GGATG(9/13)"
lib["FseI"]      = "GG_CCGG^CC"
lib["FspEI"]     = "CC(12/16)"
lib["FspI"]      = "TGC^_GCA"
lib["HaeII"]     = "R_GCGC^Y"
lib["HaeIII"]    = "GG^_CC"
lib["HgaI"]      = "GACGC(5/10)"
lib["HhaI"]      = "G^CG_C"
lib["HincII"]    = "GTY^_RAC"
lib["HindIII"]   = "A^AGCT_T"
lib["HinfI"]     = "G^ANT_C"
lib["HinP1I"]    = "G^CG_C"
lib["HpaI"]      = "GTT^_AAC"
lib["HphI"]      = "GGTGA(8/7)"
lib["HpyAV"]     = "CCTTC(6/5)"
lib["HpyCH4III"] = "AC^N_GT"
lib["HpyCH4IV"]  = "A^CG_T"
lib["HpyCH4V"]   = "TG^_CA"
lib["Hpy99I"]    = "^CGWCG_"
lib["Hpy188I"]   = "TC^N_GA"
lib["Hpy166II"]  = "GTN^_NAC"
lib["Hpy188III"] = "TC^NN_GA"
lib["I-CeuI"]    = "TAACTATAACGGTCCTAAGGTAGCGAA(-9/-13)"
lib["I-SceI"]    = "TAGGGATAACAGGGTAAT(-9/-13)"
lib["KasI"]      = "G^GCGC_C"
lib["KflI"]      = "GG^GWC_CC"
lib["KpnI"]      = "G_GTAC^C"
lib["LguI"]      = "GCTCTTC(1/4)"
lib["LpnPI"]     = "CCDG(10/14)"
lib["MboI"]      = "^_GATC"
lib["MboII"]     = "GAAGA(8/7)"
lib["MfeI"]      = "C^AATT_G"
lib["MluCI"]     = "^_AATT"
lib["MluI"]      = "A^CGCG_T"
lib["MlyI"]      = "GAGTC(5/5)"
lib["MmeI"]      = "TCCRAC(20/18)"
lib["MnlI"]      = "CCTC(7/6)"
lib["MscI"]      = "TGG^_CCA"
lib["MseI"]      = "T^TA_A"
lib["MslI"]      = "CAYNN^_NNRTG"
lib["MspA1I"]    = "CMG^_CKG"
lib["MspI"]      = "C^CG_G"
lib["MspJI"]     = "CNNR(9/13)"
lib["MwoI"]      = "GCNN^NNN_NNGC"
lib["NaeI"]      = "GCC^_GGC"
lib["NarI"]      = "GG^CG_CC"
lib["NciI"]      = "CC^S_GG"
lib["NcoI"]      = "C^CATG_G"
lib["NdeI"]      = "CA^TA_TG"
lib["NgoMIV"]    = "G^CCGG_C"
lib["NheI"]      = "G^CTAG_C"
lib["NlaIII"]    = "_CATG^"
lib["NlaIV"]     = "GGN^_NCC"
lib["NmeAIII"]   = "GCCGAG(21/19)"
lib["NotI"]      = "GC^GGCC_GC"
lib["NruI"]      = "TCG^_CGA"
lib["NsiI"]      = "A_TGCA^T"
lib["NspI"]      = "R_CATG^Y"
#lib["Nt.AlwI"]   = "GGATC(4/-5)"
#lib["Nt.BbvCI"]  = "CCTCAGC(-5/-7)"
#lib["Nt.BsmAI"]  = "GTCTC(1/-5)"
#lib["Nt.BspQI"]  = "GCTCTTC(1/-7)"
#lib["Nt.BstNBI"] = "GAGTC(4/-5)"
#lib["Nt.CviPII"] = "(0/-1)CCD"
lib["PacI"]      = "TTA_AT^TAA"
lib["PaqCI"]     = "CACCTGC(4/8)"
lib["PciI"]      = "A^CATG_T"
lib["PflMI"]     = "CCAN_NNN^NTGG"
lib["PI-PspI"]   = "TGGCAAACAGCTATTATGGGTATTATGGGT(-13/-17)"
lib["PI-SceI"]   = "ATCTATGTCGGGTGCGGAGAAAGAGGTAAT(-15/-19)"
lib["PleI"]      = "GAGTC(4/5)"
lib["PluTI"]     = "G_GCGC^C"
lib["PmeI"]      = "GTTT^_AAAC"
lib["PmlI"]      = "CAC^_GTG"
lib["PpuMI"]     = "RG^GWC_CY"
lib["PshAI"]     = "GACNN^_NNGTC"
lib["PsiI"]      = "TTA^_TAA"
lib["PspGI"]     = "^_CCWGG"
lib["PspOMI"]    = "G^GGCC_C"
lib["PspXI"]     = "VC^TCGA_GB"
lib["PstI"]      = "C_TGCA^G"
lib["PvuI"]      = "CG_AT^CG"
lib["PvuII"]     = "CAG^_CTG"
lib["RsaI"]      = "GT^_AC"
lib["RsrII"]     = "CG^GWC_CG"
lib["SacI"]      = "G_AGCT^C"
lib["SacII"]     = "CC_GC^GG"
lib["SalI"]      = "G^TCGA_C"
lib["SanDI"]     = "GG^GWC_CC" 
lib["SapI"]      = "GCTCTTC(1/4)"
lib["Sau96I"]    = "G^GNC_C"
lib["SbfI"]      = "CC_TGCA^GG"
lib["ScaI"]      = "AGT^_ACT"
lib["ScrFI"]     = "CC^N_GG"
lib["SexAI"]     = "A^CCWGG_T"
lib["SfaNI"]     = "GCATC(5/9)"
lib["SfcI"]      = "C^TRYA_G"
lib["SfiI"]      = "GGCCN_NNN^NGGCC"
lib["SfoI"]      = "GGC^_GCC"
lib["SgrAI"]     = "CR^CCGG_YG"
lib["SmaI"]      = "CCC^_GGG"
lib["SmlI"]      = "C^TYRA_G"
lib["SnaBI"]     = "TAC^_GTA"
lib["SpeI"]      = "A^CTAG_T"
lib["SphI"]      = "G_CATG^C"
lib["SrfI"]      = "GCCC^_GGGC"
lib["SspI"]      = "AAT^_ATT"
lib["StuI"]      = "AGG^_CCT"
lib["StyD4I"]    = "^_CCNGG"
lib["StyI"]      = "C^CWWG_G"
lib["SwaI"]      = "ATTT^_AAAT"
lib["TaqI"]      = "T^CG_A"
lib["TfiI"]      = "G^AWT_C"
lib["TseI"]      = "G^CWG_C"
lib["Tsp45I"]    = "^_GTSAC"
lib["TspRI"]     = "_NNCASTGNN^"
lib["Tth111I"]   = "GACN^N_NGTC"
lib["XbaI"]      = "T^CTAG_A"
lib["XcmI"]      = "CCANNNN_N^NNNNTGG"
lib["XhoI"]      = "C^TCGA_G"
lib["XmaI"]      = "C^CCGG_G"
lib["XmnI"]      = "GAANN^_NNTTC"
lib["ZraI"]      = "GAC^_GTC"

defaultkeys  = frozenset(lib.__dict__.keys())
new_cutsites = set([]) 
