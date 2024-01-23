from resources import Resources
from typing import Dict
from variant import VariantInfo


class Classify:
    """Runs the rule checks and calculates the final verdict for the variant"""

    def __init__(self):
        self.resources = Resources()
        self.resources.read_datasets()
        self.truncating = ["frameshift", "splice_acceptor", "splice_donor", "stop_gained", "start_lost",
                     "missense_variant&splic"]

    def check_PVS1(self, variant) -> int:
        """
        Certain types of variants (e.g., nonsense, frameshift, canonical
        +- 1 or 2 splice sites, initiation codon, single exon or multiexon
        deletion) in a gene where LOF is a known mechanism of disease
        SO-terms: +-2 bases from splice site = splice acceptor or donor variant. SO-terms has also splice region variants that are wider region
        """

        funcs_tmp = ["frameshift", "splice_acceptor", "splice_donor", "stop_gained", "start_lost",
                     "missense_variant&splice_acceptor", "missense_variant&splice_donor"]  
        funcs_tmp2 = "inframe"
        funcs_tmp3 = "splic"

        PVS = 0
        PVS_t1 = 0
        PVS_t2 = 0
        PVS_t3 = 0

        consequence = variant.csq['Consequence']
        for fc in funcs_tmp:
            if consequence.find(fc) >= 0 and consequence.find(funcs_tmp2) < 0:
                PVS_t1 = 1
                break
        # wait to check LOF genes use the LoFtool_percentile,but  how to know is the disease mechanism
        try:
            if self.resources.lof_genes[variant.csq['SYMBOL']] == '1':
                PVS_t2 = 1
        except KeyError:
            PVS_t2 = 0
        else:
            pass
        # print("PVSt1= %d PVSt2= %d" % (PVS_t1,PVS_t2) )

        # begin check the site is really affect the splicing

        try:
            if float(variant.csq['rf_score']) > 0.515 or float(variant.csq['ada_score']) > 0.708:
                PVS_t3 = 1
        except ValueError:
            pass
        else:
            pass

        if PVS_t1 != 0 and PVS_t2 != 0:
            PVS = 1
            if consequence.find(funcs_tmp3) >= 0 and PVS_t3 != 1:
                PVS = 0

        # Checking first/last exon from the variant.csq['EXON'] !

        """
        #begin check it in the AAChange.knownGene for the major/Canonical isoform, not 1/last exon
        #SUFU:uc001kvy.2:exon6:c.G716A:p.R239Q
        line_tmp2=variant.csq['AAChange.knownGene']
        #for cls0 in line_tmp2.split(','):
        for cls0 in re.split("[,;]",line_tmp2):
            cls0_1=cls0.split(':')
            if len(cls0_1)>1:
                trans_id=cls0_1[1]
                exon=cls0_1[2]
                try:
                    exon_lth="exon"+knownGeneCanonical_dict[trans_id]
                    #if exon==exon_lth or exon =="exon1": # not 1 or last exon
                    if exon==exon_lth: # relax for only last exon
                        PVS=0

                    try:
                        if (float(knownGeneCanonical_ed_dict[trans_id])-float( cls[Allels_flgs['Start']]  ))<50: # means close  3' of gene 50 bp.
                            PVS=0
                    except ValueError:
                        pass
                    else:
                        pass


                except KeyError:
                    pass
                else:
                    pass
         """

        return (PVS)

    def check_PS1(self,variant) -> int:
        """
        PS1 Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
        Example: Val->Leu caused by either G>C or G>T in the same codon
        variant.csq['Amino_acids']  : 'T/M'
        intervar:NOD2:NM_001293557:exon3:c.C2023T:p.R675W,NOD2:NM_022162:exon4:c.C2104T:p.R702W
        """

        PS1 = 0
        PS1_t1 = 0
        PS1_t2 = 0
        PS1_t3 = 0

        funcs_tmp = ["missense"]
        ACGTs = ["A", "C", "G", "T"]
        consequence = variant.csq['Consequence']
        for fc in funcs_tmp:
            if consequence.find(fc) >= 0:
                PS1_t1 = 1
                # need to wait to check Same amino acid change as a previously pathogenic variant

                # change info in: 'Amino_acids': 'T/M'
                aa_change = variant.csq['Amino_acids']

                aa = aa_change.split('/')[1]
                k = variant.chrom
                k = k.split('chr')[-1]
                keys_tmp2 = k + "_" + str(variant.start) + "_" + str(variant.end) + "_" + variant.alt
                try:
                    if self.resources.aa_changes[keys_tmp2] == aa:
                        PS1_t2 = 1
                except KeyError:
                    for nt in ACGTs:
                        if nt != variant.alt and nt != variant.ref:
                            keys_tmp3 = variant.chrom + "_" + str(variant.start) + "_" + str(variant.end) + "_" + nt
                            try:
                                if self.resources.aa_changes[keys_tmp3] == aa:
                                    PS1_t2 = 1
                            except KeyError:
                                pass
                            else:
                                pass

                else:
                    pass

        # Check the splicing score thing? Do I have the relevant info in the VEP annotated file already?

        try:
            if float(variant.csq['rf_score']) > 0.515 or float(
                    variant.csq['ada_score']) > 0.708:  # means alter the splicing
                PS1_t3 = 1
            if variant.csq['rf_score'] == "" or variant.csq['ada_score'] == "":  # absent also means not in splicing
                PS1_t3 = 0
        except ValueError:
            pass
        else:
            pass

        if PS1_t1 != 0 and PS1_t2 != 0:
            PS1 = 1
            if PS1_t3 == 1:  # remove the splicing affect
                PS1 = 0
        return (PS1)

    def check_PM1(self, variant) -> int:
        """
        Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of
        an enzyme) without benign variation
        """
        PM1 = 0
        PM1_t1 = 0
        PM1_t2 = 0

        funcs_tmp = ["missense"]
        consequence = variant.csq['Consequence']
        for fc in funcs_tmp:
            if consequence.find(fc) >= 0:
                PM1_t1 = 1;
            # need to wait to check whether in hot spot  or  functional domain/without benign variation

            ### Check the interpro domain thing, do we need new custom annotation to the VEP annotation?
            if variant.csq['Interpro_domain'] != '':
                k = variant.chrom
                k = k.split('chr')[-1]
                keys_tmp2 = k + "_" + variant.csq['SYMBOL']
                try:
                    found = self.resources.PM1_domains[keys_tmp2]
                    f = found.split("|")
                    for k in f:
                        if k in variant.csq['Interpro_domain'].replace("_", " "):
                            PM1_t2 = 0
                except KeyError:
                    PM1_t2 = 1
                else:
                    pass

        if PM1_t1 == 1 and PM1_t2 == 1:
            PM1 = 1

        return (PM1)

    def check_PM2(self, variant) -> int:
        """
        Absent from controls (or at extremely low frequency if recessive) (Table 6) in Exome Sequencing Project,
        1000 Genomes Project, or Exome Aggregation Consortium


        Varsome: if AF not found in GnomAD with VALID coverage or
        For AD+X+AD/AR genes AC<5
        For AR if homozygous AC < 3 OR AF <0.0001
        """

        PM2 = 0
        PM2_s = 0
        
        # Check if absent
        freqs = ['gnomADexomes_AF', 'ExAC_AF', '1000Gp3_AF']

        not_absent = False
        for key in freqs:
            if (variant.csq[key] != ''):  # absent in all 3 controls

                try:
                    if float(variant.csq[key]) > 0:
                        not_absent = True
                except ValueError:
                    pass

        # if not absent, check the ACs
        if not_absent:
            try:
                CGDinfo = self.resources.CGD[variant.csq['SYMBOL']]
                #print(CGDinfo)
            except KeyError:
                pass
            else:
                if  CGDinfo['INHERITANCE'] in ['AR', 'AD/AR']:
                    # AR: homozygous no. < 3
                    try:
                        if int(variant.csq['gnomADexomes_nhomalt']) < 3 and (
                                float(variant.csq['gnomADexomes_AF']) < 0.01):
                            PM2 = 1
                    except ValueError:
                        # nhomalt is absent
                        PM2 = 1

                elif CGDinfo['INHERITANCE'] in ['AD', 'XL']:

                    # AD: (XL, AD/AR ) AC < 5
                    try:
                        if int(variant.csq['gnomADexomes_AF']) < 0.0001: 
                            # 0.001 to 0.0001
                            PM2 = 1
                    except ValueError:
                        # AC value is absent ''
                        PM2 = 1

        else:
            if "splic" in variant.csq['Consequence']:
                try:
                    if float(variant.csq['gnomADgenomes_AF']) > 0.001:
                        PM2 = 0
                    else:
                        PM2 = 1
                except ValueError:
                    PM2 = 1
            else:
                PM2 = 1

        if PM2 == 1:
            try:
                if float(variant.csq['phyloP100way_vertebrate']) > 7.2:
                    PM2_s = 1
                    PM2 = 0
            except ValueError:
                pass

        return (PM2, PM2_s)

    def check_PM4(self, variant) -> int:
        '''
        Protein length changes as a result of in-frame deletions/insertions in a nonrepeat region or stop-loss variants

        Varsome: don't count if PVS1 was triggered?
        '''
        PM4 = 0
        PM4_t1 = 0
        PM4_t2 = 0

        # Are these the right VEP gene change terms?, checked ok
        funcs_tmp = ["inframe_insertion", "inframe_deletion", "stop_lost"]

        for fc in funcs_tmp:
            if variant.csq['Consequence'].find(fc) >= 0:
                PM4_t1 = 1

            # need to wait to check  in a nonrepeat region
        if variant.csq['rmsk'] == '':
            PM4_t2 = 1

        if variant.csq['rmsk'] != '' and variant.csq['Consequence'].find("stop_lost") >= 0:
            PM4_t2 = 1

        if PM4_t1 != 0 and PM4_t2 != 0:
            PM4 = 1

        return (PM4)

    def check_PM5(self, variant) -> int:
        """
        Novel missense change at an amino acid residue where a different missense change determined to be
        pathogenic has been seen before;Example: Arg156His is pathogenic; now you observe Arg156Cys
        NOD2:NM_001293557:exon3:c.C2023T:p.R675W,NOD2:NM_022162:exon4:c.C2104T:p.R702W
        """
        PM5 = 0
        PM5_t1 = 0
        PM5_t2 = 0
        PM5_t3 = 0

        funcs_tmp = ["missense"]  # removed 'non-synonymous' because SO-terms used by VEP doesn't include it
        ACGTs = ["A", "C", "G", "T"]

        for fc in funcs_tmp:
            if variant.csq['Consequence'].find(fc) >= 0:
                PM5_t1 = 1

                # need to wait to check no-Same amino acid change as a previously pathogenic variant
                aa_change = variant.csq['Amino_acids']
                aa = aa_change.split('/')[1]
                k = variant.chrom
                k = k.split('chr')[-1]
                keys_tmp2 = k + "_" + str(variant.start) + "_" + str(variant.end) + "_" + variant.alt
                try:
                    if self.resources.aa_changes[keys_tmp2]:
                        PM5_t2 = 0

                except KeyError:
                    PM5_t2 = 1
                    PM5_t3 = 0
                    for nt in ACGTs:
                        if nt != variant.alt and nt != variant.ref:
                            keys_tmp3 = k + "_" + str(variant.start) + "_" + str(variant.end) + "_" + nt
                            try:
                                if self.resources.aa_changes[keys_tmp3]:
                                    PM5_t3 = 1
                                if self.resources.aa_changes[keys_tmp3] == aa:
                                    PM5_t2 = 0 * PM5_t2
                            except KeyError:
                                pass
                            else:
                                pass

                else:
                    pass

        if PM5_t1 != 0 and PM5_t2 != 0 and PM5_t3 != 0:
            PM5 = 1
        return (PM5)

    def check_PP2(self, variant) -> int:
        """
        Missense variant in a gene that has a low rate of benign missense variation and in which
        missense variants are a common mechanism of disease
        """
        PP2 = 0

        funcs_tmp = ["missense"]

        for fc in funcs_tmp:
            if variant.csq['Consequence'].find(fc) >= 0:
                # need to check whether gene has a low rate of benign missense variation.....
                try:
                    if self.resources.PP2_genes[variant.csq['SYMBOL']] == '1':
                        PP2 = 1
                except KeyError:
                    PP2 = 0
                else:
                    pass

        return (PP2)

    def check_PP3(self, variant) -> int:
        """
        Multiple lines of computational evidence support a deleterious effect on the gene or gene product
        (conservation, evolutionary, splicing impact, etc.)
        sfit for conservation, GERP++_RS for evolutionary, splicing impact from dbNSFP
        """

        PP3 = 0

        # cutoff_conserv=2 # for GERP++_RS  !!!! Why this is so low???, even 3.597 was lower than for conservation

        # Try finding all the values, search the score values separately and evaluate to a threshold
        P = 0  # add +1 with every pathogenic or benign/non-path prediction
        B = 0

        preds = [
            'BayesDel_noAF_pred', 'DEOGEN2_pred',
            'SIFT_pred', 'Polyphen2_HDIV_pred',
            'PROVEAN_pred', 'M-CAP_pred', 'MetaSVM_pred',
            'MutationTaster_pred', 'MutationAssessor_pred',  # mutationAssessor: H/M
            'MetaRNN_pred', 'LIST-S2_pred', 'PrimateAI_pred'
        ]
        th_preds = [
            'REVEL', 'CADD_PHRED', 'FATHMM_MKL_C',
            'MVP_score', 'DANN_score'
        ]
        # what's the threshold for 'Eigen-phred_coding'?????
        ths = [0.5, 15, 0.5, 0.75, 0.93]
        # 'GERP++_RS'

        # Go through the D/T predictions
        for p in preds:
            try:
                if any(ex in variant.csq[p] for ex in ('D', 'H')):
                    #variant.csq[p] in ('D', 'H'):
                    P += 1
                elif any(ex in variant.csq[p] for ex in ('T', 'M', 'L', 'N')):
                    #variant.csq[p] in ('T', 'M', 'L', 'N'):
                    B += 1
            except ValueError:
                pass

        for p, s in zip(th_preds, ths):
            try:
                if float(variant.csq[p]) > s:
                    P += 1
                elif float(variant.csq[p]) <= s:
                    B += 1
            except ValueError:
                pass

        # GERP as a fall-back in the absence of any other predictions

        if P == 0 and B == 0:
            try:
                #Changed GERP++3.597 to PhyloP
                if float(variant.csq['phyloP100way_vertebrate']) > 3.81:
                    PP3 = 1
            except ValueError:
                pass
        # Prediction ratio for pathogenic prediction
        else:
            ratio = P / (P + B)

            if ratio > 0.53:
                PP3 = 1

        ###########################################
        # is this splice check needed??
        """
        try:
            if float(variant.csq['rf_score'])>0.515 or float(variant.csq['ada_score'])>0.708:
                PP3_t3=1
        except ValueError:
            pass
        else:
            pass

        if (PP3+PP3_t3)>=1:
            PP3=1
        """

        return (PP3)

    def check_PP5(self, variant) -> int:
        """
        Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory
        to perform an independent evaluation
        """
        PP5 = 0

        line_tmp2 = variant.csq['ClinVar_CLNSIG']
        if line_tmp2 != '':

            if line_tmp2.find("ikely pathogenic") >= 0 or line_tmp2.find("athogenic") >= 0:
                if line_tmp2.find("onflicting") == -1:
                    PP5 = 1
        return (PP5)

    def check_BA1(self, variant) -> int:
        """
        BA1 Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
        """
        BA1 = 0

        # Freqs_3pops={'1000g2015aug_all':0,'esp6500siv2_all':0,'gnomAD_genome_ALL':0}
        freqs = ['gnomADexomes_AF', 'ExAC_AF', '1000Gp3_AF', 'gnomADgenomes_AF']
        for key in freqs:
            try:
                if float(variant.csq[key]) > 0.05: BA1 = 1

            except ValueError:
                pass
            else:
                pass

        return (BA1)

    def check_BS1(self, variant) -> int:
        """
        Allele frequency is greater than expected for disorder (see Table 6)
        > 1% in ESP6500all ExAc? need to check more
        """
        BS1 = 0
        cutoff = 0.005 

        freqs = ['gnomADexomes_AF', 'ExAC_AF', '1000Gp3_AF', 'gnomADgenomes_AF']
        for key in freqs:
            try:
                if variant.csq[key] != '':
                    if float(variant.csq[key]) >= cutoff: BS1 = 1
            except ValueError:
                pass
            else:
                pass
        if "splic" in variant.csq['Consequence']:
            try:
                if variant.csq['gnomADgenomes_AF'] != '':
                    if float(variant.csq['gnomADgenomes_AF']) >= cutoff: BS1 = 1
            except ValueError:
                pass

        return (BS1)

    def check_BS2(self, variant) -> int:
        """
        Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked
        (hemizygous) disorder, with full penetrance expected at an early age
        check gnomAD_genome_ALL
        """
        BS2 = 0

        try:
            CGDinfo = self.resources.CGD[variant.csq['SYMBOL']]
        except KeyError:
            # Gene not in CGD / no onset info
            BS2 = 0
        else:
            if 'Adult' in CGDinfo['ONSET'] or 'N/A' in CGDinfo['ONSET']:
                BS2 = 0
            else:
                try:
                    AC = int(variant.csq['gnomADexomes_AC'])
                except ValueError:
                    #print("Not found in gnomad")
                    BS2 = 0
                else:

                    if CGDinfo['INHERITANCE'] in ['AR', 'XL']:
                        # AR or XL: AC > 3
                        if float(variant.csq['gnomADexomes_AF']) > 0.01: 
                            #old:AC > 3:
                            BS2 = 1

                    elif CGDinfo['INHERITANCE'] in ['AD']:
                        # AD: AC > 5
                        if float(variant.csq['gnomADexomes_AF']) > 0.0001:
                            #old: AC > 5: 0.001 to 0.0001
                            BS2 = 1
                    else:
                        # no inheritance info or not in AR, XL or AD
                        BS2 = 0

        return (BS2)

    def check_BP1(self, variant) -> int:
        """
        Missense variant in a gene for which primarily truncating variants are known to cause disease
        truncating:  stop_gain / frameshift deletion/  nonframshift deletion
        We defined Protein truncating variants  (4) (table S1) as single-nucleotide variants (SNVs) predicted to introduce a premature stop codon or to disrupt a splice site, small insertions or deletions (indels) predicted to disrupt a transcript reading frame, and larger deletions
        """
        BP1 = 0

        funcs_tmp = ["missense"]
        funcs_tmp2 = "splic"
        line_tmp = variant.csq['Consequence']
        for fc in funcs_tmp:
            if line_tmp.find(fc) >= 0 and line_tmp.find(funcs_tmp2) < 0:
                # need to wait to check whether truncating is the only cause disease
                try:
                    if self.resources.BP1_genes[variant.csq['SYMBOL']] == '1':
                        BP1 = 1
                except KeyError:
                    BP1 = 0
                else:
                    pass
        return (BP1)

    def check_BP3(self, variant) -> int:
        """
        In-frame deletions/insertions in a repetitive region without a known function
        if the repetitive region is in the domain, this BP3 should not be applied.
        """
        BP3 = 0
        BP3_t1 = 0
        BP3_t2 = 0

        funcs_tmp = ["inframe_deletion", "inframe_insertion"]
        line_tmp = variant.csq['Consequence']
        for fc in funcs_tmp:
            if line_tmp.find(fc) >= 0:
                BP3_t1 = 1;
            # need to wait to check  in a repeat region
        if variant.csq['rmsk'] != '' and variant.csq['Interpro_domain'] == '':  # repeat and not in domain
            BP3_t2 = 1

        if BP3_t1 != 0 and BP3_t2 != 0:
            BP3 = 1
        return (BP3)

    def check_BP4(self, variant) -> int:
        """
        Multiple lines of computational evidence suggest no impact on gene or gene product (conservation,
        evolutionary,splicing impact, etc.)
        """
        BP4 = 0
        BP4_t1 = 0
        BP4_t2 = 0

        # Try finding all the values, search the score values separately and evaluate to a threshold
        P = 0  # add +1 with every pathogenic or benign/non-path prediction
        B = 0

        preds = [
            'BayesDel_noAF_pred', 'DEOGEN2_pred',
            'SIFT_pred', 'Polyphen2_HDIV_pred',
            'PROVEAN_pred', 'M-CAP_pred', 'MetaSVM_pred',
            'MutationTaster_pred', 'MutationAssessor_pred',  # mutationAssessor: H/M
            'MetaRNN_pred', 'LIST-S2_pred', 'PrimateAI_pred'
        ]
        th_preds = [
            'REVEL', 'CADD_PHRED', 'FATHMM_MKL_C',
            'MVP_score', 'DANN_score'
        ]

        ths = [0.5, 15, 0.5, 0.75, 0.93]


        # Go through the D/T predictions
        for p in preds:
            try:
                if any(ex in variant.csq[p] for ex in ('D', 'H')):
                    #variant.csq[p] in ('D', 'H'):
                    P += 1
                elif any(ex in variant.csq[p] for ex in ('T', 'M', 'L', 'N')):
                    #variant.csq[p] in ('T', 'M', 'L', 'N'):
                    B += 1
            except ValueError:
                pass

        for p, s in zip(th_preds, ths):
            try:
                if float(variant.csq[p]) > s:
                    P += 1
                elif float(variant.csq[p]) <= s:
                    B += 1
            except ValueError:
                pass

        # GERP as a fall-back in the absence of any other predictions

        if P == 0 and B == 0:
            try:
                #Changed GERP++ < 3.561to 'phyloP100way_vertebrate'
                #and non-truncating!
                if float(variant.csq['phyloP100way_vertebrate']) < 1.4 and variant.csq['Consequence'] not in self.truncating:
                    BP4_t1 = 1
            except ValueError:
                pass
        # Prediction ratio for benign prediction
        else:
            ratio = B / (P + B)

            if ratio > 0.53:
                BP4_t1 = 1

        # remove splicing effect
        try:
            if float(variant.csq['rf_score']) <= 0.515 and float(variant.csq['ada_score']) <= 0.708:
                BP4_t2 = 1
        except ValueError:
            BP4_t2 = 1  # means absent, this site is not in splicing consensus regions
        else:
            pass

        if (BP4_t1 + BP4_t2) == 2:
            BP4 = 1

        return (BP4)

    def check_BP6(self, variant) -> int:
        """
        Reputable source recently reports variant as benign, but the evidence is not available to the
        laboratory to perform an independent evaluation; Check the ClinVar column to see whether this
        is "benign".
        """
        BP6 = 0

        line_tmp2 = variant.csq['ClinVar_CLNSIG']
        if line_tmp2 != '':
            cls3 = line_tmp2.split(';')
            clinvar_bp = cls3[0]
            if clinvar_bp.find("ikely benign") >= 0 or clinvar_bp.find("enign") >= 0:
                BP6 = 1

        return (BP6)

    def check_BP7(self, variant) -> int:
        """
        A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the
        splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly
        conserved
        """
        BP7 = 0
        BP7_t1 = 0
        BP7_t2 = 0
        cutoff_conserv = 2  # for GERP++

        funcs_tmp = ["synon"]

        line_tmp = variant.csq['Consequence']
        for fc in funcs_tmp:
            if line_tmp.find(fc) >= 0:

                # need to wait to check the  impact to the splice from dbscSNV
                # either score(ada and rf) >0.6 as splicealtering
                if variant.csq['rf_score'] == "" or variant.csq['ada_score'] == "":
                    BP7_t1 = 1  # absent means it is not in the  splice consensus sequence
                else:
                    if float(variant.csq['rf_score']) < 0.515 and float(variant.csq['ada_score']) < 0.708:
                        BP7_t1 = 1
        # check the conservation score of gerp++ > 2
        try:
            if float(variant.csq['GERP++_RS']) <= float(cutoff_conserv) or variant.csq['GERP++_RS'] == '':
                BP7_t2 = 1
        except ValueError:
            # absent means there are gaps in the multiple alignment,so cannot have the score,not conserved
            BP7_t2 = 1
        else:
            pass

        if BP7_t1 != 0 and BP7_t2 != 0:
            BP7 = 1
        return (BP7)

    def verdict(self, scores: Dict[str, int]) -> str:
        """Do the verdict by checking the sums of each rule sets"""
        verdicts = ["Pathogenic", "Likely pathogenic", "Benign", "Likely benign", "Uncertain significance"]
        PAS_out = -1
        BES_out = -1
        BPS_out = 4

        # Likely pathogenic
        # 2. 1 PS1-PS4 AND 1-2 PM1-PM6 OR
        if scores['PS'] == 1:
            if scores['PM'] == 1 or scores['PM'] == 2: PAS_out = 1
            # 1. 1 PVS1 AND 1 PM1-PM6 OR
        if scores['PVS1'] == 1:
            if scores['PM'] >= 1: PAS_out = 1
        # 3. 1 PS1-PS4 AND 2>= PP1-PP5 OR
        if scores['PS'] == 1 and scores['PP'] >= 2: PAS_out = 1
        # 4. 3>= PM1-PM6 OR
        if scores['PM'] >= 3: PAS_out = 1
        # 5. 2 PM1-PM6 AND 2>= PP1-PP5
        if scores['PM'] == 2 and scores['PP'] >= 2: PAS_out = 1
        # 6.     1 PM1-PM6 AND 3>= PP1-PP5, huom varsome like deduction!
        if scores['PM'] == 1 and scores['PP'] >= 3: PAS_out = 1

        # Pathogenic
        # 1. PVS1 AND
        if scores['PVS1'] == 1:
            # a. 1>= (PS1-PS4) OR
            if scores['PS'] >= 1: PAS_out = 0
            # b. 2>= PM1-PM6 OR
            if scores['PM'] >= 2: PAS_out = 0
            # c. 1 PM1-PM6 and 1 PP1-PP5 OR
            if scores['PM'] == 1 and scores['PP'] == 1: PAS_out = 0
            # d. 2>= PP1-PP5
            if scores['PP'] >= 2: PAS_out = 0
        # 2. 2>= PS1-PS4
        if scores['PS'] >= 2: PAS_out = 0
        # 3. 1 PS1-PS4 AND
        if scores['PS'] == 1:
            # a. 3>= PM1-PM6 OR
            if scores['PM'] >= 3: PAS_out = 0
            # b. 2 PM1-PM6 AND 2>= PP1-PP5 OR
            if scores['PM'] == 2 and scores['PP'] >= 2: PAS_out = 0
            # c. 1 PM1-PM6 AND 4>= PP1-PP5
            if scores['PM'] == 1 and scores['PP'] >= 4: PAS_out = 0

        # Likely benign    !!!!CHANGE?? varsome: any strong benign rule (BS1,BS2) is sufficient to trigger a likely benign verdict
        # https://pubmed.ncbi.nlm.nih.gov/29300386/
        # 1. OLD: 1 BS1-BS4 AND 1 BP1-BP7 OR, NEW: 1 BS1/2
        # if scores['BS'] ==1 and scores['BP']==1: BES_out = 3
        if scores['BS'] == 1: BES_out = 3
        # 2. 2>= BP1-BP7
        if scores['BP'] >= 2: BES_out = 3
        # benign
        # 1. BA1 OR 2. 2>= BS1-BS4
        if scores['BA1'] == 1 or scores['BS'] >= 2: BES_out = 2

        # VUS and final check
        if PAS_out != -1 and BES_out == -1: BPS_out = PAS_out
        if PAS_out == -1 and BES_out != -1: BPS_out = BES_out
        if PAS_out == -1 and BES_out == -1: BPS_out = 4
        if PAS_out != -1 and BES_out != -1: BPS_out = 4

        return (verdicts[BPS_out])

    def triggered(self, *args) -> str:
        rules = ['PVS1', 'PS1', 'PM2_s', 'PM1', 'PM2', 'PM4', 'PM5', 'PP2', 'PP3', 'PP5', 'BA1', 'BS1', 'BS2', 'BP1', 'BP3',
                 'BP4', 'BP6', 'BP7']
        l = []
        for a in args:
            l += a

        trig = ','.join([x for x, y in zip(rules, l) if y == 1])

        return (trig)

    def predict(self, variant: VariantInfo) -> str:
        #typehint: : VariantInfo
        """
        Predictor calls rule check methods and collects the scores. Then calls method that calculated the final verdict.
        """
        #Checking the strengthened PM2
        PM2, PM2_s = self.check_PM2(variant)
        if PM2 == 1 or PM2_s == 1:
            #skip BS2 check when PM2 is triggered
            BS2 = 0
        else:
            BS2 = self.check_BS2(variant)
        
        # Call Rules
        PVS1 = [self.check_PVS1(variant)]
        PS = [self.check_PS1(variant), PM2_s]
        PM = [self.check_PM1(variant), PM2, self.check_PM4(variant), self.check_PM5(variant)]
        PP = [self.check_PP2(variant), self.check_PP3(variant), self.check_PP5(variant)]
        BA1 = [self.check_BA1(variant)]
        BS = [self.check_BS1(variant), BS2]
        BP = [self.check_BP1(variant), self.check_BP3(variant), self.check_BP4(variant), self.check_BP6(variant),self. check_BP7(variant)]

        # PS, PM, PP, BS, BP, PVS1, BA1,
        # save the names and the sum of the rules PS1-4, PM1-6 etc
        scores = {'PVS1': sum(PVS1), 'PS': sum(PS), 'PM': sum(PM), 'PP': sum(PP), 'BA1': sum(BA1), 'BS': sum(BS),
                  'BP': sum(BP)}

        # Final verdict
        result = self.verdict(scores)

        # return the verdict AND the triggered rules!
        rules = self.triggered(PVS1, PS, PM, PP, BA1, BS, BP)

        return ((result, rules))

    
