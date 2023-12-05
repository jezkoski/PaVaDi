import gzip
import sys

class Resources():
    """Reads and keeps track of the resources used for the annotation"""

    def __init__(self):
        self.lof_genes = {}
        self.aa_changes = {}
        self.PP2_genes = {}
        self.BP1_genes = {}
        self.PS4_snps = {}
        self.PM1_domains = {}
        self.CGD = {}

    def read_datasets(self):
        #LOF
        try:
            with open('resources/PVS1.LOF.genes.hg19') as f:
                for line in f:
                    self.lof_genes[line.strip()] = '1'
        except IOError:
            print("Error: no LOFgenes file")
            sys.exit()

        #AA changes,
        try:
            with open('resources/PS1.AA.change.patho.hg19') as f:
                for line in f:
                    values=line.split('\t')
                    key=values[0]+"_"+values[1]+"_"+values[2]+"_"+values[4]
                    self.aa_changes[key]=values[6]

        except IOError:
            print("Error: no AAchange file")
            sys.exit()

        #Domain benign
    
        try:
            with open('resources/PM1_domains_with_benigns.hg19') as f:
                for line in f:
                    l=line.split('\t')
                    key=l[0]+"_"+l[1]
                    self.PM1_domains[key]=l[2]

        except IOError:
            print('Error: no domains_benign file')
            sys.exit()

        #PP2 genes
        try:
            with open('resources/PP2.genes.hg19') as f:
                for line in f:
                    l=line.strip()
                    if len(l)>1:
                        key=l
                        self.PP2_genes[key]='1'
        except IOError:
            print("Error: no PP2 genes file")
            sys.exit()
 
        #PS4 snps
        try:
            with open('resources/PS4.variants.hg19') as f:
                for line in f:
                    l=line.split('\t')
                    if len(l[0])>=1:
                        key=l[0]+"_"+l[1]+"_"+l[1]+"_"+l[3]+"_"+l[4]
                        self.PS4_snps[key]='1'
        except IOError:
            print("Error: no PS4 variants file")
            sys.exit()      

        #BP1 genes 
        try:
            with open('resources/BP1.genes.hg19') as f:
                for line in f:
                    l=line.strip()
                    if len(l)>1:
                        self.BP1_genes[l]='1'
        except IOError:
            print("Error: no BP1 genes file")
            sys.exit()

        #CGD
        try:
            with open('resources/CGD022022.txt') as f:
                for line in f:
                    l=line.strip().split('\t')
                    if len(l[0])>=1:
                        # Use gene name or HGNC id as a key?? VEP CSQ: HGNC_ID
                        key=l[0]
                        self.CGD[key] = {'HGNC ID': l[1], 'INHERITANCE': l[2], 'ONSET': l[3]}
        except IOError:
            print("Error: no CGD file")
            sys.exit()      


    def flip_ACGT(self, acgt):
        nt="";
        if acgt=="A":
            nt="T"
        if acgt=="T":
            nt="A"
        if acgt=="C":
            nt="G"
        if acgt=="G":
            nt="C"
        if acgt=="N":
            nt="N"
        if acgt=="X":
            nt="X"
        return(nt)


