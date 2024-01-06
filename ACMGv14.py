#!/usr/bin/env python
import optparse, sys
from classify import Classify
from variant import VariantInfo
from datetime import date
from contextlib import closing
from cyvcf2 import VCF,Writer

usage = """Usage: %prog [OPTION] -i  INPUT -o  OUTPUT ...
       %prog  --config=config.ini ...
"""

version = """ACMG annotator, based on InterVar, version 1.4.
"""

description = """=============================================================================
                                                                       
Interpretation of Pathogenic/Benign for variants 
"""



def Annotator(File,outfile):

    #Does the rule checks and reads in resources
    classify = Classify()
    today = date.today().strftime("%d/%m/%Y")
    
    print(f"Starting ACMG annotation with {version} on {today}")
    verdict_counts = {'Pathogenic': 0, 'Likely pathogenic': 0, 'Uncertain significance': 0, 'Likely benign': 0, 'Benign': 0}

    with closing(VCF(str(File))) as vcf:
        #ACMG info tag and the new outfile
        vcf.add_info_to_header({'ID': 'ACMG verdict', 'Description': 'Variant prediction according to ACMG rules', 'Type':'Character','Number': '1'})
        vcf.add_info_to_header({'ID': 'Triggered ACMG rules', 'Description': 'Triggered ACMG rules', 'Type':'Character','Number': '1'})
        # Add version info etc to header
        info = f"ACMGannotator version {version}, annotation date: {today}"
        vcf.add_info_to_header({'ID': 'Version', 'Description': f"{version}", 'Type': 'Character','Number': '1'})
        vcf.add_info_to_header({'ID': 'Annotation date', 'Description': f"File was annotated {today}", 'Type':'Character','Number':'1'})
        
        fname=outfile
        w=Writer(fname,vcf)

        #setup CSQ fields
        raw = vcf.raw_header.splitlines()
        csq_fields = VariantInfo.get_vep_csq_fields(raw)

        for v in vcf:

            #create new variant object with only INFO to send off
            variant=VariantInfo.get_variantinfo(v)
            variant.parse_csq(csq_fields)

            #send new variant object to predictor/assigner and get back the verdict and triggered rules
            result,rules = classify.predict(variant)
            verdict_counts[result] += 1

            #Write ACMG verdict to the INFO tag
            v.INFO["ACMG verdict"] = result
            v.INFO["Triggered ACMG rules"] = rules
            w.write_record(v)

    #Print out some overall stats from the annotation
    print(f"Total number of variants annotated: {sum(verdict_counts.values())}")
    print(verdict_counts)
    

def main():
        
    parser = optparse.OptionParser(usage=usage, version=version, description=description)
    
    parser.add_option("-?", action="help", help=optparse.SUPPRESS_HELP, dest="help")
    parser.add_option("-v", action="version", help=optparse.SUPPRESS_HELP, dest="version")
    
    #parser.add_option("-c", "--config", dest="config", action="store",
    #              help="The config file of all options. it is for your own configure file.You can edit all the options in the configure and if you use this options,you can ignore all the other options bellow", metavar="config.ini")
    
    parser.add_option("-i", "--input", dest="input", action="store",
                  help="The input file contains your variants", metavar="example/ex1.vcf.gz")

    parser.add_option("-o", "--output", dest="output", action="store",
                  help="The prefix of output file which contains the results, the file of results will be as [$$prefix].vcf.gz ", metavar="example/myanno")
    
    (options,args) = parser.parse_args()
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
        
    print("%s" %description)
    print("%s" %version)        

    if options.input != None:
        File = options.input

    if options.output != None:
        outfile = options.output

    Annotator(File,outfile)




if __name__ =="__main__":
        main()

