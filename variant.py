from contextlib import closing
from cyvcf2 import VCF, Writer
from cyvcf2 import Variant as cyVCF2variant
from typing import Any, Type, TypeVar, Optional, Dict, List
import re
import attr


V = TypeVar("V", bound="VariantInfo")

@attr.s(auto_attribs=True, repr=False, eq=False, order=False, slots=True)
class VariantInfo:
    """
    Creates Variant object with only INFO fields
    """
    
        
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    info: Dict[str, Any] = attr.Factory(dict)
    csq: Optional[Dict[str,Any]] = attr.Factory(dict)
        
    def __attrs_post_init__(self):
        if self.ref == ".":
            raise ValueError(
                "ref_allele cannot be missing ('.'). Try normalize the variant first."
            )
        if self.alt is None or self.alt == ".":
            raise ValueError(
                "alt_allele cannot be missing ('.' or None). Try normalize the variant first."
            )
            
    @classmethod
    def get_variantinfo(cls: Type[V], variant: cyVCF2variant) -> V:
        """
        Create VariantInfo object based on cyvcf2.Variant
        """
        
        return cls(chrom=variant.CHROM,
                  start=variant.start+1,
                  end=variant.end,
                  ref=variant.REF,
                  alt=variant.ALT[0],
                  info=dict(variant.INFO)
                  )
    
    def parse_csq(self, csq_fields: List[str]):
        """Parse the CSQ fields from variant"""
        
        csqvalues = self.info['CSQ'].split("|")
        self.csq = dict(zip(csq_fields, csqvalues))
        self.info["CSQ"] = self.csq

    @classmethod
    def get_vep_csq_fields(cls: Type[V], vcf_raw_headers: List[str]) -> List[str]:
        """"Extract CSQ fields from VCF header"""
    
        try:
            csq_info_header = next(l for l in reversed(vcf_raw_headers) if l.startswith("##INFO=<ID=CSQ,"))
        except StopIteration:
            raise ValueError(f"Cannot find CSQ format in the VCF header")
        m = re.search(r"Format: ([\+\-\w\|]+)", csq_info_header)
        if m:
            csq_format = m.group(1)
        else:
            raise ValueError(f"Cannot parse the CSQ field format from its INFO VCF header: {csq_info_header}")
        csq_fields = csq_format.split("|") 
    
        return csq_fields

    @classmethod
    def vcf_reader(cls: Type[V], File: str):
    
        with closing(VCF(File)) as vcf:
         #ACMG info tag and the new outfile
            vcf.add_info_to_header({'ID': 'ACMG verdict', 'Description': 'Variant prediction according to ACMG rules', 'Type':'Character','Number': '1'})
            vcf.add_info_to_header({'ID': 'Triggered ACMG rules', 'Description': 'Triggered ACMG rules', 'Type':'Character','Number': '1'})
        
            for variant in vcf:
                yield variant
