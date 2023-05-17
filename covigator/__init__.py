VERSION = "v2.1.0"
ANALYSIS_PIPELINE_VERSION = "v0.15.0"

MISSENSE_VARIANT = "missense_variant"
SYNONYMOUS_VARIANT = "synonymous_variant"
INFRAME_DELETION = "inframe_deletion"
INFRAME_INSERTION = "inframe_insertion"
DISRUPTIVE_INFRAME_DELETION = "disruptive_inframe_deletion"
CONSERVATIVE_INFRAME_DELETION = "conservative_inframe_deletion"
DISRUPTIVE_INFRAME_INSERTION = "disruptive_inframe_insertion"
CONSERVATIVE_INFRAME_INSERTION = "conservative_inframe_insertion"
INTERGENIC_VARIANT = "intergenic_variant"
INTERGENIC_DELETION = "intergenic_deletion"
FRAMESHIFT_VARIANT = "frameshift_variant"

CANONICAL_PROTEIN_MAPPING = {
            "S": ["S", "s", "spike"],
            "N": ["N", "n"],
            "M": ["M", "m"],
            "E": ["E", "e"],
            "ORF1ab": ["ORF1ab", "1ab", "orf1ab", "ORF1a", "ORF1b", "orf1a", "orf1b", "nsp3", "nsp6", "nsp5", "nsp4",
                       "nsp15", "nsp13", "nsp7", "nsp12", "NSP2"],
            "ORF3a": ["ORF3a", "orf3a"],
            "ORF8": ["ORF8", "8"],
            "ORF6": ["ORF6", "orf6"]
        }
