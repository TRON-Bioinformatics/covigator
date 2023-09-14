VERSION = "v2.2.2"
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
            "ORF1ab": ["ORF1ab", "1ab", "orf1ab", "ORF1a", "ORF1b", "orf1a", "orf1b", "nsp1", "NSP1", "nsp2", "NSP2",
                       "nsp3", "NSP3", "nsp4", "NSP4", "nsp5", "NSP5", "nsp6", "NSP6", "nsp7", "NSP7", "nsp8", "NSP8",
                       "nsp9", "NSP9", "nsp10", "NSP10", "nsp11", "NSP11", "nsp12", "NSP12", "nsp13", "NSP13", "nsp14",
                       "NSP14", "nsp15", "NSP15", "nsp16", "NSP16"],
            "ORF3a": ["ORF3a", "orf3a"],
            "ORF8": ["ORF8", "orf8", "8"],
            "ORF6": ["ORF6", "orf6"],
            "ORF7a": ["ORF7a", "orf7a"],
            "ORF7b": ["ORF7b", "orf7b"],
            "ORF10": ["ORF10", "orf10", "10"],

        }
