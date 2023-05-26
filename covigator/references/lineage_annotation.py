import os
import json
import re
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from sqlalchemy.orm import Session
from covigator import SYNONYMOUS_VARIANT, MISSENSE_VARIANT, INTERGENIC_VARIANT, INTERGENIC_DELETION, \
    CONSERVATIVE_INFRAME_DELETION, DISRUPTIVE_INFRAME_DELETION, FRAMESHIFT_VARIANT, CANONICAL_PROTEIN_MAPPING
from covigator.database.model import Lineages, LineageDefiningVariants, LineageVariant, Gene, VariantType, VariantLevel
from covigator.database.queries import Queries
from logzero import logger
from datetime import datetime
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

@dataclass
class Hgvs:
    hgvs_p: str
    reference: str
    alternate: str
    position: int
    annotation: str
    protein: str


@dataclass
class MutationInfo:
    site: str
    variant_id: str
    variant_type: str
    position: int
    reference: str
    alternate: str
    length: int
    protein: str
    level: str
    hgvs_p: str
    annotation: str
    ambiguous_alternate: bool = False


class LineageAnnotationsLoader:

    LINEAGE_CONSTELLATION_DIRECTORY = "constellations/constellations/definitions"
    LINEAGE_CONSTELLATION_GENOME = "constellations/constellations/data/SARS-CoV-2.json"
    protein_mapping = {synonym: canonical for canonical, list_of_synonyms in CANONICAL_PROTEIN_MAPPING.items()
                       for synonym in list_of_synonyms}

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=session)

        self.lineage_constellation_directory = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.LINEAGE_CONSTELLATION_DIRECTORY)
        self.gene_df = pd.read_sql(self.session.query(Gene).order_by(Gene.start).statement, self.session.bind)
        self.genome = self.load_genome()

    def load_genome(self):
        """
        Load reference genome file from pangolin constellation repository
        """
        genome_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), self.LINEAGE_CONSTELLATION_GENOME)
        with open(genome_file, "r") as file_handle:
            constellation_genome = json.load(file_handle)
        return constellation_genome.get("genome")

    def find_gene(self, position: int) -> Tuple:
        """
        For a given nucleotide position find the overlapping gene and return the canonical name with start and end
        coordinates.
        """
        if position is None:
            return None, None, None
        gene_coordinates = self.gene_df.query("start <= @position & end > @position").get(
            ["name", "start", "end"])
        if gene_coordinates.shape[0] == 0:
            logger.info("Position does not overlap any gene...")
            return None, None, None
        gene_name = gene_coordinates.get("name").item()
        gene_start = gene_coordinates.get("start").item()
        gene_end = gene_coordinates.get("end").item()
        gene_name = LineageAnnotationsLoader.protein_mapping.get(gene_name)
        return gene_name, gene_start, gene_end

    def _build_hgvs_annotation(self, frameshift, delins, ):
        pass

    def _get_protein_position(self, position_in_cds):
        """
        Calculate position in amino acid sequence
        """
        if position_in_cds % 3 != 0:
            return position_in_cds // 3 + 1
        return position_in_cds // 3


    def _is_frameshift(self, deletion_length):
        """
        Given a deletion length, check if is frameshift deletion or not
        """
        return deletion_length % 3 != 0

    def _get_last_position_of_deletion(self, protein_position, deletion_length):
        """
        Given the deletion_length (nucleotides) and the protein position return the last position affected by the
        deletion
        """
        if self._is_frameshift(deletion_length):
            return protein_position

        return protein_position + deletion_length // 3


    def get_hgvs_from_nuc_deletion(self, variant: dict) -> Hgvs:
        """
        Translate pangolin nucleotide level deletion description into hgvs protein annotation. Returns the deletion
        as hgvs string and the annotation of this variant

        Args:
            variant (dict): A dictionary with keys for genomic position and length of deletion
        Returns:
            annotation (dict): A dictionary with the annotation and HGVS notation of the deletion
        """
        genomic_position = variant.get("genomic_position")
        deletion_length = variant.get("length")
        gene, gene_start, gene_end = self.find_gene(genomic_position)
        # If deletions falls into intergenic space return nucleotide level reference, alternate and position
        if gene is None:
            return Hgvs(
                hgvs_p=None,
                reference=self.genome[genomic_position - 1: genomic_position + deletion_length],
                alternate=self.genome[genomic_position - 1],
                position=genomic_position,
                annotation=INTERGENIC_DELETION,
                protein=gene)

        cds_start = gene_start - 1
        # Given a nucleotide position calculate position relative to CDS coordinates
        # minus one as positions given in df are inclusive on both sides. Deletion positions
        # are non-inclusive => the next genomic position is affected by the deletion
        position_in_cds = genomic_position - cds_start
        # Get position in the protein and calculate the coordinates of the corresponding codon
        protein_position = self._get_protein_position(position_in_cds)
        codon = [protein_position * 3 - 2, protein_position * 3 - 1, protein_position * 3]
        delins = False
        first_changed = False
        position_last = 0
        codon_last = []
        # Annotate deletions according to SO terms.
        annotation = CONSERVATIVE_INFRAME_DELETION

        # Example: ORF1ab p.S3675_F3677del
        # Position is last base of the codon (3674). The next AA position(s) are affected by the deletion
        if position_in_cds % 3 == 0:
            first_changed = True
            protein_position += 1
            # If the deletion is longer than one codon get last codon affected by deletion
            position_last = self._get_last_position_of_deletion(protein_position - 1,
                                                                deletion_length)

        # Deletion starts at first base of codon (second and third position of codon part of the deletion)
        elif position_in_cds % 3 == 1:
            # Created deleted and WT sequence and check if AA and the sequence of the codons has changed -->
            # If the AA at that position has not changed the following (deletion_length -1) AAs are deleted
            wt_codon = Seq(self.genome[cds_start + codon[0] - 1:
                                       cds_start + codon[2]])
            wt_aa = str(wt_codon.translate())
            mt_end = cds_start + codon[0] + deletion_length
            mt_codon = Seq(self.genome[cds_start + codon[0] - 1] +  # First base of WT codon
                           self.genome[mt_end: mt_end + 2])         # First and second base after deletion
            mt_aa = str(mt_codon.translate())
            # If the codon sequence (nucleotide) has not changed the deletion is annotated as conservative
            # else it is classified as disruptive. We follow here the classification of SNPEff to keep compatibility
            # with mutations called by the pipeline.
            annotation = CONSERVATIVE_INFRAME_DELETION if str(wt_codon) == str(mt_codon) \
                else DISRUPTIVE_INFRAME_DELETION
            if wt_aa == mt_aa:
                first_changed = True
                protein_position += 1
                # Get last codon affected by deletion
                position_last = self._get_last_position_of_deletion(protein_position - 1, deletion_length)
            else:
                # If the AA has changed get the last codon affected by the deletion and set inserted AA
                # Example: p.D119_I121delinsV
                position_last = self._get_last_position_of_deletion(protein_position, deletion_length)
                delins = True
                inserted_aa = mt_aa

        # Deletion starts at second base of codon (third base of codon part of deletion)
        elif position_in_cds % 3 == 2:
            # Created deleted and WT sequence and check if AA has changed -->
            # If the AA at that position has not changed the following (deletion_length - 1) AA are deleted
            wt_codon = Seq(self.genome[cds_start + codon[0] - 1:
                                       cds_start + codon[2]])
            wt_aa = str(wt_codon.translate())
            mt_end = cds_start + codon[1] + deletion_length
            mt_codon = Seq(self.genome[cds_start + codon[0] - 1] +  # first position of WT codon
                           self.genome[cds_start + codon[0]] +      # second position of WT codon
                           self.genome[mt_end])                     # first base after deletion
            mt_aa = str(mt_codon.translate())
            # If the codon sequence (nucleotide) has not changed the deletion is annotated as conservative
            # else it is classified as disruptive. We follow here the classification of SNPEff to keep compatibility
            # with mutations called by the pipeline.
            annotation = CONSERVATIVE_INFRAME_DELETION if str(wt_codon) == str(mt_codon) \
                else DISRUPTIVE_INFRAME_DELETION
            if mt_aa == wt_aa:
                first_changed = True
                protein_position += 1
                # If the deletion is longer than one codon get last codon affected by deletion
                position_last = self._get_last_position_of_deletion(protein_position - 1, deletion_length)
            else:
                # If the AA has changed get the last codon affected by the deletion and set inserted AA
                # Example: p.D119_I121delinsV
                position_last = self._get_last_position_of_deletion(protein_position, deletion_length)
                delins = True
                inserted_aa = mt_aa

        # The first and last codon position shown the HGVS annotation depend on the effect of the mutation
        # Calculate codon for the first affected amino acid
        if first_changed:
            codon = [protein_position * 3 - 2, protein_position * 3 - 1, protein_position * 3]
        first_aa = str(Seq(self.genome[cds_start + codon[0] - 1: cds_start + codon[2]]).translate())
        # Calculate codon for the last affected amino acid if deletion is longer than one codon
        # Check consists of comparing position_last with protein_position. If they are the same only a single codon
        # is affected and therefore the codon and the AA not required
        if position_last != protein_position:
            codon_last = [position_last * 3 - 2, position_last * 3 - 1, position_last * 3]
            last_aa = str(Seq(self.genome[cds_start + codon_last[0] - 1: cds_start + codon_last[2]]).translate())

        # Create HGVS notations for different effects of the mutation
        # HGVS for frameshift deletion
        if self._is_frameshift(deletion_length):
            # Redefine in case we have frameshift variant
            annotation = FRAMESHIFT_VARIANT
            hgvs = "p.{aa}{pos}fs".format(
                aa=first_aa,
                pos=protein_position
            )
            ref_aa = first_aa
        # HGVS for a site where the deletion leads to an insertion of a new AA
        elif delins:
            hgvs = "p.{first_aa}{first_pos}_{last_aa}{last_pos}delins{ins}".format(
                first_aa=first_aa,
                first_pos=protein_position,
                last_aa=last_aa,
                last_pos=position_last,
                ins=inserted_aa
            )
            # Return all AA that were deleted
            ref_aa = str(Seq(self.genome[cds_start + codon[0] - 1:
                                         cds_start + codon_last[2]]).translate())
        # HGVS for deletion of one of multiple codons
        else:
            # Deletion spans multiple codons --> hgvs consists of first and last affected AA
            if position_last - protein_position > 0:
                hgvs = "p.{first_aa}{first_pos}_{last_aa}{last_pos}del".format(
                    first_aa=first_aa,
                    first_pos=protein_position,
                    last_aa=last_aa,
                    last_pos=position_last
                )
                ref_aa = str(Seq(self.genome[cds_start + codon[0] - 1:
                                             cds_start + codon_last[2]]).translate())
            # Single codon deleted
            else:
                hgvs = "p.{aa}{pos}del".format(
                    aa=first_aa,
                    pos=protein_position
                )
                ref_aa = first_aa

        return Hgvs(
            hgvs_p=hgvs,
            reference=ref_aa,
            alternate="del",
            position=protein_position,
            annotation=annotation,
            protein=gene)

    def get_hgvs_from_nuc_snp(self, variant) -> Hgvs:
        """
        Translate pangolin nucleotide level SNPs description into hgvs protein annotation. Returns the variant
        as hgvs string and the annotation of this variant

        Args:
            variant (dict): A dictionary with keys for genomic position, reference and alternate base

        Returns:
            annotation (dict): A dictionary with the annotation and HGVS notation of the SNP
        """
        position = variant.get("genomic_position")
        ref_base = variant.get("reference")
        alt_base = variant.get("alternate")
        annotation = SYNONYMOUS_VARIANT
        gene, gene_start, gene_end = self.find_gene(position)
        # If mutation does not overlap any gene --> Intergenic
        if gene is None:
            return Hgvs(
                hgvs_p=None,
                position=position,
                reference=ref_base,
                alternate=alt_base,
                annotation=INTERGENIC_VARIANT,
                protein=gene)

        position_in_cds = position - (gene_start - 1)
        protein_position = self._get_protein_position(position_in_cds)
        # Get wild type and mutated codon sequence
        if position_in_cds % 3 == 0:
            wt_codon = self.genome[(position - 2) - 1: position]
            mt_codon = wt_codon[0:2] + alt_base
        elif position_in_cds % 3 == 1:
            wt_codon = self.genome[position - 1: position + 2]
            mt_codon = alt_base + wt_codon[1::]
        else:
            wt_codon = self.genome[(position - 1) - 1: position + 1]
            mt_codon = wt_codon[0] + alt_base + wt_codon[2]
        # Translate and check determine type of mutation
        wt_aa = str(Seq(wt_codon).translate())
        mt_aa = str(Seq(mt_codon).translate())
        if wt_aa != mt_aa:
            annotation = MISSENSE_VARIANT
        hgvs_p = "p.{ref}{aa_pos}{alt}".format(
            ref=wt_aa,
            aa_pos=protein_position,
            alt=mt_aa
        )
        return Hgvs(
            hgvs_p=hgvs_p,
            position=protein_position,
            reference=wt_aa,
            alternate=mt_aa,
            annotation=annotation,
            protein=gene)

    def _generate_variant_id(self, hgvs: Hgvs) -> str:
        """
        Generate variant id based on hgvs object
        """
        if hgvs.protein is None:
            return "{}:{}>{}".format(hgvs.position, hgvs.reference, hgvs.alternate)

        return "{}:{}{}{}".format(hgvs.protein, hgvs.reference, hgvs.position, hgvs.alternate)


    def _parse_mutation_sites(self, variant_string: str) -> List[Dict]:
        """
        Parse mutation site definition from pango constellation files. This code
        is based on the parser found in the utility script of scorpio. Mutations
        on nucleotide level are parsed and translated to protein level.
        """
        snp_pattern = re.compile('([ACTG]+)([0-9]+)([ACTG]+)')
        insertion_pattern = re.compile(r'(\w+):(\d+)\+([a-zA-Z]+)')
        aa_mutation = re.compile(r'([a-zA-Z-*]+)(\d+)([a-zA-Z-*]*)')
        mutation_list = []
        this_mut = variant_string.split(":")
        # Insertions -> eiter nucleotide or protein coordinates
        if "+" in variant_string:
            logger.warning("Skipping nucleotide insertion {}...".format(variant_string))
        # Deletions (nucleotide level) => translated to protein level
        elif this_mut[0] == "del":
            hgvs = self.get_hgvs_from_nuc_deletion({"genomic_position": int(this_mut[1]), "length": int(this_mut[2])})
            variant_id = self._generate_variant_id(hgvs)
            mutation_info = MutationInfo(
                site=variant_string,
                variant_id=variant_id,
                variant_type=VariantType.DELETION.name,
                position=hgvs.position,
                reference=hgvs.reference,
                alternate=hgvs.alternate,
                length=int(this_mut[2]),
                protein=hgvs.protein,
                level=VariantLevel.PROTEOMIC.name if hgvs.protein is not None else VariantLevel.NUCLEOTIDE.name,
                hgvs_p=hgvs.hgvs_p,
                annotation=hgvs.annotation
            )
            mutation_list.append(mutation_info)
        # SNPs (nucleotide level, can be non-synonymous, missense and intergenic) => translated to protein level
        elif this_mut[0] in ["snp", "nuc"]:
            match = re.match(snp_pattern, this_mut[1])
            if not match:
                logger.warning("Could not parse the following site definition: {}".format(variant_string))
                return [None]
            hgvs = self.get_hgvs_from_nuc_snp({"genomic_position": int(match[2]), "reference": match[1],
                                               "alternate": match[3]})
            variant_id = self._generate_variant_id(hgvs)
            mutation_info = MutationInfo(
                site=variant_string,
                variant_id=variant_id,
                variant_type=VariantType.SNV.name,
                position=hgvs.position,
                reference=hgvs.reference,
                alternate=hgvs.alternate,
                protein=hgvs.protein,
                level=VariantLevel.PROTEOMIC.name if hgvs.protein is not None else VariantLevel.NUCLEOTIDE.name,
                hgvs_p=hgvs.hgvs_p,
                length=1,
                annotation=hgvs.annotation
            )
            mutation_list.append(mutation_info)
        # Amino acid mutations => SNVs, MNVs, Deletions
        else:
            match = re.match(aa_mutation, this_mut[1])
            if not match:
                logger.warning("Could not parse the following site definition: {}".format(variant_string))
                return [None]
            protein = LineageAnnotationsLoader.protein_mapping.get(this_mut[0])
            reference = match[1]
            alternate = match[3]
            position = int(match[2])
            if alternate not in ['-', 'del']:
                # Split grouped protein mutations e.g. KR204RG into K204R & R205G
                if len(reference) == len(alternate) and len(reference) > 1:
                    for protein_pos, aa_mutations in enumerate(zip(reference, alternate), start=position):
                        variant_id = "{}:{}{}{}".format(protein, aa_mutations[0], protein_pos, aa_mutations[1])
                        annotation = MISSENSE_VARIANT if aa_mutations[0] != aa_mutations[1] else SYNONYMOUS_VARIANT
                        mutation_info = MutationInfo(
                            site=variant_id,
                            variant_id=variant_id,
                            hgvs_p="p.{}{}{}".format(aa_mutations[0], protein_pos, aa_mutations[1]),
                            variant_type=VariantType.SNV.name,
                            position=protein_pos,
                            reference=aa_mutations[0],
                            alternate=aa_mutations[1],
                            protein=protein,
                            length=len(aa_mutations[0]),
                            level=VariantLevel.PROTEOMIC.name,
                            annotation=annotation
                        )
                        mutation_list.append(mutation_info)
                else:
                    variant_id = "{}:{}{}{}".format(protein, reference, position, alternate)
                    annotation = MISSENSE_VARIANT if reference != alternate else SYNONYMOUS_VARIANT
                    mutation_info = MutationInfo(
                        site=variant_string,
                        variant_id=variant_id,
                        hgvs_p="p.{}{}{}".format(reference, position, alternate),
                        variant_type=VariantType.SNV.name,
                        position=position,
                        reference=reference,
                        alternate=alternate,
                        protein=protein,
                        ambiguous_alternate=True if alternate == "" else False,
                        length=len(reference),
                        level=VariantLevel.PROTEOMIC.name,
                        annotation=annotation
                    )
                    mutation_list.append(mutation_info)
            else:
                alternate = "del" if alternate == "-" else "del"
                variant_id = "{}:{}{}{}".format(protein, reference, position, alternate)
                # Generate HGVS notation for deletion
                if len(reference) > 1:
                    hgvs_p = "p.{}{}_{}{}del".format(reference[0], position, reference[-1],
                                                     position + len(reference)-1)
                else:
                    hgvs_p = "p.{}{}del".format(reference, position)
                mutation_info = MutationInfo(
                    site=variant_string,
                    variant_id=variant_id,
                    hgvs_p=hgvs_p,
                    variant_type=VariantType.DELETION.name,
                    position=position,
                    reference=reference,
                    alternate=alternate,
                    protein=protein,
                    length=len(reference),
                    level=VariantLevel.PROTEOMIC.name,
                    annotation=CONSERVATIVE_INFRAME_DELETION
                )
                mutation_list.append(mutation_info)
        return mutation_list

    @staticmethod
    def _create_constellation_pango_mapping(lineage_constellation: dict) -> Dict:
        """
        Creates a dictionary, mapping pango lineage identifiers to constellation labels. This mapping
        is required in find_parent_sites to find the constellation label for a parent lineage as
        parent lineages are given as pango ids.
        """
        lineage_mapping = {}
        for constellation, annotation in lineage_constellation.items():
            for pango_id in annotation["pangolin_lineage_list"]:
                lineage_mapping[pango_id] = constellation
        return lineage_mapping

    def _find_parent_sites(self, lineage: str, lineage_constellation: dict, pango_constellation_mapping: dict) -> List:
        """
        Include for each constellation also the mutations present up in the phylogenetic tree.
        """
        # Returns parent as pangolin identifier
        parent_pango = lineage_constellation[lineage].get("parent_lineage_id", None)
        # Get corresponding constellation label to look up mutations
        parent = pango_constellation_mapping.get(parent_pango, None)
        sites = lineage_constellation[lineage].get("lineage_mutations", None)
        if parent is None:
            return sites
        else:
            sites.extend([x for x in self._find_parent_sites(parent, lineage_constellation, pango_constellation_mapping)
                          if x not in sites])
        return sites

    def _find_constellation_files(self) -> List:
        """
        Glob all lineage constellation files found in constellation repository
        """
        json_files = list(Path(self.lineage_constellation_directory).glob('*.json'))
        return json_files

    def _read_constellation_files(self):
        """
        Parse lineage defining constellation JSON files and store information
        """
        lineage_constellation = {}
        # Iterate over all json constellation files
        for this_file in self._find_constellation_files():
            with open(this_file, "r") as fd:
                data = json.load(fd)
            # Get lineage annotation information
            pangolin_lineage_list = set(data["variant"].get("Pango_lineages", []))
            if not pangolin_lineage_list:
                # MRCA lineage info is either not present or literally an empty string ""
                lineage_name = data["variant"].get("mrca_lineage", "")
                if not lineage_name:
                    lineage_name = data["variant"].get("lineage_name", "")
                pangolin_lineage_list = set([lineage_name])
            constellation_label = data["label"]
            # Optional records not present in all lineage constellations
            who_label = data["variant"].get("WHO_label", None)
            phe_label = data["variant"].get("PHE_label", None)
            voc_date = None
            vui_date = None
            variant_of_concern = False
            variant_under_investigation = False
            tags = "|".join(data["tags"])
            # Extract VOC/VOI/VUI information from tags field
            for this_tag in data["tags"]:
                if this_tag.startswith("VOC-") or this_tag.startswith("VOC "):
                    variant_of_concern = True
                    if this_tag.startswith("VOC-"):
                        voc_date = datetime.strptime(this_tag.lstrip("VOC-"), "%y%b-%d")
                        voc_date = voc_date.strftime("%Y-%m-%d")
                    # edge case for B1.1.7 constellation file
                    else:
                        voc_date = datetime.strptime(this_tag.lstrip("VOC "), "%Y%m/%d")
                        voc_date = voc_date.strftime("%Y-%m-%d")
                # Some constellations bring both VUI and V tags, others only one of the two.
                # Therefore, we have to check both cases, even if this means that we overwrite the information once.
                if this_tag.startswith("V-"):
                    variant_under_investigation = True
                    vui_date = datetime.strptime(this_tag.lstrip("V-"), "%y%b-%d")
                    vui_date = vui_date.strftime("%Y-%m-%d")
                if this_tag.startswith("VUI-"):
                    variant_under_investigation = True
                    vui_date = datetime.strptime(this_tag.lstrip("VUI-"), "%y%b-%d")
                    vui_date = vui_date.strftime("%Y-%m-%d")

            parent_lineage_id = data["variant"].get("parent_lineage", None)
            # Drop incompatible lineage calls from pangolin identifier list
            # These calls are incompatible with the scorpio constellation
            incompatible_lineage_calls = set(data["variant"].get("incompatible_lineage_calls", set()))
            pangolin_lineage_list = pangolin_lineage_list - incompatible_lineage_calls
            lineage_mutations = [self._parse_mutation_sites(x) for x in data["sites"]]
            lineage_mutations = [mut for x in lineage_mutations for mut in x if mut is not None]
            lineage_constellation[constellation_label] = {
                "pangolin_lineage_list": pangolin_lineage_list,
                "who_label": who_label,
                "phe_label": phe_label,
                "voc_date": voc_date,
                "vui_date": vui_date,
                "variant_of_concern": variant_of_concern,
                "variant_under_investigation": variant_under_investigation,
                "tags": tags,
                "parent_lineage_id": parent_lineage_id,
                "lineage_mutations": lineage_mutations
            }
        return lineage_constellation

    def _fill_lineage_mutation_table(self, lineage_constellation: dict):
        """
        Store constellation information in database
        """
        count_lineages = 0
        count_mutations = 0
        all_mutations = []
        seen_mutations = set()

        # Get a list of pangolin lineage ids with a PHE label
        lineages_with_phe = [y.get("pangolin_lineage_list") for x, y in lineage_constellation.items()
                             if y.get("phe_label")]
        lineages_with_phe = {x for lineage in lineages_with_phe for x in lineage}
        for constellation, annotation in lineage_constellation.items():
            # Drop "duplicated" lineages, as we do not need to annotate them twice in the database.
            # This is the case when a separate constellation file was created for an already existing lineage,
            # just to store it with an additional mutation, e.g. B.1.1.7+E484K.json and B.1.617.2+K417N.
            # Because it is still the same lineage, we store in the table only the constellation with a given PHE label
            pango_lineages = annotation["pangolin_lineage_list"]
            if not annotation["phe_label"]:
                pango_lineages = annotation["pangolin_lineage_list"] - lineages_with_phe

            for pangolin_lineage_id in pango_lineages:
                # Multiple pangolin identifiers can be presented in a constellation.
                # That is the case, e.g., for the delta sub-lineages AY.1 and AY.2, that were created to track local
                # spread. For now, we store the pango-lineages together with the constellation label as unique keys
                lineage = Lineages(
                    pango_lineage_id=pangolin_lineage_id,
                    constellation_id=constellation,
                    who_label=annotation["who_label"],
                    phe_label=annotation["phe_label"],
                    parent_lineage_id=annotation["parent_lineage_id"],
                    voc_date=annotation["voc_date"],
                    vui_date=annotation["vui_date"],
                    variant_under_investigation=annotation["variant_under_investigation"],
                    variant_of_concern=annotation["variant_of_concern"],
                    tags=annotation["tags"]
                )
                self.session.add(lineage)
                self.session.commit()
                count_lineages += 1
            # Append all lineage defining mutations to a list of seen mutations
            for this_mut in annotation["lineage_mutations"]:
                if this_mut.variant_id not in seen_mutations:
                    seen_mutations.add(this_mut.variant_id)
                    all_mutations.append(this_mut)
        logger.info("Loaded into the database {} lineages".format(count_lineages))
        
        # Store all observed lineage defining mutations in database
        lineage_mutations = []
        for this_mut in all_mutations:
            mutation = LineageDefiningVariants(
                variant_id=this_mut.variant_id,
                hgvs=this_mut.hgvs_p,
                variant_type=this_mut.variant_type,
                protein=this_mut.protein,
                position=this_mut.position,
                reference=this_mut.reference,
                alternate=this_mut.alternate,
                ambiguous_alternate=this_mut.ambiguous_alternate,
                annotation=this_mut.annotation,
                variant_level=this_mut.level
            )
            lineage_mutations.append(mutation)
            count_mutations += 1
        # Add all mutations at one for performance reasons
        self.session.add_all(lineage_mutations)
        self.session.commit()
        logger.info("Loaded into the database {} mutations".format(count_mutations))

    def _fill_relation_ship_table(self, lineage_constellation):
        """
        Store lineage/mutation N:M relationships.
        """
        count_relationship = 0
        relationship = []
        # Collect mutations and lineages including it into a dictionary
        mutation_lineage_mapping = {}
        pango_constellation_mapping = self._create_constellation_pango_mapping(lineage_constellation)
        for constellation, annotation in lineage_constellation.items():
            # Find mutations from parental lineage
            all_mutations_of_lineage = self._find_parent_sites(constellation, lineage_constellation,
                                                               pango_constellation_mapping)
            for this_mut in all_mutations_of_lineage:
                var_id = this_mut.variant_id
                if var_id not in mutation_lineage_mapping.keys():
                    mutation_lineage_mapping[var_id] = annotation["pangolin_lineage_list"]
                else:
                    # Create union of pango lineages with the mutation
                    mutation_lineage_mapping[var_id] = mutation_lineage_mapping[var_id] | annotation["pangolin_lineage_list"]
        # Store lineage / mutation combinations in database
        for variant_id, lineages in mutation_lineage_mapping.items():
            for pango_id in lineages:
                lineage_variant = LineageVariant(
                    pango_lineage_id=pango_id,
                    variant_id=variant_id
                )
                relationship.append(lineage_variant)
                count_relationship += 1

        self.session.add_all(relationship)
        self.session.commit()
        logger.info("Loaded into the database {} lineage-variant relationships".format(count_relationship))

    def _update_sublineages_who(self):
        """
        Update lineage table to annotate sublineages of a VOC with the WHO label
        """
        update_counter = 0
        # Get lineage annotations from database
        query = self.session.query(Lineages.pango_lineage_id.label("pangolin_lineage"), Lineages.who_label,
                                   Lineages.parent_lineage_id)
        lineages = pd.read_sql(query.statement, self.session.bind)
        lineages_without_who = lineages[lineages.who_label.isna()]
        # Find and update sublineages of VOC
        for this_lineage in lineages_without_who.itertuples():
            # skip lineages without a parent as these are not sublineages from a VOC
            if pd.isnull(this_lineage.parent_lineage_id):
                continue
            who_label = self.queries.find_parent_who_label(this_lineage.pangolin_lineage, lineages)
            self.session.query(Lineages).filter(Lineages.pango_lineage_id == this_lineage.pangolin_lineage).\
                update({'who_label': who_label}, synchronize_session=False)
            self.session.commit()
            update_counter += 1
        logger.info("Updated {} sublineages to include WHO label from parental VOC".format(update_counter))

    def load_data(self):
        # Fill lineage annotation table
        lineage_constellation = self._read_constellation_files()

        self._fill_lineage_mutation_table(lineage_constellation)
        self._update_sublineages_who()
        self._fill_relation_ship_table(lineage_constellation)
