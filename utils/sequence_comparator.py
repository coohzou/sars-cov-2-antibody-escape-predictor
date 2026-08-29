"""
Sequence comparison module
Compares user sequences against SARS-CoV-2 reference genomes.
"""

import json
import logging
import os

from Bio import Align
from Bio import SeqIO

from utils.paths import MANIFEST_PATH, PREDICTION_DATA_DIR, prediction_fasta

logger = logging.getLogger(__name__)

SPIKE_START = 21563
SPIKE_END = 25384

# Broad WHO-style categories shown to users (not sub-lineages like BA.4 / JN.1).
MAJOR_VARIANT_NAMES = (
    "Wild Type",
    "Alpha",
    "Beta",
    "Gamma",
    "Delta",
    "Omicron",
)

# Map manifest / sub-lineage labels to a major category for display and metrics.
VARIANT_TO_MAJOR = {
    "Wild Type": "Wild Type",
    "Alpha": "Alpha",
    "Beta": "Beta",
    "Gamma": "Gamma",
    "Delta": "Delta",
    "Omicron": "Omicron",
    "BA.1": "Omicron",
    "BA.2": "Omicron",
    "BA.4": "Omicron",
    "BA.5": "Omicron",
    "BA.2.75.2": "Omicron",
    "BQ.1.1": "Omicron",
    "XBB.1.5": "Omicron",
    "JN.1": "Omicron",
}


def to_major_variant(variant_name):
    """Collapse sub-lineage names to Alpha–Omicron major categories."""
    if not variant_name:
        return None
    return VARIANT_TO_MAJOR.get(variant_name, variant_name)


class SequenceComparator:
    def __init__(self):
        self.covid_threshold = 0.70
        self.references = {}
        self.spike_references = {}
        self.variant_meta = {}
        self.prediction_dir = PREDICTION_DATA_DIR
        self._load_all_references()

        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = "global"
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -1

    def _resolve_path(self, item):
        return prediction_fasta(item["file"])

    def _load_manifest(self):
        if not os.path.exists(MANIFEST_PATH):
            return None
        with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
            return json.load(f)

    @staticmethod
    def _extract_spike_region(sequence):
        seq = sequence.replace("\n", "").replace(" ", "").upper()
        if len(seq) < SPIKE_END:
            return seq
        return seq[SPIKE_START - 1:SPIKE_END]

    def _load_all_references(self):
        manifest = self._load_manifest()
        if manifest and manifest.get("variants"):
            entries = manifest["variants"]
        else:
            entries = [
                {"name": "Wild Type", "file": "wuhan_hu1_complete.fasta", "in_variants_dir": False},
                {"name": "Alpha", "file": "alpha_complete.fasta", "in_variants_dir": True},
                {"name": "Beta", "file": "beta_complete.fasta", "in_variants_dir": True},
                {"name": "Gamma", "file": "gamma_complete.fasta", "in_variants_dir": True},
                {"name": "Delta", "file": "delta_complete.fasta", "in_variants_dir": True},
                {"name": "Omicron", "file": "omicron_complete.fasta", "in_variants_dir": True},
            ]

        for item in entries:
            name = item["name"]
            path = self._resolve_path(item)
            if not os.path.exists(path):
                logger.warning("%s reference not found at %s", name, path)
                continue
            try:
                records = list(SeqIO.parse(path, "fasta"))
                if records:
                    full_seq = str(records[0].seq).upper()
                    self.references[name] = full_seq
                    self.spike_references[name] = self._extract_spike_region(full_seq)
                    self.variant_meta[name] = {
                        "lineage": item.get("lineage", name),
                        "accession": item.get("accession", ""),
                        "file": item.get("file", os.path.basename(path)),
                    }
                    logger.info(
                        "Loaded %s reference (%d bp, spike %d bp)",
                        name,
                        len(full_seq),
                        len(self.spike_references[name]),
                    )
            except Exception as exc:
                logger.error("Failed to load %s: %s", name, exc)

    def calculate_identity(self, seq1, seq2):
        if not seq1 or not seq2:
            return 0

        seq1 = self._extract_spike_region(seq1)
        seq2 = self._extract_spike_region(seq2)

        alignments = self.aligner.align(seq1, seq2)
        if not alignments:
            return 0

        best_alignment = alignments[0]
        aligned_seq1 = str(best_alignment[0])
        aligned_seq2 = str(best_alignment[1])

        matches = 0
        total = 0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a != "-" and b != "-":
                total += 1
                if a == b:
                    matches += 1

        if total == 0:
            return 0

        return round((matches / total) * 100, 2)

    def is_sars_cov2(self, query_seq):
        query_clean = query_seq.replace("\n", "").replace(" ", "").upper()

        if not query_clean or not self.spike_references:
            return False, 0, None, "Reference not available"

        best_similarity = 0
        best_variant = None

        for variant_name, ref_spike in self.spike_references.items():
            if variant_name not in MAJOR_VARIANT_NAMES:
                continue
            similarity = self.calculate_identity(query_clean, ref_spike)
            if similarity > best_similarity:
                best_similarity = similarity
                best_variant = variant_name

        major_variant = to_major_variant(best_variant)
        if best_similarity >= self.covid_threshold * 100:
            return (
                True,
                best_similarity,
                major_variant,
                f"Confirmed SARS-CoV-2 ({major_variant}, {best_similarity}%)",
            )
        return False, best_similarity, major_variant, f"Not SARS-CoV-2 ({best_similarity}%)"

    def get_all_references(self):
        return self.references

    def get_variant_meta(self):
        return self.variant_meta


comparator = SequenceComparator()
