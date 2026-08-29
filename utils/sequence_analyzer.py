import os
import logging
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from .sequence_comparator import comparator
from .paths import prediction_fasta

logger = logging.getLogger(__name__)


class SequenceAnalyzer:
    def __init__(self):
        self.mutations_to_detect = {
            'D614G': {'protein_pos': 614, 'ref_aa': 'D', 'mut_aa': 'G'},
            'N501Y': {'protein_pos': 501, 'ref_aa': 'N', 'mut_aa': 'Y'},
            'E484K': {'protein_pos': 484, 'ref_aa': 'E', 'mut_aa': 'K'},
            'K417N': {'protein_pos': 417, 'ref_aa': 'K', 'mut_aa': 'N'},
            'K417T': {'protein_pos': 417, 'ref_aa': 'K', 'mut_aa': 'T'},
            'L452R': {'protein_pos': 452, 'ref_aa': 'L', 'mut_aa': 'R'},
            'T478K': {'protein_pos': 478, 'ref_aa': 'T', 'mut_aa': 'K'},
            'E484Q': {'protein_pos': 484, 'ref_aa': 'E', 'mut_aa': 'Q'},
            'P681R': {'protein_pos': 681, 'ref_aa': 'P', 'mut_aa': 'R'},
            'P681H': {'protein_pos': 681, 'ref_aa': 'P', 'mut_aa': 'H'},
            'H655Y': {'protein_pos': 655, 'ref_aa': 'H', 'mut_aa': 'Y'},
            'G339D': {'protein_pos': 339, 'ref_aa': 'G', 'mut_aa': 'D'},
            'S371L': {'protein_pos': 371, 'ref_aa': 'S', 'mut_aa': 'L'},
            'S371F': {'protein_pos': 371, 'ref_aa': 'S', 'mut_aa': 'F'},
            'S373P': {'protein_pos': 373, 'ref_aa': 'S', 'mut_aa': 'P'},
            'S375F': {'protein_pos': 375, 'ref_aa': 'S', 'mut_aa': 'F'},
            'N440K': {'protein_pos': 440, 'ref_aa': 'N', 'mut_aa': 'K'},
            'G446S': {'protein_pos': 446, 'ref_aa': 'G', 'mut_aa': 'S'},
            'S477N': {'protein_pos': 477, 'ref_aa': 'S', 'mut_aa': 'N'},
            'Q493R': {'protein_pos': 493, 'ref_aa': 'Q', 'mut_aa': 'R'},
            'Q498R': {'protein_pos': 498, 'ref_aa': 'Q', 'mut_aa': 'R'},
            'Y505H': {'protein_pos': 505, 'ref_aa': 'Y', 'mut_aa': 'H'},
            'N679K': {'protein_pos': 679, 'ref_aa': 'N', 'mut_aa': 'K'},
            'F486V': {'protein_pos': 486, 'ref_aa': 'F', 'mut_aa': 'V'},
            'R493Q': {'protein_pos': 493, 'ref_aa': 'R', 'mut_aa': 'Q'},
            'R346T': {'protein_pos': 346, 'ref_aa': 'R', 'mut_aa': 'T'},
            'K444T': {'protein_pos': 444, 'ref_aa': 'K', 'mut_aa': 'T'},
            'N460K': {'protein_pos': 460, 'ref_aa': 'N', 'mut_aa': 'K'},
            'F486S': {'protein_pos': 486, 'ref_aa': 'F', 'mut_aa': 'S'},
            'F486P': {'protein_pos': 486, 'ref_aa': 'F', 'mut_aa': 'P'},
            'L455F': {'protein_pos': 455, 'ref_aa': 'L', 'mut_aa': 'F'},
            'F456L': {'protein_pos': 456, 'ref_aa': 'F', 'mut_aa': 'L'},
            'Q493E': {'protein_pos': 493, 'ref_aa': 'Q', 'mut_aa': 'E'},
            'T19R': {'protein_pos': 19, 'ref_aa': 'T', 'mut_aa': 'R'},
            'G142D': {'protein_pos': 142, 'ref_aa': 'G', 'mut_aa': 'D'},
            'D796Y': {'protein_pos': 796, 'ref_aa': 'D', 'mut_aa': 'Y'},
            'Q954H': {'protein_pos': 954, 'ref_aa': 'Q', 'mut_aa': 'H'},
            'N969K': {'protein_pos': 969, 'ref_aa': 'N', 'mut_aa': 'K'},
        }

        self.ref_path = prediction_fasta("wuhan_hu1_complete.fasta")
        self.reference_genome = self._load_reference_genome()
        self.SPIKE_START = 21563
        self.SPIKE_END = 25384

        self.global_aligner = Align.PairwiseAligner()
        self.global_aligner.mode = 'global'
        self.global_aligner.match_score = 2
        self.global_aligner.mismatch_score = -1
        self.global_aligner.open_gap_score = -2
        self.global_aligner.extend_gap_score = -1

        self.protein_aligner = Align.PairwiseAligner()
        self.protein_aligner.mode = 'global'
        self.protein_aligner.match_score = 2
        self.protein_aligner.mismatch_score = -1
        self.protein_aligner.open_gap_score = -10
        self.protein_aligner.extend_gap_score = -0.5

        self.reference_spike_aa = self._build_reference_spike_aa()
        self.anchor_mutations = {"D614G", "N501Y", "E484K", "L452R", "K417N", "K417T", "P681R", "P681H"}

    def _build_reference_spike_aa(self):
        if not self.reference_genome:
            return ""
        spike_nt = self.reference_genome[self.SPIKE_START - 1:self.SPIKE_END]
        if len(spike_nt) < 3000:
            return ""
        spike_aa = str(Seq(spike_nt).translate(to_stop=False))
        if "MFV" in spike_aa[:20]:
            return spike_aa[spike_aa.find("MFV"):]
        return spike_aa

    def _translate_spike_aa(self, spike_nt):
        spike_aa = ""
        for frame in range(3):
            trial_nt = spike_nt[frame:(len(spike_nt[frame:]) // 3) * 3 + frame]
            trial_aa = str(Seq(trial_nt).translate(to_stop=False))
            if "MFV" in trial_aa[:20]:
                spike_aa = trial_aa[trial_aa.find("MFV"):]
                break
        if not spike_aa:
            spike_aa = str(Seq(spike_nt).translate(to_stop=False))
        return spike_aa

    def _load_reference_genome(self):
        if os.path.exists(self.ref_path):
            try:
                with open(self.ref_path, encoding="utf-8") as handle:
                    rec = next(SeqIO.parse(handle, "fasta"))
                return str(rec.seq).upper()
            except Exception as exc:
                logger.error("Failed to load reference genome: %s", exc)
        return None

    def read_fasta(self, file_path):
        try:
            return list(SeqIO.parse(file_path, 'fasta'))
        except Exception:
            return []

    def calculate_gc_content(self, sequence):
        if not sequence:
            return 0
        seq = sequence.upper()
        return round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2) if seq else 0

    @staticmethod
    def _is_valid_sequence(sequence):
        if len(sequence) < 1000:
            return False
        allowed = set("ACGTUN-")
        invalid = sum(1 for base in sequence if base not in allowed)
        return (invalid / len(sequence)) <= 0.05

    def _find_best_reference_match(self, query_seq):
        if not comparator.spike_references:
            return None

        is_covid, similarity, best_variant, _message = comparator.is_sars_cov2(query_seq)
        return {
            'is_covid': is_covid,
            'similarity': similarity,
            'variant_name': best_variant,
        }

    def _extract_query_spike_nt(self, query_genome):
        if not self.reference_genome:
            return None
        try:
            alignment = self.global_aligner.align(self.reference_genome, query_genome)[0]
        except Exception:
            return None

        aligned_ref, aligned_query = str(alignment[0]), str(alignment[1])
        ref_coord_to_align_index = {}
        ref_coord = 1

        for i, base in enumerate(aligned_ref):
            if base != '-':
                ref_coord_to_align_index[ref_coord] = i
                ref_coord += 1

        if self.SPIKE_START not in ref_coord_to_align_index or self.SPIKE_END not in ref_coord_to_align_index:
            if len(query_genome) >= self.SPIKE_END:
                return query_genome[self.SPIKE_START - 1:self.SPIKE_END].replace("-", "").upper()
            return None

        start_idx = ref_coord_to_align_index[self.SPIKE_START]
        end_idx = ref_coord_to_align_index[self.SPIKE_END]

        spike_nt = aligned_query[start_idx:end_idx + 1].replace('-', '').upper()
        return spike_nt if len(spike_nt) >= 3000 else None

    def _build_ref_to_query_pos_map(self, query_spike_aa):
        """Map Wuhan spike positions (1-indexed) to query sequence indices (0-indexed)."""
        pos_map = {}
        if not self.reference_spike_aa or not query_spike_aa:
            return pos_map

        try:
            alignment = self.protein_aligner.align(self.reference_spike_aa, query_spike_aa)[0]
        except Exception as exc:
            logger.debug("Protein alignment failed: %s", exc)
            return pos_map

        aligned_ref = str(alignment[0])
        aligned_query = str(alignment[1])
        ref_pos = 0
        query_idx = -1

        for ref_char, query_char in zip(aligned_ref, aligned_query):
            if ref_char != '-':
                ref_pos += 1
            if query_char != '-':
                query_idx += 1
            if ref_char != '-' and query_char != '-':
                pos_map[ref_pos] = query_idx

        return pos_map

    def _detect_mutations_via_alignment(self, spike_aa):
        """Detect mutations using reference-aligned coordinates (handles indels)."""
        actual_mutations = {}
        pos_map = self._build_ref_to_query_pos_map(spike_aa)
        if len(pos_map) < 100:
            return actual_mutations

        for name, info in self.mutations_to_detect.items():
            ref_pos = info['protein_pos']
            ref_idx = ref_pos - 1
            if ref_idx < 0 or ref_idx >= len(self.reference_spike_aa):
                continue
            if self.reference_spike_aa[ref_idx] != info['ref_aa']:
                continue

            query_idx = pos_map.get(ref_pos)
            if query_idx is None or query_idx >= len(spike_aa):
                continue
            if spike_aa[query_idx] != info['mut_aa']:
                continue

            actual_mutations[name] = info

        return actual_mutations

    def _detect_mutations_via_offset_voting(self, spike_aa):
        """Fallback detection for samples where alignment mapping is unreliable."""
        actual_mutations = {}
        candidate_details = []
        offsets_count = {}

        for name, info in self.mutations_to_detect.items():
            ref_pos = info['protein_pos']
            mut_aa = info['mut_aa']

            search_range = 8
            start = max(0, ref_pos - search_range - 1)
            end = min(len(spike_aa), ref_pos + search_range)
            local_aa = spike_aa[start:end]

            found_idx = local_aa.find(mut_aa)
            if found_idx != -1:
                actual_pos = found_idx + start + 1
                offset = actual_pos - ref_pos
                candidate_details.append({'name': name, 'offset': offset, 'info': info})
                if -10 <= offset <= 2:
                    weight = 3 if name in self.anchor_mutations else 1
                    offsets_count[offset] = offsets_count.get(offset, 0) + weight

        if not offsets_count:
            return actual_mutations

        main_offset = max(offsets_count, key=offsets_count.get)
        if offsets_count[main_offset] < 2:
            return actual_mutations

        logger.debug("Locked main offset: %s", main_offset)

        for cand in candidate_details:
            if abs(cand['offset'] - main_offset) > 1:
                continue

            info = cand['info']
            ref_pos = info['protein_pos']
            ref_idx = ref_pos - 1
            aa_index = ref_pos + main_offset - 1
            if aa_index < 0 or aa_index >= len(spike_aa):
                continue
            if spike_aa[aa_index] != info['mut_aa']:
                continue
            if ref_idx < 0 or ref_idx >= len(self.reference_spike_aa):
                continue
            if self.reference_spike_aa[ref_idx] != info['ref_aa']:
                continue

            actual_mutations[cand['name']] = cand['info']

        return actual_mutations

    def detect_spike_mutations(self, query_genome):
        spike_nt = self._extract_query_spike_nt(query_genome)
        if not spike_nt:
            return {}

        spike_aa = self._translate_spike_aa(spike_nt)
        pos_map = self._build_ref_to_query_pos_map(spike_aa)

        if len(pos_map) >= 100:
            return self._detect_mutations_via_alignment(spike_aa)

        return self._detect_mutations_via_offset_voting(spike_aa)

    def analyze_sequence_file(self, file_path):
        result = {'success': False, 'error': None}

        records = self.read_fasta(file_path)
        if not records:
            result['error'] = "Cannot read file"
            return result

        if len(records) > 1:
            result['warning'] = (
                f"Multiple sequences found; analyzed first record only ({records[0].id or 'record 1'})"
            )

        dna_seq = str(records[0].seq).upper()
        if not self._is_valid_sequence(dna_seq):
            result['error'] = "Invalid sequence content (expected DNA: A/C/G/T/N)"
            return result

        result['details'] = {
            'sequence_length': len(dna_seq),
            'gc_content': self.calculate_gc_content(dna_seq),
        }

        best_match = self._find_best_reference_match(dna_seq)
        if not best_match or not best_match['is_covid']:
            result['error'] = "Not recognized as SARS-CoV-2 (Similarity < 70%)"
            result['similarity_score'] = best_match['similarity'] if best_match else 0
            return result

        result['similarity_score'] = best_match['similarity']
        result['variant'] = best_match['variant_name']
        result['variant_confidence'] = "High" if best_match['similarity'] >= 95 else "Medium"

        detected_mutations = self.detect_spike_mutations(dna_seq)
        result['detected_mutations'] = detected_mutations
        result['has_d614g'] = 'D614G' in detected_mutations
        result['success'] = True

        return result
