import customtkinter as ctk
from customtkinter import DrawEngine
DrawEngine.preferred_drawing_method = "polygon_shapes"  # Best for smooth corners
from tkinter import ttk, scrolledtext, filedialog, messagebox
from pathlib import Path
import json
import re
import copy
import os
import csv
from typing import List, Tuple, Dict
from collections import Counter, defaultdict
import primer3
import networkx as nx
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
import webbrowser


class PrimerGenerator:
    """Site-directed mutagenesis primer generator with improved 3'-extension handling."""

    def __init__(self, flank_size_bases=30):
        # Standard genetic code
        self.codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        self.stop_codons = {"TAA", "TAG", "TGA"}

        # Reverse codon table
        self.aa_to_codons = defaultdict(list)
        for codon, aa in self.codon_table.items():
            self.aa_to_codons[aa].append(codon)

        # Constants
        self.flank_size_bases = flank_size_bases
        self.flank_size_codons = flank_size_bases // 3
        self.MIN_PRIMER_LEN = 30
        self.MIN_EXTENSION_LEN = 15  # Minimum 3'-extension length
        self.TARGET_TM = 68.0
        self.MIN_TM = 55.0
        self.MAX_TM_DIFF = 3.0
        self.MAX_OVERLAP_LEN = 15  # Hard cap on overlap length

    def check_sequence_validity(self, dna_sequence: str, flank_size: int = 30) -> bool:
        """Validate DNA sequence structure and reading frame."""
        if len(dna_sequence) % 3 != 0:
            print(f"Error: Sequence length {len(dna_sequence)} is not a multiple of 3.")
            return False

        if len(dna_sequence) < (flank_size * 2 + 6):
            print("Error: Sequence too short for flanking regions.")
            return False

        # Check start codon
        start_codon_pos = dna_sequence[flank_size:flank_size + 3]
        if start_codon_pos != "ATG":
            print(f"Error: No ATG start codon after first {flank_size} bases, found '{start_codon_pos}'.")
            return False

        # Check stop codon
        stop_codon_pos = dna_sequence[-(flank_size + 3):-flank_size]
        if stop_codon_pos not in self.stop_codons:
            print(f"Error: No valid stop codon before last {flank_size} bases, found '{stop_codon_pos}'.")
            return False

        print("Sequence check passed ✓: correct ATG start and stop codon positions.")
        return True

    @staticmethod
    def _parse_mutation(mut_str: str) -> Tuple[str, int, str]:
        """Parse mutation string like 'A123B' into (original_aa, position, new_aa)."""
        match = re.match(r'([A-Z])(\d+)([A-Z])', mut_str.strip())
        if not match:
            raise ValueError(f"❌ Invalid mutation format: {mut_str}")
        return match.group(1), int(match.group(2)), match.group(3)

    def validate_mutations_against_sequence(self, dna_sequence: str, variant_list: List[List[str]]) -> None:
        """Validate mutations against wildtype sequence."""
        for variant in variant_list:
            for mut in variant:
                original_aa, pos, _ = self._parse_mutation(mut)

                # Adjust position for flanking regions
                adj_pos_codons = pos + self.flank_size_codons
                codon_start = (adj_pos_codons - 1) * 3
                codon_seq = dna_sequence[codon_start:codon_start + 3]

                wt_aa = self.codon_table.get(codon_seq, 'X')
                if wt_aa != original_aa:
                    raise ValueError(
                        f"❌ Mutation {mut} does not match WT residue '{wt_aa}' "
                        f"at coding AA position {pos} (codon {codon_seq})"
                    )

    def reverse_complement(self, sequence: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement[base] for base in sequence[::-1])

    def calculate_tm(self, sequence: str) -> float:
        """Calculate melting temperature using primer3."""
        if not sequence:
            return 0.0
        try:
            return round(primer3.calc_tm(sequence), 2)
        except Exception as e:
            print(f"Error calculating Tm for sequence '{sequence}': {e}")
            return 0.0

    def get_best_codon(self, amino_acid: str, original_codon: str) -> str:
        """Get optimal codon for amino acid, preferring original if available."""
        codons = self.aa_to_codons.get(amino_acid, ["NNN"])
        return original_codon if original_codon in codons else codons[0]

    @staticmethod
    def flatten_mutations(mutation_groups):
        flat = []
        for group in mutation_groups:
            if isinstance(group, (list, tuple)) and all(isinstance(i, tuple) and len(i) == 3 for i in group):
                flat.extend(group)
            else:
                flat.append(group)
        return flat

    def merge_close_mutations(self, mutations: List[Tuple[str, int, str]], max_distance_nt: int = 21) -> List[List[Tuple[str, int, str]]]:
        sorted_muts = sorted(mutations, key=lambda x: x[1])
        if not sorted_muts:
            return []

        merged = []
        current_group = [sorted_muts[0]]

        for mut in sorted_muts[1:]:
            last_position = current_group[-1][1]
            if (mut[1] - last_position) * 3 <= max_distance_nt:
                current_group.append(mut)
            else:
                merged.append(current_group)
                current_group = [mut]
        merged.append(current_group)
        return merged

    def _create_mutated_sequence(self, dna_sequence: str, mutations: List[Tuple[str, int, str]]) -> str:
        """Create mutated DNA sequence from list of mutations."""
        mutated_seq = list(dna_sequence)

        for aa, pos, new_aa in mutations:
            idx = ((pos - 1) + self.flank_size_codons) * 3
            original_codon = dna_sequence[idx:idx + 3]
            new_codon = self.get_best_codon(new_aa, original_codon)

            # Verify original amino acid
            if self.codon_table.get(original_codon, 'X') != aa:
                print(f"Warning: Original amino acid mismatch at position {pos}")

            mutated_seq[idx:idx + 3] = list(new_codon)

        return ''.join(mutated_seq)

    def _find_optimal_overlap(self, mutated_seq, codon_starts, codon_ends):
        min_inside_flank = 4
        crit_flank_limit = 2
        reverse_flank_weight = 2.0
        forward_flank_weight = 1.5
        wt_inside_weight = 0.5
        penalty_per_missing_flank_nt = 50
        upstream_bias_weight = 0.2

        mutation_positions = sorted(set(codon_starts))
        seq_len = len(mutated_seq)
        target_tm = self.TARGET_TM

        def search_len(curr_len):
            best_score = -1e9
            best_start, best_end, best_tm = None, None, None
            best_flanks = (None, None)
            for start_pos in range(0, seq_len - curr_len + 1):
                end_pos = start_pos + curr_len
                muts_inside = [m for m in mutation_positions if start_pos <= m < end_pos]
                if not muts_inside:
                    continue

                first_inside = min(muts_inside)
                last_inside = max(muts_inside)
                reverse_flank = first_inside - start_pos
                forward_flank = end_pos - (last_inside + 3)

                # reject critical short flanks
                if reverse_flank < crit_flank_limit or forward_flank < crit_flank_limit:
                    continue

                count = len(muts_inside)
                ext_span = self._extension_mutation_span(mutation_positions, start_pos, end_pos)
                overlap_seq = mutated_seq[start_pos:end_pos]
                overlap_tm = self.calculate_tm(overlap_seq)

                total_mut_nt = count * 3
                total_wt_nt = curr_len - total_mut_nt

                flank_penalty = 0
                if forward_flank < min_inside_flank:
                    flank_penalty += (min_inside_flank - forward_flank) * penalty_per_missing_flank_nt
                if reverse_flank < min_inside_flank:
                    flank_penalty += (min_inside_flank - reverse_flank) * penalty_per_missing_flank_nt

                score = (
                    count * 1000
                    + (forward_flank * forward_flank_weight)
                    + (reverse_flank * reverse_flank_weight)
                    + (total_wt_nt * wt_inside_weight)
                    - (ext_span * 10)
                    - flank_penalty
                    - abs(overlap_tm - target_tm)
                    - (start_pos * upstream_bias_weight)
                )

                if score > best_score:
                    best_score = score
                    best_start, best_end, best_tm = start_pos, end_pos, overlap_tm
                    best_flanks = (forward_flank, reverse_flank)

            return best_start, best_end, best_tm, best_flanks, best_score

        # --- Stage 1: try exact MAX_OVERLAP_LEN ---
        best_start, best_end, best_tm, (f_flank, r_flank), _ = search_len(self.MAX_OVERLAP_LEN)
        if best_start is not None and f_flank >= min_inside_flank and r_flank >= min_inside_flank:
            return best_start, best_end, best_tm

        # --- Stage 2: allow +1 and +2 bases, pick best scoring ---
        best_len_candidate = None
        best_len_score = -1e9
        for curr_len in [self.MAX_OVERLAP_LEN + 1, self.MAX_OVERLAP_LEN + 2]:
            cand_start, cand_end, cand_tm, _, cand_score = search_len(curr_len)
            if cand_start is not None and cand_score > best_len_score:
                best_len_candidate = (cand_start, cand_end, cand_tm)
                best_len_score = cand_score

        if best_len_candidate is not None:
            return best_len_candidate

        # --- Fallback: center between first/last mutation ---
        first_mut = min(mutation_positions)
        last_mut = max(mutation_positions)
        center = (first_mut + last_mut) // 2
        best_start = max(0, center - self.MAX_OVERLAP_LEN // 2)
        best_end = min(seq_len, best_start + self.MAX_OVERLAP_LEN)
        best_tm = self.calculate_tm(mutated_seq[best_start:best_end])
        return best_start, best_end, best_tm



    def _extension_mutation_span(self, mutation_positions, overlap_start, overlap_end):
        """
        Helper: compute span of mutations outside the overlap.
        Smaller = better.
        """
        outside = [pos for pos in mutation_positions if pos < overlap_start or pos >= overlap_end]
        if not outside:
            return 0
        return max(outside) - min(outside)

    def _build_primer(self, mutated_seq: str, overlap_seq: str, extension_start: int,
                      extension_direction: int, is_reverse: bool = False) -> Tuple[str, float]:
        """Build primer ensuring min 3'-extension length."""

        if is_reverse:
            extension_end = extension_start
            extension_start = max(0, extension_end - self.MIN_EXTENSION_LEN)
            extension_template = mutated_seq[extension_start:extension_end]
            extension = self.reverse_complement(extension_template)
            primer = self.reverse_complement(overlap_seq) + extension
        else:
            extension_end = min(len(mutated_seq), extension_start + self.MIN_EXTENSION_LEN)
            extension = mutated_seq[extension_start:extension_end]
            primer = overlap_seq + extension

        tm = self.calculate_tm(primer)
        return primer, tm

    def _extend_primer_one_base(self, mutated_seq: str, current_primer: str, overlap_seq: str,
                              extension_start: int, extension_end: int,
                              is_reverse: bool) -> Tuple[str, float, int, int]:

        MAX_PRIMER_LEN = 45  # define your max length here or use class attribute

        # Calculate current primer length
        current_len = len(current_primer)

        # Prevent extension if max length reached
        if current_len >= MAX_PRIMER_LEN:
            return current_primer, self.calculate_tm(current_primer), extension_start, extension_end

        if is_reverse:
            if extension_start == 0:
                # can't extend beyond start of sequence
                return current_primer, self.calculate_tm(current_primer), extension_start, extension_end

            extension_start = max(0, extension_start - 1)
            extension_template = mutated_seq[extension_start:extension_end]
            extension = self.reverse_complement(extension_template)
            primer = self.reverse_complement(overlap_seq) + extension
        else:
            if extension_end >= len(mutated_seq):
                # can't extend beyond end of sequence
                return current_primer, self.calculate_tm(current_primer), extension_start, extension_end

            extension_end = min(len(mutated_seq), extension_end + 1)
            extension = mutated_seq[extension_start:extension_end]
            primer = overlap_seq + extension

        tm = self.calculate_tm(primer)
        return primer, tm, extension_start, extension_end

    def count_nt_changes(self, wt_seq: str, mutations: List[Tuple[str, int, str]]) -> int:
        """Return the count of nucleotide changes vs wildtype after applying `mutations`."""
        mutated_seq = self._create_mutated_sequence(wt_seq, mutations)
        return sum(1 for a,b in zip(wt_seq, mutated_seq) if a != b)

    def split_mutations_by_max_nt_changes(
        self, wt_seq: str, mutation_group: List[Tuple[str, int, str]], max_nt_changes: int = 6
    ) -> List[List[Tuple[str, int, str]]]:
        """
        Split mutation group so that each group does not exceed max_nt_changes.
        Ensures codon-aware splits. Mutations are grouped as required to stay under the threshold.
        """
        sorted_mutations = sorted(mutation_group, key=lambda m: m[1])
        split_groups = []
        current_group = []
        for mut in sorted_mutations:
            test_list = current_group + [mut]
            if self.count_nt_changes(wt_seq, test_list) > max_nt_changes:
                if current_group:
                    split_groups.append(current_group)
                current_group = [mut]
            else:
                current_group.append(mut)
        if current_group:
            split_groups.append(current_group)
        return split_groups

    def generate_primers_for_construct(self, dna_sequence: str, mutations: List[Tuple[str, int, str]]) -> List[Dict]:
        mutations = self.flatten_mutations(mutations)
        primer_sets = []
        mutation_groups = self.merge_close_mutations(mutations)

        MAX_PRIMER_LEN = 45
        MAX_MUT_NT = 6  # Maximum nucleotide changes per primer pair

        for group in mutation_groups:
            # Split groups if nucleotide changes exceed max allowed per primer
            split_groups = self.split_mutations_by_max_nt_changes(dna_sequence, group, MAX_MUT_NT)
            previous_sub_group = []

            for i, sub_group in enumerate(split_groups):
                # Determine codon starts/ends for current mutations
                codon_starts_sub = [((m[1] - 1) + self.flank_size_codons) * 3 for m in sub_group]
                codon_ends_sub = [start + 3 for start in codon_starts_sub]

                # --- Accumulate previous mutations ---
                all_applied_muts = previous_sub_group + sub_group
                # Create mutated sequence including all previous mutations
                mutated_seq = self._create_mutated_sequence(dna_sequence, all_applied_muts)

                # Find optimal overlap for these mutations
                overlap_start, overlap_end, _ = self._find_optimal_overlap(mutated_seq, codon_starts_sub, codon_ends_sub)

                overlap_seq = mutated_seq[overlap_start:overlap_end]

                # Forward primer initial design
                forward_start = overlap_end
                forward_end = overlap_end + self.MIN_EXTENSION_LEN
                forward_primer = overlap_seq + mutated_seq[forward_start:forward_end]
                forward_tm = self.calculate_tm(forward_primer)

                # Reverse primer initial design
                reverse_end = overlap_start
                reverse_start = max(0, reverse_end - self.MIN_EXTENSION_LEN)
                reverse_template = mutated_seq[reverse_start:reverse_end]
                reverse_primer = self.reverse_complement(overlap_seq) + self.reverse_complement(reverse_template)
                reverse_tm = self.calculate_tm(reverse_primer)

                # Balance primer Tm with incremental extension
                tm_diff = abs(forward_tm - reverse_tm)
                attempts = 0
                max_attempts = 30
                while (
                    (tm_diff > self.MAX_TM_DIFF) or
                    (forward_tm < self.TARGET_TM) or
                    (reverse_tm < self.TARGET_TM)
                ) and attempts < max_attempts:
                    if forward_tm < reverse_tm and len(forward_primer) < MAX_PRIMER_LEN:
                        forward_primer, forward_tm, forward_start, forward_end = self._extend_primer_one_base(
                            mutated_seq, forward_primer, overlap_seq, forward_start, forward_end, is_reverse=False)
                    elif len(reverse_primer) < MAX_PRIMER_LEN:
                        reverse_primer, reverse_tm, reverse_start, reverse_end = self._extend_primer_one_base(
                            mutated_seq, reverse_primer, overlap_seq, reverse_start, reverse_end, is_reverse=True)
                    else:
                        break
                    tm_diff = abs(forward_tm - reverse_tm)
                    attempts += 1

                # Combine previous_sub_group + sub_group uniquely by mutation identity (aa, pos, new_aa)
                all_mut_set = { (m[0], m[1], m[2]) for m in previous_sub_group + sub_group }
                all_mut_list = sorted(all_mut_set, key=lambda x: x[1])

                # Create a unique ID string for this primer set from all covered mutations
                set_id = "_".join([f"{m[0]}{m[1]}{m[2]}" for m in all_mut_list])

                # Determine dependency: for the first primer set, no dependencies
                if len(primer_sets) == 0:
                    depends_on = []
                else:
                    # Depend on the immediate previous primer set
                    depends_on = [primer_sets[-1]['set_id']]

                primer_sets.append({
                    'set_id': set_id,
                    'depends_on': depends_on,
                    'mutations': [f"{m[0]}{m[1]}{m[2]}" for m in sub_group],  # current subgroup mutations
                    'all_covered_mutations': [f"{m[0]}{m[1]}{m[2]}" for m in all_mut_list],  # total mutations in primer sequence
                    'overlap_sequence': overlap_seq,
                    'overlap_tm': self.calculate_tm(overlap_seq),
                    'forward_primer': forward_primer,
                    'reverse_primer': reverse_primer,
                    'forward_tm': forward_tm,
                    'reverse_tm': reverse_tm,
                    'forward_length': len(forward_primer),
                    'reverse_length': len(reverse_primer),
                    'overlap_length': len(overlap_seq),
                    'tm_difference': tm_diff
                })

                previous_sub_group += sub_group

        return primer_sets



    def save_primer_list(self, primer_data: List[Dict], output_dir: str = None, filename: str = "primer_list.txt") -> None:
        """Save primer data to text and CSV files."""
        if output_dir is None:
            output_dir = os.path.expanduser("~/Mutagenesis")

        os.makedirs(output_dir, exist_ok=True)

        # Save text format
        filepath_txt = os.path.join(output_dir, filename)
        with open(filepath_txt, 'w') as f:
            f.write("Site-Directed Mutagenesis Primer List\n")
            f.write("=" * 120 + "\n\n")

            header = (f"{'Primer Name':<20}{'Primer Sequence':<60}{'Len':>6}"
                      f"{'Tm (°C)':>9}{'OverlapLen':>11}{'OverlapTm':>10}{'Notes':>30}\n")
            f.write(header)
            f.write("-" * 120 + "\n")

            for entry in primer_data:
                primer_base = "_".join(entry['all_covered_mutations'])
                notes = ",".join(entry['all_covered_mutations'])

                # Forward primer
                f.write(f"{primer_base + '_for':<20}")
                f.write(f"{entry['forward_primer']:<60}")
                f.write(f"{entry['forward_length']:>6}")
                f.write(f"{entry['forward_tm']:>9.1f}")
                f.write(f"{entry['overlap_length']:>11}")
                f.write(f"{entry['overlap_tm']:>10.1f}")
                f.write(f"{notes:>30}\n")

                # Reverse primer
                f.write(f"{primer_base + '_rev':<20}")
                f.write(f"{entry['reverse_primer']:<60}")
                f.write(f"{entry['reverse_length']:>6}")
                f.write(f"{entry['reverse_tm']:>9.1f}")
                f.write(f"{entry['overlap_length']:>11}")
                f.write(f"{entry['overlap_tm']:>10.1f}")
                f.write(f"{notes:>30}\n")
                f.write("-" * 120 + "\n")

            f.write(f"\nSummary:\nTotal primer pairs: {len(primer_data)}\n")

        # Save CSV format
        filepath_csv = os.path.join(output_dir, "primer_list.csv")
        with open(filepath_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "Primer Name", "Primer Sequence", "Length", "Tm (°C)",
                "Overlap Length", "Overlap Tm", "Mutations"
            ])

            for entry in primer_data:
                all_muts = ",".join(entry['all_covered_mutations'])
                base_name = "_".join(entry['all_covered_mutations'])
                writer.writerow([
                    f"{base_name}_for", entry['forward_primer'], entry['forward_length'],
                    entry['forward_tm'], entry['overlap_length'], entry['overlap_tm'], all_muts
                ])
                writer.writerow([
                    f"{base_name}_rev", entry['reverse_primer'], entry['reverse_length'],
                    entry['reverse_tm'], entry['overlap_length'], entry['overlap_tm'], all_muts
                ])

        print(f"Primer files saved:\n  - {filepath_txt}\n  - {filepath_csv}")

    def save_primer_data_json(self, primer_data: List[Dict], output_dir: str = None, filename: str = "primer_list.json") -> None:
        """Save primer data as JSON."""
        if output_dir is None:
            output_dir = os.path.expanduser("~/Mutagenesis")

        os.makedirs(output_dir, exist_ok=True)
        filepath_json = os.path.join(output_dir, filename)

        with open(filepath_json, 'w') as f:
            json.dump(primer_data, f, indent=2)

        print(f"Primer JSON saved to {filepath_json}")


# MutagenesisProtocol...

class MutagenesisProtocol:
    def __init__(
        self,
        variant_input,
        max_mutations_per_step,
        variant_prefix,
        output_dir_path,
        start_col=2,
        start_row=16,
        list_existing_as_steps=False,
        undo_stack=None
    ):
        self.variant_input = variant_input
        self.max_mutations_per_step = max_mutations_per_step
        self.variant_prefix = variant_prefix
        self.start_col = start_col
        self.start_row = start_row
        self.mutagenesis_dir = Path(output_dir_path)
        self.primer_json_path = self.mutagenesis_dir / "primer_list.json"
        self.output_dir = self.mutagenesis_dir / "protocols"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pdf_path = self.output_dir / "mutagenesis_protocol.pdf"
        self.databank_file = self.output_dir / "variant_databank.json"
        self.list_existing_as_steps = list_existing_as_steps
        self.undo_stack = undo_stack if undo_stack is not None else []

    @staticmethod
    def mutation_pos(mutation):
        matches = re.findall(r'\d+', mutation)
        return int(matches[0]) if matches else 0

    @staticmethod
    def variant_key(mutations):
        return ','.join(sorted(mutations, key=MutagenesisProtocol.mutation_pos))

    @staticmethod
    def extract_variant_number(name):
        match = re.search(r'variant (\d+)', name)
        return int(match.group(1)) if match else -1

    @staticmethod
    def normalize_mut(m):
        return m.strip().upper()

    def get_base_variant(self, current, existing_variants):
        """Find the best existing variant to use as a starting point"""
        best = ()
        for var in existing_variants:
            if set(var).issubset(current) and len(var) > len(best):
                best = var
        return best

    def design_primer_pathways(self, primer_list, target_variants):
        """
        Design optimal mutagenesis pathways using available primers.
        Only selects primers whose 'all_covered_mutations' are a subset of target mutations.

        primer_list: list of primers loaded from primer_list.json
        target_variants: list of target mutation sets, e.g., [["S19V"], ["P423F", "Q424V"]]
        Returns: list of workflows for each target variant
        """
        all_workflows = []

        for variant_idx, target_mutations in enumerate(target_variants):
            print(f"\nProcessing variant {variant_idx + 1}: {', '.join(target_mutations)}")

            steps = []
            target_muts = set(self.normalize_mut(m) for m in target_mutations)
            primer_map = {p['set_id']: p for p in primer_list}
            used_primers = []
            remaining_muts = target_muts.copy()

            def add_primer_with_deps(primer_id):
                # Recursively add dependencies first
                for dep in primer_map[primer_id].get('depends_on', []):
                    if dep not in used_primers:
                        add_primer_with_deps(dep)
                if primer_id not in used_primers:
                    used_primers.append(primer_id)
                    primer_info = primer_map[primer_id]
                    step_muts = [self.normalize_mut(m) for m in primer_info['all_covered_mutations']]
                    mutations_added = [self.normalize_mut(m) for m in primer_info['mutations']]
                    steps.append({
                        'primer_id': primer_id,
                        'primer_info': primer_info,
                        'all_covered_mutations': step_muts,
                        'mutations_added': mutations_added
                    })
                    remaining_muts.difference_update(mutations_added)

            while remaining_muts:
                # Find primers that only introduce mutations we want (subset of target)
                # and cover the most remaining mutations
                candidates = []
                for p in primer_list:
                    primer_all_muts = set(self.normalize_mut(m) for m in p['all_covered_mutations'])

                    # CRITICAL: Only consider primers whose all_covered_mutations are subset of target
                    if not primer_all_muts.issubset(target_muts):
                        continue

                    new_coverage = primer_all_muts & remaining_muts
                    if new_coverage:
                        candidates.append((p['set_id'], len(new_coverage), primer_all_muts, len(primer_all_muts)))

                if not candidates:
                    print(f"❌ Cannot cover remaining mutations: {remaining_muts}")
                    print("❌ No primers available that only introduce target mutations")
                    raise RuntimeError("Cannot cover all mutations with available primers that don't introduce unwanted mutations.")

                # Sort by: coverage of remaining mutations (desc), then total mutations (desc)
                candidates.sort(key=lambda x: (x[1], x[3]), reverse=True)
                best_primer_id = candidates[0][0]

                print(f"  Selected primer {best_primer_id} (covers {candidates[0][1]} remaining mutations)")
                add_primer_with_deps(best_primer_id)

            # Print the workflow
            print(f"\nWorkflow for variant {variant_idx + 1}:")
            print("Step  |  Primer Set  |  Mutations Added  |  All Covered")
            print("-" * 70)
            for step_num, step in enumerate(steps, 1):
                added_str = ', '.join(step['mutations_added'])
                covered_str = ', '.join(step['all_covered_mutations'])
                print(f"{step_num:2d}    |  {step['primer_id']:15s}  |  {added_str:15s}  |  {covered_str}")

            all_workflows.append(steps)

        return all_workflows
    
    def save_protocol_json(self, protocol_by_round, final_variants, existing_input_variants, 
                          variant_to_label_map, variant_to_final_variant, used_variants):
        """Save protocol data to JSON for GUI display"""
        protocol_data = {
            "existing_input_variants": existing_input_variants,
            "final_variants_overview": {},
            "protocol_steps": {},
            "all_variants": {}
        }
        
        # Final variants overview
        for label in sorted(variant_to_label_map, key=self.extract_variant_number):
            sorted_mutations = sorted(variant_to_label_map[label], key=self.mutation_pos)
            protocol_data["final_variants_overview"][label] = {
                "final_variant": variant_to_final_variant[label],
                "mutations": sorted_mutations
            }
        
        # Protocol steps
        for round_num in sorted(protocol_by_round):
            steps = sorted(protocol_by_round[round_num], key=lambda r: self.extract_variant_number(r[1]))
            step_data = []
            for row in steps:
                sorted_muts = sorted(row[2].split(', '), key=self.mutation_pos)
                step_data.append({
                    "new_variant": row[0],
                    "parent_variant": row[1], 
                    "mutations_added": sorted_muts,
                    "mutations_added_str": ', '.join(sorted_muts)
                })
            protocol_data["protocol_steps"][round_num] = step_data
        
        # All variants
        for muts, name in used_variants.items():
            mut_text = ', '.join(muts) if muts else '(none)'
            protocol_data["all_variants"][name] = mut_text
            
        # Save to JSON file
        protocol_json_path = self.output_dir / "protocol_data.json"
        with open(protocol_json_path, "w") as f:
            json.dump(protocol_data, f, indent=2)
        
        print(f"✅ Protocol data saved: {protocol_json_path}")

    def generate_pdf(
        self,
        protocol_by_round,
        final_variants,
        filename,
        existing_input_variants,
        variant_to_label_map,
        variant_to_final_variant,
        used_variants
    ):
        doc = SimpleDocTemplate(filename, pagesize=A4)
        elements = []
        styles = getSampleStyleSheet()
        cell_style = ParagraphStyle('mut_cell', fontSize=9, leading=11)

        if existing_input_variants:
            elements.append(Paragraph("<b>Pre-existing Variants in Databank</b>", styles['Heading2']))
            elements.append(Spacer(1, 12))
            data = [['Input Label', 'Existing Variant Label']]
            for input_label, existing_label in sorted(existing_input_variants.items()):
                data.append([input_label, existing_label or '(unknown)'])
            table = Table(data, colWidths=[150, 250])
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
            ]))
            elements.append(table)
            elements.append(PageBreak())

        if not protocol_by_round:
            elements.append(Paragraph("<b>No new variants generated.</b>", styles['Heading2']))
            elements.append(Spacer(1, 12))
            if not existing_input_variants:
                elements.append(Paragraph("No variants were generated because no input variants were provided.", styles['Normal']))
            else:
                elements.append(Paragraph("All input variants already exist in the databank.", styles['Normal']))
            doc.build(elements)
            return

        elements.append(Paragraph("<b>Final Variants Overview</b>", styles['Heading2']))
        final_table_data = [['Input Label', 'Final Variant', 'Mutations']]
        for label in sorted(variant_to_label_map, key=self.extract_variant_number):
            sorted_mutations = sorted(variant_to_label_map[label], key=self.mutation_pos)
            mut_paragraph = Paragraph(', '.join(sorted_mutations), cell_style)
            final_table_data.append([label, variant_to_final_variant[label], mut_paragraph])
        table = Table(final_table_data, colWidths=[100, 100, 260], repeatRows=1)
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
        ]))
        elements.append(table)
        elements.append(Spacer(1, 12))

        elements.append(Paragraph("<b>Step-by-Step Mutagenesis Protocol (Grouped by Round)</b>", styles['Heading2']))
        for i, round_num in enumerate(sorted(protocol_by_round), start=1):
            steps = sorted(protocol_by_round[round_num], key=lambda r: self.extract_variant_number(r[1]))
            elements.append(Paragraph(f"<b>Step {i}</b>", styles['Heading3']))
            table_data = [['New Variant', 'Parent Variant', 'Mutations Added']]
            for row in steps:
                sorted_muts = ', '.join(sorted(row[2].split(', '), key=self.mutation_pos))
                table_data.append([row[0], row[1], Paragraph(sorted_muts, cell_style)])
            step_table = Table(table_data, colWidths=[100, 120, 160], repeatRows=1)
            step_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
            ]))
            elements.append(step_table)
            elements.append(Spacer(1, 10))

        elements.append(PageBreak())
        elements.append(Paragraph("<b>All Produced Variants</b>", styles['Heading2']))
        all_variants_data = [['Variant Name', 'Mutations']]
        for muts, name in used_variants.items():
            mut_text = ', '.join(muts) if muts else '(none)'
            mut_paragraph = Paragraph(mut_text, cell_style)
            all_variants_data.append([name, mut_paragraph])
        table = Table(all_variants_data, colWidths=[140, 320], repeatRows=1)
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LEFTPADDING', (0, 0), (-1, -1), 5),
            ('RIGHTPADDING', (0, 0), (-1, -1), 5),
        ]))
        elements.append(table)

        doc.build(elements)

    def generate_label_pdf(self, variants, start_col, start_row, base_dir):
        base_dir = Path(base_dir)
        variant_json_file = base_dir / "protocols" / "variant_databank.json"
        pdf_output_file = base_dir / "protocols" / "herma_10900_labels.pdf"

        page_width, page_height = A4
        label_width = 25.3 * mm
        label_height = 9.9 * mm
        horizontal_spacing = 2.5 * mm
        top_margin = 13.5 * mm
        left_margin = 10.3 * mm
        columns = 7
        rows = 27
        font_size = 5.5

        def variant_to_key(mutations):
            return ','.join(sorted(mutations, key=self.mutation_pos))

        with open(variant_json_file, "r") as f:
            variant_data = json.load(f)

        def mutations_str_to_set(mutation_str):
            if mutation_str == '(none)':
                return set()
            return set(mutation_str.split(','))

        input_mutation_sets = [set(v) for v in variants]
        selected_variants = []
        for label, mutation_str in variant_data.items():
            variant_mut_set = mutations_str_to_set(mutation_str)
            if any(variant_mut_set.issubset(inp_set) for inp_set in input_mutation_sets):
                selected_variants.append(label)

        def variant_sort_key(label):
            match = re.search(r'\D+(\d+)$', label)
            return int(match.group(1)) if match else 0

        selected_variants = sorted(selected_variants, key=variant_sort_key)
        c = canvas.Canvas(str(pdf_output_file), pagesize=A4)
        c.setFont("Helvetica", font_size)

        label_index = 0
        total_labels = len(selected_variants)
        page_number = 1

        current_start_col = start_col
        current_start_row = start_row

        while label_index < total_labels:
            for row in range(rows):
                for col in range(columns):
                    if page_number == 1:
                        if row < current_start_row or (row == current_start_row and col < current_start_col):
                            continue
                    if label_index >= total_labels:
                        break

                    x = left_margin + col * (label_width + horizontal_spacing)
                    y = page_height - top_margin - (row + 1) * label_height

                    label = selected_variants[label_index]
                    text_x = x + label_width / 2
                    text_y = y + label_height / 2 - font_size / 2

                    c.drawCentredString(text_x, text_y, label)
                    label_index += 1

                if label_index >= total_labels:
                    break

            if label_index < total_labels:
                c.showPage()
                c.setFont("Helvetica", font_size)
                page_number += 1
                current_start_row = 0
                current_start_col = 0

        c.save()
        print(f"✅ PDF with labels generated: {pdf_output_file}")

    def run(self):
        # Check primer JSON exists
        if not self.primer_json_path.exists():
            raise FileNotFoundError(f"Primer JSON not found: {self.primer_json_path}")

        with open(self.primer_json_path, "r") as f:
            primer_data = json.load(f)

        mutation_to_groups = {}
        primer_groups = {}
        for idx, primer_set in enumerate(primer_data):
            group_muts = [self.normalize_mut(m) for m in primer_set["mutations"]]
            primer_groups[idx] = sorted(group_muts, key=self.mutation_pos)
            for mut in group_muts:
                mutation_to_groups.setdefault(mut, []).append(idx)

        if self.list_existing_as_steps:
            # SKIP DATABANK: Start fresh, only wildtype is "existing"
            variant_databank = {}
            existing_variants = {(): 'wildtype'}
        else:
            if self.databank_file.exists():
                try:
                    with open(self.databank_file, "r") as f:
                        content = f.read().strip()
                        variant_databank = json.loads(content) if content else {}
                except json.JSONDecodeError:
                    print("⚠️ Warning: Databank file invalid. Starting fresh.")
                    variant_databank = {}
            else:
                variant_databank = {}
            existing_variants = {(): 'wildtype'}
            for key, mutation_str in variant_databank.items():
                muts = tuple(mutation_str.split(',')) if mutation_str != '(none)' else ()
                existing_variants[muts] = key

        all_mutations = [mut for v in self.variant_input for mut in v]
        mutation_freq = Counter(all_mutations)
        variants = [sorted(v, key=lambda x: -mutation_freq[x]) for v in self.variant_input]

        existing_variants = {(): 'wildtype'}
        for key, mutation_str in variant_databank.items():
            muts = tuple(mutation_str.split(',')) if mutation_str != '(none)' else ()
            existing_variants[muts] = key

        final_variants = defaultdict(list)
        variant_to_label_map = {}
        variant_to_final_variant = {}
        protocol_by_round = defaultdict(list)
        used_variants = {}
        existing_input_variants = {}

        for idx, target_mutations in enumerate(variants):
            target = tuple(target_mutations)
            variant_key_str = self.variant_key(target)
            variant_label = f"variant {idx+1}"

            if variant_key_str in variant_databank.values():
                if not self.list_existing_as_steps:
                    # Just note it in the pre-existing list and skip real planning
                    for muts, name in existing_variants.items():
                        if variant_key_str == self.variant_key(muts):
                            existing_input_variants[variant_label] = name
                            break
                    continue
                else:
                    # Force regenerate a full protocol path as if it's new
                    # Remove the exact match from existing_variants so planning will run
                    for muts in list(existing_variants.keys()):
                        if self.variant_key(muts) == variant_key_str:
                            del existing_variants[muts]
                    # do NOT 'continue' here — fall through to full mutation batching logic

            base_variant = self.get_base_variant(target, existing_variants)
            current_mutations = list(base_variant)
            base_name = existing_variants[base_variant]

            target_group_ids = []
            for mut in target:
                mut_norm = self.normalize_mut(mut)
                possible_gids = mutation_to_groups.get(mut_norm, [])
                if possible_gids:
                    filtered_gids = [
                        gid for gid in possible_gids
                        if all(m in map(self.normalize_mut, target) for m in primer_groups[gid])
                        and any(m not in map(self.normalize_mut, current_mutations) for m in primer_groups[gid])
                    ]
                    if filtered_gids:
                        best_gid = max(
                            filtered_gids,
                            key=lambda gid: (
                                len([m for m in primer_groups[gid] if self.normalize_mut(m) in map(self.normalize_mut, target)]),
                                len(primer_groups[gid])
                            )
                        )
                        if best_gid not in target_group_ids:
                            target_group_ids.append(best_gid)

            needed_group_ids = []
            for gid in target_group_ids:
                group_muts = primer_groups[gid]
                if not all(m in map(self.normalize_mut, current_mutations) for m in group_muts):
                    needed_group_ids.append(gid)

            while needed_group_ids:
                batch_ids = needed_group_ids[:self.max_mutations_per_step]
                batch_muts = []
                for gid in batch_ids:
                    batch_muts.extend(primer_groups[gid])

                test_mutations = sorted(
                    set(map(self.normalize_mut, current_mutations)) | set(batch_muts),
                    key=self.mutation_pos
                )
                test_sorted = tuple(test_mutations)

                if test_sorted not in existing_variants:
                    mutation_count = len(test_sorted)
                    new_name = f"{variant_label}.{mutation_count}"
                    existing_variants[test_sorted] = new_name
                    used_variants[test_sorted] = new_name
                    protocol_by_round[mutation_count].append([new_name, base_name, ', '.join(batch_muts)])
                    base_name = new_name
                    current_mutations = test_mutations

                for gid in batch_ids:
                    needed_group_ids.remove(gid)

            final_variants[base_name].append(target)
            variant_to_label_map[variant_label] = target
            variant_to_final_variant[variant_label] = base_name
            used_variants[tuple(sorted(map(self.normalize_mut, target), key=self.mutation_pos))] = base_name

        self.generate_pdf(
            protocol_by_round, final_variants, str(self.pdf_path),
            existing_input_variants, variant_to_label_map,
            variant_to_final_variant, used_variants
        )

        self.save_protocol_json(
            protocol_by_round, final_variants, existing_input_variants,
            variant_to_label_map, variant_to_final_variant, used_variants
        )

        existing_ids = [int(v_id[-2:]) for v_id in variant_databank if v_id.startswith(self.variant_prefix)]
        next_id_num = max(existing_ids, default=0) + 1
        for muts, name in existing_variants.items():
            mutation_str = ','.join(sorted(muts, key=self.mutation_pos))
            if mutation_str not in variant_databank.values():
                new_id = f"{self.variant_prefix}{next_id_num:02d}"
                variant_databank[new_id] = mutation_str
                next_id_num += 1

        databank_file = self.output_dir / "variant_databank.json"
        if databank_file.exists():
            with open(databank_file, 'r') as f:
                try:
                    current_state = json.load(f)
                except json.JSONDecodeError:
                    current_state = {}
            self.undo_stack.append(copy.deepcopy(current_state))
        else:
            self.undo_stack.append({})  # empty dict if file does not exist

        if not self.list_existing_as_steps:
            with open(self.databank_file, "w") as f:
                json.dump(variant_databank, f, indent=2)

        # Now save new databank
        with open(databank_file, 'w') as f:
            json.dump(variant_databank, f, indent=2)

        print(f"✅ PDF generated: {self.pdf_path}")
        print(f"📁 Databank updated: {self.databank_file}")

        self.generate_label_pdf(self.variant_input, self.start_col, self.start_row, self.mutagenesis_dir)



#Tkinter GUI
# Set appearance mode and theme
ctk.set_appearance_mode("dark")  # "dark", "light", or "system"
ctk.set_default_color_theme("blue")  # "blue", "dark-blue", or "green"

class MutagenesisApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Site-Directed Mutagenesis Primer Designer")
        self.geometry("1050x500")
        self.minsize(1200, 700)  # Set minimum size instead of disabling resize
        
        # Configure main window grid
        self.grid_columnconfigure(0, weight=0, minsize=200)  # Navigation column
        self.grid_columnconfigure(1, weight=1)  # Main content column
        self.grid_rowconfigure(0, weight=1)

        self.LARGEFONT = ctk.CTkFont(family="Verdana", size=18, weight="bold")
        self.MEDIUMFONT = ctk.CTkFont(family="Verdana", size=12)
        self.SMALLFONT = ctk.CTkFont(family="Verdana", size=10)

        self.output_dir = ctk.StringVar(value=str(Path.cwd() / "Mutagenesis"))

        # Left navigation frame
        nav_frame = ctk.CTkFrame(self, width=200)
        nav_frame.grid(row=0, column=0, sticky="nsew", padx=(5,2), pady=5)
        nav_frame.grid_propagate(False)  # Maintain minimum width
        nav_frame.grid_columnconfigure(0, weight=1)

        # Navigation buttons
        nav_buttons = [
            ("User Input", InputPage),
            ("Variant Databank", DatabankPage),
            ("Wildtype Protein", WildtypeProteinPage),
            ("Variant Proteins", VariantProteinPage),  # Add this line
            ("Protocol Results", ProtocolResultsPage),
            ("Primer", PrimerPage)
        ]

        num_nav_buttons = len(nav_buttons)  # Your nav_buttons list
        nav_frame.grid_rowconfigure(num_nav_buttons, weight=1)  # Makes this row "flex" and fill space

        for i, (text, page_class) in enumerate(nav_buttons):
            btn = ctk.CTkButton(nav_frame, text=text, command=lambda p=page_class: self.show_frame(p))
            btn.grid(row=i, column=0, sticky="ew", pady=(15 if i == 0 else 5, 5), padx=15)

        # After creating nav_frame
        github_button = GithubButton(nav_frame, "https://github.com/xKosinus/Mutagenesis_Protocol/tree/windows_gui", text="GitHub Source")
        github_button.grid(row=num_nav_buttons + 1, column=0, sticky="ew", padx=15, pady=20)

        # Main container
        container = ctk.CTkFrame(self)
        container.grid(row=0, column=1, sticky="nsew", padx=(2,5), pady=5)
        container.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(0, weight=1)

        self.frames = {}
        for F in (InputPage, DatabankPage, WildtypeProteinPage, VariantProteinPage, ProtocolResultsPage, PrimerPage):  # Add VariantProteinPage
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.undo_stack = []
        self.show_frame(InputPage)

    def show_frame(self, page_class):
        frame = self.frames[page_class]
        frame.tkraise()
        # Update protein page when switching to it
        if page_class == WildtypeProteinPage:
            frame.update_protein_display()
        elif page_class == VariantProteinPage:  # Add this line
            frame.update_variant_display()      # Add this line  
        elif page_class == PrimerPage:
            frame.load_primers()

    def get_databank(self):
        path = Path(self.output_dir.get()) / "protocols" / "variant_databank.json"
        if path.exists():
            with open(path, "r") as f:
                return json.load(f)
        return {}

    def save_databank(self, databank):
        path = Path(self.output_dir.get()) / "protocols" / "variant_databank.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(databank, f, indent=2)

class PlaceholderTextbox(ctk.CTkTextbox):
    def __init__(self, master, placeholder, placeholder_fg="gray", text_fg="white", **kwargs):
        super().__init__(master, **kwargs)
        self.placeholder = placeholder
        self.placeholder_fg = placeholder_fg
        self.text_fg = text_fg
        self.insert("1.0", self.placeholder)
        self.configure(text_color=self.placeholder_fg)
        self.bind("<FocusIn>", self._clear_placeholder)
        self.bind("<FocusOut>", self._show_placeholder)

    def _clear_placeholder(self, event=None):
        text = self.get("1.0", "end").strip()
        if text == self.placeholder:
            self.delete("1.0", "end")
            self.configure(text_color=self.text_fg)

    def _show_placeholder(self, event=None):
        text = self.get("1.0", "end").strip()
        if text == "":
            self.insert("1.0", self.placeholder)
            self.configure(text_color=self.placeholder_fg)

class InputPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid weights for responsive layout
        self.grid_columnconfigure(0, weight=0, minsize=200)  # Label column
        self.grid_columnconfigure(1, weight=1)  # Input column
        self.grid_columnconfigure(2, weight=0)  # Button column
        self.grid_rowconfigure(10, weight=1)  # Output textbox row

        # Title label
        title = ctk.CTkLabel(self,
                             text="Mutagenesis Input",
                             corner_radius=10,
                             fg_color="#4a90e2",
                             text_color="white",
                             font=controller.LARGEFONT,
                             height=50)
        title.grid(row=0, column=0, padx=10, pady=20, columnspan=3, sticky="ew")

        # Output directory input
        ctk.CTkLabel(self, text="Output Directory:", font=controller.MEDIUMFONT).grid(
            row=1, column=0, sticky="w", padx=10, pady=5)
        self.output_dir_entry = ctk.CTkEntry(self, textvariable=controller.output_dir)
        self.output_dir_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=5)
        select_dir_btn = ctk.CTkButton(self, text="Choose...", command=self.select_output_dir, width=70)
        select_dir_btn.grid(row=1, column=2, sticky="w", padx=5, pady=5)

        # Label for DNA sequence input
        ctk.CTkLabel(self, text="Wildtype DNA Sequence:", font=controller.MEDIUMFONT).grid(
            row=2, column=0, sticky="nw", padx=10, pady=(15, 0))

        # Frame for DNA sequence input and info button (like var_frame)
        seq_frame = ctk.CTkFrame(self)
        seq_frame.grid(row=2, column=1, columnspan=2, pady=5, padx=5, sticky="ew")
        seq_frame.grid_columnconfigure(0, weight=1)

        # DNA sequence input textbox inside frame
        self.seq_text = PlaceholderTextbox(
            seq_frame,
            placeholder="Enter wildtype DNA sequence here...",
            placeholder_fg="gray",
            text_fg="white",
            height=120)
        self.seq_text.grid(row=0, column=0, sticky="ew", padx=(0, 5))

        # Info button inside frame, same row as textbox, next column
        question_seq_btn = ctk.CTkButton(seq_frame, text="?", width=30, height=30, command=self.show_info_sequence)
        question_seq_btn.grid(row=0, column=1, sticky="ne", pady=5)

        # Variant Mutations input
        ctk.CTkLabel(self, text="Mutations per Variant:", 
                    font=controller.MEDIUMFONT).grid(
            row=3, column=0, sticky="nw", padx=10, pady=(10, 0))
        
        # Frame for variant text and info button
        var_frame = ctk.CTkFrame(self)
        var_frame.grid(row=3, column=1, columnspan=2, pady=5, padx=5, sticky="ew")
        var_frame.grid_columnconfigure(0, weight=1)
        
        self.var_text = ctk.CTkTextbox(var_frame, height=120, font=controller.SMALLFONT)
        self.var_text.grid(row=0, column=0, sticky="ew", padx=(0, 5))
        self.var_text.insert("1.0", "S19V, S30V, Q424V, A431E\nK3Q, K36R, D76Q, L167A, Q424V\n")

        question_btn = ctk.CTkButton(var_frame, text="?", width=30, height=30, command=self.show_info_mutations)
        question_btn.grid(row=0, column=1, sticky="ne", pady=5)

        # Input parameters in a grid
        params_frame = ctk.CTkFrame(self)
        params_frame.grid(row=4, column=0, columnspan=3, sticky="ew", padx=10, pady=10)
        params_frame.grid_columnconfigure(1, weight=1)
        params_frame.grid_columnconfigure(3, weight=1)

        # Parameters
        param_configs = [
            ("Max mutations per step:", self.create_max_mut_combo, 0, 0),
            ("Variant prefix:", self.create_prefix_entry, 0, 2),
            ("Label Start Col:", self.create_start_col_entry, 1, 0),
            ("Label Start Row:", self.create_start_row_entry, 1, 2)
        ]

        for text, create_func, row, col in param_configs:
            ctk.CTkLabel(params_frame, text=text, font=controller.MEDIUMFONT).grid(
                row=row, column=col, sticky="w", padx=10, pady=5)
            widget = create_func(params_frame)
            widget.grid(row=row, column=col+1, sticky="w", padx=5, pady=5)

        # Checkbox
        self.skip_db = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(params_frame, text="Skip variant databank", variable=self.skip_db).grid(
            row=1, column=3, columnspan=2, pady=10, padx=10, sticky="e")

        # Action buttons
        button_frame = ctk.CTkFrame(self)
        button_frame.grid(row=9, column=0, columnspan=3, sticky="ew", padx=10, pady=10)
        
        run_btn = ctk.CTkButton(button_frame, text="Run Complete Workflow", 
                               font=controller.MEDIUMFONT, command=self.run_workflow)
        run_btn.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        
        undo_btn = ctk.CTkButton(button_frame, text="Undo Last Change", 
                                font=controller.MEDIUMFONT, command=self.undo_change)
        undo_btn.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        # Output textbox (expandable)
        self.output_text = ctk.CTkTextbox(self, font=controller.SMALLFONT)
        self.output_text.grid(row=10, column=0, columnspan=3, padx=10, pady=10, sticky="nsew")
        self.output_text.configure(state="disabled")

    def create_max_mut_combo(self, parent):
        self.max_mut_var = ctk.IntVar(value=1)
        return ctk.CTkComboBox(parent, values=["1", "2"], variable=self.max_mut_var, width=80)

    def create_prefix_entry(self, parent):
        self.var_prefix = ctk.StringVar(value="KWE_TA_A")
        return ctk.CTkEntry(parent, textvariable=self.var_prefix, width=150)

    def create_start_col_entry(self, parent):
        self.start_col = ctk.IntVar(value=0)
        return ctk.CTkEntry(parent, textvariable=self.start_col, width=80)

    def create_start_row_entry(self, parent):
        self.start_row = ctk.IntVar(value=0)
        return ctk.CTkEntry(parent, textvariable=self.start_row, width=80)

    def get_dna_sequence(self):
        """Get the DNA sequence from the input field, cleaned and uppercase"""
        seq = self.seq_text.get("1.0", "end").strip()
        if seq == "Enter wildtype DNA sequence here...":
            return ""
        return seq.replace('\n', '').upper()

    def show_info_mutations(self):
        messagebox.showinfo("Info",
                            "Enter mutations separated by commas, one variant per line.\nExample: S19V, S30V\nMultiple variants separated by new lines.")
        
    def show_info_sequence(self):
        messagebox.showinfo("Info",
                            "Include 30 bases at the flanking sites!")

    def select_output_dir(self):
        new_dir = filedialog.askdirectory(title="Select Output Directory", initialdir=self.controller.output_dir.get())
        if new_dir:
            self.controller.output_dir.set(new_dir)

    def run_workflow(self):
        # Enable textbox before writing
        self.output_text.configure(state="normal")
        self.output_text.delete("1.0", "end")
        
        seq = self.get_dna_sequence()
        variants_raw = self.var_text.get("1.0", "end").strip().split('\n')
        variant_list = [list(map(str.strip, line.split(','))) for line in variants_raw if line]
        max_mutations = self.max_mut_var.get()
        vp = self.var_prefix.get().strip()
        start_col = self.start_col.get()
        start_row = self.start_row.get()
        skip_db = self.skip_db.get()

        output_dir = Path(self.controller.output_dir.get())
        self.output_text.insert("end", f"Output directory: {output_dir}\n")
        self.output_text.insert("end", "=== RUNNING MUTAGENESIS WORKFLOW ===\n")
        self.output_text.insert("end", f"Sequence length: {len(seq)}\n")
        self.output_text.insert("end", f"Variants: {variant_list}\n")
        self.output_text.insert("end", f"Max mutations/step: {max_mutations}\n")
        self.output_text.insert("end", f"Variant prefix: {vp}\n")
        self.output_text.insert("end", f"Label START_COL: {start_col}, START_ROW: {start_row}\n\n")

        # ... rest of your workflow inserts
        self.output_text.insert("end", "✅ Primer design complete.\n\n")
        self.output_text.insert("end", "\n🎉 Workflow finished! Check your output directory for outputs.\n")

        # Lock textbox again
        self.output_text.configure(state="disabled")


        # ---- Primer design ----
        try:
            pg = PrimerGenerator(flank_size_bases=30)
            if not pg.check_sequence_validity(seq, flank_size=30):
                self.output_text.insert("end", "❌ Sequence failed validation.\n")
                return
            pg.validate_mutations_against_sequence(seq, variant_list)
        except Exception as e:
            self.output_text.insert("end", str(e) + "\n")
            return

        all_primer_data = []
        existing_pairs = set()
        for muts in variant_list:
            tuples = []
            for mut_str in muts:
                m = re.match(r'([A-Z])(\d+)([A-Z])', mut_str.strip())
                if m:
                    tuples.append((m.group(1), int(m.group(2)), m.group(3)))
            for primer_set in pg.generate_primers_for_construct(seq, tuples):
                key = (primer_set['forward_primer'], primer_set['reverse_primer'])
                if key not in existing_pairs:
                    all_primer_data.append(primer_set)
                    existing_pairs.add(key)
        pg.save_primer_list(all_primer_data, output_dir=output_dir)
        pg.save_primer_data_json(all_primer_data, output_dir=output_dir)
        self.output_text.insert("end", "✅ Primer design complete.\n\n")

        # ---- Protocol generation ----
        try:
            protocol = MutagenesisProtocol(
                variant_list,
                max_mutations,
                vp,
                output_dir,
                start_col,
                start_row,
                list_existing_as_steps=skip_db,
                undo_stack=self.controller.undo_stack
            )
            protocol.run()
            self.output_text.insert("end", "\n🎉 Workflow finished! Check your output directory for outputs.\n")
            self.controller.frames[DatabankPage].refresh_variants()
        except Exception as e:
            self.output_text.insert("end", str(e) + "\n")

    def undo_change(self):
        databank_path = Path(self.controller.output_dir.get()) / "protocols" / "variant_databank.json"
        self.output_text.configure(state="normal")
        try:
            if self.controller.undo_stack:
                last_state = self.controller.undo_stack.pop()
                with open(databank_path, "w") as f:
                    json.dump(last_state, f, indent=2)
                self.output_text.insert("end", "✅ Undid last change, variant databank restored.\n")
            else:
                self.output_text.insert("end", "No undo steps available.\n")
        except Exception as e:
            self.output_text.insert("end", "Undo failed: " + str(e) + "\n")
        self.output_text.configure(state="disabled")

class WildtypeProteinPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)  # Protein display area

        # Genetic code dictionary
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        # Header frame
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame,
                             text="Wildtype Protein Sequence",
                             corner_radius=10,
                             fg_color="#4a90e2",
                             text_color="white",
                             font=controller.LARGEFONT,
                             height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, padx=10, pady=10, sticky="e")

        # Info labels frame
        info_frame = ctk.CTkFrame(self)
        info_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        info_frame.grid_columnconfigure((0,1,2), weight=1)

        self.dna_length_label = ctk.CTkLabel(info_frame, text="DNA Length: N/A", font=controller.MEDIUMFONT)
        self.dna_length_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        self.coding_length_label = ctk.CTkLabel(info_frame, text="Coding Length: N/A", font=controller.MEDIUMFONT)
        self.coding_length_label.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        self.protein_length_label = ctk.CTkLabel(info_frame, text="Protein Length: N/A", font=controller.MEDIUMFONT)
        self.protein_length_label.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        # Format selection frame
        format_frame = ctk.CTkFrame(self)
        format_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=5)

        ctk.CTkLabel(format_frame, text="Display Format:", font=controller.MEDIUMFONT).grid(
            row=0, column=0, padx=10, pady=10, sticky="w")

        self.format_var = ctk.StringVar(value="1-letter")
        format_menu = ctk.CTkOptionMenu(format_frame, 
                                       values=["1-letter", "3-letter", "Both", "FASTA"],
                                       variable=self.format_var,
                                       command=self.update_display_format,
                                       width=120)
        format_menu.grid(row=0, column=1, padx=5, pady=10, sticky="w")

        refresh_btn = ctk.CTkButton(format_frame, text="Refresh", command=self.update_protein_display, width=80)
        refresh_btn.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        copy_btn = ctk.CTkButton(format_frame, text="Copy to Clipboard", command=self.copy_to_clipboard, width=120)
        copy_btn.grid(row=0, column=3, padx=10, pady=10, sticky="w")

        # Protein sequence display (expandable)
        self.protein_text = ctk.CTkTextbox(self, font=ctk.CTkFont(family="Courier", size=11))
        self.protein_text.grid(row=3, column=0, pady=10, padx=10, sticky="nsew")
        self.protein_text.configure(state="disabled")

        # Status label
        self.status_label = ctk.CTkLabel(self, text="", text_color="gray", font=controller.SMALLFONT)
        self.status_label.grid(row=4, column=0, sticky="w", padx=10, pady=5)

        # Store the current protein sequence
        self.current_protein_1letter = ""
        self.current_protein_3letter = ""

    def translate_dna_to_protein(self, dna_sequence):
        """Translate DNA sequence to protein sequence"""
        protein = []
        for i in range(0, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i:i+3]
            if codon in self.genetic_code:
                protein.append(self.genetic_code[codon])
            else:
                protein.append('X')
        return ''.join(protein)

    def get_three_letter_code(self, one_letter):
        """Convert one-letter amino acid code to three-letter code"""
        aa_map = {
            'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
            'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
            'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
            'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
            '*': 'Ter', 'X': 'Xxx'
        }
        return aa_map.get(one_letter, 'Xxx')

    def format_protein_sequence(self, sequence, format_type="1-letter", residues_per_line=60):
        """Format protein sequence for display"""
        if format_type == "1-letter":
            lines = []
            for i in range(0, len(sequence), residues_per_line):
                line_seq = sequence[i:i+residues_per_line]
                formatted_line = ""
                for j in range(0, len(line_seq), 10):
                    formatted_line += line_seq[j:j+10] + " "
                lines.append(f"{i+1:>4}: {formatted_line.strip()}")
            return '\n'.join(lines)

        elif format_type == "3-letter":
            three_letter = [self.get_three_letter_code(aa) for aa in sequence]
            lines = []
            residues_per_line_3letter = 20
            for i in range(0, len(three_letter), residues_per_line_3letter):
                line_seq = three_letter[i:i+residues_per_line_3letter]
                lines.append(f"{i+1:>4}: {'-'.join(line_seq)}")
            return '\n'.join(lines)

        elif format_type == "Both":
            result = "=== One-Letter Code ===\n"
            result += self.format_protein_sequence(sequence, "1-letter") + "\n\n"
            result += "=== Three-Letter Code ===\n"
            result += self.format_protein_sequence(sequence, "3-letter")
            return result

        elif format_type == "FASTA":
            header = ">Wildtype_Protein|Length=" + str(len(sequence))
            lines = [header]
            for i in range(0, len(sequence), 60):
                lines.append(sequence[i:i+60])
            return '\n'.join(lines)

        return sequence

    def update_protein_display(self):
        """Update the protein display based on current DNA input"""
        input_page = self.controller.frames[InputPage]
        dna_seq = input_page.get_dna_sequence()

        if not dna_seq or dna_seq == "":
            self.protein_text.configure(state="normal")
            self.protein_text.delete("1.0", "end")
            self.protein_text.insert("1.0", "No DNA sequence provided. Please enter a sequence in the Input page.")
            self.protein_text.configure(state="disabled")
            self.status_label.configure(text="No DNA sequence available")
            self.dna_length_label.configure(text="DNA Length: N/A")
            self.coding_length_label.configure(text="Coding Length: N/A")
            self.protein_length_label.configure(text="Protein Length: N/A")
            return

        # Remove 30bp flanking regions from both ends
        if len(dna_seq) > 60:
            coding_seq = dna_seq[30:-30]
            self.status_label.configure(text="✓ Removed 30bp flanking regions from both ends")
        else:
            self.protein_text.configure(state="normal")
            self.protein_text.delete("1.0", "end")
            self.protein_text.insert("1.0", "Error: DNA sequence is too short (must be >60bp with flanking regions)")
            self.protein_text.configure(state="disabled")
            self.status_label.configure(text="Error: Sequence too short", text_color="red")
            return

        # Translate to protein
        self.current_protein_1letter = self.translate_dna_to_protein(coding_seq)

        # Update info labels
        self.dna_length_label.configure(text=f"DNA Length: {len(dna_seq)}bp")
        self.coding_length_label.configure(text=f"Coding Length: {len(coding_seq)}bp")
        self.protein_length_label.configure(text=f"Protein Length: {len(self.current_protein_1letter)}aa")

        # Display formatted protein
        self.update_display_format()

    def update_display_format(self, *args):
        """Update the display format of the protein sequence"""
        if not self.current_protein_1letter:
            return

        format_type = self.format_var.get()
        formatted = self.format_protein_sequence(self.current_protein_1letter, format_type)

        self.protein_text.configure(state="normal")
        self.protein_text.delete("1.0", "end")
        self.protein_text.insert("1.0", formatted)
        self.protein_text.configure(state="disabled")

    def copy_to_clipboard(self):
        """Copy the current protein sequence to clipboard"""
        if self.current_protein_1letter:
            content = self.protein_text.get("1.0", "end-1c")
            self.clipboard_clear()
            self.clipboard_append(content)
            self.status_label.configure(text="✓ Copied to clipboard", text_color="green")
        else:
            self.status_label.configure(text="No protein sequence to copy", text_color="orange")

class VariantProteinPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)  # Protein display area

        # Genetic code dictionary for translation
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        # Header frame with title and back button
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)
        title = ctk.CTkLabel(header_frame,
                             text="Variant Protein Sequences",
                             corner_radius=10,
                             fg_color="#4a90e2",
                             text_color="white",
                             font=controller.LARGEFONT,
                             height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        back_btn = ctk.CTkButton(header_frame, text="Back to Databank", 
                                 command=lambda: controller.show_frame(controller.frames[DatabankPage]))
        back_btn.grid(row=0, column=1, padx=10, pady=10, sticky="e")

        # Info labels frame
        info_frame = ctk.CTkFrame(self)
        info_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        info_frame.grid_columnconfigure((0,1,2), weight=1)
        self.variant_count_label = ctk.CTkLabel(info_frame, text="Selected Variants: 0", font=controller.MEDIUMFONT)
        self.variant_count_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.total_sequences_label = ctk.CTkLabel(info_frame, text="Total Sequences: 0", font=controller.MEDIUMFONT)
        self.total_sequences_label.grid(row=0, column=1, padx=10, pady=10, sticky="w")
        self.avg_length_label = ctk.CTkLabel(info_frame, text="Avg Length: N/A", font=controller.MEDIUMFONT)
        self.avg_length_label.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        # Control buttons frame
        control_frame = ctk.CTkFrame(self)
        control_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=5)
        refresh_btn = ctk.CTkButton(control_frame, text="Refresh from Databank", 
                                   command=self.update_variant_display, width=150)
        refresh_btn.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        copy_btn = ctk.CTkButton(control_frame, text="Copy All to Clipboard", 
                                 command=self.copy_to_clipboard, width=150)
        copy_btn.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        # Protein sequences display (expandable)
        self.protein_text = ctk.CTkTextbox(self, font=ctk.CTkFont(family="Courier", size=11))
        self.protein_text.grid(row=3, column=0, pady=10, padx=10, sticky="nsew")
        self.protein_text.configure(state="disabled")

        # Status label
        self.status_label = ctk.CTkLabel(self, text="", text_color="gray", font=controller.SMALLFONT)
        self.status_label.grid(row=4, column=0, sticky="w", padx=10, pady=5)

        # Internal state
        self.current_variants = {}

    def translate_dna_to_protein(self, dna_sequence):
        protein = []
        for i in range(0, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i:i+3]
            protein.append(self.genetic_code.get(codon, 'X'))
        return "".join(protein)

    def apply_mutations_to_protein(self, wt_protein, mutation_str):
        protein_list = list(wt_protein)
        mutations = mutation_str.split(',') if mutation_str else []
        for mut in mutations:
            mut = mut.strip()
            if not mut:
                continue
            m = re.match(r'([A-Z])(\d+)([A-Z])', mut)
            if m:
                orig_aa, pos, new_aa = m.group(1), int(m.group(2)), m.group(3)
                idx = pos - 1
                if 0 <= idx < len(protein_list):
                    # Optional warning check
                    if protein_list[idx] != orig_aa:
                        print(f"Warning: WT AA at position {pos} is {protein_list[idx]}, expected {orig_aa}")
                    protein_list[idx] = new_aa
        return "".join(protein_list)

    def get_selected_variants(self):
        databank_page = self.controller.frames.get(DatabankPage)
        if not databank_page:
            return {}
        selected = {}
        for var_id, check_var in databank_page.check_vars:
            if check_var.get() == 1:
                selected[var_id] = True
        return selected

    def generate_variant_protein_sequences(self):
        input_page = self.controller.frames[InputPage]
        wt_dna = input_page.get_dna_sequence()
        if not wt_dna or len(wt_dna) <= 60:
            self.status_label.configure(text="Invalid or missing wildtype DNA sequence", text_color="red")
            return {}

        coding_dna = wt_dna[30:-30] if len(wt_dna) > 60 else wt_dna
        wt_protein = self.translate_dna_to_protein(coding_dna)

        databank = self.controller.get_databank()

        selected_variants = self.get_selected_variants()
        if not selected_variants:
            self.status_label.configure(text="No variants selected", text_color="orange")
            return {}

        variant_proteins = {}
        for var_id in selected_variants:
            muts_str = databank.get(var_id, "")
            variant_seq = self.apply_mutations_to_protein(wt_protein, muts_str)
            variant_proteins[var_id] = {
                'protein': variant_seq,
                'mutations': muts_str if muts_str else "(none)",
                'length': len(variant_seq)
            }

        return variant_proteins

    def format_variants_as_fasta(self, variant_proteins):
        if not variant_proteins:
            return "No variant sequences to display."
        lines = []
        for var_id, data in sorted(variant_proteins.items()):
            header = f">{var_id}|Mutations={data['mutations']}|Length={data['length']}aa"
            lines.append(header)
            seq = data['protein']
            for i in range(0, len(seq), 60):
                lines.append(seq[i:i+60])
            lines.append("")  # blank line between sequences
        return "\n".join(lines)

    def update_variant_display(self):
        self.status_label.configure(text="Generating variant protein sequences...", text_color="blue")
        variant_proteins = self.generate_variant_protein_sequences()
        self.current_variants = variant_proteins

        self.protein_text.configure(state="normal")
        self.protein_text.delete("1.0", "end")

        if not variant_proteins:
            self.protein_text.insert("1.0", "No variant sequences to display.\n\nPlease:\n1. Provide WT DNA in Input page.\n2. Select variants in Databank page.\n3. Click 'Refresh from Databank'.")
            self.variant_count_label.configure(text="Selected Variants: 0")
            self.total_sequences_label.configure(text="Total Sequences: 0")
            self.avg_length_label.configure(text="Avg Length: N/A")
            self.status_label.configure(text="No variant sequences", text_color="orange")
            self.protein_text.configure(state="disabled")
            return

        fasta_content = self.format_variants_as_fasta(variant_proteins)
        self.protein_text.insert("1.0", fasta_content)
        self.protein_text.configure(state="disabled")

        # Update info labels
        cnt = len(variant_proteins)
        lengths = [data["length"] for data in variant_proteins.values()]
        avg_len = sum(lengths) / len(lengths) if lengths else 0

        self.variant_count_label.configure(text=f"Selected Variants: {cnt}")
        self.total_sequences_label.configure(text=f"Total Sequences: {cnt}")
        self.avg_length_label.configure(text=f"Avg Length: {avg_len:.1f} aa")

        self.status_label.configure(text=f"✓ Generated {cnt} variant protein sequences", text_color="green")

    def copy_to_clipboard(self):
        if not self.current_variants:
            self.status_label.configure(text="No variant sequences to copy", text_color="orange")
            return
        content = self.protein_text.get("1.0", "end-1c")
        self.clipboard_clear()
        self.clipboard_append(content)
        self.status_label.configure(text="✓ Copied all sequences to clipboard", text_color="green")

class DatabankPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.last_deleted = {}

        # Configure main grid for responsive layout
        self.grid_columnconfigure(0, weight=1, minsize=400)  # Variant list
        self.grid_columnconfigure(1, weight=1, minsize=400)  # Search controls
        self.grid_rowconfigure(1, weight=1)  # Main content row

        # Internal state for variants and filtering
        self.check_vars = []
        self.checkboxes = []
        self.all_variants = {}
        self.filtered_variants = {}

        # Header frame
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Variant Databank",
                             corner_radius=15, fg_color="#4a90e2", text_color="white",
                             font=controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, sticky="e", padx=10, pady=10)

        # Main content frame
        content_frame = ctk.CTkFrame(self)
        content_frame.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=10, pady=5)
        content_frame.grid_columnconfigure(0, weight=1)  # Variant list
        content_frame.grid_columnconfigure(1, weight=1)  # Search section
        content_frame.grid_rowconfigure(0, weight=1)

        # Variant list frame (left side)
        variant_frame = ctk.CTkFrame(content_frame)
        variant_frame.grid(row=0, column=0, sticky="nsew", padx=(5, 2), pady=5)
        variant_frame.grid_columnconfigure(0, weight=1)
        variant_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(variant_frame, text="Variants:", font=controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        self.variant_listbox = ctk.CTkScrollableFrame(variant_frame)
        self.variant_listbox.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)

        # Search section (right side)
        search_main_frame = ctk.CTkFrame(content_frame)
        search_main_frame.grid(row=0, column=1, sticky="nsew", padx=(2, 5), pady=5)
        search_main_frame.grid_columnconfigure(0, weight=1)
        search_main_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(search_main_frame, text="Search & Filter:", font=controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        self.search_frame = ctk.CTkScrollableFrame(search_main_frame)
        self.search_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        self.search_frame.grid_columnconfigure(0, weight=1)

        self.search_fields = []
        self.search_vars = []

        # Search options frame - REORGANIZED LAYOUT
        options_frame = ctk.CTkFrame(self.search_frame)
        options_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        options_frame.grid_columnconfigure(5, weight=1)  # Make column 5 expandable

        # Row 0: Mutation count filter (now first)
        ctk.CTkLabel(options_frame, text="Mutations:", font=controller.SMALLFONT).grid(
            row=0, column=0, padx=5, pady=5, sticky="w")

        # Min mutations entry
        self.min_mutations_var = ctk.StringVar(value="")
        min_entry = ctk.CTkEntry(options_frame, textvariable=self.min_mutations_var, 
                               placeholder_text="Min", width=50)
        min_entry.grid(row=0, column=1, padx=2, pady=5, sticky="w")
        self.min_mutations_var.trace("w", self.on_search_change)

        ctk.CTkLabel(options_frame, text="to", font=controller.SMALLFONT).grid(
            row=0, column=2, padx=2, pady=5, sticky="w")

        # Max mutations entry
        self.max_mutations_var = ctk.StringVar(value="")
        max_entry = ctk.CTkEntry(options_frame, textvariable=self.max_mutations_var, 
                               placeholder_text="Max", width=50)
        max_entry.grid(row=0, column=3, padx=2, pady=5, sticky="w")
        self.max_mutations_var.trace("w", self.on_search_change)

        # Clear mutation filter button
        clear_mut_btn = ctk.CTkButton(options_frame, text="Clear", 
                                     command=self.clear_mutation_filter, width=60)
        clear_mut_btn.grid(row=0, column=4, padx=10, pady=5, sticky="w")

        # Row 1: Logic controls grouped together (now second)
        ctk.CTkLabel(options_frame, text="Logic:", font=controller.SMALLFONT).grid(
            row=1, column=0, padx=5, pady=5, sticky="w")

        self.search_logic_var = ctk.StringVar(value="AND")
        logic_menu = ctk.CTkOptionMenu(options_frame, values=["AND", "OR"], 
                                      variable=self.search_logic_var,
                                      command=self.on_search_change, width=70)
        logic_menu.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        # Add/Remove buttons next to Logic
        self.add_btn = ctk.CTkButton(options_frame, text="+", 
                                    command=self.add_search_field, width=30, height=30)
        self.add_btn.grid(row=1, column=2, padx=2, pady=5, sticky="w")

        self.remove_btn = ctk.CTkButton(options_frame, text="−", 
                                       command=self.remove_search_field,
                                       width=30, height=30, state="disabled")
        self.remove_btn.grid(row=1, column=3, padx=2, pady=5, sticky="w")

        # Clear All button next to +/- buttons
        clear_all_btn = ctk.CTkButton(options_frame, text="Clear All", 
                                     command=self.clear_all_search, width=80)
        clear_all_btn.grid(row=1, column=4, padx=10, pady=5, sticky="w")

        # Row 2: Help text
        help_label = ctk.CTkLabel(options_frame, text="Use AND to match all terms, OR for any. Filter by mutation count range (e.g., 2 to 4).",
                                 font=controller.SMALLFONT, text_color="gray")
        help_label.grid(row=2, column=0, columnspan=6, padx=5, pady=(0, 5), sticky="w")

        # Results label
        self.results_label = ctk.CTkLabel(self, text="", text_color="gray")
        self.results_label.grid(row=2, column=0, columnspan=2, sticky="w", padx=10, pady=(5, 0))

        # Action buttons frame
        button_frame = ctk.CTkFrame(self)
        button_frame.grid(row=3, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        
        buttons = [
            ("Delete Selected Variants", self.delete_variants, 0),
            ("Select All", self.select_all, 1),
            ("Deselect All", self.deselect_all, 2),
            ("Undo Delete", self.undo_delete, 3)
        ]
        
        for text, command, col in buttons:
            btn = ctk.CTkButton(button_frame, text=text, command=command)
            btn.grid(row=0, column=col, sticky="w", padx=10, pady=10)

        refresh_btn = ctk.CTkButton(button_frame, text="Refresh", command=self.refresh_variants)
        refresh_btn.grid(row=0, column=4, sticky="w", padx=10, pady=10)

        view_proteins_btn = ctk.CTkButton(button_frame, text="View Proteins", 
                                         command=lambda: controller.show_frame(VariantProteinPage))
        view_proteins_btn.grid(row=0, column=5, sticky="w", padx=10, pady=10)

        # Status label
        self.status_box = ctk.CTkLabel(self, text="", text_color="green")
        self.status_box.grid(row=4, column=0, columnspan=2, sticky="w", padx=10, pady=(0, 10))

        # Initialize with first search field
        self.add_search_field(is_first=True)
        self.refresh_variants()

    def add_search_field(self, is_first=False):
        field_num = len(self.search_fields)
        search_var = ctk.StringVar()
        search_var.trace("w", self.on_search_change)
        self.search_vars.append(search_var)
        
        field_frame = ctk.CTkFrame(self.search_frame)
        field_frame.grid(row=field_num + 2, column=0, sticky="ew", padx=5, pady=2)
        field_frame.grid_columnconfigure(1, weight=1)
        
        label_text = "Search:" if is_first else f"Search:"
        label = ctk.CTkLabel(field_frame, text=label_text, font=self.controller.SMALLFONT)
        label.grid(row=0, column=0, sticky="w", padx=5, pady=5)
        
        search_entry = ctk.CTkEntry(field_frame, textvariable=search_var,
                                   placeholder_text="Search by variant ID or mutations...")
        search_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
        
        field_data = {
            'frame': field_frame, 
            'label': label, 
            'entry': search_entry, 
            'var': search_var
        }
        self.search_fields.append(field_data)
        
        self.remove_btn.configure(state="normal" if len(self.search_fields) > 1 else "disabled")
        self.filter_variants()

    def remove_search_field(self):
        if len(self.search_fields) <= 1:
            return
        
        last_field = self.search_fields.pop()
        self.search_vars.pop()
        last_field['frame'].destroy()
        
        self.remove_btn.configure(state="normal" if len(self.search_fields) > 1 else "disabled")
        self.filter_variants()

    def clear_all_search(self):
        for search_var in self.search_vars:
            search_var.set("")

    def clear_mutation_filter(self):
        self.min_mutations_var.set("")
        self.max_mutations_var.set("")

    def count_mutations(self, mutation_string):
        """Count the number of mutations in a mutation string"""
        if not mutation_string or mutation_string == "(none)":
            return 0
        return len([m for m in mutation_string.split(',') if m.strip()])

    def on_search_change(self, *args):
        self.filter_variants()

    def get_active_search_terms(self):
        return [var.get().lower().strip() for var in self.search_vars if var.get().strip()]

    def get_mutation_count_filter(self):
        """Get the mutation count filter range"""
        min_str = self.min_mutations_var.get().strip()
        max_str = self.max_mutations_var.get().strip()
        
        min_val = None
        max_val = None
        
        try:
            if min_str:
                min_val = int(min_str)
        except ValueError:
            pass
            
        try:
            if max_str:
                max_val = int(max_str)
        except ValueError:
            pass
            
        return min_val, max_val

    def filter_variants(self):
        search_terms = self.get_active_search_terms()
        min_muts, max_muts = self.get_mutation_count_filter()
        
        if not search_terms and min_muts is None and max_muts is None:
            self.filtered_variants = self.all_variants.copy()
        else:
            self.filtered_variants = {}
            logic = self.search_logic_var.get()
            
            for var_id, mutations in self.all_variants.items():
                var_text = f"{var_id} {mutations}".lower()
                mutation_count = self.count_mutations(mutations)
                
                # Check text search criteria
                text_match = True
                if search_terms:
                    if logic == "AND":
                        text_match = all(term in var_text for term in search_terms)
                    else:  # OR logic
                        text_match = any(term in var_text for term in search_terms)
                
                # Check mutation count criteria
                count_match = True
                if min_muts is not None and mutation_count < min_muts:
                    count_match = False
                if max_muts is not None and mutation_count > max_muts:
                    count_match = False
                
                # Include variant if both criteria match
                if text_match and count_match:
                    self.filtered_variants[var_id] = mutations
        
        self.update_variant_display()

    def update_variant_display(self):
        # Clear existing widgets
        for widget in self.variant_listbox.winfo_children():
            widget.destroy()
        self.check_vars.clear()
        self.checkboxes.clear()

        # Update results label
        total = len(self.all_variants)
        filtered = len(self.filtered_variants)
        terms = self.get_active_search_terms()
        min_muts, max_muts = self.get_mutation_count_filter()
        
        filter_desc = []
        if terms:
            logic = self.search_logic_var.get()
            terms_txt = f" {logic} ".join([f'"{x}"' for x in terms])
            filter_desc.append(f"text: {terms_txt}")
        
        if min_muts is not None or max_muts is not None:
            if min_muts is not None and max_muts is not None:
                filter_desc.append(f"mutations: {min_muts}-{max_muts}")
            elif min_muts is not None:
                filter_desc.append(f"mutations: ≥{min_muts}")
            elif max_muts is not None:
                filter_desc.append(f"mutations: ≤{max_muts}")
        
        if not filter_desc:
            self.results_label.configure(text=f"Showing all {total} variants")
        else:
            filters = "; ".join(filter_desc)
            self.results_label.configure(text=f"Showing {filtered} of {total} matching [{filters}]")

        # Create checkboxes for filtered variants
        for var_id, muts in sorted(self.filtered_variants.items()):
            var_text = f"{var_id}: {muts if muts != '(none)' else 'no mutations'}"
            check_var = ctk.IntVar(value=0)
            cb = ctk.CTkCheckBox(self.variant_listbox, text=var_text, variable=check_var)
            cb.pack(fill="x", padx=10, pady=2, anchor="w")
            self.check_vars.append((var_id, check_var))
            self.checkboxes.append(cb)

    def select_all(self):
        for _, var_check_var in self.check_vars:
            var_check_var.set(1)

    def deselect_all(self):
        for _, var_check_var in self.check_vars:
            var_check_var.set(0)

    def refresh_variants(self):
        self.all_variants = self.controller.get_databank()
        self.filter_variants()

    def delete_variants(self):
        to_delete = [var_id for var_id, var_check_var in self.check_vars if var_check_var.get() == 1]
        if not to_delete:
            self.status_box.configure(text="No variants selected for deletion.", text_color="red")
            return

        databank = self.controller.get_databank()
        self.last_deleted = {var_id: databank[var_id] for var_id in to_delete if var_id in databank}

        for var_id in to_delete:
            databank.pop(var_id, None)

        self.controller.save_databank(databank)
        self.refresh_variants()
        self.status_box.configure(text=f"Deleted {len(to_delete)} variant(s). You can undo this action.", 
                                 text_color="green")

    def undo_delete(self):
        if not self.last_deleted:
            self.status_box.configure(text="Nothing to undo.", text_color="red")
            return

        databank = self.controller.get_databank()
        databank.update(self.last_deleted)
        self.controller.save_databank(databank)
        
        restored_count = len(self.last_deleted)
        self.last_deleted.clear()
        self.refresh_variants()
        self.status_box.configure(text=f"Restored {restored_count} variant(s).", text_color="green")

class ProtocolResultsPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)  # Results display area

        # Header frame
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Protocol Results",
                             corner_radius=10, fg_color="#4a90e2", text_color="white",
                             font=controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, sticky="e", padx=10, pady=10)

        # Control buttons frame
        controls_frame = ctk.CTkFrame(self)
        controls_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=5)
        
        refresh_btn = ctk.CTkButton(controls_frame, text="Refresh Results", command=self.refresh_results)
        refresh_btn.grid(row=0, column=0, sticky="w", padx=10, pady=10)

        # Results display area (expandable)
        self.results_text = ctk.CTkTextbox(self, font=controller.SMALLFONT)
        self.results_text.grid(row=2, column=0, sticky="nsew", padx=10, pady=10)

        # Load results on initialization
        self.refresh_results()

    def refresh_results(self):
        self.results_text.configure(state="normal")
        self.results_text.delete("1.0", "end")

        output_dir = Path(self.controller.output_dir.get())
        protocols_dir = output_dir / "protocols"
        protocol_file = protocols_dir / "protocol_data.json"

        if not protocol_file.exists():
            self.results_text.insert("end", "No protocol_data.json found. Run workflow first.\n")
            self.results_text.configure(state="disabled")
            return

        try:
            with open(protocol_file, "r") as f:
                protocol_data = json.load(f)
        except Exception as e:
            self.results_text.insert("end", f"Failed loading protocol_data.json:\n{e}\n")
            self.results_text.configure(state="disabled")
            return

        self.results_text.insert("end", "=== STEP-BY-STEP MUTAGENESIS PROTOCOL ===\n\n")

        # Step-by-step display
        protocol_steps = protocol_data.get("protocol_steps", {})
        if not protocol_steps:
            self.results_text.insert("end", "No protocol steps found in protocol_data.json.\n\n")
        else:
            # Final variants overview
            final_variants = protocol_data.get("final_variants_overview", {})
            if final_variants:
                self.results_text.insert("end", "FINAL VARIANTS OVERVIEW:\n")
                self.results_text.insert("end", "=" * 80 + "\n")
                for label, info in final_variants.items():
                    final_var = info.get("final_variant", "")
                    muts = ", ".join(info.get("mutations", []))
                    self.results_text.insert("end", f"{label:<15} -> {final_var:<15} | {muts}\n")
                self.results_text.insert("end", "\n")

            # Detailed protocol steps
            self.results_text.insert("end", "STEP-BY-STEP PROTOCOL:\n")
            self.results_text.insert("end", "=" * 80 + "\n")

            for step_num in sorted(map(int, protocol_steps.keys())):
                variants_in_step = protocol_steps.get(str(step_num), [])
                self.results_text.insert("end", f"\nSTEP {step_num}:\n")
                self.results_text.insert("end", "-" * 60 + "\n")
                for variant_info in variants_in_step:
                    new_var = variant_info.get("new_variant", "Unknown")
                    parent_var = variant_info.get("parent_variant", "Unknown")
                    muts_added = variant_info.get("mutations_added_str", "")
                    self.results_text.insert("end", f"Create: {new_var:<15} from {parent_var:<15}\n")
                    self.results_text.insert("end", f"       Add mutations: {muts_added}\n\n")

        # Display all variants summary
        all_variants = protocol_data.get("all_variants", {})
        self.results_text.insert("end", "\nALL VARIANTS IN PROTOCOL DATA:\n")
        self.results_text.insert("end", "=" * 80 + "\n")

        for variant, muts in sorted(all_variants.items()):
            self.results_text.insert("end", f"{variant:<15} | {muts}\n")

        self.results_text.insert("end", f"\nTotal variants: {len(all_variants)}\n")

        self.results_text.configure(state="disabled")

class PrimerPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)  # Table area

        # Header frame
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Primer List", 
                            corner_radius=10, fg_color="#4a90e2", 
                            text_color="white", font=controller.LARGEFONT,
                            height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, sticky="e", padx=10, pady=10)

        # Controls frame
        controls_frame = ctk.CTkFrame(self)
        controls_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=5)

        self.copy_btn = ctk.CTkButton(controls_frame, text="Copy Name + Sequence", command=self.copy_selected)
        self.copy_btn.grid(row=0, column=0, sticky="w", padx=10, pady=10)

        refresh_btn = ctk.CTkButton(controls_frame, text="Refresh", command=self.load_primers)
        refresh_btn.grid(row=0, column=1, sticky="w", padx=10, pady=10)

        # Scrollable frame for table (expandable)
        self.scroll_frame = ctk.CTkScrollableFrame(self)
        self.scroll_frame.grid(row=2, column=0, sticky="nsew", padx=10, pady=10)

        # Keep references to rows
        self.table_labels = []

    def load_primers(self):
        """Load primer_list.csv and display in table"""
        # Clear old table
        for widget in self.scroll_frame.winfo_children():
            widget.destroy()
        self.table_labels.clear()

        # Path to file
        path = Path(self.controller.output_dir.get()) / "primer_list.csv"
        if not path.exists():
            lbl = ctk.CTkLabel(self.scroll_frame, text="❌ No primer_list.csv found.", text_color="red")
            lbl.pack(pady=20)
            return

        try:
            # Read CSV
            with open(path, newline="", encoding='latin-1') as csvfile:
                reader = csv.DictReader(csvfile)
                primers = list(reader)

            if not primers:
                lbl = ctk.CTkLabel(self.scroll_frame, text="No primers found in file.", text_color="orange")
                lbl.pack(pady=20)
                return

            # Configure grid columns for proper spacing
            headers = ["Primer Name", "Primer Sequence", "Tm (°C)", "Overlap Tm (°C)"]
            col_weights = [23, 40, 12, 12]  # Roughly proportional to your original widths
            min_widths = [150, 250, 80, 80]   # Minimum widths to prevent too small columns

            for i, (weight, min_width) in enumerate(zip(col_weights, min_widths)):
                self.scroll_frame.grid_columnconfigure(i, weight=weight, minsize=min_width)

            # Table header
            for i, header in enumerate(headers):
                lbl = ctk.CTkLabel(self.scroll_frame, text=header, 
                                  font=ctk.CTkFont(family="Verdana", size=15, weight="bold"))
                lbl.grid(row=0, column=i, padx=5, pady=10, sticky="w")

            # Table rows
            for r, primer in enumerate(primers, start=1):
                row_labels = []
                cols = [
                    primer.get("Primer Name", ""),
                    primer.get("Primer Sequence", ""),
                    primer.get("Tm (°C)", ""),
                    primer.get("Overlap Tm", primer.get("Overlap Tm (°C)", ""))
                ]
                for c, value in enumerate(cols):
                    font = ctk.CTkFont(family="Verdana", size=13) if c == 1 else ctk.CTkFont(family="Verdana", size=13)
                    lbl = ctk.CTkLabel(self.scroll_frame, text=str(value), font=font)
                    lbl.grid(row=r, column=c, padx=5, pady=2, sticky="w")
                    row_labels.append(lbl)
                self.table_labels.append(row_labels)

        except Exception as e:
            error_lbl = ctk.CTkLabel(self.scroll_frame, text=f"Error loading primers: {str(e)}", text_color="red")
            error_lbl.pack(pady=20)

    def copy_selected(self):
        """Copy only first two columns (Name + Sequence) to clipboard"""
        if not self.table_labels:
            return
        
        lines = ["Primer Name\tPrimer Sequence"]  # Header
        for row in self.table_labels:
            name = row[0].cget("text")
            seq = row[1].cget("text")
            lines.append(f"{name}\t{seq}")
        
        text = "\n".join(lines)
        self.clipboard_clear()
        self.clipboard_append(text)

class GithubButton(ctk.CTkFrame):
    def __init__(self, master, url, text="GitHub", **kwargs):
        super().__init__(master, **kwargs)
        self.url = url
        self.button = ctk.CTkButton(self, text=text, command=self.open_link, fg_color="#222", text_color="white")
        self.button.pack(fill="x", padx=2, pady=2)

    def open_link(self):
        webbrowser.open(self.url)

if __name__ == "__main__":
    app = MutagenesisApp()
    app.mainloop()
