import customtkinter as ctk
from customtkinter import DrawEngine
DrawEngine.preferred_drawing_method = "polygon_shapes"
from customtkinter import CTkFont
from tkinter import simpledialog, filedialog, messagebox
from pathlib import Path
import json
import re
import copy
import os
import csv
from typing import List, Tuple, Dict
from collections import Counter, defaultdict
import primer3
import webbrowser
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm


# Constants and Shared Data
GENETIC_CODE = {
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

AMINO_ACID_MAP = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Xxx'
}

# Utility Functions
def translate_dna_to_protein(dna_sequence):
    """Translate DNA sequence to protein sequence"""
    protein = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        protein.append(GENETIC_CODE.get(codon, 'X'))
    return ''.join(protein)

def get_three_letter_code(one_letter):
    """Convert one-letter amino acid code to three-letter code"""
    return AMINO_ACID_MAP.get(one_letter, 'Xxx')

def parse_mutation(mut_str):
    """Parse mutation string like 'A123B' into (original_aa, position, new_aa)"""
    match = re.match(r'([A-Z])(\d+)([A-Z])', mut_str.strip())
    if not match:
        raise ValueError(f"Invalid mutation format: {mut_str}")
    return match.group(1), int(match.group(2)), match.group(3)

def mutation_position(mutation):
    """Extract position from mutation string for sorting"""
    matches = re.findall(r'\d+', mutation)
    return int(matches[0]) if matches else 0

class ScrollableFrameWithWheel(ctk.CTkScrollableFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        def on_mousewheel(event):
            self._parent_canvas.yview_scroll(int(-1 * (event.delta * 0.2)), "units")

        def on_mousewheel_linux_up(event):
            self._parent_canvas.yview_scroll(-1, "units")

        def on_mousewheel_linux_down(event):
            self._parent_canvas.yview_scroll(1, "units")

        self.bind("<Enter>", lambda e: [
            self.bind_all("<MouseWheel>", on_mousewheel),
            self.bind_all("<Button-4>", on_mousewheel_linux_up),
            self.bind_all("<Button-5>", on_mousewheel_linux_down),
        ])

        self.bind("<Leave>", lambda e: [
            self.unbind_all("<MouseWheel>"),
            self.unbind_all("<Button-4>"),
            self.unbind_all("<Button-5>"),
        ])

class ConfigManager:
    """Manage application configuration persistence"""
    
    def __init__(self, config_file="mutagenesis_config.json"):
        self.config_file = Path.home() / ".mutagenesis" / config_file
        self.config_file.parent.mkdir(parents=True, exist_ok=True)
        self.config = self.load_config()
    
    def load_config(self):
        """Load configuration from file"""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, Exception):
                return {}
        return {}
    
    def save_config(self):
        """Save configuration to file"""
        try:
            with open(self.config_file, 'w') as f:
                json.dump(self.config, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save config: {e}")
    
    def get(self, key, default=None):
        """Get configuration value"""
        return self.config.get(key, default)
    
    def set(self, key, value):
        """Set configuration value"""
        self.config[key] = value
        self.save_config()

class PrimerGenerator:
    """Site-directed mutagenesis primer generator with improved 3'-extension handling."""

    def __init__(self, flank_size_bases=30):
        self.stop_codons = {"TAA", "TAG", "TGA"}
        
        # Reverse codon table
        self.aa_to_codons = defaultdict(list)
        for codon, aa in GENETIC_CODE.items():
            self.aa_to_codons[aa].append(codon)

        # Constants
        self.flank_size_bases = flank_size_bases
        self.flank_size_codons = flank_size_bases // 3
        self.MIN_PRIMER_LEN = 30
        self.MIN_EXTENSION_LEN = 15
        self.TARGET_TM = 68.0
        self.MIN_TM = 55.0
        self.MAX_TM_DIFF = 3.0
        self.MAX_OVERLAP_LEN = 15

    def check_sequence_validity(self, dna_sequence, flank_size=30):
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

        print("Sequence check passed: correct ATG start and stop codon positions.")
        return True

    def validate_mutations_against_sequence(self, dna_sequence, variant_list):
        """Validate mutations against wildtype sequence."""
        for variant in variant_list:
            for mut in variant:
                original_aa, pos, _ = parse_mutation(mut)
                adj_pos_codons = pos + self.flank_size_codons
                codon_start = (adj_pos_codons - 1) * 3
                codon_seq = dna_sequence[codon_start:codon_start + 3]
                wt_aa = GENETIC_CODE.get(codon_seq, 'X')
                if wt_aa != original_aa:
                    raise ValueError(
                        f"Mutation {mut} does not match WT residue '{wt_aa}' "
                        f"at coding AA position {pos} (codon {codon_seq})"
                    )

    def reverse_complement(self, sequence):
        """Return reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement[base] for base in sequence[::-1])

    def calculate_tm(self, sequence):
        """Calculate melting temperature using primer3."""
        if not sequence:
            return 0.0
        try:
            return round(primer3.calc_tm(sequence), 2)
        except Exception as e:
            print(f"Error calculating Tm for sequence '{sequence}': {e}")
            return 0.0

    def get_best_codon(self, amino_acid, original_codon):
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

    def merge_close_mutations(self, mutations, max_distance_nt=21):
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

    def _create_mutated_sequence(self, dna_sequence, mutations):
        """Create mutated DNA sequence from list of mutations."""
        mutated_seq = list(dna_sequence)

        for aa, pos, new_aa in mutations:
            idx = ((pos - 1) + self.flank_size_codons) * 3
            original_codon = dna_sequence[idx:idx + 3]
            new_codon = self.get_best_codon(new_aa, original_codon)

            if GENETIC_CODE.get(original_codon, 'X') != aa:
                print(f"Warning: Original amino acid mismatch at position {pos}")

            mutated_seq[idx:idx + 3] = list(new_codon)

        return ''.join(mutated_seq)

    def _find_optimal_overlap(self, mutated_seq, codon_starts, codon_ends):
        # Implementation of overlap optimization logic
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

        # Stage 1: try exact MAX_OVERLAP_LEN
        best_start, best_end, best_tm, (f_flank, r_flank), _ = search_len(self.MAX_OVERLAP_LEN)
        if best_start is not None and f_flank >= min_inside_flank and r_flank >= min_inside_flank:
            return best_start, best_end, best_tm

        # Stage 2: allow +1 and +2 bases
        best_len_candidate = None
        best_len_score = -1e9
        for curr_len in [self.MAX_OVERLAP_LEN + 1, self.MAX_OVERLAP_LEN + 2]:
            cand_start, cand_end, cand_tm, _, cand_score = search_len(curr_len)
            if cand_start is not None and cand_score > best_len_score:
                best_len_candidate = (cand_start, cand_end, cand_tm)
                best_len_score = cand_score

        if best_len_candidate is not None:
            return best_len_candidate

        # Fallback: center between first/last mutation
        first_mut = min(mutation_positions)
        last_mut = max(mutation_positions)
        center = (first_mut + last_mut) // 2
        best_start = max(0, center - self.MAX_OVERLAP_LEN // 2)
        best_end = min(seq_len, best_start + self.MAX_OVERLAP_LEN)
        best_tm = self.calculate_tm(mutated_seq[best_start:best_end])
        return best_start, best_end, best_tm

    def _extension_mutation_span(self, mutation_positions, overlap_start, overlap_end):
        """Helper: compute span of mutations outside the overlap."""
        outside = [pos for pos in mutation_positions if pos < overlap_start or pos >= overlap_end]
        if not outside:
            return 0
        return max(outside) - min(outside)

    def _extend_primer_one_base(self, mutated_seq, current_primer, overlap_seq,
                              extension_start, extension_end, is_reverse):
        MAX_PRIMER_LEN = 45

        if len(current_primer) >= MAX_PRIMER_LEN:
            return current_primer, self.calculate_tm(current_primer), extension_start, extension_end

        if is_reverse:
            if extension_start == 0:
                return current_primer, self.calculate_tm(current_primer), extension_start, extension_end
            extension_start = max(0, extension_start - 1)
            extension_template = mutated_seq[extension_start:extension_end]
            extension = self.reverse_complement(extension_template)
            primer = self.reverse_complement(overlap_seq) + extension
        else:
            if extension_end >= len(mutated_seq):
                return current_primer, self.calculate_tm(current_primer), extension_start, extension_end
            extension_end = min(len(mutated_seq), extension_end + 1)
            extension = mutated_seq[extension_start:extension_end]
            primer = overlap_seq + extension

        tm = self.calculate_tm(primer)
        return primer, tm, extension_start, extension_end

    def count_nt_changes(self, wt_seq, mutations):
        """Return the count of nucleotide changes vs wildtype after applying mutations."""
        mutated_seq = self._create_mutated_sequence(wt_seq, mutations)
        return sum(1 for a, b in zip(wt_seq, mutated_seq) if a != b)

    def split_mutations_by_max_nt_changes(self, wt_seq, mutation_group, max_nt_changes=6):
        """Split mutation group so that each group does not exceed max_nt_changes."""
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

    def generate_primers_for_construct(self, dna_sequence, mutations):
        mutations = self.flatten_mutations(mutations)
        primer_sets = []
        mutation_groups = self.merge_close_mutations(mutations)

        MAX_PRIMER_LEN = 45
        MAX_MUT_NT = 6

        for group in mutation_groups:
            split_groups = self.split_mutations_by_max_nt_changes(dna_sequence, group, MAX_MUT_NT)
            previous_sub_group = []

            for i, sub_group in enumerate(split_groups):
                codon_starts_sub = [((m[1] - 1) + self.flank_size_codons) * 3 for m in sub_group]
                codon_ends_sub = [start + 3 for start in codon_starts_sub]

                all_applied_muts = previous_sub_group + sub_group
                mutated_seq = self._create_mutated_sequence(dna_sequence, all_applied_muts)

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

                # Balance primer Tm
                tm_diff = abs(forward_tm - reverse_tm)
                attempts = 0
                max_attempts = 30
                
                while ((tm_diff > self.MAX_TM_DIFF) or (forward_tm < self.TARGET_TM) or 
                       (reverse_tm < self.TARGET_TM)) and attempts < max_attempts:
                    
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

                all_mut_set = {(m[0], m[1], m[2]) for m in previous_sub_group + sub_group}
                all_mut_list = sorted(all_mut_set, key=lambda x: x[1])
                set_id = "_".join([f"{m[0]}{m[1]}{m[2]}" for m in all_mut_list])

                depends_on = [] if len(primer_sets) == 0 else [primer_sets[-1]['set_id']]

                primer_sets.append({
                    'set_id': set_id,
                    'depends_on': depends_on,
                    'mutations': [f"{m[0]}{m[1]}{m[2]}" for m in sub_group],
                    'all_covered_mutations': [f"{m[0]}{m[1]}{m[2]}" for m in all_mut_list],
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

    def save_primer_list(self, primer_data, output_dir=None, filename="primer_list.txt"):
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
                      f"{'Tm (C)':>9}{'OverlapLen':>11}{'OverlapTm':>10}{'Notes':>30}\n")
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
                "Primer Name", "Primer Sequence", "Length", "Tm (C)",
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

    def save_primer_data_json(self, primer_data, output_dir=None, filename="primer_list.json"):
        """Save primer data as JSON."""
        if output_dir is None:
            output_dir = os.path.expanduser("~/Mutagenesis")

        os.makedirs(output_dir, exist_ok=True)
        filepath_json = os.path.join(output_dir, filename)

        with open(filepath_json, 'w') as f:
            json.dump(primer_data, f, indent=2)

        print(f"Primer JSON saved to {filepath_json}")


class LabelPrintingSystem:
    """Complete label printing system for Herma 10900 label sheets"""
    def __init__(self):
        # Herma 10900 label specifications
        self.columns = 7
        self.rows = 27
        self.label_width = 25.3 * mm
        self.label_height = 9.9 * mm
        self.horizontal_spacing = 2.5 * mm
        self.top_margin = 13.5 * mm
        self.left_margin = 10.3 * mm
        self.font_size = 5.5

    def generate_labels_pdf(self, selected_variants, start_row, start_col, output_path):
        """Generate PDF with selected variants starting from specified position"""
        page_width, page_height = A4
        c = canvas.Canvas(str(output_path), pagesize=A4)
        c.setFont("Helvetica", self.font_size)

        # Sort variants for consistent ordering
        sorted_variants = sorted(selected_variants)
        
        variant_index = 0
        total_variants = len(sorted_variants)
        current_row = start_row
        current_col = start_col
        page_number = 1

        while variant_index < total_variants:
            # Calculate position
            x = self.left_margin + current_col * (self.label_width + self.horizontal_spacing)
            y = page_height - self.top_margin - (current_row + 1) * self.label_height

            # Draw label
            label_text = sorted_variants[variant_index]
            text_x = x + self.label_width / 2
            text_y = y + self.label_height / 2 - self.font_size / 2
            c.drawCentredString(text_x, text_y, label_text)
            
            variant_index += 1

            # Move to next position (left to right, top to bottom)
            current_col += 1
            if current_col >= self.columns:
                current_col = 0
                current_row += 1
                if current_row >= self.rows:
                    # Start new page
                    if variant_index < total_variants:
                        c.showPage()
                        c.setFont("Helvetica", self.font_size)
                        current_row = 0
                        current_col = 0
                        page_number += 1

        c.save()
        return True


class LabelSelectionWindow:
    """Interactive window for selecting starting position on label sheet"""
    def __init__(self, parent, selected_variants, output_dir, controller):
        self.parent = parent
        self.selected_variants = selected_variants
        self.output_dir = output_dir
        self.controller = controller
        self.label_system = LabelPrintingSystem()
        
        # Create window
        self.window = ctk.CTkToplevel(parent)
        self.window.title("Label Sheet Position Selection")
        self.window.geometry("800x600")
        self.window.transient(parent)
        self.window.after(100, self.window.grab_set)

        # Configure grid
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_rowconfigure(2, weight=1)

        # Selected cell tracking
        self.selected_cell = None
        self.cell_buttons = {}

        self.setup_ui()

    def setup_ui(self):
        # Title
        title_label = ctk.CTkLabel(self.window, 
                                  text=f"Select Starting Position for {len(self.selected_variants)} Labels",
                                  font=self.controller.LARGEFONT)
        title_label.grid(row=0, column=0, pady=20, sticky="ew")

        # Info label
        info_label = ctk.CTkLabel(self.window, 
                                 text="Click on a cell to select the starting position. Labels will fill left-to-right, top-to-bottom.",
                                 font=self.controller.MEDIUMFONT)
        info_label.grid(row=1, column=0, pady=10, sticky="ew")

        # Grid frame with proper configuration
        grid_frame = ScrollableFrameWithWheel(self.window, height=400)
        grid_frame.grid(row=2, column=0, sticky="nsew", padx=20, pady=10)
        
        # Configure grid frame columns to be uniform
        for col in range(self.label_system.columns):
            grid_frame.grid_columnconfigure(col, weight=1, uniform="col")

        # Create grid of buttons (7 columns x 27 rows)
        for row in range(self.label_system.rows):
            for col in range(self.label_system.columns):
                btn = ctk.CTkButton(grid_frame, 
                                   text=f"R{row+1}\nC{col+1}",
                                   width=60, height=40,
                                   font=ctk.CTkFont(size=11),
                                   command=lambda r=row, c=col: self.select_cell(r, c))
                btn.grid(row=row, column=col, padx=1, pady=1, sticky="nsew")
                self.cell_buttons[(row, col)] = btn

        # Buttons frame
        button_frame = ctk.CTkFrame(self.window)
        button_frame.grid(row=3, column=0, sticky="ew", padx=20, pady=10)

        # Status label
        self.status_label = ctk.CTkLabel(button_frame, text="Select a starting position", 
                                        font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=0, column=0, columnspan=3, pady=10)

        # Action buttons
        self.generate_btn = ctk.CTkButton(button_frame, text="Generate Labels PDF", 
                                         command=self.generate_pdf, state="disabled")
        self.generate_btn.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        cancel_btn = ctk.CTkButton(button_frame, text="Cancel", command=self.window.destroy)
        cancel_btn.grid(row=1, column=1, padx=10, pady=10, sticky="w")

        # Preview button
        self.preview_btn = ctk.CTkButton(button_frame, text="Preview Selection", 
                                   command=self.preview_selection, state="disabled")
        self.preview_btn.grid(row=1, column=2, padx=10, pady=10, sticky="w")

    def select_cell(self, row, col):
        """Handle cell selection"""
        # Reset previous selection
        if self.selected_cell:
            prev_row, prev_col = self.selected_cell
            self.cell_buttons[(prev_row, prev_col)].configure(fg_color=("gray75", "gray25"))

        # Highlight new selection
        self.selected_cell = (row, col)
        self.cell_buttons[(row, col)].configure(fg_color="green")
        
        # Update status and enable buttons
        self.status_label.configure(text=f"Selected: Row {row+1}, Column {col+1}")
        self.generate_btn.configure(state="normal")
        self.preview_btn.configure(state="normal")

    def preview_selection(self):
        """Preview which cells will be filled"""
        if not self.selected_cell:
            return

        # Reset all buttons to default color first
        for (row, col), btn in self.cell_buttons.items():
            if (row, col) != self.selected_cell:
                btn.configure(fg_color=("gray75", "gray25"))

        # Highlight cells that will be filled
        start_row, start_col = self.selected_cell
        current_row, current_col = start_row, start_col
        labels_placed = 0
        total_labels = len(self.selected_variants)

        while labels_placed < total_labels and current_row < self.label_system.rows:
            # Highlight current cell
            if (current_row, current_col) in self.cell_buttons:
                if (current_row, current_col) == self.selected_cell:
                    self.cell_buttons[(current_row, current_col)].configure(fg_color="green")
                else:
                    self.cell_buttons[(current_row, current_col)].configure(fg_color="blue")
            
            labels_placed += 1
            
            # Move to next position
            current_col += 1
            if current_col >= self.label_system.columns:
                current_col = 0
                current_row += 1

        # Update status
        if labels_placed < total_labels:
            self.status_label.configure(
                text=f"Preview: {labels_placed} labels fit on this page, {total_labels - labels_placed} continue on next page"
            )
        else:
            self.status_label.configure(text=f"Preview: All {total_labels} labels fit on this page")

    def generate_pdf(self):
        """Generate the labels PDF"""
        if not self.selected_cell:
            messagebox.showwarning("Warning", "Please select a starting position")
            return

        try:
            start_row, start_col = self.selected_cell
            output_path = Path(self.output_dir) / "protocols" / "selected_variant_labels.pdf"
            output_path.parent.mkdir(parents=True, exist_ok=True)

            success = self.label_system.generate_labels_pdf(
                self.selected_variants, start_row, start_col, output_path
            )

            if success:
                messagebox.showinfo("Success", f"Labels PDF generated successfully!\nSaved to: {output_path}")
                self.window.destroy()
            else:
                messagebox.showerror("Error", "Failed to generate labels PDF")

        except Exception as e:
            messagebox.showerror("Error", f"Error generating PDF: {str(e)}")


class MutagenesisProtocol:
    def __init__(self, variant_input, max_mutations_per_step, variant_prefix, output_dir_path,
                 start_col=2, start_row=16, list_existing_as_steps=False, undo_stack=None):
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
    def variant_key(mutations):
        return ','.join(sorted(mutations, key=mutation_position))

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
            sorted_mutations = sorted(variant_to_label_map[label], key=mutation_position)
            protocol_data["final_variants_overview"][label] = {
                "final_variant": variant_to_final_variant[label],
                "mutations": sorted_mutations
            }
        
        # Protocol steps
        for round_num in sorted(protocol_by_round):
            steps = sorted(protocol_by_round[round_num], key=lambda r: self.extract_variant_number(r[1]))
            step_data = []
            for row in steps:
                sorted_muts = sorted(row[2].split(', '), key=mutation_position)
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
        
        print(f"Protocol data saved: {protocol_json_path}")

    def generate_pdf(self, protocol_by_round, final_variants, filename, existing_input_variants,
                    variant_to_label_map, variant_to_final_variant, used_variants):
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
        final_table_data = [['Final Variant', 'Mutations']]
        for label in sorted(variant_to_label_map, key=self.extract_variant_number):
            sorted_mutations = sorted(variant_to_label_map[label], key=mutation_position)
            mut_paragraph = Paragraph(', '.join(sorted_mutations), cell_style)
            final_table_data.append([variant_to_final_variant[label], mut_paragraph])
        
        table = Table(final_table_data, colWidths=[150, 310], repeatRows=1)
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
                sorted_muts = ', '.join(sorted(row[2].split(', '), key=mutation_position))
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
            primer_groups[idx] = sorted(group_muts, key=mutation_position)
            for mut in group_muts:
                mutation_to_groups.setdefault(mut, []).append(idx)

        if self.list_existing_as_steps:
            variant_databank = {}
            existing_variants = {(): 'wildtype'}
        else:
            if self.databank_file.exists():
                try:
                    with open(self.databank_file, "r") as f:
                        content = f.read().strip()
                        variant_databank = json.loads(content) if content else {}
                except json.JSONDecodeError:
                    print("Warning: Databank file invalid. Starting fresh.")
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

        final_variants = defaultdict(list)
        variant_to_label_map = {}
        variant_to_final_variant = {}
        protocol_by_round = defaultdict(list)
        used_variants = {}
        existing_input_variants = {}

        existing_ids = [int(v_id[-2:]) for v_id in variant_databank if v_id.startswith(self.variant_prefix)]
        next_id_num = max(existing_ids, default=0) + 1
        for muts, name in existing_variants.items():
            mutation_str = ','.join(sorted(muts, key=mutation_position))
            if mutation_str not in variant_databank.values():
                new_id = f"{self.variant_prefix}{next_id_num:02d}"
                variant_databank[new_id] = mutation_str
                next_id_num += 1

        # If Skip Variant Databank is active, ignore existing variants during planning
        if self.list_existing_as_steps:
            planning_variants = {(): 'wildtype'}  # only start from wildtype
        else:
            planning_variants = dict(existing_variants)  # normal mode


        for idx, target_mutations in enumerate(variants):
            target = tuple(target_mutations)
            variant_key_str = self.variant_key(target)
            variant_label = f"{self.variant_prefix}{next_id_num:02d}"

            if variant_key_str in variant_databank.values():
                if not self.list_existing_as_steps:
                    # Just note it in the pre-existing list and skip real planning
                    for muts, name in planning_variants.items():
                        if variant_key_str == self.variant_key(muts):
                            existing_input_variants[variant_label] = name
                            break
                    continue
                else:
                    # Skip using databank variants — do nothing special
                    pass

            base_variant = self.get_base_variant(target, planning_variants)
            current_mutations = list(base_variant)
            base_name = planning_variants[base_variant]

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
                    key=mutation_position
                )
                test_sorted = tuple(test_mutations)

                if test_sorted not in planning_variants:
                    mutation_count = len(test_sorted)
                    new_name = f"{self.variant_prefix}{next_id_num:02d}"
                    planning_variants[test_sorted] = new_name
                    used_variants[test_sorted] = new_name
                    protocol_by_round[mutation_count].append([new_name, base_name, ', '.join(batch_muts)])
                    base_name = new_name
                    current_mutations = test_mutations
                    variant_databank[new_name] = ','.join(test_sorted)
                    next_id_num += 1

                for gid in batch_ids:
                    needed_group_ids.remove(gid)

            final_variants[base_name].append(target)
            variant_to_label_map[variant_label] = target
            variant_to_final_variant[variant_label] = base_name
            used_variants[tuple(sorted(map(self.normalize_mut, target), key=mutation_position))] = base_name

        # If skipping databank use, merge new variants into the existing databank
        if self.list_existing_as_steps:
            existing_variants.update({
                muts: name for muts, name in planning_variants.items() if muts not in existing_variants
            })


        self.generate_pdf(
            protocol_by_round, final_variants, str(self.pdf_path),
            existing_input_variants, variant_to_label_map,
            variant_to_final_variant, used_variants
        )

        self.save_protocol_json(
            protocol_by_round, final_variants, existing_input_variants,
            variant_to_label_map, variant_to_final_variant, used_variants
        )

        databank_file = self.output_dir / "variant_databank.json"
        if databank_file.exists():
            with open(databank_file, 'r') as f:
                try:
                    current_state = json.load(f)
                except json.JSONDecodeError:
                    current_state = {}
            self.undo_stack.append(copy.deepcopy(current_state))
        else:
            self.undo_stack.append({})

        if self.list_existing_as_steps:
            # Load existing databank
            if databank_file.exists():
                with open(databank_file, "r") as f:
                    try:
                        existing_json = json.load(f)
                    except json.JSONDecodeError:
                        existing_json = {}
            else:
                existing_json = {}

            # Merge new variants
            existing_json.update(variant_databank)

            # Save merged databank
            with open(databank_file, "w") as f:
                json.dump(existing_json, f, indent=2)
        else:
            # Normal behavior
            with open(databank_file, "w") as f:
                json.dump(variant_databank, f, indent=2)

        print(f"PDF generated: {self.pdf_path}")
        print(f"Databank updated: {self.databank_file}")

    def undo_change(self):
        databank_path = Path(self.controller.output_dir.get()) / "protocols" / "variant_databank.json"
        output_dir = Path(self.controller.output_dir.get())
        protocols_dir = output_dir / "protocols"
        
        try:
            if self.controller.undo_stack:
                last_state = self.controller.undo_stack.pop()
                
                with open(databank_path, "w") as f:
                    json.dump(last_state, f, indent=2)
                
                primer_files = [
                    output_dir / "primer_list.txt",
                    output_dir / "primer_list.csv", 
                    output_dir / "primer_list.json"
                ]
                
                for primer_file in primer_files:
                    if primer_file.exists():
                        primer_file.unlink()
                        print(f"Removed primer file: {primer_file}")
                
                protocol_files = [
                    protocols_dir / "mutagenesis_protocol.pdf",
                    protocols_dir / "protocol_data.json",
                    protocols_dir / "selected_variant_labels.pdf"
                ]
                
                for protocol_file in protocol_files:
                    if protocol_file.exists():
                        protocol_file.unlink()
                        print(f"Removed protocol file: {protocol_file}")
                
                self.status_label.configure(text="Last change undone successfully - all files cleaned up", text_color="green")
            else:
                self.status_label.configure(text="No changes to undo", text_color="orange")
        except Exception as e:
            self.status_label.configure(text=f"Undo failed: {str(e)}", text_color="red")


# Tkinter GUI Classes
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


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


class GithubButton(ctk.CTkFrame):
    def __init__(self, master, url, text="GitHub", **kwargs):
        super().__init__(master, **kwargs)
        self.url = url
        self.button = ctk.CTkButton(self, text=text, command=self.open_link, fg_color="#222", text_color="white")
        self.button.pack(fill="x", padx=2, pady=2)

    def open_link(self):
        webbrowser.open(self.url)


class BaseProteinPage:
    """Base class for protein-related pages to avoid code duplication"""
    
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
            three_letter = [get_three_letter_code(aa) for aa in sequence]
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

    def apply_mutations_to_protein(self, wt_protein, mutation_str):
        """Apply mutations to protein sequence"""
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
                    if protein_list[idx] != orig_aa:
                        print(f"Warning: WT AA at position {pos} is {protein_list[idx]}, expected {orig_aa}")
                    protein_list[idx] = new_aa
        return "".join(protein_list)


class MutagenesisApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Site-Directed Mutagenesis Primer Designer")
        self.geometry("1050x500")
        self.minsize(1200, 700)
        self.config_manager = ConfigManager()
        
        # Configure main window grid
        self.grid_columnconfigure(0, weight=0, minsize=200)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.LARGEFONT = ctk.CTkFont(family="Verdana", size=18, weight="bold")
        self.MEDIUMFONT = ctk.CTkFont(family="Verdana", size=12)
        self.SMALLFONT = ctk.CTkFont(family="Verdana", size=10)

        # Load last used output directory or use default
        last_output_dir = self.config_manager.get('last_output_dir', str(Path.cwd() / "Mutagenesis"))
        self.output_dir = ctk.StringVar(value=last_output_dir)
        
        # Trace changes to output_dir to save them
        self.output_dir.trace_add('write', lambda *args: self.save_output_dir())

        # Left navigation frame
        nav_frame = ctk.CTkFrame(self)
        nav_frame.grid(row=0, column=0, sticky="nsew", padx=(5,2), pady=5)
        nav_frame.grid_propagate(False)
        nav_frame.grid_columnconfigure(0, weight=1)

        # Navigation buttons
        nav_buttons = [
            ("User Input", InputPage),
            ("Variant Databank", DatabankPage),
            ("Wildtype Protein", WildtypeProteinPage),
            ("Variant Proteins", VariantProteinPage),
            ("Protocol Results", ProtocolResultsPage),
            ("Primer", PrimerPage),
            ("Mutation Extractor", MutationExtractorPage)
        ]

        num_nav_buttons = len(nav_buttons)
        nav_frame.grid_rowconfigure(num_nav_buttons, weight=1)

        for i, (text, page_class) in enumerate(nav_buttons):
            btn = ctk.CTkButton(nav_frame, text=text, command=lambda p=page_class: self.show_frame(p))
            btn.grid(row=i, column=0, sticky="ew", pady=(15 if i == 0 else 5, 5), padx=15)

        # GitHub button
        github_button = GithubButton(nav_frame, "https://github.com/xKosinus/Mutagenesis_Protocol/tree/windows_gui", text="GitHub Source")
        github_button.grid(row=num_nav_buttons + 1, column=0, sticky="ew", padx=15, pady=20)

        # Main container
        container = ctk.CTkFrame(self)
        container.grid(row=0, column=1, sticky="nsew", padx=(2,5), pady=5)
        container.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(0, weight=1)

        # Initialize frames
        self.frames = {}
        page_classes = [InputPage, MutationExtractorPage, DatabankPage, WildtypeProteinPage, 
                       VariantProteinPage, ProtocolResultsPage, PrimerPage]
        
        for F in page_classes:
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.undo_stack = []
        self.show_frame(InputPage)

    def show_frame(self, page_class):
        frame = self.frames[page_class]
        frame.tkraise()
        # Update specific pages when switching to them
        if page_class == WildtypeProteinPage:
            frame.update_protein_display()
        elif page_class == VariantProteinPage:
            frame.update_variant_display()
        elif page_class == PrimerPage:
            frame.load_primers()
        elif page_class == MutationExtractorPage:
            frame.refresh_reference_sequences()

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

    def save_output_dir(self):
        """Save output directory when it changes"""
        self.config_manager.set('last_output_dir', self.output_dir.get())

class InputPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid weights for responsive layout
        self.grid_columnconfigure(0, weight=0, minsize=200)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=0)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(4, weight=2)  
        self.grid_rowconfigure(5, weight=0)

        # Initialize sequence management variables
        self.sequence_undo_stack = []

        # Title label
        self.create_title_section()
        self.create_output_directory_section()
        self.create_dna_sequence_section()
        self.create_sequence_management_buttons()
        self.create_mutations_section()
        self.create_parameters_section()
        self.create_action_buttons()
        self.create_status_section()

    def create_title_section(self):
        title = ctk.CTkLabel(self, text="Mutagenesis Input", corner_radius=10,
                             fg_color="#4a90e2", text_color="white", 
                             font=self.controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=20, columnspan=3, sticky="ew")

    def create_output_directory_section(self):
        ctk.CTkLabel(self, text="Output Directory:", font=self.controller.MEDIUMFONT).grid(
            row=1, column=0, sticky="w", padx=10, pady=5)
        
        self.output_dir_entry = ctk.CTkEntry(self, textvariable=self.controller.output_dir)
        self.output_dir_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=5)
        
        select_dir_btn = ctk.CTkButton(self, text="Choose...", command=self.select_output_dir, width=70)
        select_dir_btn.grid(row=1, column=2, sticky="w", padx=5, pady=5)

    def create_dna_sequence_section(self):
        ctk.CTkLabel(self, text="Wildtype DNA Sequence:", font=self.controller.MEDIUMFONT).grid(
            row=2, column=0, sticky="nw", padx=10, pady=(15, 0))

        # Frame for DNA sequence input and info button
        seq_frame = ctk.CTkFrame(self)
        seq_frame.grid(row=2, column=1, columnspan=2, pady=5, padx=5, sticky="nsew")
        seq_frame.grid_columnconfigure(0, weight=1)
        seq_frame.grid_rowconfigure(0, weight=1)

        self.seq_text = PlaceholderTextbox(seq_frame, placeholder="Enter wildtype DNA sequence here...",
                                          placeholder_fg="gray", text_fg="white", height=150)
        self.seq_text.grid(row=0, column=0, sticky="nsew", padx=(0, 5))

        question_seq_btn = ctk.CTkButton(seq_frame, text="?", width=30, height=30, 
                                        command=self.show_info_sequence)
        question_seq_btn.grid(row=0, column=1, sticky="ne", pady=5)

    def create_sequence_management_buttons(self):
        seq_buttons_frame = ctk.CTkFrame(self)
        seq_buttons_frame.grid(row=3, column=1, columnspan=3, sticky="ew", padx=5, pady=5)
        
        buttons = [
            ("Save Sequence", self.save_sequence),
            ("Load/Manage Sequences", self.load_manage_sequences),
            ("Undo Sequence Changes", self.undo_sequence_changes)
        ]
        
        for i, (text, command) in enumerate(buttons):
            btn = ctk.CTkButton(seq_buttons_frame, text=text, command=command)
            btn.grid(row=0, column=i, padx=10, pady=5, sticky="w")

    def create_mutations_section(self):
        ctk.CTkLabel(self, text="Mutations per Variant:", font=self.controller.MEDIUMFONT).grid(
            row=4, column=0, sticky="nw", padx=10, pady=(10, 0))
        
        # Frame for variant text and info button
        var_frame = ctk.CTkFrame(self)
        var_frame.grid(row=4, column=1, columnspan=2, pady=5, padx=5, sticky="nsew")
        var_frame.grid_columnconfigure(0, weight=1)
        var_frame.grid_rowconfigure(0, weight=1)
        
        self.var_text = ctk.CTkTextbox(var_frame, height=150, font=self.controller.SMALLFONT)
        self.var_text.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        self.var_text.insert("1.0", "S19V, S30V, Q424V, A431E\nK3Q, K36R, D76Q, L167A, Q424V\n")

        question_btn = ctk.CTkButton(var_frame, text="?", width=30, height=30, 
                                    command=self.show_info_mutations)
        question_btn.grid(row=0, column=1, sticky="ne", pady=5)

    def create_parameters_section(self):
        params_frame = ctk.CTkFrame(self)
        params_frame.grid(row=5, column=0, columnspan=3, sticky="ew", padx=10, pady=10)
        params_frame.grid_columnconfigure(1, weight=1)
        params_frame.grid_columnconfigure(3, weight=1)

        # Max mutations per step
        ctk.CTkLabel(params_frame, text="Max mutations per step:", font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=5)
        
        self.max_mut_var = ctk.IntVar(value=1)
        max_mut_combo = ctk.CTkComboBox(params_frame, values=["1", "2"], variable=self.max_mut_var, width=80)
        max_mut_combo.grid(row=0, column=1, sticky="w", padx=5, pady=5)

        # Variant prefix
        ctk.CTkLabel(params_frame, text="Variant prefix:", font=self.controller.MEDIUMFONT).grid(
            row=0, column=2, sticky="w", padx=10, pady=5)
        
        self.var_prefix = ctk.StringVar(value="KWE_TA_A")
        prefix_entry = ctk.CTkEntry(params_frame, textvariable=self.var_prefix, width=150)
        prefix_entry.grid(row=0, column=3, sticky="w", padx=5, pady=5)

        # Skip databank checkbox
        self.skip_db = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(params_frame, text="Skip variant databank", variable=self.skip_db).grid(
            row=0, column=4, pady=10, padx=10, sticky="e")

    def create_action_buttons(self):
        button_frame = ctk.CTkFrame(self)
        button_frame.grid(row=6, column=0, columnspan=3, sticky="ew", padx=10, pady=10)
        
        run_btn = ctk.CTkButton(button_frame, text="Run Complete Workflow", 
                                font=self.controller.MEDIUMFONT, command=self.run_workflow)
        run_btn.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        
        undo_btn = ctk.CTkButton(button_frame, text="Undo Last Change", 
                                font=self.controller.MEDIUMFONT, command=self.undo_change)
        undo_btn.grid(row=0, column=1, padx=10, pady=10, sticky="w")

    def create_status_section(self):
        self.status_font = CTkFont(family="Verdana", size=14)
        self.status_label = ctk.CTkLabel(self, text="Ready to run workflow", 
                                        font=self.status_font, text_color="gray")
        self.status_label.grid(row=7, column=0, columnspan=3, sticky="w", padx=10,
                               pady=(0, 10))

    def get_dna_sequence(self):
        """Get the DNA sequence from the input field, cleaned and uppercase"""
        seq = self.seq_text.get("1.0", "end").strip()
        if seq == "Enter wildtype DNA sequence here...":
            return ""
        return seq.replace('\n', '').upper()

    def get_sequences_file_path(self):
        """Get the path to the wildtype_sequences.json file"""
        output_dir = Path(self.controller.output_dir.get())
        return output_dir / "wildtype_sequences.json"

    def load_sequences_database(self):
        file_path = self.get_sequences_file_path()
        if file_path.exists():
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    sequences = json.load(f)
                return sequences
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load sequences database: {str(e)}")
                return {}
        return {}

    def save_sequences_database(self, sequences_db):
        """Save sequences to JSON file"""
        file_path = self.get_sequences_file_path()
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create backup for undo functionality
        if file_path.exists():
            try:
                with open(file_path, "r") as f:
                    current_state = json.load(f)
                self.sequence_undo_stack.append(copy.deepcopy(current_state))
            except:
                self.sequence_undo_stack.append({})
        else:
            self.sequence_undo_stack.append({})
        
        # Keep only last 10 undo states to prevent memory issues
        if len(self.sequence_undo_stack) > 10:
            self.sequence_undo_stack.pop(0)
        
        with open(file_path, "w") as f:
            json.dump(sequences_db, f, indent=2)

    def save_sequence(self):
        """Save current sequence with user-provided name"""
        current_seq = self.get_dna_sequence()
        if not current_seq:
            messagebox.showwarning("Warning", "Please enter a DNA sequence before saving.")
            return

        # Get unique name from user
        while True:
            name = simpledialog.askstring("Save Sequence", "Enter a unique name for this sequence:")
            if name is None:  # User cancelled
                return
            
            name = name.strip()
            if not name:
                messagebox.showwarning("Warning", "Please enter a non-empty name.")
                continue
            
            # Check if name already exists
            sequences_db = self.load_sequences_database()
            if name in sequences_db:
                retry = messagebox.askyesno("Name Exists", 
                    f"A sequence named '{name}' already exists. Do you want to overwrite it?")
                if not retry:
                    continue
            
            # Save the sequence
            sequences_db[name] = current_seq
            self.save_sequences_database(sequences_db)
            self.status_label.configure(text=f"Sequence '{name}' saved successfully", text_color="green")
            break

    def load_manage_sequences(self):
        """Open sequence management window"""
        sequences_db = self.load_sequences_database()
        if not sequences_db:
            messagebox.showinfo("Info", "No saved sequences found.")
            return

        # Create management window
        manage_window = ctk.CTkToplevel(self)
        manage_window.title("Manage Saved Sequences")
        manage_window.geometry("800x500")
        manage_window.transient(self)
        manage_window.after(100, manage_window.grab_set)

        # Configure grid
        manage_window.grid_columnconfigure(0, weight=1)
        manage_window.grid_rowconfigure(1, weight=1)

        # Title
        title_label = ctk.CTkLabel(manage_window, text="Saved Sequences", 
                                  font=self.controller.LARGEFONT)
        title_label.grid(row=0, column=0, pady=20, sticky="ew")

        # Sequences list frame
        list_frame = ScrollableFrameWithWheel(manage_window)
        list_frame.grid(row=1, column=0, sticky="nsew", padx=20, pady=10)
        list_frame.grid_columnconfigure(1, weight=1)

        # Track checkboxes and radio buttons
        sequence_vars = {}
        load_var = ctk.StringVar(value="")

        # Headers
        headers = [("Load", 0), ("Delete", 1), ("Name", 2), ("Sequence (first 50 chars)", 3)]
        for text, col in headers:
            ctk.CTkLabel(list_frame, text=text, font=self.controller.MEDIUMFONT).grid(
                row=0, column=col, padx=10, pady=5, sticky="w")

        # Create list items
        for i, (name, sequence) in enumerate(sorted(sequences_db.items()), start=1):
            # Radio button for loading
            load_radio = ctk.CTkRadioButton(list_frame, text="", variable=load_var, value=name)
            load_radio.grid(row=i, column=0, padx=10, pady=2, sticky="w")
            
            # Checkbox for deletion
            delete_var = ctk.BooleanVar(value=False)
            sequence_vars[name] = delete_var
            delete_check = ctk.CTkCheckBox(list_frame, text="", variable=delete_var)
            delete_check.grid(row=i, column=1, padx=10, pady=2, sticky="w")
            
            # Name label
            name_label = ctk.CTkLabel(list_frame, text=name, font=self.controller.SMALLFONT)
            name_label.grid(row=i, column=2, padx=10, pady=2, sticky="w")
            
            # Sequence preview
            seq_preview = sequence[:50] + ("..." if len(sequence) > 50 else "")
            seq_label = ctk.CTkLabel(list_frame, text=seq_preview, 
                                   font=ctk.CTkFont(family="Courier", size=10))
            seq_label.grid(row=i, column=3, padx=10, pady=2, sticky="w")

        # Buttons frame
        button_frame = ctk.CTkFrame(manage_window)
        button_frame.grid(row=2, column=0, sticky="ew", padx=20, pady=10)

        def load_selected():
            selected_name = load_var.get()
            if not selected_name:
                messagebox.showwarning("Warning", "Please select a sequence to load.")
                return
            
            selected_sequence = sequences_db[selected_name]
            self.seq_text._clear_placeholder()
            self.seq_text.delete("1.0", "end")
            self.seq_text.insert("1.0", selected_sequence)
            self.status_label.configure(text=f"Loaded sequence '{selected_name}'", text_color="green")
            manage_window.destroy()

        def delete_selected():
            to_delete = [name for name, var in sequence_vars.items() if var.get()]
            if not to_delete:
                messagebox.showwarning("Warning", "Please select sequences to delete.")
                return
            
            confirm = messagebox.askyesno("Confirm Deletion", 
                f"Are you sure you want to delete {len(to_delete)} sequence(s)?")
            if not confirm:
                return
            
            # Remove selected sequences
            updated_db = {name: seq for name, seq in sequences_db.items() if name not in to_delete}
            self.save_sequences_database(updated_db)
            self.status_label.configure(text=f"Deleted {len(to_delete)} sequence(s)", text_color="orange")
            manage_window.destroy()

        buttons = [
            ("Load Selected", load_selected, 0),
            ("Delete Selected", delete_selected, 1),
            ("Cancel", manage_window.destroy, 2)
        ]
        
        for text, command, col in buttons:
            btn = ctk.CTkButton(button_frame, text=text, command=command)
            btn.grid(row=0, column=col, padx=10, pady=10, sticky="w")

    def undo_sequence_changes(self):
        """Undo last sequence database change"""
        if not self.sequence_undo_stack:
            self.status_label.configure(text="No sequence changes to undo", text_color="orange")
            return

        try:
            last_state = self.sequence_undo_stack.pop()
            file_path = self.get_sequences_file_path()
            
            with open(file_path, "w") as f:
                json.dump(last_state, f, indent=2)
            
            self.status_label.configure(text="Sequence changes undone successfully", text_color="green")
        except Exception as e:
            self.status_label.configure(text=f"Undo failed: {str(e)}", text_color="red")

    def show_info_mutations(self):
        messagebox.showinfo("Info",
                            "Enter mutations separated by commas, one variant per line.\nExample: S19V, S30V\nMultiple variants separated by new lines.")
        
    def show_info_sequence(self):
        messagebox.showinfo("Info", "Include 30 bases at the flanking sites!")

    def select_output_dir(self):
        new_dir = filedialog.askdirectory(title="Select Output Directory", 
                                         initialdir=self.controller.output_dir.get())
        if new_dir:
            self.controller.output_dir.set(new_dir)

    def run_workflow(self):
        self.status_label.configure(text="Running workflow...", text_color="blue")

        # Add validation early here
        if not self.validate_inputs():
            return
        
        seq = self.get_dna_sequence()
        variants_raw = self.var_text.get("1.0", "end").strip().split('\n')
        variant_list = [list(map(str.strip, line.split(','))) for line in variants_raw if line]
        max_mutations = self.max_mut_var.get()
        vp = self.var_prefix.get().strip()
        skip_db = self.skip_db.get()

        output_dir = Path(self.controller.output_dir.get())

        # Primer design
        try:
            pg = PrimerGenerator(flank_size_bases=30)
            if not pg.check_sequence_validity(seq, flank_size=30):
                self.status_label.configure(text="Sequence failed validation", text_color="red")
                return
            pg.validate_mutations_against_sequence(seq, variant_list)
        except Exception as e:
            self.status_label.configure(text=f"Error: {str(e)}", text_color="red")
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

        # Protocol generation
        try:
            protocol = MutagenesisProtocol(
                variant_list, max_mutations, vp, output_dir, 2, 16,
                list_existing_as_steps=skip_db, undo_stack=self.controller.undo_stack
            )
            protocol.run()
            self.status_label.configure(text="Workflow completed successfully! Check output directory.", 
                                       text_color="green")
            self.controller.frames[DatabankPage].refresh_variants()
        except Exception as e:
            self.status_label.configure(text=f"Error: {str(e)}", text_color="red")

    def undo_change(self):
        databank_path = Path(self.controller.output_dir.get()) / "protocols" / "variant_databank.json"
        output_dir = Path(self.controller.output_dir.get())
        protocols_dir = output_dir / "protocols"
        
        try:
            if self.controller.undo_stack:
                last_state = self.controller.undo_stack.pop()
                
                with open(databank_path, "w") as f:
                    json.dump(last_state, f, indent=2)
                
                primer_files = [
                    output_dir / "primer_list.txt",
                    output_dir / "primer_list.csv", 
                    output_dir / "primer_list.json"
                ]
                
                for primer_file in primer_files:
                    if primer_file.exists():
                        primer_file.unlink()
                        print(f"Removed primer file: {primer_file}")
                
                protocol_files = [
                    protocols_dir / "mutagenesis_protocol.pdf",
                    protocols_dir / "protocol_data.json",
                    protocols_dir / "selected_variant_labels.pdf"
                ]
                
                for protocol_file in protocol_files:
                    if protocol_file.exists():
                        protocol_file.unlink()
                        print(f"Removed protocol file: {protocol_file}")
                
                self.status_label.configure(text="Last change undone successfully - all files cleaned up", text_color="green")
            else:
                self.status_label.configure(text="No changes to undo", text_color="orange")
        except Exception as e:
            self.status_label.configure(text=f"Undo failed: {str(e)}", text_color="red")

    def validate_inputs(self):
        dna_seq = self.get_dna_sequence()
        mutations_raw = self.var_text.get("1.0", "end").strip()
        mutations_lines = [line.strip() for line in mutations_raw.splitlines() if line.strip()]
        mutations = []
        for line in mutations_lines:
            muts = [m.strip() for m in line.split(',') if m.strip()]
            mutations.append(muts)

        # DNA Sequence checks (same as before)...
        if not dna_seq:
            messagebox.showerror("Input Error", "DNA sequence is empty. Please enter a valid sequence.")
            return False
        if not re.fullmatch(r"[ATGC]+", dna_seq):
            messagebox.showerror("Input Error", "DNA sequence must only contain characters: A, T, G, C.")
            return False
        if len(dna_seq) % 3 != 0:
            messagebox.showerror("Input Error", "DNA sequence length must be divisible by 3.")
            return False
        if len(dna_seq) < 66:
            messagebox.showerror("Input Error", "DNA sequence is too short for required flanking regions (minimum 66 bases).")
            return False
        if dna_seq[30:33] != "ATG":
            messagebox.showerror("Input Error", f"Expected start codon ATG at bases 31-33, found '{dna_seq[30:33]}'.")
            return False
        stop_codon = dna_seq[-33:-30]
        if stop_codon not in {"TAA", "TAG", "TGA"}:
            messagebox.showerror("Input Error", f"Expected stop codon TAA, TAG or TGA before last 30 bases, found '{stop_codon}'.")
            return False

        # Mutation validations per variant
        if not mutations:
            messagebox.showerror("Input Error", "No mutations entered. Please enter at least one mutation.")
            return False

        mutation_pattern = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$")

        for idx, variant in enumerate(mutations):
            seen_in_variant = set()
            for mut in variant:
                if not mutation_pattern.match(mut):
                    messagebox.showerror("Input Error", f"Invalid mutation format: '{mut}' in line {idx + 1}. Expected format e.g. K49R.")
                    return False
                if mut in seen_in_variant:
                    messagebox.showerror("Input Error", f"Duplicate mutation '{mut}' found in variant line {idx + 1}.")
                    return False
                seen_in_variant.add(mut)

        # If all validations pass
        self.status_label.configure(text="Validation successful.", text_color="green")
        return True
    
    def on_output_dir_changed(self, *args):
        """Save output directory when it changes"""
        self.config_manager.set('last_output_dir', self.output_dir.get())

    def select_output_dir(self):
        # Get current directory (will be the last used one)
        current_dir = self.controller.output_dir.get()
        
        new_dir = filedialog.askdirectory(
            title="Select Output Directory", 
            initialdir=current_dir
        )
        
        if new_dir:
            self.controller.output_dir.set(new_dir)
            # The trace callback will automatically save it
            self.status_label.configure(
                text=f"Output directory changed to: {new_dir}", 
                text_color="green"
            )


class WildtypeProteinPage(ctk.CTkFrame, BaseProteinPage):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)

        self.create_header()
        self.create_info_section()
        self.create_format_controls()
        self.create_display_area()
        self.create_status_section()

        # Store the current protein sequence
        self.current_protein_1letter = ""

    def create_header(self):
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Wildtype Protein Sequence",
                             corner_radius=10, fg_color="#4a90e2", text_color="white",
                             font=self.controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: self.controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, padx=10, pady=10, sticky="e")

    def create_info_section(self):
        info_frame = ctk.CTkFrame(self)
        info_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        info_frame.grid_columnconfigure((0,1,2), weight=1)

        self.dna_length_label = ctk.CTkLabel(info_frame, text="DNA Length: N/A", 
                                           font=self.controller.MEDIUMFONT)
        self.dna_length_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        self.coding_length_label = ctk.CTkLabel(info_frame, text="Coding Length: N/A", 
                                              font=self.controller.MEDIUMFONT)
        self.coding_length_label.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        self.protein_length_label = ctk.CTkLabel(info_frame, text="Protein Length: N/A", 
                                               font=self.controller.MEDIUMFONT)
        self.protein_length_label.grid(row=0, column=2, padx=10, pady=10, sticky="w")

    def create_format_controls(self):
        format_frame = ctk.CTkFrame(self)
        format_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=5)

        ctk.CTkLabel(format_frame, text="Display Format:", font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, padx=10, pady=10, sticky="w")

        self.format_var = ctk.StringVar(value="1-letter")
        format_menu = ctk.CTkOptionMenu(format_frame, 
                                       values=["1-letter", "3-letter", "Both", "FASTA"],
                                       variable=self.format_var,
                                       command=self.update_display_format,
                                       width=120)
        format_menu.grid(row=0, column=1, padx=5, pady=10, sticky="w")

        refresh_btn = ctk.CTkButton(format_frame, text="Refresh", 
                                   command=self.update_protein_display, width=80)
        refresh_btn.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        copy_btn = ctk.CTkButton(format_frame, text="Copy to Clipboard", 
                                command=self.copy_to_clipboard, width=120)
        copy_btn.grid(row=0, column=3, padx=10, pady=10, sticky="w")

    def create_display_area(self):
        self.protein_text = ctk.CTkTextbox(self, font=ctk.CTkFont(family="Courier", size=11))
        self.protein_text.grid(row=3, column=0, pady=10, padx=10, sticky="nsew")
        self.protein_text.configure(state="disabled")

    def create_status_section(self):
        self.status_label = ctk.CTkLabel(self, text="", text_color="gray", 
                                        font=self.controller.SMALLFONT)
        self.status_label.grid(row=4, column=0, sticky="w", padx=10, pady=5)

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
        self.current_protein_1letter = translate_dna_to_protein(coding_seq)

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


class VariantProteinPage(ctk.CTkFrame, BaseProteinPage):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)

        self.create_header()
        self.create_info_section()
        self.create_controls()
        self.create_display_area()
        self.create_status_section()

        # Internal state
        self.current_variants = {}

    def create_header(self):
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)
        
        title = ctk.CTkLabel(header_frame, text="Variant Protein Sequences",
                             corner_radius=10, fg_color="#4a90e2", text_color="white",
                             font=self.controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        
        back_btn = ctk.CTkButton(header_frame, text="Back to Databank", 
                                 command=lambda: self.controller.show_frame(DatabankPage))
        back_btn.grid(row=0, column=1, padx=10, pady=10, sticky="e")

    def create_info_section(self):
        info_frame = ctk.CTkFrame(self)
        info_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        info_frame.grid_columnconfigure((0,1,2), weight=1)
        
        self.variant_count_label = ctk.CTkLabel(info_frame, text="Selected Variants: 0", 
                                              font=self.controller.MEDIUMFONT)
        self.variant_count_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        
        self.total_sequences_label = ctk.CTkLabel(info_frame, text="Total Sequences: 0", 
                                                font=self.controller.MEDIUMFONT)
        self.total_sequences_label.grid(row=0, column=1, padx=10, pady=10, sticky="w")
        
        self.avg_length_label = ctk.CTkLabel(info_frame, text="Avg Length: N/A", 
                                           font=self.controller.MEDIUMFONT)
        self.avg_length_label.grid(row=0, column=2, padx=10, pady=10, sticky="w")

    def create_controls(self):
        control_frame = ctk.CTkFrame(self)
        control_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=5)
        
        refresh_btn = ctk.CTkButton(control_frame, text="Refresh from Databank", 
                                   command=self.update_variant_display, width=150)
        refresh_btn.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        
        copy_btn = ctk.CTkButton(control_frame, text="Copy All to Clipboard", 
                                 command=self.copy_to_clipboard, width=150)
        copy_btn.grid(row=0, column=1, padx=10, pady=10, sticky="w")

    def create_display_area(self):
        self.protein_text = ctk.CTkTextbox(self, font=ctk.CTkFont(family="Courier", size=11))
        self.protein_text.grid(row=3, column=0, pady=10, padx=10, sticky="nsew")
        self.protein_text.configure(state="disabled")

    def create_status_section(self):
        self.status_label = ctk.CTkLabel(self, text="", text_color="gray", 
                                        font=self.controller.SMALLFONT)
        self.status_label.grid(row=4, column=0, sticky="w", padx=10, pady=5)

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
        wt_protein = translate_dna_to_protein(coding_dna)

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
        self.grid_columnconfigure(0, weight=1, minsize=400)
        self.grid_columnconfigure(1, weight=1, minsize=400)
        self.grid_rowconfigure(1, weight=1)

        # Internal state for variants and filtering
        self.check_vars = []
        self.checkboxes = []
        self.all_variants = {}
        self.filtered_variants = {}

        self.create_header()
        self.create_main_content()
        self.create_action_buttons()
        self.create_status_section()

        # Initialize with first search field
        self.add_search_field(is_first=True)
        self.refresh_variants()

    def create_header(self):
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Variant Databank",
                             corner_radius=15, fg_color="#4a90e2", text_color="white",
                             font=self.controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: self.controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, sticky="e", padx=10, pady=10)

    def create_main_content(self):
        content_frame = ctk.CTkFrame(self)
        content_frame.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=10, pady=5)
        content_frame.grid_columnconfigure(0, weight=1)
        content_frame.grid_columnconfigure(1, weight=1)
        content_frame.grid_rowconfigure(0, weight=1)

        # Variant list frame (left side)
        self.create_variant_list(content_frame)
        
        # Search section (right side)
        self.create_search_section(content_frame)

    def create_variant_list(self, parent):
        variant_frame = ctk.CTkFrame(parent)
        variant_frame.grid(row=0, column=0, sticky="nsew", padx=(5, 2), pady=5)
        variant_frame.grid_columnconfigure(0, weight=1)
        variant_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(variant_frame, text="Variants:", font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        self.variant_listbox = ScrollableFrameWithWheel(variant_frame)
        self.variant_listbox.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)

    def create_search_section(self, parent):
        search_main_frame = ctk.CTkFrame(parent)
        search_main_frame.grid(row=0, column=1, sticky="nsew", padx=(2, 5), pady=5)
        search_main_frame.grid_columnconfigure(0, weight=1)
        search_main_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(search_main_frame, text="Search & Filter:", font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        self.search_frame = ScrollableFrameWithWheel(search_main_frame)
        self.search_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        self.search_frame.grid_columnconfigure(0, weight=1)

        self.search_fields = []
        self.search_vars = []

        # Create search options
        self.create_search_options()

        # Results label
        self.results_label = ctk.CTkLabel(self, text="", text_color="gray")
        self.results_label.grid(row=2, column=0, columnspan=2, sticky="w", padx=10, pady=(5, 0))

    def create_search_options(self):
        options_frame = ctk.CTkFrame(self.search_frame)
        options_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        options_frame.grid_columnconfigure(5, weight=1)

        # Mutation count filter
        ctk.CTkLabel(options_frame, text="Mutations:", font=self.controller.SMALLFONT).grid(
            row=0, column=0, padx=5, pady=5, sticky="w")

        self.min_mutations_var = ctk.StringVar(value="")
        min_entry = ctk.CTkEntry(options_frame, textvariable=self.min_mutations_var, 
                               placeholder_text="Min", width=50)
        min_entry.grid(row=0, column=1, padx=2, pady=5, sticky="w")
        self.min_mutations_var.trace("w", self.on_search_change)

        ctk.CTkLabel(options_frame, text="to", font=self.controller.SMALLFONT).grid(
            row=0, column=2, padx=2, pady=5, sticky="w")

        self.max_mutations_var = ctk.StringVar(value="")
        max_entry = ctk.CTkEntry(options_frame, textvariable=self.max_mutations_var, 
                               placeholder_text="Max", width=50)
        max_entry.grid(row=0, column=3, padx=2, pady=5, sticky="w")
        self.max_mutations_var.trace("w", self.on_search_change)

        clear_mut_btn = ctk.CTkButton(options_frame, text="Clear", 
                                     command=self.clear_mutation_filter, width=60)
        clear_mut_btn.grid(row=0, column=4, padx=10, pady=5, sticky="w")

        # Logic controls
        ctk.CTkLabel(options_frame, text="Logic:", font=self.controller.SMALLFONT).grid(
            row=1, column=0, padx=5, pady=5, sticky="w")

        self.search_logic_var = ctk.StringVar(value="AND")
        logic_menu = ctk.CTkOptionMenu(options_frame, values=["AND", "OR"], 
                                      variable=self.search_logic_var,
                                      command=self.on_search_change, width=70)
        logic_menu.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.add_btn = ctk.CTkButton(options_frame, text="+", 
                                    command=self.add_search_field, width=30, height=30)
        self.add_btn.grid(row=1, column=2, padx=2, pady=5, sticky="w")

        self.remove_btn = ctk.CTkButton(options_frame, text="−", 
                                       command=self.remove_search_field,
                                       width=30, height=30, state="disabled")
        self.remove_btn.grid(row=1, column=3, padx=2, pady=5, sticky="w")

        clear_all_btn = ctk.CTkButton(options_frame, text="Clear All", 
                                     command=self.clear_all_search, width=80)
        clear_all_btn.grid(row=1, column=4, padx=10, pady=5, sticky="w")

        # Help text
        help_label = ctk.CTkLabel(options_frame, 
                                 text="Use AND to match all terms, OR for any. Filter by mutation count range (e.g., 2 to 4).",
                                 font=self.controller.SMALLFONT, text_color="gray")
        help_label.grid(row=2, column=0, columnspan=6, padx=5, pady=(0, 5), sticky="w")

    def create_action_buttons(self):
        button_frame = ctk.CTkFrame(self)
        button_frame.grid(row=3, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        
        # First row buttons
        first_row_buttons = [
            ("Select All", self.select_all, 0),
            ("Deselect All", self.deselect_all, 1),
            ("Delete Selected Variants", self.delete_variants, 2),
            ("Undo Delete", self.undo_delete, 3)
        ]
        
        # Second row buttons
        second_row_buttons = [
            ("Print Labels", self.print_selected_labels, 0),
            ("Refresh", self.refresh_variants, 1),
            ("View Proteins", lambda: self.controller.show_frame(VariantProteinPage), 2),
            ("Add Variants", self.add_variants_window, 3),
            ("Edit Variants", self.open_variant_editor, 4)
        ]
        
        # Create first row
        for text, command, col in first_row_buttons:
            btn = ctk.CTkButton(button_frame, text=text, command=command, width=160)
            btn.grid(row=0, column=col, sticky="w", padx=10, pady=5)
        
        # Create second row
        for text, command, col in second_row_buttons:
            btn = ctk.CTkButton(button_frame, text=text, command=command, width=160)
            btn.grid(row=1, column=col, sticky="w", padx=10, pady=5)

    def create_status_section(self):
        self.status_box = ctk.CTkLabel(self, text="", text_color="green")
        self.status_box.grid(row=4, column=0, columnspan=2, sticky="w", padx=10, pady=(0, 10))

    def add_search_field(self, is_first=False):
        field_num = len(self.search_fields)
        search_var = ctk.StringVar()
        search_var.trace("w", self.on_search_change)
        self.search_vars.append(search_var)
        
        field_frame = ctk.CTkFrame(self.search_frame)
        field_frame.grid(row=field_num + 2, column=0, sticky="ew", padx=5, pady=2)
        field_frame.grid_columnconfigure(1, weight=1)
        
        label_text = "Search:" if is_first else "Search:"
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

    def add_variants_window(self):
        """Open window for adding new variants to the databank"""
        output_dir = self.controller.output_dir.get()
        add_window = AddVariantsWindow(self, output_dir, self.controller)

    def open_variant_editor(self):
            """Open the variant editor window"""
            editor_window = VariantEditorWindow(self, self.controller)

    def print_selected_labels(self):
        """Print labels for selected variants"""
        selected_variants = []
        for var_id, check_var in self.check_vars:
            if check_var.get() == 1:
                selected_variants.append(var_id)
        
        if not selected_variants:
            self.status_box.configure(text="Please select variants to print labels for.", text_color="red")
            return
        
        # Open label selection window
        output_dir = self.controller.output_dir.get()
        label_window = LabelSelectionWindow(self, selected_variants, output_dir, self.controller)

class AddVariantsWindow:
    def __init__(self, parent, output_dir, controller):
        self.parent = parent
        self.output_dir = output_dir
        self.controller = controller
        
        # Create window
        self.window = ctk.CTkToplevel(parent)
        self.window.title("Add New Variants to Databank")
        self.window.geometry("900x600")
        self.window.transient(parent)
        
        # Configure main grid weights
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_rowconfigure(2, weight=1)

        # Track variant entries
        self.variant_entries = []

        # Setup UI first, then grab focus
        self.setup_ui()
        
        # Delay the grab_set to ensure window is fully rendered
        self.window.after(100, self.window.grab_set)

    def setup_ui(self):
        # Title
        title_label = ctk.CTkLabel(self.window, 
                                text="Add New Variants to Databank",
                                font=self.controller.LARGEFONT)
        title_label.grid(row=0, column=0, pady=20, padx=20, sticky="ew")

        # Prefix input frame
        prefix_frame = ctk.CTkFrame(self.window)
        prefix_frame.grid(row=1, column=0, sticky="ew", padx=20, pady=10)
        
        # Configure prefix frame grid
        prefix_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(prefix_frame, text="Variant Prefix:", 
                    font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=10)

        self.prefix_var = ctk.StringVar(value="KWE_TA_A")
        self.prefix_entry = ctk.CTkEntry(prefix_frame, textvariable=self.prefix_var, width=150)
        self.prefix_entry.grid(row=0, column=1, sticky="w", padx=10, pady=10)

        ctk.CTkLabel(prefix_frame, text="(Auto-generated: PREFIX01, PREFIX02, etc.)", 
                    font=self.controller.SMALLFONT, text_color="gray").grid(
            row=0, column=2, sticky="w", padx=10, pady=10)

        # Variants table frame
        table_frame = ctk.CTkFrame(self.window)
        table_frame.grid(row=2, column=0, sticky="nsew", padx=20, pady=10)
        
        # Configure table frame grid
        table_frame.grid_columnconfigure(0, weight=1)
        table_frame.grid_rowconfigure(1, weight=1)

        # Table header
        header_frame = ctk.CTkFrame(table_frame)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=(10, 5))
        
        # Configure header frame grid (single column now)
        header_frame.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(header_frame, text="Mutations (comma-separated)", 
                    font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, padx=10, pady=10, sticky="w")

        # Scrollable variants table
        self.variants_scroll_frame = ScrollableFrameWithWheel(table_frame)
        self.variants_scroll_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        
        # Configure scrollable frame grid (single column now)
        self.variants_scroll_frame.grid_columnconfigure(0, weight=1)

        # Initialize with some empty rows
        self.create_variant_rows(10)

        # Buttons frame
        button_frame = ctk.CTkFrame(self.window)
        button_frame.grid(row=3, column=0, sticky="ew", padx=20, pady=10)

        # Status label
        self.status_label = ctk.CTkLabel(button_frame, 
                                        text="Enter variant information and click Add to save", 
                                        font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=0, column=0, columnspan=4, pady=10, padx=10)

        # Action buttons
        add_btn = ctk.CTkButton(button_frame, text="Add Variants", 
                               command=self.add_variants, font=self.controller.MEDIUMFONT)
        add_btn.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        cancel_btn = ctk.CTkButton(button_frame, text="Cancel", command=self.window.destroy)
        cancel_btn.grid(row=1, column=1, padx=10, pady=10, sticky="w")

        # Utility buttons
        clear_btn = ctk.CTkButton(button_frame, text="Clear All", command=self.clear_all_entries)
        clear_btn.grid(row=1, column=2, padx=10, pady=10, sticky="w")

        add_rows_btn = ctk.CTkButton(button_frame, text="Add More Rows", command=self.add_more_rows)
        add_rows_btn.grid(row=1, column=3, padx=10, pady=10, sticky="w")

    def create_variant_rows(self, num_rows=10):
        """Create empty variant input rows"""
        for i in range(num_rows):
            row_num = len(self.variant_entries)
            
            # Only mutations entry (no name entry anymore)
            mutations_var = ctk.StringVar()
            mutations_entry = ctk.CTkEntry(self.variants_scroll_frame, textvariable=mutations_var,
                                        placeholder_text="e.g., S19V, Q424V, A431E")
            mutations_entry.grid(row=row_num, column=0, padx=5, pady=2, sticky="ew")
            
            # Store references (simplified structure)
            self.variant_entries.append({
                'mutations_var': mutations_var,
                'mutations_entry': mutations_entry
            })

    def clear_all_entries(self):
        """Clear all entry fields"""
        for entry_dict in self.variant_entries:
            entry_dict['mutations_var'].set("")
        self.status_label.configure(text="All entries cleared")

    def add_variants(self):
        """Add variants to the databank"""
        try:
            prefix = self.prefix_var.get().strip()
            if not prefix:
                messagebox.showwarning("Warning", "Please enter a variant prefix")
                return

            # Collect non-empty entries
            variants_to_add = []
            validation_errors = []
            
            for i, entry_dict in enumerate(self.variant_entries):
                mutations = entry_dict['mutations_var'].get().strip()
                
                if mutations:  # Only process rows with mutations
                    # Validate mutations format
                    is_valid, error_msg = self.validate_mutations(mutations)
                    if not is_valid:
                        validation_errors.append(f"Row {i+1}: {error_msg}")
                        continue
                    
                    variants_to_add.append({
                        'mutations': mutations,
                        'row': i+1
                    })

            if validation_errors:
                error_msg = "Validation errors found:\n\n" + "\n".join(validation_errors)
                messagebox.showerror("Validation Error", error_msg)
                return

            if not variants_to_add:
                messagebox.showwarning("Warning", "Please enter at least one variant with mutations")
                return

            # Load current databank and prepare backup for undo
            databank = self.parent.controller.get_databank()
            self.parent.controller.undo_stack.append(copy.deepcopy(databank))

            # Add variants with auto-numbering
            next_number = self.get_next_available_number()
            added_variants = []
            
            for variant_data in variants_to_add:
                # Generate variant ID automatically
                variant_id = f"{prefix}{next_number:02d}"
                next_number += 1
                
                # Check for duplicates
                if variant_id in databank:
                    overwrite = messagebox.askyesno("Duplicate Found", 
                        f"Variant '{variant_id}' already exists. Overwrite?")
                    if not overwrite:
                        continue
                
                # Add to databank
                databank[variant_id] = variant_data['mutations']
                added_variants.append(variant_id)

            if added_variants:
                # Save databank
                self.parent.controller.save_databank(databank)
                
                # Refresh parent databank display
                self.parent.refresh_variants()
                
                # Show success message
                success_msg = f"Successfully added {len(added_variants)} variant(s):\n\n"
                success_msg += "\n".join([f"• {vid}: {databank[vid]}" for vid in added_variants])
                messagebox.showinfo("Success", success_msg)
                
                self.status_label.configure(text=f"Added {len(added_variants)} variant(s) successfully")
                
                # Close window after successful addition
                self.window.destroy()
            else:
                self.status_label.configure(text="No variants were added")

        except Exception as e:
            messagebox.showerror("Error", f"Error adding variants: {str(e)}")
            self.status_label.configure(text=f"Error: {str(e)}")

class VariantEditorWindow:
    """Window for editing variant mutations while keeping labels unchanged"""
    def __init__(self, parent, controller):
        self.parent = parent
        self.controller = controller
        self.original_databank = {}
        self.modified_variants = {}
        
        # Create window
        self.window = ctk.CTkToplevel(parent)
        self.window.title("Edit Variant Mutations")
        self.window.geometry("1100x600")
        self.window.transient(parent)

        # Configure grid
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_columnconfigure(1, weight=1)
        self.window.grid_rowconfigure(1, weight=1)

        # Track entry widgets and search variables
        self.variant_entries = {}
        self.search_vars = []
        self.search_fields = []
        self.all_variants = {}
        self.filtered_variants = {}

        # Setup UI first
        self.setup_ui()
        
        # Then load variants
        self.load_variants()
        
        # Force update to ensure display
        self.window.update_idletasks()
        
        # Grab set after everything is ready
        self.window.after(100, self.window.grab_set)

    def setup_ui(self):
        # Title
        title_label = ctk.CTkLabel(self.window, 
                                  text="Edit Variant Mutations",
                                  font=self.controller.LARGEFONT)
        title_label.grid(row=0, column=0, columnspan=2, pady=20, padx=20, sticky="ew")

        # Left side - Variant list with editable mutations
        self.create_variant_list()
        
        # Right side - Search and filter
        self.create_search_section()
        
        # Bottom buttons
        self.create_action_buttons()

    def create_variant_list(self):
        variant_frame = ctk.CTkFrame(self.window)
        variant_frame.grid(row=1, column=0, sticky="nsew", padx=(20, 10), pady=10)
        variant_frame.grid_columnconfigure(0, weight=1)
        variant_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(variant_frame, text="Variants (Label: Mutations)", 
                    font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        # Scrollable frame for variants
        self.variant_scroll = ScrollableFrameWithWheel(variant_frame)
        self.variant_scroll.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        self.variant_scroll.grid_columnconfigure(0, weight=0)  # Label column
        self.variant_scroll.grid_columnconfigure(1, weight=1)  # Entry column

    def create_search_section(self):
        search_frame = ctk.CTkFrame(self.window)
        search_frame.grid(row=1, column=1, sticky="nsew", padx=(10, 20), pady=10)
        search_frame.grid_columnconfigure(0, weight=1)
        search_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(search_frame, text="Search & Filter:", 
                    font=self.controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        # Scrollable frame for search fields
        self.search_scroll = ScrollableFrameWithWheel(search_frame)
        self.search_scroll.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        self.search_scroll.grid_columnconfigure(0, weight=1)

        # Results label - CREATE THIS FIRST before search options
        self.results_label = ctk.CTkLabel(search_frame, text="", text_color="gray")
        self.results_label.grid(row=2, column=0, sticky="w", padx=10, pady=5)

        # Search options
        self.create_search_options()
        
        # Initial search field (this will trigger filter_variants which needs results_label)
        self.add_search_field(is_first=True)

    def create_search_options(self):
        options_frame = ctk.CTkFrame(self.search_scroll)
        options_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        options_frame.grid_columnconfigure(4, weight=1)

        # Logic controls
        ctk.CTkLabel(options_frame, text="Logic:", font=self.controller.SMALLFONT).grid(
            row=0, column=0, padx=5, pady=5, sticky="w")

        self.search_logic_var = ctk.StringVar(value="AND")
        logic_menu = ctk.CTkOptionMenu(options_frame, values=["AND", "OR"], 
                                      variable=self.search_logic_var,
                                      command=self.on_search_change, width=70)
        logic_menu.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        self.add_btn = ctk.CTkButton(options_frame, text="+", 
                                    command=self.add_search_field, width=30, height=30)
        self.add_btn.grid(row=0, column=2, padx=2, pady=5, sticky="w")

        self.remove_btn = ctk.CTkButton(options_frame, text="−", 
                                       command=self.remove_search_field,
                                       width=30, height=30, state="disabled")
        self.remove_btn.grid(row=0, column=3, padx=2, pady=5, sticky="w")

        clear_all_btn = ctk.CTkButton(options_frame, text="Clear All", 
                                     command=self.clear_all_search, width=80)
        clear_all_btn.grid(row=0, column=4, padx=10, pady=5, sticky="w")

    def add_search_field(self, is_first=False):
        field_num = len(self.search_fields)
        search_var = ctk.StringVar()
        search_var.trace("w", self.on_search_change)
        self.search_vars.append(search_var)
        
        field_frame = ctk.CTkFrame(self.search_scroll)
        field_frame.grid(row=field_num + 1, column=0, sticky="ew", padx=5, pady=2)
        field_frame.grid_columnconfigure(1, weight=1)
        
        label_text = "Search:" if is_first else "Search:"
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

    def on_search_change(self, *args):
        self.filter_variants()

    def get_active_search_terms(self):
        return [var.get().lower().strip() for var in self.search_vars if var.get().strip()]

    def filter_variants(self):
        search_terms = self.get_active_search_terms()
        
        if not search_terms:
            self.filtered_variants = self.all_variants.copy()
        else:
            self.filtered_variants = {}
            logic = self.search_logic_var.get()
            
            for var_id, mutations in self.all_variants.items():
                var_text = f"{var_id} {mutations}".lower()
                
                if logic == "AND":
                    if all(term in var_text for term in search_terms):
                        self.filtered_variants[var_id] = mutations
                else:  # OR logic
                    if any(term in var_text for term in search_terms):
                        self.filtered_variants[var_id] = mutations
        
        # Debug: print filtering results
        print(f"Filtered to {len(self.filtered_variants)} variants")
        
        self.update_variant_display()

    def update_variant_display(self):
        # Clear existing entries
        for widget in self.variant_scroll.winfo_children():
            widget.destroy()
        self.variant_entries.clear()

        # Update results label
        total = len(self.all_variants)
        filtered = len(self.filtered_variants)
        terms = self.get_active_search_terms()
        
        if not terms:
            self.results_label.configure(text=f"Showing all {total} variants")
        else:
            logic = self.search_logic_var.get()
            terms_txt = f" {logic} ".join([f'"{x}"' for x in terms])
            self.results_label.configure(text=f"Showing {filtered} of {total} matching [{terms_txt}]")

        # Debug: Check if we have variants to display
        print(f"Displaying {len(self.filtered_variants)} variants")
        
        # If no variants, show a message
        if not self.filtered_variants:
            no_data_label = ctk.CTkLabel(self.variant_scroll, 
                                        text="No variants found in databank.\nPlease run the workflow first to create variants.",
                                        font=self.controller.MEDIUMFONT,
                                        text_color="orange")
            no_data_label.grid(row=0, column=0, columnspan=2, pady=20, padx=20)
            return

        # Create rows for filtered variants
        for i, (var_id, mutations) in enumerate(sorted(self.filtered_variants.items())):
            # Label (read-only)
            label = ctk.CTkLabel(self.variant_scroll, text=f"{var_id}:", 
                               font=self.controller.MEDIUMFONT)
            label.grid(row=i, column=0, sticky="w", padx=(10, 5), pady=5)

            # Entry for mutations (editable)
            entry_var = ctk.StringVar(value=mutations)
            entry = ctk.CTkEntry(self.variant_scroll, textvariable=entry_var)
            entry.grid(row=i, column=1, sticky="ew", padx=(5, 10), pady=5)
            
            # Track changes
            entry_var.trace("w", lambda *args, vid=var_id: self.mark_modified(vid))
            
            self.variant_entries[var_id] = entry_var
            
        print(f"Created {len(self.variant_entries)} entry widgets")

    def mark_modified(self, var_id):
        """Track which variants have been modified"""
        if var_id in self.variant_entries:
            current_value = self.variant_entries[var_id].get()
            original_value = self.original_databank.get(var_id, "")
            
            if current_value != original_value:
                self.modified_variants[var_id] = current_value
            elif var_id in self.modified_variants:
                del self.modified_variants[var_id]

    def create_action_buttons(self):
        button_frame = ctk.CTkFrame(self.window)
        button_frame.grid(row=2, column=0, columnspan=2, sticky="ew", padx=20, pady=10)

        # Status label
        self.status_label = ctk.CTkLabel(button_frame, 
                                        text="Modify mutations as needed, then click Save Changes", 
                                        font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=0, column=0, columnspan=3, pady=10, padx=10)

        # Action buttons
        save_btn = ctk.CTkButton(button_frame, text="Save Changes", 
                                command=self.save_changes, font=self.controller.MEDIUMFONT)
        save_btn.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        cancel_btn = ctk.CTkButton(button_frame, text="Cancel", command=self.window.destroy)
        cancel_btn.grid(row=1, column=1, padx=10, pady=10, sticky="w")

        reset_btn = ctk.CTkButton(button_frame, text="Reset All Changes", 
                                 command=self.reset_changes)
        reset_btn.grid(row=1, column=2, padx=10, pady=10, sticky="w")

    def load_variants(self):
        """Load variants from databank"""
        self.all_variants = self.controller.get_databank()
        self.original_databank = copy.deepcopy(self.all_variants)
        
        # Debug: print what we loaded
        print(f"Loaded {len(self.all_variants)} variants from databank")
        
        self.filter_variants()

    def normalize_mutations(self, mutation_string):
        """Clean up mutation string - remove spaces and ensure consistency"""
        if not mutation_string or mutation_string == "(none)":
            return "(none)"
        
        # Split by comma, strip spaces from each mutation, filter empty strings
        mutations = [m.strip() for m in mutation_string.split(',') if m.strip()]
        
        # Join back without spaces after commas
        return ','.join(mutations)

    def validate_mutation_format(self, mutation_string):
        """Validate mutation format"""
        if mutation_string == "(none)":
            return True, ""
        
        mutations = [m.strip() for m in mutation_string.split(',') if m.strip()]
        
        if not mutations:
            return True, ""  # Empty is valid
        
        mutation_pattern = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$")
        
        for mut in mutations:
            if not mutation_pattern.match(mut):
                return False, f"Invalid mutation format: '{mut}'. Expected format: e.g., K49R"
        
        return True, ""

    def reset_changes(self):
        """Reset all changes to original values"""
        for var_id, entry_var in self.variant_entries.items():
            if var_id in self.original_databank:
                entry_var.set(self.original_databank[var_id])
        
        self.modified_variants.clear()
        self.status_label.configure(text="All changes reset to original values", text_color="orange")

    def save_changes(self):
        """Save modified variants to databank"""
        if not self.modified_variants:
            messagebox.showinfo("Info", "No changes to save")
            return

        # Validate and normalize all modified variants
        normalized_variants = {}
        validation_errors = []
        
        for var_id, mutation_string in self.modified_variants.items():
            # Normalize (remove spaces)
            normalized = self.normalize_mutations(mutation_string)
            
            # Validate format
            is_valid, error_msg = self.validate_mutation_format(normalized)
            if not is_valid:
                validation_errors.append(f"{var_id}: {error_msg}")
            else:
                normalized_variants[var_id] = normalized

        # Show validation errors if any
        if validation_errors:
            error_msg = "Validation errors found:\n\n" + "\n".join(validation_errors)
            messagebox.showerror("Validation Error", error_msg)
            return

        # Confirm changes
        change_summary = "\n".join([
            f"{vid}:\n  Old: {self.original_databank.get(vid, '')}\n  New: {new_val}"
            for vid, new_val in sorted(normalized_variants.items())
        ])
        
        confirm = messagebox.askyesno(
            "Confirm Changes",
            f"Save changes to {len(normalized_variants)} variant(s)?\n\n{change_summary}"
        )
        
        if not confirm:
            return

        try:
            # Backup for undo
            databank = self.controller.get_databank()
            self.controller.undo_stack.append(copy.deepcopy(databank))

            # Apply changes
            for var_id, new_mutations in normalized_variants.items():
                databank[var_id] = new_mutations

            # Save to file
            self.controller.save_databank(databank)
            
            # Refresh parent databank display
            self.parent.refresh_variants()
            
            messagebox.showinfo("Success", f"Successfully updated {len(normalized_variants)} variant(s)")
            self.window.destroy()
            
        except Exception as e:
            messagebox.showerror("Error", f"Error saving changes: {str(e)}")


class ProtocolResultsPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)

        self.create_header()
        self.create_controls()
        self.create_display_area()

        # Load results on initialization
        self.refresh_results()

    def create_header(self):
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Protocol Results",
                             corner_radius=10, fg_color="#4a90e2", text_color="white",
                             font=self.controller.LARGEFONT, height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: self.controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, sticky="e", padx=10, pady=10)

    def create_controls(self):
        controls_frame = ctk.CTkFrame(self)
        controls_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=5)

        refresh_btn = ctk.CTkButton(controls_frame, text="Refresh Results", command=self.refresh_results)
        refresh_btn.grid(row=0, column=0, sticky="w", padx=10, pady=10)

    def create_display_area(self):
        self.results_text = ctk.CTkTextbox(self, font=self.controller.SMALLFONT)
        self.results_text.grid(row=2, column=0, sticky="nsew", padx=10, pady=10)

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

class PrimerEditorWindow:
    """Window for editing individual primer sequences with +/- buttons"""
    
    def __init__(self, parent, controller, primer_data):
        self.parent = parent
        self.controller = controller
        self.primer_data = copy.deepcopy(primer_data)  # Work with a copy
        self.current_index = 0
        self.modified = False
        
        # Create window
        self.window = ctk.CTkToplevel(parent)
        self.window.title("Primer Editor")
        self.window.geometry("900x600")
        self.window.transient(parent)
        
        # Configure grid
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_rowconfigure(2, weight=1)
        
        # Get wildtype sequence for Tm calculation
        self.wt_sequence = self.get_wildtype_sequence()
        
        self.setup_ui()
        self.load_primer(0)
        
        # Grab focus after setup
        self.window.after(100, self.window.grab_set)
    
    def get_wildtype_sequence(self):
        """Get wildtype DNA sequence from input page"""
        try:
            input_page = self.controller.frames[InputPage]
            return input_page.get_dna_sequence()
        except:
            return ""
    
    def reverse_complement(self, sequence):
        """Return reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in sequence[::-1])
    
    def is_reverse_primer(self, primer_name):
        """Check if this is a reverse primer based on naming convention"""
        return '_rev' in primer_name.lower() or '_reverse' in primer_name.lower()
    
    def find_primer_position_in_wt(self, primer_seq, is_reverse=False):
        """Find where this primer sequence appears in the wildtype sequence"""
        if not self.wt_sequence or not primer_seq:
            return -1
        
        # For reverse primers, we need to work with reverse complement of WT
        search_sequence = self.wt_sequence if not is_reverse else self.reverse_complement(self.wt_sequence)
        
        # Try exact match
        pos = search_sequence.find(primer_seq)
        if pos != -1:
            return pos
        
        # Try finding a substring (primer might be modified)
        # Check overlapping regions with 80% match threshold
        for i in range(len(search_sequence) - len(primer_seq) + 1):
            matches = sum(1 for a, b in zip(search_sequence[i:i+len(primer_seq)], primer_seq) if a == b)
            if matches / len(primer_seq) > 0.8:  # 80% match
                return i
        
        return -1
    
    def get_next_base_from_wt(self, current_seq, direction='3prime', is_reverse=False):
        """Get the next base from wildtype sequence"""
        if not self.wt_sequence or not current_seq:
            return None
        
        # For reverse primers, use reverse complement of WT sequence
        search_sequence = self.wt_sequence if not is_reverse else self.reverse_complement(self.wt_sequence)
        
        pos = self.find_primer_position_in_wt(current_seq, is_reverse)
        if pos == -1:
            return None
        
        if direction == '5prime':
            # Add to 5' end - get base before current position
            if pos > 0:
                return search_sequence[pos - 1]
        else:  # 3prime
            # Add to 3' end - get base after current position
            end_pos = pos + len(current_seq)
            if end_pos < len(search_sequence):
                return search_sequence[end_pos]
        
        return None
    
    def calculate_tm(self, sequence):
        """Calculate Tm for a sequence"""
        if not sequence:
            return 0.0
        
        # Try using primer3 if available
        try:
            import primer3
            tm = primer3.calc_tm(sequence)
            if tm > 0:
                return round(tm, 2)
        except:
            pass
        
        # Fallback calculation
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        length = len(sequence)
        
        if length < 14:
            tm = at_count * 2 + gc_count * 4
        else:
            gc_percent = (gc_count / length) * 100 if length > 0 else 0
            tm = 81.5 + 0.41 * gc_percent - (675 / length) + 16.6 * (-1.3)
        
        return round(max(0.0, tm), 2)
    
    def setup_ui(self):
        # Title and navigation
        header_frame = ctk.CTkFrame(self.window)
        header_frame.grid(row=0, column=0, sticky="ew", padx=20, pady=20)
        header_frame.grid_columnconfigure(1, weight=1)
        
        ctk.CTkLabel(header_frame, text="Primer Editor", 
                    font=self.controller.LARGEFONT).grid(row=0, column=0, columnspan=3, pady=10)
        
        # Navigation buttons
        nav_frame = ctk.CTkFrame(header_frame)
        nav_frame.grid(row=1, column=0, columnspan=3, pady=10)
        
        self.prev_btn = ctk.CTkButton(nav_frame, text="◄ Previous", 
                                      command=self.prev_primer, width=100)
        self.prev_btn.grid(row=0, column=0, padx=5)
        
        self.primer_label = ctk.CTkLabel(nav_frame, text="Primer 1 of X", 
                                        font=self.controller.MEDIUMFONT)
        self.primer_label.grid(row=0, column=1, padx=20)
        
        self.next_btn = ctk.CTkButton(nav_frame, text="Next ►", 
                                      command=self.next_primer, width=100)
        self.next_btn.grid(row=0, column=2, padx=5)
        
        # Primer info
        info_frame = ctk.CTkFrame(self.window)
        info_frame.grid(row=1, column=0, sticky="ew", padx=20, pady=10)
        info_frame.grid_columnconfigure(5, weight=1)
        
        ctk.CTkLabel(info_frame, text="Primer Name:", 
                    font=self.controller.MEDIUMFONT).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.name_label = ctk.CTkLabel(info_frame, text="", 
                                      font=self.controller.MEDIUMFONT, text_color="#4a90e2")
        self.name_label.grid(row=0, column=1, sticky="w", padx=10, pady=5)
        
        ctk.CTkLabel(info_frame, text="Current Tm:", 
                    font=self.controller.MEDIUMFONT).grid(row=0, column=2, sticky="w", padx=10, pady=5)
        self.tm_label = ctk.CTkLabel(info_frame, text="0.0 C", 
                                     font=self.controller.MEDIUMFONT, text_color="green")
        self.tm_label.grid(row=0, column=3, sticky="w", padx=10, pady=5)
        
        ctk.CTkLabel(info_frame, text="Length:", 
                    font=self.controller.MEDIUMFONT).grid(row=0, column=4, sticky="w", padx=10, pady=5)
        self.length_label = ctk.CTkLabel(info_frame, text="0 bp", 
                                        font=self.controller.MEDIUMFONT)
        self.length_label.grid(row=0, column=5, sticky="w", padx=10, pady=5)
        
        # Sequence editor frame
        editor_frame = ctk.CTkFrame(self.window)
        editor_frame.grid(row=2, column=0, sticky="nsew", padx=20, pady=10)
        editor_frame.grid_columnconfigure(1, weight=1)
        editor_frame.grid_rowconfigure(0, weight=1)
        
        # Left controls (5' end)
        left_frame = ctk.CTkFrame(editor_frame)
        left_frame.grid(row=0, column=0, sticky="ns", padx=5, pady=5)
        
        ctk.CTkLabel(left_frame, text="5' End", font=self.controller.MEDIUMFONT).pack(pady=10)
        
        self.left_minus_btn = ctk.CTkButton(left_frame, text="− Remove", width=100, height=40,
                                           font=ctk.CTkFont(size=14),
                                           command=self.remove_left_base)
        self.left_minus_btn.pack(pady=10)
        
        # Show next base that will be added
        self.left_next_base_label = ctk.CTkLabel(left_frame, text="Next: ?", 
                                                 font=ctk.CTkFont(size=16, weight="bold"),
                                                 text_color="#4a90e2")
        self.left_next_base_label.pack(pady=10)
        
        self.left_plus_btn = ctk.CTkButton(left_frame, text="+ Add", width=100, height=40,
                                          font=ctk.CTkFont(size=14),
                                          command=self.add_left_base)
        self.left_plus_btn.pack(pady=10)
        
        # Sequence display (center)
        seq_frame = ctk.CTkFrame(editor_frame)
        seq_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
        seq_frame.grid_rowconfigure(1, weight=1)
        seq_frame.grid_columnconfigure(0, weight=1)
        
        ctk.CTkLabel(seq_frame, text="Primer Sequence", 
                    font=self.controller.MEDIUMFONT).grid(row=0, column=0, pady=10)
        
        self.sequence_text = ctk.CTkTextbox(seq_frame, font=ctk.CTkFont(family="Courier", size=14),
                                           wrap="char")
        self.sequence_text.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)
        
        # Right controls (3' end)
        right_frame = ctk.CTkFrame(editor_frame)
        right_frame.grid(row=0, column=2, sticky="ns", padx=5, pady=5)
        
        ctk.CTkLabel(right_frame, text="3' End", font=self.controller.MEDIUMFONT).pack(pady=10)
        
        self.right_minus_btn = ctk.CTkButton(right_frame, text="− Remove", width=100, height=40,
                                            font=ctk.CTkFont(size=14),
                                            command=self.remove_right_base)
        self.right_minus_btn.pack(pady=10)
        
        # Show next base that will be added
        self.right_next_base_label = ctk.CTkLabel(right_frame, text="Next: ?", 
                                                  font=ctk.CTkFont(size=16, weight="bold"),
                                                  text_color="#4a90e2")
        self.right_next_base_label.pack(pady=10)
        
        self.right_plus_btn = ctk.CTkButton(right_frame, text="+ Add", width=100, height=40,
                                           font=ctk.CTkFont(size=14),
                                           command=self.add_right_base)
        self.right_plus_btn.pack(pady=10)
        
        # Action buttons
        button_frame = ctk.CTkFrame(self.window)
        button_frame.grid(row=3, column=0, sticky="ew", padx=20, pady=10)
        
        self.status_label = ctk.CTkLabel(button_frame, text="Make changes and click Save", 
                                        font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=0, column=0, columnspan=4, pady=10)
        
        reset_btn = ctk.CTkButton(button_frame, text="Reset Current", 
                                 command=self.reset_current_primer)
        reset_btn.grid(row=1, column=0, padx=10, pady=10)
        
        save_btn = ctk.CTkButton(button_frame, text="Save All Changes", 
                                command=self.save_all_changes, fg_color="green")
        save_btn.grid(row=1, column=1, padx=10, pady=10)
        
        cancel_btn = ctk.CTkButton(button_frame, text="Cancel", command=self.window.destroy)
        cancel_btn.grid(row=1, column=2, padx=10, pady=10)
    
    def load_primer(self, index):
        """Load primer at given index"""
        if 0 <= index < len(self.primer_data):
            self.current_index = index
            primer = self.primer_data[index]
            
            # Check if this is a reverse primer
            primer_name = primer.get("Primer Name", "")
            self.current_is_reverse = self.is_reverse_primer(primer_name)
            
            # Update labels
            self.name_label.configure(text=primer_name)
            primer_type = " (Reverse)" if self.current_is_reverse else " (Forward)"
            self.primer_label.configure(text=f"Primer {index + 1} of {len(self.primer_data)}{primer_type}")
            
            # Update sequence display
            sequence = primer.get("Primer Sequence", "")
            self.sequence_text.delete("1.0", "end")
            
            # Format sequence with line breaks every 60 bases
            formatted = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
            self.sequence_text.insert("1.0", formatted)
            
            # Update Tm and length
            self.update_primer_info()
            
            # Update next base labels
            self.update_next_base_labels()
            
            # Update navigation buttons
            self.prev_btn.configure(state="normal" if index > 0 else "disabled")
            self.next_btn.configure(state="normal" if index < len(self.primer_data) - 1 else "disabled")
    
    def update_next_base_labels(self):
        """Update the labels showing which bases will be added next"""
        current_seq = self.get_current_sequence()
        
        # Use reverse flag for current primer
        is_reverse = getattr(self, 'current_is_reverse', False)
        
        # 5' end (left)
        next_base_5 = self.get_next_base_from_wt(current_seq, '5prime', is_reverse)
        if next_base_5:
            self.left_next_base_label.configure(text=f"Next: {next_base_5}")
            self.left_plus_btn.configure(state="normal")
        else:
            self.left_next_base_label.configure(text="Next: N/A")
            self.left_plus_btn.configure(state="disabled")
        
        # 3' end (right)
        next_base_3 = self.get_next_base_from_wt(current_seq, '3prime', is_reverse)
        if next_base_3:
            self.right_next_base_label.configure(text=f"Next: {next_base_3}")
            self.right_plus_btn.configure(state="normal")
        else:
            self.right_next_base_label.configure(text="Next: N/A")
            self.right_plus_btn.configure(state="disabled")
    
    def update_primer_info(self):
        """Update Tm and length display"""
        sequence = self.get_current_sequence()
        length = len(sequence)
        tm = self.calculate_tm(sequence)
        
        self.length_label.configure(text=f"{length} bp")
        self.tm_label.configure(text=f"{tm} C")
        
        # Update stored data
        self.primer_data[self.current_index]["Primer Sequence"] = sequence
        self.primer_data[self.current_index]["Length"] = str(length)
        self.primer_data[self.current_index]["Tm (C)"] = str(tm)
        
        # Update next base labels
        self.update_next_base_labels()
    
    def get_current_sequence(self):
        """Get current sequence from text widget"""
        text = self.sequence_text.get("1.0", "end-1c")
        # Remove whitespace and newlines
        return ''.join(text.split())
    
    def add_left_base(self):
        """Add a base to the 5' end from wildtype sequence"""
        current_seq = self.get_current_sequence()
        is_reverse = getattr(self, 'current_is_reverse', False)
        next_base = self.get_next_base_from_wt(current_seq, '5prime', is_reverse)
        
        if not next_base:
            self.status_label.configure(text="Cannot determine next base from wildtype sequence", text_color="red")
            return
        
        new_seq = next_base + current_seq
        
        self.sequence_text.delete("1.0", "end")
        formatted = '\n'.join([new_seq[i:i+60] for i in range(0, len(new_seq), 60)])
        self.sequence_text.insert("1.0", formatted)
        
        self.update_primer_info()
        self.modified = True
        primer_type = "reverse" if is_reverse else "forward"
        self.status_label.configure(text=f"Added {next_base} to 5' end ({primer_type} primer)", text_color="green")
    
    def remove_left_base(self):
        """Remove a base from the 5' end"""
        current_seq = self.get_current_sequence()
        if len(current_seq) > 1:
            removed = current_seq[0]
            new_seq = current_seq[1:]
            
            self.sequence_text.delete("1.0", "end")
            formatted = '\n'.join([new_seq[i:i+60] for i in range(0, len(new_seq), 60)])
            self.sequence_text.insert("1.0", formatted)
            
            self.update_primer_info()
            self.modified = True
            self.status_label.configure(text=f"Removed {removed} from 5' end", text_color="orange")
        else:
            self.status_label.configure(text="Cannot remove - sequence too short!", text_color="red")
    
    def add_right_base(self):
        """Add a base to the 3' end from wildtype sequence"""
        current_seq = self.get_current_sequence()
        is_reverse = getattr(self, 'current_is_reverse', False)
        next_base = self.get_next_base_from_wt(current_seq, '3prime', is_reverse)
        
        if not next_base:
            self.status_label.configure(text="Cannot determine next base from wildtype sequence", text_color="red")
            return
        
        new_seq = current_seq + next_base
        
        self.sequence_text.delete("1.0", "end")
        formatted = '\n'.join([new_seq[i:i+60] for i in range(0, len(new_seq), 60)])
        self.sequence_text.insert("1.0", formatted)
        
        self.update_primer_info()
        self.modified = True
        primer_type = "reverse" if is_reverse else "forward"
        self.status_label.configure(text=f"Added {next_base} to 3' end ({primer_type} primer)", text_color="green")
    
    def remove_right_base(self):
        """Remove a base from the 3' end"""
        current_seq = self.get_current_sequence()
        if len(current_seq) > 1:
            removed = current_seq[-1]
            new_seq = current_seq[:-1]
            
            self.sequence_text.delete("1.0", "end")
            formatted = '\n'.join([new_seq[i:i+60] for i in range(0, len(new_seq), 60)])
            self.sequence_text.insert("1.0", formatted)
            
            self.update_primer_info()
            self.modified = True
            self.status_label.configure(text=f"Removed {removed} from 3' end", text_color="orange")
        else:
            self.status_label.configure(text="Cannot remove - sequence too short!", text_color="red")
    
    def reset_current_primer(self):
        """Reset current primer to original sequence"""
        # Get original from parent
        original_data = self.parent.primer_data[self.current_index]
        self.primer_data[self.current_index] = copy.deepcopy(original_data)
        self.load_primer(self.current_index)
        self.status_label.configure(text="Primer reset to original", text_color="blue")
    
    def prev_primer(self):
        """Navigate to previous primer"""
        if self.current_index > 0:
            self.load_primer(self.current_index - 1)
    
    def next_primer(self):
        """Navigate to next primer"""
        if self.current_index < len(self.primer_data) - 1:
            self.load_primer(self.current_index + 1)
    
    def save_all_changes(self):
        """Save all modified primers back to CSV"""
        if not self.modified:
            messagebox.showinfo("No Changes", "No modifications were made.")
            return
        
        confirm = messagebox.askyesno("Confirm Save", 
                                      "Save all changes to primer_list.csv?\n\nThis will overwrite the existing file.")
        if not confirm:
            return
        
        try:
            # Save to CSV
            output_dir = Path(self.controller.output_dir.get())
            csv_path = output_dir / "primer_list.csv"
            
            # Backup original
            backup_path = output_dir / "primer_list.csv.backup"
            if csv_path.exists():
                import shutil
                shutil.copy(csv_path, backup_path)
            
            # Write updated data
            with open(csv_path, 'w', newline='', encoding='utf-8') as f:
                fieldnames = ["Primer Name", "Primer Sequence", "Length", "Tm (C)", 
                            "Overlap Length", "Overlap Tm", "Mutations"]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(self.primer_data)
            
            # Update parent's data
            self.parent.primer_data = copy.deepcopy(self.primer_data)
            self.parent.load_primers()
            
            messagebox.showinfo("Success", f"Primers saved successfully!\n\nBackup created: {backup_path.name}")
            self.window.destroy()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save primers:\n{str(e)}")

class PrimerPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)

        self.create_header()
        self.create_controls()
        self.create_table_area()

        # Keep references to rows and primer data
        self.table_labels = []
        self.primer_data = []
        
        # Cache for wildtype sequences
        self.wt_sequence = ""
        self.wt_sequence_rc = ""

    def create_header(self):
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame, text="Primer List", 
                            corner_radius=10, fg_color="#4a90e2", 
                            text_color="white", font=self.controller.LARGEFONT,
                            height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: self.controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, sticky="e", padx=10, pady=10)

    def create_controls(self):
        controls_frame = ctk.CTkFrame(self)
        controls_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=5)

        self.copy_btn = ctk.CTkButton(controls_frame, text="Copy Name + Sequence", command=self.copy_selected)
        self.copy_btn.grid(row=0, column=0, sticky="w", padx=10, pady=10)

        refresh_btn = ctk.CTkButton(controls_frame, text="Refresh", command=self.load_primers)
        refresh_btn.grid(row=0, column=1, sticky="w", padx=10, pady=10)

        edit_btn = ctk.CTkButton(controls_frame, text="Edit Primers", command=self.open_primer_editor)
        edit_btn.grid(row=0, column=2, sticky="w", padx=10, pady=10)

    def create_table_area(self):
        self.scroll_frame = ScrollableFrameWithWheel(self)
        self.scroll_frame.grid(row=2, column=0, sticky="nsew", padx=10, pady=10)

    def reverse_complement(self, sequence):
        """Return the reverse complement of a DNA sequence"""
        if not sequence:
            return ""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        rc = ''.join([complement.get(base, base) for base in sequence[::-1]])
        return rc.upper()

    def get_wildtype_sequence(self):
        """Get wildtype DNA sequence from input page"""
        try:
            input_page = self.controller.frames[InputPage]
            return input_page.get_dna_sequence()
        except:
            return ""

    def is_reverse_primer(self, primer_name):
        """Determine if primer is a reverse primer based on name"""
        name_lower = primer_name.lower()
        return ('_rev' in name_lower or 
                '_reverse' in name_lower or 
                name_lower.endswith('rev'))

    def get_reference_sequence(self, primer_name):
        """Get the appropriate reference sequence (forward or reverse complement)"""
        if self.is_reverse_primer(primer_name):
            return self.wt_sequence_rc
        else:
            return self.wt_sequence

    def find_primer_position(self, primer_seq, ref_seq):
        """Find where primer sequence appears in reference sequence"""
        if not ref_seq or not primer_seq:
            return -1
        
        # Try exact match first
        pos = ref_seq.find(primer_seq)
        if pos != -1:
            return pos
        
        # Try approximate match (75% threshold)
        best_match_pos = -1
        best_match_score = 0
        
        for i in range(len(ref_seq) - len(primer_seq) + 1):
            matches = sum(1 for a, b in zip(ref_seq[i:i+len(primer_seq)], primer_seq) if a == b)
            match_score = matches / len(primer_seq)
            
            if match_score > best_match_score and match_score > 0.75:
                best_match_score = match_score
                best_match_pos = i
        
        return best_match_pos if best_match_score > 0.75 else -1

    def get_mutated_positions(self, primer_name, primer_seq):
        """Find positions in primer that differ from wildtype sequence"""
        ref_seq = self.get_reference_sequence(primer_name)
        
        if not ref_seq or not primer_seq:
            return set()
        
        pos = self.find_primer_position(primer_seq, ref_seq)
        if pos == -1:
            return set()
        
        mutated_positions = set()
        wt_region = ref_seq[pos:pos + len(primer_seq)]
        
        for i, (primer_base, wt_base) in enumerate(zip(primer_seq, wt_region)):
            if primer_base != wt_base:
                mutated_positions.add(i)
        
        return mutated_positions

    def create_highlighted_textbox(self, parent, sequence, primer_name):
        """Create a textbox with mutation highlighting"""
        textbox = ctk.CTkTextbox(parent, height=30, width=400, 
                                font=ctk.CTkFont(family="Courier", size=12))
        textbox.grid_propagate(False)
        
        # Configure tag for red text
        textbox.tag_config("mutation", foreground="red")
        
        # Get mutated positions
        mutated_positions = self.get_mutated_positions(primer_name, sequence)
        
        # Insert sequence with highlighting
        for i, base in enumerate(sequence):
            if i in mutated_positions:
                textbox.insert("end", base, "mutation")
            else:
                textbox.insert("end", base)
        
        # Make read-only
        textbox.configure(state="disabled")
        
        return textbox

    def load_primers(self):
        """Load primer_list.csv and display in table with mutation highlighting"""
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
            # Load wildtype sequence for comparison
            self.wt_sequence = self.get_wildtype_sequence()
            if self.wt_sequence:
                self.wt_sequence_rc = self.reverse_complement(self.wt_sequence)
            else:
                # If no wildtype sequence, show warning but still display primers
                warning_lbl = ctk.CTkLabel(self.scroll_frame, 
                    text="⚠️ No wildtype sequence found - mutations won't be highlighted", 
                    text_color="orange")
                warning_lbl.grid(row=0, column=0, columnspan=4, pady=10, sticky="ew")

            # Read CSV
            with open(path, newline="", encoding='utf-8') as csvfile:
                reader = csv.DictReader(csvfile)
                primers = list(reader)

            if not primers:
                lbl = ctk.CTkLabel(self.scroll_frame, text="No primers found in file.", text_color="orange")
                lbl.pack(pady=20)
                return
            
            # Store primer data for editing
            self.primer_data = primers

            # Configure grid columns
            headers = ["Primer Name", "Primer Sequence", "Tm (C)", "Overlap Tm (C)"]
            
            # Make sequence column wider to accommodate textboxes
            self.scroll_frame.grid_columnconfigure(0, weight=0, minsize=180)  # Name
            self.scroll_frame.grid_columnconfigure(1, weight=1, minsize=420)  # Sequence (wider)
            self.scroll_frame.grid_columnconfigure(2, weight=0, minsize=80)   # Tm
            self.scroll_frame.grid_columnconfigure(3, weight=0, minsize=100)  # Overlap Tm

            # Table header row
            header_row = 0 if self.wt_sequence else 1
            
            for i, header in enumerate(headers):
                lbl = ctk.CTkLabel(self.scroll_frame, text=header, 
                                  font=ctk.CTkFont(family="Verdana", size=14, weight="bold"))
                lbl.grid(row=header_row, column=i, padx=5, pady=10, sticky="w")

            # Table rows
            start_row = header_row + 1
            for r, primer in enumerate(primers, start=start_row):
                row_widgets = []
                
                primer_name = primer.get("Primer Name", "")
                primer_seq = primer.get("Primer Sequence", "")
                
                # Column 0: Primer Name
                name_lbl = ctk.CTkLabel(self.scroll_frame, text=primer_name,
                                       font=ctk.CTkFont(family="Verdana", size=12))
                name_lbl.grid(row=r, column=0, padx=5, pady=2, sticky="w")
                row_widgets.append(name_lbl)
                
                # Column 1: Primer Sequence (with highlighting if wildtype available)
                if self.wt_sequence:
                    seq_textbox = self.create_highlighted_textbox(self.scroll_frame, primer_seq, primer_name)
                    seq_textbox.grid(row=r, column=1, padx=5, pady=2, sticky="ew")
                    row_widgets.append(seq_textbox)
                else:
                    # Fallback to regular label if no wildtype
                    seq_lbl = ctk.CTkLabel(self.scroll_frame, text=primer_seq,
                                          font=ctk.CTkFont(family="Courier", size=12))
                    seq_lbl.grid(row=r, column=1, padx=5, pady=2, sticky="w")
                    row_widgets.append(seq_lbl)
                
                # Column 2: Tm
                tm_lbl = ctk.CTkLabel(self.scroll_frame, 
                                     text=primer.get("Tm (C)", ""),
                                     font=ctk.CTkFont(family="Verdana", size=12))
                tm_lbl.grid(row=r, column=2, padx=5, pady=2, sticky="w")
                row_widgets.append(tm_lbl)
                
                # Column 3: Overlap Tm
                overlap_tm_lbl = ctk.CTkLabel(self.scroll_frame,
                                             text=primer.get("Overlap Tm", primer.get("Overlap Tm (C)", "")),
                                             font=ctk.CTkFont(family="Verdana", size=12))
                overlap_tm_lbl.grid(row=r, column=3, padx=5, pady=2, sticky="w")
                row_widgets.append(overlap_tm_lbl)
                
                self.table_labels.append(row_widgets)

        except Exception as e:
            error_lbl = ctk.CTkLabel(self.scroll_frame, text=f"Error loading primers: {str(e)}", text_color="red")
            error_lbl.pack(pady=20)

    def copy_selected(self):
        """Copy only first two columns (Name + Sequence) to clipboard"""
        if not self.table_labels:
            return
        
        lines = ["Primer Name\tPrimer Sequence"]  # Header
        for row in self.table_labels:
            # First widget is the name label
            name = row[0].cget("text")
            
            # Second widget is either a textbox or label
            seq_widget = row[1]
            if isinstance(seq_widget, ctk.CTkTextbox):
                # Get text from textbox (strip any formatting)
                seq = seq_widget.get("1.0", "end-1c").replace('\n', '')
            else:
                # Get text from label
                seq = seq_widget.cget("text")
            
            lines.append(f"{name}\t{seq}")
        
        text = "\n".join(lines)
        self.clipboard_clear()
        self.clipboard_append(text)

    def open_primer_editor(self):
        """Open the primer editor window"""
        if not self.primer_data:
            messagebox.showwarning("No Primers", "Please load primers first by clicking Refresh.")
            return
        
        PrimerEditorWindow(self, self.controller, self.primer_data)


class MutationExtractorPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Configure grid for responsive layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)  # FASTA content area

        # Header frame
        header_frame = ctk.CTkFrame(self)
        header_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        header_frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(header_frame,
                             text="Mutation Extractor",
                             corner_radius=10,
                             fg_color="#4a90e2",
                             text_color="white",
                             font=controller.LARGEFONT,
                             height=50)
        title.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        back_btn = ctk.CTkButton(header_frame, text="Back to Input", 
                                command=lambda: controller.show_frame(InputPage))
        back_btn.grid(row=0, column=1, padx=10, pady=10, sticky="e")

        # Reference sequence selection frame
        ref_frame = ctk.CTkFrame(self)
        ref_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        ref_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(ref_frame, text="Reference Sequence:", font=controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=10)

        self.ref_sequence_var = ctk.StringVar(value="")
        self.ref_dropdown = ctk.CTkOptionMenu(ref_frame, values=["No sequences available"], 
                                             variable=self.ref_sequence_var, width=200)
        self.ref_dropdown.grid(row=0, column=1, sticky="w", padx=5, pady=10)

        refresh_btn = ctk.CTkButton(ref_frame, text="Refresh", command=self.refresh_reference_sequences, width=80)
        refresh_btn.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        # FASTA file selection frame
        fasta_frame = ctk.CTkFrame(self)
        fasta_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=10)
        fasta_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(fasta_frame, text="FASTA File:", font=controller.MEDIUMFONT).grid(
            row=0, column=0, sticky="w", padx=10, pady=10)

        self.fasta_path_var = ctk.StringVar(value="")
        self.fasta_entry = ctk.CTkEntry(fasta_frame, textvariable=self.fasta_path_var, state="readonly")
        self.fasta_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=10)

        browse_btn = ctk.CTkButton(fasta_frame, text="Browse...", command=self.browse_fasta_file, width=80)
        browse_btn.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        # FASTA preview area (expandable)
        preview_label = ctk.CTkLabel(self, text="FASTA File Preview:", font=controller.MEDIUMFONT)
        preview_label.grid(row=3, column=0, sticky="nw", padx=10, pady=(10, 5))

        self.fasta_preview = ctk.CTkTextbox(self, font=ctk.CTkFont(family="Courier", size=11))
        self.fasta_preview.grid(row=3, column=0, sticky="nsew", padx=10, pady=(0, 10))
        self.fasta_preview.configure(state="disabled")

        # Action buttons frame
        action_frame = ctk.CTkFrame(self)
        action_frame.grid(row=4, column=0, sticky="ew", padx=10, pady=10)

        extract_btn = ctk.CTkButton(action_frame, text="Extract Mutations", 
                                   command=self.extract_mutations, font=controller.MEDIUMFONT)
        extract_btn.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        # Status label
        self.status_label = ctk.CTkLabel(self, text="Select reference sequence and FASTA file to begin", 
                                        font=controller.SMALLFONT, text_color="gray")
        self.status_label.grid(row=5, column=0, sticky="w", padx=10, pady=5)

        # Initialize
        self.refresh_reference_sequences()

    def get_sequences_file_path(self):
        """Get the path to the wildtype_sequences.json file"""
        output_dir = Path(self.controller.output_dir.get())
        return output_dir / "wildtype_sequences.json"

    def load_saved_sequences(self):
        """Load sequences from JSON file"""
        file_path = self.get_sequences_file_path()
        if file_path.exists():
            try:
                with open(file_path, "r") as f:
                    return json.load(f)
            except (json.JSONDecodeError, Exception):
                return {}
        return {}

    def refresh_reference_sequences(self):
        """Refresh the dropdown with available reference sequences"""
        sequences_db = self.load_saved_sequences()
        if sequences_db:
            sequence_names = list(sequences_db.keys())
            self.ref_dropdown.configure(values=sequence_names)
            if sequence_names:
                self.ref_sequence_var.set(sequence_names[0])
                self.status_label.configure(text=f"Loaded {len(sequence_names)} reference sequence(s)", text_color="green")
            else:
                self.ref_dropdown.configure(values=["No sequences available"])
                self.ref_sequence_var.set("")
        else:
            self.ref_dropdown.configure(values=["No sequences available"])
            self.ref_sequence_var.set("")
            self.status_label.configure(text="No saved sequences found. Save sequences in Input tab first.", text_color="orange")

    def browse_fasta_file(self):
        """Browse and select FASTA file"""
        from tkinter import filedialog
        
        file_path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fas *.fa"), ("All files", "*.*")]
        )
        
        if file_path:
            self.fasta_path_var.set(file_path)
            self.preview_fasta_file(file_path)

    def preview_fasta_file(self, file_path):
        """Preview the selected FASTA file"""
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Show first 1000 characters as preview
            preview_content = content[:1000]
            if len(content) > 1000:
                preview_content += "\n\n... (file truncated for preview)"
            
            self.fasta_preview.configure(state="normal")
            self.fasta_preview.delete("1.0", "end")
            self.fasta_preview.insert("1.0", preview_content)
            self.fasta_preview.configure(state="disabled")
            
            # Count sequences
            seq_count = content.count('>')
            self.status_label.configure(text=f"FASTA file loaded with {seq_count} sequence(s)", text_color="green")
            
        except Exception as e:
            self.status_label.configure(text=f"Error reading FASTA file: {str(e)}", text_color="red")

    def parse_fasta(self, file_path):
        """Parse FASTA file and return list of (header, sequence) tuples"""
        sequences = []
        try:
            with open(file_path, 'r') as f:
                current_header = None
                current_sequence = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_header is not None:
                            sequences.append((current_header, ''.join(current_sequence)))
                        # Start new sequence
                        current_header = line[1:]  # Remove '>'
                        current_sequence = []
                    elif line and current_header is not None:
                        current_sequence.append(line.upper())
                
                # Don't forget the last sequence
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_sequence)))
                    
        except Exception as e:
            raise Exception(f"Error parsing FASTA file: {str(e)}")
        
        return sequences

    def get_reference_protein_sequence(self):
        """Get the reference protein sequence from selected DNA sequence"""
        ref_name = self.ref_sequence_var.get()
        if not ref_name or ref_name == "No sequences available":
            raise Exception("Please select a reference sequence")
        
        sequences_db = self.load_saved_sequences()
        if ref_name not in sequences_db:
            raise Exception(f"Reference sequence '{ref_name}' not found")
        
        dna_sequence = sequences_db[ref_name]
        
        # Remove flanking regions (30bp from each end)
        if len(dna_sequence) <= 60:
            raise Exception("Reference DNA sequence is too short (must be >60bp with flanking regions)")
        
        coding_sequence = dna_sequence[30:-30]
        protein_sequence = translate_dna_to_protein(coding_sequence)
        
        return protein_sequence

    def simple_pairwise_alignment(self, seq1, seq2):
        """
        Simple pairwise alignment using dynamic programming (Needleman-Wunsch-like)
        Returns aligned sequences and alignment score
        """
        # Scoring parameters
        match_score = 2
        mismatch_penalty = -1
        gap_penalty = -2
        
        len1, len2 = len(seq1), len(seq2)
        
        # Initialize scoring matrix
        score_matrix = [[0 for _ in range(len2 + 1)] for _ in range(len1 + 1)]
        
        # Initialize first row and column with gap penalties
        for i in range(1, len1 + 1):
            score_matrix[i][0] = i * gap_penalty
        for j in range(1, len2 + 1):
            score_matrix[0][j] = j * gap_penalty
        
        # Fill the scoring matrix
        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                if seq1[i-1] == seq2[j-1]:
                    diagonal_score = score_matrix[i-1][j-1] + match_score
                else:
                    diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty
                
                up_score = score_matrix[i-1][j] + gap_penalty
                left_score = score_matrix[i][j-1] + gap_penalty
                
                score_matrix[i][j] = max(diagonal_score, up_score, left_score)
        
        # Traceback to get alignment
        aligned_seq1, aligned_seq2 = [], []
        i, j = len1, len2
        
        while i > 0 or j > 0:
            current_score = score_matrix[i][j] if i > 0 and j > 0 else float('-inf')
            up_score = score_matrix[i-1][j] if i > 0 else float('-inf')
            left_score = score_matrix[i][j-1] if j > 0 else float('-inf')
            diagonal_score = score_matrix[i-1][j-1] if i > 0 and j > 0 else float('-inf')
            
            if i > 0 and j > 0 and current_score == diagonal_score + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif i > 0 and current_score == up_score + gap_penalty:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            elif j > 0 and current_score == left_score + gap_penalty:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
            else:
                # Fallback - should not happen in a correct implementation
                if i > 0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')
                if j > 0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')
        
        # Reverse the sequences (traceback gives them backwards)
        aligned_seq1.reverse()
        aligned_seq2.reverse()
        
        return ''.join(aligned_seq1), ''.join(aligned_seq2), score_matrix[len1][len2]

    def find_mutations(self, reference_protein, variant_protein, variant_name):
        """Find mutations between reference and variant protein sequences using alignment"""
        try:
            # Perform pairwise alignment
            aligned_ref, aligned_var, alignment_score = self.simple_pairwise_alignment(reference_protein, variant_protein)
            
            # Calculate alignment quality (rough estimate)
            alignment_length = len(aligned_ref)
            if alignment_length == 0:
                raise Exception("Alignment failed - no aligned positions")
            
            # Count matches (excluding gaps)
            matches = sum(1 for a, b in zip(aligned_ref, aligned_var) if a == b and a != '-' and b != '-')
            aligned_positions = sum(1 for a, b in zip(aligned_ref, aligned_var) if a != '-' and b != '-')
            
            if aligned_positions == 0:
                raise Exception("No aligned positions found")
            
            identity = matches / aligned_positions
            if identity < 0.3:  # Less than 30% identity might indicate poor alignment
                raise Exception(f"Low sequence identity ({identity:.1%}) - possible alignment issue")
            
            # Extract mutations from aligned sequences
            mutations = []
            ref_position = 0  # Track position in original reference sequence (without gaps)
            
            for i, (ref_aa, var_aa) in enumerate(zip(aligned_ref, aligned_var)):
                if ref_aa != '-':  # Only count positions that exist in reference
                    ref_position += 1
                    
                    # Point mutation: both positions have amino acids and they're different
                    if var_aa != '-' and ref_aa != var_aa:
                        mutation = f"{ref_aa}{ref_position}{var_aa}"
                        mutations.append(mutation)
                    # Skip gaps and insertions in variant (we only want point mutations)
            
            return mutations, identity, alignment_length
            
        except Exception as e:
            raise Exception(f"Alignment error: {str(e)}")

    def extract_mutations(self):
        """Main function to extract mutations from FASTA sequences using alignment"""
        try:
            # Validate inputs
            if not self.fasta_path_var.get():
                raise Exception("Please select a FASTA file")
            
            # Get reference protein sequence
            self.status_label.configure(text="Processing reference sequence...", text_color="blue")
            reference_protein = self.get_reference_protein_sequence()
            
            # Parse FASTA file
            self.status_label.configure(text="Parsing FASTA file...", text_color="blue")
            fasta_sequences = self.parse_fasta(self.fasta_path_var.get())
            
            if not fasta_sequences:
                raise Exception("No sequences found in FASTA file")
            
            # Extract mutations for each sequence using alignment
            mutation_lines = []
            alignment_errors = []
            alignment_info = []
            
            self.status_label.configure(text=f"Aligning and extracting mutations from {len(fasta_sequences)} sequence(s)...", text_color="blue")
            
            for i, (header, sequence) in enumerate(fasta_sequences):
                try:
                    mutations, identity, alignment_length = self.find_mutations(reference_protein, sequence, header)
                    
                    if mutations:
                        mutation_line = ", ".join(mutations)
                        mutation_lines.append(mutation_line)
                    else:
                        # No mutations found - sequence is identical to reference in aligned regions
                        mutation_lines.append("")  # Empty line for identical sequences
                    
                    # Store alignment info for summary
                    alignment_info.append(f"{header}: {len(mutations)} mutations, {identity:.1%} identity, {alignment_length} aligned positions")
                        
                except Exception as seq_error:
                    alignment_errors.append(f"Sequence '{header}': {str(seq_error)}")
            
            # Show alignment errors if any
            if alignment_errors:
                error_msg = "The following sequences couldn't be aligned:\n\n" + "\n".join(alignment_errors)
                
                # Also show successful alignments for context
                if alignment_info:
                    error_msg += "\n\nSuccessful alignments:\n" + "\n".join(alignment_info)
                
                messagebox.showerror("Alignment Errors", error_msg)
                
                # Ask if user wants to continue with successfully processed sequences
                if mutation_lines:
                    continue_anyway = messagebox.askyesno(
                        "Continue?", 
                        f"Found {len(alignment_errors)} alignment error(s) and {len(mutation_lines)} successful extraction(s).\n\nDo you want to continue with the successfully extracted mutations?"
                    )
                    if not continue_anyway:
                        return
                else:
                    return
            
            # Show alignment summary for successful extractions
            if alignment_info and not alignment_errors:
                summary_msg = f"Extracted Sequences ({len(alignment_info)} sequences):\n"
                messagebox.showinfo("Extracted Sequences", summary_msg)
            
            # Transfer mutations to Input tab
            if mutation_lines:
                # Filter out empty lines (sequences with no mutations)
                non_empty_lines = [line for line in mutation_lines if line.strip()]
                
                if non_empty_lines:
                    input_page = self.controller.frames[InputPage]
                    input_page.var_text.delete("1.0", "end")  # Clear existing content
                    mutations_text = "\n".join(non_empty_lines)
                    input_page.var_text.insert("1.0", mutations_text)
                    
                    total_mutations = sum(len(line.split(", ")) for line in non_empty_lines if line)
                    self.status_label.configure(
                        text=f"Extracted {total_mutations} mutations from {len(non_empty_lines)} variant(s) and transferred to Input tab", 
                        text_color="green"
                    )
                    
                    # Switch to Input tab automatically
                    self.controller.show_frame(InputPage)
                else:
                    self.status_label.configure(text="No mutations found - all sequences are identical to reference in aligned regions", text_color="orange")
            else:
                self.status_label.configure(text="No mutations extracted", text_color="orange")
                
        except Exception as e:
            self.status_label.configure(text=f"Error: {str(e)}", text_color="red")


if __name__ == "__main__":
    app = MutagenesisApp()
    app.mainloop()
