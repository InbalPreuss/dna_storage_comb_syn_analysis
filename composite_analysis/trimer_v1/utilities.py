import csv
import math
import os
from pathlib import Path
from typing import Union, List, Dict
from itertools import product

import Levenshtein
from Bio import SeqIO, Align, pairwise2


# def trim_adapters_from_read(read: str, seq_adapter1: str, seq_adapter2: str, start_pos_bc: int, end_pos_bc: int) -> str:
#     """
#     Trims adapter sequences from the left and right sides of a read if they overlap
#     by at least 2 nucleotides.
#
#     :param read: The original sequence read.
#     :param seq_adapter1: The left adapter sequence.
#     :param seq_adapter2: The right adapter sequence.
#     :return: The trimmed sequence.
#     """
#
#     def find_overlap(sequence, adapter, from_end=True):
#         """Finds the longest overlap of at least 2 nucleotides between sequence and adapter."""
#         min_overlap = 2  # Minimum required overlap to trim
#         max_overlap = min(len(sequence), len(adapter))
#         temp_i = 0
#
#         for i in range(max_overlap, min_overlap - 1, -1):
#             if from_end:
#                 if sequence[-i:] == adapter[:i]:
#                     temp_i = i
#             else:
#                 if sequence[:i] == adapter[-i:]:
#                     temp_i = i
#
#         return temp_i
#
#     # Trim left adapter
#     left_trim = find_overlap(read, str(seq_adapter1), from_end=False)
#     if left_trim > 0:
#         read = read[left_trim:]
#
#     # Trim right adapter
#     right_trim = find_overlap(read, str(seq_adapter2), from_end=True)
#     if right_trim > 0:
#         read = read[:-right_trim]
#
#     return read

def trim_adapters_from_read(read: str, seq_adapter1: str, seq_adapter2: str, start_pos_bc: int, end_pos_bc: int) -> str:
    """
    Trims adapter sequences from the left and right sides of a read if they overlap
    by at least 2 nucleotides, ensuring that trimming does not affect the region
    between start_pos_bc and end_pos_bc.

    :param read: The original sequence read.
    :param seq_adapter1: The left adapter sequence.
    :param seq_adapter2: The right adapter sequence.
    :param start_pos_bc: Start position of the barcode, which should not be trimmed.
    :param end_pos_bc: End position of the barcode, which should not be trimmed.
    :return: The trimmed sequence.
    """

    def find_overlap(sequence, adapter, from_end=True, limit_start=0, limit_end=None):
        """Finds the longest overlap of at least 2 nucleotides between sequence and adapter within allowed boundaries."""
        min_overlap = 2  # Minimum required overlap to trim
        max_overlap = min(len(sequence), len(adapter))
        temp_i = 0

        for i in range(max_overlap, min_overlap - 1, -1):
            if from_end:
                if limit_start <= len(sequence) - i and sequence[-i:] == adapter[:i]:
                    temp_i = i
            else:
                if limit_end is None or i <= limit_end:
                    if sequence[:i] == adapter[-i:]:
                        temp_i = i

        return temp_i

    # Trim left adapter only before start_pos_bc
    left_trim = find_overlap(read[:start_pos_bc], seq_adapter1, from_end=False)
    if left_trim > 0:
        read = read[left_trim:]
        start_pos_bc -= left_trim
        end_pos_bc -= left_trim

    # Trim right adapter only after end_pos_bc
    right_trim = find_overlap(read[end_pos_bc:], seq_adapter2, from_end=True)
    if right_trim > 0:
        read = read[:-right_trim]

    return read


def is_within_levenshtein_distance(bc_seq: str, seq_design:str, read: str, max_distance: int,
                                   seq_adapter1: str, seq_adapter2: str, is_trim=True) -> tuple[int, int, int, List[str], str]:
    """
        Calculates Levenshtein distance between the barcode and a read.
        Returns distance, start positions, and annotations.
    """

    distance, start_pos_bc_min_dist, end_pos_bc_min_dist, start_pos_read = math.inf, math.inf, math.inf, math.inf
    annotations = []
    bc_seq_len = len(bc_seq)
    seq_design_len = len(seq_design)
    read_len = len(read)

    trimmed_read = read

    for start_pos in range(read_len - bc_seq_len + 1):
        end_pos = start_pos + bc_seq_len
        sub_read = read[start_pos:end_pos]
        distance = Levenshtein.distance(bc_seq, sub_read)
        if distance <= max_distance:
            max_distance = distance
            start_pos_bc_min_dist = start_pos
            end_pos_bc_min_dist = end_pos + 1

            # If is trimmed, then trim the adapters got in the design from the read
            if is_trim:
                trimmed_read = trim_adapters_from_read(read, seq_adapter1, seq_adapter2, start_pos_bc_min_dist, end_pos_bc_min_dist)

            # Compute the edit annotations of read compared to seq_design
            # start_pos_read = max(0, end_pos - seq_design_len)
            # annotations = compute_edit_annotations(seq_design, read[:end_pos])
            annotations = compute_edit_annotations(seq_design, trimmed_read[start_pos:start_pos+seq_design_len])


    # return max_distance, start_pos_bc_min_dist, end_pos_bc_min_dist, annotations, start_pos_read
    return max_distance, start_pos_bc_min_dist, end_pos_bc_min_dist, annotations, trimmed_read


def get_amount_of_reads_from_file(file_path: Union[Path, str]):
    # Convert the string to a Path object
    file_path = Path(file_path)
    # Determine the file format based on the file extension
    if file_path.suffix == '.fastq':
        file_format = 'fastq'
    elif file_path.suffix == '.fasta':
        file_format = 'fasta'
    else:
        raise ValueError("Unsupported file format. The file must be either .fastq or .fasta")

    # Open the file and create an index
    index = SeqIO.index(str(file_path), file_format)

    # Get the number of sequences in the index
    num_sequences = len(index)
    return num_sequences

def compute_edit_annotations(ref_seq, comp_seq):
    """
    Returns a list where each position in the reference sequence is annotated with:
    - '0'  (if the nucleotide is identical)
    - '1S' (if there is a substitution)
    - '1D' (if a deletion occurs at that position)
    - '1I' (if an insertion occurs at that position)
    """
    ref_len = len(ref_seq)
    comp_len = len(comp_seq)

    # Compute Levenshtein edit operations
    edit_operations = Levenshtein.editops(ref_seq, comp_seq)

    # Initialize annotation list
    annotations = ["0"] * ref_len  # Default all positions to "0"

    for op, ref_pos, comp_pos in edit_operations:
        if op == "replace":
            annotations[ref_pos] = f"S"  # Substitution
        elif op == "delete":
            annotations[ref_pos] = f"D"  # Deletion
        elif op == "insert":
            if ref_pos < ref_len:  # Handle insertions inside the sequence
                annotations.insert(ref_pos, f"I")  # Insertion

    return annotations[:ref_len]  # Ensure length matches ref sequence


def write_list_to_csv(data: List, file_name: Union[Path, str]) -> None:
    with open(file_name, 'a', newline='', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(data)


def open_fastq_yield(input_file: Union[str, Path]):
    with open(str(input_file), 'r') as inf:
        for record in SeqIO.parse(inf, 'fastq'):
            yield record

def open_fasta_yield(input_file: Union[str, Path]):
    with open(str(input_file), 'r') as inf:
        for record in SeqIO.parse(inf, 'fasta'):
            yield record


def generate_all_combination_oligos(input_file: Union[str, Path], output_file: Union[str, Path], alphabet: Dict):
    # Define the replacement letters
    replacement_letters = 'ACGT'

    # Create a list of all possible combinations
    combinations = list(product(replacement_letters, repeat=10))

    # Read the input CSV and write the output CSV
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        for row in reader:
            new_row = []
            for cell in row:
                if 'N' in cell:
                    for combo in combinations:
                        new_cell = cell.replace('N', ''.join(combo))
                        new_row.append(new_cell)
                else:
                    new_row.append(cell)
            writer.writerow(new_row)

    print(f"New CSV file '{output_file}' created with all combinations.")


# def find_best_alignment(seq, target_seq):
#     alignment = pairwise2.align.localms(seq, target_seq)
#     print(pairwise2.format_alignment(*alignment[0]))
#
#     x=4





# # def find_sequence_in_target(seq, target_seq, match_score=1, mismatch_score=-1, gap_open=-0.5, gap_extend=-0.1):
# #     aligner = Align.PairwiseAligner()
# #     aligner.mode = 'local'  # Use local alignment to find the best matching subsequence
# #
# #     # Set custom scoring parameters
# #     aligner.match_score = match_score
# #     aligner.mismatch_score = mismatch_score
# #     aligner.open_gap_score = gap_open
# #     aligner.extend_gap_score = gap_extend
# #
# #     # Perform the alignment
# #     alignments = aligner.align(target_seq, seq)
# #
# #     # Get the best alignment (the first one in the sorted list by default)
# #     best_alignment = alignments[0]
# #
# #     # Extract the score of the best alignment
# #     best_score = best_alignment.score
# #
# #     # Extract the start and end positions of seq in target_seq
# #     target_start = best_alignment.aligned[0][0][0]
# #     target_end = best_alignment.aligned[0][-1][1]
# #     query_start = best_alignment.aligned[1][0][0]
# #     query_end = best_alignment.aligned[1][-1][1]
# #
# #     return best_score, target_start, target_end, query_start, query_end, best_alignment
#
def find_best_alignment(seq, target_seq, output_file: Union[str, Path]):
    # seq, target_seq = "ACT", "ACG"
    # Create a pairwise aligner object
    aligner = Align.PairwiseAligner()

    # Set parameters (adjust as needed for your specific requirements)
    aligner.mode = 'local'  # local alignment for finding the best matching subsequence
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1

    # Perform the alignment
    alignments = aligner.align(seq, target_seq)

    # Get the best alignment (highest score)
    best_alignment = alignments[0]

    # Extract alignment details
    aligned_seq1 = best_alignment.aligned[0]
    aligned_seq2 = best_alignment.aligned[1]
    score = best_alignment.score
    start = aligned_seq2[0][0]
    end = aligned_seq2[-1][1]

    result= {
        'aligned_seq': str(seq[aligned_seq1[0][0]:aligned_seq1[-1][1]]),
        'aligned_target_seq': str(target_seq[aligned_seq2[0][0]:aligned_seq2[-1][1]]),
        'score': int(score),
        'start': start,
        'end': end
    }
    # print(best_alignment)
    # print(seq)
    # print(target_seq)
    with open(output_file, "a") as file:
        # Write the string into the file
        file.write(str(best_alignment))

    return result['score'], result['start'], result['end']
#
#
#
# def get_best_alignment(seq: str, target_seq: str):
#     aligner = Align.PairwiseAligner()
#     aligner.mode = 'local'  # Use local alignment to find the best matching subsequence
#
#     # Perform the alignment
#     alignments = aligner.align(target_seq, seq)
#
#     # Get the best alignment (the first one in the sorted list by default)
#     best_alignment = alignments[0]
#
#     # Extract the score of the best alignment
#     best_score = best_alignment.score
#
#     # Extract the start position of seq in target_seq
#     start_position = best_alignment.aligned[0][0][0]
#     end_position = best_alignment.aligned[0][0][1]
#
#     return best_score, start_position, end_position, best_alignment

# Ensure the directory exists before writing any file
def ensure_directory_exists(file_path):
    """Creates the directory for the given file if it doesn't exist."""
    dir_path = os.path.dirname(file_path)
    os.makedirs(dir_path, exist_ok=True)