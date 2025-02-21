import csv
import math
from pathlib import Path
from typing import Union, List, Dict
from itertools import product

import Levenshtein
from Bio import SeqIO, Align, pairwise2
import numpy as np
from scipy.special import rel_entr
from scipy.stats import entropy


def is_within_levenshtein_distance(seq: str, target_seq: str, max_distance=4) -> tuple[int, int, int]:
    distance, start_pos_temp, end_pos = math.inf, math.inf, math.inf
    seq_len = len(seq)
    target_seq_len = len(target_seq)

    for start_pos in range(target_seq_len - seq_len + 1):
        subseq = target_seq[start_pos:start_pos + seq_len]
        distance = Levenshtein.distance(seq, subseq)
        if distance <= max_distance:
            max_distance = distance
            start_pos_temp = start_pos
        if max_distance == 0:  # Early exit on perfect match
            break

    return max_distance, start_pos_temp, start_pos_temp + seq_len


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

def kl_divergence(p, q, epsilon=1e-10) -> float:
    """
    Compute KL divergence safely, avoiding division by zero issues.

    :param p: First probability distribution (observed/nucleotide_percentages).
    :param q: Second probability distribution (expected/design_distribution).
    :param epsilon: Small smoothing value to avoid log(0) issues.
    :return: KL divergence value.
    """
    p = np.array(p) / np.sum(p)
    q = np.array(q) / np.sum(q)

    # Apply smoothing to avoid zero probabilities
    q = q + epsilon
    q = q / np.sum(q)  # Re-normalize after smoothing

    return entropy(p, q)


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
    print(best_alignment)
    print(seq)
    print(target_seq)
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