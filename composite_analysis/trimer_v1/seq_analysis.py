import csv
import gzip
import itertools
import math
import os
from pathlib import Path
from typing import Union, List, Dict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt

import utilities as uts

CHUNK_SIZE = 1000

# Define chunk sizes to prevent memory overload
BARCODE_CHUNK_SIZE = 1000
SEQUENCE_CHUNK_SIZE = 1000


# Safe FASTQ file reader (Handles .gz and format issues)
# FASTQ File Reader with Reverse & Reverse Complement Processing
def read_fastq_in_chunks(file_path, reverse_sequences=False, reverse_complement_sequences=False):
    """
    Reads a FASTQ file in chunks to optimize memory usage.
    Can apply reverse and/or reverse complement transformations.
    """
    is_gzip = file_path.endswith(".gz")

    try:
        with (gzip.open(file_path, "rt") if is_gzip else open(file_path, "r")) as handle:
            chunk = []
            for record in SeqIO.parse(handle, "fastq"):
                sequence = str(record.seq)

                # Apply transformations based on config
                if reverse_sequences:
                    sequence = sequence[::-1]  # Reverse the sequence
                    chunk.append(sequence)
                if reverse_complement_sequences:
                    sequence = str(Seq(sequence).reverse_complement())  # Reverse complement
                    chunk.append(sequence)
                if not reverse_sequences and not reverse_complement_sequences:
                    chunk.append(sequence)

                # Yield chunk when full
                if len(chunk) >= SEQUENCE_CHUNK_SIZE:
                    yield chunk
                    chunk = []

            if chunk:
                yield chunk

    except Exception as e:
        print(f"ERROR: Failed to read FASTQ file '{file_path}': {e}")


class SeqAnalysis:
    def __init__(self, data_path: Union[str, Path], output_path: Union[str, Path], plot_path: Union[str, Path],
                 sequences_fastq_file: Union[str, Path], sequences_file: Union[str, Path],
                 sequences_unassembled_file: Union[str, Path],
                 sequences_assembled_file: Union[str, Path],
                 adapter_start_location: List, combinatorial_location: List, design_file: Union[str, Path],
                 barcode_location: List,
                 adapter_end_location: List, total_sequence_length: int, total_sequence_length_with_adapters: int,
                 alphabet: Dict, max_bc_distance: int,
                 combinatorial_letters_length: int,
                 blast_db_bc: Union[str, Path], blast_db_bc_identifier: Union[str, Path],
                 blast_database_path: Union[str, Path], blast_db_bc_fasta: Union[str, Path], blast_db_seq_design_fasta: Union[str, Path],
                 blast_database_bc_path: Union[str, Path], blast_database_seq_design_path: Union[str, Path],
                 general_information_file: Union[str, Path], th_minimum_len_reads_to_analyse: int,
                 csv_path: Union[str, Path],
                 reads_chunk_to_fasta_file: Union[str, Path], fasta_path: Union[str, Path], amount_of_bc: int,
                 blast_bc_results: Union[str, Path], blast_db_bc_identifier_fasta: Union[str, Path],
                 barcode_length: int, blast_bc_identifier_results: Union[str, Path],
                 query_results_path: Union[str, Path], barcode_with_sequences_distance_dict_file: Union[str, Path],
                 blast_database_adapter_path: Union[str, Path], blast_db_adapter1_fasta: Union[str, Path], blast_db_adapter2_fasta: Union[str, Path],
                 reverse_sequences: bool = False, reverse_complement_sequences: bool = False):
        self.data_path = data_path
        self.output_path = output_path
        self.plot_path = plot_path
        self.csv_path = csv_path
        self.sequences_fastq_file = sequences_fastq_file
        self.sequences_file = sequences_file
        self.sequences_unassembled_file = sequences_unassembled_file
        self.sequences_assembled_file = sequences_assembled_file
        self.adapter_start_location = adapter_start_location
        self.combinatorial_location = combinatorial_location
        self.design_file = design_file
        self.barcode_location = barcode_location
        self.adapter_end_location = adapter_end_location
        self.total_sequence_length = total_sequence_length
        self.total_sequence_length_with_adapters = total_sequence_length_with_adapters
        self.alphabet = alphabet
        self.max_bc_distance = max_bc_distance
        self.combinatorial_letters_length = combinatorial_letters_length
        self.blast_db_bc = blast_db_bc
        self.blast_db_bc_identifier = blast_db_bc_identifier
        self.blast_database_path = blast_database_path
        self.blast_database_seq_design_path = blast_database_seq_design_path
        self.blast_database_bc_path = blast_database_bc_path
        self.blast_db_bc_fasta = blast_db_bc_fasta
        self.blast_db_seq_design_fasta = blast_db_seq_design_fasta
        self.general_information_file = general_information_file
        self.th_minimum_len_reads_to_analyse = th_minimum_len_reads_to_analyse
        self.reads_chunk_to_fasta_file = reads_chunk_to_fasta_file
        self.fasta_path = fasta_path
        self.amount_of_bc = amount_of_bc
        self.blast_bc_results = blast_bc_results
        self.blast_db_bc_identifier_fasta = blast_db_bc_identifier_fasta
        self.barcode_length = barcode_length
        self.blast_bc_identifier_results = blast_bc_identifier_results
        self.query_results_path = query_results_path
        self.barcode_with_sequences_distance_dict_file = barcode_with_sequences_distance_dict_file
        self.reverse_sequences = reverse_sequences
        self.reverse_complement_sequences = reverse_complement_sequences
        self.blast_database_adapter_path = blast_database_adapter_path
        self.blast_db_adapter1_fasta = blast_db_adapter1_fasta
        self.blast_db_adapter2_fasta = blast_db_adapter2_fasta

        self.create_output_dirs()

    def run(self):
        self.extract_sequences()
        self.extract_part_of_seqs()

        # Levenstein distance
        self.find_seqs_per_barcodes_using_chosen_method(blast_db_bc_fasta=self.blast_db_bc_fasta,
                                                        blast_db_seq_design_fasta=self.blast_db_seq_design_fasta,
                                                        distance_method='levenshtein',
                                                        output_file=self.barcode_with_sequences_distance_dict_file)
        # self.find_seqs_per_barcodes_using(blast_db_fasta=self.blast_db_bc_fasta,
        #                                   distance_method='alignment',
        #                                   output_file=self.barcode_with_sequences_distance_dict_file)

    def create_output_dirs(self):
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.blast_database_path, exist_ok=True)
        os.makedirs(self.blast_database_seq_design_path, exist_ok=True)
        os.makedirs(self.blast_database_bc_path, exist_ok=True)
        os.makedirs(self.blast_database_adapter_path, exist_ok=True)
        os.makedirs(self.csv_path, exist_ok=True)
        os.makedirs(self.fasta_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)

    def extract_sequences(self):
        # with open(self.sequences_fastq_file, "r") as handle, open(self.sequences_file, "w") as output_handle:
        #     for record in SeqIO.parse(handle, "fastq"):
        #         output_handle.write(str(record.seq) + "\n")

        """Extracts sequences from FASTQ and saves them to a text file."""

        uts.ensure_directory_exists(self.sequences_file)  # Ensure output folder exists

        with open(self.sequences_file, "w") as output_handle:
            for chunk in read_fastq_in_chunks(self.sequences_fastq_file, reverse_sequences=self.reverse_sequences,
                                                   reverse_complement_sequences=self.reverse_complement_sequences):
                output_handle.write("\n".join(chunk) + "\n")  # Write batch at once

    def extract_part_of_seq(self, seq_part_list: List[str], blast_db_fasta: Union[str, Path]):
        with open(blast_db_fasta, 'w') as fasta_file:
            chunk_iter = pd.read_csv(self.design_file, chunksize=CHUNK_SIZE)
            for chunk in chunk_iter:
                # chunk['barcode'] = chunk['barcode']

                for index, row in chunk.iterrows():
                    seq_part_from_row = ''
                    for seq_part_i in seq_part_list:
                        seq_part_from_row += row[seq_part_i]
                    fasta_file.write(f">{int(row['barcode_id'])}\n{seq_part_from_row}\n")

    def find_seqs_per_barcodes_using_chosen_method(self, blast_db_bc_fasta: Union[str, Path],
                                                   blast_db_seq_design_fasta: Union[str, Path],
                                                   distance_method: str,
                                                   output_file: Union[str, Path]) -> None:
        output_csv_file = self.output_path + output_file

        with open(output_csv_file, "w", newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["barcode_id", "barcode", "sequence", "distance", "start_pos", "end_pos", "number_of_sequences",
                 "sequence_length"])

            bc_iter = uts.open_fasta_yield(blast_db_bc_fasta)
            bc_chunk = list(itertools.islice(bc_iter, CHUNK_SIZE))

            seq_design_iter = uts.open_fasta_yield(blast_db_seq_design_fasta)
            seq_design_chunk = list(itertools.islice(seq_design_iter, CHUNK_SIZE))

            adapter1_iter = uts.open_fasta_yield(self.blast_db_adapter1_fasta)
            adapter1_chunk = list(itertools.islice(adapter1_iter, CHUNK_SIZE))

            adapter2_iter = uts.open_fasta_yield(self.blast_db_adapter2_fasta)
            adapter2_chunk = list(itertools.islice(adapter2_iter, CHUNK_SIZE))


            for barcode_info, seq_design_info, adapter1_info, adapter2_info in zip(bc_chunk,seq_design_chunk, adapter1_chunk, adapter2_chunk):
                sequence_count = 0
                matching_sequences = []
                scores = []

                bc_seq = barcode_info.seq
                barcode_i = int(barcode_info.id)

                seq_design = seq_design_info.seq
                seq_adapter1 = adapter1_info.seq
                seq_adapter2 = adapter2_info.seq
                seq_design_i = int(seq_design_info.id)

                # Header: First column is the full sequence, then each nucleotide
                header = ["Compared Sequence","start_pos","end_pos"] + list(seq_design)
                annotations_data = []

                with open(self.sequences_file, "r") as input_handle:
                    for read in input_handle:
                        read = read.strip()
                        sequence_length = len(read)


                        if distance_method == 'levenshtein':
                            distance, start_pos_bc, end_pos_bc, annotations, trimmed_read = uts.is_within_levenshtein_distance(
                                bc_seq=bc_seq,
                                seq_design=seq_design,
                                read=read,
                                max_distance=self.max_bc_distance,
                                seq_adapter1=seq_adapter1,
                                seq_adapter2=seq_adapter2)
                            if (distance <= self.max_bc_distance) and (start_pos_bc != math.inf):
                                matching_sequences.append((trimmed_read, distance, start_pos_bc, end_pos_bc, sequence_length))
                                sequence_count += 1
                                annotations_data.append([trimmed_read,start_pos_bc,end_pos_bc] + annotations)
                        elif distance_method == 'alignment':
                            score, start_pos_bc, end_pos_bc = uts.find_best_alignment(seq=bc_seq, target_seq=read,
                                                                                output_file=self.output_path + '/alignment.txt')
                            if score >= 20:
                                scores.append(score)
                                matching_sequences.append((read, score, start_pos_bc, end_pos_bc, sequence_length))
                                sequence_count += 1
                        # elif distance_method == 'levenshtein_':
                        #     pass

                for read, distance, start_pos_bc, end_pos_bc, sequence_length in matching_sequences:
                    csv_writer.writerow(
                        [barcode_i, bc_seq, read, distance, start_pos_bc, end_pos_bc, sequence_count, sequence_length])

                # Save to CSV
                df = pd.DataFrame(annotations_data, columns=header)
                df.to_csv(self.csv_path+f"/{barcode_i}/annotations_{distance_method}.csv", index=False)
                if distance_method == 'alignment':
                    # plot histogram of scores
                    plt.hist(scores)
                    plt.xlabel('Score')
                    plt.ylabel('Frequency')
                    plt.title(f'Bc={barcode_i + 1}, Score Distribution, \n Total number of scores: {len(scores)}')
                    plt.savefig(self.plot_path + f'bc={barcode_i + 1}_score_distribution.png')
                    plt.show()
                    plt.close()

    def extract_part_of_seqs(self):
        self.extract_part_of_seq(seq_part_list=['barcode'], blast_db_fasta=self.blast_db_bc_fasta)
        self.extract_part_of_seq(seq_part_list=['sequence'], blast_db_fasta=self.blast_db_seq_design_fasta)
        self.extract_part_of_seq(seq_part_list=['Adapter1'], blast_db_fasta=self.blast_db_adapter1_fasta)
        self.extract_part_of_seq(seq_part_list=['Adapter2'], blast_db_fasta=self.blast_db_adapter2_fasta)
