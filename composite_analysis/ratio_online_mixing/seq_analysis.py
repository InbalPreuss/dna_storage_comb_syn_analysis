import csv
import gzip
import itertools
import math
import os
import shutil
from pathlib import Path
from typing import Union, List, Dict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt

import utilities as uts

BARCODE_CHUNK_SIZE = 1000
SEQUENCE_CHUNK_SIZE = 1000


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
                 blast_db_bc: Union[str, Path], blast_db_bc_identifier: Union[str, Path], blast_database_path: Union[str, Path], blast_db_bc_fasta: Union[str, Path],
                 general_information_file: Union[str, Path], th_minimum_len_reads_to_analyse: int,
                 csv_path: Union[str, Path],
                 reads_chunk_to_fasta_file: Union[str, Path], fasta_path: Union[str, Path], amount_of_bc: int,
                 blast_bc_results: Union[str, Path], blast_db_bc_identifier_fasta: Union[str, Path],
                 barcode_length: int, blast_bc_identifier_results: Union[str, Path], query_results_path: Union[str, Path],
                 barcode_with_sequences_distance_dict_file: Union[str, Path],
                 sequences_R1_fastq_file: Union[str, Path], sequences_R2_fastq_file: Union[str, Path]):
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
        self.blast_db_bc_fasta = blast_db_bc_fasta
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
        self.sequences_R1_fastq_file = sequences_R1_fastq_file
        self.sequences_R2_fastq_file = sequences_R2_fastq_file

        self.create_output_dirs()

    def run(self):
        # self.get_sequences_from_file(input_file=self.sequences_fastq_file, output_file=self.sequences_file)
        # self.get_sequences_from_file(input_file=self.sequences_R1_fastq_file, output_file=self.sequences_file, is_reversed=True)
        # self.get_sequences_from_file(input_file=self.sequences_R1_fastq_file, output_file=self.sequences_file)
        # self.get_sequences_from_file(input_file=self.sequences_R2_fastq_file, output_file=self.sequences_file)

        # TODO: Uncomment this - important
        # # Get sequences from file
        # self.get_sequences_from_file(input_file=self.sequences_R1_fastq_file, output_file=self.sequences_file,
        #                              is_reversed=True)

        # TODO: Uncomment this - important
        # # Get barcode sequences
        # self.get_part_of_seq(seq_part_list=['barcode'], blast_db_fasta=self.blast_db_bc_fasta)

        # Levenstein distance
        self.find_seqs_per_barcodes_using(blast_db_fasta=self.blast_db_bc_fasta,
                                          distance_method='levenshtein',
                                          output_file=self.barcode_with_sequences_distance_dict_file)

    def create_output_dirs(self):
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.blast_database_path, exist_ok=True)
        os.makedirs(self.csv_path, exist_ok=True)
        os.makedirs(self.fasta_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)

    # def get_sequences_from_file(self):
    #     with open(self.sequences_fastq_file, "r") as handle, open(self.sequences_file, "w") as output_handle:
    #         for record in SeqIO.parse(handle, "fastq"):
    #             output_handle.write(str(record.seq) + "\n")

    def get_sequences_from_file(self, input_file: Union[str, Path],
                                output_file: Union[str, Path],
                                is_reversed=False) -> None:
        is_gzip = input_file.endswith(".gz")

        # Handle gzip file extraction if needed
        if is_gzip:
            extracted_file = input_file[:-3]  # Remove .gz extension
            with gzip.open(input_file, "rt") as gz_handle, open(extracted_file, "w") as temp_handle:
                shutil.copyfileobj(gz_handle, temp_handle)
        else:
            extracted_file = input_file

        # Process sequences
        with open(extracted_file, "r") as handle, open(output_file, "w") as output_handle:
            for record in SeqIO.parse(handle, "fastq"):
                if is_reversed:
                    sequence = str(record.seq[::-1])
                    output_handle.write(sequence + "\n")
                else:
                    sequence = str(record.seq)
                    output_handle.write(sequence + "\n")

                # Write also the reverse complement
                sequence_rev_comp = str(Seq(sequence).reverse_complement())
                output_handle.write(sequence_rev_comp + "\n")




        # Recompress the original input file if it was originally a gzip file
        if is_gzip:
            with open(extracted_file, "rb") as f_in, gzip.open(input_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

            # Cleanup temporary extracted file
            os.remove(extracted_file)

    def get_part_of_seq(self, seq_part_list: List[str], blast_db_fasta: Union[str, Path]):
        with open(blast_db_fasta, 'w') as fasta_file:
            chunk_iter = pd.read_csv(self.design_file, chunksize=BARCODE_CHUNK_SIZE)
            for chunk in chunk_iter:
                # chunk['barcode'] = chunk['barcode']

                for index, row in chunk.iterrows():
                    seq_part_from_row = ''
                    for seq_part_i in seq_part_list:
                        seq_part_from_row += row[seq_part_i]
                    fasta_file.write(f">{int(row['barcode_id'])}\n{seq_part_from_row}\n")

    # def find_seqs_per_barcodes_using(self, blast_db_fasta: Union[str, Path], distance_method: str, output_file: Union[str, Path]) -> None:
    #     output_csv_file = self.output_path + output_file
    #
    #     with open(output_csv_file, "w", newline='') as csv_file:
    #         csv_writer = csv.writer(csv_file)
    #         csv_writer.writerow(
    #             ["barcode_id", "barcode", "sequence", "distance", "start_pos", "end_pos", "number_of_sequences", "sequence_length"])
    #
    #         reads_iter = uts.open_fasta_yield(blast_db_fasta)
    #         reads_chunk = list(itertools.islice(reads_iter, CHUNK_SIZE))
    #
    #         for barcode_info in reads_chunk:
    #         # for barcode_info in reversed(reads_chunk): #TODO: delete
    #             barcode = barcode_info.seq
    #             sequence_count = 0
    #             matching_sequences = []
    #             scores = []
    #
    #             barcode_i = int(barcode_info.id)
    #
    #             with open(self.sequences_file, "r") as input_handle:
    #                 for seq in input_handle:
    #                     seq = seq.strip()
    #                     sequence_length = len(seq)
    #                     if distance_method == 'levenshtein':
    #                         distance, start_pos, end_pos = uts.is_within_levenshtein_distance(seq=barcode,
    #                                                                                           target_seq=seq,
    #                                                                                           max_distance=self.max_bc_distance)
    #                         if (distance <= self.max_bc_distance) and (start_pos != math.inf):
    #                             matching_sequences.append((seq, distance, start_pos, end_pos, sequence_length))
    #                             sequence_count += 1
    #                     elif distance_method == 'alignment':
    #                         score, start_pos, end_pos = uts.find_best_alignment(seq=barcode, target_seq=seq, output_file=self.output_path + '/alignment.txt')
    #                         if score >= 20:
    #                             scores.append(score)
    #                             matching_sequences.append((seq, score, start_pos, end_pos, sequence_length))
    #                             sequence_count += 1
    #
    #
    #             for seq, distance, start_pos, end_pos, sequence_length in matching_sequences:
    #                     csv_writer.writerow([barcode_i, barcode, seq, distance, start_pos, end_pos, sequence_count, sequence_length])
    #             if distance_method == 'alignment':
    #                 # plot histogram of scores
    #                 plt.hist(scores)
    #                 plt.xlabel('Score')
    #                 plt.ylabel('Frequency')
    #                 plt.title(f'Bc={barcode_i+1}, Score Distribution, \n Total number of scores: {len(scores)}')
    #                 plt.savefig(self.plot_path + f'bc={barcode_i+1}_score_distribution.png')
    #                 plt.show()
    #                 plt.close()


    def find_seqs_per_barcodes_using(self, blast_db_fasta: Union[str, Path], distance_method: str,
                                     output_file: Union[str, Path]) -> None:
        output_csv_file = self.output_path + output_file

        with open(output_csv_file, "w", newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["barcode_id", "barcode", "sequence", "distance", "start_pos", "end_pos", "number_of_sequences",
                 "sequence_length"]
            )

            # Read barcode sequences in chunks
            reads_iter = uts.open_fasta_yield(blast_db_fasta)
            while True:
                reads_chunk = list(itertools.islice(reads_iter, BARCODE_CHUNK_SIZE))
                if not reads_chunk:
                    break  # Stop when no more sequences to process

                for barcode_info in reads_chunk:
                    barcode = barcode_info.seq
                    barcode_i = int(barcode_info.id)

                    matching_sequences = []
                    scores = []

                    # Open sequence file ONCE and read in chunks
                    with open(self.sequences_file, "r") as input_handle:
                        while True:
                            sequence_batch = list(itertools.islice(input_handle, SEQUENCE_CHUNK_SIZE))
                            if not sequence_batch:
                                break  # Stop when no more sequences to process

                            # Process each sequence in the batch
                            for seq in sequence_batch:
                                seq = seq.strip()
                                sequence_length = len(seq)

                                if distance_method == 'levenshtein':
                                    distance, start_pos, end_pos = uts.is_within_levenshtein_distance(
                                        seq=barcode, target_seq=seq, max_distance=self.max_bc_distance
                                    )
                                    if distance <= self.max_bc_distance and start_pos != math.inf:
                                        matching_sequences.append((seq, distance, start_pos, end_pos, sequence_length))

                                elif distance_method == 'alignment':
                                    score, start_pos, end_pos = uts.find_best_alignment(
                                        seq=barcode, target_seq=seq, output_file=self.output_path + '/alignment.txt'
                                    )
                                    if score >= 20:
                                        scores.append(score)
                                        matching_sequences.append((seq, score, start_pos, end_pos, sequence_length))

                    sequence_count = len(matching_sequences)  # More efficient than incrementing in loop

                    # Write results in bulk to reduce I/O operations
                    csv_writer.writerows(
                        [[barcode_i, barcode, seq, distance, start_pos, end_pos, sequence_count, sequence_length]
                         for seq, distance, start_pos, end_pos, sequence_length in matching_sequences]
                    )

                    # Plot histogram if using alignment method
                    if distance_method == 'alignment' and scores:
                        plt.hist(scores, bins=np.linspace(min(scores), max(scores), 20))  # Use bins for efficiency
                        plt.xlabel('Score')
                        plt.ylabel('Frequency')
                        plt.title(f'Bc={barcode_i + 1}, Score Distribution, \n Total number of scores: {len(scores)}')
                        plt.savefig(f"{self.plot_path}bc={barcode_i + 1}_score_distribution.png")
                        plt.close()  # Avoid displaying plots for performance