import csv
import math
import os
from pathlib import Path
from typing import Dict, Union, List

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
import numpy as np

from utilities import is_within_levenshtein_distance


class SeqAnalysis:
    def __init__(self, data_path: Union[str, Path], output_path: Union[str, Path], plot_path: Union[str, Path],
                 sequences_fastq_file: Union[str, Path], sequences_file: Union[str, Path],
                 sequences_unassembled_file: Union[str, Path],
                 sequences_assembled_file: Union[str, Path],
                 adapter_start_location: List, combinatorial_location: List, design_file: Union[str, Path],
                 identifier_location: List, barcode_location: List,
                 adapter_end_location: List, total_sequence_length: int, total_sequence_length_with_adapters: int,
                 alphabet: Dict, max_bc_distance: int,
                 combinatorial_letters_length: int, identifier_length: int):
        self.data_path = data_path
        self.output_path = output_path
        self.plot_path = plot_path
        self.sequences_fastq_file = sequences_fastq_file
        self.sequences_file = sequences_file
        self.sequences_unassembled_file = sequences_unassembled_file
        self.sequences_assembled_file = sequences_assembled_file
        self.adapter_start_location = adapter_start_location
        self.combinatorial_location = combinatorial_location
        self.design_file = design_file
        self.identifier_location = identifier_location
        self.barcode_location = barcode_location
        self.adapter_end_location = adapter_end_location
        self.total_sequence_length = total_sequence_length
        self.total_sequence_length_with_adapters = total_sequence_length_with_adapters
        self.alphabet = alphabet
        self.max_bc_distance = max_bc_distance
        self.combinatorial_letters_length = combinatorial_letters_length
        self.identifier_length = identifier_length

        self.sequences = []
        self.barcodes = []
        self.barcode_with_sequences_dict = {}

        self.create_output_dirs()

    def run(self):
        self.get_sequences_from_file()
        self.get_barcodes()
        self.find_seqs_per_barcodes_using_distance()

    def create_output_dirs(self):
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)

    def get_sequences_from_file(self):
        with open(self.sequences_fastq_file, "r") as handle, open(self.sequences_file, "w") as output_handle:
            for record in SeqIO.parse(handle, "fastq"):
                output_handle.write(str(record.seq) + "\n")

    def get_barcodes(self) -> pd.Series:
        df = pd.read_csv(self.design_file)
        self.barcodes = df['barcode'].str.upper()
        return self.barcodes

    def find_seqs_per_barcodes_using_distance(self) -> None:
        output_csv_file = self.output_path + f'/barcode_with_sequences_distance_dict.csv'
        with open(output_csv_file, "w", newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["barcode idx", "barcode", "sequence", "distance", "start_pos", "end_pos", "number_of_sequences"])

            for barcode_i, barcode in enumerate(self.barcodes):
                sequence_count = 0
                matching_sequences = []

                with open(self.sequences_file, "r") as input_handle:
                    for seq in input_handle:
                        seq = seq.strip()
                        distance, start_pos, end_pos = is_within_levenshtein_distance(seq=barcode.upper(),
                                                                                      target_seq=seq,
                                                                                      max_distance=self.max_bc_distance)
                        if (distance <= self.max_bc_distance) and (start_pos != math.inf):
                            matching_sequences.append((seq, distance, start_pos, end_pos))
                            sequence_count += 1

                for seq, distance, start_pos, end_pos in matching_sequences:
                    csv_writer.writerow([barcode_i, barcode, seq, distance, start_pos, end_pos, sequence_count])
