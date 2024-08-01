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

        self.plot_sequence_length_distribution()
        self.plot_all_nucleotide_distribution()

    def create_output_dirs(self):
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)

    def get_sequences_from_file(self):
        with open(self.sequences_fastq_file, "r") as handle, open(self.sequences_file, "w") as output_handle:
            for record in SeqIO.parse(handle, "fastq"):
                output_handle.write(str(record.seq) + "\n")

    def plot_sequence_length_distribution(self):
        sequence_lengths = []

        with open(self.sequences_file, "r") as input_handle:
            for line in input_handle:
                sequence_lengths.append(len(line.strip()))

        plt.hist(sequence_lengths, bins=40)
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')
        plt.title(f'Sequence Length Distribution, \n Total number of sequences: {len(sequence_lengths)}')
        plt.savefig(self.plot_path + 'sequence_length_distribution.png')
        plt.show()

    def plot_all_nucleotide_distribution(self):
        self.plot_nucleotide_distribution(plot_parts_of_seq=0)
        self.plot_nucleotide_distribution(plot_parts_of_seq=1)
        self.plot_nucleotide_distribution(plot_parts_of_seq=2)
        self.plot_nucleotide_distribution(plot_parts_of_seq=3)
        self.plot_nucleotide_distribution(plot_parts_of_seq=0, stacked=False)
        self.plot_nucleotide_distribution(plot_parts_of_seq=1, stacked=False)
        self.plot_nucleotide_distribution(plot_parts_of_seq=2, stacked=False)
        self.plot_nucleotide_distribution(plot_parts_of_seq=3, stacked=False)

    def plot_nucleotide_distribution(self, plot_parts_of_seq: int, stacked=True):
        with open(self.output_path + '/barcode_with_sequences_distance_dict.csv', 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)

            # Initialize a dictionary to store sequences for each barcode
            barcode_sequences = {}

            for row in csv_reader:
                barcode_idx = row['barcode idx']
                barcode = row['barcode']
                sequence = row['sequence']
                distance = int(row['distance'])
                start_pos = int(row['start_pos'])
                end_pos = int(row['end_pos'])

                if barcode_idx not in barcode_sequences:
                    barcode_sequences[barcode_idx] = {'sequences': []}

                barcode_sequences[barcode_idx]['sequences'].append(
                    {sequence: {'distance': distance, 'start_pos': start_pos, 'end_pos': end_pos}}
                )

                # Process each barcode once it reaches a certain threshold
                if len(barcode_sequences[barcode_idx]['sequences']) >= 1000:  # Adjust threshold as needed
                    self._plot_nucleotide_distribution_for_barcode(barcode, barcode_sequences[barcode],
                                                                   plot_parts_of_seq, stacked)
                    barcode_sequences[barcode_idx] = {'sequences': []}

            # Process remaining sequences for each barcode
            for barcode_idx, data in barcode_sequences.items():
                if data['sequences']:
                    self._plot_nucleotide_distribution_for_barcode(barcode_idx, barcode, data, plot_parts_of_seq,
                                                                   stacked)

    def _plot_nucleotide_distribution_for_barcode(self, barcode_idx, barcode, data, plot_parts_of_seq: int,
                                                  stacked=True):
        os.makedirs(self.plot_path + f'{barcode_idx}/', exist_ok=True)
        title = f'nuc dist for bc {barcode_idx}, stacked={stacked}, '
        sequences = []

        if plot_parts_of_seq == 0:  # entire seq
            title += 'entire seq'
            sequences = [list(seq_info) for seq_dict in data['sequences'] for seq_info in seq_dict]
        elif plot_parts_of_seq == 1:  # without bc
            title += 'comb letters and identifier'
            for seq_dict in data['sequences']:
                for seq_info, info in seq_dict.items():
                    start_pos = info['start_pos']
                    sequences.append(list(seq_info[:start_pos]))
        elif plot_parts_of_seq == 2:  # only bc
            title += 'only bc'
            for seq_dict in data['sequences']:
                for seq_info, info in seq_dict.items():
                    start_pos = info['start_pos']
                    end_pos = info['end_pos']
                    sequences.append(list(seq_info[start_pos:end_pos]))
        elif plot_parts_of_seq == 3:  # only 16 nuc before th BC
            title += 'comb letters and identifier only 16nuc'
            for seq_dict in data['sequences']:
                for seq_info, info in seq_dict.items():
                    start_pos = info['start_pos']
                    sequences.append(list(
                        seq_info[(start_pos - self.identifier_length - self.combinatorial_letters_length):start_pos]))

        max_length = max(len(seq) for seq in sequences)

        # Initialize counts
        counts = {nuc: np.zeros(max_length) for nuc in 'ACGT'}

        # Count nucleotides at each position
        for seq in sequences:
            for i, nuc in enumerate(seq):
                if nuc in counts:
                    counts[nuc][i] += 1

        # Calculate the average counts at each position
        total_sequences = len(sequences)
        avg_counts = {nuc: counts[nuc] / total_sequences for nuc in 'ACGT'}

        # Create a DataFrame for easier plotting
        df = pd.DataFrame(avg_counts)
        df.index.name = 'Position'

        # Plot the histogram
        fig, ax = plt.subplots(figsize=(40, 6))  # Increase the figure size to make the x-axis wider
        df.plot(kind='bar', stacked=stacked, ax=ax, width=1.0)
        plt.title(f'Nucleotide Distribution for Barcode: {barcode_idx}')
        plt.xlabel('Position')
        plt.ylabel('Frequency')
        plt.xticks(rotation=45)
        plt.legend(title='Nucleotide')
        plt.savefig(self.plot_path + f'{barcode_idx}/' + title + '.png')
        plt.show()

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

    def barcode_sequences_dict_to_csv(self) -> None:
        output_csv_file = self.output_path + f'/barcode_with_sequences_distance_dict.csv'
        with open(output_csv_file, "w") as csv_file:
            csv_file.write("barcode,sequence,distance,start_pos,end_pos,number_of_sequences\n")

            for barcode, data in self.barcode_with_sequences_dict.items():
                num_sequences = len(data['sequences'])
                for seq_dict in data['sequences']:
                    for seq, info in seq_dict.items():
                        row = f"{barcode},{seq},{info['distance']},{info['start_pos']},{info['end_pos']},{num_sequences}\n"
                        csv_file.write(row)
