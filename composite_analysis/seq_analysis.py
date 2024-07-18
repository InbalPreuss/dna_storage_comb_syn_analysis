import math
import os
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
import Levenshtein
import numpy as np

import config


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

    return max_distance, start_pos_temp, start_pos_temp + seq_len


class SeqAnalysis:
    def __init__(self, config: dict):
        self.data_path = config['data_path']
        self.output_path = config['output_path']
        self.plot_path = config['plot_path']
        self.sequences_file = config['sequences_file']
        self.sequences_unassembled_file = config['sequences.unassembled_file']
        self.sequences_assembled_file = config['sequences.assembled_file']
        self.adapter_start_location = config['adapter_start_location']
        self.combinatorial_location = config['combinatorial_location']
        self.design_file = config['design_file']
        self.identifier_location = config['identifier_location']
        self.barcode_location = config['barcode_location']
        self.adapter_end_location = config['adapter_end_location']
        self.total_sequence_length = config['total_sequence_length']
        self.total_sequence_length_with_adapters = config['total_sequence_length_with_adapters']
        self.ratio = config['ratio']
        self.max_bc_distance = config['max_bc_distance']
        self.combinatorial_letters_length = config['combinatorial_letters_length']
        self.identifier_length = config['identifier_length']

        self.sequences = []
        self.barcodes = []
        self.barcode_with_sequences_dict = {}

        self.create_output_dirs()

    def run(self):
        self.get_sequences_from_file(self.sequences_file)
        self.get_barcodes()

        # self.plot_sequence_length_distribution()
        self.find_seqs_per_barcodes()
        self.plot_all_nucleotide_distribution()

    def create_output_dirs(self):
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)

    def get_sequences_from_file(self, file_name: str) -> list:
        with open(file_name, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                self.sequences.append(str(record.seq))
        return self.sequences

    def plot_sequence_length_distribution(self):
        sequence_lengths = [len(seq) for seq in self.sequences]
        plt.hist(sequence_lengths, bins=40)
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')
        plt.title(f'Sequence Length Distribution, \n Total number of sequences: {len(self.sequences)}')
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
        for bc_i, (barcode, data) in enumerate(self.barcode_with_sequences_dict.items()):
            os.makedirs(self.plot_path + f'{bc_i}/', exist_ok=True)
            title = f'nuc dist for bc {bc_i}, stacked={stacked}, '
            sequences = []
            if plot_parts_of_seq == 0:  # entire seq
                title += 'entire seq'
                sequences = [list(seq_info) for seq_dict in data['sequences'] for seq_info in seq_dict]
            elif plot_parts_of_seq == 1:  # without bc
                title += 'comb letters and identifier'
                # Extract sequences up to start_pos - 1
                for seq_dict in data['sequences']:
                    for seq_info, info in seq_dict.items():
                        start_pos = info['start_pos']
                        sequences.append(list(seq_info[:start_pos]))
            elif plot_parts_of_seq == 2:  # only bc
                title += f'only bc'
                # Extract sequences up to start_pos - 1
                for seq_dict in data['sequences']:
                    for seq_info, info in seq_dict.items():
                        start_pos = info['start_pos']
                        end_pos = info['end_pos']
                        sequences.append(list(seq_info[start_pos:end_pos]))
            elif plot_parts_of_seq == 3:  # only 16 nuc before th BC
                title += 'comb letters and identifier only 16nuc'
                # Extract sequences only in the length of comb letters + identifier, up to start_pos - 1
                for seq_dict in data['sequences']:
                    for seq_info, info in seq_dict.items():
                        start_pos = info['start_pos']
                        sequences.append(list(seq_info[(start_pos-self.identifier_length-self.combinatorial_letters_length):start_pos]))

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
            plt.title(f'Nucleotide Distribution for Barcode: {barcode}')
            plt.xlabel('Position')
            plt.ylabel('Average Count')
            plt.xticks(rotation=45)
            # plt.tight_layout()
            plt.legend(title='Nucleotide')
            plt.savefig(self.plot_path + f'{bc_i}/' + title + '.png')
            plt.show()

    def get_barcodes(self) -> pd.Series:
        df = pd.read_csv(self.design_file)
        self.barcodes = df['barcode'].str.upper()
        return self.barcodes

    def find_seqs_per_barcodes(self) -> dict:
        for barcode in self.barcodes:
            self.barcode_with_sequences_dict[barcode] = {'sequences': []}
            for seq in self.sequences:
                distance, start_pos, end_pos = is_within_levenshtein_distance(seq=barcode.upper(), target_seq=seq)
                if (distance <= self.max_bc_distance) and (start_pos != math.inf):
                    self.barcode_with_sequences_dict[barcode]['sequences'].append(
                        {seq: {'distance': distance, 'start_pos': start_pos, 'end_pos': end_pos}})
        self.barcode_sequences_dict_to_csv()
        return self.barcode_with_sequences_dict

    def barcode_sequences_dict_to_csv(self):
        rows = []

        for barcode, data in self.barcode_with_sequences_dict.items():
            num_sequences = len(data['sequences'])
            for seq_dict in data['sequences']:
                for seq, info in seq_dict.items():
                    row = {
                        'barcode': barcode,
                        'sequence': seq,
                        'distance': info['distance'],
                        'start_pos': info['start_pos'],
                        'end_pos': info['end_pos'],
                        'number_of_sequences': num_sequences
                    }
                    rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(self.plot_path + f'/barcode_with_sequences_dict.csv', index=False)


if __name__ == '__main__':
    # call SeqAnalysis class
    seq_analysis = SeqAnalysis(config.config)
    # get sequences from file
    seq_analysis.run()
