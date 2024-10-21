import csv
import os
from collections import defaultdict
from pathlib import Path
from typing import Union, List, Dict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


class Plot():
    def __init__(self, sequences_file: Union[str, Path], plot_path: Union[str, Path], output_path: Union[str, Path],
                 data_path: Union[str, Path], sequences_fastq_file: Union[str, Path],
                 sequences_unassembled_file: Union[str, Path], sequences_assembled_file: Union[str, Path],
                 adapter_start_location: List, combinatorial_location: List, design_file: Union[str, Path],
                 barcode_location: List, adapter_end_location: List,
                 total_sequence_length: int, total_sequence_length_with_adapters: int, alphabet: Dict,
                 max_bc_distance: int, combinatorial_letters_length: int):
        self.sequences_file = sequences_file
        self.plot_path = plot_path
        self.output_path = output_path
        self.data_path = data_path
        self.sequences_fastq_file = sequences_fastq_file
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

    def run(self):
        # self.create_output_dirs()

        # self.plot_sequence_length_distribution()
        # self.plot_all_nucleotide_distribution(input_file='/barcode_with_sequences_distance_dict.csv')
        self.calculate_nucleotide_percentage_single_row(input_file='/barcode_with_sequences_distance_dict.csv')

    def create_output_dirs(self):
        os.makedirs(self.plot_path, exist_ok=True)

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

    def plot_all_nucleotide_distribution(self, input_file: Union[str, Path]):
        self.plot_nucleotide_distribution(plot_parts_of_seq=0, input_file=input_file)
        # self.plot_nucleotide_distribution(plot_parts_of_seq=1, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=2, input_file=input_file)
        # self.plot_nucleotide_distribution(plot_parts_of_seq=3, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=0, stacked=False, input_file=input_file)
        # self.plot_nucleotide_distribution(plot_parts_of_seq=1, stacked=False, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=2, stacked=False, input_file=input_file)
        # self.plot_nucleotide_distribution(plot_parts_of_seq=3, stacked=False, input_file=input_file)

    def plot_nucleotide_distribution(self, plot_parts_of_seq: int, input_file: Union[str, Path], stacked=True):
        with open(self.output_path + input_file, 'r') as csv_file:
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

                if distance > 0:
                    continue

                if barcode_idx not in barcode_sequences:
                    barcode_sequences[barcode_idx] = {'sequences': []}

                barcode_sequences[barcode_idx]['sequences'].append(
                    {sequence: {'distance': distance, 'start_pos': start_pos, 'end_pos': end_pos}}
                )

                # Process each barcode once it reaches a certain threshold
                if len(barcode_sequences[barcode_idx]['sequences']) >= 1000:  # Adjust threshold as needed
                    title = self._plot_nucleotide_distribution_for_barcode(barcode_idx, barcode, barcode_sequences[barcode_idx],
                                                                           plot_parts_of_seq, stacked)
                    barcode_sequences[barcode_idx] = {'sequences': []}

            # Process remaining sequences for each barcode
            for barcode_idx, data in barcode_sequences.items():
                if data['sequences']:
                    title = self._plot_nucleotide_distribution_for_barcode(barcode_idx, barcode, data,
                                                                           plot_parts_of_seq, stacked)
            plt.savefig(self.plot_path + f'{barcode_idx}/' + title + '.png')
            plt.close()

    def _plot_nucleotide_distribution_for_barcode(self, barcode_idx, barcode, data, plot_parts_of_seq: int,
                                                  stacked=True, ):
        os.makedirs(self.plot_path + f'{barcode_idx}/', exist_ok=True)
        title = f'nuc dist for bc {barcode_idx}, stacked={stacked}, '
        sequences = []
        start_positions = []
        end_positions = []

        if plot_parts_of_seq == 0:  # entire seq
            title += 'entire seq'
            sequences = [list(seq_info) for seq_dict in data['sequences'] for seq_info in seq_dict]
        elif plot_parts_of_seq == 1:  # without bc
            title += 'comb letters and bc'
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
                    start_positions.append(start_pos)
                    end_positions.append(end_pos)
        elif plot_parts_of_seq == 3:  # only 16 nuc before th BC
            title += 'comb letters and identifier only 16nuc'
            for seq_dict in data['sequences']:
                for seq_info, info in seq_dict.items():
                    start_pos = info['start_pos']
                    seq_start_pos = (start_pos - self.identifier_length - self.combinatorial_letters_length)
                    sequences.append(list(
                        seq_info[seq_start_pos:start_pos]))
                    start_positions.append(-16)
                    end_positions.append(-1)

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
        fig, ax = plt.subplots(figsize=(40, 6))
        df.plot(kind='bar', stacked=stacked, ax=ax, width=1.0)
        plt.title(f'Nucleotide Distribution for Barcode: {barcode_idx}')
        plt.xlabel('Position')
        plt.ylabel('Frequency')
        if len(start_positions) != 0 and len(end_positions) != 0:
            # Generate x-tick labels based on start and end positions
            xticks = np.arange(max_length)
            xtick_labels = [f"{start_positions[0] + i}" for i in xticks]  # Modify this logic as per the requirement

            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels, rotation=45)
        else:
            plt.xticks(rotation=45)
        plt.legend(title='Nucleotide')
        plt.savefig(self.plot_path + f'{barcode_idx}/' + title + '.png')
        # plt.show()

        return title

    from collections import defaultdict

    # # Function to calculate nucleotide percentage for each barcode ID as a single row
    # def calculate_nucleotide_percentage_single_row(self, input_file: Union[str, Path]):
    #
    #     df = pd.read_csv(self.output_path + input_file)
    #
    #     # Filter the dataframe where the 'distance' column is 0
    #     df_filtered = df[df['distance'] == 0]
    #
    #     # Dictionary to store results for each barcode
    #     result_rows = []
    #
    #     # Group the dataframe by 'barcode idx'
    #     grouped = df_filtered.groupby('barcode idx')
    #
    #     # For each group (barcode id)
    #     for barcode_id, group in grouped:
    #         max_len = group['sequence'].apply(len).max()
    #         nucleotide_count = defaultdict(lambda: [0] * max_len)
    #         total_sequences = [0] * max_len
    #
    #         # Process each sequence in the group
    #         for _, row in group.iterrows():
    #             sequence = row['sequence']
    #             start_pos = row['start_pos']
    #             for i, nucleotide in enumerate(sequence[start_pos:], start=start_pos):
    #                 if i < max_len:
    #                     nucleotide_count[nucleotide][i] += 1
    #                     total_sequences[i] += 1
    #
    #         # Calculate the percentage for each nucleotide at each position
    #         percentage_row = []
    #         for i in range(max_len):
    #             if total_sequences[i] > 0:
    #                 for nucleotide in 'ACGT':
    #                     percentage_row.append((nucleotide_count[nucleotide][i] / total_sequences[i]) * 100)
    #             else:
    #                 percentage_row.extend([0, 0, 0, 0])  # Append 0 if no sequence at this position
    #
    #         result_rows.append([barcode_id] + percentage_row)
    #
    #     # Dynamically adjust the column names based on the length of the nucleotide percentage rows
    #     num_columns = len(result_rows[0]) - 1  # Exclude the barcode_id
    #     columns_filtered = ['barcode_id'] + [f'{nucleotide}_{i}' for i in range(num_columns // 4) for nucleotide in
    #                                          'ACGT']
    #
    #     # Create a DataFrame for the filtered results
    #     percentage_df_filtered = pd.DataFrame(result_rows, columns=columns_filtered)
    #
    #     # Save the result to a CSV file
    #     output_csv_path_filtered = self.output_path + '/csv/nucleotide_percentage_single_row_filtered_fixed.csv'
    #     percentage_df_filtered.to_csv(output_csv_path_filtered, index=False)

    def calculate_nucleotide_percentage_single_row(self, input_file: Union[str, Path]):

        df = pd.read_csv(self.output_path + input_file)

        # Filter the dataframe where the 'distance' column is 0
        df_filtered = df[df['distance'] == 0]

        # Dictionary to store results for each barcode
        result_rows = []

        # Group the dataframe by 'barcode idx'
        grouped = df_filtered.groupby('barcode idx')

        # For each group (barcode id)
        for barcode_id, group in grouped:
            max_len = group['sequence'].apply(len).max()
            nucleotide_count = defaultdict(lambda: [0] * max_len)
            total_sequences = [0] * max_len

            # Process each sequence in the group
            for _, row in group.iterrows():
                sequence = row['sequence']
                start_pos = row['start_pos']
                # for i, nucleotide in enumerate(sequence[start_pos:], start=start_pos):
                for i, nucleotide in enumerate(sequence[start_pos:]):
                    if i < max_len - start_pos:
                        nucleotide_count[nucleotide][i] += 1
                        total_sequences[i] += 1

            # Calculate the percentage for each nucleotide at each position
            percentage_row = []
            for i in range(max_len):
                if total_sequences[i] > 0:
                    for nucleotide in 'ACGT':
                        percentage_row.append((nucleotide_count[nucleotide][i] / total_sequences[i]) * 100)
                else:
                    percentage_row.extend([0, 0, 0, 0])  # Append 0 if no sequence at this position

            result_rows.append([barcode_id] + percentage_row)

        # Dynamically adjust the column names based on the length of the nucleotide percentage rows
        num_columns = len(result_rows[0]) - 1  # Exclude the barcode_id
        columns_filtered = ['barcode_id'] + [f'{nucleotide}_{i}' for i in range(num_columns // 4) for nucleotide in
                                             'ACGT']

        # Create a DataFrame for the filtered results
        percentage_df_filtered = pd.DataFrame(result_rows, columns=columns_filtered)

        # Save the result to a CSV file
        output_csv_path_filtered = self.output_path + '/csv/nucleotide_percentage_single_row_filtered_fixed.csv'
        percentage_df_filtered.to_csv(output_csv_path_filtered, index=False)

        # Now, plot the stacked bar graph for each barcode id
        self.plot_stacked_bars(percentage_df_filtered)

    def plot_stacked_bars(self, df):
        # For each unique barcode id in the dataframe, plot the stacked bars
        for barcode_id in df['barcode_id'].unique():
            self.plot_stacked_bar_for_barcode(df, barcode_id, 0, 43)

    def plot_stacked_bar_for_barcode(self, df, barcode_id, start_pos, end_pos):
        # Filter data for the given barcode id
        df_barcode = df[df['barcode_id'] == barcode_id]

        # Extract percentage columns for A, C, G, T
        total_positions = len(df_barcode.columns[1:]) // 4  # Total number of positions
        positions = range(start_pos, min(end_pos + 1, total_positions))  # Restrict positions to the start and end range

        A_percentages = df_barcode[[col for col in df_barcode.columns if 'A_' in col]].values.flatten()[
                        start_pos:end_pos + 1]
        C_percentages = df_barcode[[col for col in df_barcode.columns if 'C_' in col]].values.flatten()[
                        start_pos:end_pos + 1]
        G_percentages = df_barcode[[col for col in df_barcode.columns if 'G_' in col]].values.flatten()[
                        start_pos:end_pos + 1]
        T_percentages = df_barcode[[col for col in df_barcode.columns if 'T_' in col]].values.flatten()[
                        start_pos:end_pos + 1]

        # Prepare data for the stacked bar plot
        fig, ax = plt.subplots(figsize=(45, 15))

        # Plot the bars
        bar_width = 0.8
        A_bars = ax.bar(positions, A_percentages, bar_width, label='A')
        C_bars = ax.bar(positions, C_percentages, bar_width, bottom=A_percentages, label='C')
        G_bars = ax.bar(positions, G_percentages, bar_width, bottom=A_percentages + C_percentages, label='G')
        T_bars = ax.bar(positions, T_percentages, bar_width, bottom=A_percentages + C_percentages + G_percentages,
                        label='T')

        # Annotate bars with the actual percentages
        for bars, nucleotide in zip([A_bars, C_bars, G_bars, T_bars], ['A', 'C', 'G', 'T']):
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    if nucleotide in ['A', 'G']:
                        text_color = 'white'  # White text for A and G
                    else:
                        text_color = 'black'  # Black text for C and T
                    ax.annotate(f'{height:.1f}%',
                                xy=(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2),
                                xytext=(0, 25),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='center', fontsize=20, rotation=90, color=text_color)

        # Set labels and title
        ax.set_xlabel('Position in Sequence', fontsize=20)
        ax.set_ylabel('Percentage (%)', fontsize=20)
        ax.set_title(f'Nucleotide Percentage at Each Position (Barcode ID {barcode_id})', fontsize=20)
        ax.set_xticks(positions)
        ax.set_xticklabels([f'{i}' for i in positions], fontsize=14)  # Set the font size of x-tick labels here
        ax.tick_params(axis='x', labelsize=20)  # Set x-tick label size
        ax.tick_params(axis='y', labelsize=20)  # Adjust y-tick label size if needed
        ax.legend(fontsize=20)
        # ax.set_xticklabels([f'Pos {i}' for i in positions])
        ax.legend(fontsize=20)

        # Save the plot
        plt.savefig(self.plot_path + f'{barcode_id}/' + f'plot_stacked_bar_for_barcode_{start_pos}_{end_pos}.png')

        # Show the plot
        plt.tight_layout()
        plt.show()






