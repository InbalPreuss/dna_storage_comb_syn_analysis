import csv
import math
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
                 max_bc_distance: int, combinatorial_letters_length: int, csv_path: Union[str, Path]):
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
        self.csv_path = csv_path

    def run(self):
        # self.create_output_dirs()

        # self.plot_sequence_length_distribution()
        # self.plot_sequence_length_distribution(sequence_range=(0, 48))
        # self.plot_sequence_length_distribution(sequence_range=(108, 108))
        # self.plot_sequence_length_distribution(sequence_range=(45, 45))
        # self.plot_sequence_length_distribution(sequence_range=(44, 44))
        # self.plot_sequence_length_distribution(sequence_range=(43, 43))
        # self.plot_sequence_length_distribution(sequence_range=(42, 42))
        # self.plot_sequence_length_distribution(sequence_range=(41, 41))
        # self.plot_all_nucleotide_distribution(input_file='/barcode_with_sequences_distance_dict.csv')
        # calculate_nucleotide_percentage_per_bc_file = self.calculate_nucleotide_percentage_per_bc(input_file='/barcode_with_sequences_distance_dict.csv')
        calculate_nucleotide_percentage_per_bc_file,  calculate_nucleotide_count_per_bc_file= self.calculate_nucleotide_count_and_percentage_per_bc(input_file='/barcode_with_sequences_distance_dict.csv')
        # self.calculate_sequence_percentage_design(self.design_file)

        # Example usage
        # self.foreach_bc_and_each_comb_letter_analysis_graph(nucleotide_csv=calculate_nucleotide_percentage_per_bc_file) #TODO: uncomment this
        # self.foreach_bc_and_each_comb_letter_analysis_graph(nucleotide_csv='output/csv/nucleotide_percentage_per_bc_seq_length_44_distance_0.csv') #TODO: delete this

        self.calculate_foreach_bc_the_max_likelihood_letter_per_position(count_csv=calculate_nucleotide_count_per_bc_file)

    def create_output_dirs(self):
        os.makedirs(self.plot_path, exist_ok=True)

    # def plot_sequence_length_distribution(self):
    #     sequence_lengths = []
    #
    #     with open(self.sequences_file, "r") as input_handle:
    #         for line in input_handle:
    #             sequence_lengths.append(len(line.strip()))
    #
    #     plt.hist(sequence_lengths, bins=40)
    #     plt.xlabel('Sequence Length')
    #     plt.ylabel('Frequency')
    #     plt.title(f'Sequence Length Distribution, \n Total number of sequences: {len(sequence_lengths)}')
    #     plt.savefig(self.plot_path + 'sequence_length_distribution.png')
    #     plt.show()

    def plot_sequence_length_distribution(self, sequence_range=None):
        """
        Plots the sequence length distribution. If sequence_range is provided, only sequences
        within that length range will be considered for the plot.

        :param sequence_range: A tuple (min_length, max_length) to filter sequences by length.
                               If None, the full distribution is plotted.
        """
        sequence_lengths = []

        # Read the sequence lengths from the file
        with open(self.sequences_file, "r") as input_handle:
            for line in input_handle:
                seq_length = len(line.strip())
                sequence_lengths.append(seq_length)

        # Filter by the specified range if provided
        if sequence_range:
            min_length, max_length = sequence_range
            sequence_lengths = [length for length in sequence_lengths if min_length <= length <= max_length]

        # Plot the histogram
        plt.hist(sequence_lengths, bins=40)
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')

        # Modify the title based on the filtering
        if sequence_range:
            plt.title(
                f'Sequence Length Distribution for Range ({min_length}-{max_length}),\n Total: {len(sequence_lengths)}')
            fig_name = f'sequence_length_distribution_{min_length}_{max_length}.png'
        else:
            plt.title(f'Sequence Length Distribution, \n Total: {len(sequence_lengths)}')
            fig_name = 'sequence_length_distribution.png'

        # Save and display the plot
        plt.savefig(self.plot_path + fig_name)
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
                barcode_idx = row['barcode_id']
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

    def calculate_nucleotide_count_and_percentage_per_bc(self, input_file: Union[str, Path]):
        DISTANCE = 0
        SEQUENCE_LENGTH = 44
        PARAMS_STR = 'seq_length_' + str(SEQUENCE_LENGTH) + '_distance_' + str(DISTANCE)

        df = pd.read_csv(self.output_path + input_file)

        # Filter the dataframe where the 'distance' column is DISTANCE and 'sequence' has SEQUENCE_LENGTH
        df_filtered = df[(df['distance'] == DISTANCE) & (df['sequence'].str.len() == SEQUENCE_LENGTH)]

        # Dictionary to store the results for each barcode
        result_rows_percentage = []
        result_rows_count = []

        # Group the dataframe by 'barcode id'
        grouped = df_filtered.groupby('barcode_id')

        # For each group (barcode id)
        for barcode_id, group in grouped:
            max_len = group['sequence'].apply(len).max()
            nucleotide_count = defaultdict(lambda: [0] * max_len)
            total_sequences = [0] * max_len

            # Process each sequence in the group
            for _, row in group.iterrows():
                sequence = row['sequence']
                start_pos = row['start_pos']

                # Count nucleotides at each position
                for i, nucleotide in enumerate(sequence[start_pos:]):
                    if i < max_len - start_pos:
                        nucleotide_count[nucleotide][i] += 1
                        total_sequences[i] += 1

            # Save the count of each nucleotide at each position
            count_row = []
            percentage_row = []

            for i in range(max_len):
                if total_sequences[i] > 0:
                    for nucleotide in 'ACGT':
                        # Add counts to count_row
                        count_row.append(nucleotide_count[nucleotide][i])
                        # Add percentages to percentage_row
                        percentage_row.append((nucleotide_count[nucleotide][i] / total_sequences[i]) * 100)
                else:
                    # Append 0 for both counts and percentages if no sequence at this position
                    count_row.extend([0, 0, 0, 0])
                    percentage_row.extend([0, 0, 0, 0])

            # Append rows for both count and percentage
            result_rows_count.append([barcode_id] + count_row)
            result_rows_percentage.append([barcode_id] + percentage_row)

        # Dynamically adjust the column names based on the length of the nucleotide count/percentage rows
        num_columns = len(result_rows_count[0]) - 1  # Exclude the barcode_id
        columns_filtered = ['barcode_id'] + [f'{nucleotide}_{i}' for i in range(num_columns // 4) for nucleotide in
                                             'ACGT']

        # Create DataFrames for the filtered results
        count_df_filtered = pd.DataFrame(result_rows_count, columns=columns_filtered)
        percentage_df_filtered = pd.DataFrame(result_rows_percentage, columns=columns_filtered)

        # Save the count CSV file
        output_csv_path_count = self.csv_path + f'nucleotide_count_per_bc_{PARAMS_STR}_{",".join(self.alphabet.keys())}.csv'
        count_df_filtered.to_csv(output_csv_path_count, index=False)

        # Save the percentage CSV file
        output_csv_path_percentage = self.csv_path + f'nucleotide_percentage_per_bc_{PARAMS_STR}_{",".join(self.alphabet.keys())}.csv'
        percentage_df_filtered.to_csv(output_csv_path_percentage, index=False)

        # Optionally plot the stacked bar graph for each barcode id (if needed for percentages)
        self.plot_stacked_bars(df=percentage_df_filtered, params_str=PARAMS_STR)

        return output_csv_path_count, output_csv_path_percentage
    def calculate_nucleotide_percentage_per_bc(self, input_file: Union[str, Path]):
        DISTANCE = 0
        SEQUENCE_LENGTH = 44
        PARAMS_STR = 'seq_length_'+str(SEQUENCE_LENGTH)+'_distance_'+str(DISTANCE)

        df = pd.read_csv(self.output_path + input_file)

        # Filter the dataframe where the 'distance' column is DISTANCE
        df_filtered = df[(df['distance'] == DISTANCE) & (df['sequence'].str.len() == SEQUENCE_LENGTH)]

        # Dictionary to store results for each barcode
        result_rows = []

        # Group the dataframe by 'barcode idx'
        grouped = df_filtered.groupby('barcode_id')

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
        output_csv_path_filtered = self.csv_path + 'nucleotide_percentage_single_row_filtered_fixed_'+PARAMS_STR+'.csv'
        percentage_df_filtered.to_csv(output_csv_path_filtered, index=False)

        # Now, plot the stacked bar graph for each barcode id
        self.plot_stacked_bars(df=percentage_df_filtered, params_str=PARAMS_STR)

        return output_csv_path_filtered

    def plot_stacked_bars(self, df, params_str: str):
        START_POS = 0
        END_POS = 43

        # For each unique barcode id in the dataframe, plot the stacked bars
        for barcode_id in df['barcode_id'].unique():
            self.plot_stacked_bar_for_barcode(df, barcode_id, START_POS, END_POS, params_str)

    def plot_stacked_bar_for_barcode(self, df, barcode_id, start_pos, end_pos, params_str: str):
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
        # Set consistent y-axis limits
        ax.set_ylim(0, 100.5)  # This ensures all percentages are within the range 0 to 100

        # Save the plot
        plt.savefig(self.plot_path + f'{barcode_id}/' + f'plot_stacked_bar_for_barcode_{start_pos}_{end_pos}_'+params_str+'.png')

        # Show the plot
        plt.tight_layout()
        plt.show()

    def calculate_sequence_percentage_design(self, input_file: Union[str, Path]):
        """
        Calculate the nucleotide percentage for each position in the 'sequence without adapters' column based on the
        config['alphabet'] mapping, and generate a stacked bar plot.
        """
        df = pd.read_csv(input_file)

        # Extract the 'sequence without adapters' column
        sequences = df['sequence without adapters']
        barcode_idx = df['Seq Name']

        # Process each sequence
        for sequence, barcode_idx in zip(sequences, barcode_idx):
            # Dictionary to store percentages for each position
            max_len = max(sequences.str.len())  # Find the maximum length of sequences
            nucleotide_count = defaultdict(lambda: {'A': 0, 'C': 0, 'G': 0, 'T': 0})
            total_sequences = [0] * max_len

            for i, char in enumerate(sequence):
                if char in self.alphabet:
                    nucleotide_dist = self.alphabet[char]
                    for nucleotide in 'ACGT':
                        nucleotide_count[i][nucleotide] += nucleotide_dist[nucleotide]
                    total_sequences[i] += 1

            # Calculate the percentage for each nucleotide at each position
            percentage_data = {'Position': [], 'A': [], 'C': [], 'G': [], 'T': []}
            for i in range(max_len):
                if total_sequences[i] > 0:
                    percentage_data['Position'].append(i)
                    for nucleotide in 'ACGT':
                        percentage_data[nucleotide].append((nucleotide_count[i][nucleotide] / total_sequences[i]) * 100)

            # Create DataFrame for percentages
            percentage_df = pd.DataFrame(percentage_data)

            # Save the result to a CSV file
            output_csv_path = self.csv_path + f'/{barcode_idx}/sequence_percentage_distribution_design.csv'
            os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
            percentage_df.to_csv(output_csv_path, index=False)

            # Plot the stacked bar graph
            self.plot_stacked_bars_design(percentage_df, barcode_idx)

    def plot_stacked_bars_design(self, df, barcode_idx):
        """
        Plot a stacked bar graph for nucleotide percentages at each position.
        """
        positions = df['Position']
        A_percentages = df['A']
        C_percentages = df['C']
        G_percentages = df['G']
        T_percentages = df['T']

        # Create the plot
        fig, ax = plt.subplots(figsize=(45, 15))
        bar_width = 0.8

        # Plot the stacked bars
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
                    text_color = 'white' if nucleotide in ['A', 'G'] else 'black'
                    ax.annotate(f'{height:.1f}%',
                                xy=(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2),
                                xytext=(0, 25),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='center', fontsize=20, rotation=90, color=text_color)


        # Set labels and title
        ax.set_xlabel('Position in Sequence', fontsize=20)
        ax.set_ylabel('Percentage (%)', fontsize=20)
        ax.set_title(f'Nucleotide Percentage at Each Position (Barcode ID {barcode_idx})', fontsize=20)
        ax.set_xticks(positions)
        ax.set_xticklabels([f'{i}' for i in positions], fontsize=14)  # Set the font size of x-tick labels here
        ax.tick_params(axis='x', labelsize=20)  # Set x-tick label size
        ax.tick_params(axis='y', labelsize=20)  # Adjust y-tick label size if needed
        ax.legend(fontsize=20)
        # ax.set_xticklabels([f'Pos {i}' for i in positions])
        ax.legend(fontsize=20)
        # Set consistent y-axis limits
        ax.set_ylim(0, 100.5)  # This ensures all percentages are within the range 0 to 100

        # Save and show the plot
        output_plot_path = self.plot_path + f'{barcode_idx}/sequence_percentage_distribution_design.png'
        plt.savefig(output_plot_path)
        plt.tight_layout()
        plt.show()

    import pandas as pd

    def foreach_bc_and_each_comb_letter_analysis_graph(self, nucleotide_csv,
                                                       combinatorial_letters=['L', 'N', 'I', 'B']):
        # Load the design CSV and nucleotide percentage CSV
        design_df = pd.read_csv(self.design_file)
        nucleotide_df = pd.read_csv(nucleotide_csv)

        # Iterate over each row in the Design CSV (Seq Name and combinatorial letters column)
        for index, row in design_df.iterrows():
            seq_name = row['barcode_id']
            # combinatorial_seq = row['combinatorial letters']
            combinatorial_seq = row['sequence without adapters']
            nucleotide_df_barcode_id = nucleotide_df[nucleotide_df['barcode_id'] == seq_name]

            # Create a dictionary to store the positions and percentages for each combinatorial letter found
            letter_percentages = {letter: {'positions': [], 'percentages': []} for letter in combinatorial_letters}

            # Iterate over each letter in the combinatorial sequence
            for pos, letter in enumerate(combinatorial_seq):
                if letter in combinatorial_letters:
                    # For each letter (A, C, G, T), extract its percentage from the nucleotide_df
                    nucleotide_columns = [f"{base}_{pos}" for base in ['A', 'C', 'G', 'T']]
                    nucleotide_percentages = nucleotide_df_barcode_id[nucleotide_columns].iloc[0].values

                    # Store the percentages
                    letter_percentages[letter]['positions'].append(pos)
                    letter_percentages[letter]['percentages'].append(nucleotide_percentages)

                    # Save results to CSV for each combinatorial letter
                    for letter, data in letter_percentages.items():
                        if data['positions']:
                            # Create a DataFrame with positions and percentages for each base
                            percentage_df = pd.DataFrame(
                                data['percentages'],
                                columns=[f"A", f"C", f"G", f"T"]
                            )
                            percentage_df['Position'] = data['positions']

                            # Save the DataFrame to CSV for this combinatorial letter
                            percentage_df.to_csv(f"{self.csv_path}{seq_name}/{letter}_percentages.csv", index=False)

                            # Now generate the graph for each letter
                            plt.figure(figsize=(8, 6))
                            plt.plot(percentage_df['Position'], percentage_df['C'], label='C', color='blue', marker='o',
                                     linestyle='-')
                            plt.plot(percentage_df['Position'], percentage_df['T'], label='T', color='orange',
                                     marker='o', linestyle='-')

                            # Add horizontal lines at each data point for C and T
                            plt.axhline(y=self.alphabet[letter]['C']*100, color='blue', linestyle='--', linewidth=0.5)  # Horizontal lines for C
                            plt.axhline(y=self.alphabet[letter]['T']*100, color='orange', linestyle='--', linewidth=0.5)  # Horizontal lines for C
                            # for y in percentage_df['C']:
                            #     plt.axhline(y=y, color='blue', linestyle='--', linewidth=0.5)  # Horizontal lines for C
                            # for y in percentage_df['T']:
                            #     plt.axhline(y=y, color='orange', linestyle='--',
                            #                 linewidth=0.5)  # Horizontal lines for T

                            # Add labels and title
                            plt.xlabel('Position')
                            plt.ylabel('Percentage')
                            plt.ylim(0, 100)
                            plt.yticks(list(range(0, 100, 10)))
                            plt.title(f"{letter} - {seq_name}, C:{int(self.alphabet[letter]['C']*100)}%, T:{int(self.alphabet[letter]['T']*100)}%")
                            plt.legend()

                            # Save the plot
                            graph_path = os.path.join(self.plot_path, f"{seq_name}/{letter}_graph.png")
                            plt.savefig(graph_path)
                            plt.close()


    # Function to adjust probabilities with epsilon and normalize
    def adjust_probs_with_epsilon(self, probs):
        # Define epsilon value
        eps = 0.1

        # Add epsilon to non-zero probabilities
        adjusted_probs = {nuc: prob + eps if prob > 0 else eps for nuc, prob in probs.items()}
        # Calculate normalization factor
        normalization_factor = sum(adjusted_probs.values())
        # Normalize probabilities
        normalized_probs = {nuc: adjusted_probs[nuc] / normalization_factor for nuc in adjusted_probs}
        return normalized_probs

    def calculate_foreach_bc_the_max_likelihood_letter_per_position(self, count_csv: str):
        # Load the nucleotide count CSV
        df = pd.read_csv(count_csv)

        # List of nucleotides
        nucleotides = ['A', 'C', 'G', 'T']

        # List to store the result for each barcode_id and position
        result_rows = []
        score_rows = []  # List to store the score for each alphabet letter at each position for each barcode_id

        # Iterate over each barcode_id in the count CSV
        for index, row in df.iterrows():
            barcode_id = row['barcode_id']
            max_len = (len(row) - 1) // 4  # Exclude barcode_id, and divide by 4 for ACGT columns

            # Store the selected σ for this barcode
            selected_sigmas = []
            position_scores = {'barcode_id': barcode_id}  # Dictionary to store scores for each position

            # For each position in the sequence
            for pos in range(max_len):
                # Extract the counts for A, C, G, T at the current position
                counts = {nuc: row[f"{nuc}_{pos}"] for nuc in nucleotides}

                # Calculate the score for each sigma (from the adjusted probability dictionary)
                max_score = -np.inf
                best_sigma = None
                scores_for_position = {}

                for sigma, probs in self.alphabet.items():
                    # Adjust probabilities using epsilon
                    adjusted_probs = self.adjust_probs_with_epsilon(probs)
                    score = 0

                    # Compute Σ_(i∈{A,C,G,T}) n_i log(p_i) with adjusted probabilities
                    for nuc in nucleotides:
                        n_i = counts[nuc]
                        p_i = adjusted_probs[nuc]

                        # Avoid log(0) by skipping terms where p_i is 0
                        if p_i > 0 and n_i > 0:
                            score += n_i * math.log(p_i)

                    # Update the best sigma if the current score is higher
                    if score > max_score:
                        max_score = score
                        best_sigma = sigma

                    # Store the score for the current sigma at this position
                    scores_for_position[sigma] = score

                # Append the best sigma for this position
                selected_sigmas.append(best_sigma)

                # Save scores for each sigma at this position
                for sigma, score in scores_for_position.items():
                    position_scores[f"{sigma}_{pos}"] = score

            # Append the results for this barcode_id
            result_rows.append([barcode_id] + selected_sigmas)
            score_rows.append(position_scores)  # Append the score row for the current barcode_id

        # Generate dynamic column names based on the number of positions
        columns = ['barcode_id'] + [f"Position_{i}" for i in range(max_len)]

        # Create a DataFrame to store the results
        result_df = pd.DataFrame(result_rows, columns=columns)
        score_df = pd.DataFrame(score_rows)

        # Save the results and scores to new CSVs
        result_output_csv_path = self.csv_path + f'selected_sigma_per_position_{",".join(self.alphabet.keys())}.csv'
        scores_output_csv_path = self.csv_path + f'sigma_scores_per_position_{",".join(self.alphabet.keys())}.csv'

        result_df.to_csv(result_output_csv_path, index=False)
        score_df.to_csv(scores_output_csv_path, index=False)

        return result_output_csv_path, scores_output_csv_path







