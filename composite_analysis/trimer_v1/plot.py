import csv
import itertools
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Union, List, Dict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

import utilities as uts

CHUNK_SIZE = 1000


class Plot():
    def __init__(self, sequences_file: Union[str, Path], plot_path: Union[str, Path], output_path: Union[str, Path],
                 data_path: Union[str, Path], sequences_fastq_file: Union[str, Path],
                 sequences_unassembled_file: Union[str, Path], sequences_assembled_file: Union[str, Path],
                 adapter_start_location: List, combinatorial_location: List, design_file: Union[str, Path],
                 barcode_location: List, adapter_end_location: List,
                 total_sequence_length: int, total_sequence_length_with_adapters: int, alphabet: Dict,
                 max_bc_distance: int, combinatorial_letters_length: int, csv_path: Union[str, Path],
                 barcode_with_sequences_distance_dict_file: Union[str, Path], barcode_length: int, blast_db_bc_fasta: Union[str, Path],
                 blast_db_seq_design_fasta: Union[str, Path]) :
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
        self.barcode_with_sequences_distance_dict_file = barcode_with_sequences_distance_dict_file
        self.barcode_length = barcode_length
        self.blast_db_bc_fasta = blast_db_bc_fasta
        self.blast_db_seq_design_fasta = blast_db_seq_design_fasta

    def run(self):
        self.create_output_dirs()

        # self.plot_sequence_length_distribution()
        # self.plot_sequence_length_distribution(sequence_range=(20, 40))
        # self.plot_all_nucleotide_distribution(input_file=self.barcode_with_sequences_distance_dict_file)
        # ## calculate_nucleotide_percentage_per_bc_file = self.calculate_nucleotide_percentage_per_bc(input_file=self.barcode_with_sequences_distance_dict_file)
        # calculate_nucleotide_percentage_per_bc_file,  calculate_nucleotide_count_per_bc_file= self.calculate_nucleotide_count_and_percentage_per_bc(input_file=f'/barcode_with_sequences_distance_dict.csv')
        # self.calculate_sequence_percentage_design(self.design_file)
        # #
        # # # Example usage
        # # # self.foreach_bc_and_each_comb_letter_analysis_graph(nucleotide_csv=calculate_nucleotide_percentage_per_bc_file) #TODO: uncomment this
        # # # self.foreach_bc_and_each_comb_letter_analysis_graph(nucleotide_csv='output/csv/nucleotide_percentage_per_bc_seq_length_44_distance_0.csv') #TODO: delete this
        # #
        # # self.calculate_foreach_bc_the_max_likelihood_letter_per_position(count_csv=calculate_nucleotide_count_per_bc_file)
        #
        self.plot_error_distribution()
        self.plot_long_deletion_error_distribution(long_d_th=1)
        self.plot_long_deletion_error_distribution(long_d_th=2)
        self.plot_long_deletion_error_distribution(long_d_th=3)
        self.plot_long_deletion_error_distribution(long_d_th=4)
        self.plot_long_deletion_error_distribution(long_d_th=5)
        self.plot_long_deletion_error_distribution(long_d_th=6)
        self.plot_long_deletion_error_distribution(long_d_th=7)
        # self.plot_long_del_error_distribution(long_d_th=1)
        # self.plot_long_del_error_distribution(long_d_th=2)
        # self.plot_long_del_error_distribution(long_d_th=3)
        # self.plot_long_del_error_distribution(long_d_th=4)
        # self.plot_long_del_error_distribution(long_d_th=5)

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
        self.plot_nucleotide_distribution(plot_parts_of_seq=1, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=2, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=3, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=0, stacked=False, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=1, stacked=False, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=2, stacked=False, input_file=input_file)
        self.plot_nucleotide_distribution(plot_parts_of_seq=3, stacked=False, input_file=input_file)

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
        elif plot_parts_of_seq == 1:  # comb letters and bc
            title += 'comb letters and bc'
            for seq_dict in data['sequences']:
                for seq_info, info in seq_dict.items():
                    end_pos = info['end_pos']
                    sequences.append(list(seq_info[:end_pos]))

                # Adjust sequence alignment to the right
            max_length = max(len(seq) for seq in sequences)
            aligned_sequences = []
            for seq in sequences:
                shift = max_length - len(seq)  # Shift to align right
                aligned_sequences.append(['-'] * shift + seq)  # Add '-' padding at the start

            sequences = aligned_sequences  # Replace original sequences

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
                    seq_start_pos = (start_pos - self.combinatorial_letters_length)
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
        fig, ax = plt.subplots(figsize=(40, 10))
        df.plot(kind='bar', stacked=stacked, ax=ax, width=1.0)
        plt.title(f'Nucleotide Distribution for Barcode: {barcode_idx}', fontsize=35)
        plt.xlabel('Position', fontsize=30)
        plt.ylabel('Frequency', fontsize=30)
        if len(start_positions) != 0 and len(end_positions) != 0:
            # Generate x-tick labels based on start and end positions
            xticks = np.arange(max_length)
            xtick_labels = [f"{start_positions[0] + i}" for i in xticks]  # Modify this logic as per the requirement

            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels, rotation=45)
        else:
            plt.xticks(rotation=45, fontsize=25)
            plt.yticks(fontsize=25)
        plt.legend(title='Nucleotide', fontsize=25, title_fontsize=25)
        plt.tight_layout()
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
        SEQUENCE_LENGTH = 31
        PARAMS_STR = 'seq_length_' + str(SEQUENCE_LENGTH) + '_distance_' + str(DISTANCE)

        count_sequences_df = pd.DataFrame(columns=["barcode_id", "count"])

        df = pd.read_csv(self.csv_path + input_file)

        # Filter the dataframe where the 'distance' column is DISTANCE and 'sequence' has SEQUENCE_LENGTH
        df_filtered = df[(df['distance'] == DISTANCE) & (df['sequence'].str.len() == SEQUENCE_LENGTH)]

        # Dictionary to store the results for each barcode
        result_rows_percentage = []
        result_rows_count = []

        # Group the dataframe by 'barcode id'
        grouped = df_filtered.groupby('barcode_id')

        # For each group (barcode id)
        for barcode_id, group in grouped:

            if barcode_id not in count_sequences_df['barcode_id'].values:
                # Add a new row to the DataFrame
                count_sequences_df = pd.concat([
                    count_sequences_df,
                    pd.DataFrame({"barcode_id": [barcode_id], "count": [len(group['sequence'])]})
                ], ignore_index=True)
            else:
                # Update the existing count (if needed)
                count_sequences_df.loc[count_sequences_df['barcode_id'] == barcode_id, 'count'] = len(group['sequence'])

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

        # Save the DataFrame to a CSV file
        csv_file_path = self.csv_path + f'count_sequences_{PARAMS_STR}_{",".join(self.alphabet.keys())}.csv'
        count_sequences_df.to_csv(csv_file_path, index=False)

        # Optionally plot the stacked bar graph for each barcode id (if needed for percentages)
        self.plot_stacked_bars(df=percentage_df_filtered, params_str=PARAMS_STR)

        return output_csv_path_count, output_csv_path_percentage
    def calculate_nucleotide_percentage_per_bc(self, input_file: Union[str, Path]):
        DISTANCE = 0
        SEQUENCE_LENGTH = 31
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
        # barcode_id = 13
        print(f'barcode_id={barcode_id}')

        # Load the sequence from Design.csv
        design_df = pd.read_csv(self.design_file)
        sequence = design_df[design_df['barcode_id'] == barcode_id]['sequence without adapters'].values[0]

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
        fig, ax = plt.subplots(figsize=(45, 20))

        # Plot the bars
        bar_width = 0.8
        A_bars = ax.bar(positions, A_percentages, bar_width, label='A')
        C_bars = ax.bar(positions, C_percentages, bar_width, bottom=A_percentages, label='C')
        G_bars = ax.bar(positions, G_percentages, bar_width, bottom=A_percentages + C_percentages, label='G')
        T_bars = ax.bar(positions, T_percentages, bar_width, bottom=A_percentages + C_percentages + G_percentages,
                        label='T')

        # Add yellow lines based on self.alphabet
        for pos, letter in zip(positions, sequence):
            if letter in self.alphabet:
                props = self.alphabet[letter]
                cumulative_height = 0
                for nucleotide, proportion in props.items():
                    if proportion > 0:
                        height = proportion * 100  # Convert proportion to percentage
                        bar_start = pos - 0.4  # Center the yellow line within the bar width
                        bar_end = pos + 0.4
                        ax.plot([bar_start, bar_end], [cumulative_height + height, cumulative_height + height],
                                color='yellow', linewidth=6)
                        cumulative_height += height

        # Annotate bars with the actual percentages
        for bars, nucleotide in zip([A_bars, C_bars, G_bars, T_bars], ['A', 'C', 'G', 'T']):
            for bar, letter in zip(bars, sequence):
                height = bar.get_height()
                if height > 0:
                    text_color = 'white' if nucleotide in ['A', 'G'] else 'black'
                    ax.annotate(f'{height:.1f}%',
                                xy=(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2),
                                xytext=(0, 25),
                                textcoords="offset points",
                                ha='center', va='center', fontsize=20, rotation=90, color=text_color)

                # Annotate with sequence letter
                ax.annotate(letter,
                            xy=(bar.get_x() + bar.get_width() / 2, 103),
                            ha='center', va='bottom', fontsize=45, color='black')

        # Draw a hollow square over the first 18 bars
        first_bar_position = positions[0] - bar_width / 2
        last_bar_position = positions[self.barcode_length - 1] + bar_width / 2
        rectangle_height = 105
        rectangle = Rectangle((first_bar_position, 0), last_bar_position - first_bar_position, rectangle_height,
                              fill=False, edgecolor='purple', linewidth=13)
        ax.add_patch(rectangle)

        # Place a text label "BC" next to the square
        ax.text(last_bar_position - 0.5 * (last_bar_position - first_bar_position), -5, 'BC', fontsize=45,
                va='top', ha='center', color='purple')

        # Add table for self.alphabet
        cell_text = []
        row_labels = []
        for key, values in self.alphabet.items():
            row_labels.append(key)
            cell_text.append([f'{int(values[nuc] * 100)}' for nuc in ['A', 'C', 'G', 'T']])


        table = ax.table(cellText=cell_text,
                         rowLabels=row_labels,
                         colLabels=['A (%)', 'C (%)', 'G (%)', 'T (%)'],
                         bbox=[1.05, 0.1, 0.2, 0.8])
        table.auto_set_font_size(False)
        table.set_fontsize(30)  # Increased font size for better visibility
        table.scale(0.5, 3)  # Adjusted scale for tighter fit

        # Set labels and title
        ax.set_xlabel('Position in Sequence', fontsize=50)
        ax.set_ylabel('Percentage (%)', fontsize=50)
        ax.set_title(f'Nucleotide Percentage at Each Position (Barcode ID {barcode_id})', fontsize=50, y=1.08)
        ax.set_xticks(positions)
        ax.set_xticklabels([f'{i}' for i in positions], fontsize=14)
        ax.tick_params(axis='x', labelsize=35)
        ax.tick_params(axis='y', labelsize=35)
        ax.legend(fontsize=35)
        ax.set_ylim(0, 103.5)
        ax.set_xlim(positions[0] - 5, positions[-1] + 5)
        plt.subplots_adjust(top=0.85)
        plt.tight_layout()

        # Save the plot
        os.makedirs(os.path.dirname(self.plot_path + f'{barcode_id}/'), exist_ok=True)
        plt.savefig(
            self.plot_path + f'{barcode_id}/' + f'plot_stacked_bar_for_barcode_{start_pos}_{end_pos}_' + params_str + '.png')
        plt.close()


    def calculate_sequence_percentage_design(self, input_file: Union[str, Path]):
        """
        Calculate the nucleotide percentage for each position in the 'sequence without adapters' column based on the
        config['alphabet'] mapping, and generate a stacked bar plot.
        """
        df = pd.read_csv(input_file)

        # Extract the 'sequence without adapters' column
        sequences = df['sequence without adapters']
        barcode_idx = df['barcode_id']

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

        # Load the sequence from Design.csv
        design_df = pd.read_csv(self.design_file)
        sequence = design_df[design_df['barcode_id'] == barcode_idx]['sequence without adapters'].values[0]

        positions = df['Position']
        A_percentages = df['A']
        C_percentages = df['C']
        G_percentages = df['G']
        T_percentages = df['T']

        # Create the plot
        fig, ax = plt.subplots(figsize=(45, 20))
        bar_width = 0.8

        # Plot the stacked bars
        A_bars = ax.bar(positions, A_percentages, bar_width, label='A')
        C_bars = ax.bar(positions, C_percentages, bar_width, bottom=A_percentages, label='C')
        G_bars = ax.bar(positions, G_percentages, bar_width, bottom=A_percentages + C_percentages, label='G')
        T_bars = ax.bar(positions, T_percentages, bar_width, bottom=A_percentages + C_percentages + G_percentages,
                        label='T')

        # Annotate bars with the actual percentages and sequence letters
        for bars, nucleotide in zip([A_bars, C_bars, G_bars, T_bars], ['A', 'C', 'G', 'T']):
            for bar, letter in zip(bars, sequence):
                height = bar.get_height()
                # Annotate with percentage
                if height > 0:
                    text_color = 'white' if nucleotide in ['A', 'G'] else 'black'
                    ax.annotate(f'{height:.1f}%',
                                xy=(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2),
                                xytext=(0, 25),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='center', fontsize=20, rotation=90, color=text_color)
                # Annotate with sequence letter
                ax.annotate(letter,
                            xy=(bar.get_x() + bar.get_width() / 2, 100),  # place at top of plot area
                            ha='center', va='bottom', fontsize=45, color='black')

        # Draw a hollow square over the first 18 bars
        first_bar_position = positions[0] - bar_width / 2  # start slightly before the first bar
        last_bar_position = positions[self.barcode_length-1] + bar_width / 2  # end slightly after the 18th bar
        rectangle_height = 105  # slightly higher than the axis limit to cover full bars
        rectangle = Rectangle((first_bar_position, 0), last_bar_position - first_bar_position, rectangle_height,
                              fill=False, edgecolor='purple', linewidth=13)
        ax.add_patch(rectangle)

        # Place a text label "BC" next to the square
        ax.text(last_bar_position - 0.5 * (last_bar_position - first_bar_position), -5, 'BC', fontsize=45, va='top', ha='center', color='purple')

        # Add table for self.alphabet
        cell_text = []
        row_labels = []
        for key, values in self.alphabet.items():
            row_labels.append(key)
            cell_text.append([f'{int(values[nuc] * 100)}' for nuc in ['A', 'C', 'G', 'T']])

        table = ax.table(cellText=cell_text,
                         rowLabels=row_labels,
                         colLabels=['A (%)', 'C (%)', 'G (%)', 'T (%)'],
                         bbox=[1.05, 0.1, 0.2, 0.8])
        table.auto_set_font_size(False)
        table.set_fontsize(30)  # Increased font size for better visibility
        table.scale(0.5, 3)  # Adjusted scale for tighter fit

        # Set labels and title
        ax.set_xlabel('Position in Sequence', fontsize=50)
        ax.set_ylabel('Percentage (%)', fontsize=50)
        ax.set_title(f'Nucleotide Percentage at Each Position (Barcode ID {barcode_idx})', fontsize=50,  y=1.08)
        ax.set_xticks(positions)
        ax.set_xticklabels([f'{i}' for i in positions], fontsize=14)  # Set the font size of x-tick labels here
        ax.tick_params(axis='x', labelsize=35)  # Set x-tick label size
        ax.tick_params(axis='y', labelsize=35)  # Adjust y-tick label size if needed
        ax.legend(fontsize=30)

        # Set consistent y-axis limits
        ax.set_ylim(0, 100.5)  # This ensures all percentages are within the range 0 to 100
        ax.set_xlim(positions.iloc[0] - 5, positions.iloc[-1] + 5)  # Ensure "BC" is within the visible range
        plt.subplots_adjust(top=0.85)  # Fine-tuned top margin
        plt.tight_layout()


        # Save and show the plot
        output_plot_path = self.plot_path + f'{barcode_idx}/sequence_percentage_distribution_design.png'
        os.makedirs(os.path.dirname(self.plot_path+f'{barcode_idx}/'), exist_ok=True)
        plt.savefig(output_plot_path)
        # plt.show()

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

    def plot_error_distribution(self, distance_method: str = "levenshtein"):
        """
        Reads the annotations_levenshtein.csv file, processes errors (S, D, I),
        and generates three bar plots for substitutions, deletions, and insertions.
        """

        seq_design_iter = uts.open_fasta_yield(self.blast_db_seq_design_fasta)
        seq_design_chunk = list(itertools.islice(seq_design_iter, CHUNK_SIZE))
        
        for seq_design_info in seq_design_chunk:
            seq_design = seq_design_info.seq
            barcode_i = int(seq_design_info.id)

            # Read CSV file
            df = pd.read_csv(self.csv_path+f"/{barcode_i}/annotations_{distance_method}.csv", sep=",")
            print(f"df={df.columns}")
            # Extract column names representing sequence positions
            df.columns = df.columns.str.strip()
            print(f"df.columns={df.columns}")
            sequence_positions = df.columns[3:]  # Skip first three columns: Compared Sequence, start_pos, end_pos
            print(f"sequence_positions={sequence_positions}")

            # Fix duplicate column names using numerical indices
            renamed_positions = [f"{i + 1}{nuc}" for i, nuc in enumerate(seq_design)]  # Use design_seq as labels

            # Rename columns to match sequence positions uniquely
            df.columns = ["Compared Sequence", "start_pos", "end_pos"] + renamed_positions

            # Initialize error counts
            substitution_counts = {pos: 0 for pos in renamed_positions}
            deletion_counts = {pos: 0 for pos in renamed_positions}
            insertion_counts = {pos: 0 for pos in renamed_positions}

            # Count errors in each position
            for pos in renamed_positions:
                substitution_counts[pos] = (df[pos] == "S").sum()
                deletion_counts[pos] = (df[pos] == "D").sum()
                insertion_counts[pos] = (df[pos] == "I").sum()

            # Convert to lists for plotting
            positions = [f"{nuc}" for i, nuc in enumerate(renamed_positions)]
            sub_values = list(substitution_counts.values())
            del_values = list(deletion_counts.values())
            ins_values = list(insertion_counts.values())

            # plt.xticks(list(design_seq))

            # Plot Substitutions
            plt.figure(figsize=(12, 6))
            plt.tight_layout()
            plt.bar(positions, sub_values)
            plt.xlabel("Position (Nucleotide)")
            plt.ylabel("Count of Substitutions (S)")
            plt.title("Distribution of Substitution Errors")
            plt.xticks(rotation=90)
            plt.savefig(self.plot_path + f'{barcode_i}/' + f"{distance_method}_Substitution" + '.png')

            # Plot Deletions
            plt.figure(figsize=(12, 6))
            plt.tight_layout()
            plt.bar(positions, del_values)
            plt.xlabel("Position (Nucleotide)")
            plt.ylabel("Count of Deletions (D)")
            plt.title("Distribution of Deletion Errors")
            plt.xticks(rotation=90)
            plt.savefig(self.plot_path + f'{barcode_i}/' + f"{distance_method}_Deletions" + '.png')

            # Plot Insertions
            plt.figure(figsize=(12, 6))
            plt.tight_layout()
            plt.bar(positions, ins_values)
            plt.xlabel("Position (Nucleotide)")
            plt.ylabel("Count of Insertions (I)")
            plt.title("Distribution of Insertion Errors")
            plt.xticks(rotation=90)
            plt.savefig(self.plot_path + f'{barcode_i}/' + f"{distance_method}_Insertion" + '.png')

    def plot_long_del_error_distribution(self, long_d_th: int, distance_method: str = "levenshtein"):
        """
        Reads the annotations_levenshtein.csv file, processes long deletions (D),
        and generates a histogram of their start positions based on a given threshold.

        - `long_d_th`: Minimum consecutive deletions required to be counted.
        """

        reads_iter = uts.open_fasta_yield(self.blast_db_bc_fasta)
        reads_chunk = list(itertools.islice(reads_iter, CHUNK_SIZE))

        for barcode_info in reads_chunk:
            barcode_i = int(barcode_info.id)
            design_seq = barcode_info.seq

            # Read CSV file
            df = pd.read_csv(self.csv_path + f"/{barcode_i}/annotations_{distance_method}.csv", sep=",")
            df.columns = df.columns.str.strip()  # Remove any whitespace from column names

            # Extract nucleotide positions (skip 'Compared Sequence', 'start_pos', 'end_pos')
            sequence_positions = df.columns[3:]

            # Initialize list to store positions of long deletions
            long_deletion_positions = []

            # Iterate over rows
            for _, row in df.iterrows():
                start_pos = row["start_pos"]  # Base starting position

                deletion_streak = 0  # Track consecutive deletions
                for i, col in enumerate(sequence_positions):
                    if row[col] == "D":  # If a deletion is found
                        deletion_streak += 1
                    else:
                        # If streak meets threshold, store the start position
                        if deletion_streak >= long_d_th:
                            long_deletion_positions.append(start_pos + i - deletion_streak + 1)
                        deletion_streak = 0  # Reset streak

                    # If streak meets threshold, store the start position
                    if deletion_streak >= long_d_th:
                        long_deletion_positions.append(start_pos + i - deletion_streak + 1)
                        deletion_streak = 0  # Reset streak

                # If a streak ends at the last position, check and store it
                if deletion_streak >= long_d_th:
                    long_deletion_positions.append(start_pos + len(sequence_positions) - deletion_streak)

            renamed_positions = [f"{i + 1}{nuc}" for i, nuc in enumerate(design_seq)]  # Use design_seq as labels
            # Plot histogram
            plt.figure(figsize=(10, 5))
            plt.hist(long_deletion_positions,
                     bins=range(min(long_deletion_positions, default=0), max(long_deletion_positions, default=1) + 1),
                     edgecolor='black', alpha=0.7)
            plt.xlabel("Start Position of Long Deletions")
            plt.ylabel("Count of Long Deletions")
            plt.title(f"Histogram of Long Deletions (Threshold = {long_d_th})")
            plt.grid(axis="y", linestyle="--", alpha=0.7)

            # Save the plot
            plt.savefig(self.plot_path + f'{barcode_i}/' + f"{distance_method}_LongDeletions_{long_d_th}.png")
            plt.close()



    def plot_long_deletion_error_distribution(self, long_d_th: int, distance_method: str = "levenshtein"):
        """
        Reads the annotations_levenshtein.csv file, processes long deletions (D),
        and generates a bar plot similar to the other error distributions.

        - `long_d_th`: Minimum consecutive deletions required to be counted as a long deletion.
        """

        seq_design_iter = uts.open_fasta_yield(self.blast_db_seq_design_fasta)
        seq_design_chunk = list(itertools.islice(seq_design_iter, CHUNK_SIZE))

        for seq_design_info in seq_design_chunk:
            seq_design = seq_design_info.seq
            barcode_i = int(seq_design_info.id)

            # Read CSV file
            df = pd.read_csv(self.csv_path + f"/{barcode_i}/annotations_{distance_method}.csv", sep=",")
            df.columns = df.columns.str.strip()  # Remove any whitespace from column names

            # Extract column names representing sequence positions
            sequence_positions = df.columns[3:]  # Skip first three columns: Compared Sequence, start_pos, end_pos

            # Generate x-axis labels (1C, 2C, 3A, ...)
            renamed_positions = [f"{i + 1}{nuc}" for i, nuc in enumerate(seq_design)]

            # Rename DataFrame columns accordingly
            df.columns = ["Compared Sequence", "start_pos", "end_pos"] + renamed_positions

            # Initialize long deletion counts
            long_deletion_counts = {pos: 0 for pos in renamed_positions}

            # Iterate over rows to track long deletions
            for _, row in df.iterrows():
                deletion_streak = 0  # Track consecutive deletions
                for i, col in enumerate(renamed_positions):
                    # if row[col] == "D":  # If a deletion is found
                    #     deletion_streak += 1
                    # else:
                    #     # If streak meets threshold, add to count
                    #     if deletion_streak >= long_d_th:
                    #         long_deletion_counts[renamed_positions[i - deletion_streak]] += 1  # Assign to start of streak
                    #     deletion_streak = 0  # Reset streak
                    x=3
                    if row[col] == "D":  # If a deletion is found
                        deletion_streak += 1
                    else:
                        # If streak meets threshold, add to count
                        if deletion_streak >= long_d_th:
                            long_deletion_counts[renamed_positions[i - deletion_streak]] += 1  # Assign to start of streak
                        deletion_streak = 0  # Reset streak
                    if deletion_streak >= long_d_th:
                        long_deletion_counts[renamed_positions[i]] += 1  # Assign to start of streak
                        deletion_streak = 0  # Reset streak

                # If a streak ends at the last position, check and store it
                if deletion_streak >= long_d_th:
                    long_deletion_counts[renamed_positions[-deletion_streak]] += 1  # Assign to start of streak

            # Convert to lists for plotting
            positions = renamed_positions
            long_del_values = list(long_deletion_counts.values())

            # Plot Long Deletions
            plt.figure(figsize=(12, 6))
            plt.tight_layout()
            plt.bar(positions, long_del_values)
            plt.xlabel("Position (Nucleotide)")
            plt.ylabel(f"Count of Long Deletions (≥{long_d_th})")
            plt.title(f"Distribution of Long Deletions (Threshold = {long_d_th})")
            plt.xticks(rotation=90)
            plt.savefig(self.plot_path + f'{barcode_i}/' + f"{distance_method}_LongDeletions_{long_d_th}.png")
            plt.close()







