from seq_analysis import SeqAnalysis

import os
import numpy as np
from matplotlib import pyplot as plt


class Plot(SeqAnalysis):
    def __init__(self, data_path, output_path, plot_path, sequences_fastq_file, sequences_file, sequences_unassembled_file,
                 sequences_assembled_file, adapter_start_location, combinatorial_location, design_file,
                 identifier_location,
                 barcode_location, adapter_end_location, total_sequence_length, total_sequence_length_with_adapters,
                 alphabet,
                 max_bc_distance, combinatorial_letters_length, identifier_length):
        super().__init__(data_path, output_path, plot_path, sequences_fastq_file, sequences_file, sequences_unassembled_file,
                         sequences_assembled_file,
                         adapter_start_location, combinatorial_location, design_file, identifier_location,
                         barcode_location,
                         adapter_end_location, total_sequence_length, total_sequence_length_with_adapters, alphabet,
                         max_bc_distance, combinatorial_letters_length, identifier_length)

    def run(self):
        self.plot_sequence_length_distribution()
        self.plot_all_nucleotide_distribution()
