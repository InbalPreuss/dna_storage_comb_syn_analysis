from config import build_config
from seq_analysis import SeqAnalysis

if __name__ == '__main__':
    # call SeqAnalysis class
    config = build_config()
    seq_analysis = SeqAnalysis(data_path=config['data_path'],
                               output_path=config['output_path'],
                               plot_path=config['plot_path'],
                               sequences_fastq_file=config['sequences_fastq_file'],
                               sequences_file=config['sequences_file'],
                               sequences_unassembled_file=config['sequences.unassembled_file'],
                               sequences_assembled_file=config['sequences.assembled_file'],
                               adapter_start_location=config['adapter_start_location'],
                               combinatorial_location=config['combinatorial_location'],
                               design_file=config['design_file'],
                               identifier_location=config['identifier_location'],
                               barcode_location=config['barcode_location'],
                               adapter_end_location=config['adapter_end_location'],
                               total_sequence_length=config['total_sequence_length'],
                               total_sequence_length_with_adapters=config['total_sequence_length_with_adapters'],
                               alphabet=config['alphabet'],
                               max_bc_distance=config['max_bc_distance'],
                               combinatorial_letters_length=config['combinatorial_letters_length'],
                               identifier_length=config['identifier_length'])
    # get sequences from file
    seq_analysis.run()

