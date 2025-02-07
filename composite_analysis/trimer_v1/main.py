from config import build_config
from seq_analysis import SeqAnalysis
from plot import Plot

if __name__ == '__main__':
    # call SeqAnalysis class
    config = build_config()
    if(True):
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
                                   barcode_location=config['barcode_location'],
                                   adapter_end_location=config['adapter_end_location'],
                                   total_sequence_length=config['total_sequence_length'],
                                   total_sequence_length_with_adapters=config['total_sequence_length_with_adapters'],
                                   alphabet=config['alphabet'],
                                   max_bc_distance=config['max_bc_distance'],
                                   combinatorial_letters_length=config['combinatorial_letters_length'],
                                   blast_db_bc=config['blast_db_bc'],
                                   blast_db_bc_identifier=config['blast_db_bc_identifier'],
                                   blast_database_path=config['blast_database_path'],
                                   blast_db_bc_fasta=config['blast_db_bc_fasta'],
                                   blast_db_seq_design_fasta=config['blast_db_seq_design_fasta'],
                                   blast_database_seq_design_path=config['blast_database_seq_design_path'],
                                   blast_database_bc_path=config['blast_database_bc_path'],
                                   general_information_file=config['general_information_file'],
                                   th_minimum_len_reads_to_analyse=config['th_minimum_len_reads_to_analyse'],
                                   csv_path=config['csv_path'],
                                   reads_chunk_to_fasta_file=config['reads_chunk_to_fasta_file'],
                                   fasta_path=config['fasta_path'],
                                   amount_of_bc=config['amount_of_bc'],
                                   blast_bc_results=config['blast_bc_results'],
                                   blast_db_bc_identifier_fasta=config['blast_db_bc_identifier_fasta'],
                                   barcode_length=config['barcode_length'],
                                   blast_bc_identifier_results=config['blast_bc_identifier_results'],
                                   query_results_path=config['query_results_path'],
                                   barcode_with_sequences_distance_dict_file=config['barcode_with_sequences_distance_dict_file'],)
        # get sequences from file
        seq_analysis.run()

    if (True):
        # plot
        plot = Plot(data_path=config['data_path'],
                    output_path=config['output_path'],
                    plot_path=config['plot_path'],
                    sequences_fastq_file=config['sequences_fastq_file'],
                    sequences_file=config['sequences_file'],
                    sequences_unassembled_file=config['sequences.unassembled_file'],
                    sequences_assembled_file=config['sequences.assembled_file'],
                    adapter_start_location=config['adapter_start_location'],
                    combinatorial_location=config['combinatorial_location'],
                    design_file=config['design_file'],
                    barcode_location=config['barcode_location'],
                    adapter_end_location=config['adapter_end_location'],
                    total_sequence_length=config['total_sequence_length'],
                    total_sequence_length_with_adapters=config['total_sequence_length_with_adapters'],
                    alphabet=config['alphabet'],
                    max_bc_distance=config['max_bc_distance'],
                    combinatorial_letters_length=config['combinatorial_letters_length'],
                    csv_path=config['csv_path'],
                    barcode_with_sequences_distance_dict_file=config['barcode_with_sequences_distance_dict_file'],
                    barcode_length=config['barcode_length'],
                    blast_db_bc_fasta=config['blast_db_bc_fasta'],
                    blast_db_seq_design_fasta=config['blast_db_seq_design_fasta'],)
        plot.run()
