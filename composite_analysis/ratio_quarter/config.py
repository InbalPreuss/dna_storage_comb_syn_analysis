# Description: Configuration file for composite_analysis package.
# Sequence info:
# 1-33 adapter
# 34-43 Combinatorial letters (in this case all 10 N are combination of A,C,G and T, with 1:1:1:1 ratio_c_t.)
# 44-49 identifier
# 50-79 barcode
# 80-113 adapter
# TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNAACCTTagttctaacaaccaaaatatgattagcttaCTGTCTCTTATACACATCTCCGAGCCCACGAGAC
# adapter: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG, length: 33
# combinatorial letters: NNNNNNNNNN (10 Ns), length: 10
# identifier: AACCTT, length: 6
# barcode: agttctaacaaccaaaatatgattagctta, length: 30
# adapter: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC, length: 34
# length of the sequence with adapters: 113
# length of the sequence without adapters: 47

def build_config():
    config = {
        'data_path': 'data/',
        'output_path': 'output/',
        'plot_path': 'output/plots/',
        'csv_path': 'output/csv/',
        'fasta_path': 'output/fasta/',

        'general_information_file': 'output/csv/general_information.csv',
        'reads_chunk_to_fasta_file': "output/fasta/reads_chunk_to_fasta",

        'sequences_fastq_file': 'data/merged.assembled.fastq',
        'sequences.unassembled_file': 'data/merged.unassembled.fastq',
        'sequences.assembled_file': 'data/merged.assembled.fastq',
        'design_file': 'data/design.csv',

        'sequences_file': 'output/sequences_file.txt',

        'adapter_start_location': [1,33],
        'combinatorial_location': [34,43],
        'identifier_location': [44,49],
        'barcode_location': [50,79],
        'adapter_end_location': [80,113],
        'alphabet': {'N': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}},
        'amount_of_bc': 10,

        'max_bc_distance': 10,
        'th_minimum_len_reads_to_analyse': 1,

        # Blast database bc
        'blast_database_path': "output/blast_database/",
        'blast_db_bc': "output/blast_database/blast_db_bc",
        'blast_db_bc_identifier': "output/blast_database/blast_db_bc_identifier",
        'blast_db_bc_fasta': "output/blast_database/blast_db_bc.fasta",
        'blast_db_bc_identifier_fasta': "output/blast_database/blast_db_bc_identifier_fasta.fasta",
        'blast_bc_results': "output/blast_database/blast_bc_results.csv",
        'blast_bc_identifier_results': "output/blast_database/blast_bc_identifier_results.csv",
        'query_results_path': "output/blast_database/query_results.txt",

        # Blast database bc + identifier
        'blast_database_path': "output/blast_database/",
        'blast_db_bc': "output/blast_database/blast_db_bc",
        'blast_db_bc_identifier': "output/blast_database/blast_db_bc_identifier",
        'blast_db_bc_fasta': "output/blast_database/blast_db_bc.fasta",
        'blast_db_bc_identifier_fasta': "output/blast_database/blast_db_bc_identifier_fasta.fasta",
        'blast_bc_results': "output/blast_database/blast_bc_results.csv",
        'blast_bc_identifier_results': "output/blast_database/blast_bc_identifier_results.csv",
        'query_results_path': "output/blast_database/query_results.txt",
    }
    config['barcode_length'] = config['barcode_location'][1] - config['barcode_location'][0] + 1 # 30
    config['identifier_length'] = config['identifier_location'][1] - config['identifier_location'][0] + 1  # 6
    config['combinatorial_letters_length'] = config['combinatorial_location'][1] - config['combinatorial_location'][0] + 1 # 10
    config['total_sequence_length'] = config['barcode_length'] + config['combinatorial_letters_length'] + config['identifier_length'] # 46
    config['adapter_start_location_length'] = config['adapter_start_location'][1] - config['adapter_start_location'][0] + 1 # 33
    config['adapter_end_location_length'] = config['adapter_end_location'][1] - config['adapter_end_location'][0] + 1 # 34
    config['total_sequence_length_with_adapters'] = config['total_sequence_length'] + config['adapter_start_location'][1] + config['adapter_end_location'][1] # 113

    return config
