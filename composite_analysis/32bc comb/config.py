

def build_config():
    config = {
        'data_path': 'data/',
        'output_path': 'output/',
        'plot_path': 'output/plots/',
        'csv_path': 'output/csv/',
        'fasta_path': 'output/fasta/',

        'general_information_file': 'output/csv/general_information.csv',
        'reads_chunk_to_fasta_file': "output/fasta/reads_chunk_to_fasta",

        'sequences_fastq_file': 'data/merged_output.fastq',
        'sequences.unassembled_file': 'data/merged.unassembled.fastq',
        'sequences.assembled_file': 'data/merged.assembled.fastq',
        'design_file': 'data/design.csv',

        'sequences_file': 'output/sequences_file.txt',

        'adapter_start_location': [1,20],
        'barcode_location': [21,28],
        'combinatorial_location': [29, 64],
        'adapter_end_location': [65,85],
        'alphabet': {'0': {'A': 0.5, 'C': 0.5, 'G': 0, 'T': 0},
                     '1': {'A': 0, 'C': 0, 'G': 0.5, 'T': 0.5},
                     '2': {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5},
                     '3': {'A': 0.5, 'C': 0, 'G': 0.5, 'T': 0},
                     '4': {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5},
                     '5': {'A': 0, 'C': 0.5, 'G': 0.5, 'T': 0},
                     '6': {'A': 0.33, 'C': 0.33, 'G': 0.33, 'T': 0},
                     '7': {'A': 0, 'C': 0.33, 'G': 0.33, 'T': 0.33},
                     '8': {'A': 0.33, 'C': 0, 'G': 0.33, 'T': 0.33},
                     '9': {'A': 0.33, 'C': 0.33, 'G': 0, 'T': 0.33},
                     'I': {'A': 0, 'C': 0.6, 'G': 0, 'T': 0.4},
                     'N': {'A': 0, 'C': 0.7, 'G': 0, 'T': 0.3},
                     'B': {'A': 0, 'C': 0.8, 'G': 0, 'T': 0.2},
                     'L': {'A': 0, 'C': 0.9, 'G': 0, 'T': 0.1},
                     'M': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},},
        'amount_of_bc': 42,

        'max_bc_distance': 2,
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
    config['barcode_length'] = config['barcode_location'][1] - config['barcode_location'][0] + 1 # 8
    config['combinatorial_letters_length'] = config['combinatorial_location'][1] - config['combinatorial_location'][0] + 1 # 36
    config['total_sequence_length'] = config['barcode_length'] + config['combinatorial_letters_length'] # 44
    config['adapter_start_location_length'] = config['adapter_start_location'][1] - config['adapter_start_location'][0] + 1 # 20
    config['adapter_end_location_length'] = config['adapter_end_location'][1] - config['adapter_end_location'][0] + 1 # 21
    config['total_sequence_length_with_adapters'] = config['total_sequence_length'] + config['adapter_start_location'][1] + config['adapter_end_location'][1] # 85

    return config
