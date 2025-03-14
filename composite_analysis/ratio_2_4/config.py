

def build_config():
    config = {
        'data_path': 'data/',
        'output_path': 'output/',
        'plot_path': 'output/plots/',
        'csv_path': 'output/csv/',
        'fasta_path': 'output/fasta/',

        'general_information_file': 'output/csv/general_information.csv',
        'reads_chunk_to_fasta_file': "output/fasta/reads_chunk_to_fasta",

        'sequences_fastq_file': 'data/315_merged_output.fastq',
        'sequences.unassembled_file': 'data/315_unassembled_forward.fastq',
        'sequences.assembled_file': 'data/315_unassembled_reverse.fastq',
        'design_file': 'data/design.csv',

        'sequences_file': 'output/sequences_file.txt',

        'adapter_start_location': [1,33], # 33
        'barcode_location': [34,51], # 18
        'combinatorial_location': [52, 64], # 13
        'adapter_end_location': [65,98], # 34
        # 'alphabet': {'M': {'A': 0.5, 'C': 0.5, 'G': 0, 'T': 0},
        #              'K': {'A': 0, 'C': 0, 'G': 0.5, 'T': 0.5},
        #              'Y': {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5},
        #              'R': {'A': 0.5, 'C': 0, 'G': 0.5, 'T': 0},
        #              'W': {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5},
        #              'S': {'A': 0, 'C': 0.5, 'G': 0.5, 'T': 0},
        #              'V': {'A': 0.33, 'C': 0.33, 'G': 0.33, 'T': 0},
        #              'B': {'A': 0, 'C': 0.33, 'G': 0.33, 'T': 0.33},
        #              'D': {'A': 0.33, 'C': 0, 'G': 0.33, 'T': 0.33},
        #              'H': {'A': 0.33, 'C': 0.33, 'G': 0, 'T': 0.33},
        #              '4': {'A': 0, 'C': 0.6, 'G': 0, 'T': 0.4},
        #              '3': {'A': 0, 'C': 0.7, 'G': 0, 'T': 0.3},
        #              '2': {'A': 0, 'C': 0.8, 'G': 0, 'T': 0.2},
        #              '1': {'A': 0, 'C': 0.9, 'G': 0, 'T': 0.1},
        #              'N': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
        #              'a': {'A': 1, 'C': 0, 'G': 0, 'T': 0},
        #              'c': {'A': 0, 'C': 1, 'G': 0, 'T': 0},
        #              'g': {'A': 0, 'C': 0, 'G': 1, 'T': 0},
        #              't': {'A': 0, 'C': 0, 'G': 0, 'T': 1},},
        'alphabet': {'M': {'A': 0.5, 'C': 0.5, 'G': 0, 'T': 0},
                     'K': {'A': 0, 'C': 0, 'G': 0.5, 'T': 0.5},
                     'Y': {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5},
                     'R': {'A': 0.5, 'C': 0, 'G': 0.5, 'T': 0},
                     'W': {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5},
                     'S': {'A': 0, 'C': 0.5, 'G': 0.5, 'T': 0},
                     'V': {'A': 0.33, 'C': 0.33, 'G': 0.33, 'T': 0},
                     'B': {'A': 0, 'C': 0.33, 'G': 0.33, 'T': 0.33},
                     'D': {'A': 0.33, 'C': 0, 'G': 0.33, 'T': 0.33},
                     'H': {'A': 0.33, 'C': 0.33, 'G': 0, 'T': 0.33},
                     '4': {'A': 0, 'C': 0.6, 'G': 0, 'T': 0.4},
                     '3': {'A': 0, 'C': 0.7, 'G': 0, 'T': 0.3},
                     '2': {'A': 0, 'C': 0.8, 'G': 0, 'T': 0.2},
                     '1': {'A': 0, 'C': 0.9, 'G': 0, 'T': 0.1},
                     'N': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
                     'A': {'A': 1, 'C': 0, 'G': 0, 'T': 0},
                     'C': {'A': 0, 'C': 1, 'G': 0, 'T': 0},
                     'G': {'A': 0, 'C': 0, 'G': 1, 'T': 0},
                     'T': {'A': 0, 'C': 0, 'G': 0, 'T': 1}, },

        'amount_of_bc': 42,

        'max_bc_distance': 2,
        'th_minimum_len_reads_to_analyse': 1,

        # CSV output
        'barcode_with_sequences_distance_dict_file': 'csv/barcode_with_sequences_distance_dict.csv',


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
    config['barcode_length'] = config['barcode_location'][1] - config['barcode_location'][0] + 1 # 18
    config['combinatorial_letters_length'] = config['combinatorial_location'][1] - config['combinatorial_location'][0] + 1 # 13
    config['total_sequence_length'] = config['barcode_length'] + config['combinatorial_letters_length'] # 31
    config['adapter_start_location_length'] = config['adapter_start_location'][1] - config['adapter_start_location'][0] + 1 # 33
    config['adapter_end_location_length'] = config['adapter_end_location'][1] - config['adapter_end_location'][0] + 1 # 34
    config['total_sequence_length_with_adapters'] = config['total_sequence_length'] + config['adapter_start_location'][1] + config['adapter_end_location'][1] # 98

    return config
