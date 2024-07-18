# Description: Configuration file for composite_analysis package.
# Sequence info:
# 1-33 adapter
# 34-43 Combinatorial letters (in this case all 10 N are combination of A,C,G and T, with 1:1:1:1 ratio.)
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


config = {
    'data_path': 'data/',
    'output_path': 'output/',
    'plot_path': 'output/plots/',

    'sequences_file': 'data/merged.assembled.fastq',
    'sequences.unassembled_file': 'data/merged.unassembled.fastq',
    'sequences.assembled_file': 'data/merged.assembled.fastq',
    'design_file': 'data/design.csv',

    'adapter_start_location': [1,33],
    'combinatorial_location': [34,43],
    'identifier_location': [44,49],
    'barcode_location': [50,79],
    'adapter_end_location': [80,113],
    'total_sequence_length_with_adapters': 113,
    'total_sequence_length': 47,
    'ratio': [1, 1, 1, 1],
    'combinatorial_letters_length': 8,
    'identifier_length': 6,

    'max_bc_distance': 10,

}