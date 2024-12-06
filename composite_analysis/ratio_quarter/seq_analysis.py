import collections
import csv
import itertools
import math
import os
import subprocess
from pathlib import Path
from typing import Dict, Union, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO, Align
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import utilities as uts

CHUNK_SIZE = 1000


class SeqAnalysis:
    def __init__(self, data_path: Union[str, Path], output_path: Union[str, Path], plot_path: Union[str, Path],
                 sequences_fastq_file: Union[str, Path], sequences_file: Union[str, Path],
                 sequences_unassembled_file: Union[str, Path],
                 sequences_assembled_file: Union[str, Path],
                 adapter_start_location: List, combinatorial_location: List, design_file: Union[str, Path],
                 identifier_location: List, barcode_location: List,
                 adapter_end_location: List, total_sequence_length: int, total_sequence_length_with_adapters: int,
                 alphabet: Dict, max_bc_distance: int,
                 combinatorial_letters_length: int, identifier_length: int,
                 blast_db_bc: Union[str, Path], blast_db_bc_identifier: Union[str, Path], blast_database_path: Union[str, Path], blast_db_bc_fasta: Union[str, Path],
                 general_information_file: Union[str, Path], th_minimum_len_reads_to_analyse: int,
                 csv_path: Union[str, Path],
                 reads_chunk_to_fasta_file: Union[str, Path], fasta_path: Union[str, Path], amount_of_bc: int,
                 blast_bc_results: Union[str, Path], blast_db_bc_identifier_fasta: Union[str, Path],
                 barcode_length: int, blast_bc_identifier_results: Union[str, Path], query_results_path: Union[str, Path]):
        self.data_path = data_path
        self.output_path = output_path
        self.plot_path = plot_path
        self.csv_path = csv_path
        self.sequences_fastq_file = sequences_fastq_file
        self.sequences_file = sequences_file
        self.sequences_unassembled_file = sequences_unassembled_file
        self.sequences_assembled_file = sequences_assembled_file
        self.adapter_start_location = adapter_start_location
        self.combinatorial_location = combinatorial_location
        self.design_file = design_file
        self.identifier_location = identifier_location
        self.barcode_location = barcode_location
        self.adapter_end_location = adapter_end_location
        self.total_sequence_length = total_sequence_length
        self.total_sequence_length_with_adapters = total_sequence_length_with_adapters
        self.alphabet = alphabet
        self.max_bc_distance = max_bc_distance
        self.combinatorial_letters_length = combinatorial_letters_length
        self.identifier_length = identifier_length
        self.blast_db_bc = blast_db_bc
        self.blast_db_bc_identifier = blast_db_bc_identifier
        self.blast_database_path = blast_database_path
        self.blast_db_bc_fasta = blast_db_bc_fasta
        self.general_information_file = general_information_file
        self.th_minimum_len_reads_to_analyse = th_minimum_len_reads_to_analyse
        self.reads_chunk_to_fasta_file = reads_chunk_to_fasta_file
        self.fasta_path = fasta_path
        self.amount_of_bc = amount_of_bc
        self.blast_bc_results = blast_bc_results
        self.blast_db_bc_identifier_fasta = blast_db_bc_identifier_fasta
        self.barcode_length = barcode_length
        self.blast_bc_identifier_results = blast_bc_identifier_results
        self.query_results_path = query_results_path

        self.create_output_dirs()

    def run(self):
        self.get_sequences_from_file()
        self.get_part_of_seq(seq_part_list=['barcode'], blast_db_fasta=self.blast_db_bc_fasta)
        # self.find_seq_per_barcode_identifier_using_blast(fasta_file=self.blast_db_bc_fasta, blast_db=self.blast_db_bc, query_results_path=self.query_results_path)

        # Levenstein distance
        self.find_seqs_per_barcodes_using(blast_db_fasta=self.blast_db_bc_fasta,
                                                   distance_method='levenshtein',
                                                   output_file='/barcode_with_sequences_distance_dict.csv')
        # self.find_seqs_per_barcodes_using(blast_db_fasta=self.blast_db_bc_fasta,
        #                                   distance_method='alignment',
        #                                   output_file='/barcode_with_sequences_alignment_dict.csv')

        # # find seq using the bc+identifier
        # self.get_part_of_seq(seq_part_list=['identifier', 'barcode'], blast_db_fasta=self.blast_db_bc_identifier_fasta)
        # self.find_seq_per_barcode_identifier_using_blast(fasta_file=self.blast_db_bc_identifier_fasta, blast_db=self.blast_db_bc_identifier, query_results_path=self.query_results_path)

    def create_output_dirs(self):
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.blast_database_path, exist_ok=True)
        os.makedirs(self.csv_path, exist_ok=True)
        os.makedirs(self.fasta_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)
        # os.makedirs(self.blast_db, exist_ok=True)

    def get_sequences_from_file(self):
        with open(self.sequences_fastq_file, "r") as handle, open(self.sequences_file, "w") as output_handle:
            for record in SeqIO.parse(handle, "fastq"):
                output_handle.write(str(record.seq) + "\n")

    def get_barcodes(self):
        with open(self.blast_db_bc_fasta, 'w') as fasta_file:
            chunk_iter = pd.read_csv(self.design_file, chunksize=CHUNK_SIZE)
            for chunk in chunk_iter:
                chunk['barcode'] = chunk['barcode'].str.upper()

                for index, barcode in chunk['barcode'].items():
                    fasta_file.write(f">{index}\n{barcode}\n")

    def get_part_of_seq(self, seq_part_list: List[str], blast_db_fasta: Union[str, Path]):
        with open(blast_db_fasta, 'w') as fasta_file:
            chunk_iter = pd.read_csv(self.design_file, chunksize=CHUNK_SIZE)
            for chunk in chunk_iter:
                chunk['barcode'] = chunk['barcode'].str.upper()

                for index, row in chunk.iterrows():
                    seq_part_from_row = ''
                    for seq_part_i in seq_part_list:
                        seq_part_from_row += row[seq_part_i]
                    fasta_file.write(f">{index}\n{seq_part_from_row}\n")

    def find_seqs_per_barcodes_using(self, blast_db_fasta: Union[str, Path], distance_method: str, output_file: Union[str, Path]) -> None:
        output_csv_file = self.output_path + output_file

        with open(output_csv_file, "w", newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["barcode idx", "barcode", "sequence", "distance", "start_pos", "end_pos", "number_of_sequences"])

            reads_iter = uts.open_fasta_yield(blast_db_fasta)
            reads_chunk = list(itertools.islice(reads_iter, CHUNK_SIZE))

            for barcode_i, barcode_info in enumerate(reads_chunk):
                barcode = barcode_info.seq
                sequence_count = 0
                matching_sequences = []
                scores = []

                with open(self.sequences_file, "r") as input_handle:
                    for seq in input_handle:
                        seq = seq.strip()
                        if distance_method == 'levenshtein':
                            distance, start_pos, end_pos = uts.is_within_levenshtein_distance(seq=barcode.upper(),
                                                                                              target_seq=seq,
                                                                                              max_distance=self.max_bc_distance)
                            if (distance <= self.max_bc_distance) and (start_pos != math.inf):
                                matching_sequences.append((seq, distance, start_pos, end_pos))
                                sequence_count += 1
                        elif distance_method == 'alignment':
                            score, start_pos, end_pos = uts.find_best_alignment(seq=barcode.upper(), target_seq=seq, output_file=self.output_path + '/alignment.txt')
                            if score >= 20:
                                scores.append(score)
                                matching_sequences.append((seq, score, start_pos, end_pos))
                                sequence_count += 1

                for seq, distance, start_pos, end_pos in matching_sequences:
                    csv_writer.writerow([barcode_i, barcode, seq, distance, start_pos, end_pos, sequence_count])
                # plot histogram of scores
                plt.hist(scores)
                plt.xlabel('Score')
                plt.ylabel('Frequency')
                plt.title(f'Bc={barcode_i+1}, Score Distribution, \n Total number of scores: {len(scores)}')
                plt.savefig(self.plot_path + f'bc={barcode_i+1}_score_distribution.png')
                plt.show()
                plt.close()

    def create_blast_database(self, fasta_file: Union[str, Path], blast_db: Union[str, Path]):
        makeblastdb_cmd = f"makeblastdb -in {fasta_file} -dbtype nucl -out {blast_db}"
        try:
            subprocess.run(makeblastdb_cmd, shell=True, check=True)
            print(f"BLAST database '{blast_db}' created successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error creating BLAST database: {e}")

    def extract_bc_and_reads_results_to_csv(self, blast_db: Union[str, Path], word_size: int):

        chunk_size = 1000
        chunk_end = chunk_size
        reads_chunk_i = 1
        res = list()
        failed = 0
        number_of_reads = 0
        is_run_loop = True
        reads_length_counts = {'length': 'count'}

        number_of_reads += uts.get_amount_of_reads_from_file(file_path=self.sequences_fastq_file)
        uts.write_list_to_csv(['# reads', number_of_reads], self.general_information_file)

        reads_iter = uts.open_fastq_yield(self.sequences_fastq_file)

        ''' Usning blastn to analyze the reads '''
        while is_run_loop:
            if number_of_reads < chunk_end:
                chunk_size = number_of_reads % chunk_size
                is_run_loop = False
            print(chunk_size)
            reads_chunk = list(itertools.islice(reads_iter, chunk_size))
            reads_specific_len, reads_length_counts = self.retrieve_reads_in_specific_len_at_least_dict(
                reads=reads_chunk,
                length=self.th_minimum_len_reads_to_analyse, reads_length_counts=reads_length_counts)
            if len(reads_specific_len) == 0:
                if not is_run_loop:
                    break
                chunk_end += chunk_size
                reads_chunk_i += 1
                continue
            fasta_chunk_idx = chunk_end - (chunk_end - (reads_chunk_i * chunk_size)) + chunk_size
            output_fasta_file = self.reads_chunk_to_fasta_format(reads=reads_specific_len, idx=fasta_chunk_idx)

            query_results = self.run_blastn_to_find_location_of_oligo_in_reads_chunks(reads_chunk_i=reads_chunk_i,
                                                                                      output_fasta_file=output_fasta_file,
                                                                                      blast_db=blast_db,
                                                                                      word_size=word_size)

            for read_idx, (read_id, bc_pos_in_read) in enumerate(query_results.items()):
                read = reads_specific_len[read_id]
                try:
                    bc_pos_dict = self.find_pos_in_read_list(bc_pos_in_read=bc_pos_in_read)
                except:
                    print(f'read_id={read_id}')
                    print(read)
                    continue

                res.append(self.identify_oligo(bc_pos_dict=bc_pos_dict,
                                               read=read.seq))
                if res[-1].__contains__(0):
                    failed += 1

            uts.write_list_to_csv(
                ['processed reads:', chunk_end, 'failed reads:', failed, str(100 * failed / (chunk_end + 1)) + '%',
                 'failed.', 'query_results: ', len(query_results)],
                self.general_information_file)
            with open(self.blast_bc_results, "ab") as f:
                np.savetxt(f, res, fmt='%i', delimiter=",")
            res = list()

            if not is_run_loop:
                break
            chunk_end += chunk_size
            reads_chunk_i += 1

        uts.write_dict_to_csv_as_dict(reads_length_counts, self.count_reads_len_file)
        self.hist_length_counts_reads()

    def extract_bc_identifier_and_reads_results_to_csv(self, word_size: int, blast_db: Union[str, Path], blast_results: Union[str, Path], query_results_path: Union[str, Path]):

        chunk_size = 1000
        chunk_end = chunk_size
        reads_chunk_i = 1
        res = list()
        failed = 0
        number_of_reads = 0
        is_run_loop = True
        reads_length_counts = {'length': 'count'}

        number_of_reads += uts.get_amount_of_reads_from_file(file_path=self.sequences_fastq_file)
        uts.write_list_to_csv(['# reads', number_of_reads], self.general_information_file)

        reads_iter = uts.open_fastq_yield(self.sequences_fastq_file)

        ''' Usning blastn to analyze the reads '''
        while is_run_loop:
            if number_of_reads < chunk_end:
                chunk_size = number_of_reads % chunk_size
                is_run_loop = False
            print(chunk_size)
            reads_chunk = list(itertools.islice(reads_iter, chunk_size))
            reads_specific_len, reads_length_counts = self.retrieve_reads_in_specific_len_at_least_dict(
                reads=reads_chunk,
                length=self.th_minimum_len_reads_to_analyse, reads_length_counts=reads_length_counts)
            if len(reads_specific_len) == 0:
                if not is_run_loop:
                    break
                chunk_end += chunk_size
                reads_chunk_i += 1
                continue
            fasta_chunk_idx = chunk_end - (chunk_end - (reads_chunk_i * chunk_size)) + chunk_size
            output_fasta_file = self.reads_chunk_to_fasta_format(reads=reads_specific_len, idx=fasta_chunk_idx)

            query_results = self.run_blastn_to_find_location_of_oligo_in_reads_chunks(reads_chunk_i=reads_chunk_i,
                                                                                      output_fasta_file=output_fasta_file,
                                                                                      blast_db=blast_db,
                                                                                      word_size=word_size)

            for read_idx, (read_id, bc_pos_in_read) in enumerate(query_results.items()):
                read = reads_specific_len[read_id]
                try:
                    bc_idx, bc_pos = self.find_pos_in_read_list(bc_pos_in_read=bc_pos_in_read)
                except:
                    print(f'read_id={read_id}')
                    print(read)
                    continue

                res = [bc_idx, self.identify_oligo(bc_pos=bc_pos,
                                               read=read.seq)]
                uts.write_list_to_csv(res, blast_results)

                self.write_query_results_to_file(query_results=query_results, file_name=query_results_path)

            uts.write_list_to_csv(
                ['processed reads:', chunk_end, 'failed reads:', failed, str(100 * failed / (chunk_end + 1)) + '%',
                 'failed.', 'query_results: ', len(query_results)],
                self.general_information_file)


            if not is_run_loop:
                break
            chunk_end += chunk_size
            reads_chunk_i += 1

        # uts.write_dict_to_csv_as_dict(reads_length_counts, self.count_reads_len_file)
        # self.hist_length_counts_reads()

    def retrieve_reads_in_specific_len_at_least_dict(self, reads: List[SeqIO.SeqRecord], length: int,
                                                     reads_length_counts: Dict) -> Tuple[Dict, Dict]:
        reads_specific_len = {}
        for r in reads:
            read_length = len(r)
            if read_length >= length:
                reads_specific_len[r.id] = r

                if read_length in reads_length_counts:
                    reads_length_counts[read_length] += 1
                else:
                    reads_length_counts[read_length] = 1
        print(f'{len(reads_specific_len)} reads in len at least {length}')

        return reads_specific_len, reads_length_counts

    def find_seq_per_barcode_using_blast(self, fasta_file: Union[str, Path], blast_db: Union[str, Path], word_size: int):
        self.create_blast_database(fasta_file=fasta_file, blast_db=blast_db)
        self.extract_bc_and_reads_results_to_csv(blast_db=blast_db, word_size=word_size)
    def find_seq_per_barcode_identifier_using_blast(self, fasta_file: Union[str, Path], blast_db: Union[str, Path], query_results_path: Union[str, Path]):
        self.create_blast_database(fasta_file=fasta_file, blast_db=blast_db)
        self.extract_bc_identifier_and_reads_results_to_csv(blast_db=blast_db,
                                                            word_size=self.barcode_length+self.identifier_length,
                                                            blast_results=self.blast_bc_identifier_results,
                                                            query_results_path=query_results_path)

    def reads_chunk_to_fasta_format(self, reads, idx) -> str:
        # Write the sequences and their reverse complements to the output file in FASTA format
        output_fasta_file = self.reads_chunk_to_fasta_file + str(idx) + '.fasta'
        with open(output_fasta_file, "w") as output_handle:
            for read_id, read_seq in reads.items():
                # Write the original sequence
                original_seq_record = SeqRecord(read_seq.seq,
                                                id=read_id,
                                                description="original")
                SeqIO.write(original_seq_record, output_handle, "fasta")

                # Generate and write the reverse complement
                reverse_complement_seq = str(read_seq.seq.reverse_complement())
                reverse_complement_record = SeqRecord(Seq(reverse_complement_seq),
                                                      id=read_id,
                                                      description="reverse_complement")
                SeqIO.write(reverse_complement_record, output_handle, "fasta")

        return output_fasta_file

    def run_blastn_to_find_location_of_oligo_in_reads_chunks(self, reads_chunk_i: int,
                                                             output_fasta_file: Union[str, Path],
                                                             blast_db: Union[str, Path],
                                                             word_size: int) -> Dict:
        """qseqid: Query Seq-id
            sseqid: Subject Seq-id
            pident: Percentage of identical matches
            length: Alignment length
            mismatch: Number of mismatches
            gapopen: Number of gap openings
            qstart: Start of alignment in query
            qend: End of alignment in query
            sstart: Start of alignment in subject
            send: End of alignment in subject
            evalue: Expect value
            bitscore: Bit score
            btop: Blast traceback operations"""


        # blastn_cmd = (
        #     f'blastn -query {output_fasta_file} '
        #     f'-db {blast_db} '
        #     '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qseq sseq" '
        #     f'-word_size {max(word_size - 10, 20)} -gapopen 2 -gapextend 1 -reward 1 -penalty -2 '
        #     f'-min_raw_gapped_score {max(word_size - 10, 20)} -num_threads 6'
        # )

        blastn_cmd = (
            f'blastn -query {output_fasta_file} '
            f'-db {blast_db} '
            '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qseq sseq" '
            f'-gapopen 2 -gapextend 1 -reward 1 -penalty -2 '
            f'-min_raw_gapped_score {max(word_size - 10, 20)} -num_threads 6'
        )

        # Run the command and capture the output
        try:
            query_results = subprocess.check_output(blastn_cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(e.output.decode())

        blastn_output_idx_file = self.blast_database_path + 'blastn_output_' + str(reads_chunk_i) + '.txt'
        with open(blastn_output_idx_file, "wb") as binary_file:

            # Write bytes to file
            binary_file.write(query_results)
        query_results = {}
        with open(blastn_output_idx_file) as blastn_output_file:
            for line in blastn_output_file:
                fields = line.strip().split("\t")
                read_id = fields[0]
                uni_id = fields[1]
                # percentage_of_identical_matches = fields[2]
                alignment_length = fields[3]
                number_of_mismatches = fields[4]
                number_of_gap_openings = fields[5]
                start_of_alignment_in_read = fields[6]
                end_of_alignment_in_read = fields[7]
                # start_of_alignment_in_uni = fields[8]
                # end_of_alignment_in_uni = fields[9]
                # expect_value = fields[10]
                bit_score = fields[11]
                # blast_traceback_operations = fields[12]
                # aligned_part_of_query_sequence = fields[13]
                # aligned_part_of_subject_sequence = fields[14]
                # alignment_location = f"{uni_id}:{start_of_alignment_in_read}-{end_of_alignment_in_read}"
                if read_id not in query_results:
                    query_results[read_id] = {}
                    query_results[read_id][uni_id] = []
                if uni_id not in query_results[read_id]:
                    query_results[read_id][uni_id] = []
                query_results[read_id][uni_id].append((start_of_alignment_in_read, end_of_alignment_in_read, bit_score,
                                                       alignment_length, number_of_mismatches, number_of_gap_openings))

        return query_results

    def find_pos_in_read_list(self, bc_pos_in_read: Dict) -> List[int]:
        experiment_bc = self.convert_bc_pos_in_read_to_list(bc_pos_in_read=bc_pos_in_read)
        # merged_design_uni_pos_and_pos_in_read_list = self.merge_design_pos_and_pos_in_read_list(
        #     experiment_uni_pos=experiment_uni_pos_list)
        # return merged_design_uni_pos_and_pos_in_read_list
        return experiment_bc

    def convert_bc_pos_in_read_to_list(self, bc_pos_in_read: Dict) -> Tuple[int]:
        bc_pos_in_read_sorted = collections.OrderedDict(
            sorted(bc_pos_in_read.items(), key=lambda t: t[0], reverse=False))
        experiment_bc_pos = [0] * self.amount_of_bc
        experiment_bc = Tuple[int, int]
        best_bc_info_prev = -1
        for bc_name, bc_info in bc_pos_in_read_sorted.items():
            # experiment_uni_pos[]
            bc_idx = int(bc_name)
            best_bit_score = -1
            best_bc_info = -1
            for info in bc_info:
                if float(info[2]) > best_bit_score:
                    best_bit_score = float(info[2])
                    best_bc_info = info

            if best_bc_info_prev < int(best_bc_info[0]):
                best_bc_info_prev = int(best_bc_info[0]) - 1
                experiment_bc_pos[bc_idx - 1] = best_bc_info_prev
                experiment_bc = (bc_idx, best_bc_info_prev)

        # blast give the positions with +1, so we fix it by subtracting it by 1
        # return experiment_bc_pos
        return experiment_bc

    def identify_oligo(self, bc_pos: int, read: SeqIO.SeqRecord) -> SeqIO.SeqRecord:

        # self.identify_identifier(bc_pos_dict=bc_pos_dict,
        #                          read=read)
        payload = self.identify_payload(bc_pos=bc_pos,
                              read=read)
        return payload

    def identify_identifier(self, bc_pos_dict: Dict, read: SeqIO.SeqRecord):
        pass

    def identify_payload(self, bc_pos: int, read: SeqIO.SeqRecord) -> SeqIO.SeqRecord:
        return str(read[:bc_pos])

    def write_query_results_to_file(self, query_results: Dict, file_name: Union[str, Path]) -> None:
        with open(file_name, 'a') as f:
            for query_id, alignments in query_results.items():
                f.write(f"Query ID: {query_id}\n")
                for subject_id, matches in alignments.items():
                    f.write(f"  Subject ID: {subject_id}\n")
                    for match in matches:
                        f.write(f"    {match}\n")
                f.write("\n")  # Add a newline between each query result for readability
