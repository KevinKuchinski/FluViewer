import sys
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from collections import Counter
from math import ceil, log10


def main():
    version = '0.1.9'
    args = parse_args(sys.argv, version)
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer/\n')
    print('Job name:', args['-n'])
    print('Fwd reads:', args['-f'])
    print('Rev reads:', args['-r'])
    print('Reference sequences:', args['-d'])
    print()
    check_input_files(args['-f'], args['-r'], args['-d'])
    check_database(args['-d'])
    make_output_dir(args['-n'])
    normalize_depth(args['-n'], args['-f'], args['-r'], args['-N'], args['-g'])
    assemble_contigs(args['-n'], args['-g'])
    blast_results = blast_contigs(args['-d'], args['-n'], args['-g'],
                                  args['-T'], args['-i'], args['-l'])
    blast_results = filter_contig_blast_results(blast_results,
                                                args['-n'], args['-g'])
    make_scaffold_seqs(blast_results, args['-n'], args['-g'])
    blast_results = blast_scaffolds(args['-d'], args['-n'],
                                    args['-g'], args['-T'])
    blast_results = filter_scaffold_blast_results(blast_results)
    make_mapping_refs(blast_results, args['-d'], args['-n'])
    map_reads(args['-n'], args['-g'], args['-q'])
    call_variants(args['-n'], args['-q'], args['-L'], args['-g'])
    mask_ambig_low_cov(args['-n'], args['-D'], args['-v'], args['-V'],
                       args['-q'], args['-g'])
    make_consensus_seqs(args['-n'], args['-g'])
    write_report(args['-n'], args['-g'])
    make_plots(args['-n'], args['-g'])
    if args['-g']:
        garbage_collection(args['-n'])
    print('\nDone.\n')
    exit(0)


def parse_args(args, version):
    ''' Check that each argument has been set only once. '''
    arg_names = set(arg for arg in args if arg[0] == '-')
    multi_set_args = set()
    for arg in arg_names:
        if args.count(arg) > 1:
            multi_set_args.add(arg)
    if multi_set_args != set():
        multi_set_args = ', '.join(sorted(multi_set_args))
        print('\nERROR: The following arguments have been set more than once: '
              f'{multi_set_args}')
        print_usage(version)
        exit(1)
    ''' Create dict for runtime arguments and their provided values. '''
    arg_values = {}
    for arg_1, arg_2 in zip(args[1:-1], args[2:]):
        if arg_1[0] == '-':
            if arg_2[0] != '-':
                arg_values[arg_1] = arg_2
            else:
                arg_values[arg_1] = ''
    if args[-1][0] == '-':
        arg_values[args[-1]] = ''
    ''' Set defaults, mins, and maxs. '''
    required_args = {'-n', '-f', '-r', '-d'}
    arg_value_types = {'-n': str, '-f': str, '-r': str, '-d': str, '-i': float,
                       '-l': int, '-D': int, '-q': int, '-v': float,
                       '-V': float, '-N': int, '-L': int, '-T': int, '-g': bool}
    min_arg_values = {'-i': 0, '-l': 32, '-D': 1, '-q': 0,
                      '-v': 0, '-V': 0, '-N': 1, '-L':1, '-T': 0}
    max_arg_values = {'-i': 100, '-v': 1, '-V': 1}
    default_arg_values = {'-i': 90, '-l': 50, '-D': 20, '-q': 20, '-v': 0.75,
                          '-V': 0.25, '-N': 200, '-L': 200, '-g': True}
    ''' Set garbage flag to boolean value if not provided. '''
    if '-g' not in arg_values:
        arg_values['-g'] = True
    elif '-g' in arg_values and arg_values['-g'] == '':
        arg_values['-g'] = False
    elif '-g' in arg_values and arg_values['-g'] != '':
        print('\nERROR: Argument -g should be provided without a value.\n')
        print_usage(version)
        exit(1)
    ''' Check if all required arguments were provided. '''
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys():
            missing_args.add(required_arg)
    if missing_args != set():
        missing_args = ', '.join(sorted(missing_args))
        print(f'\nERROR: Values must be provided for the argument following '
              f'arguments: {missing_args}')
        print_usage(version)
        exit(1)
    ''' Check if unrecognized arguments were provided. '''
    recognized_args = set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args.add(provided_arg)
    if unrecognized_args != set():
        unrecognized_args = ', '.join(sorted(unrecognized_args))
        print(f'\nERROR: The following provided arguments are not recognized: '
              f'{unrecognized_args}')
        print_usage(version)
        exit(1)
    ''' Check if provided values are of the correct type. '''
    for arg, value in arg_values.items():
        bad_types = set()
        for arg, value in arg_values.items():
            try:
                arg_values[arg] = arg_value_types[arg](value)
            except ValueError:
                value_type = str(arg_value_types[arg])
                value_type = value_type.split('<')[1].split('>')[0]
                if arg_value_types[arg] != bool:
                    print(f'\nERROR: Value for argument {arg} must be of type '
                          f'{value_type}')
                else:
                    print(f'\nERROR: Argument {arg} must be provided without a '
                          f'value') 
                print_usage(version)
                exit(1)
    ''' Check if provided values are within the correct range. '''
    for arg, value in arg_values.items():
        bad_values = set()
        if arg in min_arg_values.keys() and value < min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be greater than '
                  f'{min_arg_values[arg]}')
            print_usage(version)
            exit(1)
        if arg in max_arg_values.keys() and value > max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be less than '
                  f'{max_arg_values[arg]}')
            print_usage(version)
            exit(1) 
    ''' Check if arguments were provided without values. '''
    empty_args = set()
    for arg, value in arg_values.items():
        if value == '':
            empty_args.add(arg)
    if empty_args != set():
        empty_args = ', '.join(sorted(empty_args))
        print(f'\nERROR: The following arguments were provided without values: '
              f'{empty_args}')
        print_usage(version)
        exit(1)
    ''' Check if variant masking threshold is lower than variant calling
    threshold. '''
    if '-v' in arg_values and '-V' in arg_values:
        if arg_values['-V'] >= arg_values['-v']:
            print(f'\nERROR: Variant allele fraction for masking ambiguous '
                  f'variants must be lower than variant allele fraction for '
                  f'calling variants.')
            print_usage(version)
            exit(1)
    ''' Assign default values to unspecified arguments. '''
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value
    ''' Return keyword args and their values. '''
    return arg_values


def print_usage(version):
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer\n')
    print('Usage: FluViewer -n <output_name> -f <path_to_fwd_reads> '
          '-r <path_to_rev_reads> -d <path_to_db_file> [ <optional_args> ]\n')
    print('Required arguments:')
    print(' -n : output name')
    print(' -f : path to FASTQ file containing forward reads')
    print(' -r : path to FASTQ file containing reverse reads')
    print(' -d : path to FASTA file containing FluViewer database')
    print('Optional arguments:')
    print(' -i : Minimum sequence identity between database reference'
          ' sequences and contigs (percentage, default = 90, min = 0,'
          ' max = 100)')
    print(' -l : Minimum length of alignment between database reference'
          ' sequences and contigs (int, default = 50, min = 32)')
    print(' -D : minimum read depth for base calling (int, default = 20, '
          ' min = 1)')
    print(' -q : Minimum PHRED score for mapping quality and base quality'
          ' during variant calling (int, default = 20, min = 0)')
    print(' -v : Variant allele fraction threshold for calling variants'
          ' (float, default = 0.75, min = 0, max = 1)')
    print(' -V : Variant allele fraction threshold for masking ambiguous variants'
          ' (float, default = 0.25, min = 0, max = 1')
    print(' -N : Target depth for pre-normalization of reads'
          ' (int, default = 200, min = 1)')
    print(' -L : Coverage depth limit for variant calling'
          ' (int, default = 200, min = 1)')
    print(' -T : Threads used for contig/scaffold alignments'
          ' (int, default = 1, min = 1)')
    print(' -g : Disable garbage collection and retain intermediate analysis'
          ' files (no value provided)')
    print()


def check_input_files(fwd_reads, rev_reads, db):
    ''' Check that all provided input files exist. '''
    for file in [fwd_reads, rev_reads, db]:
        if not os.path.isfile(file):
            print('\nERROR: Input file does not exist.\n')
            error_code = 1
            exit(error_code)


def check_database(db):
    ''' Checks the contents of the provided reference sequence database to
    ensure proper header formatting and unambiguous sequences. ''' 
    print('Checking reference sequence database...')
    total_entries = 0
    unique_headers = list()
    with open(db, 'r') as input_file:
        for line in input_file:
            line = line.strip()
            if line[0] == '>':
                total_entries += 1
                header = line
                unique_headers.append(header)
                if len(line.split('|')) != 4:
                    print(f'\nERROR: The header for the following database entry '
                          f'does not contain the expected number of |-delimited '
                          f'fields:\n{header}\n')
                    exit(1)
                else:
                    accession, name, segment, subtype = line.split('|')
                    if any([name.count('(') != 1, name.count(')') != 1,
                            name[-1] != ')', name.index('(') > name.index(')')]):
                        print(f'\nERROR: The strain name (strain subtype) for '
                              f'the following database entry is improperly'
                              f' formatted:\n{header}\n')
                        exit(1)
                    if segment not in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
                        print(f'\nERROR: The segment indicated for the following '
                              f'database entry is not recognized:\n{header}\n')
                        exit(1)
            else:
                if len(line.strip()) != sum(line.count(base) for base in 'ATGC'):
                        print(f'\nERROR: The sequence provided for the '
                              f'following database entry contains ambiguous or '
                              f'lower-case nucleotides:\n{header}\n')
                        exit(1)
    unique_headers = Counter(unique_headers)
    for header in unique_headers:
        if unique_headers[header] > 1:
            print(f'\nERROR: The following database entry does not have a unique '
                  f'header:\n{header}\n')
            exit(1)


def make_output_dir(out_name):
    ''' Create a directory for output after checking that it does not
    already exist. '''
    if os.path.exists(out_name) == False:
        os.mkdir(out_name)
        logs_path = os.path.join(out_name, 'logs')
        os.mkdir(logs_path)
    else:
        print('\nERROR: Output directory already exists! Aborting analysis.\n')
        error_code = 1
        exit(error_code)


def run(terminal_command, out_name, process_name, error_code, collect_garbage):
    ''' A generalized function for running subprocesses, logging their output, and
    trapping erroneous exit statuses. '''
    stdout_file = os.path.join(out_name, 'logs',
                               f'{process_name}_stdout.txt')
    stderr_file = os.path.join(out_name, 'logs',
                               f'{process_name}_stderr.txt')
    for file in [stdout_file, stderr_file]:
        if file != None:
            with open(file, 'w') as log_file:
                log_file.write('*' * 80 + '\n')
                log_file.write('Terminal command:' + '\n')
                log_file.write(terminal_command + '\n')
                log_file.write('*' * 80 + '\n')
    stdout_file = open(stdout_file, 'w')
    stderr_file = open(stderr_file, 'w')
    complete_process = subprocess.run(terminal_command, stdout=stdout_file,
                                      stderr=stderr_file, shell=True)
    stdout_file.close()
    stderr_file.close()
    return_code = complete_process.returncode
    if return_code != 0:
        print(f'\nERROR: Subprocess {process_name} failed (Exit status: '
              f'{return_code})\n')
        if collect_garbage:
            garbage_collection(out_name)
        exit(error_code)


def garbage_collection(out_name):
    ''' Clean up unneccessary intermediate files (many of which occupy
    substantial amounts of storage). '''
    if os.path.isdir(os.path.join(out_name, f'spades_output')):
        shutil.rmtree(os.path.join(out_name, f'spades_output'))
    files = []
    files += [f'R1.fq', f'R2.fq', 'contigs_blast.tsv', 'scaffolds.fa',
              'scaffolds_blast.tsv', 'alignment.sam', 'pileup.vcf',
              'variants.bcf', 'variants.bcf.csi', 'low_cov.tsv', 'ambig.tsv',
              'variants.tsv', 'masked.bed', 'reads_mapped.tsv',
              'depth_of_cov_samtools.tsv', 'depth_of_cov_freebayes.tsv']
    files += [f'{out_name}_mapping_refs.fa.' + suffix
              for suffix in ['amb', 'ann', 'bwt', 'fai', 'pac', 'sa']]
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    files += [f'{segment}_contigs.fa' for segment in segments]
    files += [f'{segment}_contigs.afa' for segment in segments]
    files += [f'{segment}_contigs.dnd' for segment in segments]
    for file in files:
        file = os.path.join(out_name, file)
        if os.path.isfile(file):
            os.remove(file)
    


def normalize_depth(out_name, fwd_reads_raw, rev_reads_raw, depth,
                    collect_garbage):
    ''' BBNorm is run on the input reads to downsample regions of deep coverage
    (using a k-mer frequency approach). This balances coverage, increases
    analysis speed, and limits the impacts of artefactual reads. '''
    print('Normalizing depth of coverage and subsampling reads...')
    fwd_reads = os.path.join(out_name, f'R1.fq')
    rev_reads = os.path.join(out_name, f'R2.fq')
    terminal_command = (f'bbnorm.sh in={fwd_reads_raw} in2={rev_reads_raw} '
                        f'out={fwd_reads} out2={rev_reads} target={depth}')
    process_name = 'bbnorm'
    error_code = 2
    run(terminal_command, out_name, process_name, error_code, collect_garbage)

    
def assemble_contigs(out_name, collect_garbage):
    ''' Normalized, downsampled reads are assembled de novo into contigs
    using SPAdes. '''
    print('Assembling reads into contigs...')
    spades_output = os.path.join(out_name, 'spades_output')
    fwd_reads = os.path.join(out_name, f'R1.fq')
    rev_reads = os.path.join(out_name, f'R2.fq')
    terminal_command = (f'spades.py --rnaviral --isolate -1 {fwd_reads} '
                        f'-2 {rev_reads} -o {spades_output}')
    process_name = 'spades'
    error_code = 3
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    if not os.path.isfile(os.path.join(spades_output, 'contigs.fasta')):
        print('\nERROR: No contigs assembled! Aborting analysis.\n')
        if collect_garbage:
            garbage_collection(out_name)
        error_code = 4
        exit(error_code)


def blast_contigs(db, out_name, collect_garbage, threads, identity, length):
    ''' Contigs are aligned to reference sequences using BLASTn. '''
    print('Aligning contigs to reference sequences...')
    if any([os.path.exists(db + '.' + suffix) == False
            for suffix in ['nhr', 'nin' , 'nsq']]):
        terminal_command = (f'makeblastdb -in {db} -dbtype nucl')
        process_name = 'makeblastdb_contigs'
        error_code = 5
        run(terminal_command, out_name, process_name, error_code,
            collect_garbage)
    blast_output = os.path.join(out_name, 'contigs_blast.tsv')
    spades_output = os.path.join(out_name, 'spades_output')
    contigs = os.path.join(spades_output, 'contigs.fasta')
    cols = 'qseqid sseqid pident length bitscore sstart send qseq sseq slen'
    terminal_command = (f'blastn -query {contigs} -db {db} '
                        f'-num_threads {threads} -outfmt "6 {cols}" '
                        f'> {blast_output}')
    process_name = 'blastn_contigs'
    error_code = 6
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    blast_results = pd.read_csv(blast_output, names=cols.split(' '), sep='\t')
    blast_results = blast_results[blast_results['pident']>=identity]
    blast_results = blast_results[blast_results['length']>=length]
    if len(blast_results) == 0:
        print(f'\nERROR: No contigs aligned to reference sequences! '
              f'Aborting analysis.\n')
        if collect_garbage:
            garbage_collection(out_name)
        error_code = 7
        exit(error_code)
    return blast_results


def filter_contig_blast_results(blast_results, out_name, collect_garbage):
    ''' Contigs alignments are filtered to discard spurious alignments (the
    length and sequence identity of each alignment must exceed certain
    thresholds). Afterwards, a single reference sequence is selected for each
    genome segment: the reference sequence with the most positions covered by
    contigs, with ties being broken by the highest sequence identity, then the
    longest reference sequence length, then the first alphabetically. Once a
    reference sequence has been chosen for each segment, only contig alignments
    to those reference sequences are retained. '''
    print('Filtering contig alignments...')
    ''' Annotate each ref seq with its segment and subtype. '''
    subject_annots = blast_results[['sseqid']].drop_duplicates()
    get_segment = lambda row: row['sseqid'].split('|')[2]
    subject_annots['segment'] = subject_annots.apply(get_segment, axis=1)
    get_subtype = lambda row: row['sseqid'].split('|')[3]
    subject_annots['subtype'] = subject_annots.apply(get_subtype, axis=1)
    blast_results = pd.merge(blast_results, subject_annots, on='sseqid')
    ''' Check for evidence of mixed infections. First, find best alignments(s)
    for each contig. Next, check if each segment is represented by only one
    subtype. '''
    cols = ['qseqid', 'bitscore']
    max_bitscores = blast_results[cols].groupby('qseqid').max().reset_index()
    best_results = pd.merge(blast_results, max_bitscores, on=cols)
    for segment in best_results['segment'].unique():
        segment_results = best_results[best_results['segment']==segment]
        segment_subtypes = segment_results['subtype'].unique()
        if len(segment_subtypes) > 1:
            segment_subtypes = ', '.join(segment_subtypes)
            print(f'\nERROR: Multiple subtypes detected for segment {segment} '
                  f'({segment_subtypes})! Aborting analysis.\n')
            if collect_garbage:
                garbage_collection(out_name)
            error_code = 8
            exit(error_code) 
    ''' Find ref seq(s) most covered by contigs. '''
    def count_cov_pos(data_frame):
        cov_positions = set()
        for index, row in data_frame.iterrows():
            start = min([row['sstart'], row['send']])
            end = max([row['sstart'], row['send']])
            cov_positions = cov_positions.union(set(range(start, end + 1)))
        return len(cov_positions)
    cols = ['sseqid', 'segment', 'subtype', 'sstart', 'send']
    group_cols = ['sseqid', 'segment', 'subtype']
    cov_pos = blast_results[cols].drop_duplicates()
    cov_pos = cov_pos.groupby(group_cols).apply(count_cov_pos).reset_index()
    cov_pos.columns = ['sseqid', 'segment', 'subtype', 'covered_positions']
    cols = ['segment', 'subtype', 'covered_positions']
    group_cols = ['segment', 'subtype']
    max_cov_pos = cov_pos[cols].drop_duplicates()
    max_cov_pos = max_cov_pos.groupby(group_cols).max().reset_index()
    merge_cols = ['segment', 'subtype', 'covered_positions']
    max_cov_pos = pd.merge(cov_pos, max_cov_pos, on=merge_cols)
    max_cov_pos = max_cov_pos[['sseqid', 'covered_positions']]
    blast_results = pd.merge(blast_results, max_cov_pos, on='sseqid')
    ''' Find remaining ref seq(s) with most identical positions. '''
    def count_id_pos(data_frame):
        identical_positions = set()
        for index, row in data_frame.iterrows():
            start = min([row['sstart'], row['send']])
            increment = 1 if row['sstart'] <= row['send'] else -1
            subject_position = start
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if sbase in 'ATGC' and qbase == sbase:
                    identical_positions.add(subject_position)
                if sbase != '-':
                    subject_position += increment
        return len(identical_positions)
    cols = ['sseqid', 'segment', 'subtype', 'sstart', 'send', 'qseq', 'sseq']
    group_cols = ['sseqid', 'segment', 'subtype']
    ident_pos = blast_results[cols].drop_duplicates()
    ident_pos = ident_pos.groupby(group_cols).apply(count_id_pos).reset_index()
    ident_pos.columns = ['sseqid', 'segment', 'subtype', 'identical_positions']
    cols = ['segment', 'subtype', 'identical_positions']
    group_cols = ['segment', 'subtype']
    max_ident_pos = ident_pos[cols].drop_duplicates()
    max_ident_pos = max_ident_pos.groupby(group_cols).max().reset_index()
    merge_cols = ['segment', 'subtype', 'identical_positions']
    max_ident_pos = pd.merge(ident_pos, max_ident_pos, on=merge_cols)
    cols = ['sseqid', 'identical_positions']
    max_ident_pos = max_ident_pos[cols].drop_duplicates()
    blast_results = pd.merge(blast_results, max_ident_pos, on='sseqid')
    ''' Take longest remaining ref seq for each segment/subtype. '''
    cols = ['segment', 'subtype', 'slen']
    group_cols = ['segment', 'subtype']
    longest_sseq = blast_results[cols].drop_duplicates()
    longest_sseq = longest_sseq.groupby(group_cols).min().reset_index()
    blast_results = pd.merge(blast_results, longest_sseq, on=cols)
    ''' Take first alphabetical remaining ref seq for each segment/subtype. '''
    cols = ['segment', 'subtype', 'sseqid']
    group_cols = ['segment', 'subtype']
    first_alpha_sseq = blast_results[cols].drop_duplicates()
    first_alpha_sseq = first_alpha_sseq.groupby(group_cols).min().reset_index()
    blast_results = pd.merge(blast_results, first_alpha_sseq, on=cols)
    return blast_results


def make_scaffold_seqs(blast_results, out_name, collect_garbage):
    ''' A scaffold sequence is created for each genome segment by joining and
    collapsing all the contigs describing that segment. First, the contigs are
    aligned to reference sequences, and a single reference sequence is chosen
    for each reference sequence (this occurs in the filter_contigs_blast_results
    func). Next, unaligned leading and trailing sequences are trimmed from the
    contigs. Next, leading and trailing Ns are added to the contig so that it is
    properly positioned within the segment (based on the subject-start and
    subject-end coordinates of its alignment to the selected reference sequence).
    Next, clustalW is used to generate a multiple sequence alignment of the
    trimmed, positioned contigs. This multiple sequence alignment is used to
    generate a consensus sequence of the regions of the segment covered by
    contigs. '''
    print('Creating scaffolds...')
    ''' Make sure contigs are all in the forward orientation. '''
    rev_comp_bases = {'A': 'T',
                      'T': 'A',
                      'G': 'C',
                      'C': 'G',
                      'N': 'N',
                      '-': '-',
                      'W': 'W',
                      'S': 'S',
                      'M': 'K',
                      'K': 'M',
                      'R': 'Y',
                      'Y': 'R',
                      'B': 'V',
                      'D': 'H',
                      'H': 'D',
                      'V': 'B'}
    rev_comp_seq = lambda seq: ''.join(rev_comp_bases[base]
                                       for base in seq[::-1])
    get_start = lambda row: min([row['sstart'], row['send']])
    blast_results['start'] = blast_results.apply(get_start, axis=1)
    get_end = lambda row: max([row['sstart'], row['send']])
    blast_results['end'] = blast_results.apply(get_end, axis=1)
    def flip_qseq(row):
        if row['sstart'] > row['send']:
            return rev_comp_seq(row['qseq'])
        else:
            return row['qseq']
    blast_results['qseq'] = blast_results.apply(flip_qseq, axis=1)
    def flip_sseq(row):
        if row['sstart'] > row['send']:
            return rev_comp_seq(row['sseq'])
        else:
            return row['sseq']
    blast_results['sseq'] = blast_results.apply(flip_sseq, axis=1)
    ''' Trim contigs based on their alignments to reference sequences. Also
    add leading and trailing Ns to contig so that it is properly positioned
    within the genome segment. '''
    segments = blast_results['segment'].unique()
    contig_counter = {segment: 0 for segment in segments}
    scaffold_seqs = {}
    for segment in segments:
        contigs = os.path.join(out_name, f'{segment}_contigs.fa')
        with open(contigs, 'w') as output_file:
            contig_results = blast_results[blast_results['segment']==segment]
            for index, row in contig_results.iterrows():
                header = f'>{segment}_contig_{contig_counter[segment]}\n'
                output_file.write(header)
                seq = 'N' * (row['start'] - 1)
                seq += row['qseq'].replace('-', '')
                seq += ('N' * (row['slen'] - row['end']))
                output_file.write(seq + '\n')
                contig_counter[segment] += 1
        ''' Generate multiple sequence alignments of trimmed/positioned
        contigs. '''
        aligned_contigs = os.path.join(out_name, f'{segment}_contigs.afa')
        if contig_counter[segment] > 1:
            terminal_command = (f'clustalw -INFILE={contigs} '
                                f'-OUTFILE={aligned_contigs} -OUTPUT=FASTA')
            process_name = f'clustalw_{segment}'
            error_code = 9
            run(terminal_command, out_name, process_name, error_code,
                collect_garbage)
        else:
            shutil.copyfile(contigs, aligned_contigs)
        ''' Replace leading and trailing Ns with dots so that they are ignored
        when determining consensus bases. '''
        seqs = {}
        with open(aligned_contigs, 'r') as input_file:
            for line in input_file:
                if line[0] == '>':
                    header = line.strip()
                    seqs[header] = ''
                else:
                    seqs[header] += line.strip()
        clean_seqs = []
        for seq in seqs.values():
            head_len = len(seq) - len(seq.lstrip('N-'))
            tail_len = len(seq) - len(seq.rstrip('N-'))
            seq = seq.strip('N-')
            seq = ('.' * head_len) + seq
            seq += ('.' * tail_len)
            clean_seqs.append(seq)
        ''' Check that all seqs in the multiple seq alignment are the
        same length. '''
        alignment_lengths = set(len(seq) for seq in clean_seqs)
        if len(alignment_lengths) > 1:
            print(f'\nERROR: Multiple sequence alignment for {segment} '
                  f'generated unequal alignment lengths! Aborting analysis.\n')
            if collect_garbage:
                garbage_collection(out_name)
            error_code = 10
            exit(error_code)
        ''' Create consensus sequence of multiple seq alignments, i.e. the
        scaffolds. '''
        alignment_length = list(alignment_lengths)[0]
        scaffold_seq = ''
        for i in range(alignment_length):
            bases = Counter(seq[i] for seq in clean_seqs if seq[i] not in '.')
            if bases.most_common(1) == []:
                scaffold_seq += 'N'
            elif bases.most_common(1)[0][1] / len(bases) > 0.5:
                scaffold_seq += bases.most_common(1)[0][0]
            else:
                scaffold_seq += 'N'
        scaffold_seqs[segment] = scaffold_seq
    ''' Write out scaffolds. '''
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    segments = [segment for segment in segments if segment in scaffold_seqs]
    scaffold_seqs = {segment: scaffold_seqs[segment] for segment in segments}
    scaffolds = os.path.join(out_name, 'scaffolds.fa')
    with open(scaffolds, 'w') as output_file:
        for segment, seq in scaffold_seqs.items():
            header = f'>{out_name}|{segment}_scaffold'
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def blast_scaffolds(db, out_name, collect_garbage, threads):
    ''' Scaffold sequences are aligned to reference sequences using BLASTn. '''
    print('Aligning scaffolds to reference sequences...')
    db_seqs = sum(line[0] == '>' for line in open(db, 'r').readlines())
    if any([os.path.exists(db + '.' + suffix) == False
            for suffix in ['nhr', 'nin' , 'nsq']]):
        terminal_command = (f'makeblastdb -in {db} -dbtype nucl')
        process_name = 'makeblastdb_scaffolds'
        error_code = 11
        run(terminal_command, out_name, process_name, error_code,
            collect_garbage)
    scaffolds = os.path.join(out_name, 'scaffolds.fa')
    blast_output = os.path.join(out_name, 'scaffolds_blast.tsv')
    cols = 'qseqid sseqid bitscore sstart send qseq'
    terminal_command = (f'blastn -query {scaffolds} -db {db} '
                        f'-num_threads {threads} -max_target_seqs '
                        f'{db_seqs} -outfmt "6 {cols}" > {blast_output}')
    process_name = 'blastn_scaffolds'
    error_code = 12
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    blast_results = pd.read_csv(blast_output, names=cols.split(' '), sep='\t')
    if len(blast_results) == 0:
        print(f'\nERROR: No scaffolds aligned to reference sequences! '
              f'Aborting analysis.\n')
        if collect_garbage:
            garbage_collection(out_name)
        error_code = 13
        exit(error_code)
    return blast_results


def filter_scaffold_blast_results(blast_results):
    ''' A single reference sequence is chosen for each segment scaffold. First,
    the bitscores of all alignments between a scaffold and a reference sequence
    are summed. The bitscore sum is calculated for each scaffold-reference
    sequence pairing. The reference sequence giving the highest bitscore sum is
    selected for each scaffold, with ties being broken by using the first
    alphabetically. Once the best-matching reference sequence has been selected for
    each segment scaffold, all other alignments are discarded. '''
    print('Filtering scaffold alignments...')
    ''' Annotate scaffold seqs with segment. '''
    query_annots = blast_results[['qseqid']].drop_duplicates()
    get_segment = lambda row: row['qseqid'].split('|')[1].split('_')[0]
    query_annots['segment'] = query_annots.apply(get_segment, axis=1)
    blast_results = pd.merge(blast_results, query_annots, on='qseqid')
    ''' Find best-matching reference sequence for each segment. '''
    cols = ['segment', 'sseqid', 'bitscore']
    group_cols = ['segment', 'sseqid']
    combo_scores = blast_results[cols].groupby(group_cols).sum().reset_index()
    cols = ['segment', 'bitscore']
    group_cols = ['segment']
    max_scores = combo_scores[cols].groupby(group_cols).max().reset_index()
    merge_cols = ['segment', 'bitscore']
    max_scores = pd.merge(max_scores, combo_scores, on=merge_cols)
    cols = ['segment', 'sseqid']
    group_cols = ['segment']
    first_alpha = max_scores[cols].groupby(group_cols).min().reset_index()
    merge_cols = ['segment', 'sseqid']
    blast_results = pd.merge(blast_results, first_alpha, on=merge_cols)
    get_start = lambda row: min([row['sstart'], row['send']])
    blast_results['start'] = blast_results.apply(get_start, axis=1)
    get_end = lambda row: max([row['sstart'], row['send']])
    blast_results['end'] = blast_results.apply(get_end, axis=1)
    return blast_results


def make_mapping_refs(blast_results, db, out_name):
    ''' Mapping references are created for each genome segment. These consist of
    the scaffold for that segment, with all Ns in the scaffold filled-in using
    the corresponding positions from that scaffold's best-matching reference
    sequence. '''
    print('Creating mapping references...')
    ''' Create dict with best-matching ref seq for each segment. '''
    sseqids = blast_results['sseqid'].unique()
    best_ref_seqs = {seq_name: '' for seq_name in sseqids}
    with open(db, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip()[1:]
            elif header in best_ref_seqs:
                best_ref_seqs[header] += line.strip()
    ''' Create mapping ref for each segment. '''
    def make_map_ref(data_frame):
        data_frame = data_frame.sort_values(by='start')
        sseq = best_ref_seqs[data_frame['sseqid'].min()]
        last_position = 0
        seq = ''
        for index, row in data_frame.iterrows():
            seq += sseq[last_position:row['start'] - 1]
            seq += row['qseq'].upper()
            last_position = row['end']
        seq += sseq[last_position:].lower()
        seq = seq.replace('-', '')
        return seq
    cols = ['sseqid', 'start', 'end', 'qseq']
    group_cols = ['sseqid']
    blast_results = blast_results[cols]
    blast_results = blast_results.groupby(group_cols).apply(make_map_ref)
    blast_results = blast_results.reset_index()
    blast_results.columns = ['sseqid', 'mapping_seq']
    ''' Annotate segment and subtype. '''
    get_segment = lambda row: row['sseqid'].split('|')[2].split('_')[0]
    blast_results['segment'] = blast_results.apply(get_segment, axis=1)
    get_subtype = lambda row: row['sseqid'].split('|')[3].split('_')[0]
    blast_results['subtype'] = blast_results.apply(get_subtype, axis=1)
    ''' Write mapping refs to FASTA. '''
    segment_order = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    get_segment_order = lambda row: segment_order.index(row['segment'])
    blast_results['sort'] = blast_results.apply(get_segment_order, axis=1)
    blast_results = blast_results.sort_values(by='sort')
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs, 'w') as output_file:
        for index, row in blast_results.iterrows():
            accession, ref_name, segment, subtype = row['sseqid'].split('|')[:4]
            accession = accession.lstrip('>')
            ref_name = ref_name.replace('(', '|').replace(')', '')
            header = f'>{out_name}|{segment}|{subtype}|{accession}|{ref_name}'
            seq = row['mapping_seq']
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def map_reads(out_name, collect_garbage, min_qual):
    ''' Normalized, downsampled reads (normalize_depth func) are mapped to the
    mapping references (make_mapping_refs func) using BWA mem. The alignment
    is filtered to retain only paired reads, then sorted and indexed. '''
    print('Mapping reads to mapping references...')
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    terminal_command = (f'bwa index {mapping_refs}')
    process_name = 'bwa_index'
    error_code = 14
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    fwd_reads = os.path.join(out_name, 'R1.fq')
    rev_reads = os.path.join(out_name, 'R2.fq')
    alignment = os.path.join(out_name, 'alignment.sam')
    terminal_command = (f'bwa mem {mapping_refs} {fwd_reads} {rev_reads} '
                        f'> {alignment}')
    process_name = 'bwa_mem'
    error_code = 15
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    filtered_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    terminal_command = (f'samtools view -f 1 -F 2828 -q {min_qual} '
                        f'-h {alignment} | samtools sort -o {filtered_alignment}')
    process_name = 'samtools_view'
    error_code = 16
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    terminal_command = (f'samtools index {filtered_alignment}')
    process_name = 'samtools_index'
    error_code = 17
    run(terminal_command, out_name, process_name, error_code, collect_garbage)


def call_variants(out_name, min_qual, max_depth, collect_garbage):
    ''' FreeBayes is used to create a pileup and call variants from the
    BAM file output (map_reads func). '''
    print('Calling variants...')
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    filtered_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    pileup = os.path.join(out_name, 'pileup.vcf')
    terminal_command = (f'freebayes -f {mapping_refs} {filtered_alignment} -p 1 '
                        f'--limit-coverage {max_depth} '
                        f'--min-mapping-quality {min_qual} '
                        f'--min-base-quality {min_qual} --pooled-continuous '
                        f'--report-monomorphic --haplotype-length 0 '
                        f'--min-alternate-count 1 --min-alternate-fraction 0 '
                        f'> {pileup}')
    process_name = 'freebayes'
    error_code = 18
    run(terminal_command, out_name, process_name, error_code, collect_garbage)


def mask_ambig_low_cov(out_name, min_depth, vaf_call, vaf_ambig,
                       min_qual, collect_garbage):
    ''' The FreeBayes VCF output is parsed, analyzing total read depth and
    variant read depth at each position. This allows positions to be masked
    for low coverage based on the read depth considered by FreeBayes (which
    could be lower than the read depth in the BAM depending on how FreeBayes
    applies it mapping quality and base quality filters). This also allows
    positions to be masked as ambiguous when the number of reads differing
    from the reference exceeds a threshold, but is not sufficient enough to
    confidently call as a specific variant. '''
    print('Masking ambiguous and low coverage positions...')
    ''' Open input/output files and initialize dicts. '''
    pileup = open(os.path.join(out_name, 'pileup.vcf'), 'r')
    variants = open(os.path.join(out_name, f'{out_name}_variants.vcf'), 'w')
    depth_of_cov = os.path.join(out_name, f'depth_of_cov_freebayes.tsv')
    depth_of_cov = open(depth_of_cov, 'w')
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    low_cov_pos = {segment: set() for segment in segments}
    ambig_pos = {segment: set() for segment in segments}
    variant_pos = {segment: set() for segment in segments}
    segment_name, segment_length = dict(), {None: 0}
    segment_length[None] = 0
    ''' Parse header '''
    line = pileup.readline()
    while line != '' and line[0] == '#':
        variants.write(line)
        if line[:10] == '##contig=<':
            name = line.strip().split('<ID=')[1].split(',length=')[0]
            segment = name.split('|')[1]
            length = int(line.strip().split(',length=')[1].split('>')[0])
            segment_name[segment] = name
            segment_length[segment] = length
        line = pileup.readline()
    ''' Parse body '''
    last_segment = None
    last_position = 0
    while line != '':
        fields = line.strip().split('\t')
        fields = (fields[0], fields[1], fields[3], fields[4], fields[5],
                  fields[8], fields[9])
        name, position, ref, alt, qual, keys, values = fields
        segment = name.split('|')[1]
        if segment != last_segment:
            if last_position < segment_length[last_segment]:
                for p in range(last_position + 1,
                               segment_length[last_segment] + 1):
                    low_cov_pos[last_segment].add(p)
                    depth_of_cov_line = [name, str(p), '0']
                    depth_of_cov_line = '\t'.join(depth_of_cov_line)
                    depth_of_cov.write(depth_of_cov_line + '\n')
            last_position = 0
        last_segment = segment
        position = int(position)
        if position != last_position + 1:
            for p in range(last_position + 1, position):
                low_cov_pos[segment].add(p)
                depth_of_cov_line = [name, str(p), '0']
                depth_of_cov_line = '\t'.join(depth_of_cov_line)
                depth_of_cov.write(depth_of_cov_line + '\n')
        qual = float(qual)
        info = {k: v for k, v in zip(keys.split(':'), values.split(':'))}
        if 'DP' in info and info['DP'].isnumeric():
            total_depth = int(info['DP'])
        else:
            total_depth = 0
        if 'AO' in info:
            alt_depths = tuple(int(i) if i.isnumeric() else 0
                               for i in info['AO'].split(','))
        else:
            alt_depths = (0, )
        max_alt_depth = max(alt_depths)
        total_alt_depth = sum(alt_depths)
        max_vaf = max_alt_depth / total_depth if total_depth > 0 else 0
        total_vaf = total_alt_depth / total_depth if total_depth > 0 else 0
        if all([qual >= min_qual, max_vaf >= vaf_call,
                total_depth >= min_depth]):
            variants.write(line)
        position -= 1
        for p in ref:
            position += 1
            depth_of_cov_line = [name, str(position), str(total_depth)]
            depth_of_cov_line = '\t'.join(depth_of_cov_line)
            depth_of_cov.write(depth_of_cov_line + '\n')
            if total_depth < min_depth:
                low_cov_pos[segment].add(position)
            elif total_vaf >= vaf_ambig and max_vaf < vaf_call:
                ambig_pos[segment].add(position)
        last_position = position
        line = pileup.readline()
    if last_position < segment_length[last_segment]:
        for p in range(last_position + 1, segment_length[last_segment] + 1):
            low_cov_pos[last_segment].add(p)
            depth_of_cov_line = [name, str(p), '0']
            depth_of_cov_line = '\t'.join(depth_of_cov_line)
            depth_of_cov.write(depth_of_cov_line + '\n')
    ''' Close input/output files '''
    pileup.close()
    variants.close()
    depth_of_cov.close()
    ''' Convert sets of low cov positions into tuples representing zero-indexed
    spans of masked positions (start, end).'''
    masked_pos = dict()
    for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
        masked_pos[segment] = low_cov_pos[segment].union(ambig_pos[segment])
        masked_pos[segment] = sorted(masked_pos[segment])
    spans = {segment: set() for segment in segments}
    segments = [segment for segment in segments
                if masked_pos[segment] != list()]
    for segment in segments:
        span_start = masked_pos[segment][0]
        for pos_A, pos_B in zip(masked_pos[segment][:-1],
                                masked_pos[segment][1:]):
            if pos_B != pos_A + 1:
                span_end = pos_A
                spans[segment].add((span_start - 1, span_end - 1))
                span_start = pos_B
        span_end = masked_pos[segment][-1]
        spans[segment].add((span_start - 1, span_end - 1))
    spans = {segment: sorted(spans[segment]) for segment in segments}
    ''' Write spans of low cov positions to TSV file for depth of coverage
    plots. '''
    with open(os.path.join(out_name, 'low_cov.tsv'), 'w') as output_file:
        for segment, segment_spans in spans.items():
            for (start, end) in segment_spans:
                line = [segment_name[segment], start, end]
                line = '\t'.join(str(i) for i in line)
                output_file.write(line + '\n')
    ''' Write ambiguous positions to TSV file. '''
    with open(os.path.join(out_name, 'ambig.tsv'), 'w') as output_file:
        for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
            for position in ambig_pos[segment]:
                line = [segment_name[segment], position]
                line = '\t'.join(str(i) for i in line)
                output_file.write(line + '\n')
    ''' Write variant positions to TSV file. '''
    with open(os.path.join(out_name, 'variants.tsv'), 'w') as output_file:
        for segment in 'PB2 PB1 PA HA NP NA M NS'.split(' '):
            for position in variant_pos[segment]:
                line = [segment_name[segment], position]
                line = '\t'.join(str(i) for i in line)
                output_file.write(line + '\n')
    ''' Write spans of masked positions to BED file in BedGraph format. '''
    with open(os.path.join(out_name, 'masked.bed'), 'w') as output_file:
        for segment, segment_spans in spans.items():
            for (start, end) in segment_spans:
                line = [segment_name[segment], start, end + 1, 0]
                line = '\t'.join(str(i) for i in line)
                output_file.write(line + '\n')


def make_consensus_seqs(out_name, collect_garbage):
    ''' High quality variants and masked positions (mask_ambig_low_cov func) are
    applied to the mapping references (make_mapping_refs) to generate the final
    consensus sequences for each segment. '''
    print('Generating consensus sequences...')
    ''' Zip and index VCF. '''
    variants = os.path.join(out_name, f'{out_name}_variants.vcf')
    zipped_variants = os.path.join(out_name, 'variants.bcf')
    terminal_command = (f'bcftools view {variants} -Ob -o {zipped_variants}')
    process_name = 'bcftools_view'
    error_code = 19
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    terminal_command = (f'bcftools index {zipped_variants}')
    process_name = 'bcftools_index'
    error_code = 20
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    ''' Apply variants to mapping refs. '''
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    masked = os.path.join(out_name, 'masked.bed')
    consensus_seqs = os.path.join(out_name, f'{out_name}_consensus_seqs.fa')
    terminal_command = (f'cat {mapping_refs} | bcftools consensus -m {masked} '
                        f'{zipped_variants} > {consensus_seqs}')
    process_name = 'bcftools_consensus'
    error_code = 21
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    ''' Reformat FASTA headers and remove whitespace. '''
    clean_seqs = {}
    with open(consensus_seqs, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip()
                clean_seqs[header] = ''
            else:
                clean_seqs[header] += line.strip().upper()
    with open(consensus_seqs, 'w') as output_file:
        for header, seq in clean_seqs.items():
            header = '|'.join(header.split('|')[:3])
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def write_report(out_name, collect_garbage):
    ''' Generate a report for each segment. '''
    print('Writing report...')
    ''' Count reads mapped to each segment and add to report. '''
    filtered_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    reads_mapped = os.path.join(out_name, 'reads_mapped.tsv')
    terminal_command = (f'samtools idxstats {filtered_alignment} > '
                        f'{reads_mapped}')
    process_name = 'samtools_idxstats'
    error_code = 22
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    cols = 'seq_name seq_length reads_mapped reads_unmapped'.split(' ')
    reads_mapped = pd.read_csv(reads_mapped, sep='\t', names=cols)
    reads_mapped = reads_mapped.replace('*', np.nan).dropna()
    get_seq_name = lambda row: '|'.join(row['seq_name'].split('|')[:3])
    reads_mapped['seq_name'] = reads_mapped.apply(get_seq_name, axis=1) 
    cols = ['seq_name', 'reads_mapped', 'seq_length']
    report = reads_mapped[cols].drop_duplicates()
    ''' Annotate segment and subtype. '''
    get_segment = lambda row: row['seq_name'].split('|')[1]
    report['segment'] = report.apply(get_segment, axis=1)
    get_subtype = lambda row: row['seq_name'].split('|')[2]
    report['subtype'] = report.apply(get_subtype, axis=1) 
    ''' Add scaffold completeness to report. '''
    seqs = {}
    scaffolds = os.path.join(out_name, 'scaffolds.fa')
    with open(scaffolds, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                seq_name = line[1:].strip()
                seqs[seq_name] = ''
            else:
                seqs[seq_name] += line.strip()
    completeness = {}
    for seq_name, seq in seqs.items():
        segment = seq_name.split('|')[1].split('_')[0]
        perc = sum(seq.count(base) for base in 'ATGC') * 100 / len(seq)
        perc = round(perc, 2)
        completeness[segment] = perc
    report['scaffold_completeness'] = report['segment'].map(completeness)
    ''' Add consensus completeness to report. '''                                                                                                                                                                     
    seqs = {}
    consensus_seqs = os.path.join(out_name, f'{out_name}_consensus_seqs.fa')
    with open(consensus_seqs, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                seq_name = line[1:].strip()
                seqs[seq_name] = ''
            else:
                seqs[seq_name] += line.strip()
    completeness = {}
    for seq_name, seq in seqs.items():
        segment = seq_name.split('|')[1]
        perc = sum(seq.count(base) for base in 'ATGC') * 100 / len(seq)
        perc = round(perc, 2)
        completeness[segment] = perc
    report['consensus_completeness'] = report['segment'].map(completeness)
    ''' Add best ref seq to report. '''                                                                                                                                                                               
    ref_seqs_used = {}
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                line = line[1:].strip().split('|')
                seq_name, segment, subtype = line[:3]
                accession, ref_name, ref_subtype = line[3:]
                seq_name = f'{seq_name}|{segment}|{subtype}'
                ref_seqs_used[seq_name] = (f'{accession}|{ref_name}'
                                           f'({ref_subtype})')
    report['ref_seq_used'] = report['seq_name'].map(ref_seqs_used)
    ''' Write report to TSV file. '''                                                                                                                                
    segment_order ='PB2 PB1 PA HA NP NA M NS'.split(' ')
    get_sort_value = lambda row: segment_order.index(row['segment'])
    report['sort'] = report.apply(get_sort_value, axis=1)                                                                                                                              
    report = report.sort_values(by='sort')                                                                                                                                                                              
    cols = ['seq_name', 'segment', 'subtype', 'reads_mapped', 'seq_length',
            'scaffold_completeness', 'consensus_completeness', 'ref_seq_used']
    report = report[cols]
    report.to_csv(os.path.join(out_name, f'{out_name}_report.tsv'),
                  index=False, sep='\t')


def make_plots(out_name, collect_garbage):
    ''' Generate depth of coverage plots for each segment. '''
    print('Making depth of coverage plots...')
    ''' Get depth of coverage using samtools. '''
    filtered_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    depth_of_cov = os.path.join(out_name, 'depth_of_cov_samtools.tsv')
    terminal_command = (f'samtools depth {filtered_alignment} > {depth_of_cov}')
    process_name = 'samtools_depth'
    error_code = 23
    run(terminal_command, out_name, process_name, error_code, collect_garbage)
    cols = ['seq_name', 'position', 'depth_samtools']
    samtools_data = pd.read_csv(depth_of_cov, sep='\t', names=cols)
    ''' Get depth of coverage from FreeBayes. '''
    depth_of_cov = os.path.join(out_name, 'depth_of_cov_freebayes.tsv')
    cols = ['seq_name', 'position', 'depth_freebayes']
    freebayes_data = pd.read_csv(depth_of_cov, sep='\t', names=cols)
    ''' Merge samtools and freebayes data. '''
    data = pd.merge(samtools_data, freebayes_data, on=['seq_name', 'position'],
                    how='left')
    ''' Annotate with segment. '''
    get_segment = lambda row: row['seq_name'].split('|')[1]
    data['segment'] = data.apply(get_segment, axis=1)
    ''' Make plots. '''
    sb.set_style('whitegrid')
    segments = data['segment'].unique()
    fig_size = (8.5, 2 * len(segments))
    fig, axs = plt.subplots(len(segments), 1, sharex=True, figsize=fig_size)
    max_position = data['position'].max()
    x_ticks = [100 * i for i in range(1, ceil(max_position / 100))]
    max_depth = max([data['depth_samtools'].max(), data['depth_freebayes'].max()])
    y_max = 10 ** ceil(log10(max_depth))
    y_ticks = [10 ** i for i in range(ceil(log10(max_depth)))]
    for segment, ax in zip(segments, axs):
        segment_data = data[data['segment']==segment]
        sb.lineplot(x='position', y='depth_samtools', ax=ax, data=segment_data,
                    color='black')
        sb.lineplot(x='position', y='depth_freebayes', ax=ax, data=segment_data,
                    color='blue')
        ax.set_xlim(1, max_position)
        ax.set_xlabel('Position')
        ax.set_xticks(x_ticks)
        ax.set_ylim(1, y_max)
        ax.set_ylabel('Read depth')
        ax.set_yscale('log')
        ax.set_yticks(y_ticks)
        ax.set_title(segment)
        ax.axvspan(segment_data['position'].max(), max_position, color='grey')
        with open(os.path.join(out_name, 'low_cov.tsv'), 'r') as input_file:
            for line in input_file:
                seq_name, start, stop = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvspan(int(start), int(stop), color='red', alpha=0.1)
        with open(os.path.join(out_name, 'ambig.tsv'), 'r') as input_file:
            for line in input_file:
                seq_name, position = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvline(position, color='orange')
        with open(os.path.join(out_name, 'variants.tsv'), 'r') as input_file:
            for line in input_file:
                seq_name, position = line.strip().split('\t')
                if seq_name.split('|')[1] == segment:
                    ax.axvline(position, color='blue')        
    plt.tight_layout()
    plots = os.path.join(out_name, f'{out_name}_depth_of_cov.png')
    plt.savefig(plots, dpi=400)


if __name__ == '__main__':
    main()
