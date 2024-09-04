import sys as sys
import os as os
import shutil as shutil
import subprocess as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb


from collections import Counter
from math import log10


def main():
    version = '0.2.1'
    args = parse_args(sys.argv, version)
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer/\n')
    print('Job name:', args['-n'])
    print('Fwd reads:', args['-f'])
    print('Rev reads:', args['-r'])
    print('Reference sequences:', args['-d'])
    print()
    check_input_files(args['-f'], args['-r'], args['-d'])
    make_output_dir(args['-n'])
    length_ranges = check_database(args['-d'])
    assemble_contigs(args['-n'], args['-f'], args['-r'])
    blast_results = align_contigs_to_ref_seqs(args['-n'], args['-d'],
                                              args['-T'], args['-i'],
                                              args['-l'])
    blast_results = filter_contig_alignments(blast_results)
    blast_results = mixed_infection_check(blast_results, args['-m'])
    make_scaffolds(args['-n'], blast_results)
    blast_results = align_scaffolds_to_ref_seqs(args['-n'], args['-d'],
                                                args['-T'])
    blast_results = filter_scaffold_alignments(blast_results)
    make_mapping_refs(args['-n'], blast_results, args['-d'])
    map_reads_to_mapping_refs(args['-n'], args['-f'], args['-r'], args['-q'])
    make_consensus_seqs(args['-n'], args['-L'], args['-q'], args['-D'],
                        args['-V'], args['-v'], length_ranges, args['-t'])
    write_report(args['-n'], args['-D'], args['-V'], args['-v'])
    make_plot(args['-n'], args['-D'], args['-V'], args['-v'])
    if args['-g']:
        collect_garbage(args['-n'])
    print('\nDone.\n')


def print_usage(version: str):
    ''' Display usage message and condensed manual for FluViewer. '''
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
    print(' -m : allow analysis of mixed infections'
          ' (no value provided)')
    print(' -D : minimum read depth for base calling (int, default = 20, '
          ' min = 1)')
    print(' -q : Minimum PHRED score for mapping quality and base quality'
          ' during variant calling (int, default = 20, min = 0)')
    print(' -v : Variant allele fraction threshold for calling variants'
          ' (float, default = 0.75, min = 0, max = 1)')
    print(' -V : Variant allele fraction threshold for masking ambiguous variants'
          ' (float, default = 0.25, min = 0, max = 1')
    print(' -L : Coverage depth limit for variant calling'
          ' (int, default = 100, min = 1)')
    print(' -t : Length range tolerance for consensus sequences'
          ' (percentage, default=1, min=0, max=100)')
    print(' -T : Threads used for contig/scaffold alignments'
          ' (int, default = 1, min = 1)')
    print(' -g : Disable garbage collection and retain intermediate analysis files'
          ' (no value provided)')
    print()


def parse_args(args: list, version: str):
    ''' Parse command line input, ensure that required arguments have been
    provided with values, ensure that provided values are of the correct
    type and fall within accepted parameters. '''
    ''' Check that each argument has been set only once. '''
    arg_names = set(arg for arg in args if arg[0] == '-')
    multi_set_args = set()
    for arg in arg_names:
        if args.count(arg) > 1:
            multi_set_args.add(arg)
    if multi_set_args != set():
        multi_set_args = ', '.join(sorted(multi_set_args))
        print('\nERROR: The following arguments have been set more than once: '
              f'{multi_set_args}\nAnalysis aborted.\n')
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
    ''' Set expected types, default values, min values, and max values. '''
    required_args = {'-n', '-f', '-r', '-d'}
    arg_value_types = {'-n': str, '-f': str, '-r': str, '-d': str, '-i': float,
                       '-l': int, '-m': bool, '-D': int, '-q': int, '-v': float,
                       '-V': float, '-L': int, '-t': float, '-T': int,
                       '-g': bool}
    min_arg_values = {'-i': 0, '-l': 32, '-D': 1, '-q': 0, '-v': 0, '-V': 0,
                      '-L':1, '-t': 0, '-T': 1}
    max_arg_values = {'-i': 100, '-v': 1, '-V': 1, '-t': 100}
    default_arg_values = {'-i': 90, '-l': 50, '-m': False, '-D': 20, '-q': 20,
                          '-v': 0.75, '-V': 0.25, '-L': 100, '-t': 1, '-T': 1,
                          '-g': True}
    ''' Set garbage collection flag to boolean value if not provided. '''
    if '-g' not in arg_values:
        arg_values['-g'] = True
    elif '-g' in arg_values and arg_values['-g'] == '':
        arg_values['-g'] = False
    elif '-g' in arg_values and arg_values['-g'] != '':
        print('\nERROR: Argument -g should be provided without a value.'
              '\nAnalysis aborted.\n')
        print_usage(version)
        exit(2)
    ''' Set allow mixtures flag to boolen value if not provided. '''
    if '-m' not in arg_values:
        arg_values['-m'] = False
    elif '-m' in arg_values and arg_values['-m'] == '':
        arg_values['-m'] = True
    elif '-m' in arg_values and arg_values['-m'] != '':
        print('\nERROR: Argument -m should be provided without a value.'
              '\nAnalysis aborted.\n')
        print_usage(version)
        exit(3)
    ''' Check if all required arguments were provided. '''
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys():
            missing_args.add(required_arg)
    if missing_args != set():
        missing_args = ', '.join(sorted(missing_args))
        print(f'\nERROR: Values must be provided for the argument following '
              f'arguments: {missing_args}\nAnalysis aborted.\n')
        print_usage(version)
        exit(4)
    ''' Check if unrecognized arguments were provided. '''
    recognized_args = set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args.add(provided_arg)
    if unrecognized_args != set():
        unrecognized_args = ', '.join(sorted(unrecognized_args))
        print(f'\nERROR: The following provided arguments are not recognized: '
              f'{unrecognized_args}\nAnalysis aborted.\n')
        print_usage(version)
        exit(5)
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
                          f'value.\nAnalysis aborted.\n') 
                print_usage(version)
                exit(6)
    ''' Check if provided values are within the correct range. '''
    for arg, value in arg_values.items():
        bad_values = set()
        if arg in min_arg_values.keys() and value < min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must not be less than '
                  f'{min_arg_values[arg]}\nAnalysis aborted.\n')
            print_usage(version)
            exit(7)
        if arg in max_arg_values.keys() and value > max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must not exceed '
                  f'{max_arg_values[arg]}\nAnalysis aborted.\n')
            print_usage(version)
            exit(8) 
    ''' Check if arguments were provided without values. '''
    empty_args = set()
    for arg, value in arg_values.items():
        if value == '':
            empty_args.add(arg)
    if empty_args != set():
        empty_args = ', '.join(sorted(empty_args))
        print(f'\nERROR: The following arguments were provided without values: '
              f'{empty_args}\nAnalysis aborted.\n')
        print_usage(version)
        exit(9)
    ''' Check if variant masking threshold is lower than variant calling
    threshold. '''
    if '-v' in arg_values and '-V' in arg_values:
        if arg_values['-V'] >= arg_values['-v']:
            print(f'\nERROR: Variant allele fraction for masking ambiguous '
                  f'variants must be lower than variant allele fraction for '
                  f'calling variants.\nAnalysis aborted.\n')
            print_usage(version)
            exit(10)
    ''' Assign default values to unspecified arguments. '''
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value
    ''' Check provided argument values for reserved characters. '''
    reserved_chars = '*;\&|()<>'
    if any(char in str(value) for char in reserved_chars
           for value in arg_values.values()):
        print('\nERROR: Argument values cannot contain any of the following '
              f'reserved characters: *;\&|()<>')
        exit(11)
    ''' Return keyword args and their values. '''
    return arg_values


def check_input_files(fwd_reads: os.path, rev_reads: os.path, db: os.path):
    ''' Check that all provided input files exist. '''
    for file in [fwd_reads, rev_reads, db]:
        if not os.path.isfile(file):
            print(f'\nERROR: Input file {file} does not exist.'
                  '\nAnalysis aborted.\n')
            exit(12)


def make_output_dir(out_name: str):
    ''' Create a directory for output after checking that it does not
    already exist. '''
    if os.path.exists(out_name) == False:
        os.mkdir(out_name)
        logs_path = os.path.join(out_name, 'logs')
        os.mkdir(logs_path)
    else:
        print(f'\nERROR: Output directory {out_name} already exists!'
              '\nAnalysis aborted.\n')
        exit(13)


def check_database(db: os.path):
    ''' Checks the contents of the provided reference sequence database to
    ensure proper header formatting and unambiguous sequences. Returns a
    DataFrame with the min and max lengths for each species/segment/subtype. '''
    print('Checking reference sequence database...')
    seqs = read_fasta(db)
    species = []
    segments = []
    subtypes = []
    lengths = []
    for header, seq in seqs.items():
        if len(header.split('|')) < 5:
            print(f'\nERROR: The following header does not contain the minimum '
                  f'number of |-delimited fields:\n{header}'
                  f'\nAnalysis aborted.\n')
            exit(14)
        seq_species, segment, subtype = header.split('|')[2:5] 
        species.append(seq_species)
        segments.append(segment)
        subtypes.append(subtype)
        lengths.append(len(seq))
    db_data = pd.DataFrame()
    db_data['species'] = species
    db_data['segment'] = segments
    db_data['subtype'] = subtypes
    db_data['length'] = lengths
    cols = ['species', 'segment', 'subtype', 'length']
    group_cols = ['species', 'segment', 'subtype']
    min_lengths = db_data[cols].groupby(group_cols).min().reset_index()
    min_lengths.columns = ['species', 'segment', 'subtype', 'min_length']
    max_lengths = db_data[cols].groupby(group_cols).max().reset_index()
    max_lengths.columns = ['species', 'segment', 'subtype', 'max_length']
    merge_cols = ['species', 'segment', 'subtype']
    length_ranges = pd.merge(min_lengths, max_lengths, on=merge_cols)
    indices = ['species', 'segment', 'subtype']
    length_ranges = length_ranges.set_index(indices)
    return length_ranges


def read_fasta(file_path: os.path):
    ''' A generalized function for reading the contents of a FASTA file into
    a dict of header: sequence. '''
    if not os.path.exists(file_path):
        print('\nERROR: File {file_path} does not exist!\nAnalysis aborted.\n')
        exit(15)
    with open(file_path, 'r') as input_file:
        seqs = {}
        header = ''
        seqs[header] = ''
        headers = []
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                seqs[header] = ''
                headers.append(header)
            else:
                seqs[header] += line.strip()
    headers = Counter(headers)
    for header in headers:
        if headers[header] > 1:
            print(f'\nERROR: The following header is used for multiple sequences:'
                  f'\n{header}\nAnalysis aborted.\n')
            exit(16)
    if seqs[''] != '':
        print(f'\nERROR: FASTA file {file_path} Contains sequence before first '
              f'header!\nAnalysis aborted.\n')
        exit(17)
    seqs = {header: seq for header, seq in seqs.items() if header != ''}
    return seqs


def run(out_name: str, command: str, process_name: str, error_code: int):
    ''' A generalized function for running subprocesses, logging their output,
    and trapping erroneous exit statuses. '''
    stdout_file = os.path.join(out_name, 'logs', f'{process_name}_stdout.txt')
    stderr_file = os.path.join(out_name, 'logs', f'{process_name}_stderr.txt')
    stdout_file = open(stdout_file, 'w')
    stderr_file = open(stderr_file, 'w')
    complete_process = sp.run(command, stdout=stdout_file, stderr=stderr_file,
                              shell=True)
    stdout_file.close()
    stderr_file.close()
    return_code = complete_process.returncode
    if return_code != 0:
        print(f'\nERROR: Subprocess {process_name} failed '
              f'(Exit status: {return_code})\nAnalysis aborted.\n')
        exit(error_code)


def assemble_contigs(out_name: str, fwd_reads: os.path, rev_reads: os.path):
    ''' Reads are assembled de novo into contigs using SPAdes. '''
    print('Assembling reads into contigs...')
    spades_output = os.path.join(out_name, 'spades_output')
    command = (f'spades.py --rnaviral --isolate -1 {fwd_reads} -2 {rev_reads} '
               f'-o {spades_output}')
    process_name = 'spades'
    error_code = 18
    run(out_name, command, process_name, error_code)
    if not os.path.isfile(os.path.join(spades_output, 'contigs.fasta')):
        print('\nERROR: No contigs assembled!\nAnalysis aborted.\n')
        error_code = 19
        exit(error_code)


def align_contigs_to_ref_seqs(out_name: str, db: os.path, num_threads: int,
                              identity: float, length: int):
    ''' Contigs are aligned to reference sequences using BLASTn. Initial
    alignment filtering based on identity and length of alignment. '''
    print('Aligning contigs to reference sequences...')
    if any([os.path.exists(db + '.' + suffix) == False
            for suffix in ['nhr', 'nin' , 'nsq']]):
        command = f'makeblastdb -in {db} -dbtype nucl'
        process_name = 'makeblastdb_contigs'
        error_code = 20
        run(out_name, command, process_name, error_code)
    blast_output = os.path.join(out_name, 'contigs_blast.tsv')
    contigs = os.path.join(out_name, 'spades_output', 'contigs.fasta')
    cols = 'qseqid sseqid pident length bitscore sstart send qseq sseq slen'
    command = (f'blastn -query {contigs} -db {db} -num_threads {num_threads} '
               f'-outfmt "6 {cols}" > {blast_output}')
    process_name = 'blastn_contigs'
    error_code = 21
    run(out_name, command, process_name, error_code)
    blast_results = pd.read_csv(blast_output, names=cols.split(' '), sep='\t')
    blast_results = blast_results[blast_results['pident']>=identity]
    blast_results = blast_results[blast_results['length']>=length]
    if len(blast_results) == 0:
        print(f'\nERROR: No contigs aligned to reference sequences!'
              f'\nAnalysis aborted.\n')
        exit(22)
    return blast_results


def filter_contig_alignments(blast_results: pd.DataFrame):
    print('Filtering contig alignments...')
    ''' Alignments of contigs to reference sequences are filtered to discard
    contigs whose best matches (by bitscore) are not composed of multiple
    species, segments, or subtypes (other than "none"). '''
    ''' Get best matches for each contig. '''
    cols = ['qseqid', 'sseqid', 'bitscore']
    group_cols = ['qseqid', 'sseqid']
    combined_bitscores = blast_results[cols].groupby(group_cols).sum()
    combined_bitscores = combined_bitscores.reset_index()
    cols = ['qseqid', 'bitscore']
    group_cols = ['qseqid']
    max_combined_bitscores = combined_bitscores[cols].groupby(group_cols).max()
    max_combined_bitscores = max_combined_bitscores.reset_index()
    merge_cols = ['qseqid', 'bitscore']
    blast_results = pd.merge(blast_results, max_combined_bitscores, on=merge_cols)
    ''' Annotate subjects with species, segment, and subtype. '''
    annots = blast_results[['sseqid']].drop_duplicates()
    get_species = lambda row: row['sseqid'].split('|')[2]
    annots['species'] = annots.apply(get_species, axis=1)
    get_segment = lambda row: row['sseqid'].split('|')[3]
    annots['segment'] = annots.apply(get_segment, axis=1)
    get_subtype = lambda row: row['sseqid'].split('|')[4]
    annots['subtype'] = annots.apply(get_subtype, axis=1)
    blast_results = pd.merge(blast_results, annots, on='sseqid')
    ''' Discard contigs whose best matches are multiple species. '''
    cols = ['qseqid', 'species']
    num_species = blast_results[cols].drop_duplicates()
    group_cols = ['qseqid']
    num_species = num_species.groupby(group_cols).count().reset_index()
    num_species = num_species[num_species['species']==1][['qseqid']]
    blast_results = pd.merge(blast_results, num_species, on='qseqid')
    ''' Discard contigs whose best matches are multiple segments. '''
    cols = ['qseqid', 'segment']
    num_segments = blast_results[cols].drop_duplicates()
    group_cols = ['qseqid']
    num_segments = num_segments.groupby(group_cols).count().reset_index()
    num_segments = num_segments[num_segments['segment']==1][['qseqid']]
    blast_results = pd.merge(blast_results, num_segments, on='qseqid')
    ''' Discard contigs whose best matches are multiple subtypes. '''
    cols = ['qseqid', 'subtype']
    num_subtypes = blast_results[cols].drop_duplicates()
    num_subtypes['subtype'] = num_subtypes['subtype'].replace('none', np.nan)
    group_cols = ['qseqid']
    num_subtypes = num_subtypes.groupby(group_cols).count().reset_index()
    num_subtypes = num_subtypes[num_subtypes['subtype']<=1][['qseqid']]
    blast_results = pd.merge(blast_results, num_subtypes, on='qseqid')
    if len(blast_results) == 0:
        print(f'\nERROR: No contigs aligned to reference sequences!'
              f'\nAnalysis aborted.\n')
        exit(22)
    return blast_results


def mixed_infection_check(blast_results: pd.DataFrame, allow_mixtures: bool):
    ''' Filtered alignments of contigs to reference sequences are analyzed
    to determine if multiple species or subtypes are present. If analysis of
    mixed infections is allowed (-m), scaffols/mapping refs/consensus seqs are
    generated only for species/segments with subtype other than "none". If
    analysis of mixed infections is not allowed, an error is raised and analysis
    aborts. '''
    print('Checking for mixed infection...')
    ''' Check for mixture of species. '''
    species = blast_results['species'].unique()
    num_species = len(species)
    if num_species > 1:
        species = ', '.join(sorted(species))
        if not allow_mixtures:
            print(f'\nERROR: Multiple species detected ({species})!'
                  f'\nAnalysis aborted.\n')
            exit(23)
        print(f'\nWARNING: Multiple species detected ({species})!'
              f'\nAnalysis continuing...\n')
    ''' Check if any species/segment combination has mixed subtypes. '''
    cols = ['species', 'segment', 'subtype']
    mixed_segments = blast_results[cols].drop_duplicates()
    mixed_segments['subtype'] = mixed_segments['subtype'].replace('none', np.nan)
    group_cols = ['species', 'segment']
    mixed_segments = mixed_segments.groupby(group_cols).count().reset_index()
    mixed_segments = mixed_segments[mixed_segments['subtype']>1]
    mixed_segments = mixed_segments[['species', 'segment']]
    merge_cols = ['species', 'segment']
    mixed_segments = pd.merge(blast_results, mixed_segments, on=merge_cols)
    if len(mixed_segments) > 0:
        if not allow_mixtures:
            print('\nERROR: Multiple subtypes detected for the following '
                  f'species/segments:')
        else:
            print('\nWARNING: Multiple subtypes detected for the following '
                  f'species/segments:')            
        for species in sorted(mixed_segments['species'].unique()):
            species_results = mixed_segments[mixed_segments['species']==species]
            segment_order = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
            for segment in sorted(species_results['segment'].unique(),
                                  key=lambda seg: segment_order.index(seg)):
                segment_results = species_results[species_results['segment']==segment]
                subtypes = ', '.join(sorted(segment_results['subtype'].unique()))
                print(f' {species}/{segment}: ({subtypes})')
        if not allow_mixtures:
            print('\nAnalysis aborted.\n')
            exit(24)
        else:
            print('\nAnalysis continuing...\n')
    ''' Discard alignments for species/segments without defined subtype if mixed
    infection.  '''
    if num_species > 1 or len(mixed_segments) > 0:
        blast_results = blast_results[blast_results['subtype']!='none']
    if len(blast_results) == 0:
        print(f'\nERROR: No contigs aligned to reference sequences!'
              f'\nAnalysis aborted.\n')
        exit(22)
    return blast_results


def make_scaffolds(out_name: str, blast_results: pd.DataFrame):
    ''' Generate a scaffold sequence for each species/segment/subtype by
    collapsing all contigs describing that species/segment/subtype into a single
    sequence. '''
    print('Making scaffolds...')
    scaffolds = os.path.join(out_name, f'{out_name}_scaffolds.fa')
    with open(scaffolds, 'w') as output_file:
        cols = ['species', 'segment', 'subtype']
        scaffold_types = blast_results[cols].drop_duplicates()
        for index, row in scaffold_types.iterrows():
            species = row['species']
            scaffold_data = blast_results[blast_results['species']==species]
            segment = row['segment']
            scaffold_data = scaffold_data[scaffold_data['segment']==segment]
            subtype = row['subtype']
            scaffold_data = scaffold_data[scaffold_data['subtype']==subtype]
            print(f' {species}/{segment} ({subtype})')
            scaffold_data = pick_scaffolding_ref(scaffold_data)
            scaffold_seq = make_scaffold_seq(out_name, species, segment,
                                             subtype, scaffold_data)
            header = f'>{out_name}|{species}|{segment}|{subtype}|'
            output_file.write(header + '\n')
            output_file.write(scaffold_seq + '\n')
            

def pick_scaffolding_ref(scaffold_data: pd.DataFrame):
    ''' Pick a reference sequence for a species/segment/subtype to use for
    scaffolding. Contig alignments against this sequence are used to roughly
    position the contig inside the segment before conducting a multiple sequence
    alignment on the contigs. The scaffolding reference is chosen by selecting
    the reference seq with the most positions covered by contigs, followed by the
    reference seq with the most identity to its contigs. If there are multiple
    options remaining, the longest is chosen, followed by the first
    alphabetically. '''
    print('  Picking reference sequence for scaffolding...')
    ''' Find ref seq(s) most covered by contigs. '''
    def count_cov_pos(data_frame):
        cov_positions = set()
        for index, row in data_frame.iterrows():
            start = min([row['sstart'], row['send']])
            end = max([row['sstart'], row['send']])
            cov_positions = cov_positions.union(set(range(start, end + 1)))
        return len(cov_positions)
    cols = ['sseqid', 'sstart', 'send']
    cov_pos = scaffold_data[cols].drop_duplicates()
    group_cols = ['sseqid']
    cov_pos = cov_pos.groupby(group_cols).apply(count_cov_pos).reset_index()
    cov_pos.columns = ['sseqid', 'covered_pos']
    max_cov_pos = cov_pos['covered_pos'].max()
    max_cov_pos = cov_pos[cov_pos['covered_pos']==max_cov_pos][['sseqid']]
    scaffold_data = pd.merge(scaffold_data, max_cov_pos, on='sseqid')
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
    cols = ['sseqid', 'sstart', 'send', 'qseq', 'sseq']
    ident_pos = scaffold_data[cols].drop_duplicates()
    group_cols = ['sseqid']
    ident_pos = ident_pos.groupby(group_cols).apply(count_id_pos).reset_index()
    ident_pos.columns = ['sseqid', 'identical_pos']
    max_ident_pos = ident_pos['identical_pos'].max()
    max_ident_pos = ident_pos[ident_pos['identical_pos']==max_ident_pos]
    max_ident_pos = max_ident_pos[['sseqid']]
    scaffold_data = pd.merge(scaffold_data, max_cov_pos, on='sseqid')
    ''' Take longest remaining ref seq(s) for each segment/subtype. '''
    max_slen = scaffold_data['slen'].max()
    scaffold_data = scaffold_data[scaffold_data['slen']==max_slen]
    ''' Take first alphabetical remaining ref seq for each segment/subtype. '''
    first_alpha = scaffold_data['sseqid'].min()
    scaffold_data = scaffold_data[scaffold_data['sseqid']==first_alpha]
    return scaffold_data


def make_scaffold_seq(out_name: str, species: str, segment: str, subtype: str,
                      scaffold_results: pd.DataFrame):
    ''' Generate the scaffold sequence for a species/segment/subtype. Contig
    alignments against the chosen reference seq are used to roughly position
    each contig with the segment (a number of terminal Ns are added to each
    end of each contig based on the start and end coordinates of the contig's
    alignment). At this time, any part of a contig extending past the scaffolding
    reference is trimmed. A mutliple sequence alignment is then performed on the
    trimmed, positioned contigs. The scaffold sequence is generated from the
    consensus of that multiple sequence alignment. Ambiguous positions (the
    dominant base does not account for at least 50% of bases) are replaced with
    Ns. Positions not covered by any contigs are also replaced with Ns in the
    scaffold. '''
    print('  Trimming and positioning contigs...')
    ''' Get alignment start and end coordinates. '''
    get_start = lambda row: min([row['sstart'], row['send']])
    scaffold_results['start'] = scaffold_results.apply(get_start, axis=1)
    get_end = lambda row: max([row['sstart'], row['send']])
    scaffold_results['end'] = scaffold_results.apply(get_end, axis=1)
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
    rev_comp_seq = lambda seq: ''.join(rev_comp_bases[base] for base in seq[::-1])
    def flip_qseq(row):
        if row['sstart'] > row['send']:
            return rev_comp_seq(row['qseq'])
        else:
            return row['qseq']
    scaffold_results['qseq'] = scaffold_results.apply(flip_qseq, axis=1)
    def flip_sseq(row):
        if row['sstart'] > row['send']:
            return rev_comp_seq(row['sseq'])
        else:
            return row['sseq']
    scaffold_results['sseq'] = scaffold_results.apply(flip_sseq, axis=1)
    ''' Write all contigs to a FASTA file. Trim contigs based on their alignments
    to reference sequences. Also add leading and trailing Ns to contig so that
    it is properly positioned within the reference sequence. '''
    file_name = f'{species}_{segment}_{subtype}_contigs.fa'
    contigs = os.path.join(out_name, file_name)
    with open(contigs, 'w') as contig_file:
        contig_counter = 1
        for index, row in scaffold_results.iterrows():
            header = f'>{species}_{segment}_{subtype}_contig_{contig_counter}\n'
            contig_file.write(header)
            seq = 'N' * (row['start'] - 1)
            seq += row['qseq'].replace('-', '')
            seq += ('N' * (row['slen'] - row['end']))
            contig_file.write(seq + '\n')
            contig_counter += 1
    ''' Generate multiple sequence alignments of trimmed/positioned contigs. '''
    print('  Collapsing contigs...')
    file_name = f'{species}_{segment}_{subtype}_contigs.afa'
    aligned_contigs = os.path.join(out_name, file_name)
    if contig_counter > 2:
        command = (f'clustalw -INFILE={contigs} -OUTFILE={aligned_contigs} '
                   f'-OUTPUT=FASTA')
        process_name = f'clustalw_{segment}'
        error_code = 25
        run(out_name, command, process_name, error_code)
    else:
        shutil.copyfile(contigs, aligned_contigs)
    ''' Load aligned contigs and replace leading and trailing Ns and gaps with
    dots so that they are ignored when determining consensus bases. '''
    seqs = read_fasta(aligned_contigs)
    clean_seqs = []
    for seq in seqs.values():
        head_len = len(seq) - len(seq.lstrip('N-'))
        tail_len = len(seq) - len(seq.rstrip('N-'))
        seq = seq.strip('N-')
        seq = ('.' * head_len) + seq
        seq += ('.' * tail_len)
        clean_seqs.append(seq)
    ''' Check that all seqs in the multiple seq alignment are the same length. '''
    alignment_lengths = set(len(seq) for seq in clean_seqs)
    if len(alignment_lengths) > 1:
        print(f'\nERROR: Multiple sequence alignment for {species}/{segment} '
              f'({subtype}) generated unequal alignment lengths! '
              f'\nAnalysis aborted.\n')
        error_code = 26
        exit(error_code)
    ''' Create consensus sequence of multiple seq alignments, ie the scaffold. '''
    print('  Generating scaffold...')
    alignment_length = list(alignment_lengths)[0]
    scaffold_seq = ''
    for i in range(alignment_length):
        bases = Counter(seq[i] for seq in clean_seqs if seq[i] != '.')
        if bases.most_common(1) == []:
            scaffold_seq += 'N'
        elif bases.most_common(1)[0][1] / len(bases) > 0.5:
            scaffold_seq += bases.most_common(1)[0][0]
        else:
            scaffold_seq += 'N'
    return scaffold_seq


def align_scaffolds_to_ref_seqs(out_name: str, db: os.path, num_threads: int):
    ''' Scaffolds are aligned to reference sequences using BLASTn. '''
    print('Aligning scaffolds to reference sequences...')
    if any([os.path.exists(db + '.' + suffix) == False
            for suffix in ['nhr', 'nin' , 'nsq']]):
        command = f'makeblastdb -in {db} -dbtype nucl'
        process_name = 'makeblastdb_scaffolds'
        error_code = 27
        run(out_name, command, process_name, error_code)
    with open(db) as input_file:
        num_db_seqs = sum(line[0] == '>' for line in input_file)
    blast_output = os.path.join(out_name, 'scaffolds_blast.tsv')
    scaffolds = os.path.join(out_name, f'{out_name}_scaffolds.fa')
    cols = 'qseqid sseqid bitscore sstart send qseq sseq'
    command = (f'blastn -query {scaffolds} -db {db} -num_threads {num_threads} '
               f'-max_target_seqs {num_db_seqs} -outfmt "6 {cols}" '
               f'> {blast_output}')
    process_name = 'blastn_scaffolds'
    error_code = 28
    run(out_name, command, process_name, error_code)
    blast_results = pd.read_csv(blast_output, names=cols.split(' '), sep='\t')
    if len(blast_results) == 0:
        print(f'\nERROR: No scaffolds aligned to reference sequences! '
              f'\nAnalysis aborted.\n')
        exit(29)
    return blast_results


def filter_scaffold_alignments(blast_results: pd.DataFrame):
    ''' For each scaffold, choose a reference sequence that will be used to fill
    in Ns. This scaffold with filled-in Ns is the mapping reference. The chosen
    reference sequences are each scaffold's best match (based on summed
    bitscores). If a scaffold has multiple best matches based on summed bitscores,
    choose the first alphabetically. '''
    print('Choosing reference sequences used for infilling scaffolds...')
    ''' Remove reversed alignments (they should be artefactual at this point). '''
    blast_results = blast_results[blast_results['send'] > blast_results['sstart']]
    ''' Get best matches for each scaffold. '''
    cols = ['qseqid', 'sseqid', 'bitscore']
    group_cols = ['qseqid', 'sseqid']
    combined_bitscores = blast_results[cols].groupby(group_cols).sum()
    combined_bitscores = combined_bitscores.reset_index()
    cols = ['qseqid', 'bitscore']
    group_cols = ['qseqid']
    max_combined_bitscores = combined_bitscores[cols].groupby(group_cols).max()
    max_combined_bitscores = max_combined_bitscores.reset_index()
    merge_cols = ['qseqid', 'bitscore']
    best_matches = pd.merge(combined_bitscores, max_combined_bitscores, on=merge_cols)
    best_matches = best_matches[['qseqid', 'sseqid']].drop_duplicates()
    blast_results = pd.merge(blast_results, best_matches, on=['qseqid', 'sseqid'])
    ''' Get first best-matching ref seq (alphabetical) in case of ties. '''
    cols = ['qseqid', 'sseqid']
    group_cols = ['qseqid']
    first_alpha = blast_results[cols].groupby(group_cols).min().reset_index()
    blast_results = pd.merge(blast_results, first_alpha, on=['qseqid', 'sseqid'])
    if len(blast_results) == 0:
        print(f'\nERROR: No scaffolds aligned to reference sequences! '
              f'\nAnalysis aborted.\n')
        exit(29)
    return blast_results


def make_mapping_refs(out_name: str, blast_results: pd.DataFrame, db: os.path):
    ''' Mapping references are created for each species/segment/subtype detected.
    These consist of the scaffold with all Ns in the scaffold filled-in using
    the corresponding positions from that scaffold's best-matching reference
    sequence. '''
    print('Creating mapping references...')
    ''' Create dict with best-matching ref seq for each segment. '''
    sseqids = blast_results['sseqid'].unique()
    best_ref_seqs = {seq_name: '' for seq_name in sseqids}
    with open(db, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
            elif header in best_ref_seqs:
                best_ref_seqs[header] += line.strip()
    ''' Create mapping ref for each species/segment/subtype. '''
    def make_map_ref(data_frame):
        data_frame = data_frame.sort_values(by=['sstart', 'send'],
                                            ascending=[True, False])
        ref_seq = best_ref_seqs[data_frame['sseqid'].min()]
        last_position = 0
        seq = ''
        for index, row in data_frame.iterrows():
            if row['sstart'] > last_position:
                seq += ref_seq[last_position:row['sstart'] - 1]
            if row['send'] > last_position:
                qseq = row['qseq'].upper()
                sseq = row['sseq'].upper()
                if row['sstart'] <= last_position:
                    start = (last_position - row['sstart']) + 1
                    qseq = qseq[start:]
                    sseq = sseq[start:]
                for qbase, sbase in zip(qseq, sseq):
                    if qbase in 'ATGC':
                        seq += qbase
                    else:
                        seq += sbase
                last_position = row['send']
        seq += ref_seq[last_position:].upper()
        seq = seq.replace('-', '')
        return seq
    cols = ['sseqid', 'sstart', 'send', 'qseq', 'sseq']
    group_cols = ['sseqid']
    blast_results = blast_results[cols]
    blast_results = blast_results.groupby(group_cols).apply(make_map_ref)
    blast_results = blast_results.reset_index()
    blast_results.columns = ['sseqid', 'mapping_seq']
    ''' Annotate sequences. '''
    get_species = lambda row: row['sseqid'].split('|')[2].split('_')[0]
    blast_results['species'] = blast_results.apply(get_species, axis=1)
    get_segment = lambda row: row['sseqid'].split('|')[3].split('_')[0]
    blast_results['segment'] = blast_results.apply(get_segment, axis=1)
    get_subtype = lambda row: row['sseqid'].split('|')[4].split('_')[0]
    blast_results['subtype'] = blast_results.apply(get_subtype, axis=1)
    ''' Write mapping refs to FASTA. '''
    segment_order = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    get_segment_order = lambda row: segment_order.index(row['segment'])
    blast_results['segment_sort'] = blast_results.apply(get_segment_order,
                                                        axis=1)
    sort_cols = ['species', 'segment_sort', 'subtype']
    blast_results = blast_results.sort_values(by=sort_cols)
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs, 'w') as output_file:
        for index, row in blast_results.iterrows():
            sseqid = row['sseqid'].split('|')[:5]
            accession, ref_name, species, segment, subtype = sseqid
            accession = accession.lstrip('>')
            header = [out_name, species, segment, subtype, accession, ref_name]
            header = '|'.join(header)
            header = f'>{header}|'
            seq = row['mapping_seq']
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def map_reads_to_mapping_refs(out_name: str, fwd_reads: os.path,
                              rev_reads: os.path, min_qual: float):
    ''' Input reads are mapped to the mapping references (make_mapping_refs func)
    using BWA mem. The alignment is filtered to retain only paired reads in their
    primary alignments, then the alignment is sorted and indexed. '''
    print('Mapping reads to mapping references...')
    ''' Index combined mapping references. '''
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    command = f'bwa index {mapping_refs}'
    process_name = 'bwa_index'
    error_code = 30
    run(out_name, command, process_name, error_code)
    ''' Align input reads to combined mapping references. '''
    alignment = os.path.join(out_name, 'alignment.sam')
    command = f'bwa mem {mapping_refs} {fwd_reads} {rev_reads} > {alignment}'
    process_name = 'bwa_mem'
    error_code = 31
    run(out_name, command, process_name, error_code)
    ''' Filter, sort, and index combined alignment. '''
    filtered_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    command = (f'samtools view -f 1 -F 2316 -q {min_qual} -h {alignment} '
               f'| samtools sort -o {filtered_alignment}')
    process_name = 'samtools_view_combined'
    error_code = 32
    run(out_name, command, process_name, error_code)
    command = f'samtools index {filtered_alignment}'
    process_name = 'samtools_index_combined'
    error_code = 33
    run(out_name, command, process_name, error_code)


def make_consensus_seqs(out_name: str, max_depth: int, min_qual: float,
                        min_depth: int, vaf_ambig: float, vaf_call: float,
                        length_ranges: pd.DataFrame, length_tolerance: float):
    ''' A consensus sequence is generated for each species/segment/subtype.
    First, variants are called from the mapped reads using FreeBayes, then
    ambiguous and low coverage positions are masked, then variant calls and
    masking are applied to the mapping reference. These operations are performed
    separately on each species/segment/subtype; FreeBayes has bizarre behaviour
    when assessing coverage/support for the 5' end of new chromosomes after the
    first chromosome, occasionally resulting in extensive masking of the 5' end
    of segments for "low coverage" despite adequate read depth. This behaviour
    has been reported to the FreeBayes developers (reported as issue #781 on the
    FreeBayes GitHub), so this is an admittedly hacky solution to this issue for
    the moment. Following individual consensus sequence generation, all
    consensus sequences are collected into a single file, and each consensus
    sequence is checked to make sure its length falls within the tolerated range
    (the length of the shortest/longest example of that species/segment/subtype
    in the provided database +/- the length_tolerance percentage). '''
    print('Making consensus sequences...')
    ''' Load headers of mapping references. '''
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs, 'r') as input_file:
        seq_names = set()
        for line in input_file:
            if line[0] == '>':
                seq_name = line.strip().lstrip('>')
                seq_names.add(seq_name)
    segments = 'PB2 PB1 PA HA NP NA M NS'.split(' ')
    get_order = lambda seq_name: (seq_name.split('|')[1],
                                  segments.index(seq_name.split('|')[2]),
                                  seq_name.split('|')[3])
    seq_names = sorted(seq_names, key=get_order)
    for seq_name in seq_names:
        species, segment, subtype = seq_name.split('|')[1:4]
        print(f' {species}/{segment} ({subtype})')
        subset_mapped_reads_and_mapping_ref(out_name, seq_name)
        make_freebayes_pileup(out_name, species, segment, subtype, max_depth,
                              min_qual)
        make_bed_of_masked_pos(out_name, seq_name, min_depth, vaf_ambig,
                               vaf_call)
        make_vcf_of_variant_pos(out_name, species, segment, subtype, min_depth,
                                vaf_call)
        make_consensus_seq(out_name, species, segment, subtype)
    ''' Gather all consensus seqs into one FASTA. '''
    seqs = {}
    consensus_seqs = os.path.join(out_name, f'{out_name}_consensus_seqs.fa')
    with open(consensus_seqs, 'w') as output_file:
        for seq_name in seq_names:
            species, segment, subtype = seq_name.split('|')[1:4]
            file_name = f'{species}_{segment}_{subtype}_consensus_seq.fa'
            consensus_seq = os.path.join(out_name, file_name)
            for header, seq in read_fasta(consensus_seq).items():
                header = '|'.join(header.split('|')[:4]) + '|'
                seqs[header] = seq
                output_file.write('>' + header + '\n')
                output_file.write(seq + '\n')
    ''' Check that consensus seq lengths are within tolerated range. '''
    for header, seq in seqs.items():
        species, segment, subtype = header.split('|')[1:4]
        index_values = species, segment, subtype
        min_length = length_ranges.loc[index_values]['min_length']
        max_length = length_ranges.loc[index_values]['max_length']
        min_length = round(min_length * ((100 - length_tolerance) / 100), 0)
        max_length = round(max_length * ((100 + length_tolerance) / 100), 0)
        if not (min_length <= len(seq) <= max_length):
            print(f'\nERROR: The consensus sequence generated for '
                  f'species/segment (subtype) {species}/{segment} ({subtype}) '
                  f'is not within the tolerated length range ({int(min_length)} '
                  f'to {int(max_length)} bases).\nAnalysis aborted.\n')
            exit(34)


def subset_mapped_reads_and_mapping_ref(out_name: str, seq_name: str):
    ''' Create new FASTA and BAM files for an individual species/segment/subtype
    by subseting the mapping references FASTA and alignment BAM. '''
    species, segment, subtype = seq_name.split('|')[1:4]
    ''' Subset mapped reads. '''
    full_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    file_name = f'{species}_{segment}_{subtype}_alignment.bam'
    subset_alignment = os.path.join(out_name, file_name)
    command = (f'samtools view -h {full_alignment} "{seq_name}" '
               f'| samtools sort -o {subset_alignment}')
    process_name = f'samtools_view_{species}_{segment}_{subtype}_subset'
    error_code = 35
    run(out_name, command, process_name, error_code)
    command = (f'samtools index {subset_alignment}')
    process_name = f'samtools_index_{species}_{segment}_{subtype}_subset'
    error_code = 36
    run(out_name, command, process_name, error_code)
    ''' Subset mapping reference FASTA. '''
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    mapping_refs = read_fasta(mapping_refs)
    file_name = f'{species}_{segment}_{subtype}_mapping_ref.fa'
    mapping_ref = os.path.join(out_name, file_name)
    with open(mapping_ref, 'w') as output_file:
        for header, seq in mapping_refs.items():
            if tuple(header.split('|')[1:4]) == (species, segment, subtype):
                output_file.write('>' + header + '\n')
                output_file.write(seq + '\n')


def make_freebayes_pileup(out_name: str, species: str, segment:str, subtype:str,
                          max_depth: int, min_qual: float):
    ''' Run FreeBayes on an individual species/segment/subtype. ''' 
    print(f'  Generating pileup...')
    file_name = f'{species}_{segment}_{subtype}_mapping_ref.fa'
    mapping_ref = os.path.join(out_name, file_name)
    file_name = f'{species}_{segment}_{subtype}_alignment.bam'
    subset_alignment = os.path.join(out_name, file_name)
    file_name = f'{species}_{segment}_{subtype}_pileup.vcf'
    pileup = os.path.join(out_name, file_name)
    command = (f'freebayes -f {mapping_ref} {subset_alignment} -p 1 '
               f'--limit-coverage {max_depth} --min-mapping-quality {min_qual} '
               f'--min-base-quality {min_qual} --pooled-continuous '
               f'--report-monomorphic --haplotype-length 0 '
               f'--use-best-n-alleles 4 --min-alternate-count 1 '
               f'--min-alternate-fraction 0 > {pileup}')
    process_name = f'freebayes_{species}_{segment}_{subtype}'
    error_code = 37
    run(out_name, command, process_name, error_code)


def parse_freebayes_pileup(out_name: str, species: str, segment: str,
                           subtype: str):
    ''' Parse the VCF output of FreeBayes, returning for each position the
    read depth, number of reads with the dominant variant allele, and total
    number of reads with variant alleles. '''
    file_name = f'{species}_{segment}_{subtype}_pileup.vcf'
    pileup = open(os.path.join(out_name, file_name), 'r')
    line = pileup.readline()
    while line != '' and line[0] == '#':
        if f'contig=<ID={out_name}|{species}|{segment}|{subtype}|' in line:
            seq_length = int(line.split(',length=')[1].split('>')[0])
        line = pileup.readline()
    last_position = 0
    last_depth = 0
    last_max_alt_depth = 0
    last_total_alt_depth = 0
    while line != '':
        fields = line.strip().split('\t')
        fields = (fields[0], int(fields[1]), fields[3],
                  fields[4], fields[8], fields[9])
        seq_name, position, ref, alt, keys, values = fields
        if position != last_position + 1:
            for int_pos in range(last_position + 1, position):
                yield int_pos, last_depth, last_max_alt_depth, last_total_alt_depth
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
        yield position, total_depth, max_alt_depth, total_alt_depth
        last_position = position
        last_depth = total_depth
        last_max_alt_depth = max_alt_depth
        last_total_alt_depth = total_alt_depth
        line = pileup.readline()
    if last_position < seq_length:
        for int_pos in range(last_position + 1, seq_length + 1):
            yield int_pos, last_depth, last_max_alt_depth, last_total_alt_depth
    pileup.close()


def get_low_cov_pos(out_name: str, species: str, segment: str, subtype: str,
                    min_depth: int):
    ''' Returns a set containing all the positions with low coverage
    (total read coverage less than min_depth). '''
    low_cov_pos = set()
    for line in parse_freebayes_pileup(out_name, species, segment, subtype):
        position, total_depth, max_vaf, total_vaf = line
        if total_depth < min_depth:
            low_cov_pos.add(position)
    return low_cov_pos


def get_ambig_pos(out_name: str, species: str, segment: str, subtype: str,
                  min_depth: int, vaf_ambig: float, vaf_call: float):
    ''' Returns a set containing all the ambiguous positions (proportion of
    reads with variant alleles exceeds vaf_ambig, but the proportion of reads
    with the dominant variant allele does not exceed vaf_call). '''
    ambig_pos = set()
    for line in parse_freebayes_pileup(out_name, species, segment, subtype):
        position, total_depth, max_vaf, total_vaf = line
        if (total_depth >= min_depth and
            total_vaf / total_depth >= vaf_ambig and
            max_vaf / total_depth < vaf_call):
            ambig_pos.add(position)
    return ambig_pos


def get_variant_pos(out_name: str, species: str, segment: str, subtype: str,
                    min_depth: int, vaf_ambig: float, vaf_call: float):
    ''' Returns a set containing all the variant positions (proportion of
    reads with the dominant variant allele exceeds vaf_call). '''
    variant_pos = set()
    for line in parse_freebayes_pileup(out_name, species, segment, subtype):
        position, total_depth, max_vaf, total_vaf = line
        if total_depth >= min_depth and max_vaf / total_depth >= vaf_call:
            variant_pos.add(position)
    return variant_pos


def get_spans(positions):
    ''' Takes a set of positions, identifies all spans of consecutive
    positions. Spans are returned as tuples of (first position in span, last
    position in span). '''
    spans = set()
    positions = sorted(positions)
    if len(positions) > 0:
        span_start = positions[0]
        last_position = positions[0]
        for position in positions[1:]:
            if position != last_position + 1:
                span_end = last_position
                spans.add((span_start, span_end))
                span_start = position
            last_position = position
        span_end = last_position
        spans.add((span_start, span_end))
    spans = sorted(spans)
    return spans


def make_bed_of_masked_pos(out_name: str, seq_name: str, min_depth: int,
                           vaf_ambig: float, vaf_call: float):
    ''' Make a BED file of positions to mask in the consensus sequence. Use
    the zero-indexed BedGraph format. '''
    print('  Masking low coverage and ambiguous positions...')
    species, segment, subtype = seq_name.split('|')[1:4]
    ''' Parse FreeBayes pileup to identify positions to be masked due to low
    coverage or ambiguous base call. '''
    low_cov_pos = get_low_cov_pos(out_name, species, segment, subtype,
                                  min_depth)
    ambig_pos = get_ambig_pos(out_name, species, segment, subtype, min_depth,
                              vaf_ambig, vaf_call)
    masked_pos = low_cov_pos.union(ambig_pos)
    masked_spans = get_spans(masked_pos)
    ''' Write spans of masked positions to BED file in BedGraph format. '''
    file_name = f'{species}_{segment}_{subtype}_masked.bed'
    masked_pos = os.path.join(out_name, file_name)
    with open(masked_pos, 'w') as output_file:
        for (start, end) in masked_spans:
            line = [seq_name, start - 1, end, 0]
            line = '\t'.join(str(i) for i in line)
            output_file.write(line + '\n')


def make_vcf_of_variant_pos(out_name: str, species: str, segment: str,
                            subtype: str, min_depth: int, vaf_call: float):
    ''' Filter FreeBayes output to only output variant positions where the
    proportion of reads with the dominant variant allele exceeds vaf_call. '''
    print('  Calling variant positions...')
    ''' Parse FreeBayes pileup and write high-confidence variant lines to
    VCF file. '''
    file_name = f'{species}_{segment}_{subtype}_pileup.vcf'
    pileup = open(os.path.join(out_name, file_name), 'r')
    file_name = f'{species}_{segment}_{subtype}_variants.vcf'
    variants = open(os.path.join(out_name, file_name), 'w')
    line = pileup.readline()
    while line != '' and line[0] == '#':
        variants.write(line)
        line = pileup.readline()
    while line != '':
        fields = line.strip().split('\t')
        fields = (fields[0], int(fields[1]), fields[3],
                  fields[4], fields[8], fields[9])
        seq_name, position, ref, alt, keys, values = fields
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
        if (total_depth >= min_depth and
            max_alt_depth / total_depth >= vaf_call):
            variants.write(line)
        line = pileup.readline()
    pileup.close()
    variants.close()
    ''' Zip and index VCF file. '''
    file_name = f'{species}_{segment}_{subtype}_variants.vcf'
    variants = os.path.join(out_name, file_name)
    file_name = f'{species}_{segment}_{subtype}_variants.bcf'
    zipped_variants = os.path.join(out_name, file_name)
    command = f'bcftools view {variants} -Ob -o {zipped_variants}'
    process_name = 'bcftools_view_{species}_{segment}_{subtype}'
    error_code = 38
    run(out_name, command, process_name, error_code)
    command = f'bcftools index {zipped_variants}'
    process_name = 'bcftools_index_{species}_{segment}_{subtype}'
    error_code = 39
    run(out_name, command, process_name, error_code)


def make_consensus_seq(out_name: str, species: str, segment: str, subtype: str):
    ''' Apply variant positions VCF and masked positions BED to mapping
    reference to generate consensus sequence for an individual
    species/segment/subtype. '''
    print('  Generating consensus sequence...')
    ''' Apply variants to mapping ref. '''
    file_name = f'{species}_{segment}_{subtype}_mapping_ref.fa'
    mapping_ref = os.path.join(out_name, file_name)
    file_name = f'{species}_{segment}_{subtype}_masked.bed'
    masked = os.path.join(out_name, file_name)
    file_name = f'{species}_{segment}_{subtype}_variants.bcf'
    zipped_variants = os.path.join(out_name, file_name)
    file_name = f'{species}_{segment}_{subtype}_consensus_seq.fa'
    consensus_seq = os.path.join(out_name, file_name)
    command = (f'cat {mapping_ref} | bcftools consensus -m {masked} '
               f'{zipped_variants} > {consensus_seq}')
    process_name = 'bcftools_consensus_{species}_{segment}_{subtype}'
    error_code = 40
    run(out_name, command, process_name, error_code)


def write_report(out_name:str, min_depth: int, vaf_ambig: float,
                 vaf_call: float):
    ''' Generate a report with metrics about each species/segment/subtype in
    the library. '''
    print('Writing report...')
    ''' Count reads mapped to each segment and add to report. '''
    filtered_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    reads_mapped = os.path.join(out_name, 'reads_mapped.tsv')
    command = f'samtools idxstats {filtered_alignment} > {reads_mapped}'
    process_name = 'samtools_idxstats'
    error_code = 41
    run(out_name, command, process_name, error_code)
    cols = 'seq_name seq_length reads_mapped reads_unmapped'.split(' ')
    reads_mapped = pd.read_csv(reads_mapped, sep='\t', names=cols)
    reads_mapped = reads_mapped.replace('*', np.nan).dropna()
    get_seq_name = lambda row: '|'.join(row['seq_name'].split('|')[:4]) + '|'
    reads_mapped['seq_name'] = reads_mapped.apply(get_seq_name, axis=1) 
    cols = ['seq_name', 'reads_mapped', 'seq_length']
    report = reads_mapped[cols].drop_duplicates()
    ''' Get scaffold seq completeness. '''
    scaffolds = os.path.join(out_name, f'{out_name}_scaffolds.fa')
    scaffolds = read_fasta(scaffolds)
    scaffolds = {'|'.join(header.split('|')[:4]) + '|': seq
                 for header, seq in scaffolds.items()}
    scaffolds = {header: sum(char in 'ATGC' for char in seq)
                 for header, seq in scaffolds.items()}
    report['scaffold_completeness'] = report['seq_name'].map(scaffolds)
    get_perc = lambda row: row['scaffold_completeness'] * 100 / row['seq_length']
    report['scaffold_completeness'] = report.apply(get_perc, axis=1)
    report['scaffold_completeness'] = round(report['scaffold_completeness'], 1)
    ''' Get consensus seq completeness. '''
    consensus_seqs = os.path.join(out_name, f'{out_name}_consensus_seqs.fa')
    consensus_seqs = read_fasta(consensus_seqs)
    consensus_seqs = {'|'.join(header.split('|')[:4]) + '|': seq
                      for header, seq in consensus_seqs.items()}
    consensus_seqs = {header: sum(char in 'ATGC' for char in seq)
                      for header, seq in consensus_seqs.items()}
    report['consensus_completeness'] = report['seq_name'].map(consensus_seqs)
    get_perc = lambda row: row['consensus_completeness'] * 100 / row['seq_length']
    report['consensus_completeness'] = report.apply(get_perc, axis=1)
    report['consensus_completeness'] = round(report['consensus_completeness'], 1)
    ''' Count low cov positions. '''
    low_cov_pos = {}
    for seq_name in report['seq_name'].unique():
        species, segment, subtype = seq_name.split('|')[1:4]
        low_cov_pos[seq_name] = get_low_cov_pos(out_name, species, segment,
                                                subtype, min_depth)
    low_cov_pos = {seq_name: len(pos) for seq_name, pos in low_cov_pos.items()}
    report['low_cov_pos'] = report['seq_name'].map(low_cov_pos)
    get_perc = lambda row: row['low_cov_pos'] * 100 / row['seq_length']
    report['low_cov_perc'] = report.apply(get_perc, axis=1)
    report['low_cov_perc'] = round(report['low_cov_perc'], 1)
    ''' Count ambig positions. '''
    ambig_pos = {}
    for seq_name in report['seq_name'].unique():
        species, segment, subtype = seq_name.split('|')[1:4]
        ambig_pos[seq_name] = get_ambig_pos(out_name, species, segment, subtype,
                                            min_depth, vaf_ambig, vaf_call)
    ambig_pos = {seq_name: len(pos) for seq_name, pos in ambig_pos.items()}
    report['ambig_pos'] = report['seq_name'].map(ambig_pos)
    get_perc = lambda row: row['ambig_pos'] * 100 / row['seq_length']
    report['ambig_perc'] = report.apply(get_perc, axis=1)
    report['ambig_perc'] = round(report['ambig_perc'], 1)
    ''' Count variant positions. '''
    variant_pos = {}
    for seq_name in report['seq_name'].unique():
        species, segment, subtype = seq_name.split('|')[1:4]
        variant_pos[seq_name] = get_variant_pos(out_name, species, segment,
                                                subtype, min_depth, vaf_ambig,
                                                vaf_call)
    variant_pos = {seq_name: len(pos) for seq_name, pos in variant_pos.items()}
    report['variant_pos'] = report['seq_name'].map(variant_pos)
    get_perc = lambda row: row['variant_pos'] * 100 / row['seq_length']
    report['variant_perc'] = report.apply(get_perc, axis=1)
    report['variant_perc'] = round(report['variant_perc'], 1)
    ''' Get ref seq used for mapping reference. '''
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    ref_seqs = read_fasta(mapping_refs)
    get_seq_name = lambda header: '|'.join(header.split('|')[:4]) + '|'
    get_ref_seq = lambda header: '|'.join(header.split('|')[-3:-1])
    ref_seqs = {get_seq_name(header): get_ref_seq(header)
                for header in ref_seqs.keys()}
    report['ref_seq'] = report['seq_name'].map(ref_seqs)
    ''' Write out report. '''
    cols = ['seq_name', 'seq_length', 'reads_mapped', 'scaffold_completeness',
            'consensus_completeness', 'low_cov_perc', 'ambig_perc',
            'variant_perc', 'ref_seq']
    report = report[cols].drop_duplicates()
    report_path = os.path.join(out_name, f'{out_name}_report.tsv')
    report.to_csv(report_path, index=False, sep='\t')


def make_plot(out_name: str, min_depth: int, vaf_ambig: float, vaf_call: float):
    ''' Generate a plot with depth of coverage and low coverage/ambiguous/variant
    positions for each species/segment/subtype in the library. '''
    print('Making plot...')
    ''' Initialize plot. '''
    sb.set_style('whitegrid')
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    mapping_refs = read_fasta(mapping_refs)
    fig_size = (8.5, 3 * len(mapping_refs.keys()))
    fig, axs = plt.subplots(len(mapping_refs.keys()), sharex=True,
                            figsize=fig_size)
    ''' Load depth of cov data. '''
    depth_of_cov = get_depth_of_cov(out_name)
    max_position = depth_of_cov['position'].max()
    max_depth = depth_of_cov['depth'].max()
    min_depth = depth_of_cov[depth_of_cov['depth']>0]['depth'].min()
    if min_depth > 0 and max_depth - min_depth > 500:
        log_scale = True
    else:
        log_scale = False
    for seq_name, ax in zip(mapping_refs.keys(), axs):
        species, segment, subtype = seq_name.split('|')[1:4]
        ax_data = depth_of_cov[depth_of_cov['seq_name']==seq_name]
        sb.lineplot(x='position', y='depth', hue='tool', data=ax_data, ax=ax,
                    legend=False, hue_order=('bedtools', 'freebayes'),
                    palette=('grey', 'black'), lw=1)
        ax.axvspan(ax_data['position'].max(), max_position, color='grey')
        ''' Format x axis. '''
        ax.set_xlabel('Position')
        ax.set_xlim(1, max_position)
        xticks = [1] + [200 * i for i in range(1, int(max_position / 200) + 1)]
        ax.set_xticks(xticks)
        ''' Format y axis. '''
        label = 'Depth of coverage'
        if log_scale:
            ax.set_yscale('log')
            label += ' (log)'
            ymin = 1
            yticks = [10 ** i for i in range(int(log10(max_depth)) + 1)]
        else:
            ymin = 0
            yticks = [20 * i for i in range(int(max_depth / 20) + 1)]
        ax.set_ylabel(label)
        ax.set_ylim(ymin, max_depth)
        ax.set_yticks(yticks)
        ''' Add spans of low cov positions. '''
        ax.axhline(min_depth, xmin=0, xmax=ax_data['position'].max() / max_position,
                   color='red', alpha = 0.5, lw=1)
        low_cov_pos = get_low_cov_pos(out_name, species, segment, subtype, min_depth)
        for start, end in get_spans(low_cov_pos):
            ax.axvspan(start, end, color='red', alpha=0.2)
        ''' Add spans of ambiguous positions. '''
        ambig_pos = get_ambig_pos(out_name, species, segment, subtype, min_depth,
                                  vaf_ambig, vaf_call)
        for start, end in get_spans(ambig_pos):
            ax.axvspan(start, end, color='blue', alpha=0.2)
        ''' Add spans of variant positions. '''
        variant_pos = get_variant_pos(out_name, species, segment, subtype, min_depth,
                                      vaf_ambig, vaf_call)
        for start, end in get_spans(variant_pos):
            ax.axvspan(start, end, color='green', alpha=0.2)
        ''' Add title. '''
        ax.set_title(f'{species}/{segment} ({subtype})')
    ''' Save fig. '''
    plt.tight_layout()
    plt.savefig(os.path.join(out_name, f'{out_name}_plot.png'), dpi=400)


def get_depth_of_cov(out_name: str):
    ''' Generate DataFrame with depth of coverage at each position in each
    sequence. Generate depth of coverage from combined alignment using bedtools
    and from each individual species/segment/subtype by parsing FreeBayes
    output.'''
    print(' Getting depth of coverage...')
    ''' Get depth of cov for combined mapping using bedtools. '''
    full_alignment = os.path.join(out_name, f'{out_name}_alignment.bam')
    depth_of_cov_bedtools = os.path.join(out_name, f'depth_of_cov_bedtools.tsv')
    command = (f'bedtools genomecov -d -ibam {full_alignment} '
               f'> {depth_of_cov_bedtools}')
    process_name = 'bedtools_genomecov'
    error_code = 42
    run(out_name, command, process_name, error_code)
    cols = ['seq_name', 'position', 'depth']
    depth_of_cov_bedtools = pd.read_csv(depth_of_cov_bedtools, sep='\t',
                                        names=cols)
    depth_of_cov_bedtools['tool'] = 'bedtools'
    ''' Get depth of cov for each species/segment/subtype by parsing FreeBayes
    output. '''
    mapping_refs = os.path.join(out_name, f'{out_name}_mapping_refs.fa')
    mapping_refs = read_fasta(mapping_refs)
    seq_names = []
    positions = []
    depths = []
    for seq_name in mapping_refs.keys():
        species, segment, subtype = seq_name.split('|')[1:4]
        for values in parse_freebayes_pileup(out_name, species, segment,
                                             subtype):
            position, depth, max_alt_depth, total_alt_depth = values
            seq_names.append(seq_name)
            positions.append(position)
            depths.append(depth)
    depth_of_cov_freebayes = pd.DataFrame()
    depth_of_cov_freebayes['seq_name'] = seq_names
    depth_of_cov_freebayes['position'] = positions
    depth_of_cov_freebayes['depth'] = depths
    depth_of_cov_freebayes['tool'] = 'freebayes'
    ''' Concat results. '''
    depth_of_cov = pd.concat([depth_of_cov_bedtools, depth_of_cov_freebayes],
                         sort=True, ignore_index=True)
    return depth_of_cov


def collect_garbage(out_name: str):
    ''' Remove space-consuming intermediate files. '''
    shutil.rmtree(os.path.join(out_name, 'spades_output'))
    files_to_keep = ['alignment.bam', 'alignment.bam.bai', 'consensus_seqs.fa',
                     'report.tsv', 'plot.png']
    files_to_keep = [f'{out_name}_{file}' for file in files_to_keep]
    for file in os.listdir(out_name):
        if os.path.isfile(os.path.join(out_name, file)):
            if file not in files_to_keep:
                os.remove(os.path.join(out_name, file))


if __name__ == '__main__':
    main()
