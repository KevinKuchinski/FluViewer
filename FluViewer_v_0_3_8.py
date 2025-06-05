import sys as sys
import os as os
import shutil as shutil
import subprocess as sp
import numpy as np
import pandas as pd


from collections import Counter


def main():
    version = '0.3.8'
    args = parse_args(sys.argv, version)
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer/\n')
    print('Job name:', args['-n'])
    print('Output directory:', args['-o'])
    print('Fwd reads:', args['-f'])
    print('Rev reads:', args['-r'])
    print('Reference sequence DB:', args['-d'])
    print('CPUs used:', args['-T'])
    check_input_files(args['-f'], args['-r'], args['-d'])
    make_output_dir(args['-o'])
    length_ranges = check_database(args['-d'])
    assemble_contigs(args['-o'], args['-n'], args['-f'], args['-r'], args['-T'])

    blast_results = align_contigs_to_ref_seqs(args['-o'], args['-n'], args['-d'], args['-T'],
                                              args['-i'], args['-l'])

    blast_results = annotate_contig_alignments(blast_results)
    blast_results = mixed_infection_check(blast_results, args['-m'])
    make_scaffolds(args['-o'], args['-n'], blast_results)
    blast_results = align_scaffolds_to_ref_seqs(args['-o'], args['-n'], args['-d'], args['-T'])
    blast_results = filter_scaffold_alignments(blast_results)
    make_mapping_refs(args['-o'], args['-n'], blast_results, args['-d'])
    map_reads_to_mapping_refs(args['-o'], args['-n'], args['-f'], args['-r'], args['-q'])
    hi_qual_allele_calls = get_hi_qual_allele_calls(args['-o'], args['-n'], args['-D'], args['-q'])
    consensus_allele_calls = get_consensus_allele_calls(hi_qual_allele_calls, args['-D'], args['-c'],
                                                        args['-a'], args['-A'])
    make_consensus_seqs(args['-o'], args['-n'], consensus_allele_calls, length_ranges, args['-L'])
    write_mpileup_report(args['-o'], args['-n'], consensus_allele_calls, hi_qual_allele_calls)
    write_ambig_report(args['-o'], args['-n'], consensus_allele_calls, hi_qual_allele_calls)
    write_consensus_seq_report(args['-o'], args['-n'], args['-D'])
    if args['-g']:
        collect_garbage(args['-o'], args['-n'])
    print('\nDone.\n')


def print_usage(version: str):
    ''' Display usage message and condensed manual for FluViewer. '''
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer\n')
    print('Usage: FluViewer -n <output_name> -o <output_dir> -f <path_to_fwd_reads> '
          '-r <path_to_rev_reads> -d <path_to_db_file> [ <optional_args> ]\n')
    print('Required arguments:')
    print(' -n : output name')
    print(' -o : output directory')
    print(' -f : path to FASTQ file containing forward reads')
    print(' -r : path to FASTQ file containing reverse reads')
    print(' -d : path to FASTA file containing FluViewer database')
    print('Optional arguments:')
    print(' -i : Minimum sequence identity between database reference sequences and contigs '
          '(percentage, default = 90, min = 0, max = 100)')
    print(' -l : Minimum length of alignment between database reference sequences and contigs '
          '(int, default = 100, min = 32)')
    print(' -D : minimum read depth for base calling (int, default = 20, min = 1)')
    print(' -q : Minimum PHRED score for mapping quality and base quality during variant calling '
          '(int, default = 20, min = 0)')
    print(' -c : Consensus allele fraction (float, default = 0.75, min = 0, max = 1)')
    print(' -a : Alternate allele fraction (float, default = 0.05, min = 0, max = 1)')
    print(' -L : Length range tolerance for consensus sequences '
          '(percentage, default=1, min=0, max=100)')
    print(' -T : Threads used for computationally intensive steps '
          '(int, default = all available CPUs, min = 1)')
    print('Optional flags:')
    print(' -m : allow analysis of mixed infections (no value provided)')
    print(' -A : Mask ambiguous indels (no value provided)')
    print(' -g : Disable garbage collection and retain intermediate analysis files '
          '(no value provided)')
    print()


def parse_args(args: list, version: str):
    ''' Parse command line input. Ensure that required arguments have been provided with values.
    Ensure that provided values are of the correct type and fall within accepted parameters.
    Check that each argument has been set only once. '''
    arg_names = set(arg for arg in args if arg[0] == '-')
    multi_set_args = set()
    for arg in arg_names:
        if args.count(arg) > 1:
            multi_set_args.add(arg)
    if multi_set_args != set():
        multi_set_args = ', '.join(sorted(multi_set_args))
        print('\nERROR: The following arguments have been set more than once: '
              '{multi_set_args}\nAnalysis aborted.\n')
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
    ''' Check if all required arguments were provided. '''
    required_args = {'-n', '-o', '-f', '-r', '-d'}
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys():
            missing_args.add(required_arg)
    if missing_args != set():
        missing_args = ', '.join(sorted(missing_args))
        print(f'\nERROR: Values must be provided for the argument following '
              f'arguments: {missing_args}\nAnalysis aborted.\n')
        print_usage(version)
        exit(1)    
    ''' For arguments that are set without values, set to boolean values. '''
    default_arg_values = {'-m': False, '-A': False, '-g': True}
    bad_args = set()
    for arg, value in default_arg_values.items():
        if arg not in arg_values:
            arg_values[arg] = value
        elif arg in arg_values and arg_values[arg] == '':
            arg_values[arg] = not value
        elif arg in arg_values and arg_values[arg] != '':
            bad_args.add(arg)
    if bad_args != set():
        bad_args = ', '.join(sorted(bad_args))
        print(f'\nERROR: The following arguments should be provided without values: {bad_args}'
              '\nAnalysis aborted.\n')
        exit(1)
    ''' Check that all arguments were provided with values. '''
    bad_args = set()
    for arg, arg_value in arg_values.items():
        if arg_value == '':
            bad_args.add(arg)
    if bad_args != set():
        bad_args = ', '.join(sorted(bad_flag_args))
        print(f'\nERROR: The following arguments must be provided with values: {bad_args}'
              '\nAnalysis aborted.\n')
        exit(1)    
    ''' Check if unrecognized arguments were provided. '''
    arg_value_types = {'-n': str, '-o': str, '-f': str, '-r': str, '-d': str, '-i': float, '-l': int,
                       '-m': bool, '-D': int, '-q': int, '-c': float, '-a': float, '-A': bool,
                       '-L': float, '-T': int, '-g': bool}
    recognized_args = set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args.add(provided_arg)
    if unrecognized_args != set():
        unrecognized_args = ', '.join(sorted(unrecognized_args))
        print(f'\nERROR: The following provided arguments are not recognized: {unrecognized_args}'
              '\nAnalysis aborted.\n')
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
                print(f'\nERROR: Value for argument {arg} must be of type {value_type}.'
                      '\nAnalysis aborted.\n')
                print_usage(version)
                exit(1)
    ''' Check if provided values are within the correct range. '''
    min_arg_values = {'-i': 0, '-l': 32, '-D': 1, '-q': 0, '-c': 0, '-a': 0, '-L': 0, '-T': 1}
    max_arg_values = {'-i': 100, '-c': 1, '-a': 1, '-L': 100}
    for arg, value in arg_values.items():
        if arg in min_arg_values.keys() and value < min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must not be less than {min_arg_values[arg]}'
                  '\nAnalysis aborted.\n')
            print_usage(version)
            exit(1)
        if arg in max_arg_values.keys() and value > max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must not exceed {max_arg_values[arg]}'
                  '\nAnalysis aborted.\n')
            print_usage(version)
            exit(1) 
    ''' Check if alternate allele fraction is lower than consensus allele fraction. '''
    if '-c' in arg_values and '-a' in arg_values:
        if arg_values['-a'] >= arg_values['-c']:
            print('\nERROR: Alternate allele fraction (-v) must be lower than consensus allele '
                  'fraction (-c).\nAnalysis aborted.\n')
            print_usage(version)
            exit(1)
    ''' Assign default values to unspecified arguments. '''
    default_arg_values = {'-i': 90, '-l': 100, '-m': False, '-D': 20, '-q': 20, '-c': 0.75,
                          '-a': 0.05, '-A': False, '-L': 1, '-T': 1, '-g': True}
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value
    ''' Check provided argument values for reserved characters. '''
    reserved_chars = '*;\&|()<>'
    if any(char in str(value) for char in reserved_chars
           for value in arg_values.values()):
        print('\nERROR: Argument values cannot contain any of the following reserved characters: '
              '*;\&|()<>')
        exit(1)
    ''' Return keyword args and their values. '''
    return arg_values


def check_input_files(fwd_reads: os.path, rev_reads: os.path, db: os.path):
    ''' Check that all provided input files exist. '''
    for file in [fwd_reads, rev_reads, db]:
        if not os.path.isfile(file):
            print(f'\nERROR: Input file {file} does not exist.\nAnalysis aborted.\n')
            exit(2)


def make_output_dir(out_dir: str):
    ''' Create a directory for output after checking that it does not already exist. '''
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
        print(f'\nWARNING: Output directory {out_dir} already exists! '
              'Files will be overwritten.\nAnalysis continuing...\n')
    logs_path = os.path.join(out_dir, 'logs')
    if not os.path.exists(logs_path):
        os.mkdir(logs_path)
    else:
        print(f'\nWARNING: Logs directory {logs_path} already exists! '
              'Files will be overwritten.\nAnalysis continuing...\n')


def check_database(db: os.path):
    ''' Checks the contents of the provided reference sequence database to ensure proper header
    formatting and unambiguous sequences. Returns a Pandas DataFrame with the min and max lengths
    for each species/segment/subtype. '''
    print('Checking reference sequence database...')
    ''' Make sure each database entry has a header with the expected number of fields. '''
    seqs = read_fasta(db)
    bad_headers = set()
    for header in seqs.keys():
        if len(header.split('|')) < 5:
               bad_headers.add(header)
    if bad_headers != set():
        print(f'\nERROR: The following headers do not contain the minimum number of |-delimited fields:')
        for header in bad_header:
            print(f' {header}')
        print('\nAnalysis aborted.\n')
        exit(3)
    ''' Parse database entries and get sequence lengths. '''
    species = []
    segments = []
    subtypes = []
    lengths = []
    for header, seq in seqs.items():
        seq_species, segment, subtype = header.split('|')[2:5] 
        species.append(seq_species)
        segments.append(segment)
        subtypes.append(subtype)
        lengths.append(len(seq))
    ''' Create DataFrame with length ranges. '''
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
    ''' A generalized function for reading the contents of a FASTA file into a dict of
    header: sequence. '''
    if not os.path.exists(file_path):
        print('\nERROR: File {file_path} does not exist!\nAnalysis aborted.\n')
        exit(4)
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
    ''' Check that each header appears only once. '''
    headers = Counter(headers)
    bad_headers = set()
    for header in headers:
        if headers[header] > 1:
            bad_headers.add(header)
    if bad_headers != set():
        print('\nERROR: The following headers appear more than once in the provided database:')
        for header in bad_headers:
            print(f' {header}')
        print('\nAnalysis aborted.\n')
        exit(5)
    if seqs[''] != '':
        print(f'\nERROR: FASTA file {file_path} Contains sequence before first header!'
              '\nAnalysis aborted.\n')
        exit(6)
    seqs = {header: seq for header, seq in seqs.items() if header != ''}
    return seqs


def run(out_dir: str, out_name: str, command: str, process_name: str, error_code: int):
    ''' A generalized function for running subprocesses, logging their output, and trapping
    erroneous exit statuses. '''
    stdout_file = os.path.join(out_dir, 'logs', f'{process_name}_stdout.txt')
    stderr_file = os.path.join(out_dir, 'logs', f'{process_name}_stderr.txt')
    stdout_file = open(stdout_file, 'w')
    stderr_file = open(stderr_file, 'w')
    complete_process = sp.run(command, stdout=stdout_file, stderr=stderr_file, shell=True)
    stdout_file.close()
    stderr_file.close()
    return_code = complete_process.returncode
    if return_code != 0:
        print(f'\nERROR: Subprocess {process_name} failed (Exit status: {return_code})'
              '\nAnalysis aborted.\n')
        exit(error_code)


def assemble_contigs(out_dir: str, out_name: str, fwd_reads: os.path, rev_reads: os.path,
                     threads: int):
    ''' Reads are assembled de novo into contigs using SPAdes. '''
    print('Assembling reads into contigs...')
    spades_output = os.path.join(out_dir, 'spades_output')
    command = (f'spades.py --rnaviral --threads {threads} -1 {fwd_reads} -2 {rev_reads} '
               f'-o {spades_output}')
    process_name = 'spades'
    error_code = 7
    run(out_dir, out_name, command, process_name, error_code)
    if not os.path.isfile(os.path.join(spades_output, 'contigs.fasta')):
        print('\nERROR: No contigs assembled!\nAnalysis aborted.\n')
        exit(8)


def align_contigs_to_ref_seqs(out_dir: str, out_name: str, db: os.path, num_threads: int, identity: float,
                              length: int):
    ''' Contigs are aligned to reference sequences using BLASTn. Initial alignment filtering based on
    identity and length of alignment. '''
    print('Aligning contigs to reference sequences...')
    if any([os.path.exists(db + '.' + suffix) == False
            for suffix in ['nhr', 'nin' , 'nsq']]):
        command = f'makeblastdb -in {db} -dbtype nucl'
        process_name = 'makeblastdb_contigs'
        error_code = 9
        run(out_dir, out_name, command, process_name, error_code)
    blast_output = os.path.join(out_dir, 'contigs_blast.tsv')
    contigs = os.path.join(out_dir, 'spades_output', 'contigs.fasta')
    cols = 'qseqid sseqid pident length bitscore sstart send qseq sseq slen'
    command = (f'blastn -query {contigs} -db {db} -num_threads {num_threads} -outfmt "6 {cols}" '
               f'> {blast_output}')
    process_name = 'blastn_contigs'
    error_code = 10
    run(out_dir, out_name, command, process_name, error_code)
    blast_results = pd.read_csv(blast_output, names=cols.split(' '), sep='\t')
    blast_results = blast_results[blast_results['pident']>=identity]
    blast_results = blast_results[blast_results['length']>=length]
    if len(blast_results) == 0:
        print(f'\nERROR: No contigs aligned to reference sequences!\nAnalysis aborted.\n')
        exit(11)
    return blast_results


def annotate_contig_alignments(blast_results: pd.DataFrame):
    print('Filtering contig alignments...')
    cols = ['qseqid', 'sseqid', 'bitscore']
    group_cols = ['qseqid', 'sseqid']
    combined_bitscores = blast_results[cols].groupby(group_cols).sum().reset_index()
    cols = ['qseqid', 'bitscore']
    max_combined_bitscores = combined_bitscores[cols].groupby('qseqid').max().reset_index()
    best_matches = pd.merge(combined_bitscores, max_combined_bitscores, on=['qseqid', 'bitscore'])
    best_matches = best_matches[['qseqid', 'sseqid']].drop_duplicates()
    blast_results = pd.merge(blast_results, best_matches, on=['qseqid', 'sseqid'])
    print('Annotating contig alignments...')
    ''' Alignments of contigs to reference sequences are annotated with species, segment, and subtype. '''
    annots = blast_results[['sseqid']].drop_duplicates()
    annots['species'] = annots.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
    annots['segment'] = annots.apply(lambda row: row['sseqid'].split('|')[3], axis=1)
    annots['subtype'] = annots.apply(lambda row: row['sseqid'].split('|')[4], axis=1)
    blast_results = pd.merge(blast_results, annots, on='sseqid')
    return blast_results


def mixed_infection_check(blast_results: pd.DataFrame, allow_mixtures: bool):
    ''' Filtered alignments of contigs to reference sequences are analyzed to determine if multiple
    species or subtypes are present. If analysis of mixed infections is allowed (-m),
    scaffolds, mapping references, and consensus sequencess are only generated for species/segments
    with subtype other than "none". If analysis of mixed infections is not allowed, an error is raised
    and analysis aborts. '''
    print('Checking for mixed infection...')
    ''' Check for mixture of species. '''
    species = blast_results['species'].unique()
    num_species = len(species)
    if num_species > 1:
        species = ', '.join(sorted(species))
        if not allow_mixtures:
            print(f'\nERROR: Multiple species detected ({species})!\nAnalysis aborted.\n')
            exit(12)
        print(f'\nWARNING: Multiple species detected ({species})!\nAnalysis continuing...')
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
            print('\nERROR: Multiple subtypes detected for the following species/segments:')
        else:
            print('\nWARNING: Multiple subtypes detected for the following species/segments:')            
        for species in sorted(mixed_segments['species'].unique()):
            species_results = mixed_segments[mixed_segments['species']==species]
            for segment in species_results['segment'].unique():
                segment_results = species_results[species_results['segment']==segment]
                subtypes = ', '.join(sorted(segment_results['subtype'].unique()))
                print(f' {species}/{segment}: ({subtypes})')
        if not allow_mixtures:
            print('\nAnalysis aborted.\n')
            exit(13)
        else:
            print('\nAnalysis continuing...')
    ''' Discard alignments for species/segments without defined subtype if mixed
    infection.  '''
    if num_species > 1 or len(mixed_segments) > 0:
        blast_results = blast_results[blast_results['subtype']!='none']
    if len(blast_results) == 0:
        print(f'\nERROR: No contigs aligned to reference sequences!\nAnalysis aborted.\n')
        exit(11)
    return blast_results


def make_scaffolds(out_dir: str, out_name: str, blast_results: pd.DataFrame):
    ''' Generate a scaffold sequence for each species/segment/subtype by collapsing all contigs
    describing that species/segment/subtype into a single sequence. '''
    print('Making scaffolds...')
    scaffolds = os.path.join(out_dir, f'{out_name}_scaffolds.fa')
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
            scaffold_seq = make_scaffold_seq(out_dir, out_name, species, segment, subtype, scaffold_data)
            header = f'>{out_name}|{species}|{segment}|{subtype}|'
            output_file.write(header + '\n')
            output_file.write(scaffold_seq + '\n')


def pick_scaffolding_ref(scaffold_data: pd.DataFrame):
    ''' Pick a reference sequence for a species/segment/subtype to use for scaffolding. Contig
    alignments against this sequence are used to roughly position the contig inside the segment
    before conducting a multiple sequence alignment on the contigs. The scaffolding reference is
    chosen by selecting the reference seq with the most positions covered by contigs. If there are
    multiple options remaining, the longest is chosen, followed by the first alphabetically. '''
    print('  Picking reference sequence for scaffolding...')
    ''' Find ref seq(s) most covered by contigs. '''
    sseq_cov = {sseq: set() for sseq in scaffold_data['sseqid'].unique()}
    cols = ['sseqid', 'sstart', 'send']
    for index, row in scaffold_data[cols].drop_duplicates().iterrows():
        start, stop = sorted((row['sstart'], row['send']))
        for position in range(start, stop + 1):
            sseq_cov[row['sseqid']].add(position)
    sseq_cov = {sseq: len(sseq_cov[sseq]) for sseq in sseq_cov.keys()}
    scaffold_data['sseq_cov'] = scaffold_data.apply(lambda row: sseq_cov[row['sseqid']], axis=1)
    max_cov = scaffold_data['sseq_cov'].max()
    scaffold_data = scaffold_data[scaffold_data['sseq_cov']==max_cov]
    ''' Take longest remaining ref seq(s) for each segment/subtype. '''
    max_slen = scaffold_data['slen'].max()
    scaffold_data = scaffold_data[scaffold_data['slen']==max_slen]
    ''' Take first alphabetical remaining ref seq for each segment/subtype. '''
    first_alpha = scaffold_data['sseqid'].min()
    scaffold_data = scaffold_data[scaffold_data['sseqid']==first_alpha]
    return scaffold_data


def flip_reversed_alignments(blast_results: pd.DataFrame):
    ''' Helper function that takes a DataFrame of blastn results and replaces qseq and sseq strings
    with their reverse complements if the alignment is reversed. It also creates columns for start
    and end coordinates that are essentially sstart and send with the smaller of these values as
    the start. '''
    ''' Get alignment start and end coordinates. '''
    blast_results['start'] = blast_results.apply(lambda row: min((row['sstart'], row['send'])), axis=1)
    blast_results['end'] = blast_results.apply(lambda row: max((row['sstart'], row['send'])), axis=1)
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
    blast_results['qseq'] = blast_results.apply(flip_qseq, axis=1)
    def flip_sseq(row):
        if row['sstart'] > row['send']:
            return rev_comp_seq(row['sseq'])
        else:
            return row['sseq']
    blast_results['sseq'] = blast_results.apply(flip_sseq, axis=1)
    return blast_results


def make_scaffold_seq(out_dir: str, out_name: str, species: str, segment: str, subtype: str,
                      scaffold_results: pd.DataFrame):
    ''' Generate the scaffold sequence for a genome segment. Contig alignments against the chosen
    reference seq are used to roughly position each contig within the segment (a number of terminal
    Ns are added to each end of each contig based on the start and end coordinates of the contig's
    alignment). At this time, any part of a contig extending past the scaffolding reference is
    trimmed. A mutliple sequence alignment is then performed on the trimmed, positioned contigs.
    The scaffold sequence is generated from the consensus of that multiple sequence alignment.
    Ambiguous positions (the dominant base does not account for at least 50% of bases) are replaced
    with Ns. Positions not covered by any contigs are also replaced with Ns in the scaffold. '''
    print('  Trimming and positioning contigs...')
    ''' Make sure all alignments are in the forward orientation. '''
    scaffold_results = flip_reversed_alignments(scaffold_results)
    ''' Write all contigs to a FASTA file. Trim contigs based on their alignments to reference
    sequences. Also add leading and trailing Ns to contig so that it is properly positioned within
    the reference sequence. '''
    contigs = os.path.join(out_dir, f'{species}_{segment}_{subtype}_contigs.fa')
    with open(contigs, 'w') as contig_file:
        contig_counter = 1
        for index, row in scaffold_results.iterrows():
            header = f'>{species}_{segment}_{subtype}_contig_{contig_counter}'
            contig_file.write(header + '\n')
            seq = 'N' * (row['start'] - 1)
            seq += row['qseq'].replace('-', '')
            seq += ('N' * (row['slen'] - row['end']))
            contig_file.write(seq + '\n')
            contig_counter += 1
    ''' Generate multiple sequence alignments of trimmed/positioned contigs. '''
    print('  Collapsing contigs...')
    aligned_contigs = os.path.join(out_dir, f'{species}_{segment}_{subtype}_contigs.afa')
    if contig_counter > 2:
        command = (f'clustalw -INFILE={contigs} -OUTFILE={aligned_contigs} -OUTPUT=FASTA')
        process_name = f'clustalw_{segment}'
        error_code = 14
        run(out_dir, out_name, command, process_name, error_code)
    else:
        shutil.copyfile(contigs, aligned_contigs)
    ''' Load aligned contigs and replace leading and trailing Ns and gaps with dots so that they
    are ignored when determining consensus bases. '''
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
        exit(15)
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
    scaffold_seq = scaffold_seq.replace('-', '')
    return scaffold_seq


def align_scaffolds_to_ref_seqs(out_dir: str, out_name: str, db: os.path, num_threads: int):
    ''' Scaffolds are aligned to reference sequences using BLASTn. '''
    print('Aligning scaffolds to reference sequences...')
    if any([os.path.exists(db + '.' + suffix) == False
            for suffix in ['nhr', 'nin' , 'nsq']]):
        command = f'makeblastdb -in {db} -dbtype nucl'
        process_name = 'makeblastdb_scaffolds'
        error_code = 16
        run(out_dir, out_name, command, process_name, error_code)
    with open(db) as input_file:
        num_db_seqs = sum(line[0] == '>' for line in input_file)
    blast_output = os.path.join(out_dir, 'scaffolds_blast.tsv')
    scaffolds = os.path.join(out_dir, f'{out_name}_scaffolds.fa')
    cols = 'qseqid sseqid bitscore sstart send qseq sseq'
    command = (f'blastn -query {scaffolds} -db {db} -num_threads {num_threads} '
               f'-max_target_seqs {num_db_seqs} -outfmt "6 {cols}" > {blast_output}')
    process_name = 'blastn_scaffolds'
    error_code = 17
    run(out_dir, out_name, command, process_name, error_code)
    blast_results = pd.read_csv(blast_output, names=cols.split(' '), sep='\t')
    if len(blast_results) == 0:
        print(f'\nERROR: No scaffolds aligned to reference sequences! \nAnalysis aborted.\n')
        exit(18)
    return blast_results


def filter_scaffold_alignments(blast_results: pd.DataFrame):
    ''' For each scaffold, choose a reference sequence that will be used to fill in Ns. This scaffold
    with filled-in Ns is the mapping reference. The chosen reference sequences are each scaffold's
    best match (based on summed bitscores). If a scaffold has multiple best matches based on summed
    bitscores, choose the first alphabetically. '''
    print('Choosing reference sequences used for infilling scaffolds...')
    ''' Make sure all alignments are in the forward orientation. '''
    blast_results = flip_reversed_alignments(blast_results)
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
    ''' Get first best-matching ref seq (alphabetically) in case of ties. '''
    cols = ['qseqid', 'sseqid']
    group_cols = ['qseqid']
    first_alpha = blast_results[cols].groupby(group_cols).min().reset_index()
    blast_results = pd.merge(blast_results, first_alpha, on=['qseqid', 'sseqid'])
    return blast_results


def make_mapping_ref(blast_results: pd.DataFrame, species: str, segment: str, subtype: str, ref_seq: str):
    ''' This function generates a mapping reference for single genome segment. '''
    blast_results = blast_results[blast_results['species']==species]
    blast_results = blast_results[blast_results['segment']==segment]
    blast_results = blast_results[blast_results['subtype']==subtype]
    blast_results = blast_results.sort_values(by=['start', 'end'], ascending=[True, False])
    last_position = 0
    seq = ''
    for index, row in blast_results.iterrows():
        if row['start'] > last_position:
            seq += ref_seq[last_position:row['start'] - 1]
        if row['end'] > last_position:
            qseq = row['qseq'].upper()
            sseq = row['sseq'].upper()
            if row['start'] <= last_position:
                start = (last_position - row['start']) + 1
                qseq = qseq[start:]
                sseq = sseq[start:]
            for qbase, sbase in zip(qseq, sseq):
                if qbase in 'ATGC':
                    seq += qbase
                else:
                    seq += sbase
            last_position = row['end']
    seq += ref_seq[last_position:].upper()
    seq = seq.replace('-', '')
    return seq


def make_mapping_refs(out_dir: str, out_name: str, blast_results: pd.DataFrame, db: os.path):
    ''' Mapping references are created for each species/segment/subtype detected. These consist of
    the scaffold with all Ns in the scaffold filled-in using the corresponding positions from that
    scaffold's best-matching reference sequence. '''
    print('Creating mapping references...')
    ''' Create dict with best-matching ref seq for each segment. '''
    best_ref_seqs = {seq_name: '' for seq_name in blast_results['sseqid'].unique()}
    with open(db, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
            elif header in best_ref_seqs:
                best_ref_seqs[header] += line.strip()
    ''' Annotate blast results with species, segment, and subtype '''
    annots = blast_results[['sseqid']].drop_duplicates()
    annots['species'] = annots.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
    annots['segment'] = annots.apply(lambda row: row['sseqid'].split('|')[3], axis=1)
    annots['subtype'] = annots.apply(lambda row: row['sseqid'].split('|')[4], axis=1)
    blast_results = pd.merge(blast_results, annots, on='sseqid')
    cols = ['sseqid', 'species', 'segment', 'subtype', 'start', 'end', 'qseq', 'sseq']
    blast_results = blast_results[cols].drop_duplicates()
    ''' Create mapping ref for each species/segment/subtype. '''
    mapping_refs = {}
    cols = ['species', 'segment', 'subtype', 'sseqid']
    for index, row in blast_results[cols].drop_duplicates().iterrows():
        species, segment, subtype = row['species'], row['segment'], row['subtype']
        print(f' {species}/{segment} ({subtype})')
        ref_seq_name, ref_seq = row['sseqid'], best_ref_seqs[row['sseqid']]
        mapping_ref = make_mapping_ref(blast_results, species, segment, subtype, ref_seq)
        mapping_refs[(species, segment, subtype, ref_seq_name)] = mapping_ref
    ''' Write out mapping references in FASTA. '''
    mapping_refs_path = os.path.join(out_dir, f'{out_name}_mapping_refs.fa')
    with open(mapping_refs_path, 'w') as output_file:
        for (species, segment, subtype, ref_seq_name), seq in mapping_refs.items():
            accession, ref_name, species, segment, subtype = ref_seq_name.lstrip('>').split('|')[:5]
            header = [out_name, species, segment, subtype, accession, ref_name]
            header = '>' + '|'.join(header) + '|'
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def map_reads_to_mapping_refs(out_dir: str, out_name: str, fwd_reads: os.path, rev_reads: os.path,
                              min_qual: float):
    ''' Input reads are mapped to the mapping references using BWA mem. The alignment is filtered to
    retain only paired reads in their primary alignments, then the alignment is sorted and indexed. '''
    print('Mapping reads to mapping references...')
    ''' Index combined mapping references. '''
    mapping_refs = os.path.join(out_dir, f'{out_name}_mapping_refs.fa')
    command = f'bwa index {mapping_refs}'
    process_name = 'bwa_index'
    error_code = 30
    run(out_dir, out_name, command, process_name, error_code)
    ''' Align input reads to combined mapping references. '''
    alignment = os.path.join(out_dir, 'alignment.sam')
    command = f'bwa mem {mapping_refs} {fwd_reads} {rev_reads} > {alignment}'
    process_name = 'bwa_mem'
    error_code = 31
    run(out_dir, out_name, command, process_name, error_code)
    ''' Filter, sort, and index combined alignment. '''
    filtered_alignment = os.path.join(out_dir, f'{out_name}_alignment.bam')
    command = (f'samtools view -f 1 -F 2316 -q {min_qual} -h {alignment} | samtools sort '
               f'-o {filtered_alignment}')
    process_name = 'samtools_view_combined'
    error_code = 32
    run(out_dir, out_name, command, process_name, error_code)
    command = f'samtools index {filtered_alignment}'
    process_name = 'samtools_index_combined'
    error_code = 33
    run(out_dir, out_name, command, process_name, error_code)


def parse_mpileup_line(ref_base: str, allele_calls: str, allele_quals: str, min_qual: int,
                       seq_name: str, position: int):
    ''' This function parses a single line from an mpileup file and returns the high quality allele
    calls at the position described by this line. Alleles are returned in a Counter object
    {allele: count}. Only three columns from the mpileup line are needed: the reference base at that
    position, the string of allele calls piled up at that position, and the string of PHRED base
    quality scores that correspond to the allele calls. '''
    i = 0 # Initialize index counter
    parsed_allele_calls = list() # Initialize list to hold parsed allele calls
    while i < len(allele_calls):
        if allele_calls[i] == '^':
            ''' The ^ symbol in the allele call string denotes that the following allele is at the
            start of its read. The symbol following ^ is the PHRED mapping quality of that read.
            This information is not used by this parser so these positions in the allele calls
            string are skipped past. Poorly mapped reads should be filtered out of the SAM/BAM
            before creating the mpileup. '''
            i += 2
        if allele_calls[i] in '.,':
            ''' The . and , symbols denote that a read contained the reference base in that position.
            Two different symbols are used to denote forward vs reverse orientated reads, but this
            information is not used by this parser. '''
            parsed_allele_calls.append(ref_base)
            i += 1
        elif allele_calls[i].upper() in 'ATGCN':
            ''' The A, T, G, and C symbols denote that a read contained the indicated single
            nucleotide mismatch in that position. Lower case letters are used to denote those
            mismatches on a reverse orientated read, but that information is not used by this
            parser. '''
            parsed_allele_calls.append(allele_calls[i].upper())
            i += 1
        elif allele_calls[i] in '*#':
            ''' The * and # symbols denote that a read did not contain that position because of an
            upstream deletion. Two different symbols are used to denote forward vs reverse
            orientated reads, but this information is not used by this parser. '''
            parsed_allele_calls.append('*')
            i += 1
        elif allele_calls[i] == '-':
            ''' The symbol - denotes the beginning of a deletion in that read. It will be followed
            by an integer indicating the length of the deletion, which will then be followed by the
            deleted nucleotides. NB that the symbol preceding the - symbol is considered part of the
            deletion: it indicates the nucleotide at the current position, which is the position
            preceding the deletion. For example, a read which contains the reference base at the
            current position, followed by a deletion of ATGC between the current position and the
            next position, will be encoded as .-4ATGC If the reference base in the current position
            was C, this deletion allele would be parsed as C-ATGC. '''
            i += 1
            deletion_length = ''
            while allele_calls[i].isnumeric():
                deletion_length += allele_calls[i]
                i += 1
            deletion_length = int(deletion_length)
            deletion = allele_calls[i:i+deletion_length].upper()
            deletion = parsed_allele_calls[-1] + '-' + deletion
            parsed_allele_calls = parsed_allele_calls[:-1]
            parsed_allele_calls.append(deletion)
            i += deletion_length
        elif allele_calls[i] == '+':
            ''' The symbol + denotes the beginning of an insertion in that read. It will be followed
            by an integer indicating the length of the insertion, which will then be followed by the
            inserted nucleotides. NB that the symbol preceding the + symbol is considered part of
            the insertion: it indicates the nucleotide at the current position, which is the position
            preceding the insertion. For example, a read which contains the reference base at the
            current position, followed by an insertion of ATGC between the current position and the
            next position, will be encoded as .+4ATGC If the reference base in the current position
            was C, this insertion allele would be parsed as C+ATGC. '''
            i += 1
            insertion_length = ''
            while allele_calls[i].isnumeric():
                insertion_length += allele_calls[i]
                i += 1
            insertion_length = int(insertion_length)
            insertion = allele_calls[i:i+insertion_length].upper()
            insertion = parsed_allele_calls[-1] + '+' + insertion
            parsed_allele_calls = parsed_allele_calls[:-1]
            parsed_allele_calls.append(insertion)
            i += insertion_length
        else:
            ''' Any other symbol in the allele call string is skipped past. '''
            i += 1
    ''' The length of the parsed allele calls list should match the length of the allele quality
    string. '''
    if len(parsed_allele_calls) != len(allele_quals):
        print(f'\nERROR: number of allele calls does not match number of quality scores at position '
              f'{position} in sequence {seq_name}!\nAnalysis aborted.\n')
        exit()
    ''' Discard allele calls if their PHRED base quality score is below the threshold set by
    min_qual. '''
    hi_qual_allele_calls = Counter()
    for call, qual in zip(parsed_allele_calls, allele_quals):
        if ord(qual) - 33 >= min_qual:
            hi_qual_allele_calls.update((call, ))
    return hi_qual_allele_calls


def get_hi_qual_allele_calls(out_dir: str, out_name: str, min_depth: int, min_qual: int):
    ''' This function creates an mpileup from a sorted BAM file using samtools. Each line of the
    mpileup is parsed to extract all high-quality allele calls at that line's position in mpileup.
    These high-quality allele calls are returned in a dict of {seq_name containing
    species/segment/subtype: a position-by-position list of high-quality allele calls (as Counter
    objects)}. '''
    print('Generating mpileup...')
    ''' Generate mpileup '''
    filtered_alignment = os.path.join(out_dir, f'{out_name}_alignment.bam')
    mapping_refs = os.path.join(out_dir, f'{out_name}_mapping_refs.fa')
    mpileup = os.path.join(out_dir, f'mpileup.tsv')
    max_depth = max((min_depth * 10, 200))
    command = (f'samtools mpileup -B -aa -d {max_depth} -f {mapping_refs} {filtered_alignment} '
               f'> {mpileup}')
    process_name = 'samtools_mpileup'
    error_code = 20
    run(out_dir, out_name, command, process_name, error_code)
    print('Parsing mpileup...')
    ''' Get all seq names '''
    seq_names = set()
    with open(mpileup, 'r') as input_file:
        for line in input_file:
            line = line.strip().split('\t')
            seq_name = line[0]
            seq_names.add(seq_name)
    ''' Get list of high quality allele calls at each position of each seq '''
    hi_qual_allele_calls = {seq_name: list() for seq_name in seq_names}
    with open(mpileup, 'r') as input_file:
        for line in input_file:
            line = line.strip().split('\t')
            seq_name = line[0]
            position = line[1]
            ref_base = line[2]
            allele_calls = line[4]
            allele_quals = line[5]
            calls = parse_mpileup_line(ref_base, allele_calls, allele_quals, min_qual, seq_name,
                                       position)
            hi_qual_allele_calls[seq_name].append(calls)
    return hi_qual_allele_calls


def get_degen_base(alt_bases: tuple):
    ''' Takes a tuple (sorted alphabetically) of alternate nucleotides and returns the corresponding
    single-letter degenerate nucleotide code. NB tuples of nucleotides are only passed to this
    function when there is definitely a degeneracy. Thus, tuples containing a single nucleotide must
    reflect situations where that nucleotide was an alternate allele among other symbols, eg
    deletion characters like * and #. This is why tuples of single nucleotides return the ambiguity
    character N. '''
    degen_bases = {tuple(): 'N',
                   ('A',): 'N',
                   ('C',): 'N',
                   ('G',): 'N',
                   ('T',): 'N',
                   ('A', 'T'): 'W',
                   ('C', 'G'): 'S',
                   ('G', 'T'): 'K',
                   ('A', 'C'): 'M',
                   ('C', 'T'): 'Y',
                   ('A', 'G'): 'R',
                   ('A', 'C', 'G'): 'V',
                   ('A', 'C', 'T'): 'H',
                   ('A', 'G', 'T'): 'D',
                   ('C', 'G', 'T'): 'B',
                   ('A', 'C', 'G', 'T'): 'N'}
    return degen_bases[alt_bases]


def get_consensus_allele_call(allele_calls: Counter, min_consensus_frac: float, min_alt_frac: float,
                              ambig_indels: bool):
    ''' This function takes a count of allele calls and returns the consensus allele call. The
    consensus allele will be the most abundant allele if it accounts for at least min_consensus_frac
    of the total number of allele calls. It must also be the only allele that is the most abundant.
    Otherwise, the consensus allele will be masked with degeneracy/ambiguity characters. The
    degeneracy/ambiguity characters are chosen to reflect all alleles accounting for at least
    min_alt_frac of the total number of allele calls. '''
    snps = Counter()
    indels = Counter()
    ''' Begin by splitting each allele call into two components: the single nucleotide at the current
    position and the indel between the current position and the next position. For alleles that have
    no indel component, the dummy indel '' is used. '''
    for call, count in allele_calls.items():
        if '-' in call:
            snp, indel = call.split('-')
            indel = '-' + indel
        elif '+' in call:
            snp, indel = call.split('+')
            indel = '+' + indel
        else:
            snp, indel = call, ''
        snps.update({snp: count})
        indels.update({indel: count})
    ''' Count the total number of alleles at the current position. '''
    total_depth = allele_calls.total()
    ''' Check if the most common snp satisfies the consensus criteria (there are no other snps that
    are equally abundant and it accounts for at least min_consensus_frac of the snps). '''
    most_common_snps = snps.most_common(1)
    most_common_snp_fraction = most_common_snps[0][1] / total_depth
    if (len(most_common_snps) == 1 and most_common_snp_fraction >= min_consensus_frac):
        final_snp = most_common_snps[0][0]
    else:
        ''' If the most common snp does not satisfy the consensus criteria, find the degeneracy
        character representing all nucleotides accounting for at least min_alt_frac of the snps. '''
        alt_snps = set()
        for snp, count in snps.items():
            alt_snp_frac = count / total_depth
            if snp != '*' and alt_snp_frac >= min_alt_frac:
                alt_snps.add(snp)
        final_snp = get_degen_base(tuple(sorted(alt_snps)))
    ''' Check if the most common indel satisfies the consensus criteria (there are no other indels
    that are equally abundant and it accounts for at least min_consensus_frac of the indels). '''
    most_common_indels = indels.most_common(1)
    most_common_indel_frac = most_common_indels[0][1] / total_depth
    if (len(most_common_indels) == 1 and most_common_indel_frac >= min_consensus_frac):
        final_indel = most_common_indels[0][0]
    elif ambig_indels:
        ''' If the most common indel does not satisfy the consensus criteria and ambiguous indel
        masking is active, create an ambiguous indel of the same  type and length as the most common
        alternative indel. If there are multiple alternative indels that are the most common and
        they are all insertions,  create an ambiguous insertion with the length of the longest one.
        If there  are multiple alternative indels that are the most common and they are all deletions
        or a mixture of insertions and deletions, create an ambiguous  deletion with the length of
        the longest deletion. '''
        alt_indels = Counter()
        for indel, count in indels.items():
            alt_indel_frac = count / total_depth
            if indel != '' and alt_indel_frac >= min_alt_frac:
                alt_indels.update((indel,))
        most_common_alt_indels = set()
        for indel in alt_indels.most_common(1):
            most_common_alt_indels.add(indel[0])
        if len(most_common_alt_indels) == 0:
            final_indel = ''
        elif len(most_common_alt_indels) == 1:
            final_indel = list(most_common_alt_indels)[0]
            final_indel = final_indel[0] +  ('N' * (len(final_indel[1:])))
        else:
            if all(indel[0] == '+' for indel in most_common_alt_indels):
                final_indel = max(most_common_alt_indels, key=lambda i: len(i[0]))
                final_indel = final_indel[0] +  ('N' * (len(final_indel[1:])))
            else:
                most_common_alt_indels = set(indel for indel in most_common_alt_indels
                                             if indel[0] == '-')
                final_indel = max(most_common_alt_indels, key=lambda i: len(i[0]))
                final_indel = final_indel[0] +  ('N' * (len(final_indel[1:])))
    else:
        ''' If the most common indel does not satisfy the consensus criteria but ambiguous indel
        masking is not active, add no indel. '''
        final_indel = ''
    ''' Combine the best snp and best indel to create the consensus allele call. '''
    consensus_allele = final_snp + final_indel
    return consensus_allele


def get_consensus_allele_calls(hi_qual_allele_calls: dict, min_depth: int, min_consensus_frac: float,
                               min_alt_frac: float, ambig_indels: bool):
    ''' This function evaluates the high-quality allele calls at each position of each genome segment
    and returns the a dict of {seq_name containing species/segment/subtype: a position-by-position
    list of consensus allele calls}. '''
    print('Calling consensus alleles...')
    consensus_allele_calls = {seq_name: list() for seq_name in hi_qual_allele_calls.keys()}
    for seq_name, seq_calls in hi_qual_allele_calls.items():
        for calls in seq_calls:
            hi_qual_depth = calls.total()
            if hi_qual_depth >= min_depth:
                call = get_consensus_allele_call(calls, min_consensus_frac, min_alt_frac, ambig_indels)
                consensus_allele_call = call
            else:
                consensus_allele_call = 'N'
            consensus_allele_calls[seq_name].append(consensus_allele_call)
    return consensus_allele_calls


def make_consensus_seq(consensus_allele_calls: list):
    ''' This function takes a list of consensus allele calls and converts it into a string
    representing the consensus sequence. '''
    consensus_seq = ''
    i = 0
    while i < len(consensus_allele_calls):
        call = consensus_allele_calls[i]
        if '-' in call:
            ''' If there is a deletion in the consensus allele, either delete the appropriate
            number of following alleles or mask them with Ns if the deletion is ambiguous. '''
            snp, indel = call.split('-')
            consensus_seq += snp
            deletion_length = len(indel)
            if all(char == 'N' for char in indel):
                ''' Ambiguous deletion. '''
                i += 1
                consensus_seq += ('N' * deletion_length)
                i += deletion_length

            else:
                ''' Unambiguous deletion. '''
                i += 1
                i += deletion_length

        elif '+' in call:
            ''' Add insertion. '''
            snp, indel = call.split('+')
            consensus_seq += snp
            consensus_seq += indel
            i += 1
        else:
            snp, indel = call, ''
            if snp in 'ATGCWSKMYRVHDBN':
                consensus_seq += snp
            i += 1
    return consensus_seq


def make_consensus_seqs(out_dir: str, out_name: str, consensus_allele_calls: dict,
                        length_ranges: pd.DataFrame, length_tolerance: float):    
    ''' This function generates a consensus sequence for each genome segment then writes those
    sequences to a FASTA file. It also checks that these sequences are within the tolerated length
    ranges. '''
    print('Generating consensus seqs...')
    consensus_seqs = {}
    for seq_name, consensus_alleles in consensus_allele_calls.items():
        consensus_seq = make_consensus_seq(consensus_alleles)
        consensus_seqs[seq_name] = consensus_seq
    ''' Write out consensus seq '''
    print('Writing consensus seqs...')
    consensus_seqs_path = os.path.join(out_dir, f'{out_name}_consensus_seqs.fa')
    with open(consensus_seqs_path, 'w') as output_file:
        for seq_name, seq in consensus_seqs.items():
            header = '|'.join(seq_name.split('|')[:4]) + '|'
            output_file.write(f'>{header}\n')
            output_file.write(seq + '\n')
    ''' Check that consensus seq lengths are within tolerated range. '''
    print('Checking length of consensus sequences...')
    bad_seqs = set()
    for header, seq in consensus_seqs.items():
        species, segment, subtype = header.split('|')[1:4]
        index_values = species, segment, subtype
        min_length = length_ranges.loc[index_values]['min_length']
        max_length = length_ranges.loc[index_values]['max_length']
        min_length = round(min_length * ((100 - length_tolerance) / 100), 0)
        max_length = round(max_length * ((100 + length_tolerance) / 100), 0)
        if not (min_length <= len(seq) <= max_length):
            bad_seqs.add((species, segment, subtype, len(seq), min_length, max_length))
    if bad_seqs != set():
        print('\nERROR: The consensus sequence generated for the following are not within the '
              'tolerated length range:')
        for species, segment, subtype, seq_length, min_length, max_length in bad_seqs:
            print(f' {species}/{segment} ({subtype})\tActual length: {seq_length}'
                  f'\tTolerated length: {min_length} to {max_length}')
        print('\nAnalysis aborted.\n')
        exit(21)


def write_mpileup_report(out_dir: str, out_name: str, consensus_allele_calls: dict, hi_qual_allele_calls: dict):
    ''' This function generates a TSV report detailing how the mpileup was parsed at each position
    of each genome segment. '''
    print('Writing mpileup report...')
    consensus_seqs = read_fasta(os.path.join(out_dir, f'{out_name}_consensus_seqs.fa'))
    mpileup_report_path = os.path.join(out_dir, f'{out_name}_mpileup_report.tsv')
    with open(mpileup_report_path, 'w') as output_file:
        line = ('seq_name',
                'consensus_pos',
                'consensus_base',
                'mpileup_pos',
                'consensus_allele',
                'total_depth',
                'allele',
                'count',
                'frac')
        line = '\t'.join(str(s) for s in line)
        output_file.write(line + '\n')
        for seq_name in consensus_allele_calls.keys():
            mpileup_position, consensus_position = 0, 0
            while mpileup_position < len(consensus_allele_calls[seq_name]):
                consensus_allele = consensus_allele_calls[seq_name][mpileup_position]
                while consensus_allele == '*':
                    alleles = hi_qual_allele_calls[seq_name][mpileup_position]
                    total_depth = alleles.total()
                    for allele, count in alleles.most_common():
                        frac = round(count / total_depth, 3)
                        line = (seq_name,
                                '!',
                                '!',
                                mpileup_position + 1,
                                consensus_allele,
                                total_depth,
                                allele,
                                count,
                                frac)
                        line = '\t'.join(str(s) for s in line)
                        output_file.write(line + '\n')
                    mpileup_position += 1
                    consensus_allele = consensus_allele_calls[seq_name][mpileup_position]
                short_seq_name = '|'.join(seq_name.split('|')[:4]) + '|'
                consensus_base = consensus_seqs[short_seq_name][consensus_position]
                alleles = hi_qual_allele_calls[seq_name][mpileup_position]
                total_depth = alleles.total()
                if total_depth == 0:
                    line = (seq_name,
                            consensus_position + 1,
                            consensus_base,
                            mpileup_position + 1,
                            'N',
                            total_depth,
                            '!',
                            '!',
                            '!')
                    line = '\t'.join(str(s) for s in line)
                    output_file.write(line + '\n')
                for allele, count in alleles.most_common():
                    frac = round(count / total_depth, 3)
                    line = (seq_name,
                            consensus_position + 1,
                            consensus_base,
                            mpileup_position + 1,
                            consensus_allele,
                            total_depth,
                            allele,
                            count,
                            frac)
                    line = '\t'.join(str(s) for s in line)
                    output_file.write(line + '\n')
                if '+' in consensus_allele:
                    for char in consensus_allele.split('+')[1]:
                        consensus_position += 1
                        short_seq_name = '|'.join(seq_name.split('|')[:4]) + '|'
                        consensus_base = consensus_seqs[short_seq_name][consensus_position]
                        line = (seq_name,
                                consensus_position + 1,
                                consensus_base,
                                '!',
                                '!',
                                '!',
                                '!',
                                '!',
                                '!')
                        line = '\t'.join(str(s) for s in line)
                        output_file.write(line + '\n')
                elif '-' in consensus_allele:
                    if all(char != 'N' for char in consensus_allele.split('-')[1]):
                        for char in consensus_allele.split('-')[1]:
                            mpileup_position += 1
                            consensus_allele = consensus_allele_calls[seq_name][mpileup_position]
                            alleles = hi_qual_allele_calls[seq_name][mpileup_position]
                            total_depth = alleles.total()
                            for allele, count in alleles.most_common():
                                frac = round(count / total_depth, 3)
                                line = (seq_name,
                                        '!',
                                        '!',
                                        mpileup_position + 1,
                                        consensus_allele,
                                        total_depth,
                                        allele,
                                        count,
                                        frac)
                                line = '\t'.join(str(s) for s in line)
                                output_file.write(line + '\n')
                mpileup_position += 1
                consensus_position += 1


def write_ambig_report(out_dir: str, out_name: str, consensus_allele_calls: dict,
                       hi_qual_allele_calls: dict):
    ''' This function generates a TSV report detailing the high-quality allele counts at each
    consensus sequence position that was masked for low coverage or ambiguity. '''
    print('Writing ambiguity report...')
    consensus_seqs = read_fasta(os.path.join(out_dir, f'{out_name}_consensus_seqs.fa'))
    ambig_report_path = os.path.join(out_dir, f'{out_name}_ambig_report.tsv')
    with open(ambig_report_path, 'w') as output_file:
        line = ('seq_name',
                'consensus_pos',
                'consensus_base',
                'total_depth',
                'allele',
                'count',
                'frac')
        line = '\t'.join(str(s) for s in line)
        output_file.write(line + '\n')
        for seq_name in consensus_allele_calls.keys():
            mpileup_position, consensus_position = 0, 0
            while mpileup_position < len(consensus_allele_calls[seq_name]):
                consensus_allele = consensus_allele_calls[seq_name][mpileup_position]
                ''' Check if the consensus allele is ambiguous, ie a degenerate base or an indel
                composed entirely of Ns. '''
                degen_base = len(consensus_allele) == 1 and consensus_allele not in 'ATGC*'
                ambig_insertion = ('-' in consensus_allele and
                                   all(char == 'N' for char in consensus_allele.split('-')[1]))
                ambig_deletion = ('+' in consensus_allele and
                                  all(char == 'N' for char in consensus_allele.split('+')[1]))
                if any((degen_base, ambig_insertion, ambig_deletion)):
                    ambig_length = len(consensus_allele.replace('-', '').replace('+', ''))
                    short_seq_name = '|'.join(seq_name.split('|')[:4]) + '|'
                    start_pos, end_pos = consensus_position, consensus_position + ambig_length
                    consensus_motif = consensus_seqs[short_seq_name][start_pos:end_pos]
                    alleles = hi_qual_allele_calls[seq_name][mpileup_position]
                    total_depth = alleles.total()
                    for allele, count in alleles.most_common():
                        frac = round(count / total_depth, 3)
                        line = (short_seq_name,
                                consensus_position + 1,
                                consensus_motif,
                                total_depth,
                                allele,
                                count,
                                frac)
                        line = '\t'.join(str(s) for s in line)
                        output_file.write(line + '\n')
                if '+' in consensus_allele:
                    for char in consensus_allele.split('+')[1]:
                        consensus_position += 1
                elif '-' in consensus_allele:
                    if all(char != 'N' for char in consensus_allele.split('-')[1]):
                        for char in consensus_allele.split('-')[1]:
                            mpileup_position += 1
                mpileup_position += 1
                consensus_position += 1


def write_consensus_seq_report(out_dir: str, out_name: str, min_depth: int):
    ''' This function generates a TSV report detailing important metrics about the consensus
    sequences for each genome segment. '''
    print('Write consensus seqs report...')
    scaffold_seqs = read_fasta(os.path.join(out_dir, f'{out_name}_scaffolds.fa'))
    mapping_refs = read_fasta(os.path.join(out_dir, f'{out_name}_mapping_refs.fa'))
    consensus_seqs = read_fasta(os.path.join(out_dir, f'{out_name}_consensus_seqs.fa'))
    consensus_seqs_report_path = os.path.join(out_dir, f'{out_name}_consensus_seqs_report.tsv')
    with open(consensus_seqs_report_path, 'w') as output_file:
        line = ('seq_name',
                'scaffold_completeness',
                'total_positions',
                'called_positions',
                'perc_called',
                'low_cov_position',
                'perc_low_cov',
                'degen_positions',
                'perc_degen',
                'ref_seq_used')
        line = '\t'.join(str(s) for s in line)
        output_file.write(line + '\n')
        for seq_name in mapping_refs.keys():
            mpileup_report_path = os.path.join(out_dir, f'{out_name}_mpileup_report.tsv')
            with open(mpileup_report_path, 'r') as input_file:
                input_file.readline()
                total_positions = 0
                called_positions = 0
                low_cov_positions = 0
                degen_positions = 0
                last_position = 0
                last_depth = 0
                for line in input_file:
                    line = line.strip().split('\t')
                    if line[0] == seq_name:
                        consensus_position = line[1]
                        if consensus_position != '!' and consensus_position != last_position:
                            total_positions += 1
                            consensus_base = line[2]
                            total_depth = int(line[5]) if line[5] != '!' else last_depth
                            if consensus_base in 'ATGC':
                                called_positions += 1
                            elif consensus_base not in 'ATGC' and total_depth >= min_depth:
                                degen_positions += 1
                            elif total_depth < min_depth:
                                low_cov_positions += 1
                            last_position = consensus_position
                            last_depth = total_depth
            perc_called = round(called_positions * 100 / total_positions, 1)
            perc_low_cov = round(low_cov_positions * 100 / total_positions, 1)
            perc_degen = round(degen_positions * 100 / total_positions, 1)
            ref_seq_used = '|'.join(seq_name.split('|')[-3:-1]) + '|'
            seq_name = '|'.join(seq_name.split('|')[:4]) + '|'
            if any((len(consensus_seqs[seq_name]) != total_positions,
                    sum(consensus_seqs[seq_name].count(base) for base in 'ATGC') != called_positions)):
                print(f'\nERROR: Report values do not match!\nAnalysis aborted.\n')
                exit(22)
            scaffold_completeness = sum(scaffold_seqs[seq_name].count(base) for base in 'ATGC')
            scaffold_completeness = scaffold_completeness * 100 / len(scaffold_seqs[seq_name])
            scaffold_completeness = round(scaffold_completeness, 1)                          
            line = (seq_name,
                    scaffold_completeness,
                    total_positions,
                    called_positions,
                    perc_called,
                    low_cov_positions,
                    perc_low_cov,
                    degen_positions,
                    perc_degen,
                    ref_seq_used)
            line = '\t'.join(str(s) for s in line)
            output_file.write(line + '\n')


def collect_garbage(out_dir: str, out_name: str):
    ''' Remove space-consuming intermediate files. '''
    shutil.rmtree(os.path.join(out_dir, 'spades_output'))
    files_to_keep = ['alignment.bam', 'alignment.bam.bai', 'scaffolds.fa', 'mapping_refs.fa',
                     'consensus_seqs.fa', 'mpileup_report.tsv', 'ambig_report.tsv',
                     'consensus_seqs_report.tsv']
    files_to_keep = [f'{out_name}_{file}' for file in files_to_keep]
    for file in os.listdir(out_dir):
        if os.path.isfile(os.path.join(out_dir, file)):
            if file not in files_to_keep:
                os.remove(os.path.join(out_dir, file))


if __name__ == '__main__':
    main()
