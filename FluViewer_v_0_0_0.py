import argparse as ap
import os as os
import subprocess as sp
import shutil as sh
import numpy as np
import pandas as pd


def main():
    version = '0.0.0'
    print(f'\nFluViewer v{version}\n')
    args = parse_arguments()
    print(f'Analyzing library {args.output}...')
    check_input_exists([args.fwd_reads, args.rev_reads, args.db])
    check_input_empty([args.fwd_reads, args.rev_reads, args.db])
    print(f'Fwd reads: {args.fwd_reads}')
    print(f'Rev reads: {args.rev_reads}')
    print()
    make_out_dir(args.output)
    print('\nGENERATING CONSENSUS CONTIG SEQS...')
    contigs = assemble_contigs(args.output, args.fwd_reads, args.rev_reads, garbage_collection=True)
    blast_out = align_contigs_to_ref_seqs(args.output, contigs, args.db)
    best_contigs = find_best_contigs(args.output, blast_out, args.min_cov, garbage_collection=True)
    best_contigs = write_best_contigs_fasta(args.output, best_contigs, contigs, garbage_collection=True)
    bam_out = map_reads_to_best_contigs(args.output, best_contigs, args.fwd_reads, args.rev_reads, garbage_collection=True)
    vcf_out = call_variants(args.output, args.min_depth, args.min_qual, best_contigs, bam_out, garbage_collection=True)
    consensus_seqs = make_consensus_seqs(args.output, bam_out, args.min_depth, args.min_cov, best_contigs, vcf_out, garbage_collection=True)
    print('\nWRITING REPORT...')
    sequenced_bases = count_sequenced_bases_in_consensus_seqs(consensus_seqs)
    consensus_seq_lengths = get_consensus_seq_lengths(consensus_seqs)
    reads_mapped_to_contigs = count_reads_mapped_to_contigs(args.output, bam_out, garbage_collection=True)
    write_reports(args.output, sequenced_bases, consensus_seq_lengths, reads_mapped_to_contigs)
    clean_headers(consensus_seqs)
    print('\nDONE: Analysis finished succesfully.\n')
    exit(0)


def parse_arguments():
    '''Parses command line arguments.'''
    parser = ap.ArgumentParser()
    parser.add_argument('-f', '--fwd_reads', type=str, required=True, help='Fwd reads')
    parser.add_argument('-r', '--rev_reads', type=str, required=True, help='Rev reads')
    parser.add_argument('-d', '--db', type=str, required=True, help='Ref seqs DB')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output dir')
    parser.add_argument('-D', '--min_depth', type=int, required=False, default=20,
                        help='Min read depth for base calling (int, default=20)')
    parser.add_argument('-q', '--min_qual', type=int, required=False, default=30,
                        help='Min PHRED qual for base calling and read mapping (int, default=30)')
    parser.add_argument('-c', '--min_cov', type=float, required=False, default=25,
                        help='Min ref seq coverage for contigs (percentage, default=25)')
    args = parser.parse_args()
    return args


def run(terminal_command, error_msg, stdout_file, stderr_file):
    '''Runs terminal command and directs stdout and stderr into indicated files.'''
    log_files = {}
    for file, dest in zip([stdout_file, stderr_file], ['stdout', 'stderr']):
        if file != None:
            log_files[dest] = open(file, 'w')
            log_files[dest].write('*' * 80 + '\n')
            log_files[dest].write('Terminal command:\n')
            log_files[dest].write(terminal_command + '\n')
            log_files[dest].write('*' * 80 + '\n')
        else:
            log_files[dest] = None
    completed_process = sp.run(terminal_command, stdout=log_files['stdout'], stderr=log_files['stderr'], shell=True)        
    for file in zip([stdout_file, stderr_file], ['stdout', 'stderr']):
        if file != None:
            log_files[dest].close()
    if completed_process.returncode != 0:
        print('\nERROR:', error_msg)
        exit(1)
    return completed_process


def check_input_exists(file_list):
    '''Checks if input files exist.'''
    for file in file_list:
        if os.path.exists(file) == False:
            print('\nERROR: Input file does not exist:')
            print(file)
            exit(1)

def check_input_empty(file_list):
    '''Checks if input files are not empty.'''
    for file in file_list:
        if os.path.getsize(file) == 0:
            print('\nERROR: Input file is empty:')
            print(file)
            exit(1)            


def make_out_dir(out_dir):
    '''Creates output dir and a dir for logs within.'''
    print('Creating directory for output...')
    if os.path.exists(out_dir) == False:
        os.mkdir(out_dir)
        logs_path = os.path.join(out_dir, 'logs')
        os.mkdir(logs_path)
    else:
        if os.path.isdir(out_dir) == False:
            print('\nERROR: Cannot create output directory because a file with that name already exists.')
            exit(1)
        if not os.listdir(out_dir):
            print('\nWARNING: Output directory already exists but is empty. Analysis will continue.')
            os.mkdir(out_dir)
            logs_path = os.path.join(out_dir, 'logs')
            os.mkdir(logs_path)
        else:
            print('\nERROR: Output directory already exists and is not empty.')
            exit(1)


def assemble_contigs(output, fwd_reads, rev_reads, contig_type='contigs', garbage_collection=True):
    '''Assmebles contigs from fwd_reads and rev_reads FASTQ files. Sends output to spades_out_dir.
    Returns path to contigs FASTA file.'''
    print('Assembling reads into contigs...')
    spades_out = os.path.join(output, output + '_spades_results')
    terminal_command = f'spades.py --rnaviral --isolate -1 {fwd_reads} -2 {rev_reads} -o {spades_out}'
    error_msg = 'spades terminated with errors while assembling reads into contigs.'
    stdout_file = os.path.join(output, 'logs', output + '_spades_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_spades_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    old_contigs = os.path.join(spades_out, f'{contig_type}.fasta')
    if os.path.exists(old_contigs) == False:
        print(f'\nDONE: No {contig_type} assembled from reads.')
        exit(0)
    contigs = os.path.join(output, output + '_contigs.fa')
    sh.copy2(old_contigs, contigs)
    if garbage_collection == True:
        sh.rmtree(spades_out)
    return contigs


def align_contigs_to_ref_seqs(output, contigs, ref_seqs_db):
    '''Align contigs to reference sequences with BLASTn. Returns path to BLASTn results in TSV file.'''
    print(f'Aligning contigs to ref seqs in {ref_seqs_db}...')
    if any([os.path.exists(ref_seqs_db + '.' + suffix) == False for suffix in ['nhr', 'nin' , 'nsq']]):
        print('\nERROR: blastn db files do not exist for ref seqs.')
        exit(1)
    blast_out = os.path.join(output, output + '_contigs_blast_results.tsv')
    terminal_command = (f'blastn -query {contigs} -db {ref_seqs_db} -outfmt'
                        f' "6 qseqid sseqid bitscore qlen slen" > {blast_out}')
    error_msg = 'blastn terminated with errors while aligning contigs to ref seqs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_contigs_blast_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    if os.path.getsize(blast_out) == 0:
        print('\nDONE: No contigs aligned to ref seqs.')
        exit(0)        
    return blast_out


def find_best_contigs(output, blast_out, min_cov, garbage_collection=True):
    '''Find best contig for each genome segment. Returns datasheet with best contigs.'''
    print('Finding best contig for each genome segment...')
    cols = 'qseqid sseqid bitscore qlen slen'.split(' ')
    blast_results = pd.read_csv(blast_out, sep='\t', names=cols)
    if garbage_collection == True:
        os.remove(blast_out)
    # Annotate alignments with segment and subtype
    blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
    blast_results['subtype'] = blast_results.apply(lambda row: row['sseqid'].split('|')[3], axis=1)
    # Keep only best alingments for each contig (by bitscore)
    best_bitscores = blast_results[['qseqid', 'bitscore']].groupby('qseqid').max().reset_index()
    blast_results = pd.merge(blast_results, best_bitscores, on=['qseqid', 'bitscore'])
    # Discard contigs whose best alignments are to multiple segments
    segment_counts = blast_results[['qseqid', 'segment']].drop_duplicates()
    segment_counts = segment_counts.groupby('qseqid').size().reset_index()
    segment_counts = segment_counts[segment_counts[0]==1][['qseqid']]
    blast_results = pd.merge(blast_results, segment_counts, on='qseqid')
    # Discard contigs whose best alignments are to multiple subtypes
    subtype_counts = blast_results[['qseqid', 'segment']].drop_duplicates()
    subtype_counts = subtype_counts.groupby('qseqid').size().reset_index()
    subtype_counts = subtype_counts[subtype_counts[0]==1][['qseqid']]
    blast_results = pd.merge(blast_results, subtype_counts, on='qseqid')
    # Discard contigs that do not provide minimum coverage of a segment
    median_slen = blast_results[['qseqid', 'slen']].groupby('qseqid').quantile(interpolation='higher').reset_index()
    blast_results = pd.merge(blast_results, median_slen, on=['qseqid', 'slen'])
    blast_results = blast_results[blast_results['qlen'] * 100 / blast_results['slen'] >= min_cov]
    # De-duplicate rows from contigs with best alignments to multiple ref seqs
    cols = ['qseqid', 'segment', 'subtype', 'slen']
    blast_results = blast_results[cols].drop_duplicates()
    if len(blast_results) == 0:
        print('DONE: No valid contigs found.')
        exit(0)
    return blast_results


def write_best_contigs_fasta(output, best_contigs_sheet, contigs, garbage_collection=True):
    '''Looks up best contigs in contigs FASTA file and writes them to their own FASTA file.
    Returns path to best contigs FASTA file.'''
    with open(contigs, 'r') as input_file:
        contig_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                contig_seqs[header] = ''
            else:
                contig_seqs[header] += line.strip()
    if garbage_collection == True:
        os.remove(contigs)
    best_contigs = os.path.join(output, output + '_best_contigs.fa')
    with open(best_contigs, 'w') as output_file:
        contig_counter = 1
        for index in best_contigs_sheet.index:
            contig_name = best_contigs_sheet['qseqid'][index]
            segment = best_contigs_sheet['segment'][index]
            subtype = best_contigs_sheet['subtype'][index]
            segment_length = best_contigs_sheet['slen'][index]
            header = f'>{output}_contig_{contig_counter}|{segment}|{subtype}|{segment_length}'
            output_file.write(header + '\n')
            output_file.write(contig_seqs[contig_name] + '\n')
            contig_counter += 1
    return best_contigs


def map_reads_to_best_contigs(output, best_contigs, fwd_reads, rev_reads, garbage_collection=True):
    '''Maps reads to best contigs with BWA mem. Filters, sorts, and indexes mappings with samtools.
    Returns path to filtered/sorted/indexed BAM file.'''
    print('Mapping reads to best contigs...')
    terminal_command = f'bwa index {best_contigs}'
    error_msg = 'bwa index terminated with errors while indexing best contigs.'
    stdout_file = os.path.join(output, 'logs', output + '_contigs_bwa_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_contigs_bwa_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    sam_out = os.path.join(output, output + '_contigs_alignment.sam')
    terminal_command = f'bwa mem {best_contigs} {fwd_reads} {rev_reads} > {sam_out}'
    error_msg = 'bwa mem terminated with errors while mapping reads to best contigs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_contigs_bwa_mem_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    bam_out = os.path.join(output, output + '_contigs_alignment_filtered_sorted.bam')
    terminal_command = f'samtools view -f 3 -F 2828 -q 30 -h {sam_out} | samtools sort -o {bam_out}'
    error_msg = 'samtools view/sort terminated with errors while filtering and sorting read mappings to best contigs.'
    stdout_file = os.path.join(output, 'logs', output + '_contigs_samtools_view_sort_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_contigs_samtools_view_sort_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    if garbage_collection == True:
        os.remove(sam_out)
    terminal_command = f'samtools index {bam_out}'
    error_msg = 'samtools index terminated with errors while indexing read mappings.'
    stdout_file = os.path.join(output, 'logs', output + '_contigs_samtools_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_contigs_samtools_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    return bam_out


def call_variants(output, min_depth, min_qual, best_contigs, bam_out, garbage_collection=True):
    '''Call variants with bcftools. Returns path to VCF file.'''
    print('Calling variants from mapped reads...')
    vcf_out = os.path.join(output, output + '_variants.vcf.gz')
    terminal_command = (f'bcftools mpileup -q {min_qual} -Q {min_qual}'
                        f' -m {min_depth} -Ou -f {best_contigs} {bam_out}'
                        f' | bcftools call --ploidy 1 -M -mv -Oz -o {vcf_out}')
    error_msg = 'bcftools mpileup/call terminated with errors while calling variants.'
    stdout_file = os.path.join(output, 'logs', output + '_bcftools_mpileup_call_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bcftools_mpileup_call_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    terminal_command = f'bcftools index {vcf_out}'
    error_msg = 'bcftools index terminated with errors while indexing variant calls.'
    stdout_file = os.path.join(output, 'logs', output + '_bcftools_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bcftools_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    return vcf_out


def make_consensus_seqs(output, bam_out, min_depth, min_cov, best_contigs, vcf_out, garbage_collection=True):
    '''Apply variants to consensus seqs with bcftools. Returns path to consensus seqs FASTA file.'''
    print('Masking low coverage positions in best contigs...')
    low_cov = os.path.join(output, output + '_low_cov.bed')
    terminal_command = f"bedtools genomecov -bga -ibam {bam_out} | awk '$4<{min_depth} {{print}}' > {low_cov}"
    error_msg = 'bedtools genomecov terminated with errors while masking low coverage positions in best contigs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_contigs_bedtools_genomecov_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    print('Generating consensus contigs seqs...')
    consensus_seqs = os.path.join(output, output + '_consensus_contig_seqs.fa')
    terminal_command = (f'cat {best_contigs} | bcftools consensus -m {low_cov} {vcf_out}'
                        f' | seqtk seq -l 0 > {consensus_seqs}')
    error_msg = 'bcftools consensus terminated with errors while applying variants.'
    stdout_file = os.path.join(output, 'logs', output + '_bcftools_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bcftools_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    if garbage_collection == True:
        files = [best_contigs + suffix for suffix in ['', '.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']]
        files += [vcf_out + suffix for suffix in ['', '.csi']]
        files += [low_cov]
        for file in files: 
            os.remove(file)
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip()
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    segment_order = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']
    header_order = sorted(seqs.keys(), key=lambda s: segment_order.index(s.split('|')[1]))
    seqs = {header: seqs[header] for header in header_order}
    with open(consensus_seqs, 'w') as output_file:
        for header, seq in seqs.items():
            seq = seq.strip('N')
            ref_seq_length = int(header.split('|')[-1])
            if len(seq) * 100 / ref_seq_length >= min_cov:
                output_file.write(header + '\n')
                output_file.write(seq + '\n')
    return consensus_seqs


def count_sequenced_bases_in_consensus_seqs(consensus_seqs):
    print('Counting sequenced bases in consensus contig seqs...')
    check_input_exists([consensus_seqs])
    check_input_empty([consensus_seqs])
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    sequenced_bases = {header: len(seq) - seq.count('N') for header, seq in seqs.items()}
    sequenced_bases = pd.DataFrame.from_dict(sequenced_bases, orient='index').reset_index()
    sequenced_bases.columns = ['contig_name', 'sequenced_bases']
    return sequenced_bases


def get_consensus_seq_lengths(consensus_seqs):
    print('Getting lengths of consensus contig seqs...')
    check_input_exists([consensus_seqs])
    check_input_empty([consensus_seqs])
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    consensus_seq_lengths = {header: len(seq) for header, seq in seqs.items()}
    consensus_seq_lengths = pd.DataFrame.from_dict(consensus_seq_lengths, orient='index').reset_index()
    consensus_seq_lengths.columns = ['contig_name', 'contig_length']
    return consensus_seq_lengths


def count_reads_mapped_to_contigs(output, bam_out, garbage_collection=True):
    idxstats = os.path.join(output, output + '_reads_mapped_to_best_contigs.tsv')
    terminal_command = f'samtools idxstats {bam_out} > {idxstats}'
    error_msg = 'samtools idxstats terminated with errors while counting reads mapped to best contigs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_samtools_idxstats_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    cols = ['contig_name', 'contig_length', 'mapped_reads', 'unmapped_reads']
    reads_mapped_to_contigs = pd.read_csv(idxstats, sep='\t', names=cols).replace('*', np.nan).dropna()
    cols = ['contig_name', 'mapped_reads']
    reads_mapped_to_contigs = reads_mapped_to_contigs[cols]
    if garbage_collection == True:
        files = [bam_out + suffix for suffix in ['', '.bai']] + [idxstats]
        for file in files:
            os.remove(file)
    return reads_mapped_to_contigs


def write_reports(output, sequenced_bases, consensus_seq_lengths, reads_mapped_to_contigs):
    print('Compiling data for report...')
    report = pd.merge(reads_mapped_to_contigs, sequenced_bases, on='contig_name')
    report = pd.merge(report, consensus_seq_lengths, on='contig_name')
    report['segment'] = report.apply(lambda row: row['contig_name'].split('|')[1], axis=1)
    report['subtype'] = report.apply(lambda row: row['contig_name'].split('|')[2], axis=1)
    report['ref_seq_length'] = report.apply(lambda row: int(row['contig_name'].split('|')[-1]), axis=1)
    report['contig_name'] = report.apply(lambda row: row['contig_name'].split('|')[0], axis=1)
    report['segment_cov'] = round(report['sequenced_bases'] * 100 / report['ref_seq_length'], 2)
    print('Writing contig report...')
    cols = ['contig_name', 'segment', 'subtype', 'mapped_reads', 'contig_length', 'sequenced_bases', 'segment_cov']
    contig_report = report[cols].drop_duplicates()
    report_path = os.path.join(output, output + '_contigs_report.tsv')
    contig_report.to_csv(report_path, sep='\t', index=False)


def clean_headers(consensus_seqs):
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip()
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    with open(consensus_seqs, 'w') as output_file:
        for header, seq in seqs.items():
            header = '|'.join(header.split('|')[:-1])
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


if __name__ == '__main__':
    main()
