# FluViewer

FluViewer is an automated pipeline for generating influenza A and B virus (flu) genome sequences from FASTQ data. If provided with a sufficiently diverse and representative database of reference sequences, it can generate sequences regardless of host and subtype/lineage without any human intervention required.

Here is a brief description of the FluViewer process. First, the provided reads are assembled de novo into contigs. The contigs are then aligned to a database of flu reference sequences. These alignments are used to trim contigs and roughly position them within their respective genome segment. Afterwards, a multiple sequence alignment in conducted on the trimmed/positioned contigs, generating scaffold sequences for each genome segment. Next, these scaffolds are aligned to the reference sequence database to find their best matches. These best matches are used to fill in any missing regions in the scaffold, thereby creating mapping references. The provided reads are mapped to these mapping references, then variants are called, low coverage positions are masked, and final consensus sequences are generated for each genome segment. 

## Installation
1. Create a virtual environment and install the necessary dependencies using the YAML file provided in this repository. For example, if using conda:
```
conda create -n FluViewer -f FluViewer_v_0_2_x.yaml
```

2. Activate the FluViewer environment created in the previous step. For example, if using conda:
```
conda activate FluViewer
```

3. Install the latest version of FluViewer from the Python Packing Index (PyPI).
```
pip3 install FluViewer
```

4. Download and unzip the default FluViewer DB (FluViewer_db_v_0_2_0.fa.gz) provided in this repository. Custom DBs can be created and used as well (instructions below).

## Usage
```
FluViewer -f <path_to_fwd_reads> -r <path_to_rev_reads> -d <path_to_db_file> -n <output_name> [ <optional_args> ]
```

<b>Required arguments:</b>

-f : path to FASTQ file containing forward reads (trim sequencing adapters/primer before analysis)

-r : path to FASTQ file containing reverse reads (trim sequencing adapters/primer before analysis)

-d : path to FASTA file containing FluViewer database (details below)

-n : output name (creates directory with this name for output, includes this name in output files, and in consensus sequence headers)


<b>Optional arguments:</b>

-i : Minimum sequence identity between database reference sequences and contigs (percentage, default = 90, min = 0, max = 100)

-l : Minimum length of alignment between database reference sequences and contigs (int, default = 50, min = 32)

-D : minimum read depth for base calling (int, default = 20,  min = 1)

-q : Minimum PHRED score for mapping quality and base quality during variant calling (int, default = 20, min = 0)

-v : Variant allele fraction threshold for calling variants (float, default = 0.75, min = 0, max = 1)

-V : Variant allele fraction threshold for masking ambiguous variants (float, default = 0.25, min = 0, max = 1

-L : Coverage depth limit for variant calling (int, default = 100, min = 1)

-t : Length tolerance for consensus sequences (percentage, default = 1, min = 0, max = 100)

-T : Threads used for BLAST alignments (int, default = 1, min = 1)

<b>Optional flags:</b>

-m : Allow analysis of mixed infections

-g : Disable garbage collection and retain intermediate analysis files


## FluViewer Database
FluViewer requires a curated FASTA file "database" of flu reference sequences. Headers for these sequences must be formatted and annotated as follows:
```
>unique_id|strain_name(strain_subtype)|sequence_species|sequence_segment|sequence_subtype
```
Here are some example entries:
```
>CY230322|A/Washington/32/2017(H3N2)|A|PB2|none
TCAATTATATTCAGCATGGAAAGAATAAAAGAACTACGGAATCTAATGTCGCAGTCTCGCACTCGCGA...

>JX309816|A/Singapore/TT454/2010(H1N1)|A|HA|H1
CAAAAGCAACAAAAATGAAGGCAATACTAGTAGTTCTGCTATATACATTTACAACCGCAAATGCAGACA...

>MH669720|A/Iowa/52/2018(H3N2)|A|NA|N2
AGGAAAGATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATG...

>EPI_ISL_413816|B/Iowa/08/2020(yamagata)|B|PB2|none
GTTTTCAAGATGACATTGGCTAAAATTGAATTGTTAAAGCAACTGTTAAGGGACAATGAAGCCAAAACA...

>EPI_ISL_413816|B/Iowa/08/2020(yamagata)|B|HA|Yamagata
ATTTTCTAATATCCACAAAATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCAGACCG...

>EPI_ISL_413816|B/Iowa/08/2020(yamagata)|B|NA|Yamagata
ATCTTCTCAAAAACTGAGGCAAATAGGCCAAAAATGAACAATGCTACCTTCAACTATACAAACGTTAAC...

```
For influenza A viruses and influenza B viruses, strain_subtype should reflect the HA/NA subtype or lineage of the isolate (eg H1N1 or Yamagata). 
For HA segments of influenza A viruses, segment_subtype should reflect only the HA subtype of the isolate (eg H3 for the HA segment of an H3N2 virus). Similarly, for NA segments of influenza A viruses, segment_subtype should reflect only the NA subtype of the isolate (eg N2 for the NA segment of an H3N2 virus). For HA and NA segments of influenza B viruses, segment_subtype should reflect the lineage of the isolate (eg Yamagata).

For internal segments (i.e. PB2, PB1, PA, NP, M, and NS), strain_subtype should reflect the subtypes/lineage of the isolate, but 'none' should be entered for sequence_subtype.

FluViewer will only accept reference sequences composed entirely of uppercase canonical nucleotides (i.e. A, T, G, and C).

## FluViewer Output
FluViewer generates four main output files for each library:
1. A FASTA file containing consensus sequences for each genome segments
2. A sorted BAM file with reads mapped to the mapping references generated for that library
3. A report TSV file describing segment, subtype, and sequencing metrics for each consensus sequence generated
4. Depth of coverage plots for each segment

Headers in the consensus sequences FASTA file have the following format:
```
>output_name|species|segment|subject|
```

The report TSV files contain the following columns:

<b>seq_name</b> : the name of the consensus sequence described by this row

<b>seq_length</b> : the estimated length of the genome segment described by this row

<b>reads_mapped</b> : the number of sequencing reads mapped to this segment

<b>scaffold_completeness</b> : the percentage of nucleotide positions in genome segment that were present in the scaffold that was assembled from the provided reads

<b>scaffold_completeness</b> : the percentage of nucleotide positions that were sequenced to sufficient depth in the consensus sequence generated for this genome segment

<b>low_cov_perc</b> : the percentage of nucleotide positions that were masked in the consensus sequence due to insufficient sequencing depth (determined by the depth threshold set by -D)

<b>ambig_perc</b> : the percentage of nucleotide positions that were masked in the consensus sequence generated for this genome segment because of mixed base calls (determined by the VAF thresholds set by - and -V)

<b>variant_perc</b> : the percentage of nucleotide positions in the consensus sequence that were called as variants (in relation to the mapping reference)

<b>ref_seq</b> : the unique ID and strain name of the scaffold's best-matching reference sequence used for filling in missing regions in the scaffold (if the scaffold completeness was 100%, then this is provided pro forma as none of it was used to create the mapping reference)

The depth of coverage plots contains the following elements:
- A black line indicating the depth of coverage pre-variant calling
- A grey line indicating the depth of coverage post-variant calling
- Red shading covering positions where masking was applied because coverage was too low
- Blue shading covering positions where masking was applied because base calls were ambiguous
- Green shading covering positions with variants (in relation to mapping reference)
