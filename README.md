# FluViewer

FluViewer is an automated pipeline for generating influenza A virus (IAV) genome sequences from FASTQ data. If provided with a sufficiently diverse and representative database of IAV reference sequences, it can generate sequences regardless of host and subtype without any human intervention required.

Here is a brief description of the FluViewer process. First, the provided reads are normalized and downsampled using a kmer-based approach to reduce any excessive coverage of certain genome regions. Next, the normalized/downsampled reads are assembled de novo into contigs. The contigs are then aligned to a database of IAV reference sequences. These alignments are used to trim contigs and roughly position them within their respective genome segment. Afterwards, a multiple sequencing alignment in conducted on the trimmed/positioned contigs, generating scaffold sequences for each IAV genome segment. Next, these scaffolds are aligned to the IAV reference sequence database to find their best matches. These best matches are used to fill in any missing regions in the scaffold, creating mapping references. The normalized/downsampled reads are mapped to these mapping references, then variants are called and the final consensus genomes are produced. 

## Installation
1. Create a virtual environment and install the necessary dependencies using the YAML file provided in this repository. For example, if using conda:
```
conda create -n FluViewer -f FluViewer_v_0_1_x.yaml
```

2. Activate the FluViewer environment created in the previous step. For example, if using conda:
```
conda activate FluViewer
```

3. Install the latest version of FluViewer from the Python Packing Index (PyPI).
```
pip3 install FluViewer
```

4. Download and unzip the default FluViewer DB (FluViewer_db.fa.gz) provided in this repository. Custom DBs can be created and used as well (instructions below).

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

-v : Variant allele fraction threshold for calling variants (float, default = 0.95, min = 0, max = 1)

-V : Variant allele fraction threshold for masking ambiguous variants (float, default = 0.25, min = 0, max = 1

-N : Target depth for pre-normalization of reads (int, default = 200, min = 1)

-T : Threads used for BLAST alignments (int, default = 1, min = 1)


<b>Optional flags:</b>

-g : Disable garbage collection and retain intermediate analysis files


## FluViewer Database
FluViewer requires a curated FASTA file "database" of IAV reference sequences. Headers for these sequences must be formatted and annotated as follows:
```
>unique_id|strain_name(strain_subtype)|sequence_segment|sequence_subtype
```
Here are some example entries:
```
>CY230322|A/Washington/32/2017(H3N2)|PB2|none
TCAATTATATTCAGCATGGAAAGAATAAAAGAACTACGGAATCTAATGTCGCAGTCTCGCACTCGCGA...

>JX309816|A/Singapore/TT454/2010(H1N1)|HA|H1
CAAAAGCAACAAAAATGAAGGCAATACTAGTAGTTCTGCTATATACATTTACAACCGCAAATGCAGACA...

>MH669720|A/Iowa/52/2018(H3N2)|NA|N2
AGGAAAGATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATG...
```
For HA and NA segments, strain_subtype should reflect the HA and NA subtypes of the isolate (eg H1N1), but sequence_subtype should only indicate the HA or NA subtype of the segment sequence of the entry (eg H1 for an HA sequence or N1 for an NA sequence).

For internal segments (i.e. PB2, PB1, PA, NP, M, and NS), strain_subtype should reflect the HA/NA subtypes of the isolate, but 'none' should be entered for sequence_subtype. If strain_subtype is unknown, 'none' should be entered there as well.

FluViewer will only accept reference sequences composed entirely of uppercase canonical nucleotides (i.e. A, T, G, and C).

## FluViewer Output
FluViewer generates four main output files for each library:
1. A FASTA file containing consensus sequences for the IAV genome segments
2. A sorted BAM file with reads mapped to the mapping references generated for that library (the mapping reference is also retained)
3. A report TSV file describing segment, subtype, and sequencing metrics for each consensus sequence generated
4. Depth of coverage plots for each segment

Headers in the FASTA file have the following format:
```
>output_name|segment|subject
```


The report TSV files contain the following columns:

<b>seq_name</b> : the name of the consensus sequence described by this row

<b>segment</b> : IAV genome segment (PB2, PB1, PA, HA, NP, NA, M, NS)

<b>subtype</b> : HA or NA subtype ("none" for internal segments)

<b>reads_mapped</b> : the number of sequencing reads mapped to this segment (post-normalization/downsampling)

<b>seq_length</b> : the length (in nucleotides) of the consensus sequence generated by FluViewer

<b>scaffold_completeness</b> : the number of nucleotide positions in the scaffold that were assembled from the provided reads (post-normalization/downsampling)

<b>consensus_completeness</b> : the number of nucleotide positions in the consensus with a succesful base call (e.g. A, T, G, or C)

<b>ref_seq_used</b> : the unique ID and strain name of the scaffold's best-matching reference sequence used for filling in missing regions in the scaffold (if the scaffold completeness was 100%, then this is provided pro forma as none of it was used to create the mapping reference)


The depth of coverage plots contains the following elements:
- A black line indicating the depth of coverage pre-variant calling
- A blue line indicating the depth of coverage post-variant calling
- Red shading covering positions where coverage was too low for base calling
- Orange lines indicating positions where excess variation resulted in an ambiguous base call
- Blue lines indicating positions where a variant was called
