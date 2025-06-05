# FluViewer

FluViewer is an automated pipeline for generating influenza A and B virus (flu) genome sequences from FASTQ data. If provided with a sufficiently diverse and representative database of reference sequences, it can generate sequences regardless of host and subtype/lineage without any human intervention required.

Here is a brief outline of the FluViewer process:
1. Sequencing reads are assembled de novo into contigs
2. Contigs are aligned to a database of flu reference sequences. These alignments are used to trim contigs and roughly position them within their respective genome segments. 
3. A multiple sequence alignment in conducted on the trimmed/positioned contigs for each genome segment. The aligned contigs are collapsed into a single consensus scaffold sequence for each genome segment. 
4. Scaffolds are aligned to the reference sequence database to find their best matches. These best matches are used to fill in any missing regions in the scaffold, thereby creating mapping references.
5. Sequencing reads are mapped to the mapping references.
6. A consensus sequence is generated for each genome segment.

## Installation
1. Create a virtual environment and install the necessary dependencies using the YAML file provided in this repository. For example, if using conda:
```
conda create -n FluViewer -f FluViewer_v_0_3_x.yaml
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
FluViewer -f <path_to_fwd_reads> -r <path_to_rev_reads> -d <path_to_db_file> -n <output_name> -o <output_dir> [ <optional_args> ]
```

<b>Required arguments:</b>

-f : path to FASTQ file containing forward reads (trim sequencing adapters/primer before analysis)

-r : path to FASTQ file containing reverse reads (trim sequencing adapters/primer before analysis)

-d : path to FASTA file containing FluViewer database --- See FluViewer Database section below for more information

-n : output name (the name incorporated into output file names and consensus sequence headers)

-o : path to directory where output files are written


<b>Optional arguments:</b>

-i : Minimum sequence identity between database reference sequences and contigs (percentage, default = 90, min = 0, max = 100)

-l : Minimum length of alignment between database reference sequences and contigs (int, default = 100, min = 32)

-D : minimum read depth for base calling (int, default = 20,  min = 1)

-q : Minimum PHRED score for mapping quality and base quality during variant calling (int, default = 20, min = 0)

-c : Consensus allele fraction threshold (float, default = 0.75, min = 0, max = 1) --- See Variant Calling section below for more information

-a : Alternate allele fraction threshold (float, default = 0.05, min = 0, max = 1) --- See Variant Calling section below for more information

-L : Length tolerance for consensus sequences (percentage, default = 1, min = 0, max = 100) --- See Length Filtering section below for more information

-T : Threads used for BLAST alignments (int, default = 1, min = 1)

<b>Optional flags:</b>

-m : Allow analysis of mixed infections --- See Mixed Infection Analysis section below for more information

-A : Mask ambiguous indels --- See Variant Calling section below for more information

-g : Disable garbage collection and retain intermediate analysis files


## FluViewer Database
FluViewer requires a "database", which is a curated FASTA file of flu reference sequences. Headers for these sequences must be formatted and annotated as follows:
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
FluViewer generates five main output files:
1. A FASTA file containing consensus sequences generated for each genome segment
2. A sorted BAM file with reads mapped to the mapping references
3. An mpileup report TSV describing the alleles observed at each position (based on the sorted BAM file)
4. An ambiguous position report TSV describing all observed alleles at positions masked with ambiguous characters (including positions masked due to low coverage) 
5. A consensus sequence report TSV describing completeness metrics for each genome segment
FluViewer also generates output logs and saves key intermediate outputs, eg the FASTA file of scaffold sequences and the FASTA file of mapping references.

### Consensus sequence FASTA
Headers in the consensus sequences FASTA file have the following format:
```
>output_name|species|segment|subtype|
```

### mpileup report TSV

<b>seq_name</b> : the name of the mapping reference sequence described by this row

<b>consensus_pos</b> : the 1-indexed position in the consensus sequence described by this row

<b>consensus_base</b> : the allele selected for the consensus sequence at this position

<b>mpileup_pos</b> : the 1-indexed position in the mpileup described by this row (this can be different from the consensus_pos due to upstream indels)

<b>mpileup_pos</b> : the 1-indexed position in the mpileup described by this row (this can be different from the consensus_pos due to upstream indels) 

<b>consensus_allele</b> : the consensus allele selected for this position in the mpileup

<b>total_depth</b> : the number of reads covering this position that had basecalls meeting or exceeding the minimum quality threashold (-q)

<b>allele</b> : one of the alleles observed at this position in the mpileup

<b>count</b> : the number of reads in which the indicated allele appeared

<b>frac</b> : the fraction of reads in which the indicated allele appeared

### Ambiguous position report TSV

<b>seq_name</b> : the name of the consensus sequence described by this row

<b>consensus_pos</b> : the 1-indexed position in the consensus sequence described by this row

<b>consensus_base</b> : the allele selected for the consensus sequence at this position

<b>total_depth</b> : the number of reads covering this position that had basecalls meeting or exceeding the minimum quality threashold (-q)

<b>allele</b> : one of the alleles observed at this position in the mpileup

<b>count</b> : the number of reads in which the indicated allele appeared

<b>frac</b> : the fraction of reads in which the indicated allele appeared

### Consensus sequence report TSV

<b>seq_name</b> : the name of the consensus sequence described by this row

<b>scaffold_completeness</b> : the percentage of positions in the scaffold for the sequence that were covered by contigs assembled from the provided reads

<b>total_positions</b> : the number of positions in this consensus sequence

<b>called_positions</b> : the number of positions in this consensus sequence that were called as an A, T, G, or C

<b>perc_called</b> : the percentage of positions in this consensus sequence that were called as an A, T, G, or C

<b>low_cov_position</b> : the number of positions in this consensus sequence that were masked with Ns because they had insufficient read coverage (as set by -D)

<b>perc_low_cov</b> : the percentage of positions in this consensus sequence that were masked with Ns because they had insufficient read coverage (as set by -D)

<b>degen_positions</b> : the number of positions in this consensus sequence that were masked with ambiguity characters because a consensus allele could not be called (as set by -c and -a)

<b>perc_degen</b> : the percentage of positions in this consensus sequence that were masked with ambiguity characters because a consensus allele could not be called (as set by -c and -a)

<b>ref_seq_used</b> : the unique ID and strain name of the scaffold's best-matching reference sequence used for filling in missing regions in the scaffold (if the scaffold completeness was 100%, then this is provided pro forma as none of it was used to create the mapping reference)

## Length Filtering
For quality control purposes, FluViewer checks the length of each genome segment consensus sequence it generates to ensure that it falls within the expected range. The expected range is determined by the contents of the FluViewer database and a tolerance set by the user (-L). FluViewer will reject any genome segment sequence it generates that is n percent shorter than the shortest example of that species/segment/subtype in the database, where n is the tolerance percentage set by the user. Similarly, it will reject any genome segment sequence it generates that is n percent longer than the longest example of that species/segment/subtype in the database, where n is the tolerance percentage set by the user.

## Mixed Infection Analysis
By default, FluViewer will abort an analysis if it detects the presence of more than one subtype. This is to prevent the accidental creation of inaccurate, artefactual sequences combining reads originating from different viruses. This risk is particularly high for segments with highly conserved regions, eg M segment. Some influenza A hosts (particularly wild waterfowl) are often co-infected with multiple viruses, however. In these cases, the user can enable mixed infection analysis by setting -m. When mixed infection analysis is enabled, FluViewer will only generate consensus segments for species/segments that can be subtyped. In practice, this means it will only proceed with scaffolds whose best database matches were annotated with a subtype other than "none". Conventionally, only HA and NA segments have subtype annotations other than "none", ie H1-18 for influenza A HA segments, N1-11 for influenza A NA segments, and Victoria or Yamagata for influenza B HA and NA segments. If a user is creating their own database, they might attempt to categorize and annotate other segments into "subtypes" for the purposes of resolving mixed infections. This is left to the user's discretion, however.

## Variant Calling
FluViewer now implements its own variant calling. It does this by generating an mpileup TSV with samtools then parsing this output. For each position in the mpileup, the following steps are performed to choose a consensus allele:
1. Get all high-quality alleles from this mpileup position, ie PHRED scores meeting or exceeding the minimum quality threshold (-q).
2. Extract an SNP and indel from each observed allele in the mpiled. The SNP is the single nucleotide in the read that was observed at the current position in the mpileup. The indel is the insertion or deletion in the read that was observed between the current position and the next position in the mpileup. If there was no indel observed in a read, a dummy blank indel is extracted. 
3. Choose the best SNP for this position in the consensus sequence: If the fraction of reads containing the most common SNP meets or exceeds the consensus allele fraction threshold (-c): append this SNP to the consensus sequence. If the fraction of reads containing the most common SNP does not meet or exceed the consensus allele fraction threshold: identify the set of SNPs that each met or exceeed the alternate allele fraction (-a), determine which IUPAC ambiguity character describes this set of nucleotides, then append this ambiguity character to the consensus sequence.
4. (-A is not set and ambiguous indel masking is not enabled) Choose the best indel for between this position in the consensus sequence and the next position: If the fraction of reads containing the most common indel (including dummy blank indels) meets or exceeds the consensus allele fraction threshold: select this indel for the consensus allele. If the fraction of reads containing the most common indel does not meet or exceed the consensus allele fraction: do not add any indel to the consensus sequence between this position and the next.
4. (-A is set and ambiguous indel masking is enabled) Choose the best indel for between this position in the consensus sequence and the next position: If the fraction of reads containing the most common indel (including dummy blank indels) meets or exceeds the consensus allele fraction threshold: select this indel for the consensus allele. If the fraction of reads containing the most common indel does not meet or exceed the consensus allele fraction: identify the most common indel that meets or exceeeds the alternate allele fraction threshold. If that indel is an insertion, append a number of Ns to the consensus sequence equal to the length of the insertion. If that indel is a deletion, mask the next n positions in the consensus position with Ns, where n is the length of the deletion.
