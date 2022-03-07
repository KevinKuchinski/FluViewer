# FluViewer
Tool for generating influenza A virus genome sequences from FASTQ data

## Installation
1. Copy the repository:
```
git clone https://github.com/KevinKuchinski/FluViewer.git
```
2. Create an environment with following dependencies installed:
- python v3.8.5
- pandas v1.3.5
- spades v3.15.3
- blast v2.12.0
- bwa v0.7.17
- samtools v1.14
- bcftools v1.14
- bedtools v2.30.0
- seqtk v1.3
3. Unzip the FluViewer_db.fa.gz file and create a blastn database:
```
makeblastdb -in FluViewer_db.fa -dbtype nucl
```

## Usage
```
python FluViewer_v_0_0_0.py -f <path_to_fwd_reads> -r <path_to_rev_reads> -d <path_to_db_file> -o <output_name> [-D <min_depth> -q <min_qual> -c <min_cov>]
```
