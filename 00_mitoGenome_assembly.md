# Mitochondrial Diversity and Differentiation
To assess mitochondrial differentiation between Australian fairy tern and tara iti, we assembled mito-genomes using [NOVOPlastyv4.3.3](https://github.com/ndierckx/NOVOPlasty). Assemblies were seeded using the congeneric little tern ([*Sternula albifrons*](https://www.ncbi.nlm.nih.gov/nuccore/NC_028176.1)) mito-genome assembly and concatenated reads from multiple libraries.  
```
zcat ${dir}lib{1,2}/AU01_*_R1.fq.gz > little_tern_seeded_mtDNA/AU01_R1.fastq
zcat ${dir}lib{1,2}/AU01_*_R2.fq.gz > little_tern_seeded_mtDNA/AU01_R2.fastq

zcat ${dir}lib{1,2}/SP01_*_R1.fq.gz > little_tern_seeded_mtDNA/SP01_R1.fastq
zcat ${dir}lib{1,2}/SP01_*_R2.fq.gz > little_tern_seeded_mtDNA/SP01_R2.fastq
```
## NOVOPlasty Assembly
Below is an example of the config file used to run NOVOPlasty for one Australian fairy tern and one tara iti. Output directory, input files, read insert size and read length were adjusted for each sample.
```
Project:
-----------------------
Project name          = Australian_fairy_tern
Type                  = mito
Genome Range          = 15000-22000
K-mer                 = 33
Max memory            = 55
Extended log          = 1
Save assembled reads  = no
Seed Input            = reference/little_tern_mitoGenome.fasta
Extend seed directly  = no
Reference sequence    =
Variance detection    = no
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 123
Insert size           = 311
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = little_tern_seeded_mtDNA/AU01_R1.fastq
Reverse reads         = little_tern_seeded_mtDNA/AU01_R2.fastq
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              = no

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = yes
Reduce ambigious N's  = yes
Output path           = little_tern_seeded_mtDNA/AU01_mitoassembly/
```
