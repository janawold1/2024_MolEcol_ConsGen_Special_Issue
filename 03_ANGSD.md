# Population genomics with ANGSD
ANGSDvX.X and ngsTools were used to estimate summary statistics for the *de novo* tara iti assembly and scaffolded as per the common tern reference (see [01_read_QC_and_alignment.md](github.com/janawold1/2023_EVOLAPP_Special_Issue)). The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. All analyses were performed under the SAMtools genotype likelihood model `-GL 1` unless otherwise appropriate.  
```
for ref in a b c 
    do
    base=$(basename $ref .fasta)
    angsd -P 8 -b GLOBAL.list -ref $ref -out ${base}/GLOBAL.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -C 50 -baq 1 -minMap 15 -doQsdist 1 -doDepth 1 \
        -doCounts 1 -maxDepths 2200
done
```
Courtesy of scripts provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. For all subsequent analyses, a minimum mapping quality of XXX, minimum base quality of XXX, minimum depth of XXX, and a maximum depth of XXX are implemented. Depending on the analysis, a missingness threshold of 0% or 10% was used.  

## Population Structure
0% & 10% missiness. PCA. 
Considering high inbreeding for tara iti
## Genetic Differentiation
The `global` data set aligned to all three genomes was used to estimate genetic distance.  
```
for ref in a b c
    do
    base=$(basename $ref .fasta)
    angsd -P 8 -b GLOBAL.list -ref $ref -out ${base}/GLOBAL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -C 50 -baq 1 -minMapQ 15 -minQ 20 -minInd 22 \
        -setMinDepth 144 -setMaxDepth 1694 -skipTriallelic 1 \
        -SNP_pval 1e-3 -GL 1 -doCounts 1 -doMajorMinor 1 \
        -doMaf 1 -doGeno 8 -doPost 1
done
```
This generates the required `*.mafs.gz` and `*.geno.gz` files for estimating genetic distance with [ngsDistvX.X](github.com/mfumagalli/ngsTools). A `pops.label` file denoting the population of origin (one entry on a new line for each samples) is necessary for estimating genetic distance.  
```
NSITES=$(zcat GLOBAL.mafs.gz | tail -n +2 | wc -l)
echo $NSITES

ngsDist -verbose 1 -geno GLOBAL.geno.gz \
    -probs -n_ind 22 -n_sites $NSITES \
    -labels pops.label -o GLOBAL.dist -P 8
```
## Population Diversity and Inbreeding
```
for POP in AU TI GLOBAL
    do
    angsd -P 8 -b ${POP}.list -ref $ref -anc $ref \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \ 
        -trim 0 -C 50 -minMapQ 20 -minQ 20 -min Ind 11 \
        -skipTriallelic 1 -setMinDepth 66 -setMaxDepth 770 /
        -GL 1 -doMajorMinor 1 -doMaf -doCounts 1 -doSaf 1 XXXX
done
```
## Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS.  
```
realSFS AU.saf.idx > AU.sfs
realSFS TI.saf.idx > TI.sfs
realSFS AU.saf.idx TI.saf.idx > AU_TI.sfs
```