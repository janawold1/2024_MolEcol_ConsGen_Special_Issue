# Population genomics with ANGSD
[ANGSDv0.940-dirty](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) (htslib:1.18) and ngsTools were used to estimate summary statistics for the *de novo* tara iti assembly and scaffolded as per the common tern reference (see [01_read_QC_and_alignment.md](github.com/janawold1/2023_EVOLAPP_Special_Issue/blob/main/01_read_QC_and_alignment.md) for details). The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. This was repeated alignments to all three reference geomes (Common tern, Tara iti scaffolded using common tern, and the unscaffolded tara iti reference). Sex chromosomes were excluded for alignments to the scaffolded assemblies (i.e., common tern, tara iti scaffolded using common tern). All analyses were performed under the SAMtools genotype likelihood model `-GL 1` unless otherwise denoted.  
```
angsd -P 16 -b GLOBAL.list -ref $ref -out ${dir}angsd_${base}/GLOBAL.qc \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -trim 0 -C 50 -baq 1 -minMapQ 15 -doQsdist 1 -doDepth 1 \
    -doCounts 1 -maxDepth 1000
```
Courtesy of scripts provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.
|              Reference              |minimum mapQ|minimum Q|minimum depth|maximum depth|
| ----------------------------------- | ---------- | ------- | ----------- | ----------- |
|        unscaffolded tara iti        |     15     |    20   |     300     |     630     |
|tara iti scaffolded using common tern|     20     |    20   |     300     |     630     |
|             common tern             |     15     |    20   |     300     |     630     |

A missingness threshold of 0% or 10% were used when appropriate for all subsequent analyses. Below is code run for analyses performed on alignments to the tara iti reference scaffolded using the common tern assembly. All analyses were performed in a similar manner, with the variables adjusted as per the table above.  
## Population Structure
0% & 10% missiness. PCA.  
Considering high inbreeding for tara iti
## Global Genetic Differentiation and Diversity
The `global` data set aligned to all three genomes was used to estimate genetic distance.  
```
angsd -P 16 -b GLOBAL.list -ref $ref -out ${base}/GLOBAL -uniqueOnly 1 -remove_bads 1 \
    -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 38 \
    -setMinDepth 300 -setMaxDepth 630 -skipTriallelic 1 -SNP_pval 1e-3 -GL 1 \
    -doCounts 1 -doMajorMinor 1 -doMaf 1 -doGeno 8 -doPost 1
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