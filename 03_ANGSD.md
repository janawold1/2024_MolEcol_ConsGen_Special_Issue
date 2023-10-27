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
|      Population     |minimum mapQ|minimum Q|minimum depth|maximum depth|minimum individuals|
| ------------------- | ---------- | ------- | ----------- | ----------- | ----------------- |
|Australian Fairy Tern|     20     |    20   |     114     |     350     |        19         |
|      Tara iti       |     20     |    20   |     114     |     350     |        19         |
|       Global        |     20     |    20   |     300     |     630     |        38         |

A missingness threshold of 0% or 10% were used when appropriate for all subsequent analyses. Below is code run for analyses performed on alignments to the tara iti reference scaffolded using the common tern assembly. All analyses were performed in a similar manner, with the variables adjusted as per the table above.  
## Removing coding Regions
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.31.0 was used to find the complement of the UCSC annotation and exclude coding regions of the genome. An additional 1kb of sequence on either side of annotations were included to reduce linkage.  

First, the annotations for autosomal chromosomes were extracted and the window size for these regions increased by 1kb on either side.  
```
zcat reference/GCA_009819605.1/genes/GCA_009819605.1_bSteHir1.pri.xenoRefGene.gtf.gz | \
    sort | \
    grep -v WNM | \
    grep -v CM020500 | \
    grep -v CM020463 | \
    grep -v CM020462 | \
    awk '{print $1"_RagTag\t"$4-1000"\t"$5+1000"\t"$2}'> reference/GCA_009819605.1/genes/all_autosomal_annotations.bed

 zcat reference/GCA_009819605.1/genes/GCA_009819605.1_bSteHir1.pri.augustus.gtf.gz | \
    sort | \
    grep -v WNM | \
    grep -v CM020500 | \
    grep -v CM020463 | \
    grep -v CM020462 | \
    awk '{print $1"_RagTag\t"$4-1000"\t"$5+1000"\t"$2}' >> reference/GCA_009819605.1/genes/all_autosomal_annotations.bed

sort -k1,1 -k2,2n reference/GCA_009819605.1/genes/all_autosomal_annotations.bed > reference/GCA_009819605.1/genes/all_autosomal_annotations_sorted.bed
```
These files were then merged, leaving roughly 22,173 annotated genes.  
```
bedtools merge \
    -i reference/GCA_009819605.1/genes/all_autosomal_annotations_sorted.bed \
    -c 1,2,3 -o distinct > reference/GCA_009819605.1/genes/all_autosomal_annotations_merged.bed
```
The complement regions were extracted with BEDtools and used to filter regions for ANGSD.  
```
awk -v OFS='\t' {'print $1,$2'} reference/TI_scaffolded_as_CT.fasta.gz.fai > reference/TI_scaffolded_as_CT_bedtools_genome.txt

bedtools complement -i reference/GCA_009819605.1/genes/all_autosomal_annotations_sorted.bed \
    -g reference/TI_scaffolded_as_CT_bedtools_genome.txt | grep -v sf | grep -v WNM | \
    grep -v CM020500 | grep -v CM020463 | grep -v CM020462 > angsd_scaffolded/TI_scaffolded_neutral_regions.bed
```
This left 23,113 regions covering 724,015,594 bp for analyses. This corresponds to roughly 60% of the genome. 

We also extracted the autosomoal chromosomes from the common tern assembly, and renamed them with the same names in the tara iti reference assembly. This was so we could polarize the site frequency spectrum with common tern as the ancestral state. 
```
faSomeRecords common_tern.fa common_tern_autosomes.bed common_tern_autosomes.fasta

awk '{print $1"\t"$1"_RagTag"} common_tern_autosomes.bed > renames.txt

while read -r line
    do
    old=$(echo $line | awk '{print $1}')
    new=$(echo $line | awk '{print $2}')
    echo "RENAMING $old TO $new"
    sed -i "s/$old/$new/g" common_tern_autosomes.fasta
done < renames.txt
```
Once the autosomes were extracted and renamed, the common tern reference had to be subsetted to be usable for SFS analyses with ANGSD and GADMA. To ensure the appropriate regions of the genome were included for polarising the SFS, the common tern assembly was aligned to the scaffolded tara iti reference with [minimap2](https://github.com/lh3/minimap2) v2.26-r1175. This is by no means a perfect solution, but helped with maximising the likelihood that the appropriate region from the common tern assembly was used.  
```
minimap2 -x asm5 TI_scaffolded_as_CT.fasta.gz common_tern_autosomes.fasta.gz -a -o tern2tern.sam
samtools view -F tern2tern.sam | samtools fasta > common_tern_subset.fasta
```
```
bgzip common_tern_subset.fasta
samtools faidx common_tern_subset.fasta.gz
```
## Inbreeding Estimates
To estimate relative levels of inbreeding and account for the high liklihood that tara iti do not conform to Hardy-Weinberg Equilibrium (HWE), we first generated population specific genotype likelihoods for input into [ngsF](https://github.com/mfumagalli/ngsTools) v1.2.0-STD. The outputs of these analyses can also be used as a prior for populations that are not in hardy-weinburg equilibrium (HWE).  
```
angsd -P 8 -b AU.list -ref $ref -out inbreeding/AU -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 350 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 8 -b TI.list -ref $ref -out inbreeding/TI -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 8 -b GLOBAL.list -ref $ref -out inbreeding/GLOBAL -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
```
This left 3,510,734 sites in the `AU` dataset, 1,089,050 sites in the `TI` dataset, and 6,015,922 sites in the `GLOBAL` dataset upon the completion of ANGSD. Next, `ngsF` was used to estimate inbreeding. First, an initial search was performed. The below was perfomed for each of the `AU`, `TI`, and `GLOBAL` data sets.  
```
zcat inbreeding/AU.glf.gz > inbreeding/AU.glf

NSITES=$(zcat inbreeding/AU.mafs.gz | tail -n+2 | wc -l)

ngsF --n_ind 19 --n_sites $NSITES --glf inbreeding/AU.glf --out inbreeding/AU_approx_indF \
    --approx_EM --init_values u --n_threads 8

ngsF --n_ind 19 --n_sites $NSITES --glf inbreeding/AU.glf --out inbreeding/AU.indF \
    --init_values inbreeding/AU_approx_indF.pars --n_threads 8

cat inbreeding/AU.indF
```
We also estimated per-individual inbreeding with [ngsF-HMM](https://github.com/mfumagalli/ngsTools) v1.1.0.
```

```
## Population Structure
0% & 10% missiness. PCA.  
```
angsd -P 8 -b GLOBAL.bamlist -ref $ref -out structure/GLOBAL_noMiss -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1

angsd -P 8 -b GLOBAL.bamlist -ref $ref -out structure/GLOBAL_0.1miss -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1
```
Running [ngsCovar](https://github.com/mfumagalli/ngsPopGen) for the data set with no missingnes and 10% missingness.  
```
gunzip structure/GLOBAL_noMiss.geno.gz
gunzip structure/GLOBAL_0.1miss.geno.gz

NSITES=$(zcat structure/GLOBAL_noMiss.mafs.gz | tail -n+2 | wc -l)
ngsCovar -probfile structure/GLOBAL_noMiss.geno \
    -outfile structure/GLOBAL_noMiss.covar -nind 38 -nsites $NSITES -call 0 -norm 0

NSITES=$(zcat structure/GLOBAL_0.1miss.mafs.gz | tail -n+2 | wc -l)
ngsCovar -probfile structure/GLOBAL_0.1miss.geno \
    -outfile structure/GLOBAL_0.1miss.covar -nind 38 -nsites $NSITES -call 0 -norm 0

Rscript -e 'write.table(cbind(seq(1,38),rep(1,38),c(rep("AU",19),rep("TI",19))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="structure/GLOBAL_noMiss.clst", quote=F)'
Rscript plotPCA.R -i structure/GLOBAL_noMiss.covar -c 1-2 -a structure/GLOBAL_noMiss.clst -o structure/GLOBAL_noMiss.pca.pdf

Rscript -e 'write.table(cbind(seq(1,38),rep(1,38),c(rep("AU",19),rep("TI",19))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="structure/GLOBAL_0.1miss.clst", quote=F)'
Rscript plotPCA.R -i structure/GLOBAL_0.1miss.covar -c 1-2 -a structure/GLOBAL_0.1miss.clst -o structure/GLOBAL_0.1miss.pca.pdf
```

### Population Structure with Inbreeding
Considering high inbreeding for tara iti
### Global Genetic Differentiation
Only the `GLOBAL` data set was used to estimate genetic distance.  
```
angsd -P 8 -b AU.list -ref $ref -out distance/AU -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 350 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1

angsd -P 8 -b TI.list -ref $ref -out distance/TI -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1

angsd -P 8 -b GLOBAL.list -ref $ref -out distance/GLOBAL -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1
```
As with the inbreeding calculations above, this command not only generates the required `*.mafs.gz`, but also the `*.geno.gz` encoded as posterior probabilities of all possible genotypes required for estimating genetic distance with [ngsDistv1.0.10](github.com/mfumagalli/ngsTools). A `pops.label` file denoting the population of origin (one entry on a new line for each samples) is necessary for estimating genetic distance.  
```
NSITES=$(zcat GLOBAL.mafs.gz | tail -n +2 | wc -l)
echo $NSITES

ngsDist -verbose 1 -geno GLOBAL.geno.gz -probs -n_ind 38 -n_sites $NSITES -labels pops.label -o distance/GLOBAL.dist
```
## Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS. Here, we are using the common tern as the ancestral state.  
```
anc=reference/common_tern_autosomes.fasta.gz
region=TI_scaffolded_neutral_regions.bed
for POP in AU TI
    do
    angsd -P 8 -b ${POP}.list -ref $ref -anc $ref -out sfs/${POP}_fold -rf ${region} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doSaf 1
        angsd -P 8 -b ${POP}.list -ref $ref -anc $anc -out sfs/${POP}_unfold -rf ${region} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doSaf 1
done

angsd -P 8 -b GLOBAL.list -ref $ref -anc $anc -out sfs/GLOBAL -rf ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doSaf 1
```

```
realSFS sfs/AU.saf.idx > sfs/AU.sfs
realSFS sfs/TI.saf.idx > sfs/TI.sfs
realSFS sfs/AU.saf.idx sfs/TI.saf.idx > sfs/AU_TI.sfs
```