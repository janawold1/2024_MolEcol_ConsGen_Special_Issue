# Population genomics with ANGSD
[ANGSDv0.935](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) and ngsTools were used to estimate summary statistics for Australian fairy tern (*Sternual nereis nereis*) sampled from Western Australia and tara iti (*S. nereis davisae*) from Northland, NZ. The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. PCR duplicates were marked and sex chromosomes were excluded for SNP-based population analyses.  
```
ls ${DIR}*_markdup_autosomes.bam > ${ANGSD}GLOBAL.list
ls ${DIR}{SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}TI.list
ls ${DIR}AU*_markdup_autosomes.bam > ${ANGSD}AU.list

REF=${DIR}Katies_genome/Katie_5kb_ragtag.fa

for POP in GLOBAL AU TI
    do
    angsd -P 16 -b ${ANGSD}${POP}.list -ref $REF -out ${ANGSD}qc/${POP}.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 800
done
```
Courtesy of scripts provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.

|      Population     |minimum mapQ|minimum Q|minimum depth|maximum depth|minimum individuals|
| ------------------- | ---------- | ------- | ----------- | ----------- | ----------------- |
|Australian Fairy Tern|     20     |    20   |     170     |     360     |        19         |
|      Tara iti       |     20     |    20   |     100     |     280     |        18         |
|       Global        |     20     |    20   |     300     |     740     |        37         |

A missingness threshold of 0% or 10% were used when appropriate for all subsequent analyses. Below is code run for analyses performed on alignments to the tara iti reference scaffolded using the common tern assembly. All analyses for each group were performed in a similar manner. Quality thresholds were adjusted as per the table above.  

Before progressing further, the SAMtools genotype likelihood model `-GL 1` and GATK genotype likelihood model `-GL 2` were trialled to see if there were differences.

## Excluding coding Regions
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.31.0 was used to find the complement of the UCSC annotation and exclude coding regions of the genome. An additional 1kb of sequence on either side of annotations were included to reduce linkage.  

First, the annotations for autosomal chromosomes were extracted and the window size for these regions increased by 1kb on either side for each annotation file.  
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
No regions within 1kb of the start/end of the chromosome were retained to account for poor assembly at repetitive regions of the genome. This left 22,167 regions covering 615,256,912 bp for analyses. This corresponds to roughly 52% of the genome.  

This file was then indexed for use with the `angsd -sites` option.  
```
angsd sites index TI_scaffolded_neutral_regions.bed
```

We also extracted the autosomoal chromosomes from the common tern assembly, and renamed them with the same names in the tara iti reference assembly. This was so we could polarize the site frequency spectrum with common tern as the ancestral state. 

### Test for Sex Specific differences?

## Common tern ancestral alleles
It is useful to have the ancestral state for some of the analyses below. To this end, we aligned short reads from the common tern genome assembly to the tara iti genome and converted to a Fasta file. We opted this approach over using the common tern genome assembly directly as the chromosomes had different sizes and were not compatible with ANGSD methods ([brief discussion here](https://www.biostars.org/p/298013/)). This also has the benefit of ensuring that comparisons are made between the same regions, regardless of potential rearrangements. The reads for common tern were trimmed, aligned, and duplicates marked in the same manner as the fairy tern population short-reads. Would only work when using a file as input to command as below.  
```
ls ${DIR}bSteHir1_markdup_autosomes.bam > ${ANGSD}CT.list

angsd -P 16 -doFasta 1 -doCounts 1 -out ${ANGSD}bSteHir1_ancestral -b ${ANGSD}CT.list
```
## Inbreeding Estimates
### F statistics
To estimate relative levels of inbreeding and account for the high likelihood that tara iti do not conform to Hardy-Weinberg Equilibrium (HWE), we first generated population specific genotype likelihoods for input into [ngsF](https://github.com/mfumagalli/ngsTools) v1.2.0-STD. The outputs of these analyses can also be used as a prior for populations that are not in hardy-weinburg equilibrium (HWE).  
```
angsd -P 16 -b AU.list -ref $REF -out inbreeding/AU \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 16 -b TI.list -ref $REF -out inbreeding/TI -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 18 -setMinDepth 100 -setMaxDepth 280 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 16 -b GLOBAL.list -ref $REF -out inbreeding/GLOBAL -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 650 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
```
This left 980,578 sites in the `AU` dataset, 297,787 sites in the `TI` dataset, and 1,665,527 sites in the `GLOBAL` dataset upon the completion of ANGSD. Next, `ngsF` was used to estimate inbreeding. First, an initial search was performed. The below was perfomed for each of the `AU`, `TI`, and `GLOBAL` data sets upon the completion of ANGSD.  
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
### Relatedness
We estimated relatedness among individuals with [ngsRelate](https://github.com/ANGSD/NgsRelate).  
```
for POP in AU TI
    do
    zcat inbreeding/${POP}.mafs.gz | cut -f 6 | sed 1d > inbreeding/${POP}_freqs
    NSITES=$(zcat inbreeding ${POP}.mafs.gz | tail -n +2 | wc -l)
    ngsRelate -g inbreeding/${POP}.glf -n 19 -f inbreeding/${POP}_freqs -O inbreeding/${POP}_relatedness
done
```

### Runs of Homozygosity
Two methods were used to estimate runs of homozygosity, [ROHAN](https://github.com/grenaud/rohan) and an approach using ANGSD.  

For ROHAN, all PCR-duplicates were removed as recommended.  
```
samtools markdup -@16 -r --write-index *_nodup_autosomes.cram *_rohan.bam
```
And the programme was run on for each individual.
```
for bam in *_rohan.bam
    do
    base=$(basename $bam _rohan.bam)
    echo "RUNNING ROHAN FOR ${base}"
    rohan --size 50000 --rohmu 4.6e-9 -t 16 --tstv 2.36 -o output/${base} $ref $bam
done
```

### Global Heterozygosity
Here, we implemented a global (genome-wide heterozygosity) method from ANGSD. Essentially, this estimate is a proportion of heterozygous genotypes / genome size (excluding regions of the genome with low confidence). Unlike other runs of ANGSD, individual BAMs are used to estimate hetereozygosity, which is simply second value in the SFS/AFS.  

For other programs and to estimate a Ts/Tv rate for Rohan, a BCF was produced using ANGSD.  
```
angsd -P 8 -b GLOBAL.list -ref $ref -out ${ANGSD}samtools/genotypes/global_SAMtools_genotypes -uniqueOnly 1 \
    -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 -minMapQ 20 -minQ 20 -minInd 38 \
    -setMinDepth 300 -setMaxDepth 630 -doCounts 1 -skipTriallelic 1 -doBcf 1 -GL 1 \
    -doPost 1 -doMaf 1 -doGeno 10 --ignore-RG 0 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3

angsd -P 8 -b GLOBAL.list -ref $ref -out ${ANGSD}gatk/genotypes/global_GATK_genotypes -uniqueOnly 1 \
    -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 -minMapQ 20 -minQ 20 -minInd 38 \
    -setMinDepth 300 -setMaxDepth 630 -doCounts 1 -skipTriallelic 1 -doBcf 1 -GL 2 \
    -doPost 1 -doMaf 1 -doGeno 10 --ignore-RG 0 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3

vcftools --vcf ${ANGSD}gatk/genotypes/global_GATK_genotypes \
    --TsTv-sumamry \
    --out ${ANGSD}gatk/genotypes/global_GATK_genotypes

vcftools --vcf ${ANGSD}samtools/genotypes/global_SAMtools_genotypes \
    --TsTv-sumamry \
    --out ${ANGSD}samtools/genotypes/global_SAMtools_genotypes
```
```
for SAMP in ${BAM}*_markdup_autosomes.bam
    do
    BASE=$(basename $SAMP _markdup_autosomes.bam)
    angsd -i $SAMP -anc $ANC -ref $REF -out ${ANGSD}heterozygosity/${BASE} -dosaf 1 -GL 1 -doCounts 1
    realSFS ${ANGSD}heterozygosity/${BASE}.saf.idx > heterozygosity/${BASE}_est.ml
done
```
Once the SFS was estimated for each individual, the number of sites was estimated from the sum of all scaffold sizes included in the bam file.
```

```
Individual heterozygosity was then output to a file.
```
SIZE=$(awk '{sum+=$3}; END {print sum}' Katie_autosomes2.bed)
echo $SIZE



for tool in gatk samtools
    do
    for SAMP in ${ANGSD}${tool}/heterozygosity/*_est.ml
        do
        BASE=$(basename $SAMP _est.ml)
        HET=$(awk '{print $2/1088797119}' $SAMP)
        printf "$BASE\t$HET\t$tool\n" >> ${ANGSD}${tool}/heterozygosity/individual_het.tsv
    done
done

printf "Sample\tHeterozygosity\tTool\n" > ${ANGSD}heterozygosity_summary.tsv 

cat ${ANGSD}gatk/heterozygosity/individual_het.tsv >> ${ANGSD}heterozygosity_summary.tsv
cat ${ANGSD}samtools/heterozygosity/individual_het.tsv >> ${ANGSD}heterozygosity_summary.tsv
```

## Population Structure
### PCA
To visualise population structure using a PCA, we first ran ANGSD using the `-doGeno 32` and `-doPost 1` options.  
```
angsd -P 16 -b ${ANGSD}GLOBAL.list -ref $REF -out ${ANGSD}structure_PCA/GLOBAL_noMiss \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1
```
Then, [ngsCovar](https://github.com/mfumagalli/ngsPopGen) was run to estimate covariance.  
```
gunzip structure_PCA/GLOBAL_noMiss.geno.gz

NSITES=$(zcat structure_PCA/GLOBAL_noMiss.mafs.gz | tail -n+2 | wc -l)
ngsCovar -probfile structure_PCA/GLOBAL_noMiss.geno \
    -outfile structure_PCA/GLOBAL_noMiss.covar -nind 38 -nsites $NSITES -call 0 -norm 0
```
Before plotting with the Rscripts supplied by ngsTools.  
```
Rscript -e 'write.table(cbind(seq(1,38),rep(1,38),c(rep("AU",19),rep("TI",19))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="structure_PCA/GLOBAL_noMiss.clst", quote=F)'
Rscript plotPCA.R -i structure_PCA/GLOBAL_noMiss.covar -c 1-2 -a structure_PCA/GLOBAL_noMiss.clst -o structure_PCA/GLOBAL_noMiss.pca.pdf
```
### MDS
To construct a MDS for fairy tern populations, we first ran ANGSD as below. Notably, the only change is the `-doGeno 8` flag.  
```
angsd -P 16 -b GLOBAL.list -ref $REF -out ${ANGSD}structure_MDS/GLOBAL \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1
```
This command not only generates the required `*.mafs.gz`, but also the `*.geno.gz`, which contains posterior probabilities of all possible genotypes required for estimating genetic distance with [ngsDistv1.0.10](github.com/mfumagalli/ngsTools). In addition, a `pops.label` file denoting the population of origin (one entry on a new line for each samples) is necessary for estimating genetic distance.  
```
NSITES=$(zcat GLOBAL.mafs.gz | tail -n +2 | wc -l)
echo $NSITES

ngsDist -verbose 1 -geno ${ANGSD}distance/GLOBAL.geno.gz -probs \
    -n_ind 38 -n_sites $NSITES -labels ${ANGSD}pops.label -o ${ANGSD}distance/GLOBAL.dist
```
Extract and construct MDS.  
```
tail -n +3 ${ANGSD}distance/GLOBAL.dist | Rscript --vanilla --slave getMDS.R \
    --no_header --data_symm -n 4 -m "mds" -o ${ANGSD}distance/GLOBAL.mds

head distance/GLOBAL.mds
```
And finally plot the MDS.  
```
Rscript plotMDS.R -i ${ANGSD}distance/GLOBAL.mds -c 1-2 -a ${ANGSD}structure/GLOBAL_noMiss.clst -o ${ANGSD}distance/GLOBAL_mds.pdf
```
### Population Structure with Inbreeding
Initial inbreeding estimates for tara iti indicate that the population is likely out of hardy-weinburg equilibrium (HWE). To account for this, relative inbreeding levels were incorporated into assessments of population structure.  

### Geographic Population Structure (GPS)
Analyses of population structure using PCAs have been shown to exhibit bias ([Elhaik et al 2022](https://doi.org/10.1038/ncomms4513)). Admixture based appraoches, like those implemented in [GPS](https://github.com/arash-darzian/Geographic_Population_Structure_GPS/tree/main) have been proposed. This program leverages XXX may not want to use as it has a strong aDNA focus XXX.

## Summary Statistics
### Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS. Here, we are using the common tern as the ancestral state.  
```
ANC=${ANGSD}bSteHir1_ancestral.fasta
region=TI_scaffolded_neutral_regions.bed

angsd -P 16 -b ${ANGSD}AU.list -ref $REF -anc $ANC -out ${ANGSD}gatk/sfs/AU \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
    -GL 2 -doSaf 1

angsd -P 16 -b ${ANGSD}TI.list -ref $REF -anc $ANC -out ${ANGSD}gatk/sfs/TI \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ20 -minInd 18 -setMinDepth 100 -setMaxDepth 280 -doCounts 1 \
    -GL 2 -doSaf 1

angsd -P 8 -b ${ANGSD}GLOBAL.list -ref $REF -anc $ANC -out ${ANGSD}gatk/sfs/GLOBAL \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doSaf 1
```

```
realSFS -fold 0 -bootstrap 100 -tole 1e-6 ${ANGSD}gatk/sfs/AU_unfolded.saf.idx > ${ANGSD}gatk/sfs/AU_realSFS_unfolded_100boots.sfs

realSFS -fold 0 -bootstrap 100 -tole 1e-6 ${ANGSD}gatk/sfs/TI_unfolded.saf.idx > ${ANGSD}gatk/sfs/TI_realSFS_unfolded_100boots.sfs

realSFS -fold 0 -bootstrap 100 -tole 1e-6 ${ANGSD}gatk/sfs/AU_unfolded.saf.idx ${ANGSD}gatk/sfs/TI_fold.saf.idx > ${ANGSD}gatk/sfs/GLOBAL_realSFS_unfolded_100boots.sfs
```
## Theta Statistics
Using the ANGSD `-doSaf 1` and realSFS outputs, theta neutrality statistics were estimated using realSFS and associated thetaStat programmes.  
```
for POP in AU TI total
    do
    tail -n1 ${ANGSD}gatk/sfs/${POP}_unfold_realSFS_unfolded_100boots.sfs > ${ANGSD}gatk/diversity/${POP}_unfold.sfs
    realSFS saf2theta sfs/${POP}_unfold.saf.idx -outname diversity/${POP}_fold -sfs diversity/${POP}_unfold.sfs -fold 0
    thetaStat do_stat ${ANGSD}gatk/diversity/${POP}_unfold.thetas.idx
    thetaStat do_stat ${ANGSD}gatk/diversity/${POP}_unfold.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}gatk/diversity/${POP}_unfold_thetasWindow
done
```
Only the Watterson's Theta and Tajima's D can be estimated with a folded SFS. The mean value of Tajima's D was estimated by the below.  

Genetic differentiation.  
```
realSFS fst index sfs/AU.saf.idx sfs/TI.saf.idx -sfs diversity/total_fold.sfs -fstout distance/total_fst

realSFS fst stats distance/total_fst.idx \
    -tole 1e-6 -win 50000 -step 10000 -whichFst 1 > distance/total_fst_results.tsv
```
This yielded the result:  

```
-> Assuming idxname:distance/test.fst.idx
        -> Assuming .fst.gz file: distance/test.fst.gz
        -> args: tole:0.000001 nthreads:4 maxiter:100 nsites(block):0 start:(null) chr:(null) start:-1 stop:-1 fstout:(null) oldout:0 seed:-1 bootstrap:0 resample_chr:0 whichFst:1 fold:1 ref:(null) anc:reference/TI_scaffolded_as_CT.fasta.gz
        -> FST.Unweight[nObs:387393642]:0.002465 Fst.Weight:0.857685
```

We also ran:
```
realSFS fst stats2 distance/total_fst.idx \
    -tole 1e-6 -fold 1 -anc $ref -win 50000 -step 10000 -whichFST 1 > distance/total_fst2_results.tsv
```

```
angsd -P 8 -b AU.list -ref ${ref} -anc ${ref} -out diversity/AU_folded -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
    -GL 1 -doSaf 1 -doThetas 1 -pest sfs/AU_fold.sfs
```