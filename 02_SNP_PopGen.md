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
|       Global        |     20     |    20   |     300     |     630     |        37         |

A missingness threshold of 0% or 10% were used when appropriate for all subsequent analyses. Below is code run for analyses performed on alignments to the tara iti reference scaffolded using the common tern assembly. All analyses for each group were performed in a similar manner. Quality thresholds were adjusted as per the table above.  

Before progressing further, the SAMtools genotype likelihood model `-GL 1` and GATK genotype likelihood model `-GL 2` were trialled to see if there were differences.

## Excluding coding Regions
Gene prediction for the tara iti genome assembly was performed using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus/tree/master) v3.5.0. 
```
augustus --sample=100 --alternatives=false --temperature=3 --species=chicken $REF > reference/Katie_AUGUSTUS.gff
```
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
angsd -P 16 -b AU.list -ref $REF -out inbreeding/AU -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 16 -b TI.list -ref $REF -out inbreeding/TI -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 18 -setMinDepth 100 -setMaxDepth 280 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 16 -b GLOBAL.list -ref $REF -out inbreeding/GLOBAL -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 16 -b TI.list -ref $REF -out inbreeding/TI \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 18 -setMinDepth 100 -setMaxDepth 280 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
```
This left 1,839,732 sites in the `AU` dataset, 677,440 sites in the `TI` dataset, and XXX,XXX sites in the `GLOBAL` dataset upon the completion of ANGSD. Next, `ngsF` was used to estimate inbreeding. First, an initial search was performed. The below was perfomed for each of the `AU`, `TI`, and `GLOBAL` data sets upon the completion of ANGSD.  
```
for POP in AU TI
    do
    zcat ${ANGSD}inbreeding/${POP}.glf.gz > ${ANGSD}inbreeding/${POP}.glf
    NSITES=$(zcat ${ANGSD}inbreeding/${POP}.mafs.gz | tail -n+2 | wc -l)
    if [[ "${POP}" == "AU"]]
        then
        ngsF --n_ind 19 --n_sites $NSITES --glf inbreeding/${POP}.glf --out inbreeding/${POP}_approx_indF \
            --approx_EM --init_values u --n_threads 8
        ngsF --n_ind 19 --n_sites $NSITES --glf inbreeding/${POP}.glf --out inbreeding/${POP}.indF \
        --init_values inbreeding/${POP}_approx_indF.pars --n_threads 8
        else
        ngsF --n_ind 18 --n_sites $NSITES --glf inbreeding/${POP}.glf --out inbreeding/${POP}_approx_indF \
            --approx_EM --init_values u --n_threads 8
        ngsF --n_ind 18 --n_sites $NSITES --glf inbreeding/${POP}.glf --out inbreeding/${POP}.indF \
        --init_values inbreeding/${POP}_approx_indF.pars --n_threads 8
    fi
done

cat ${ANGSD}inbreeding/AU.indF
cat ${ANGSD}inbreeding/TI.indF
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
samtools markdup -@16 -r --write-index *_markdup_autosomes.bam *_rohan.bam
```
And the Ts/Tv ratio was estimated as a prior for ROHan with VCFtools v0.1.15 using the BCF outputs from ANGSD.
```
for POP in AU TI
    do
    angsd -P 16 -b GLOBAL.list -ref $REF -out ${ANGSD}samtools/genotypes/GLOBAL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 37 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
        -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 10 --ignore-RG 0
    vcftools --bcf ${ANGSD}gatk/genotypes/global_GATK_genotypes.bcf \
        --TsTv-summary --out ${ANGSD}gatk/genotypes/${POP}_GATK_genotypes
done
```
This resulted in a TsTv ratio of 1.871 for tara iti and 2.507 for Australian fairy tern.  

Finally, ROHan was run on for each individual.
```
./faSomeRecords reference/Katie_5kb_ragtag.fa Katie_autosomes2.bed reference/Katie_ragtag_autosomes2.fa
samtools faidx reference/Katie_ragtag_autosomes2.fa
REF=reference/Katie_ragtag_autosomes2.fa

for BAM in *_rohan.bam
    do
    BASE=$(basename $bam _rohan.bam)
    printf "STARTED RUNNING ROHAN FOR ${BASE} AT "
    date
    if [[ "$BASE" == "AU"*]]
        then
        rohan -t 16 --tstv 2.507 --size 50000 --rohmu 4.6e-9 -o output/${BASE} $REF $BAM
        else
        rohan -t 16 --tstv 1.871 --size 50000 --rohmu 4.6e-9 -o output/${BASE} $REF $BAM
    fi
    printf "FINISHED RUNNING ROHAN FOR ${BASE} AT "
    date
done
```

### Global Heterozygosity
Here, we implemented a global (genome-wide heterozygosity) method from ANGSD. Essentially, this estimate is a proportion of heterozygous genotypes / genome size (excluding regions of the genome with low confidence). Unlike other runs of ANGSD, individual BAMs are used to estimate hetereozygosity, which is simply second value in the SFS/AFS.  

```
angsd -P 16 -b GLOBAL.list -ref $ref -out ${ANGSD}gatk/genotypes/global_GATK_genotypes -uniqueOnly 1 \
    -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 -minMapQ 20 -minQ 20 -minInd 37 \
    -setMinDepth 300 -setMaxDepth 630 -doCounts 1 -skipTriallelic 1 -doBcf 1 -GL 2 \
    -doPost 1 -doMaf 1 -doGeno 10 --ignore-RG 0 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3

vcftools --bcf ${ANGSD}samtools/genotypes/global_SAMtools_genotypes.bcf \
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

printf "Sample\tHeterozygosity\tTool\tPopulation\n" > ${ANGSD}individual_het.tsv 

for TOOL in gatk samtools
    do
    for SAMP in ${ANGSD}${TOOL}/heterozygosity/*_est.ml
        do
        BASE=$(basename $SAMP _est.ml)
        HET=$(awk '{print $2/1088797119}' $SAMP)
        if [[ "$BASE" == *"AU"* ]]
            then
            printf "$BASE\t$HET\t$TOOL\tAU\n" >> ${ANGSD}individual_het.tsv
            else
            printf "$BASE\t$HET\t$TOOL\tNZ\n" >> ${ANGSD}individual_het.tsv
        fi
    done
done
```

## Population Structure
### PCA
To visualise population structure using a PCA, we first ran ANGSD using the `-doGeno 32` and `-doPost 1` options.  
```
angsd -P 16 -b ${ANGSD}GLOBAL.list -ref $REF -out ${ANGSD}samtools/structure_PCA/GLOBAL_noMiss \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 37 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1
```
Then, [ngsCovar](https://github.com/mfumagalli/ngsPopGen) was run to estimate covariance.  
```
gunzip structure_PCA/GLOBAL_noMiss.geno.gz

NSITES=$(zcat structure_PCA/GLOBAL_noMiss.mafs.gz | tail -n+2 | wc -l)
ngsCovar -probfile structure_PCA/GLOBAL.geno \
    -nind 37 -nsites $NSITES -outfile structure_PCA/GLOBAL_covar -call 0 -norm 0
```
Before plotting with the Rscripts supplied by ngsTools.  
```
Rscript -e 'write.table(cbind(seq(1,37),rep(1,37),c(rep("AU",19),rep("NZ",18))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="structure_PCA/GLOBAL_clst", quote=F)'

Rscript plotPCA.R -i structure_PCA/GLOBAL_covar -c 1-2 -a structure_PCA/GLOBAL_clst -o structure_PCA/GLOBAL_pca.pdf
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

cat GLOBAL.list | sed 's%/path/to/bams/%%g' | sed 's/_markdup_autosomes.bam//g' > ${ANGSD}mds_list

ngsDist -verbose 1 -geno ${ANGSD}structure_MDS/GLOBAL.geno.gz -probs \
    -n_ind 37 -n_sites $NSITES -labels ${ANGSD}mds_list -o ${ANGSD}structure_MDS/GLOBAL_dist
```
Extract and construct MDS.  
```
tail -n +3 ${ANGSD}structure_MDS/GLOBAL_dist | Rscript --vanilla --slave getMDS.R \
    --no_header --data_symm -n 4 -m "mds" -o ${ANGSD}structure_MDS/GLOBAL.mds

head structure_MDS/GLOBAL.mds
```
And finally plot the MDS.  
```
Rscript plotMDS.R -i ${ANGSD}distance/GLOBAL.mds \
    -c 1-2 -a ${ANGSD}structure_MDS/GLOBAL_clst \
    -o ${ANGSD}structure_MDS/GLOBAL_mds.pdf
```
### Population Structure with Inbreeding
Initial inbreeding estimates for tara iti indicate that the population is likely out of hardy-weinburg equilibrium (HWE). To account for this, relative inbreeding levels were incorporated into assessments of population structure.  

### Geographic Population Structure (GPS)
Analyses of population structure using PCAs have been shown to exhibit bias ([Elhaik et al 2022](https://doi.org/10.1038/ncomms4513)). Admixture based approaches, like those implemented in [GPS](https://github.com/arash-darzian/Geographic_Population_Structure_GPS/tree/main) have been proposed. This program leverages XXX may not want to use as it has a strong aDNA focus XXX.

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
```

```
realSFS -P 16 -bootstrap 100 -tole 1e-6 ${ANGSD}gatk/sfs/AU.saf.idx > ${ANGSD}gatk/sfs/AU_100boots.sfs

realSFS -P 16 -bootstrap 100 -tole 1e-6 ${ANGSD}gatk/sfs/TI.saf.idx > ${ANGSD}gatk/sfs/TI_100boots.sfs

realSFS -P 16 -bootstrap 100 -tole 1e-6 ${ANGSD}gatk/sfs/AU.saf.idx ${ANGSD}gatk/sfs/TI.saf.idx > ${ANGSD}gatk/sfs/GLOBAL_100boots.sfs
```
## Theta Statistics
Using the ANGSD `-doSaf 1` and realSFS outputs, theta neutrality statistics and nucleotide diversity were estimated using realSFS and associated thetaStat programmes.  
```
i=1
while read -r line
    do
    echo $i
    echo $line > ${ANGSD}samtools/sfs/AU_rep${i}.sfs
    i=$(( i + 1 ))
done < ${ANGSD}samtools/sfs/AU_100boots.sfs

i=1
while read -r line
    do
    echo $i
    echo $line > ${ANGSD}samtools/sfs/TI_rep${i}.sfs
    i=$(( i + 1 ))
done < ${ANGSD}samtools/sfs/TI_100boots.sfs

for SFS in ${ANGSD}samtools/sfs/AU_rep{1..100}.sfs 
    do
    BASE=$(basename $SFS .sfs)
    echo $BASE
    realSFS saf2theta ${ANGSD}samtools/sfs/AU.saf.idx -sfs ${SFS} -outname ${ANGSD}samtools/diversity/${BASE}
    thetaStat do_stat ${ANGSD}samtools/diversity/${BASE}.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}samtools/diversity/${BASE}_thetas_windows
done

for SFS in ${ANGSD}samtools/sfs/TI_rep{1..100}.sfs 
    do
    BASE=$(basename $sfs .sfs)
    echo $BASE
    realSFS saf2theta -P 26 ${ANGSD}samtools/sfs/TI.saf.idx -sfs ${SFS} -outname ${ANGSD}samtools/diversity/${BASE}
    thetaStat do_stat ${ANGSD}samtools/diversity/${BASE}.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}samtools/diversity/${BASE}_thetas_windows
done

for SFS in ${DIR}samtools/sfs/reps/GLOBAL_rep*.sfs
    do
    BASE=$(basename $SFS .sfs)
    realSFS fst index ${DIR}samtools/sfs/AU.saf.idx ${DIR}samtools/sfs/TI.saf.idx -sfs ${SFS} -fstout ${DIR}samtools/diversity/${BASE} -whichFst 1
    realSFS fst stats2 ${DIR}samtools/diversity/${BASE}.fst.idx -win 50000 -step 10000 > ${DIR}samtools/diversity/${BASE}_fst.txt
done
#thetaStat do_stat ${ANGSD}diversity/${POP}.thetas.idx -outnames 
#thetaStat do_stat ${ANGSD}diversity/${POP}.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}diversity/${POP}_thetasWindow
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
    -tole 1e-6 -anc $ANC -win 50000 -step 10000 -whichFST 1 > distance/total_fst2_results.tsv
```
```
 ##run in R
yc<-scan("GLOBAL_subset.sfs")
    source("plot2dSFS.R")

plot2<-function(s,...){
    dim(s)<-c(39,37)
    s[1]<-NA
    s[39,37]<-NA
s<-s/sum(s,na.rm=T)

    pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")
    pplot(s/sum(s,na.rm=T),pal=pal,...)
}

plot2(yc,ylab="AU",xlab="NZ")
#x11() # (not needed if you use Rstudio)
#plot2(yj,ylab="YRI",xlab="JPT")
#x11() #(not needed if you use Rstudio)
#plot2(jc,ylab="JPT",xlab="CEU")
```