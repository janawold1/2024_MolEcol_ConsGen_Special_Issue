# Population genomics with ANGSD
[ANGSDv0.935](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) and ngsTools were used to estimate summary statistics for Australian fairy tern (*Sternual nereis nereis*) sampled from Western Australia and tara iti (*S. nereis davisae*) from Northland, NZ. The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

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
Courtesy of scripts provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.

|        Population        |minimum mapQ|minimum Q|minimum depth|maximum depth|Number of individuals|
| ------------------------ | ---------- | ------- | ----------- | ----------- | ------------------- |
|Australian Fairy Tern (WA)|     20     |    20   |     200     |     420     |         19          |
|         Tara iti         |     20     |    20   |     90      |     280     |         15          |
|     Global fairy tern    |     20     |    20   |     204     |     630     |         34          |
|           Kakī           |     20     |    20   |     XXX     |     XXX     |         XX          |

 Below is code run for analyses performed on alignments of fairy terns to the tara iti reference scaffolded using the common tern assembly (see [00_genome_assembly.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/00_genome_assembly.md)) and alignments of kakī to the high-quality reference assembly for this species (see [here](https://www.genomics-aotearoa.org.nz/our-work/completed-projects/high-quality-genomes) for details). All analyses for each group were performed in a similar manner for downstream comparisons. Quality thresholds were adjusted as per the table above.  

These thresholds left X,XXX,XXX sites in the WA population of Australian fairy tern (`AU`) dataset, XXX,XXX sites in the tara iti dataset (`TI`), XXX,XXX sites in the fairy tern dataset (`GLOBAL`), and X,XXX,XXX in the kakī dataset (`KI`) for analyses.  

## Excluding putative coding Regions
Neither the tara iti or kakī reference assemblies have high-quality transcriptomes. This is significant given the potential of selection to influence some of the analyses performed here (e.g., estimates of historical N<sub>e</sub>). To assess the potential consequences of the inclusion of these regions, and to have a 'first look' at the genomic diversity around coding regions, we performed *ab initio* gene prediction using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus/tree/master) v3.5.0.  
```
augustus --sample=100 --alternatives=false --temperature=3 --species=chicken $REF > reference/Katie_AUGUSTUS.gff
```
Autosomal chromosomes were extracted from these predictions, and [BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.30.0 was then used to sort, and merge these putative gene regions. An additional 1kb of sequence on either side of annotations were included to reduce linkage.  

First, the annotations for .  
```
grep -v "#" reference/gene_predictions/Katie_AUGUSTUS.gff | \
    grep -v WNMW0 | \
    grep -v CM020459 | \
    grep -v CM020460 | \
    grep -v CM020461 | \
    grep -v CM020462 | \
    grep -v CM020463 | \
    grep -v contig_ | \
    grep -v ntLink_ | \
    bedtools sort -i - > angsd/augustus_autosomal_predictions.bed

bedtools merge -i angsd/augustus_autosomal_predictions.bed > angsd/augustus_autosomal_predictions_merged.bed
```
 To account for linkage, the window size for these regions increased by 1kb on either side. The file was adusted in cases where the addition of a 1kb buffer extended beyond the start (n = 11) or end of the chromosome (n = 7).
```
awk '{print $1"\t"$2-1000"\t"$3+1000}' angsd/augustus_autosomal_predictions_merged.bed > angsd/augustus_autosomal_predictions_merged_add1kb.bed
```
To define putatively neutral sites for analyses, we extracted the complement regions as below.  
```
awk -v OFS='\t' {'print $1,$2'} reference/Katie_5kb_ragtag.fa.fai | \
    grep -v WNMW0 | \
    grep -v CM020459 | \
    grep -v CM020460 | \
    grep -v CM020461 | \
    grep -v CM020462 | \
    grep -v CM020463 | \
    grep -v contig_ | \
    grep -v ntLink_ > reference/Katie_autosome_lengths.txt

bedtools complement \
    -i angsd/augustus_autosomal_predictions_merged_add1kb.bed \
    -g reference/Katie_autosome_lengths.txt > angsd/Katie_neutral_sites.bed
```
This left 27,337 regions covering 562,663,009 bp for analyses. This corresponds to roughly 47% of the genome.  

This file was then indexed for ANGSD with `angsd sites`.  
```
angsd sites index TI_scaffolded_neutral_regions.bed
```
## Common tern ancestral alleles
It is important to have the ancestral state for some of the analyses below. To this end, we aligned short reads from the common tern genome assembly to the tara iti genome and converted to a Fasta file. We opted this approach over using the common tern genome assembly directly as the chromosomes had different sizes and were not compatible with ANGSD methods ([brief discussion here](https://www.biostars.org/p/298013/)). This also has the benefit of ensuring that comparisons are made between the same regions, regardless of potential rearrangements. The reads for common tern were trimmed, aligned, and duplicates marked in the same manner as the fairy tern population short-reads. Would only work when using a file as input to command as below.  
```
ls ${DIR}bSteHir1_markdup_autosomes.bam > ${ANGSD}CT.list

angsd -P 16 -doFasta 1 -doCounts 1 -out ${ANGSD}bSteHir1_ancestral -b ${ANGSD}CT.list
```
## Inbreeding Estimates
### F statistics
To estimate relative levels of inbreeding and account for the high likelihood that tara iti do not conform to Hardy-Weinberg Equilibrium (HWE), we first generated population specific genotype likelihoods for input into [ngsF](https://github.com/mfumagalli/ngsTools) v1.2.0-STD. The outputs of these analyses can also be used as a prior for populations that are not in hardy-weinburg equilibrium (HWE).  
```
for POP in AU TI KI
    do
    if [[ "${POP}" == "AU" ]]
        then
        angsd -P 26 -b AU.list -ref $REF -anc $ANC -sites $SITES -out inbreeding/AU \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
            -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
    elif [[ "${POP}" == "TI" ]]
        then
        angsd -P 26 -b TI.list -ref $REF -anc $ANC -sites $SITES -out inbreeding/TI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 12 -setMinDepth 90 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
    else
        angsd -P 26 KI.list -ref $KREF -anc $KANC -sites $SITES -out inbreeding/KI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 12 -setMinDepth 90 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
    fi
done
```
Next, `ngsF` was used to estimate inbreeding. First, an initial search was performed. The below was perfomed for each of the `AU`, `TI`, and `GLOBAL` data sets upon the completion of ANGSD.  
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
Pending XXX
```
### Relatedness
We estimated relatedness, specifically the D<sub>xy</sub> metric, among individuals with [ngsRelate](https://github.com/ANGSD/NgsRelate).  
```
for POP in AU TI
    do
    zcat inbreeding/${POP}.mafs.gz | cut -f 6 | sed 1d > inbreeding/${POP}_freqs
    NSITES=$(zcat inbreeding ${POP}.mafs.gz | tail -n +2 | wc -l)
    if [[ "$POP" == "AU" ]]
        then
        ngsRelate -g inbreeding/${POP}.glf -n 19 -f inbreeding/${POP}_freqs -O inbreeding/${POP}_relatedness
        else
        ngsRelate -g inbreeding/${POP}.glf -n 15 -f inbreeding/${POP}_freqs -O inbreeding/${POP}_relatedness
    fi
done
```

### Runs of Homozygosity
Two methods were used to estimate runs of homozygosity, [ROHAN](https://github.com/grenaud/rohan) and an approach using ANGSD.  

For ROHAN (and SV discovery), all PCR-duplicates were removed as recommended.  
```
samtools markdup -@16 -r --write-index *_markdup_autosomes.bam *_rohan.bam
```
And the Ts/Tv ratio was estimated as a prior for ROHan with VCFtools v0.1.15 using the BCF outputs from ANGSD.
```
for POP in AU TI
    do
    if [[ "$POP" == "AU" ]]
        then
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $REF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 10 --ignore-RG 0
        vcftools --bcf ${ANGSD}samtools/genotypes/${POP}_genotypes.bcf \
            --TsTv-summary --out ${ANGSD}samtools/genotypes/${POP}_genotypes
        else
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $REF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 90 -setMaxDepth 280 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 10 --ignore-RG 0
        vcftools --bcf ${ANGSD}samtools/genotypes/${POP}_genotypes.bcf \
            --TsTv-summary --out ${ANGSD}samtools/genotypes/${POP}_genotypes
    fi
done
```
This resulted in a TsTv ratio of 1.871 for tara iti and 2.507 for Australian fairy tern.  

Finally, ROHan was run on for a subset of individuals. The `--rohmu` flag was varied to ensure regions containing ROHs were detected as per [this discussion](https://github.com/grenaud/ROHan/issues/12#issuecomment-1935539239). The value of 5 x 10<sup>-5</sup> and a window size of 50kb equates 2.5 heterozygous genotypes within this window. This is a potentially stringent cutoff. For comparisions with similar estimates in threatened rattlesnake populations (*Sistrurus catenatus* and *S. tergeminus*) by [Ochoa et al. 2021](https://doi.org/10.1111/mec.16147), we also performed the analyses below with a value of 5 x 10<sup>-4</sup>.  
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
        rohan -t 16 --tstv 2.507 --size 50000 --rohmu 5e-5 -o output/${BASE} $REF $BAM
        else
        rohan -t 16 --tstv 1.871 --size 50000 --rohmu 5e-5 -o output/${BASE} $REF $BAM
    fi
    printf "FINISHED RUNNING ROHAN FOR ${BASE} AT "
    date
done
```
We then prepared summary plots for ROHs, one file each for tara iti and Australian fairy tern containing all ROH sizes calculated using the mid. estimates of heterozygosity.  
```
for FILE in rohmu_5e5/AU*_subset.mid.hmmrohl.gz
    do
    BASE=$(basename $FILE _subset.mid.hmmrohl.gz)
    zcat $FILE | tail -n+2 | cut -f 5 | awk -v SAMP="$BASE" '{ if ( $1 <= 200000 ) len="short_ROH";
        else if ( $1 > 200000 && $1 <= 700000 ) len ="medium_ROH";
        else len ="long_ROH";
        print SAMP"\t"$1"\t"len"\tAU"; }' >> ROHs.tsv
done

for FILE in rohmu_5e5/{SND,SP,TI}*_subset.mid.hmmrohl.gz
    do
    BASE=$(basename $FILE _subset.mid.hmmrohl.gz)
    zcat $FILE | tail -n+2 | cut -f 5 | awk -v SAMP="$BASE" '{ if ( $1 <= 200000 ) len="short_ROH";
        else if ( $1 > 200000 && $1 <= 700000 ) len ="medium_ROH";
        else len ="long_ROH";
        print SAMP"\t"$1"\t"len"\tNZ"; }' >> ROHs.tsv
done
```
We also parsed a file with subsetting the summary files of all individuals.  
```
printf "Sample\tPerc_in_ROH\tMean_ROH_size\tPopulation\n" > ROH_summary.txt

for FILE in rohan/rohmu_5e5/*.summary.txt
    do
    BASE=$(basename $FILE .summary.txt)
    PERC=$(tail -n4 $FILE | head -n1 | awk '{print $5}')
    SIZE=$(tail -n1 $FILE | awk '{print $6}')
    if [[ "$BASE" == "AU"* ]]
        then
        printf "${BASE}\t${PERC}\t${SIZE}\tAU\n" >> ROH_summary.txt
        else
        printf "${BASE}\t${PERC}\t${SIZE}\tNZ\n" >> ROH_summary.txt
    fi
done
```

### Global Heterozygosity
Here, we implemented a global (genome-wide heterozygosity) method from ANGSD. Essentially, this estimate is a proportion of heterozygous genotypes / genome size (excluding regions of the genome with low confidence). Unlike other runs of ANGSD, individual BAMs are used to estimate hetereozygosity, which is simply second value in the SFS/AFS.  
```
for SAMP in ${BAM}*_markdup_autosomes.bam
    do
    BASE=$(basename $SAMP _markdup_autosomes.bam)
    angsd -i $SAMP -anc $ANC -ref $REF -out ${ANGSD}heterozygosity/${BASE} -dosaf 1 -GL 1 -doCounts 1
    realSFS ${ANGSD}heterozygosity/${BASE}.saf.idx > heterozygosity/${BASE}_est.ml
done
```
Once the SFS was estimated for each individual, the number of sites was estimated from the sum of all scaffold sizes included in the bam file and output to a file.
```
SIZE=$(awk '{sum+=$3}; END {print sum}' Katie_autosomes2.bed)
echo $SIZE

NEUTRAL=$(awk '{print $3 - $2}' Katie_neutral_sites.bed | awk '{sum+=$1}; END {print sum}')

printf "Sample\tHeterozygosity\tTool\tPopulation\n" > ${ANGSD}individual_het.tsv 

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
```

## Population Structure
### PCA
To visualise population structure using a PCA, we first ran ANGSD using the `-doGeno 32` and `-doPost 1` options.  
```
angsd -P 16 -b ${ANGSD}GLOBAL.list -ref $REF -out ${ANGSD}samtools/structure_PCA/GLOBAL \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 35 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
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
    -minMapQ 20 -minQ 20 -minInd 35 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
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

## Summary Statistics
### Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS. Here, we are using the common tern as the ancestral state.  
```
ANC=${ANGSD}bSteHir1_ancestral.fasta
region=TI_scaffolded_neutral_regions.bed

angsd -P 16 -b ${ANGSD}AU.list -ref $REF -anc $ANC -out ${ANGSD}samtools/sfs/AU \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
    -GL 1 -doSaf 1

angsd -P 16 -b ${ANGSD}AU.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}neutral/sfs/AU \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 170 -setMaxDepth 360 -doCounts 1 \
    -GL 1 -doSaf 1

angsd -P 16 -b ${ANGSD}TI.list -ref $REF -anc $ANC -out ${ANGSD}samtools/sfs/TI \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 100 -setMaxDepth 280 -doCounts 1 \
    -GL 1 -doSaf 1

angsd -P 16 -b ${ANGSD}TI.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}neutral/sfs/TI \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 100 -setMaxDepth 280 -doCounts 1 \
    -GL 1 -doSaf 1
```
First, bootstrapped estimates were generated under resampling for confidence intervals of nucleotide diversity and F-statistics.
```
realSFS -P 40 -bootstrap 100 -resample_chr 1 -tole 1e-6 ${ANGSD}${TOOL}/sfs/AU.saf.idx > ${ANGSD}${TOOL}/sfs/AU_100boots_resample.sfs
realSFS -P 40 -bootstrap 100 -resample_chr 1 -tole 1e-6 ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/TI_100boots_resample.sfs
realSFS -P 40 -bootstrap 100 -resample_chr 1 -tole 1e-6 ${ANGSD}${TOOL}/sfs/AU.saf.idx ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/GLOBAL_100boots_resample.sfs
```
Then a standard run of `realSFS` was performed for comparison as there remains some discussion around the performance of relative differences in the resampling methods for bootstrapping ([see here](https://github.com/ANGSD/angsd/pull/195)).  
```
realSFS -P 40 ${ANGSD}${TOOL}/sfs/AU.saf.idx ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/GLOBAL_no_boots.sfs
```
## Theta Statistics
Using the ANGSD `-doSaf 1` and realSFS outputs, theta neutrality statistics and nucleotide diversity were estimated using realSFS and associated thetaStat programmes.  
```
for POP in AU TI
    do
    i=1
    while read -r line
    do
        echo $i
        echo $line > ${ANGSD}samtools/sfs/reps/${POP}_rep${i}.sfs
        i=$(( i + 1 ))
    done < ${ANGSD}samtools/sfs/${POP}_100boots.sfs
    i=1
done

for SFS in ${ANGSD}samtools/sfs/reps/*_rep{1..100}.sfs 
    do
    BASE=$(basename $SFS .sfs)
    POP=$(echo "${BASE%_rep*}")
    echo $BASE
    if [[ "$POP" == "AU" ]]
        then
        realSFS saf2theta ${ANGSD}${TOOL}/sfs/AU.saf.idx -sfs ${SFS} -outname ${ANGSD}${TOOL}/diversity/${BASE}
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${BASE}.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}/diversity/${BASE}_thetas_10KBwindows
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${BASE}.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}/diversity/${BASE}_thetas_50KBwindows
    elif [[ "$POP" == "TI" ]]
        realSFS saf2theta ${ANGSD}${TOOL}/sfs/TI.saf.idx -sfs $SFS -outname ${ANGSD}${TOOL}/diversity/${BASE}
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${BASE}.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}/diversity/${BASE}_thetas_10KBwindows
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${BASE}.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}/diversity/${BASE}_thetas_50KBwindows
    else
        realSFS fst index ${DIR}${TOOL}/sfs/AU.saf.idx ${DIR}${TOOL}/sfs/TI.saf.idx -sfs ${SFS} -fstout ${DIR}${TOOL}/distance/${BASE} -whichFst 1
        realSFS fst stats2 ${DIR}${TOOL}/distance/${BASE}.fst.idx -win 10000 -step 1000 > ${DIR}${TOOL}/distance/${BASE}_fst_10KBwindows.txt
        realSFS fst stats2 ${DIR}${TOOL}/distance/${BASE}.fst.idx -win 50000 -step 10000 > ${DIR}${TOOL}/distance/${BASE}_fst_50KBwindows.txt
    fi
done
```
Genetic differentiation.  
```
realSFS fst index sfs/AU.saf.idx sfs/TI.saf.idx -sfs sfs/GLOBAL_no_boots.sfs -fstout distance/GLOBAL_fst

realSFS fst stats distance/GLOBAL_fst.idx \
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
## Putative Genetic Load
Given that the tara iti reference genome could not be annotated with RNA sequencing, ensure a robust and conservative assessment of genetic load that is translatable across species comparisons, we limited load estimates to highly conserved BUSCO genes in the fairy tern species complex (*Sterna nereis* spp.) and kakī (*Himantopus novazealandiae*) for comparison.  

For the fairy tern and kakī reference genomes, we concatenated single copy BUSCO sequences into species specific `GFF` files.
```
cat Katies_genome/ragtag_busco/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.gff > Katies_genome/ragtag_busco/merged_single_copy_busco_sequences.gff
cat kaki_genome/busco_output/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.gff > kaki_genome/busco_output/merged_single_copy_busco_sequences.gff
```
We then sorted and each file and converted them to `BED` format with [BEDtools](https://bedtools.readthedocs.io/en/latest/)v2.30.  
```
bedtools sort -i Katies_genome/ragtag_busco/merged_single_copy_busco_sequences.gff > angsd/single_copy_BUSCO.bed
bedtools merge -i angsd/single_copy_BUSCO.bed > angsd/single_copy_BUSCO.merged.bed
```
We then filtered the `BED` file to excluded unplaced scaffolds and sex chromosomes (consistant with above population analyses) and indexed this `BED` file with ANGSD.
```
BUSCO=${ANGSD}single_copy_BUSCO.merged.bed

for POP in AU TI
    do
    angsd -b CEBO.filelist -sites $BUSCO -out ${ANGSD}load/${POP} \
        -minMapQ 20 -minQ 20 -gl 1 -domajorminor 1 -snp_pval 1e-6 -domaf 1 \
        -doGlf 3 -doGeno 4 -doPost 1 -postCutoff 0.95
done
```