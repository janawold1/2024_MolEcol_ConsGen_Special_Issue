# Population genomics with ANGSD
[ANGSDv0.935](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) and ngsTools were used to estimate summary statistics for Australian fairy tern (*Sternual nereis nereis*) sampled from Western Australia and tara iti (*S. nereis davisae*) from Northland, NZ. The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. PCR duplicates were marked and sex chromosomes were excluded for SNP-based population analyses.  
```
ls ${DIR}{AU,SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}GLOBAL.list
ls ${DIR}{SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}TI.list
ls ${DIR}AU*_markdup_autosomes.bam > ${ANGSD}AU.list
ls ${DIR}H0*_markdup_autosomes.bam > ${ANGSD}KI.list

TREF=${DIR}SP01_genome/SP01_5kb_ragtag.fa
KREF=${DIR}kaki_genome/himNova-hic-scaff.fa

for POP in GLOBAL AU TI KI_10x
    do
    if (( $POP == KI_10x ))
    then
        angsd -P 32 -b ${ANGSD}${POP}.list -ref $KREF -out ${ANGSD}qc/${POP}.qc \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 800
    else
        angsd -P 32 -b ${ANGSD}${POP}.list -ref $REF -out ${ANGSD}qc/${POP}.qc \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 800
    fi
done

angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -out ${ANGSD}qc/KI.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 2000
```
Courtesy of scripts provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.

|        Population        |Minimum MapQ|Minimum Q|Minimum Depth|Maximum Depth|Number of Individuals|
| ------------------------ | ---------- | ------- | ----------- | ----------- | ------------------- |
|Australian Fairy Tern (WA)|     20     |    20   |     200     |     350     |         19          |
|         Tara iti         |     20     |    20   |     120     |     280     |         15          |
|     Global fairy tern    |     20     |    20   |     272     |     630     |         34          |
|         Kakī 10x         |     20     |    20   |     192     |     300     |         24          |
|         Kakī 50x         |     20     |    20   |     700     |    1,200    |         24          |

 Below is code run for analyses performed on alignments of fairy terns to the tara iti reference scaffolded using the common tern assembly (see [00_genome_assembly.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/00_genome_assembly.md)) and alignments of kakī to the high-quality reference assembly for this species (see [here](https://www.genomics-aotearoa.org.nz/our-work/completed-projects/high-quality-genomes) for details). All analyses for each group were performed in a similar manner for downstream comparisons. Quality thresholds were adjusted as per the table above.  

In addition, scaffolds that adhered to expected population sequence coverage (e.g., ~10x) and likely represented autosomal chromosomes were extracted from the reference fasta and individual BAM alignments before identifying putative coding regions as per below.  

## Excluding putative coding Regions
Neither the tara iti or kakī reference assemblies have transcriptomes. To assess the potential consequences of the inclusion of regions under selection, and to have a 'first look' at the genomic diversity around coding regions, we performed *ab initio* gene prediction using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus/tree/master) v3.5.0.  
```
augustus --sample=100 --alternatives-from-sampling=true --strand=both \
        --gff3=yes --progress=true --temperature=3 --species=chicken \
        --protein=on $TREF > ${DIR}augustus/TI_AUGUSTUS.gff

augustus --sample=100 --alternatives-from-sampling=true --strand=both \
        --gff3=yes --progress=true --temperature=3 --species=chicken \
        --protein=on $KREF > ${DIR}augustus/kaki_AUGUSTUS.gff
```
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.30.0 was then used to sort, and merge these putative gene regions. An additional 1kb of sequence on either side of annotations were included to reduce linkage. Duplicate regions were then merged with `bedtools merge`.  
```
grep -v "#" reference/gene_predictions/SP01_AUGUSTUS.gff | \
    bedtools sort -i - > angsd/augustus_autosomal_predictions.bed

bedtools merge -i angsd/augustus_autosomal_predictions.bed > angsd/augustus_autosomal_predictions_merged.bed

grep -v "#" kaki_genome/kaki_AUGUSTUS.gff | \
    grep -v scaffold_4 | \
    bedtools sort -i - > angsd/kaki_augustus_predictions.bed
```
 To account for linkage, the window size for these regions increased by 1kb on either side. The file was adusted in cases where the addition of a 1kb buffer extended beyond the start (n = 11) or end of the chromosome (n = 7).
```
awk '{print $1"\t"$2-1000"\t"$3+1000}' angsd/augustus_autosomal_predictions_merged.bed > angsd/augustus_autosomal_predictions_merged_add1kb.bed
awk '{print $1"\t"$2-1000"\t"$3+1000}' angsd/kaki_augustus_predictions_merged.bed > angsd/kaki_augustus_predictions_merged_add1kb.bed
```
To define putatively neutral sites for analyses, we extracted the complement regions as below.  
```
awk -v OFS='\t' {'print $1,$2'} ${TREF}$.fai > reference/SP01_autosome_lengths.txt

bedtools complement \
    -i angsd/augustus_autosomal_predictions_merged_add1kb.bed \
    -g reference/SP01_autosome_lengths.txt > angsd/SP01_neutral_sites.bed
```
For the fairy tern, this left 27,337 regions covering 562,663,009 bp for analyses, which roughly corresponds to 47% of the genome.  

This file was then indexed for ANGSD with `angsd sites`.  
```
angsd sites index TI_scaffolded_neutral_regions.bed
```
## Diversity Indices
### Runs of Homozygosity 
[ROHAN](https://github.com/grenaud/rohan)was used to call runs of homozygosity (ROHs). Following the recommendations, we removed all PCR-duplicates from BAM alignment files, estimated a transition/transversion ratio with [VCFtools](https://vcftools.github.io/index.html) v0.1.15, and defined an allowable level of heterozygosity within a 50kb window.  

Starting with SAMtools, we removed PCR duplicates from our markdup files.  
```
samtools markdup -@16 -r --write-index *_markdup_autosomes.bam *_rohan.bam
```
And the Ts/Tv ratio was estimated as a prior for ROHan with VCFtools using the BCF outputs from ANGSD.
```
for POP in AU TI
    do
    if [[ "$POP" == "AU" ]]
        then
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $TREF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6 -doGeno 10 --ignore-RG 0
        vcftools --bcf ${ANGSD}samtools/genotypes/${POP}_genotypes.bcf \
            --TsTv-summary --out ${ANGSD}samtools/genotypes/${POP}_genotypes
        elif [[ "$POP" == "TI" ]]
        then
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $TREF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6 -doGeno 10 --ignore-RG 0
        vcftools --bcf ${ANGSD}samtools/genotypes/${POP}_genotypes.bcf \
            --TsTv-summary --out ${ANGSD}samtools/genotypes/${POP}_genotypes
        else
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $KREF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 10 --ignore-RG 0
        vcftools --bcf ${ANGSD}samtools/genotypes/${POP}_genotypes.bcf \
            --TsTv-summary --out ${ANGSD}samtools/genotypes/${POP}_genotypes
    fi
done
```
This resulted in a TsTv ratio of 1.871 for tara iti, 2.507 for Australian fairy tern, and 2.725 for kakī.  

Finally, ROHan was run with the `--rohmu` flag varied to ensure regions containing ROHs were detected as per [this discussion](https://github.com/grenaud/ROHan/issues/12#issuecomment-1935539239). The value of 5 x 10<sup>-5</sup> and a window size of 50kb equates 2.5 heterozygous genotypes within this window. For kakī, this threshold was left at the default setting of 1 x 10<sup>-5</sup> and equated to 0.5 heterozygous sites within a 50kb window. This much lower ROH threshold for kakī is attributable to the much higher sequence depth for the kakī data (i.e, target 50x vs 10x sequence depth).  
```
./faSomeRecords reference/SP01_5kb_ragtag.fa SP01_autosomes2.bed reference/SP01_ragtag_autosomes2.fa
./faSomeRecords reference/himNova-hic-scaff.fa himNova_autsosomes.bed reference/himNova_autosomes.fa

samtools faidx reference/SP01_ragtag_autosomes2.fa
samtools faidx reference/himNova_autosomes.fa

REF=reference/SP01_ragtag_autosomes2.fa
KIREF=reference/himNova_autosomes.fa

for BAM in *_rohan.bam
    do
    BASE=$(basename $bam _rohan.bam)
    printf "STARTED RUNNING ROHAN FOR ${BASE} AT "
    date
    if [[ "$BASE" == "AU"*]]
        then
        rohan -t 16 --tstv 2.507 --size 50000 --rohmu 5e-5 -o output/${BASE} $TREF $BAM
    elif [[ "$BASE" == "H0"*]]
        then
        rohan -t 16 --tstv 2.725 --size 50000 -o output/${BASE} $KREF $BAM
    else
        rohan -t 16 --tstv 1.871 --size 50000 --rohmu 5e-5 -o output/${BASE} $TREF $BAM
    fi
    printf "FINISHED RUNNING ROHAN FOR ${BASE} AT "
    date
done
```
We then prepared summary plots for ROHs, one file each for tara iti and Australian fairy tern containing all ROH sizes calculated using the mid. estimates of heterozygosity.  
```
for FILE in fairy_tern/rohmu_5e5/AU*_subset.mid.hmmrohl.gz
    do
    BASE=$(basename $FILE _subset.mid.hmmrohl.gz)
    zcat $FILE | tail -n+2 | cut -f 2-5 | awk -v SAMP="$BASE" '{ if ( $4 <= 200000 ) len="Short ROH";
        else if ( $4 > 200000 && $4 <= 700000 ) len ="Medium ROH";
        else len ="Long ROH";
        print SAMP"\t"$0"\t"len"\tAU"; }' >> ROHs.tsv
done

for FILE in fairy_tern/rohmu_5e5/{SND,SP,TI}*_subset.mid.hmmrohl.gz
    do
    BASE=$(basename $FILE _subset.mid.hmmrohl.gz)
    zcat $FILE | tail -n+2 | cut -f 2-5 | awk -v SAMP="$BASE" '{ if ( $4 <= 200000 ) len="Short ROH";
        else if ( $4 > 200000 && $4 <= 700000 ) len ="Medium ROH";
        else len ="Long ROH";
        print SAMP"\t"$0"\t"len"\tNZ"; }' >> ROHs.tsv
done

for FILE in kaki/rohmu_default/H0*.mid.hmmrohl.gz
    do
    BASE=$(basename $FILE .mid.hmmrohl.gz)
    zcat $FILE | tail -n+2 | cut -f 2-5 | awk -v SAMP="$BASE" '{ if ( $4 <= 200000 ) len="Short ROH";
        else if ( $4 > 200000 && $4 <= 700000 ) len ="Medium ROH";
        else len ="Long ROH";
        print SAMP"\t"$0"\t"len"\tKI"; }' >> ROHs.tsv
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
    elif [[ "$BASE" == "H0"* ]]
        then
        printf "${BASE}\t${PERC}\t${SIZE}\tKI\n" >> ROH_summary.txt
    else
        printf "${BASE}\t${PERC}\t${SIZE}\tNZ\n" >> ROH_summary.txt
    fi
done
```
### Individual Heterozygosity
Here, we implemented a global (genome-wide heterozygosity) method from ANGSD. Essentially, this estimate is a proportion of heterozygous genotypes / genome size (excluding regions of the genome with low confidence). Unlike other runs of ANGSD, individual BAMs are used to estimate hetereozygosity, which is simply second value in the SFS/AFS.  
```
SAMP=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /nesi/nobackup/uc03718/angsd/GLOBAL.list)
NAME=$(basename $SAMP _autosomes_nodup.bam)

printf "\nSTARTED RUNNING ANGSD TO ESTIMATE HETEROZYGOSITY FOR $NAME AT "
date

angsd -P 32 -i ${SAMP} -ref $TREF -out ${DIR}samtools/heterozygosity/${NAME} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -angbaq 1 \
        -minMapQ 20 -minQ 20 -doCounts 1 -dosaf 1 -GL 1

realSFS ${DIR}samtools/heterozygosity/${NAME}.saf.ids > ${DIR}samtools/heterozygosity/${NAME}_est.ml
```
Once the SFS was estimated for each individual, the number of sites was estimated from the sum of all scaffold sizes included in the bam file and output to a file.
```
printf "Sample\tHeterozygosity\tTool\tPopulation\n" > ${ANGSD}individual_het.tsv 

for SAMP in ${ANGSD}${TOOL}/heterozygosity/*_est.ml
    do
    BASE=$(basename $SAMP _fold_est.ml)
    TOT=$(awk '{print $1 + $2 + $3}' $SAMP)
    HET=$(awk -v var=$TOT '{print $2/var}' $SAMP)
    if [[ "$BASE" == "AU"* ]]
        then
        printf "$BASE\t$HET\t$TOOL\tAU\n" >> ${ANGSD}individual_het.tsv
    elif [[ "$BASE" == "H0"* ]]
        then
        printf "$BASE\t$HET\t$TOOL\tKI\n" >> ${ANGSD}individual_het.tsv
    else
        printf "$BASE\t$HET\t$TOOL\tNZ\n" >> ${ANGSD}individual_het.tsv
    fi
done
```
### Local Heterozygosity 
```
for POP in AU TI KI
    do
    if [[ "$POP" == "AU" ]]
        then
        angsd -bam ${ANGSD}${POP}.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/heterozygosity/${POP}_perSite \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doHWE 1 -domajorminor 1 -doMaf 1 -SNP_pval 1e-6
        angsd -bam ${ANGSD}${POP}.list -ref $TREF -anc $TANC -sites $TSITES -out ${ANGSD}samtools/heterozygosity/${POP}_perSite \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doHWE 1 -domajorminor 1 -doMaf 1 -SNP_pval 1e-6
    elif [[ "$POP" == "TI" ]]
        then
        angsd -bam ${ANGSD}${POP}.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/heterozygosity/${POP}_perSite \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doHWE 1 -domajorminor 1 -doMaf 1 -SNP_pval 1e-6
        angsd -bam ${ANGSD}${POP}.list -ref $TREF -anc $TANC -sites $TSITES -out ${ANGSD}samtools/heterozygosity/${POP}_perSite \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doHWE 1 -domajorminor 1 -doMaf 1 -SNP_pval 1e-6
    else
        angsd -bam ${ANGSD}${POP}.list -ref $KREF -anc $KANC -out ${ANGSD}samtools/heterozygosity/${POP}_perSite \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doHWE 1 -domajorminor 1 -doMaf 1 -SNP_pval 1e-6
        angsd -bam ${ANGSD}${POP}.list -ref $KREF -anc $KANC -sites $KSITES -out ${ANGSD}neutral/heterozygosity/${POP}_perSite \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doHWE 1 -domajorminor 1 -doMaf 1 -SNP_pval 1e-6
    fi
done
```

### Relatedness
We estimated relatedness between our samples as per:
```
angsd -P 32 -b ${ANGSD}${POP}.list -ref ${REF} -out ${ANGSD}samtools/relatedness/${POP} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 -doMajorMinor 1 -GL 1 -doMaf 1 \
    -skipTriallelic 1 -SNP_pval 1e-6 -doGlf 3

zcat ${ANGSD}samtools/relatedness/${POP}.mafs.gz | cut -f6 | sed 1d > ${ANGSD}samtools/relatedness/${POP}_whole-genome_freq

ngsrelate -g ${ANGSD}samtools/relatedness/${POP}.glf.gz \
    -n 15 -f ${ANGSD}samtools/relatedness/${POP}_whole-genome_freq \
    -O ${ANGSD}samtools/relatedness/${POP}_whole-genome_relatedness
```

## Population Structure
### MDS
To construct a MDS for fairy tern populations, we first ran ANGSD as below. Notably, here we add the `-doGeno 8` flag.  
```
angsd -P 26 -b GLOBAL.list -ref $REF -anc $ANC -out ${ANGSD}structure_MDS/GLOBAL \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 34 -setMinDepth 272 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1
```
This command not only generates the required `*.mafs.gz`, but also the `*.geno.gz`, which contains posterior probabilities of all possible genotypes required for estimating genetic distance with [ngsDistv1.0.10](https://github.com/mfumagalli/ngsTools). In addition, a `pops.label` file denoting the population of origin (one entry on a new line for each samples) is necessary for estimating genetic distance.  
```
NSITES=$(zcat GLOBAL.mafs.gz | tail -n +2 | wc -l)
echo $NSITES

cat GLOBAL.list | sed 's%/path/to/bams/%%g' | sed 's/_markdup_autosomes.bam//g' > ${ANGSD}mds_list

ngsDist -verbose 1 -geno ${ANGSD}structure_MDS/GLOBAL.geno.gz -probs \
    -n_ind 34 -n_sites $NSITES -labels ${ANGSD}mds_list -o ${ANGSD}structure_MDS/GLOBAL_dist
```
Extract and construct MDS.  
```
tail -n +3 ${ANGSD}structure_MDS/GLOBAL_dist | Rscript --vanilla --slave getMDS.R \
    --no_header --data_symm -n 4 -m "mds" -o ${ANGSD}structure_MDS/GLOBAL.mds

head structure_MDS/GLOBAL.mds
```
This MDS output is plotted in [06_visualisations.ipynb].  
## Summary Statistics
### Site Frequency Spectrum
#### Inferring Ancestral alleles
It is important to have the ancestral state for some of the analyses below. To this end, we used short reads sequenced by the [B10K](https://b10k.genomics.cn/species.html) consortium and stored by the China National GeneBank DataBase (CNGBdb) for large-bill tern ([*Phaetusa simplex*](https://db.cngb.org/search/organism/297813/) NCBI ID: 297813, CNGB Sample ID: CNS0103161, CNGB Experiment ID: CNX0082959) and black skimmer ([*Rhynchops niger*](https://db.cngb.org/search/organism/227184/) NCBI ID: 227184, CNGB Sample ID: CNS0103105, CNGB Experiment ID: CNX0082867 - CNX0082869) to polarise the site frequency spectrum for fairy terns. Whereas resequence data generated by [Galla et al 2019](https://doi.org/10.3390/genes10010009) for the avocet (*Recurvirostra avosetta*) and pied stilt (*Himantopus himantopus*) were used to polarise kakī. This ([brief discussion](https://www.biostars.org/p/298013/)) outlines the benefit of using this approach for ANGSD-base analyses. The reads for each of the outgroups were trimmed, aligned, and duplicates marked in the same manner as the fairy tern and kakī population short-reads.  
```
angsd -P 16 -doFasta 2 -doCounts 1 -out ${ANGSD}tern_ancestral -i ${DIR}merged_tern_markdup_autosomes.bam
angsd -P 16 -doFasta 2 -doCounts 1 -out ${ANGSD}stilt_ancestral -i ${DIR}merged_stilt_markdup_autosomes.bam
angsd -P 16 -doFasta 2 -doCounts 1 -out ${ANGSD}pied_alleles -b pied.list
angsd -P 16 -doFasta 2 -doCounts 1 -out ${ANGSD}AFT_alleles -b AU.list
```
#### SFS Estimation
We then estimated the SFS using ANGSD and realSFS.
```
for TOOL in samtools neutral
    do
    printf "STARTED RUNNING ANGSD FOR $TOOL DATA SET AT "
    date
    if [[ "$TOOL" == "samtools" ]]
        then
        angsd -P 16 -b ${ANGSD}AU.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/sfs/AU \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}TI.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/sfs/TI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -anc $KANC -out ${ANGSD}samtools/sfs/KI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doSaf 1
    else
        angsd -P 16 -b ${ANGSD}AU.list -ref $TREF -anc $TANC -sites $TSITES -out ${ANGSD}neutral/sfs/AU \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}TI.list -ref $TREF -anc $TANC -sites $TSITES -out ${ANGSD}neutral/sfs/TI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -anc $KANC -sites $KSITES -out ${ANGSD}neutral/sfs/KI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doSaf 1
    fi
done
```
Then `realSFS` was run to estimate the site frequency spectrum and theta-based statistics.  
```
realSFS -P 40 -anc ${TANC} -ref ${TREF} ${ANGSD}${TOOL}/sfs/AU.saf.idx > ${ANGSD}${TOOL}/sfs/AU.sfs
realSFS -P 40 -anc ${TANC} -ref ${TREF} ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/TI.sfs
realSFS -P 40 -anc ${KANC} -ref ${KREF} ${ANGSD}${TOOL}sfs/KI.saf.idx > ${ANGSD}${TOOL}sfs/KI.sfs
realSFS -P 40 -anc ${TANC} -ref ${TREF} ${ANGSD}${TOOL}/sfs/AU.saf.idx ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/GLOBAL.sfs
```

### Estimating Population Differentiation (Watterhson's Theta)
Here we estimate Theta (F<sub>ST</sub>) between Australian fairy tern populations in Western Australia and tara iti. To do so, we leverage the site frequency spectrum (SFS) generated above.  
```
realSFS fst index ${DIR}${TOOL}/sfs/AU.saf.idx ${DIR}${TOOL}/sfs/TI.saf.idx -sfs ${DIR}${TOOL}/sfs/GLOBAL.sfs -fstout ${DIR}${TOOL}/distance/GLOBAL -whichFst 1
realSFS fst stats2 ${DIR}${TOOL}/distance/GLOBAL.fst.idx \
    -tole 1e-6 -ref ${TREF} -anc ${TANC} \
    -win 10000 -step 1000 > ${DIR}${TOOL}/distance/GLOBAL_stat2_10KBwindow_1KBstep.tsv

realSFS fst stats2 ${DIR}${TOOL}/distance/GLOBAL.fst.idx \
    -tole 1e-6 -ref ${TREF} -anc ${TANC} \
    -win 50000 -step 10000 > ${DIR}${TOOL}/distance/GLOBAL_stat2_50KBwindow_10KBstep.tsv

realSFS fst stats ${DIR}${TOOL}/distance/GLOBAL.fst.idx \
    -tole 1e-6 -ref ${TREF} -anc $ANC \
    -win 50000 -step 10000 -whichFst 1 > ${ANGSD}${TOOL}distance/GLOBAL_stat1_50KBwindow_10KBstep.tsv
```
The final command estimated a global weighted and unweighted F<sub>ST</sub> of `0.001981` and `0.838787` respectively for 259,924,909 putatively neutral sites and an unweighted F<sub>ST</sub> of `0.001987` and weighted F<sub>ST</sub> of `0.836593` for all 507,424,045 sites.  

### Estimating Nucleotide Diversity (π)
We estimated nucleotide diversity (π) for the West Australian population of fairy tern and tara iti independently, using the site frequency spectrum (SFS) estimated above as a prior.  
```
for POP in AU TI KI
    do
    if [[ "$POP" == "AU" ]]
        then
        angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}${TOOL}diversity/${POP}_pest \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doSaf 1 -pest ${ANGSD}${TOOL}/sfs/${POP}.sfs
        realSFS -p 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx > ${ANGSD}${TOOL}diversity/${POP}_pest.sfs
        realSFS saf2theta ${ANGSD}${TOOL}/diversity/${POP}_pest.saf.idx -sfs ${ANGSD}${TOOL}diversity/${POP}_pest.sfs -outname ${ANGSD}${TOOL}/diversity/${POP}_pest
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_10KBwindows_1KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_50KBwindows_10KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_all_Thetas
    elif [[ "$POP" == "TI" ]]
        then
        angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}${TOOL}diversity/${POP}_pest \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doSaf 1 -pest ${ANGSD}${TOOL}/sfs/${POP}.sfs
        realSFS -p 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx > ${ANGSD}${TOOL}diversity/${POP}_pest.sfs
        realSFS saf2theta ${ANGSD}${TOOL}/diversity/${POP}_pest.saf.idx -sfs ${ANGSD}${TOOL}diversity/${POP}_pest.sfs -outname ${ANGSD}${TOOL}/diversity/${POP}_pest
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_10KBwindows_1KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_50KBwindows_10KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_all_Thetas
    else
        angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}${TOOL}diversity/${POP}_pest \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doSaf 1 -pest ${ANGSD}${TOOL}/sfs/${POP}.sfs
        realSFS -p 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx > ${ANGSD}${TOOL}diversity/${POP}_pest.sfs
        realSFS saf2theta ${ANGSD}${TOOL}/diversity/${POP}_pest.saf.idx -sfs ${ANGSD}${TOOL}diversity/${POP}_pest.sfs -outname ${ANGSD}${TOOL}/diversity/${POP}_pest
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_10KBwindows_1KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_50KBwindows_10KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_all_Thetas
done
```

### Estimating D<sub>XY</sub>
When aiming to estimate D<sub>XY</sub>, we need to first get a list of sites variable in each population. This is so sites that may be fixed in one population are included in our estimates.  
```
angsd -P 24 -b ${ANGSD}AU.list -ref $TREF -anc $TANC -sites $TSITES -out ${ANGSD}neutral/diversity/AU_initialDxy \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 -doMajorMinor 1 \
    -GL 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6

angsd -P 24 -b ${ANGSD}TI.list -ref $TREF -anc $TANC -sites $TSITES -out ${ANGSD}neutral/diversity/TI_initialDxy \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 -doMajorMinor 1 \
    -GL 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6

angsd -P 24 -b ${ANGSD}AU.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/diversity/AU_initialDxy \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 -doMajorMinor 1 \
    -GL 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6

angsd -P 24 -b ${ANGSD}TI.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/diversity/TI_initialDxy \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 -doMajorMinor 1 \
    -GL 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6
```
We then extracted the chromosome and position from each population and merge them into a sites file.
```
zcat ${ANGSD}neutral/diversity/{AU,TI}_initialDxy.mafs.gz | awk '{print $1"\t"$2}' | grep -v "chromo" | sort -k1,1 -k2,2n | uniq > ${ANGSD}neutral_Dxy_sites.bed
zcat ${ANGSD}samtools/diversity/*_initialDxy.mafs.gz | awk '{print $1"\t"$2}' | grep -v "chromo" | sort -k1,1 -k2,2n | uniq > ${ANGSD}whole-genome_Dxy_sites.bed

angsd sites index neutral_Dxy_sites.bed
angsd sites index whole-genome_Dxy_sites.bed
```

And finally estimated allele frequency estimates at these overlapping sites.  
```
for POP in AU TI
    do
    angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites ${ANGSD}neutral_DXYsites.bed \
        -out ${ANGSD}neutral/diversity/${POP}_anc_DXY -GL 1 -doMajorMinor 5 -doMaf 1
    angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -sites ${ANGSD}neutral_DXYsites.bed \
        -out ${ANGSD}neutral/diversity/${POP}_ref_DXY -GL 1 -doMajorMinor 4 -doMaf 1
    angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites ${ANGSD}whole-genome_DXYsites.bed \
        -out ${ANGSD}samtools/diversity/${POP}_anc_DXY -GL 1 -doMajorMinor 5 -doMaf 1
    angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -sites ${ANGSD}whole-genome_DXYsites.bed \
        -out ${ANGSD}samtools/diversity/${POP}_ref_DXY -GL 1 -doMajorMinor 4 -doMaf 1
done
```
And finally we adapted a Rscript by Joshua Penalba [here](https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/calcDxy.R) to calculate DXY using the allele frequencies in each respective population's MAF file.

## Putative Genetic Load
Given that the tara iti reference genome could not be annotated with RNA sequencing, ensure a robust and conservative assessment of genetic load that is translatable across species comparisons, we limited load estimates to highly conserved BUSCO genes in the fairy tern species complex (*Sterna nereis* spp.) and kakī (*Himantopus novazealandiae*) for comparison.  

To start, we ran BUSCO v5.4.7 for the tara iti and kakī reference genome.
```
cd SP01_genome
busco --in SP01_5kb_ragtag.fa --out busco/ --mode genome --lineage_dataset aves_odb10 --cpu 32
cd kaki_genome
busco --in himNova-hic-scaff.fa --out busco/ --mode genome --lineage_dataset aves_odb10 --cpu 32
```
We then concatenated single copy BUSCO sequences into species specific `GFF` files.  
```
cat SP01_genome/busco/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.gff > SP01_genome/busco/merged_single_copy_busco_sequences.gff
cat kaki_genome/busco/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.gff > kaki_genome/busco/merged_single_copy_busco_sequences.gff
```
We then sorted and each file and converted them to `BED` format with [BEDtools](https://bedtools.readthedocs.io/en/latest/)v2.30.  
```
bedtools sort -i SP01_genome/busco/merged_single_copy_busco_sequences.gff > SP01_genome/fairy_single_copy_BUSCO.bed
bedtoosl merge -i SP01_genome/fairy_single_copy_BUSCO.bed > SP01_genome/fairy_single_copy_BUSCO.merged.bed

bedtools sort -i kaki_genome/busco/merged_single_copy_busco_sequences.gff > kaki_genome/kaki_single_copy_BUSCO.bed
bedtools merge -i kaki_genome/kaki_single_copy_BUSCO.bed > kaki_genome/kaki_single_copy_BUSCO.merged.bed
```
We then filtered the `BED` file to excluded unplaced scaffolds and sex chromosomes (consistant with above population analyses) and indexed this `BED` file with ANGSD.

### Sites Polarised with ANGSD
To estimate masked and realised load in each fairy tern population and kakī we used ANGSD to output the major minor alleles, the called genotype, the posterior probability of the called genotype and all possible genotypes (`-doGeno 31`) to a BCF file (`-doBcf 1`). 
```
for POP in AU TI
    do
    if [[ "$POP" == "AU" ]]
        then
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/genotypes/${POP}_whole-genome_polarized \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -doPost 1 -postCutoff 0.95 -doBcf 1 -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6 -doGeno 31 --ignore-RG 0
        elif [[ "$POP" == "TI" ]]
        then
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $TREF -anc $TANC -out ${ANGSD}samtools/genotypes/${POP}_whole-genome_polarized \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -doPost 1 -postCutoff 0.95 -doBcf 1 -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6 -doGeno 31 --ignore-RG 0
        else
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $KREF -anc $KANC -out ${ANGSD}samtools/genotypes/${POP}_whole-genome_polarized \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -doPost 1 -postCutoff 0.95 -doBcf 1 -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6 -doGeno 31 --ignore-RG 0
    fi
done
```

### SIFT
We then used [SIFT4G](https://github.com/rvaser/sift4g) to identify deleterious sites within these BUSCO genes. To do so, we leveraged the outputs from BUSCO to create a database called for our tara iti and kakī reference genomes.  

Using the same method as when we concatenated the `GFF` files used to create `BED` files for variant discovery, we also obtained the reference protein sequences for each complete single copy BUSCO gene.
```
cat SP01_genome/busco/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.faa > SP01_genome/busco/merged_single_copy_busco_prot.fa
cat kaki_genome/busco/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.faa > kaki_genome/busco/merged_single_copy_busco_prot.fa
```
The SIFT documentation suggests using `gffread`, however this is deprecated and unmaintained. We opted to convert the `GFF` files output from BUSCO to `GTF` with [agat](https://github.com/NBISweden/AGAT) v. 1.4.0. First we fixed the concatenated `GFF` from BUSCO to be more compatible with SIFT and VEP.  
```
agat_sp_manage_IDs.pl --gff references/kaki_merged_single_copy_busco.gff -o references/kaki_BUSCO_check.gff
agat_sp_manage_IDs.pl --gff references/SP01_single_copy_BUSCO.gff -o references/SP01_BUSCO_check.gff
```
Then sorted and compressed `GFF`.  
```
cat references/kaki_BUSCO_check.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > references/kaki_BUSCO_check.gff.gz
cat references/SP01_BUSCO_check.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > references/SP01_BUSCO_check.gff.gz
```
This `GFF` was then converted to `GTF` with AGAT and compressed for SIFT.  
```
agat_convert_sp_gff2gtf.pl --gff SP01_genome/busco/SP01_BUSCO_check.gff -o SP01_genome/busco/SP01_scBUSCO.gtf
agat_convert_sp_gff2gtf.pl --gff kaki_genome/busco/kaki_BUSCO_check.gff -o kaki_genome/busco/kaki_scBUSCO.gtf

bgzip SP01_genome/busco/SP01_scBUSCO.gtf
cp SP01_genome/busco/SP01_scBUSCO.gtf.gz fairy_SIFT_databases/gene-annotation-src/

bgzip kaki_genome/busco/kaki_scBUSCO.gtf
cp kaki_genome/busco/kaki_scBUSCO.gtf.gz kaki_SIFT_databases/gene-annotation-src/
```
The tara iti reference genome had each scaffold represented on a single line. This is incompatible with the notation for SIFT programmes. To address this, we used the Unix command `fold` to reduce line widths to 60 characters.  
```
fold -w 60 SP01_5kb_ragtag > SP01_5kb_ragtag_fold.fa
bgzip SP01_5kb_fold.fa
```
We then used [SIFT4G_Create_Genomic_DB](https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB) to construct a protein database from BUSCO genes identified in our tara iti and kakī reference genomes.  

First, we downloaded and configured the docker image.
```
git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git
cd SIFT4G_Create_Genomic_DB
docker build -t sift4g_db .
```
Relevant directories were then created.
```
mkdir -p fairy_SIFT_databases/{gene-annotation-src,chr-src,dbSNP}
```
Then, a config file was formatted as below, adjusting the relevant settings as necessary.  
```
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_COTE_TABLE=2
MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial

PARENT_DIR=fairy_SIFT_databases
ORG=sterna_nereis_davisae
ORG_VERSION=SP01_v1.0
DBSNP_VCF_FILE=TI_calls.vcf.gz

# Running SIFT4G, this path works for the Dockerfile
SIFT4G_PATH=/home/jana/sift4g/bin/sift4g

# POTEIN_DB needs to be uncompressed
PROTEIN_DB=fairy_SIFT_databases/uniref90.fasta

# Subdirectories, don't need to change
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=log.txt
ZLOGFILE=log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP

# Doesn't need to change
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord
```
SIFT databases were then constructed.
```
sudo docker run -it --user $(id -u):$(id -g) -v /home/jana/:/home/jana/ sift4g_db /bin/bash
perl make-SIF-db-all.pl -c fairy_terns.txt
exit
```
And finally we annotated the filtered VCF file using our SIFT database.  
```
java -jar /home/jana/SIFT_Annotator.jar -c -i ${DIR}GLOBAL_polarized_filtered.vcf \
    -d fairy_SIFT_databases/SP01_v1/ \
    -r SIFT/fairy_output -t
```

### Variant Effect Predictor
VEP was run below using the GFF file constructed above.  
```
perl vep -i GLOBAL_whole-genome_polarized.vcf \
    --custom references/SP01_BUSCO_check.gff.gz,FAIRY_GFF,gff \
    --fasta references/SP01_5kb_ragtag_fold60.fa.gz \
    --everything \
    -o vep_data/angsd_high_confidence_BUSCO_SNPs/GLOBAL_whole-genome_polarized

perl vep -i KI_whole-genome_polarized.vcf \
    --custom references/kaki_BUSCO_check.gff.gz,KAKI_GFF,gff \
    --fasta references/him_Nova-hic-scaff.fa \
    --everything \
    -o vep_data/angsd_high_confidence_BUSCO_SNPs/KI_whole-genome_polarized
```
### Intersecting calls between VEP & SIFT
To find sites that both had an impact (VEP) and were deleterious (SIFT), we first filtered the SIFT output for those sites with a SIFT score <= 0.05.
```
awk '{print $1"\t"$2"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' GLOBAL_whole-genome_SIFTannotations.txt | awk '{ if (($5 != "NA" && $6 <= 0.05 && $7 >= 2.75) || ($6 != "NA" && $6 <=0.05 && $7 >= 3.5)) print $0}' > deleterious_SIFT_fairy.tsv
awk '{print $1"\t"$2"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' GLOBAL_whole-genome_SIFTannotations.txt | awk '{ if (($5 != "NA" && $6 > 0.05 && $7 >= 2.75) || ($6 !="NA" && $6 > 0.05 && $7 <= 3.5) print $0}' > tolerant_SIFT_fairy.tsv

awk '{print $1"\t"$2"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' kaki_whole-genome_SIFTannotations.txt | awk '{ if (($5 != "NA" && $6 <= 0.05 && $7 >= 2.75) || ($5 != "NA" && $6 <=0.05 && $7 >= 3.5)) print $0}' > deleterious_SIFT_kaki.tsv
awk '{print $1"\t"$2"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' kaki_whole-genome_SIFTannotations.txt | awk '{ if (($5 != "NA" && $6 > 0.05 && $7 >= 2.75) || ($5 !="NA" && $6 > 0.05 && $7 <= 3.5) print $0}' > tolerant_SIFT_kaki.tsv

awk '{print $1"\t"$2"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' KI_10x_whole-genome_polarized_filtered_SIFTannotations.txt | awk '{ if (($5 != "NA" && $6 <= 0.05 && $7 >= 2.75) || ($5 != "NA" && $6 <=0.05 && $7 >= 3.5)) print $0}' > deleterious_SIFT_KI_10x.tsv
awk '{print $1"\t"$2"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' KI_10x_whole-genome_polarized_filtered_SIFTannotations.txt | awk '{ if (($5 != "NA" && $6 > 0.05 && $7 >= 2.75) || ($5 !="NA" && $6 > 0.05 && $7 <= 3.5) print $0}' > tolerant_SIFT_KI_10x.tsv
```
SIFT sites were then extracted, and used to determine impact from VEP.
```

grep -v "#" GLOBAL_whole-genome.tsv | grep -v IMPACT=MODIFIER > VEP_impacts_fairy.tsv
grep -v "#" KI_whole-genome.tsv | grep -v IMPACT=MODIFIER > VEP_impacts_kaki.tsv
grep -v "#" KI_10x_whole-genome.tsv | grep -v IMPACT=MODIFIER > VEP_impacts_KI_10x.tsv

for POP in fairy kaki
    do
    awk '{print $1":"$2}' deleterious_SIFT_${POP}.tsv > SIFT_sites_${POP}.txt
    while read -r line
        do
        grep "$line" VEP_impacts_${POP}.tsv >> SIFT_VEP_intersect_${POP}.txt
    done < SIFT_sites_${POP}.txt
done
```
A brief check showed that all sites marked deleterious in SIFT had some impact in VEP.  
These intersecting sites were then extracted from the filtered SNPs, and allele frequency estimated for each fairy tern population.
```
while read -r line
    do
    CHROM=$(echo $line | awk '{print $1}')
    POSIT=$(echo $line | awk '{print $2}')
    AALLE=$(bcftools query -t ${CHROM}:${POSIT} -f '[%GT\n]' GLOBAL_whole-genome_polarised_filtered_AU.vcf | sed 's%0/0%0%g' | sed 's%0/1%1%g'| sed 's%1/1%2%g' | awk '{sum += $1}; END {print sum}')
    TALLE$(bcftools query -t ${CHROM}:${POSIT} -f '[%GT\n]' GLOBAL_whole-genome_polarised_filtered_TI.vcf | sed 's%0/0%0%g' | sed 's%0/1%1%g'| sed 's%1/1%2%g' | awk '{sum += $1}; END {print sum}')
    printf "$line\t$AALLE\tAU\n" >> pop_harmful_allele_frequency.tsv
    printf "$line\t$TALLE\tTI\n" >> pop_harmful_allele_frequency.tsv
done < deleterious_SIFT_fairy.tsv
```
Then we counted the number of alleles per site, per individual for all those nonsynonymous mutations classed as either intolerant, tolerant.  
```
while read -r line
    do
    CHROM=$(echo $line | awk '{print $1}')
    POSIT=$(echo $line | awk '{print $2}')
    bcftools query -t ${CHROM}:${POSIT} -f '[%CHROM\t%POS\t$SAMPLE\t%GT\tAU\tTolerant]' GLOBAL_whole-genome_filtered_AU.vcf >> indiv_harmful_allele_frequency.tsv
    bcftools query -t ${CHROM}:${POSIT} -f '[%CHROM\t%POS\t$SAMPLE\t%GT\tTI\tTolerant]' GLOBAL_whole-genome_filtered_TI.vcf >> indiv_harmful_allele_frequency.tsv
done < tolerant_SIFT_fairy.tsv

while read -r line
    do
    CHROM=$(echo $line | awk '{print $1}')
    POSIT=$(echo $line | awk '{print $2}')
    bcftools query -t ${CHROM}:${POSIT} -f '[%CHROM\t%POS\t$SAMPLE\t%GT\tAU\tIntolerant]' GLOBAL_whole-genome_filtered_AU.vcf >> indiv_harmful_allele_frequency.tsv
    bcftools query -t ${CHROM}:${POSIT} -f '[%CHROM\t%POS\t$SAMPLE\t%GT\tTI\tIntolerant]' GLOBAL_whole-genome_filtered_TI.vcf >> indiv_harmful_allele_frequency.tsv
done < deleterious_SIFT_fairy.tsv
```