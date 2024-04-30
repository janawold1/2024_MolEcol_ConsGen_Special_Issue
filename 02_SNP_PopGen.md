# Population genomics with ANGSD
[ANGSDv0.935](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) and ngsTools were used to estimate summary statistics for Australian fairy tern (*Sternual nereis nereis*) sampled from Western Australia and tara iti (*S. nereis davisae*) from Northland, NZ. The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. PCR duplicates were marked and sex chromosomes were excluded for SNP-based population analyses.  
```
ls ${DIR}{AU,SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}GLOBAL.list
ls ${DIR}{SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}TI.list
ls ${DIR}AU*_markdup_autosomes.bam > ${ANGSD}AU.list
ls ${DIR}H0*_markdup_autosomes.bam > ${ANGSD}KI.list

TREF=${DIR}Katies_genome/Katie_5kb_ragtag.fa
KREF=${DIR}kaki_genome/himNova-hic-scaff.fa

for POP in GLOBAL AU TI
    do
    angsd -P 16 -b ${ANGSD}${POP}.list -ref $REF -out ${ANGSD}qc/${POP}.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 800
done

angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -out ${ANGSD}qc/KI.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 2000
```
Courtesy of scripts provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.

|        Population        |minimum mapQ|minimum Q|minimum depth|maximum depth|Number of individuals|
| ------------------------ | ---------- | ------- | ----------- | ----------- | ------------------- |
|Australian Fairy Tern (WA)|     20     |    20   |     200     |     420     |         19          |
|         Tara iti         |     20     |    20   |     120     |     280     |         15          |
|     Global fairy tern    |     20     |    20   |     272     |     630     |         34          |
|           Kakī           |     20     |    20   |     700     |    1,200    |         26          |

 Below is code run for analyses performed on alignments of fairy terns to the tara iti reference scaffolded using the common tern assembly (see [00_genome_assembly.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/00_genome_assembly.md)) and alignments of kakī to the high-quality reference assembly for this species (see [here](https://www.genomics-aotearoa.org.nz/our-work/completed-projects/high-quality-genomes) for details). All analyses for each group were performed in a similar manner for downstream comparisons. Quality thresholds were adjusted as per the table above.  

The below table outlines the number of sites retained with the thresholds above for the WA population of Australian fairy tern (`AU`) dataset, the tara iti dataset (`TI`), the entire the fairy tern dataset (`GLOBAL`), and the kakī dataset (`KI`) for analyses.  
|   Data Set    |     AU    |    TI   |  GLOBAL |   Kakī  |
|:-------------:|:---------:|:-------:|:-------:|:-------:|
| Neutral Sites |  842,606  | 167,370 | XXX,XXX | XXX,XXX |
| Whole Genome  | 1,625,899 | 335,979 | XXX,XXX | XXX,XXX |

## Excluding putative coding Regions
Neither the tara iti or kakī reference assemblies have high-quality transcriptomes. This is significant given the potential of selection to influence some of the analyses performed here (e.g., estimates of historical N<sub>e</sub>). To assess the potential consequences of the inclusion of these regions, and to have a 'first look' at the genomic diversity around coding regions, we performed *ab initio* gene prediction using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus/tree/master) v3.5.0.  
```
augustus --sample=100 --alternatives=false --temperature=3 \
    --species=chicken reference/Katie_ragtag.fa > reference/Katie_AUGUSTUS.gff
augustus --sample=100 --alternatives=false --temperature=3 \
    --species=chicken reference/himNova-hic-scaff.fa > reference/kaki_AUGUSTUS.gff
```
Autosomal chromosomes were extracted from these predictions, and [BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.30.0 was then used to sort, and merge these putative gene regions. An additional 1kb of sequence on either side of annotations were included to reduce linkage.  

First, the annotations for sex chromosomes, unplaced scaffolds, or those that consistently did not adhere to expected sequencing depths were excluded as they likely represented misassemblies. Duplicate regions were then merged with `bedtools merge`.  
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
For the fairy tern, this left 27,337 regions covering 562,663,009 bp for analyses. This corresponds to roughly 47% of the genome.  

This file was then indexed for ANGSD with `angsd sites`.  
```
angsd sites index TI_scaffolded_neutral_regions.bed
```
## Ancestral alleles
It is important to have the ancestral state for some of the analyses below. To this end, we used short reads from the VGP common tern genome assembly ([*Sterna hirundo*](https://www.genomeark.org/vgp-all/Sterna_hirundo.html)) and short reads from the killdeer ([*Charadrius vociferus*](https://www.ncbi.nlm.nih.gov/sra/SRX328486[accn]) accession SRR943994) to polarise the site frequency spectrum for fairy terns and kakī respectively. We opted this approach over using these genome assemblies directly as the chromosomes had different sizes and were not compatible with ANGSD methods ([brief discussion here](https://www.biostars.org/p/298013/)). The reads for both common tern and avocet were trimmed, aligned, and duplicates marked in the same manner as the fairy tern and kakī population short-reads.  
```
angsd -P 16 -doFasta 1 -doCounts 1 -out ${ANGSD}bSteHir1_ancestral -i ${DIR}bSteHir1_markdup_autosomes.bam
angsd -P 16 -doFasta 1 -doCounts 1 -out ${ANGSD}killdeer_ancestral -i ${DIR}avocet_markdup_autosomes.bam
```
## Relatedness Estimates
To estimate relatedness among individuals within each of the three populations assessed, we first generated population specific genotype likelihoods.  
```
for POP in AU TI KI
    do
    if [[ "${POP}" == "AU" ]]
        then
        angsd -P 26 -b AU.list -ref $TREF -anc $ANC -sites $SITES -out inbreeding/AU \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
    elif [[ "${POP}" == "TI" ]]
        then
        angsd -P 26 -b TI.list -ref $TREF -anc $ANC -sites $SITES -out inbreeding/TI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
    else
        angsd -P 26 KI.list -ref $KREF -anc $KANC -sites $KSITES -out inbreeding/KI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
    fi
done
```
Next, we estimated relatedness, specifically the R<sub>AB</sub> metric, among individuals with [ngsRelate](https://github.com/ANGSD/NgsRelate).  
```
for POP in AU TI
    do
    zcat inbreeding/${POP}.mafs.gz | cut -f 7 | sed 1d > inbreeding/${POP}.freqs
    if [[ "$POP" == "AU" ]]
        then
        ngsRelate -g inbreeding/${POP}.glf -n 19 -f inbreeding/${POP}.freqs -O inbreeding/${POP}_relatedness
        else
        ngsRelate -g inbreeding/${POP}.glf -n 15 -f inbreeding/${POP}.freqs -O inbreeding/${POP}_relatedness
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
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $TREF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 10 --ignore-RG 0
        vcftools --bcf ${ANGSD}samtools/genotypes/${POP}_genotypes.bcf \
            --TsTv-summary --out ${ANGSD}samtools/genotypes/${POP}_genotypes
        elif [[ "$POP" == "TI" ]]
        then
        angsd -P 24 -b ${ANGSD}${POP}.list -ref $TREF -out ${ANGSD}samtools/genotypes/${POP}_genotypes \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -doPost 1 -doBcf 1 -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 10 --ignore-RG 0
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
./faSomeRecords reference/Katie_5kb_ragtag.fa Katie_autosomes2.bed reference/Katie_ragtag_autosomes2.fa
./faSomeRecords reference/himNova-hic-scaff.fa himNova_autsosomes.bed reference/himNova_autosomes.fa

samtools faidx reference/Katie_ragtag_autosomes2.fa
samtools faidx reference/himNova_autosomes.fa

REF=reference/Katie_ragtag_autosomes2.fa
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

for FILE in rohmu_5e5/H0*_subset.mid.hmmrohl.gz
    do
    BASE=$(basename $FILE _subset.mid.hmmrohl.gz)
    zcat $FILE | tail -n+2 | cut -f 5 | awk -v SAMP="$BASE" '{ if ( $1 <= 200000 ) len="short_ROH";
        else if ( $1 > 200000 && $1 <= 700000 ) len ="medium_ROH";
        else len ="long_ROH";
        print SAMP"\t"$1"\t"len"\tKI"; }' >> ROHs.tsv
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
### Global Heterozygosity
Here, we implemented a global (genome-wide heterozygosity) method from ANGSD. Essentially, this estimate is a proportion of heterozygous genotypes / genome size (excluding regions of the genome with low confidence). Unlike other runs of ANGSD, individual BAMs are used to estimate hetereozygosity, which is simply second value in the SFS/AFS.  
```
for SAMP in ${BAM}*_markdup_autosomes.bam
    do
    BASE=$(basename $SAMP _markdup_autosomes.bam)
    angsd -i $SAMP -anc $ANC -ref $REF -out ${ANGSD}samtools/heterozygosity/${BASE} -dosaf 1 -GL 1 -doCounts 1
    angsd -i $SAMP -anc $ANC -ref $REF -sites $SITES -out ${ANGSD}neutral/heterozygosity/${BASE} -dosaf 1 -GL 1 -doCounts 1
    realSFS ${ANGSD}samtools/heterozygosity/${BASE}.saf.idx > ${ANGSD}samtools/heterozygosity/${BASE}_est.ml
    realSFS ${ANGSD}neutral/heterozygosity/${BASE}.saf.idx > ${ANGSD}neutral/heterozygosity/${BASE}_est.ml
done
```
Once the SFS was estimated for each individual, the number of sites was estimated from the sum of all scaffold sizes included in the bam file and output to a file.
```
printf "Sample\tHeterozygosity\tTool\tPopulation\n" > ${ANGSD}individual_het.tsv 

for SAMP in ${ANGSD}${TOOL}/heterozygosity/*_est.ml
    do
    BASE=$(basename $SAMP _est.ml)
    TOT=$(awk '{print $1 + $2 + $3}' $SAMP)
    HET=$(awk -v var=$TOT '{print $2/var}' $SAMP)
    if [[ "$BASE" == *"AU"* ]]
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

## Population Structure
### MDS
To construct a MDS for fairy tern populations, we first ran ANGSD as below. Notably, the only change is the `-doGeno 8` flag.  
```
angsd -P 26 -b GLOBAL.list -ref $REF -anc $ANC -out ${ANGSD}structure_MDS/GLOBAL \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 34 -setMinDepth 272 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1
```
This command not only generates the required `*.mafs.gz`, but also the `*.geno.gz`, which contains posterior probabilities of all possible genotypes required for estimating genetic distance with [ngsDistv1.0.10](github.com/mfumagalli/ngsTools). In addition, a `pops.label` file denoting the population of origin (one entry on a new line for each samples) is necessary for estimating genetic distance.  
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
And finally plot the MDS.  
```
Rscript plotMDS.R -i ${ANGSD}distance/GLOBAL.mds \
    -c 1-2 -a ${ANGSD}structure_MDS/GLOBAL_clst \
    -o ${ANGSD}structure_MDS/GLOBAL_mds.pdf
```

## Summary Statistics
### Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS. Here, we are using the common tern as the ancestral state.  
```
ANC=${ANGSD}bSteHir1_ancestral.fasta
region=TI_scaffolded_neutral_regions.bed

for TOOL in samtools neutral
    do
    printf "STARTED RUNNING ANGSD FOR $TOOL DATA SET AT "
    date
    if [[ "$TOOL" == "samtools" ]]
        then
        angsd -P 16 -b ${ANGSD}AU.list -ref $REF -anc $ANC -out ${ANGSD}samtools/sfs/AU \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}TI.list -ref $REF -anc $ANC -out ${ANGSD}samtools/sfs/TI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -anc $KANC -out ${ANGSD}samtools/sfs/KI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doSaf 1
    else
        angsd -P 16 -b ${ANGSD}AU.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}neutral/sfs/AU \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 200 -setMaxDepth 420 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}TI.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}neutral/sfs/TI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doSaf 1
        angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -anc $KANC -sites $SITES -out ${ANGSD}neutral/sfs/KI \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ20 -minInd 15 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doSaf 1
    fi
done
```
Then `realSFS` was run to estimate the site frequency spectrum and theta-based statistics.  
```
realSFS -P 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}/sfs/AU.saf.idx > ${ANGSD}${TOOL}/sfs/AU.sfs
realSFS -P 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/TI.sfs
realSFS -P 40 -anc ${KANC} -ref ${KREF} ${ANGSD}${TOOL}sfs/KI.saf.idx > ${ANGSD}${TOOL}sfs/KI.sfs
realSFS -P 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}/sfs/AU.saf.idx ${ANGSD}${TOOL}/sfs/TI.saf.idx > ${ANGSD}${TOOL}/sfs/GLOBAL.sfs
```
Finally, for estimates of contemporary N<sub>e</sub>, we generated a SFS in *dadi* format. 
```
realSFS dadi -P 40 ${ANGSD}${TOOL}/sfs/AU.saf.idx ${ANGSD}${TOOL}/sfs/TI.saf.idx \
    -sfs ${ANGSD}${TOOL}/sfs/AU.sfs -sfs ${ANGSD}${TOOL}/sfs/TI.sfs -anc $ANC \
    -ref $REF > ${ANGSD}${TOOL}/sfs/GLOBAL.dadi
```

### Estimating F<sub>ST</sub>
Here we estimate F<sub>ST</sub> between Australian fairy tern populations in Western Australia and tara iti. To do so, we leverage the site frequency spectrum (SFS) generated above. We did this for both the conventional method for estimating the SFS and the bootstrapped and resampled method. Below, we outline how we prepared the bootstrapped method to generate F<sub>ST</sub> statistics for each replicate.  
```
realSFS fst index ${DIR}${TOOL}/sfs/AU.saf.idx ${DIR}${TOOL}/sfs/TI.saf.idx -sfs ${SFS} -fstout ${DIR}${TOOL}/distance/GLOBAL -whichFst 1
realSFS fst stats2 ${DIR}${TOOL}/distance/GLOBAL.fst.idx -win 10000 -step 1000 > ${DIR}${TOOL}/distance/GLOBAL_fst2_10KBwindow_1KBstep.tsv
realSFS fst stats2 ${DIR}${TOOL}/distance/GLOBAL.fst.idx -win 50000 -step 10000 > ${DIR}${TOOL}/distance/GLOBAL_fst2_50KBwindow_10KBstep.tsv
realSFS fst stats distance/GLOBAL_fst.idx -tole 1e-6 -anc $ANC -win 50000 -step 10000 -whichFst 1 > ${ANGSD}${TOOL}distance/GLOBAL_fst1_50KBwindow_10KBstep.tsv
```
The final command estimated a global weighted and unweighted F<sub>ST</sub> of `XXX` and `XXX` respectively for putatively neutral sites and an unweighted F<sub>ST</sub> of `0.002177` and weighted F<sub>ST</sub> of `0.843221` for all sites.  

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
            -GL 1 -doSaf 1 -pest ${ANGSD}${TOOL}sfs/${POP}.sfs
        realSFS -p 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx > ${ANGSD}${TOOL}diversity/${POP}_pest.sfs
        realSFS saf2theta ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx -sfs ${ANGSD}${TOOL}diversity/${POP}_pest.sfs -outname ${ANGSD}${TOOL}/diversity/${POP}_pest
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_10KBwindows_1KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_50KBwindows_10KBstep
    elif [[ "$POP" == "TI" ]]
        then
        angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}${TOOL}diversity/${POP}_pest \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 120 -setMaxDepth 280 -doCounts 1 \
            -GL 1 -doSaf 1 -pest ${ANGSD}${TOOL}sfs/${POP}.sfs
        realSFS -p 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx > ${ANGSD}${TOOL}diversity/${POP}_pest.sfs
        realSFS saf2theta ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx -sfs ${ANGSD}${TOOL}diversity/${POP}_pest.sfs -outname ${ANGSD}${TOOL}/diversity/${POP}_pest
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_10KBwindows_1KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_50KBwindows_10KBstep
    else
        angsd -P 26 -b ${ANGSD}${POP}.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}${TOOL}diversity/${POP}_pest \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 24 -setMinDepth 700 -setMaxDepth 1200 -doCounts 1 \
            -GL 1 -doSaf 1 -pest ${ANGSD}${TOOL}sfs/${POP}.sfs
        realSFS -p 40 -anc ${ANC} -ref ${REF} ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx > ${ANGSD}${TOOL}diversity/${POP}_pest.sfs
        realSFS saf2theta ${ANGSD}${TOOL}diversity/${POP}_pest.saf.idx -sfs ${ANGSD}${TOOL}diversity/${POP}_pest.sfs -outname ${ANGSD}${TOOL}/diversity/${POP}_pest
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 10000 -step 1000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_10KBwindows_1KBstep
        thetaStat do_stat ${ANGSD}${TOOL}/diversity/${POP}_pest.thetas.idx -win 50000 -step 10000 -outnames ${ANGSD}${TOOL}diversity/${POP}_pest_thetas_50KBwindows_10KBstep
done
```

### Estimating D<sub>XY</sub>
When aiming to estimate D<sub>XY</sub>, we need to first get a list of sites common to all populations. This is so sites that may be fixed in one population are included in our estimates. To generate this list of sites we ensured that the outputs were polarised with common tern as below.  
```
angsd -P -b ${ANGSD}GLOBAL.list -ref $REF -anc $ANC -sites $SITES -out ${ANGSD}neutral/diversity/GLOBAL_initial_Dxy \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMap 20 -minQ 20 -minInd 34 -setMinDepth 272 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3

angsd -P 16 -b ${ANGSD}GLOBAL.list -ref $REF -anc $ANC -out ${ANGSD}samtools/diversity/GLOBAL_initial_Dxy \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 34 -setMinDepth 272 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 
```
The sites for use in D<sub>XY</sub> calculations were then extracted and indexed with `ANGSD`.
```
zcat ${ANGSD}neutral/diversity/GLOBAL_initial_Dxy.mafs.gz | tail -n+2 | awk '{print $1"\t"$2}' > ${ANGSD}neutral_global_DXYsites.bed
zcat ${ANGSD}samtools/diversity/GLOBAL_initial_Dxy.mafs.gz | tail -n+2 | awk '{print $1"\t"$2}' > ${ANGSD}samtools_global_DXYsites.bed

angsd sites index ${ANGSD}neutral_global_DXYsites.bed
angsd sites index ${ANGSD}samtools_global_DXYsites.bed
```
And finally, population D<sub>XY</sub> estimates were generated.  
```
for POP in AU TI
    do
    angsd -P 26 -b ${DIR}${POP}.list -ref $REF -anc $ANC -sites ${DIR}neutral_global_DXYsites.bed \
        -out ${DIR}neutral/diversity/${POP}_DXY -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -doCounts 1 -GL 1 -doMajorMinor 1 -doMaf 1
    angsd -P 26 -b ${DIR}${POP}.list -ref $REF -anc $ANC -sites ${DIR}samtools_global_DXYsites.bed \
        -out ${DIR}samtools/diversity/${POP}_DXY -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -doCounts 1 -GL 1 -doMajorMinor 1 -doMaf 1
done
```

## Historical N<sub>e</sub> Inference with PSMC
We created a consensus diploid sequence for three tara iti, Australian fairy tern and kakī samples. The conversion from FASTQ format to PSMCfa must use a `gzipped` file as input.  
```
for BAM in variants/fairy_bam/nodup/{AU01,AU27,AU28,SND06,SP01,TI37}_nodup_autosomes.bam
    do
    BASE=$(basename $BAM _nodup_autosomes.bam)
    printf "BEGAN GENERATING CONSENSUS SEQUENCE FOR $BASE AT "
    date
    bcftools mpileup --threads 26 -Ou -Q 30 -q 20 -f $REF $BAM | \
        bcftools call --threads 26 -c | \
        vcfutils.pl vcf2fq -d 8 -D 125 -Q 30 > ${PSMC}fastq/${BASE}.fq
    gzip ${PSMC}fastq/${BASE}_psmc.fq
    fq2psmcfa -q20 ${PSMC}fastq/${BASE}_psmc.fq.gz > ${PSMC}psmcfa/${BASE}.psmcfa
    printf "FINISHED GENERATING CONSENSUS SEQUENCE FOR $BASE AT "
    date
done
```
Once the file was prepared, we ran PSMC under the the following conditions... We assumed a generation time of 3 years for tara iti as the age of first breeding ranges between 2-4 years of age. For kakī, this was increased to 6 years. The avian mutation rate is estimated to range from 1.23 x 10<sup>-9</sup> - 2.21 x 10<sup>-9</sup>. Assuming the lower range of this estimate, a rate of 3.69 x 10<sup>-9</sup> and 7.38 x 10<sup>-9</sup> were used for fairy terns and kakī respectively.  
```
for FA in psmc/psmcfa/*.psmcfa
    do
    SAMP=$(basename $FA .psmcfa)
    splitfa $FA > psmc/psmcfa/${SAMP}_split.psmcfa
    printf "STARTED RUNNING PSMC FOR $SAMP AT "
    date
    psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o psmc/out/${SAMP}_diploid.psmc $FA
    seq 1 100 | parallel -j 10 psmc -N30 -t5 -r5 -b -p "4+30*2+4+6+10" -o psmc/out/${SAMP}_round-{}.psmc psmc/psmcfa/${SAMP}_split.psmcfa
    cat psmc/out/${SAMP}_*.psmc > psmc/out/${SAMP}_combined.psmc
    printf "FINISHED RUNNING PSMC FOR $SAMP AT "
    date
done

psmc_plot.pl -u 3.69e-09 -g 3 psmc/AU psmc/out/AU*_combined.psmc
psmc_plot.pl -u 3.69e-09 -g 3 psmc/TI psmc/out/TI*_combined.psmc
psmc_plot.pl -u 7.38e-09 -g 6 psmc/KK psmc/out/KK*_combined.psmc
```
## Contemporary N<sub>e</sub> Inference with StairwayPlot2


## Population Demography and Connectivity Inference with GADMA
[GADMA](https://gadma.readthedocs.io/en/latest/user_manual/input_data/snp_data_format.html) leverages the joint SFS to infer the demographic history of multiple populations. It can implement [dadi](https://bitbucket.org/gutenkunstlab/dadi/), [moments](https://github.com/MomentsLD/moments), [momi2](https://github.com/popgenmethods/momi2/), and [momentsLD](https://github.com/MomentsLD/moments).  

GADMA can take multiple input formats. Here we estimated the joint SFS in ANGSD and realSFS as above, with the exception that the `dadi` flag was on when running realSFS. This output was then converted with a perl script [realSFS2dadi.pl](https://github.com/z0on/2bRAD_denovo/blob/master/realsfs2dadi.pl).
```
realSFS dadi -ref $REF -anc $ANC AU.saf.idx TI.saf.ids -sfs AU.sfs -sfs TI.sfs > GLOBAL.dadi
realSFS2dadi.pl GLOBAL.dadi 19 15 > GLOBAL_GADMA_SNPformat.txt
```
GADMA is relatively straightforward and easy to run. A parameter file defining specific settings can be used as input. One important parameter is `sequence length`, which denotes the number of sites used to build the data (SFS in our case). Fortunately, the realSFS programme includes this (`nSites`) as part of its progress output. For the neutral dataset this was 512,061,918 sites while it was 976,823,680 sites for the whole-genome data set.  

Because the SNPs in either SFS were not filtered for linkage, we updated the paramfile to reflect this and provided a directory for bootstrapping. This has important implications for model selection methods, see [here](https://gadma.readthedocs.io/en/latest/user_manual/input_data/input_data.html#extra-information-about-data) for more information.  

As with PSMC and StairwayPlot2 above, we ran GADMA with a conservative mutation rate of 1.23e-9 and a generation time of 3 years. Finally, the Selection and Dominance options were set to true prior to running with `gadma -p GADMA_neutral.params -o GADMA_neutral_moments/`.  

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
    bcftools mpileup --threads 26 -Ou -Q30 -q20 -f $REF -R ${BUSCO} \
        -b ${ANGSD}${POP}.list -a DP,AD,ADF,ADR,SP | \
        bcftools call --threads 26 -mv -O b -o ${LOAD}TI_calls.bcf
done
```