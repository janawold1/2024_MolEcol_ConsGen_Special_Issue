# Structural Variant Discovery and Analysis
Structural variant (SV) analyses followed a similar pipeline where variant discovery, variant quality filtering, and initial genotyping was performed under recommended pipelines. High quality genotypes representing homozygous reference, heterozygous and homozygous alternate calls were then plotted with SAMplot and evaluated using PlotCritic for call refinement. Calls passing evaluation with PlotCritic were then merged using SURVIVOR and used to construct a genome graph with the VG tools suite. Detailed methods for each of these steps using Delly, Manta and Smoove are outlined below.  
## Delly Discovery
SV discovery with common tern and scaffolded tara iti assemblies
```
for BAM in nodup/*_nodup_autosomes.bam
    do
    BASE=$(basename ${BAM} _nodup_autosomes.bam)
    printf "BEGAN RUNNING DELLY FOR ${BASE} AT "
    date
    delly call -g ${REF} -o delly/raw_calls/${BASE}_nodup.bcf ${BAM}
    printf "FINISHED AT "
    date
done
```
Once initial calls were made, the file was merged and filtered for two different minimum sizes.  
```
delly merge -o delly/01_raw_merged_calls.bcf delly/raw_calls/*_nodup.bcf
```
All breakend calls were excluded as they likely indicate unresolved complex variation. The remanining SVs were required to `PASS` all Delly filters and have `PRECISE` breakpoints.  
```
bcftools view -i 'FILTER=="PASS" & INFO/PRECISE==1 & SVTYPE!="BND"' \
    -O b -o delly/02_SV_filtered.bcf delly/01_raw_merged_calls.bcf
```
Number of SVs called and number passing filtering thresholds:  
|    SV Type   | Number Called | Filtered Count |
| ------------ | ------------- | -------------- |
|  Breakends   |      748      |        0       |
|  Deletions   |     7,929     |      7,318     |
| Duplications |      599      |       190      |
|  Insertions  |      992      |       992      |
|  Inversions  |    15,640     |      6,103     |
|  **Total**   |  **25,908**   |   **14,603**   |

These final SVs were then merged with the other datasets and used as input into the VG graph as outlined below.  
## Manta Discovery
[Manta](https://github.com/Illumina/manta) v1.6.0 was used to call SVs for Australian fairy tern and tara iti. Three samples had to be excluded for Manta to run, AU13, TI06, TI34 & TI35. The errors indicated issues with the proportion of reads and read depth statistics. Each chromosome was called independently with the `--callRegions` flag to save computational resource and increase efficiency. Five samples did not pass [Manta's read-pair orientation threshold](https://github.com/Illumina/manta/issues/168) of 90% and were excluded from SV discovery. This included two Australian samples and three tara iti samples (AU08, AU13, SND11, TI35).

Running Manta is relatively simple, with the initial configuration setup as per:
```
bgzip Katie_autosomes2.bed
tabix -s 1 -b 2 -e 3 Katie_autosomes2.bed.gz 

configManta.py --referenceFasta $REF --callRegions reference/Katie_autosomes2.bed.gz --runDir $DIR --bam AU01_nodup_autosomes.bam \
    --bam AU03_nodup_autosomes.bam --bam AU04_nodup_autosomes.bam --bam AU06_nodup_autosomes.bam --bam AU09_nodup_autosomes.bam \
    --bam AU14_nodup_autosomes.bam --bam AU17_nodup_autosomes.bam --bam AU20_nodup_autosomes.bam --bam AU21_nodup_autosomes.bam \
    --bam AU23_nodup_autosomes.bam --bam AU24_nodup_autosomes.bam --bam AU25_nodup_autosomes.bam --bam AU27_nodup_autosomes.bam \
    --bam AU28_nodup_autosomes.bam --bam AU29_nodup_autosomes.bam --bam AU30_nodup_autosomes.bam --bam AU33_nodup_autosomes.bam \
    --bam SND04_nodup_autosomes.bam --bam SND05_nodup_autosomes.bam --bam SND06_nodup_autosomes.bam --bam SND15_nodup_autosomes.bam \
    --bam SP02_nodup_autosomes.bam --bam SP03_nodup_autosomes.bam --bam SP07_nodup_autosomes.bam --bam TI21_nodup_autosomes.bam \
    --bam TI36_nodup_autosomes.bam --bam TI37_nodup_autosomes.bam --bam TI38_nodup_autosomes.bam --bam TI40_nodup_autosomes.bam \
    --bam TI41_nodup_autosomes.bam
``` 
And Manta executed on the resulting `runWorkflow.py` file in the designated output directories.  

Filtering of raw Manta calls was relatively simple. Because Manta outputs inversions as breakends ([see here](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md)), Inversion calls were converted from Breakends using the `convertInversions.py` script supplied by Manta.  
```
convertInversion.py ~/anaconda3/envs/samtools/bin/samtools $REF \
    manta/results/variants/diploidSV.vcf.gz > manta/results/variants/diploid_INV_convert.vcf
```
Then all reads had to pass all 'hard' filtering thresholds and have 'PRECISE' breakpoints.  
```
bcftools view -i 'INFO/SVTYPE!="BND" & FILTER="PASS" & IMPRECISE=0' \
    -O z -o manta/results/variants/diploid_INV_filtered.vcf.gz \
    manta/results/variants/diploid_INV_convert.vcf
```

The raw calls initially comprised of:  
|    SV Type   | Number Called | Converting Inversions | Filtered Count |
| ------------ | ------------- | --------------------- | -------------- |
|  Breakends   |    20,658     |          742          |        0       |
|  Deletions   |     7,070     |         7,070         |      6,158     |
| Duplications |      663      |          663          |       420      |
|  Insertions  |     4,065     |         4,065         |      3,845     |
|  Inversions  |       0       |         9,958         |      8,208     |
|  **Total**   |  **32,456**   |      **22,498**       |   **18,631**   |  

## Smoove Discovery

```
for BAM in ${dir}*_nodup_autosomes.bam
    do
    BASE=$(basename ${BAM} _nodup_autosomes.bam)
    echo "RUNNING SMOOVE CALL FOR ${BASE}..."
    smoove call --name ${BASE} --fasta ${REF} --outdir smoove/raw_calls/ --genotype ${BAM}
done
```
Then calls for all individuals were merged into a single file.  
```
smoove merge --name 01_raw_merged --fasta ${REF} --outdir smoove/ smoove/raw_calls/*.genotyped.vcf.gz
```
After merging, `INFO/IMPRECISE` and all breakend (`SVTYPE=BND`) calls were excluded.
```
bcftools view -i 'INFO/IMPRECISE==0 & INFO/SVTYPE!="BND"' -O v -o smoove/02_smoove_precise.vcf smoove/01_raw_merged.sites.vcf
```

The raw calls initially comprised of:  
|    SV Type   | Number Called | Filtered Count |
| ------------ | ------------- | -------------- |
|  Breakends   |     5,908     |        0       |
|  Deletions   |     6,046     |       664      |
| Duplications |     1,267     |       58       |
|  Insertions  |       0       |        0       |
|  Inversions  |    11,812     |       345      |
|  **Total**   |  **25,033**   |    **1,067**   |

## Validating Filtered SV Calls
[SAMplot](https://github.com/ryanlayer/samplot) and [plotCritic](https://github.com/jbelyeu/PlotCritic) vX.X were used to evaluate SV calls from Delly, Smoove and Manta. However, SAMplot is only able to plot Deletions, Duplications and Inversions. For this step, sites for each of the tools were first extracted with BCFtools.  
```
for tool in delly manta smoove
    do
    echo "EXTRACTING SITES FROM $tool..."
    bcftools query -i 'SVTYPE!="INS"' -f '%CHROM\t%POS\t%END\t%SVTYPE' ${tool}/02_filteredSVs.vcf > ${tool}/samplot_sites.tsv
done
```
Then used as input for plotting for 6 Australian fairy terns (3 female, 3 male) and 6 tara iti (3 female, 3 male).  
```
for tool in delly manta smoove
    do
    while read -r line
        do
        printf "STARTED RUNNING SAMPLOT FOR $tool AT "
        date
```
## Merging Filtered SV Calls


# Graph Construction with VG


## Population-scale Genotyping


## Genotype Filtering


# Population Analysis
 
### SVs and ROH Correlations - May be for another day
 May need to consider multi-species comparisons. 
 Same or different patterns observed in SNPs vs SVs? Genotyping rate will likely be low - highlight as an area of growth and how their inclusion in demographic stuff will be key. 

 Similar patters off heterozygosity, Fst and structure?

 To what extent do these SVs represent functional variation? This is the value add to the Sunnocks frame - We can't get at funcgtional variation so we look at all these indirect methods of inferring adaptive potential 

 Knowing number of fixed differences.  