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
|    SV Type   | Number Called | Filtered Count | Curated SVs |
| ------------ | ------------- | -------------- | ----------- |
|  Breakends   |      748      |        0       |      0      |
|  Deletions   |     7,929     |      7,318     |    6,033    |
| Duplications |      599      |       190      |      0      |
|  Insertions  |      992      |       992      |      0      |
|  Inversions  |    15,640     |      6,103     |      0      |
|  **Total**   |  **25,908**   |   **14,603**   |  **6,033**  |

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
    -O v -o manta/results/variants/diploid_INV_convert_filtered.vcf \
    manta/results/variants/diploid_INV_convert.vcf
```

The raw calls initially comprised of:  
|    SV Type   | Number Called | Converting Inversions | Filtered Count | Curated SVs |
| ------------ | ------------- | --------------------- | -------------- | ----------- |
|  Breakends   |     4,858     |          728          |        0       |      0      |
|  Deletions   |     7,135     |         7,135         |      6,199     |    4,575    |
| Duplications |      666      |          666          |       425      |      0      |
|  Insertions  |     4,207     |         4,207         |      3,973     |      0      |
|  Inversions  |       0       |         2,065         |      1,508     |      0      |
|  **Total**   |  **16,866**   |      **14,801**       |   **12,105**   |  **4,575**  |

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
|    SV Type   | Number Called | Filtered Count | Curated SVs |
| ------------ | ------------- | -------------- | ----------- |
|  Breakends   |     5,908     |        0       |      0      |
|  Deletions   |     6,046     |       664      |     643     |
| Duplications |     1,267     |       58       |      0      |
|  Insertions  |       0       |        0       |      0      |
|  Inversions  |    11,812     |       345      |      0      |
|  **Total**   |  **25,033**   |    **1,067**   |   **643**   |

## Validating Filtered SV Calls
[SAMplot](https://github.com/ryanlayer/samplot) and [plotCritic](https://github.com/jbelyeu/PlotCritic) were used to evaluate SV calls from Delly, Smoove and Manta. However, SAMplot is only able to plot Deletions, Duplications and Inversions. In addition, except in a few instances, recovering the full sequences representing inversion and duplication haplotypes from short-read data alone is challenging. Only delection calls were curated for genome graph construction as they generally have clear support (read depth, and split reads) and obvious breakpoints. These are important aspects for accurate genotyping from graphs.  

Given the limitations of short-read data alone to resolve insertions, inversions, and duplications at the haplotype level, we limited our graphs to deletions. We required that deletion calls had to have exact breakpoints, and that they were supported by evidence from both split-read and read depth variation. There was no minimum or maximum size limitations.  

Deletions were curated for each of the three tools independently prior to merging. To generate input files representing heterozygous and homozygous alternate sites for these SVs, we focused on genotypes from each of the tools as they could provide clues as to which samples provided support for SV calls.  
```
for tool in delly manta smoove
    do
    echo "EXTRACTING SITES FROM $tool..."
    bcftools query -i 'SVTYPE!="INS" & GT="het"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | cut -f1-8 > ${tool}/samplot_het_sites.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="AA"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | cut -f1-8 > ${tool}/samplot_homozygousAlt_sites.tsv
done
```
Then used as input for plotting for 4 fairy tern samples.  
```
for tool in delly manta smoove
    do
    for geno in {homAlt,het}_n{1..4}
        do
            while read -r line
                do
                    printf "STARTED RUNNING SAMPLOT FOR $tool AT "
                    date
                    chrom=$(echo $line | awk '{print $1}')
                    start=$(echo $line | awk '{print $2}')
                    end=$(echo $line | awk '{print $3}')
                    type=$(echo $line | awk '{print $4}')
                    samp1=$(echo $line | awk '{print $5}')
                    samp2=$(echo $line | awk '{print $6}')
                    samp3=$(echo $line | awk '{print $7}')
                    samp4=$(echo $line | awk '{print $8}')
                    samplot plot -n $samp1 $samp2 $samp3 $samp4 \
                        -b align/nodup/${samp1}_nodup_autosomes.bam \
                        align/nodup/${samp2}_nodup_autosomes.bam \
                        align/nodup/${samp3}_nodup_autosomes.bam \
                        align/nodup/${samp4}_nodup_autosomes.bam \
                        -o ${tool}/samplot_outputs/${geno}/${chrom}_${start}_${end}_${type}.png \
                        -c ${chrom} -s ${start} -e ${end} -t ${type}
            done < ${tool}/samplot_${geno}_sites.tsv
    done
done
```
However, not all SV calls had support for a minimum of 4 individuals. These calls were filtered out (e.g., `awk { if (NR==5) {print $0} } > samplot_het_n1.tsv`) and plotted in a similar manner as above.

[PlotCritic](https://github.com/jbelyeu/PlotCritic) vX.X was used to evaluate whether calls had both read depth and split-read support prior to merging. But SMOOVE doesn't include the reference allele in the output VCF. To correct this, we normalised the SMOOVE output.  
```
bcftools norm --check-ref s --fasta-ref $REF -O z -o smoove/07_smoove_plotcritic_norm.vcf.gz smoove/06_smoove_plotcritic.vcf
bcftools sort -O z -o smoove/08_smoove_plotcritic_norm.sorted.vcf.gz smoove/07_smoove_plotcritic_norm.vcf.gz
tabix -p vcf smoove/08_smoove_plotcritic_norm.sorted.vcf.gz
```
Then all three call sets were merged.  
```
bcftools merge -m none -O z -o vg/DEL_calls.vcf.gz \
    delly/05_delly_plotcritic.vcf.gz \
    manta/03_manta_plotcritic.vcf.gz \
    smoove/08_smoove_plotcritic_norm.sorted.vcf.gz
tabix -p vcf DEL_calls.vcf.gz
```
This left 7,958 total deletions for genotyping.  

# Graph Construction and alignment
Graph construction, alignment and genotyping was performed using the VG toolsuite.  
```
vg autoindex \
	--workflow giraffe \
	-r ${dir}graphs/Jane_by_chromosome/Jane_chromosome_7.fa \
	-v ${dir}variants/SVs/sniffles/sniffles_chr7_PRECISE_minSUPP5.vcf.gz \
	-p ${dir}graphs/sniffles_total_giraffe -t 16
```
And finally, mapping can be done by passing the required indices.  
```
for fq2 in reads/*_R2.fq
    do
    indiv=$(basename $fq2 _R2.fq)
    printf "STARTED ALIGNING $indiv TO GRAPH AT "
    date
    vg giraffe \
	    -Z vg/fairy_tern.giraffe.gbz \
	    -m vg/fairy_tern.min \
	    -d vg/fairy_tern.dist \
	    -f reads/${indiv}_R1.fq \
	    -f reads/${indiv}_R2.fq > vg/gam/${indiv}_giraffe.gam
```
Then the output `.gbz` from `vg autoindex` was converted to `.xg` format for variant calls.
```
vg convert -x --drop-haplotypes \
	vg/fairy_tern.giraffe.gbz > vg/genotyping/fairy_tern.xg
```
And individual SV genotyping was performed.  
```
vg pack --threads 16 -x vg/fairy_tern.xg \
	-g vg/gam/${indiv}_giraffe.gam \
	-Q 5 -o vg/genotyping/${indiv}_giraffe.pack

vg call --genotype-snarls vg/genotyping/fairy_tern.xg \
	-k vg/genotyping/${indiv}_giraffe.pack > vg/genotyping/${indiv}_giraffe.vcf
```

# Population Analysis
 
### SVs and ROH Correlations - May be for another day
May need to consider multi-species comparisons.  

Same or different patterns observed in SNPs vs SVs? Genotyping rate will likely be low - highlight as an area of growth and how their inclusion in demographic stuff will be key.  

Similar patters off heterozygosity, Fst and structure?

To what extent do these SVs represent functional variation? This is the value add to the Sunnocks frame - We can't get at funcgtional variation so we look at all these indirect methods of inferring adaptive potential  

Knowing number of fixed differences.