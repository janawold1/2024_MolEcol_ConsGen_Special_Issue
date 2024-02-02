# Structural Variant Discovery and Analysis
## Delly Discovery
SV discovery with common tern and scaffolded tara iti assemblies
```
while read -r BAM
    do
    BASE=$(basename ${BAM} _markdup_autosomes.cram)
    printf "BEGAN RUNNING DELLY FOR ${BASE} AT "
    date
    delly call -g ${REF} -o delly/raw_calls/${BASE}.bcf ${BAM}
    printf "FINISHED AT "
    date
done < sample_list.txt
```
Once initial calls were made, the file was merged and filtered for two different minimum sizes.  
```
delly merge -o delly/01_raw_merged_calls.bcf delly/raw_calls/*.bcf
```
All breakend calls were excluded as they likely indicate unresolved complex variation. The remanining SVs were required to `PASS` all Delly filters and have `PRECISE` breakpoints.  
```
bcftools view -i 'FILTER=="PASS" & INFO/PRECISE==1 & SVTYPE!="BND"' \
    -O b -o delly/02_SV_filtered.bcf delly/01_raw_merged_calls.bcf
```
Number of SVs called and number passing filtering thresholds:  
|    SV Type   | Number Called | Filtered Count |
| ------------ | ------------- | -------------- |
|  Breakends   |      750      |        0       |
|  Deletions   |     7,525     |      7,313     |
| Duplications |      600      |       190      |
|  Insertions  |      991      |       991      |
|  Inversions  |    15,637     |      6,100     |
|  **Total**   |  **25,503**   |   **14,594**   |

These final SVs were then merged with the other datasets and used as input into the VG graph as outlined below.  
## Manta Discovery
[Manta](https://github.com/Illumina/manta) v1.6.0 was used to call SVs for Australian fairy tern and tara iti. Three samples had to be excluded for Manta to run, AU13, TI06, TI34 & TI35. The errors indicated issues with the proportion of reads and read depth statistics. Each chromosome was called independently with the `--callRegions` flag to save computational resource and increase efficiency. Five samples did not pass [Manta's read-pair orientation threshold](https://github.com/Illumina/manta/issues/168) of 90% and were excluded from SV discovery. This included two Australian samples and three tara iti samples (AU08, AU13, SND11, TI34, TI35).

Running Manta is relatively simple, with the initial configuration setup as per:
```
bgzip Katie_autosomes2.bed
tabix -s 1 -b 2 -e 3 Katie_autosomes2.bed.gz 

configManta.py --referenceFasta $REF --callRegions reference/Katie_autosomes2.bed.gz --runDir $DIR --bam AU01_nodup_autosomes.bam \
    --bam AU03_nodup_autosomes.bam --bam AU04_nodup_autosomes.bam --bam AU06_nodup_autosomes.bam --bam AU09_nodup_autosomes.bam \
    --bam AU14_nodup_autosomes.bam --bam AU17_nodup_autosomes.bam --bam AU20_nodup_autosomes.bam --bam AU21_nodup_autosomes.bam \
    --bam AU23_nodup_autosomes.bam --bam AU24_nodup_autosomes.bam --bam AU25_nodup_autosomes.bam --bam AU27_nodup_autosomes.bam \
    --bam AU28_nodup_autosomes.bam --bam AU29_nodup_autosomes.bam --bam AU30_nodup_autosomes.bam --bam AU33_nodup_autosomes.bam \
    --bam SND05_nodup_autosomes.bam --bam SND06_nodup_autosomes.bam --bam SND15_nodup_autosomes.bam --bam SP01_nodup_autosomes.bam \
    --bam SP02_nodup_autosomes.bam --bam SP03_nodup_autosomes.bam --bam SP07_nodup_autosomes.bam --bam TI21_nodup_autosomes.bam \
    --bam TI22_nodup_autosomes.bam --bam TI36_nodup_autosomes.bam --bam TI37_nodup_autosomes.bam --bam TI38_nodup_autosomes.bam \
    --bam TI40_nodup_autosomes.bam --bam TI41_nodup_autosomes.bam
``` 
And Manta executed on the resulting `runWorkflow.py` file in the designated output directories.  

The raw calls initially comprised of:  
|    SV Type   | Number Called | Filtered Count |
| ------------ | ------------- | -------------- |
|  Breakends   |      XXX      |        0       |
|  Deletions   |     X,XXX     |      X,XXX     |
| Duplications |      XXX      |       XXX      |
|  Insertions  |      XXX      |       XXX      |
|  Inversions  |     X,XXX     |      X,XXX     |
|  **Total**   |  **XX,XXX**   |   **XX,XXX**   |

Filtering of raw Manta calls was relatively simple. First, Inversion calls were converted from Breakends using the `convertInversions.py` script supplied by Manta. Then all reads had to pass all 'hard' filtering thresholds and have 'PRECISE' breakpoints.  

Files for individual chromsomes were then concatenated into a single file with `bcftools`.  

## Smoove Discovery

```
for BAM in ${dir}*_markdup_autosomes.cram
    do
    BASE=$(basename ${BAM} _markdup_autosomes.bam)
    echo "RUNNING SMOOVE CALL FOR ${BASE}..."
    smoove call --name ${BASE} --fasta ${REF} --outdir smoove/raw_calls/ --genotype ${BAM}
done
```
Then calls for all individuals were merged into a single file.  
```
smoove merge --name 01_raw_calls \
    --fasta ${ref} --outdir smoove/ \
    smoove/raw_calls/*.genotyped.vcf.gz
```
After merging, `INFO/IMPRECISE` calls were excluded.
```
bcftools view -i 'INFO/IMPRECISE==0 & INFO/SVTYPE!="BND"' \
    -O v -o smoove/02_smoove_precise.vcf \
    smoove/01_raw_merged.sites.vcf
```

The raw calls initially comprised of:  
|    SV Type   | Number Called | Filtered Count |
| ------------ | ------------- | -------------- |
|  Breakends   |     5,904     |        0       |
|  Deletions   |     6,023     |       666      |
| Duplications |     1,268     |       59       |
|  Insertions  |       0       |        0       |
|  Inversions  |    11,813     |       345      |
|  **Total**   |  **25,008**   |    **1,070**   |


## Merging Filtered SV Calls


## Graph Construction with VG


## Population-scale Genotyping


### Genotype Filtering


## Population Analysis
 
 ### SVs and ROH Correlations - May be for another day
 May need to consider multi-species comparisons. 
 Same or different patterns observed in SNPs vs SVs? Genotyping rate will likely be low - highlight as an area of growth and how their inclusion in demographic stuff will be key. 

 Similar patters off heterozygosity, Fst and structure?

 To what extent do these SVs represent functional variation? This is the value add to the Sunnocks frame - We can't get at funcgtional variation so we look at all these indirect methods of inferring adaptive potential 

 Knowing number of fixed differences.  