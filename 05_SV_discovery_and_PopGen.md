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

The raw calls initially comprised of:  
|    SV Type   | Total Number |
| ------------ | ------------ |
|  Breakends   |      875     |
|  Deletions   |     7,576    |
| Duplications |      611     |
|  Insertions  |      976     |
|  Inversions  |     1,558    |

### Delly Filtering
These files were used to filter different SV types for quality. This is because Delly is best calling Inversions and Duplications greater than 300bp in length and INDELs greater than 50bp in length. SVs of different types were required to `PASS` all Delly filters and have `PRECISCE` breakpoints.  
```
bcftools view -i 'FILTER=="PASS" & INFO/PRECISE==1 & SVTYPE!="BND"' \
    -O b -o delly/02_SV_filtered.bcf delly/01_raw_merged_calls.bcf
```
These final SVs were then merged with the other datasets and used as input into the VG graph as outlined below.  
## Manta Discovery
[Manta](https://github.com/Illumina/manta) v1.6.0 was used to call SVs for Australian fairy tern and tara iti. Three samples had to be excluded for Manta to run, AU13, TI06, TI34 & TI35. The errors indicated issues with the proportion of reads and read depth statistics. Each chromosome was called independently with the `--callRegions` flag to save computational resource and increase efficiency. Running Manta is relatively simple, with the initial configuration setup as per:
```
configManta.py --referenceFasta $REF --runDir $DIR \
        --bam AU01_markdup_autosomes.bam --bam AU03_markdup_autosomes.bam --bam AU04_markdup_autosomes.bam --bam AU06_markdup_autosomes.bam \
        --bam AU08_markdup_autosomes.bam --bam AU09_markdup_autosomes.bam --bam AU13_markdup_autosomes.bam --bam AU14_markdup_autosomes.bam \
        --bam AU17_markdup_autosomes.bam --bam AU20_markdup_autosomes.bam --bam AU21_markdup_autosomes.bam --bam AU23_markdup_autosomes.bam \
        --bam AU24_markdup_autosomes.bam --bam AU25_markdup_autosomes.bam --bam AU27_markdup_autosomes.bam --bam AU28_markdup_autosomes.bam \
        --bam AU29_markdup_autosomes.bam --bam AU30_markdup_autosomes.bam --bam AU33_markdup_autosomes.bam --bam SND04_markdup_autosomes.bam \
        --bam SND05_markdup_autosomes.bam --bam SND06_markdup_autosomes.bam --bam SND11_markdup_autosomes.bam --bam SND15_markdup_autosomes.bam \
        --bam SP01_markdup_autosomes.bam --bam SP02_markdup_autosomes.bam --bam SP03_markdup_autosomes.bam --bam SP07_markdup_autosomes.bam \
        --bam TI21_markdup_autosomes.bam --bam TI22_markdup_autosomes.bam --bam TI34_markdup_autosomes.bam --bam TI35_markdup_autosomes.bam \
        --bam TI36_markdup_autosomes.bam --bam TI37_markdup_autosomes.bam --bam TI38_markdup_autosomes.bam --bam TI40_markdup_autosomes.bam \
        --bam TI41_markdup_autosomes.bam
``` 
And Manta executed on the resulting `runWorkflow.py` file in the designated output directories.  

The raw calls initially comprised of:  
|    SV Type   | Total Number |
| ------------ | ------------ |
|  Breakends   |      XXX     |
|  Deletions   |     X,XXX    |
| Duplications |      XXX     |
|  Insertions  |      XXX     |
|  Inversions  |     X,XXX    |

### Manta Filtering
Filtering of raw Manta calls was relatively simple. First, Inversion calls were converted from Breakends using the `convertInversions.py` script supplied by Manta. Then all reads had to pass all 'hard' filtering thresholds and have 'PRECISE' breakpoints.  

Files for individual chromsomes were then concatenated into a single file with `bcftools`.  

## Smoove Discovery

```
for cram in ${dir}*_nodup_autosomes.cram
    do
    base=$(basename ${cram} _nodup_autosomes.cram)
    echo "RUNNING SMOOVE CALL FOR ${base}..."
    smoove call --name ${base} --fasta ${ref} --outdir smoove/raw_calls/ --genotype ${cram}
done
```
Then calls for all individuals were merged into a single file.  
```
echo "Merging all called variants..."
smoove merge --name bwa_smoove -f ${ref} --outdir ${out} ${out}SV_calls_male/*.genotyped.vcf.gz ${out}SV_calls_female/*.genotyped.vcf.gz
```
The raw calls initially comprised of:  
|    SV Type   | Total Number |
| ------------ | ------------ |
|  Breakends   |      XXX     |
|  Deletions   |     X,XXX    |
| Duplications |      XXX     |
|  Insertions  |      XXX     |
|  Inversions  |     X,XXX    |

### Smoove Filtering
```
Used bwa_smoove_annotated.vcf.gz for unfiltered SV summary stats. Created filtered file as per:
bcftools view -t ^NC_044302.2 -O v -o ${out}01_smoove_unfiltered.vcf ${out}bwa_smoove.annotated.vcf.gz
bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0-168] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0-168] > 1.3) | (SVTYPE = "INV")' \
    -O v -o ${out}02_smoove_SVfiltered.vcf ${out}01_smoove_unfiltered.vcf
bcftools view -i 'SVLEN < 50000 & SVLEN > -50000' -O v -o 03_smoove_length_filtered.vcf 02_smoove_SVfiltered.vcf
bcftools view -i '(MSHQ>=3)' -O v -o ${out}04_smoove_genofiltered.vcf ${out}03_smoove_length_filtered.vcf
```

## Merging Filtered SV Calls


## Graph Construction with VG


## Population-scale Genotyping


### Genotype Filtering


## Population Analysis
