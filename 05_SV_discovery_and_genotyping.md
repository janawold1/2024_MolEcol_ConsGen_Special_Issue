# Structural Variant Discovery and Analysis
## Delly Discovery
SV discovery with common tern and scaffolded tara iti assemblies
```
for cram in ${dir}*_nodup_autosomes.cram
    do
    base=$(basename ${cram} _nodup_autosomes.cram)
    printf "\nRUNNING DELLY ${base}\n"
    delly call -g ${ref} -o delly/raw_calls/${base}.bcf ${cram}
done
```
### Delly Filtering
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

These files were used to filter different SV types for quality. This is because Delly is best calling Inversions and Duplications greater than 300bp in length and INDELs greater than 50bp in length. SVs of different types were required to `PASS` all Delly filters and have `PRECISCE` breakpoints.  
```
bcftools view -i 'FILTER=="PASS" & INFO/PRECISE==1 & SVTYPE!="BND"' \
    -O b -o delly/02_SV_filtered.bcf delly/01_raw_merged_calls.bcf
```
These final SVs were used as input into the VG graph outlined below.  
## Manta Discovery
[Manta](https://github.com/Illumina/manta) v1.6.0 was used to call SVs for Australian fairy tern and tara iti. Three samples had to be excluded for Manta to run, AU13, TI06, TI34 & TI35. The errors indicated issues with the proportion of reads and read depth statistics. Each chromosome was called independently with the `--callRegions` flag to save computational resource and increase efficiency. Running Manta is relatively simple, with the initial configuration setup as per:
```
ref=/media/jana/BigData/tara_iti_publication/reference/TI_scaffolded_as_CT.fasta.gz
out=/media/jana/BigData/tara_iti_publication/manta/

configManta.py --referenceFasta ${ref} --runDir ${out} --callRegions ${chr} \
--bam AU01_nodup_autosomes.cram --bam AU03_nodup_autosomes.cram --bam AU04_nodup_autosomes.cram \
--bam AU06_nodup_autosomes.cram --bam AU08_nodup_autosomes.cram --bam AU09_nodup_autosomes.cram \
--bam AU14_nodup_autosomes.cram --bam AU17_nodup_autosomes.cram --bam AU20_nodup_autosomes.cram \
--bam AU21_nodup_autosomes.cram --bam AU23_nodup_autosomes.cram --bam AU24_nodup_autosomes.cram \
--bam AU25_nodup_autosomes.cram --bam AU27_nodup_autosomes.cram --bam AU28_nodup_autosomes.cram \
--bam AU29_nodup_autosomes.cram --bam AU30_nodup_autosomes.cram --bam AU33_nodup_autosomes.cram \
--bam SND04_nodup_autosomes.cram --bam SND05_nodup_autosomes.cram --bam SND06_nodup_autosomes.cram \
--bam SND11_nodup_autosomes.cram --bam SND15_nodup_autosomes.cram --bam SP01_nodup_autosomes.cram \
--bam SP02_nodup_autosomes.cram --bam SP03_nodup_autosomes.cram --bam SP07_nodup_autosomes.cram \
--bam TI21_nodup_autosomes.cram --bam TI22_nodup_autosomes.cram --bam TI36_nodup_autosomes.cram \
--bam TI37_nodup_autosomes.cram --bam TI38_nodup_autosomes.cram --bam TI39_nodup_autosomes.cram \
--bam TI40_nodup_autosomes.cram --bam TI41_nodup_autosomes.cram
``` 
And Manta executed on the resulting `runWorkflow.py` file in the designated output directories.  

### Manta Filtering
Filtering of raw Manta calls was relatively simple. First, Inversion calls were converted from Breakends using the `convertInversions.py` script supplied by Manta. Then all reads had to pass all 'hard' filtering thresholds and have 'PRECISE' breakpoints.  

Files for individual chromsomes were then concatenated into a single file with `bcftools`.  

## Smoove Discovery


### Smoove Filtering


## Merging Filtered SV Calls


# Graph Construction with VG


## Population Genotyping


### Genotype Filtering


### Population Analysis
