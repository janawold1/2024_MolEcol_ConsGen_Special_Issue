# Structural Variant Discovery and Analysis
Structural variant (SV) analyses followed a similar pipeline where variant discovery, variant quality filtering, and initial genotyping was performed under recommended pipelines. High quality genotypes representing homozygous reference, heterozygous and homozygous alternate calls were then plotted with SAMplot and evaluated using PlotCritic for call refinement. Calls passing evaluation with PlotCritic were then merged using SURVIVOR and used to construct a genome graph with the VG tools suite. Detailed methods for each of these steps using Delly, Manta and Smoove are outlined below.  

## Variant Discovery
[Manta](https://github.com/Illumina/manta) v1.6.0 was used to call SVs for Australian fairy tern and tara iti. Three samples had to be excluded for Manta to run, AU13, TI06, TI34 & TI35. The errors indicated issues with the proportion of reads and read depth statistics. Each chromosome was called independently with the `--callRegions` flag to save computational resource and increase efficiency. Five samples did not pass [Manta's read-pair orientation threshold](https://github.com/Illumina/manta/issues/168) of 90% and were excluded from SV discovery. This included two Australian samples and three tara iti samples (AU08, AU13, SND11, TI35).

Running Manta is relatively simple, with the initial configuration setup for each individual chromosome as per:
```
while read -r line
    do
    CHROM=$(echo $line | awk '{print $1}')
    REGION=$(echo $line | awk '{print $1":1-"$3}')
    mkdir manta/${CHROM}
    echo $REGION
    configManta.py --referenceFasta $REF --region $REGION --runDir manta/${CHROM}/ --bam AU01_nodup_autosomes.bam \
        --bam AU03_nodup_autosomes.bam --bam AU04_nodup_autosomes.bam --bam AU06_nodup_autosomes.bam --bam AU09_nodup_autosomes.bam \
        --bam AU14_nodup_autosomes.bam --bam AU17_nodup_autosomes.bam --bam AU20_nodup_autosomes.bam --bam AU21_nodup_autosomes.bam \
        --bam AU23_nodup_autosomes.bam --bam AU24_nodup_autosomes.bam --bam AU25_nodup_autosomes.bam --bam AU27_nodup_autosomes.bam \
        --bam AU28_nodup_autosomes.bam --bam AU29_nodup_autosomes.bam --bam AU30_nodup_autosomes.bam --bam AU33_nodup_autosomes.bam \
        --bam SND04_nodup_autosomes.bam --bam SND05_nodup_autosomes.bam --bam SND06_nodup_autosomes.bam --bam SND15_nodup_autosomes.bam \
        --bam SP02_nodup_autosomes.bam --bam SP03_nodup_autosomes.bam --bam SP07_nodup_autosomes.bam --bam TI21_nodup_autosomes.bam \
        --bam TI36_nodup_autosomes.bam --bam TI37_nodup_autosomes.bam --bam TI38_nodup_autosomes.bam --bam TI40_nodup_autosomes.bam \
        --bam TI41_nodup_autosomes.bam
done < reference/Katie_autosomes2.bed
``` 
And Manta executed on the resulting `runWorkflow.py` file in the designated output directories. Once individual chromosome runs were complete, VCFs were concatenated into a single file for SV quality filtering. 
```
echo CM0204{37..58}.1_RagTag | tr ' ' '\n' | while read chrom; do if [[ ! -d manta/${chrom} ]]; then continue; fi; find manta/${chrom}/results/variants/ -name diploidSV.vcf.gz | sort; done > manta/fairy_vcf_list

echo scaffold_{1..28} | tr ' ' '\n' | while read chrom; do if [[ ! -d manta/${chrom} ]]; then continue; fi; find manta/${chrom}/results/variants/ -name diploidSV.vcf.gz | sort; done > manta/kaki_vcf_list

bcftools concat --naive --file-list manta/fairy_vcf_list \
    -O z -o manta/01_fairy_SVcalls.vcf.gz

bcftools concat --naive --file-list manta/kaki_vcf_list \
    -O z -o manta/01_kaki_SVcalls.vcf.gz
```

## Filtering
Filtering of raw Manta calls was relatively simple. Because Manta outputs inversions as breakends ([see here](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md)), Inversion calls were converted from Breakends using the `convertInversions.py` script supplied by Manta.  
```
convertInversion.py ~/anaconda3/envs/samtools/bin/samtools $REF \
    manta/01_fairy_SVcalls.vcf.gz > manta/02_fairy_SVcalls_INVconverted.vcf

convertInversion.py ~/anaconda3/envs/samtools/bin/samtools $REF \
    manta/01_kaki_SVcalls.vcf.gz > manta/02_kaki_SVcalls_INVconverted.vcf
```
Then all reads had to pass all 'hard' filtering thresholds and have 'PRECISE' breakpoints.  
```
bcftools view -i 'FILTER=="PASS" & IMPRECISE==0' \
    -O z -o manta/03_fairy_filtered.vcf.gz manta/02_fairy_SVcalls_INVconverted.vcf
bcftools index manta/03_fairy_filtered.vcf.gz

bcftools view -i 'FILTER=="PASS" & IMPRECISE==0' \
    -O z -o manta/03_kaki_filtered.vcf.gz manta/02_kaki_SVcalls_INVconverted.vcf
bcftools index manta/03_kaki_filtered.vcf.gz
```
## Genotyping

Then performed genotyping as per:
```
while read -r line
    do
    printf "STARTED GENOTYPING $line AT "
    date
    graphtyper genotype_sv reference/Katie_5kb_ragtag.fa manta/03_fairy_filtered.vcf.gz --sams=fairy.list --threads=16 --region=$line
done < fairy_chroms.tsv
```
Genotypes for individual chromosomes were then merged and filtered to include SVs that passed in the `AGGREGATED` genotyping model for graphtyper and had a mean GQ>=25.
```
echo CM0204{37..58}.1_RagTag | tr ' ' '\n' | while read chrom; do if [[ ! -d graphtyper/sv_results/${chrom} ]]; then continue; fi; find graphtyper/sv_results/${chrom} -name "*.vcf.gz" | sort; done > graphtyper/fairy_vcf_list

bcftools concat --naive --file-list graphtyper/fairy_vcf_list \
    -O z -o graphtyper/fairy_manta_graphtyper.vcf.gz

bcftools view -i 'SVMODEL=="AGGREGATED" & mean(GQ)>=15 & N_PASS(GT=="mis")==0' \
    -O z -o graphtyper/fairy_manta_graphtyper_filtered.vcf.gz \
    graphtyper/fairy_manta_graphtyper.vcf.gz
```
This genotyping and filtering approach was also used for kakī.  

|   SV Type   | Raw FT Genotypes | Filtered FT Genotypes | KĪ Genotypes | KĪ GenoFiltered |
|:-----------:|:----------------:|:---------------------:|:------------:|:---------------:|
|  Deletion   |       5,948      |         4,245         |       XXX    |       XXX       |
| Duplication |        450       |           78          |       XXX    |       XXX       |
|  Insertion  |       3,659      |         2,653         |       XXX    |       XXX       |
|  Inversion  |        486       |           0           |       XXX    |        0        |


### Genotype Curation


## Diversity
### SV Type and Size distributions
XXX

### Individual Heterozygosity (H<sub>o</sub>)
Genotypes for each individual were estimated and then plotted in 06_visualisations.ipynb.  

```
while read -r line
    do
    echo $line
    TOT=$(bcftools query -s ${line} -i 'GT!="mis" & GQ>=15 & SVMODEL=="AGGREGATED"' -f '[%GT\n]' graphtyper/fairy_manta_graphtyper.vcf | wc -l)
    HET=$(bcftools query -s ${line} -i 'GT!="mis" & GQ>=15 & SVMODEL=="AGGREGATED"' -f '[%GT\n]' graphtyper/fairy_manta_graphtyper.vcf | sort | uniq -c | awk '{print $1}' | tr '\n' '\t' | awk -v var=$TOT '{print $2/var}')
    if [[ "$line" == "AU"* ]]
    then
        printf "$line\t$HET\tAU\n" >> graphtyper/individual_svHet.tsv
    elseif [[ "$line" == "H0"* ]]
        then
        printf "$line\t$HET\tKI\n" >> graphtyper/individual_svHet.tsv
    else
        printf "$line\t$HET\tNZ\n" >> graphtyper/individual_svHet.tsv
    fi
done < samples.txt
```

### Site Frequency Spectrum
To generate a site frequency spectrum for SVs, we generated population specific VCFs for 

## Population Structure and Differentiation
### MDS
We prepared the called SVs for ngsDist by using BCFtools to output the posterior probabilities (`FORMAT/PL`) for each individual's genotype.  
```
bcftools query -f '%CHROM\t%POS\t[%PL\t]\n' graphtyper/fairy_manta_graphtyper_filtered.vcf.gz | sed 's/,/\t/g' > graphtyper/fairy_mds/fairy_geno_svs.tsv
gzip graphtyper/fairy_mds/fairy_geno_svs.tsv

NSITES=$(zcat graphtyper/fairy_mds/fairy_geno_svs.tsv.gz | wc -l)
ngsDist --geno graphtyper/fairy_mds/fairy_geno_svs.tsv.gz --probs --log_scale \
    -n_ind 34 -n_sites $NSITES -labels graphtyper/fairy_mds/mds_list -o graphtyper/fairy_mds/GLOBAL_svDist
```

```
tail -n +3 graphtyper/fairy_mds/GLOBAL_svDist | Rscript --vanilla --slave getMDS.R \
    --no_header --data_symm -n 4 -m "mds" -o graphtyper/fairy_mds/GLOBAL_sv.mds
```
### F<sub>ST</sub>

### SVs and ROH Correlations - May be for another day
May need to consider multi-species comparisons.  

Same or different patterns observed in SNPs vs SVs? Genotyping rate will likely be low - highlight as an area of growth and how their inclusion in demographic stuff will be key.  

Similar patters off heterozygosity, Fst and structure?

To what extent do these SVs represent functional variation? This is the value add to the Sunnocks frame - We can't get at funcgtional variation so we look at all these indirect methods of inferring adaptive potential  

Knowing number of fixed differences.  