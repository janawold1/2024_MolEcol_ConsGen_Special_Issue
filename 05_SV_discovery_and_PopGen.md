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
    -O z -o graphtyper/01_fairy_genotypes.vcf.gz

bcftools view -i 'SVMODEL="AGGREGATED" & mean(GQ)>=25 & N_PASS(GT="mis")=0 & FILTER="PASS" & AC<68 & SVTYPE!="INS"' \
    -O z -o graphtyper/02_fairy_genotypes_filtered_noINS.vcf.gz \
    graphtyper/01_fairy_genotypes.vcf.gz
bcftools view -i 'SVMODEL="AGGREGATED" & mean(GQ)>=25 & N_PASS(GT="mis")=0 & FILTER="PASS" & AC<68 & SVTYPE="INS" & INFO/SEQ!="."' \
    -O z -o graphtyper/02_fairy_genotypes_filtered_INS.vcf.gz \
    graphtyper/01_fairy_genotypes.vcf.gz

bcftools index grapthyper/02_fairy_genotypes_INS.vcf.gz
bcftools index grapthyper/02_fairy_genotypes_noINS.vcf.gz

bcftools concat -a -O z -o graphtyper/02_fairy_genotypes_filtered.vcf.gz 02_fairy_genotypes_filtered_noINS.vcf.gz 02_fairy_genotypes_filtered_INS.vcf.gz

rm graphtyper/02_fairy_genotypes_*
```
This genotyping and filtering approach was also used for kakī.  

|   SV Type   | Raw FT Genotypes | Filtered FT Genotypes | KĪ Genotypes | KĪ GenoFiltered |
|:-----------:|:----------------:|:---------------------:|:------------:|:---------------:|
|  Deletion   |       5,948      |         1,087         |     6,912    |       533       |
| Duplication |        450       |           41          |       674    |        15       |
|  Insertion  |       3,659      |         1,544         |     5,304    |      1,192      |
|  Inversion  |        486       |           0           |       61     |        0        |
|    Total    |       X,XXX      |         2,672         |     X,XXX    |      1,740      |


## Diversity
### SV Type and Size distributions
XXX

### Individual Heterozygosity (H<sub>o</sub>)
Genotypes for each individual were estimated and then plotted in 06_visualisations.ipynb.  

```
while read -r line
    do
    TOT=$(bcftools query -s ${line} -i 'GT!="mis" & GQ>=25 & SVMODEL=="AGGREGATED" & FORMAT/FT=="PASS"' -f '[%GT\n]' graphtyper/02_fairy_genotypes_filtered.vcf.gz | wc -l)
    HET=$(bcftools query -s ${line} -i 'GT!="mis" & GQ>=25 & SVMODEL=="AGGREGATED" & FORMAT/FT=="PASS" & GT=="het"' -f '[%GT\n]' graphtyper/02_fairy_genotypes_filtered.vcf.gz | wc -l | awk -v var=$TOT '{print $1/var}')
    if [[ "$line" == "AU"* ]]
    then
        printf "$line\t$HET\t$TOT\tAU\n"
        printf "$line\t$HET\t$TOT\tAU\n" >> graphtyper/individual_svHet.tsv
    else
        printf "$line\t$HET\t$TOT\tTI\n"
        printf "$line\t$HET\t$TOT\tTI\n" >> graphtyper/individual_svHet.tsv
    fi
done < samples.txt
```

## Population Structure and Differentiation
### MDS
We prepared the called SVs for ngsDist by using BCFtools to output the posterior probabilities (`FORMAT/PL`) for each individual's genotype.  
```
bcftools query -f '%CHROM\t%POS\t[%PL\t]\n' graphtyper/02_fairy_genotypes_filtered.vcf.gz | sed 's/,/\t/g' | gzip -c > graphtyper/fairy_mds/fairy_SVgenos.tsv.gz

NSITES=$(zcat graphtyper/fairy_mds/fairy_SVgenos.tsv.gz | wc -l)
ngsDist --geno graphtyper/fairy_mds/fairy_SVgenos.tsv.gz --probs --log_scale \
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

## Intersections with BUSCO genes
After calling BUSCO sites for the tara iti and kakī genome assemblies, we extracted regions corresponding with complete, single-copy BUSCOs and those regions corresponding to our autosomal chromosomes.
```
cut -f2,3,4,5 run_aves_odb10/full_table.tsv | grep Complete | awk '{print $2"\t"$3"\t"$4"\t"$5}' > Katie_scBUSCOs.bed

while read -r line
    do
    CHROM=(echo $line | awk '{print $1}')
    printf "EXTRACTING CHROMOSOME $CHROM...\n"
    grep "$CHROM" Katie_scBUSCOs.bed >> Katie_autosome_scBUSCOs.bed
done < Katie_autosomes2.bed

sort -k1,1 -k2,2n Katie_autosome_scBUSCOs.bed > Katie_autosome_scBUSCOs.sorted.bed
```
After extracting BUSCOs corresponding to our autosomal scaffolds, we identified SVs overlapping these genes by extracting SV regions and finding intersecting sites with BCFtools and BEDtools.  

```
bcftools query -f '%CHROM\t%POS\t%SVLEN\n' 02_fairy_genotypes_filtered.vcf.gz | awk '{print $1"\t"$2"\t"$2+$3}' | sort -k1,1 -k2,2n > fairy_sv.bed

bedtools intersect -wa -a fairy_sv.bed  -b Katie_autsome_scBUSCOs.sorted.bed | awk '{print $1"\t"$2}' > fairy_BUSCO_svs.txt

bcftools view -T fairy_BUSCO_svs.txt -O z -o 03_fairy_BUSCO.vcf.gz 02_fairy_genotypes_filtered.vcf.gz
```
## Site Frequency Spectrum
To generate a site frequency spectrum for SVs, we generated population specific VCFs for each group and attempted to polarise the VCFs with genotypes from reads used to polarise the SNP-base SFS. However, SVs were not reliably genotyped in either tern outgroup or the pied avocet for kakī. This made comparisons of allele counts intersecting BUSCO genes with SNP-based estimates of load challenging. However, given that the BUSCO regions were identified as single-copy and complete in our reference individuals, we assumed any SV intersecting with these regions would be under strong purifying selection and therefore is the 'derived' allele. To ensure comparisons between fairy terns was comparable, we excluded any SVs private to one group as below.   

```
mkdir ./{polarise_TI,polarise_KI}

bcftools view -S TI_sampleIDs.txt 02_fairy_genotypes_filtered.vcf.gz | bcftools view -i 'GT=="alt"' -O z -o 03_TI_genotypes_filtered.vcf.gz
bcftools view -S AU_sampleIDs.txt 02_fairy_genotypes_filtered.vcf.gz | bcftools view -i 'GT=="alt"' -O z -o 03_AU_genotypes_filtered.vcf.gz

bcftools index 03_TI_genotypes_filtered.vcf.gz
bcftools index 03_AU_genotypes_filtered.vcf.gz

bcftools isec 03_TI_genotypes_filtered.vcf.gz 03_AU_gentoypes_filtered.vcf.gz -p polarise_TI/
```
Finally, we created our individuals and population frequency files.  
```
bcftools query -T fairy_BUSCOsvs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tTI\tIntersecting\n]' polarise_TI/0003.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^fairy_BUSCOsvs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tTI\tNonintersecting\n]' polarise_TI/0003.vcf >> indiv_harmful_SVallele_frequency.tsv

bcftools query -T fairy_BUSCOsvs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tAU\tIntersecting\n]' polarise_TI/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^fairy_BUSCOsvs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tAU\tNonintersecting\n]' polarise_TI/0002.vcf >> indiv_harmful_SVallele_frequency.tsv

bcftools query -T kaki_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tKI\tIntersecting\n]' polarise_KI/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^kaki_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tKI\tNonintersecting\n]' polarise_KI/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
```
