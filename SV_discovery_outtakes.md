# SNPs
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
# SVs

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

bcftools view -i 'SVTYPE=="DEL"' -O b -o delly/03_filtered_DEL.bcf delly/02_SV_filtered.bcf
```
Number of SVs called and number passing filtering thresholds:  
|    SV Type   | # FT SVs Called | # KI SVs Called | # FT Filtered | # KI Filtered | FT Curated SVs | KI Curated SVs |
| ------------ | --------------- | --------------- | ------------- | ------------- | -------------- | -------------- |
|  Breakends   |       748       |       728       |       0       |       0       |        0       |       0        |
|  Deletions   |      7,929      |      5,009      |     7,318     |     4,497     |      6,033     |     X,XXX      |
| Duplications |       599       |       371       |      190      |       83      |        0       |       0        |
|  Insertions  |       992       |       971       |      992      |      971      |        0       |       0        |
|  Inversions  |     15,640      |       867       |     6,103     |      464      |        0       |       0        |
|  **Total**   |   **25,908**    |    **7,946**    |  **14,603**   |   **6,015**   |    **6,033**   |   **X,XXX**    |

These final SVs were then merged with the other datasets and used as input into the VG graph as outlined below.  

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
bcftoo
```

The raw calls initially comprised of:  
|    SV Type   | # FT Called | # KI Called | # FT Filtered | # KI Filtered | FT Curated SVs | KI Curated SVs |
| ------------ | ----------- | ----------- | ------------- | ------------- | -------------- | -------------- |
|  Breakends   |    5,908    |    2,020    |       0       |       0       |        0       |       0        |
|  Deletions   |    6,046    |    3,907    |      664      |      319      |       643      |      319       |
| Duplications |    1,267    |     573     |      58       |       14      |        0       |       0        |
|  Insertions  |      0      |      0      |       0       |       0       |        0       |       0        |
|  Inversions  |   11,812    |     313     |      345      |       2       |        0       |       0        |
|  **Total**   | **25,033**  | **25,063**  |   **1,067**   |    **335**    |     **643**    |    **319**     |

## Validating Filtered SV Calls
[SAMplot](https://github.com/ryanlayer/samplot) and [plotCritic](https://github.com/jbelyeu/PlotCritic) were used to evaluate SV calls from Delly, Smoove and Manta. However, SAMplot is only able to plot Deletions, Duplications and Inversions. In addition, except in a few instances, recovering the full sequences representing inversion and duplication haplotypes from short-read data alone is challenging. Only delection calls were curated for genome graph construction as they generally have clear support (read depth, and split reads) and obvious breakpoints. These are important aspects for accurate genotyping from graphs.  

Given the limitations of short-read data alone to resolve insertions, inversions, and duplications at the haplotype level, we limited our graphs to deletions. We required that deletion calls had to have exact breakpoints, and that they were supported by evidence from both split-read and read depth variation. There was no minimum or maximum size limitations.  

Deletions were curated for each of the three tools independently prior to merging. To generate input files representing heterozygous and homozygous alternate sites for these SVs, we focused on genotypes from each of the tools as they could provide clues as to which samples provided support for SV calls.  
```
for tool in delly manta smoove
    do
    echo "EXTRACTING SITES FROM $tool..."
    bcftools query -i 'SVTYPE!="INS" & GT="het"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF>=8) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}' >> ${tool}/samplot_het_n4.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="het"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF==7) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;}' >> ${tool}/samplot_het_n3.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="het"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF==6) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}' >> ${tool}/samplot_het_n2.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="het"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF==5) print $1"\t"$2"\t"$3"\t"$4"\t"$5;}' >> ${tool}/samplot_het_n1.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="AA"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF>=8) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;} > ${tool}/samplot_homAlt_n4.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="AA"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF==7) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;} > ${tool}/samplot_homAlt_n3.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="AA"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF==6) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;} > ${tool}/samplot_homAlt_n2.tsv
    bcftools query -i 'SVTYPE!="INS" & GT="AA"' -f '%CHROM\t%POS\t%END\t%SVTYPE\t[%SAMPLE\t]\n' ${tool}/02_filteredSVs.vcf | awk '{if (NF==5) print $1"\t"$2"\t"$3"\t"$4"\t"$5;} > ${tool}/samplot_homAlt_n1.tsv
done
```
Then used as input for plotting for 4 fairy tern samples.  
```
for geno in {homAlt,het,ref}
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
                    align/nodup/${samp4}_nodup_autosomes.bam \-t
                    -o ${tool}/samplot_outputs/${geno}/${chrom}_${start}_${end}_${type}.png \
                    -c ${chrom} -s ${start} -e ${end} -t ${type}
        done < ${tool}/samplot_${geno}_sites.tsv
done
```
However, not all SV calls had support for a minimum of 4 individuals. These calls were filtered out (e.g., `awk { if (NR==5) {print $0} } > samplot_het_n1.tsv`) and plotted in a similar manner as above.

[PlotCritic](https://github.com/jbelyeu/PlotCritic) v1.0.1 was used to evaluate whether calls had both read depth and split-read support prior to merging.  
```
plotcritic -p $tool -i ${tool}/kaki_samplot/samplot/ -q "Is this Deletion supported?" -A "y":"Yes" "n":"No" "u":"Uncertain"
```

But SMOOVE doesn't include the reference allele in the output VCF. To correct this, we normalised the SMOOVE output.  
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
Ran into some possible hiccups with `vg giraffe` that are consistent with the kākāpō trials. For instance, `vg call` has an interesting quirk when used after alignment with `giraffe` where variants that are close to one another are 'nested', causing multiallelic calls (see [here](https://www.biostars.org/p/9578818/#9578853) for some details).

Due to this, and the poor mapping quality yielded by `vg giraffe`, we decided to use the `vg map` program instead. This approach will potentially help with nested multiallelics as it theoretically can use the VCF when the `-a` flag is on during graph construction. However, super memory and RAM hungry. To make it work, need to 1) Build graphs for individual chromosomes; 2) Coordinate the node IDs; 3) Index graphs (most memory hungry step); 4) generate chromosome specific `.fq` files for individuals; 5) map reads; and finally 6) call genotypes.

To first construct chromosomes specific graphs. Then we consolidated the node IDs. Note, my code is slightly different from what is suggested by VG Wiki. This is because my chromosome names are not purely numerical and I was RAM limited for running jobs like construction in parallel.  
```
for CHR in CM020{437..458}.1_RagTag
    do
    printf "STARTED CONSTRUCTING GRAPH FOR $CHR AT "
    date
    vg construct -f -S -a -t 16 -R $CHR -r $REF -v vg/DEL_calls.vcf.gz > vg/graphs/${CHR}_graph.vg
done

vg ids -j vg/graphs/CM020*
```
And indexed.
```
for GRAPH in vg/graphs/CM020*_RagTag_graph.vg
    do
    CHR=$(basename $GRAPH _graph.vg)
    printf "STARTED GENERATING XG INDEX FOR $CHR AT "
    date
    vg index -t 16 -x vg/graphs/${CHR}_graph.xg -g vg/graphs/${CHR}_graph.gcsa $GRAPH
done
```
Extracting chromosome specific reads for individuals, mapping and finally genotyping individuals.
```
for BAM in ${DIR}nodup/*_nodup_autosomes.bam
    do
    INDIV=$(basename $BAM _nodup_autosomes.bam)
    while read -r line
        do
        SCF=$(echo $line | awk '{print $1}')
        CHR=$(echo $line | awk '{print $2}')
        printf "STARTED EXTRACTING READS FOR $INDIV CHROMOSOME $CHR AT "
        date
        samtools view -@16 -b $BAM | samtools fastq -@16 -1 vg/reads/${indiv}_chr${CHR}_R1.fq -2 vg/reads/${INDIV}_chr${CHR}_R2.fq
        printf "STARTED ALIGNING READS FOR $INDIV CHR${CHR} AT "
        date
        vg map -t 16 -x vg/graphs/${SCF}_graph.xg -g vg/graphs/${SCF}_graph.gcsa -f vg/reads/${INDIV}_chr${CHR}_R1.fq -f vg/reads/${INDIV}_chr${CHR}_R2.fq > vg/gam/${INDIV}_chr${CHR}.gam
        printf "NOW GENOTYPING $INDIV CHR${CHR} AT "
        date
        vg pack -t 16 -x vg/graphs/${SCF}_graph.xg -g vg/gam/${INDIV}_chr${CHR}.gam -Q 5 -o vg/genotypes/${INDIV}_chr${CHR}.pack
        vg call --genotype-snarls vg/graphs/${SCF}_graph.xg -k vg/genotypes/${INDIV}_chr${CHR}.pack -t 16 -s ${INDIV} > vg/genotypes/${INDIV}_chr${CHR}.vcf
        bgzip vg/genotypes/${INDIV}_chr${CHR}.vcf
        tabix -p vcf vg/genotypes/${INDIV}_chr${CHR}.vcf.gz
    done < vg/chr_list.tsv
done
```

## Mapping quality comparisons
The below is modified from [this](https://gtpb.github.io/CPANG18/pages/toy_examples) tutorial for assessing mapping scores of reads aligned to own assemblies and to the graph.

MapQ scores for reads aligned to genome graph were extracted with:
```
for GAM in vg/gam/AU*.gam
    do
    INDIV=$(basename $GAM .gam)
    vg view -aj ${GAM} | jq -cr '[.name, .mapping_quality] | @tsv' > vg/gam/${indiv}_mapQual.tsv
    samtools view bam/${indiv}_nodup_autosomes.bam | awk '{print $1"\t"$5}' > vg/gam/${indiv}_bamQual.tsv
    join <(sort vg/gam/${indiv}_bamQual.tsv ) <(sort vg/gam/${indiv}_mapQual.tsv ) | awk -v VAR=$INDIV '{print var"\t"$0}' > vg/gam/${indiv}_qual.tsv
done

printf "Individual\tRead ID\tLinear Mapping Q Score\tGenome Graph Q Score\n" > vg/gam/fairy_mapping_scores.tsv
cat vg/gam/*_qual.tsv >> vg/gam/fairy_mapping_scores.tsv

```
Finally, the output of mapping quality for reads aligned to Katie's genome and the genome graph were merged with

```
awk '{ if ($4 < 0) print $1 }' ${indiv}_mapQual_compared.tsv | wc -l # Number of reads that aligned better to self
awk '{ if ($4 == 0) print $1 }' ${indiv}_mapQual_compared.tsv | wc -l # Number of alignments that were the same quality
awk '{ if ($4 > 0) print $1 }' ${indiv}_mapQual_compared.tsv | wc -l # Number of reads that aligned better in the graph
```
Then looked at distribution of mapping quality scores with:
```
awk '{print $2}' DEL/vg_maps/Ariki.tsv | sort -n | awk ' {print $0 (NF<1 ? OFS "NA" : "") }' | uniq -c
```
### Visualising mapping quality
First augmented the input file
```
printf "indiv\tread_id\tscore\tdata\n" > fairy_mapping_scores.tsv

awk '{print $1"\t"$2"\t"$3"\tgenome_graph"} *_mapQual.tsv >> fairy_mapping_scores.tsv
awk '{print $1"\t"$2"\t"$4"\tlinear_genome"} *_bamQual.tsv >> fairy_mapping_scores.tsv
```