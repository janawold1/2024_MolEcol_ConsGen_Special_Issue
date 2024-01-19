# Illumina Read Processing and Alignment
## Initial QC and Trimming
Initial sequencing results were visualised with FastQCvX.X and MultiQCvX.X. Reads were trimmed using [TrimGalorevX.X](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).  
```
ref=/dir/to/reference.fasta
dir=/dir/to/reads/
reports=/dir/to/fastq_out/
trim_out=/dir/to/trimmmed_reads/

for samp in ${dir}*R1.fastq.gz
do
base=$(basename ${samp} _R1.fastq.gz)
echo "Running Trim_galore for ${base}..."
trim_galore --paired \
    --nextseq 28 \
    --2colour 20 \
    --cores 8 \
    --fastqc_args "--nogroup --outdir ${reports} --threads 32" \
    --length 50 \
    --output_dir ${trim_out} \
    --clip_R1 20 \
    --clip_R2 20 \
    --three_prime_clip_R1 5 \
    --three_prime_clip_R2 5 \
    --retain_unpaired \
    --length_1 55 \
    --length_2 55 \
    ${dir}${base}_R1.fastq.gz \
    ${dir}${base}_R2.fastq.gz
done
```
## Alignment
For population analyses, Illumina short-reads were aligned, PCR-duplicates removed, and alignment statistics estimated in the same manner for both Australian fairy tern and tara iti (outlined below). All manipulation of alignment files was performed with [SAMtools v1.16](https://www.htslib.org/).  
```
#!/bin/bash -e
ref=/media/jana/BigData/tara_iti_publication/reference/Katie_5kb_ragtag.fa
reads=/media/jana/BigData/tara_iti_publication/reads/
out=/media/jana/BigData/tara_iti_publication/alignments/
for LIB in lib1 lib2 LIC001 LIC002
    for SAMP in ${INPUT}${LIB}/*_R1.fq.gz
        do
        BASE=$(basename ${SAMP} _${LIB}_R1.fq.gz)
        INFOLINE=$(zcat ${SAMP} | head -n 1)
        INSTRUMENT=`echo ${INFOLINE} | cut -d ':' -f1`
        INSTRUMENTRUN=`echo ${INFOLINE} | cut -d ':' -f2`
        FLOWCELL=`echo ${INFOLINE} | cut -d ':' -f3`
        LANE=`echo ${INFOLINE} | cut -d ':' -f4`
        INDEX=`echo ${INFOLINE} | cut -d ':' -f10`

        RGID="ID:${INSTRUMENT}_${INSTRUMENTRUN}_${FLOWCELL}_${LANE}_${INDEX}"
        RGPL="PL:Illumina"
        RGPU="PU:${FLOWCELL}.${LANE}"
        RGLB="LB:${BASE}_${LIB}"
        RGSM="SM:${BASE}"

        printf "BEGAN ALIGNING READS FOR ID:${BASE}, LIB:${LIB} AT "
        date
        bwa mem -M -R @RG'\t'$RGID'\t'$RGPL'\t'$RGPU'\t'$RGLB'\t'$RGSM -t 46 $REF $SAMP ${INPUT}${LIB}/${BASE}_${LIB}_R2.fq.gz | samtools sort -@ 46 -o ${OUTPUT}${BASE}.bam
        printf "FINISHED ALIGNING READS FOR ID:${BASE}, LIB:${LIB} AT "
        date
        wait
        done
done
```
Files were then merged for each of the datasets, the reference sample (SP01) is the only individual sequenced at two different facilites.  
```
samtools merge -@ 16 \
    -o ${out}bam/SP01_merged.bam \
    ${out}bam/SP01_lib1.bam ${out}bam/SP01_lib2.bam \
    ${out}bam/SP01_LIC001.bam \
    ${out}bam/SP01_LIC002.bam
```
 
 This `for` loop below was used to merge all other samples.  
```
for lib in lib2 LIC002
    do
    for bam in ${out}bam/*_${lib}.bam
        do
        base=$(basename $bam _lib2.bam)
        base=$(basename $base _LIC002.bam)
        echo "MERGING BAMS FOR $base..."
        samtools merge -@ 16 -o ${out}bam/${base}_merged.bam ${out}bam/${base}*.bam
    done
done
```

All alignments to all three reference assemblies were sorted and PCR duplicates removed using `SAMtools`.  
```
for bam in ${dir}bam/*_merged.bam
    do
    base=$(basename $bam _merged.bam)
    echo "NSORTING AND FIXING MATE PAIRS FOR ${base}..."
    samtools sort -@16 -n ${bam} | samtools fixmate -@16 -m -c - ${dir}bam/${base}.fixmate.bam
    samtools sort -@16 -o ${dir}bam/${base}.fixmate.sorted.bam ${dir}bam/${base}.fixmate.bam
    echo "MARKING PCR DUPLICATES FOR ${base}..."
    samtools markdup -O CRAM -@16 --reference ${ref} --write-index \
        ${dir}bam/${base}.fixmate.sorted.bam ${dir}cram/${base}_nodup.cram
done
```
Once BAMs were processed, files were converted to CRAM format and autosomal chromosomes were extracted for population analyses. Mean mapping quality and alignment depth was then estimated for comparisons between the three datasets and files converted to CRAM format to save disk space.  
```
cut -f1,2 TI_as_CT.fasta.gz.fai | grep "CM020" | grep -v CM020462.1_RagTags | grep -v CM020463.1_RagTags | awk '{print $1"\t1\t"$2} > TI_as_CT_autosomes.txt
cut -f1,2 common_tern.fasta.gz.fai grep "CM020" | grep -v CM020462.1 | grep -v CM020463.1 | awk '{print $1"\t1\t"$2} > common_tern_autosomes.txt

for samp in ${dir}nodup_bam/*_nodup.bam
    do
    base=$(basename $samp _nodup.bam)
    echo "EXTRACTING AUTOSOMES FOR ${base}..."
    samtools view -@16 -L reference/TI_as_CT_autosomes.bed -T ${ref} -O CRAM \
        --write-index -o ${dir}cram/${base}_nodup_autosomes.cram ${dir}cram/${base}_nodup.cram
    echo "FINISHED CONVERTING AND INDEXING FILES FOR ${base}, NOW ESTIMATING STATS..."
    samtools stats -@ 16 ${dir}cram/${base}_nodup.cram > ${dir}align_stats/${base}_nodup.stats
    samtools stats -@ 16 ${dir}cram/${base}_nodup_autosomes.cram > ${dir}align_stats/${base}_nodup_autosomes.stats
    echo "FINISHED PROCESSING ${base}!"
done
```