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
REF=~/reference/Katie_5kb_ragtag.fa
INPUT=~/reads/
OUTPUT=~/alignments/
for LIB in lib1 lib2 LIC001 LIC002
    do
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
    done
done
```
Files were then merged for each of the datasets, the reference sample (SP01) is the only individual sequenced at two different facilites.  
```
samtools merge -@ 16 \
    -o ${DIR}SP01_merged.bam \
    ${DIR}SP01_lib1.bam ${DIR}SP01_lib2.bam \
    ${DIR}SP01_LIC001.bam ${DIR}SP01_LIC002.bam
```
 
 This `for` loop below was used to merge all other samples. The sequence batches were run simultaneously by changing the suffix of the file (i.e., LIC002, lib2) in the initiation of the `for` loop and for the creation of the `$BASE` variable.  
```
for BAM in ${DIR}*_lib2.bam
    do
    BASE=$(basename $bam _lib2.bam)
    echo "MERGING BAMS FOR $BASE..."
    samtools merge -@64 -o ${DIR}${BASE}_merged.bam ${DIR}${BASE}*.bam
done
```
All alignments were then sorted and PCR duplicates marked and removed using `SAMtools` v.1.16.  
```
for BAM in ${DIR}*_merged.bam
    do
    BASE=$(basename $BAM _merged.bam)
    echo "NSORTING AND FIXING MATE PAIRS FOR ${BASE}\n"
    samtools sort -@64 -n ${BAM} | samtools fixmate -@64 -m -c - ${DIR}${BASE}.fixmate.bam
    samtools sort -@64 -o ${DIR}bam/${BASE}.fixmate.sorted.bam ${DIR}${BASE}.fixmate.bam
    printf "MARKING AND REMOVING PCR DUPLICATES FOR ${BASE} AT "
    date
    samtools markdup -@64 --write-index ${DIR}${BASE}.fixmate.sorted.bam ${DIR}${BASE}_markdup.bam
    samtools markdup -r -@64 --write-index ${DIR}${BASE}.fixmate.sorted.bam ${DIR}${BASE}_nodup.bam
    printf "\nFINISHED MARKING AND REMOVING DUPLICATES FOR ${BASE} AT "
    date
done
```
Once BAMs were processed, autosomal chromosomes were extracted for population analyses. Mean mapping quality and alignment depth was then estimated using [QualiMap](http://qualimap.conesalab.org/) v2.3.  
```
cut -f1,2 TKatie_5kb_ragtag.fa.fai | grep "CM020" | grep -v CM020462.1_RagTag | grep -v CM020463.1_RagTag | awk '{print $1"\t1\t"$2} > Katie_autosomes.bed

for SAMP in ${DIR}*_markdup.bam
    do
    BASE=$(basename $SAMP _markdup.bam)
    printf "EXTRACTING AUTOSOMES FOR ${BASE} AT "
    date
    samtools view -@64 -L Katies_genome/Katie_autosomes.bed \
        --write-index -o ${DIR}markdup/${BASE}_markdup_autosomes.bam ${SAMP}
    samtools view -@64 -L Katies_genome/Katie_autosomes.bed \
        --write-index -o ${DIR}nodup/${BASE}_nodup_autosomes.bam ${DIR}nodup/${BASE}_nodup.bam
    printf "FINISHED CONVERTING AND INDEXING FILES FOR ${BASE}, NOW ESTIMATING STATS AT "
    date
    mosdepth --threads 24 --fast-mode --by 50 ${DIR}mosdepth/${BASE}_markdup_autosomes ${SAMP}
    mosdepth --threads 24 --fast-mode --by 50 ${DIR}mosdepth/${BASE}_nodup_autosomes ${SAMP}
    qualimap bamqc -bam ${DIR}markdup/${BASE}_markdup_autosomes.bam -outdir ${DIR}qualimap/${BASE}_markdup_autosomes -outformat PDF:HTML --java-mem-size=3G
    qualimap bamqc -bam ${DIR}nodup/${BASE}_nodup_autosomes.bam -outdir ${DIR}qualimap/${BASE}_nodup_autosomes -outformat PDF:HTML --java-mem-size=3G
    echo "FINISHED PROCESSING ${BASE}!"
done
```
For ease of comparisons, mosdepth outputs were also plotted with:
```
python ~/anaconda3/envs/mosdepth/scripts/plot-dist.py ${data}nodup_bam_stats/*.global.dist.txt
```

Three additional scaffolds `CM020459.1`, `CM020460.1`, and `CM02061.1` were excluded from downstream analyses as these scaffolds consistently had read depths much higher than expected and likely represents improperly assembled representations of these chromosomes. A new bedfile, cleverly named `Katie_autosomes2.bed`, was used to remove these scaffolds prior to analyses.  
```
for SAMP in ${DIR}*_autosomes.bam
    do
    BASE=$(basename $SAMP _autosomes.bam)
    printf "EXTRACTING AUTOSOMES FOR ${BASE} AT "
    date
    samtools view -@64 -L Katies_genome/Katie_autosomes2.bed \
        --write-index -o ${DIR}markdup/${BASE}_markdup_autosomes2.bam ${SAMP}
    samtools view -@64 -L Katies_genome/Katie_autosomes2.bed \
        --write-index -o ${DIR}nodup/${BASE}_nodup_autosomes2.bam ${DIR}nodup/${BASE}_nodup.bam
    echo "FINISHED PROCESSING ${BASE}!"
done
```
 The BAM files originally filtered for autosomomal scaffolds were subsequently overwritten with the newly filtered bam files. Thus all files denoted as `*_{markdup,nodup}_autosomes.bam` contain only 22 scaffolds.  