# Illumina Read Processing and Alignment
## Initial QC and Trimming
Initial sequencing results were visualised with FastQCvX.X and MultiQCvX.X. Reads were trimmed using [TrimGalorevX.X]().  

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
The *de novo* tara iti assembly was highly fragmented. Because common terns and fairy terns are congeneric species, we used the common tern as a reference to scaffold the tara iti assembly in an attempt to retain information for polarizing the site frequency spectrum (SFS) and maximise our ability to call structural variants. This was accomplished by whole-genome alignment with Minimap as implemented in [RagTag](https://github.com/malonge/RagTag).

Overall, reads were aligned to three reference genomes: 1) the common tern genome assembly; 2) the *de novo* tara iti assembly; and 3) the tara iti assembly scaffolded using the common tern as reference. These three datasets were aligned, PCR-duplicates removed and alignment statistics estimated in the same manner (outlined below).  
```
ragtag.py scaffold -o reference/ragtag_TI_as_CT/ reference/common_tern.fasta reference/tara_iti_masurca_v2.fasta
```
Reads were aligned to all three genomes for the common tern, tara iti and scaffolded tara iti datasets.  
```
dir=/path/to/reads/
out=/path/to/alignments/

echo "Indexing the reference genome $ref"
time bwa index $ref

echo "Now beginning alignments...."
for samp in ${dir}*_R1.fq.gz 
    do
    #Remember to be explicit with file location
    base=$(basename $samp _R1.fq.gz)
    infoline=$(zcat ${samp} | head -n 1)
    instrument=`echo ${infoline} | cut -d ':' -f1`
    instrumentrun=`echo ${infoline} | cut -d ':' -f2`
    flowcell=`echo ${infoline} | cut -d ':' -f3`
    lane=`echo ${infoline} | cut -d ':' -f4`
    index=`echo ${infoline} | cut -d ':' -f10`
    
    #Now to incorporate this information into the alignment
    rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
    rgpl="PL:Illumina"
    rgpu="PU:${flowcell}_${lane}"
    rglb="LB:${base}_${lane}"
    rgsm="SM:${base}"

    #Be explicit with file location for read 2 and the sam file output
    echo "Aligning reads for $base" 
    time bwa mem -M \
        -R @RG'\t'${rgid}'\t'${rgpl}'\t'${rgpu}'\t'${rglb}'\t'${rgsm} \
        -t 64 ${ref} ${samp} ${dir}${base}_R2.fq.gz | samtools view -b -@ 64 > ${out}bam/${base}.bam
    time samtools sort -@ 64 -o ${out}bam/${base}.sorted.bam ${out}bam/${base}.bam
done
```
## Remove PCR duplicates and Initial alignment Stats
All alignments to all three reference assemblies were sorted and PCR duplicates removed using `SAMtools`. Mean mapping quality and alignment depth was then estimated for comparisons between the three datasets.  
```
for bam in ${dir}bam/*.bam
    do
    base=$(basename $bam .bam)
    echo "Removing PCR duplicates from ${base}..."
    samtools sort -@16 -n -o ${dir}bam/${base}.nsorted.bam ${bam}
    samtools fixmate -@16 -r -m -c ${dir}bam/${base}.nsorted.bam ${dir}bam/${base}.fixmate.bam
    samtools sort -@16 -o ${dir}bam/${base}.fixmate.sorted.bam ${dir}bam/${base}.fixmate.bam
    samtools markdup -@16 ${dir}bam/${base}.fixmate.sorted.bam ${dir}bam/${base}_nodup.bam
    samtools stats -@16 ${dir}bam/${base}_nodup.bam > ${dir}bam/${base}_nodup.stats
    qualimap bamqc -bam ${dir}bam/${base}_nodup.bam -nw 10000 -nt 32 -c -outdir ${dir}bam/${base}.graphmap --java-mem-size=8G
done
```