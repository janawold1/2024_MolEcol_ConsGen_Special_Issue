# Tara iti Genome Assembly
Two tara iti clutchmates (1 male, 1 female) passed away as a result of a storm in the 2018/2019 breeding season. The female chick was sequenced with both short- and long-read technology. Short-reads, used for assembly polishing and assessing sequencing batch-effects, were generated using an Illumina NovaSeq 6000 at the [Institute of Clinical Molecular Biology](https://www.ikmb.uni-kiel.de/resources/sequencing/whole-genome-de-novo-sequencing) (Kiel, Germany) and Livestock Improvement Corp (Hamilton, NZ). Long-read sequencing was performed using a ONT PromethION (R10 flowcells, 5000Hz, LSK114 kit) at [Bragato Research Institute](https://bri.co.nz/) (Lincoln, NZ).  

The short-read data were processed as per [01_pop_sequence_QC_and_alignment.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/01_pop_sequence_QC_and_alignment.md), while the ONT data was prepared as per below.

## Basecalling and Read Trimming
One library preparation was sufficient for all ONT sequencing. The PromethION flow cell was washed and loaded four times. All ONT `.pod5` files were basecalled using [Dorado](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) v0.5.0 and the `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` model. First, stereo basecalling was performed. Here, `$pod5_dir` represents the directory where raw `.pod5` files were stored, `$indiv` represents the individual ID of the reference genome and `$model` represents the Dorado basecalling model used.  
```
model=dna_r10.4.1_e8.2_400bps_sup@v4.3.0
#dorado download --model ${model}

dorado basecaller --min-qscore 10 --emit-moves ${model} \
  ${pod5_dir}${indiv}/${lib}/ > ${out}${indiv}/${indiv}_${lib}_R10_moves.sam --device 'cuda:all'
```

`.sam` outputs from basecalling were then indexed using [SAMtools](https://github.com/samtools/samtools) v1.16 and [duplex tools](https://github.com/nanoporetech/duplex-tools) v0.2.20 was used to identify duplex pairs.  
```
samtools index -@ 16 ${out}${indiv}/${indiv}_${lib}_R10_moves.sam

duplex_tools pair --output_dir ${out}${indiv}/${lib}_duplex/ \
  --prefix ${indiv}_${lib}_R10 --verbose --threads 32 \
  ${out}${indiv}/${indiv}_${lib}_R10_moves.sam

duplex_tools split_pairs --debug --threads 32 \
  ${out}${indiv}/${indiv}_${lib}_R10_moves.sam \
  ${pod5_dir}${indiv}/${lib} \
  ${dup_pod}${indiv}/${lib}
```
Duplex basecalling for paired and split-paired reads was performed with Dorado. Prior to performing this step, the outputs from duplex tools were renamed to have the structure of `$indiv_$lib_*`. This denotes the genome ID (`$indiv`) and ONT run (`$lib`).  
```
dorado duplex ${model} ${pod5_dir}${indiv}/${lib}/ \
  --pairs ${out}${indiv}/${lib}_pair_ids_filtered.txt > ${out}${indiv}/${indiv}_${lib}_duplex_orig.sam

cat ${split}${indiv}/${lib}/*_split_duplex_pair_ids.txt > ${out}${indiv}/${lib}_split_duplex_pair_ids.txt

dorado duplex ${model} ${split}${indiv}/${lib}/ \
  --pairs ${out}${indiv}/${lib}_split_duplex_pair_ids.txt > ${out}${indiv}/${indiv}_${lib}_duplex_splitduplex.sam
```
### Read Trimming 
Reads from both simplex and duplex calls were converted to FastQ format with SAMtools.
```
samtools fastq ${out}${indiv}/${indiv}_${lib}_R10_moves.sam > ${out}${indiv}/${indiv}.fastq
samtools fastq ${out}${indiv}/${indiv}_${lib}_duplex_orig.sam >> ${out}${indiv}/${indiv}.fastq
samtools fastq ${out}${indiv}/${indiv}_${lib}_duplex_splitduplex.sam >> ${out}${indiv}/${indiv}.fastq
```
Adapter trimming was not performed as Dorado now performs trimming as part of the basecalling process. [Chopper](https://github.com/wdecoster/chopper) v0.5.0, which is the Rust implementation of NanoFilt & NanoLyse, was used to trim reads to a minimum Q score of 20 and 3 minimum lengths, 1kb, 5kb and 10kb.  
```
cat ${out}${indiv}/${indiv}.fastq | chopper -q 20 -l 1000 > ${out}${indiv}/${indiv}_dorado_q20_1kb.fq
cat ${out}${indiv}/${indiv}.fastq | chopper -q 20 -l 5000 > ${out}${indiv}/${indiv}_dorado_q20_5kb.fq
cat ${out}${indiv}/${indiv}.fastq | chopper -q 20 -l 10000 > ${out}${indiv}/${indiv}_dorado_q20_10kb.fq
```
## Initial Genome Assembly
Initial genome assemblies using reads trimmed to a minimum Q-score of 20 and a minimium length of either 1, 5 or 10kb were performed using [FLYE v2.9.1](https://github.com/fenderglass/Flye).  
```
flye --nano-raw ${reads}${indiv}/${indiv}_q20_5kb.fastq \
  --out-dir ${out}${indiv}/dorado_q20_5kb_simple/ --genome-size 1.2g \
  --threads 24 --debug
done
```
**Table 1. Initial assembly statistics:**
| Read Inputs | Estimated Depth (FLYE) | # Scaffolds | N50 (Mbp) | L50 | Largest Scaffold Size (Mbp) |
|:-----------:|:----------------------:|:-----------:|:---------:|:---:|:---------------------------:|
|   Q20, 1kb  |           43           |     692     |  27.4     | 13  |             74.7            |
|   Q20, 5kb  |           38           |     497     |  24.3     | 13  |             77.6            |
|  Q20, 10kb  |           29           |     486     |  24.3     | 14  |             83.6            |

## Polishing and Scaffolding
A combination of [Racon v1.5.0](https://github.com/isovic/racon) and [longstitch v1.0.4](https://github.com/bcgsc/LongStitch) were used to polish and scaffold the genome. The below script aligns reads to the draft assembly using [Minimap2 v2.24](https://github.com/lh3/minimap2) for polishing with RACON and scaffolding with Longstitch. Polishing and scaffolding was repeated for two rounds.  
```
ml purge
ml load minimap2/2.24-GCC-11.3.0
ml load Racon/1.5.0-GCC-11.3.0
ml load LongStitch/1.0.4-Miniconda3

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/
reads=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_basecalled_reads/
asm=q20_5kb

for indiv in SP01
        do
        mkdir ${dir}${indiv}/${asm}_polish
        printf "\nRUNNING RACON USING Q20 5kb ASSEMBLY FOR ${indiv} AT "
        date
        minimap2 -t 32 -ax map-ont ${dir}${indiv}/dorado_${asm}_simple/${indiv}_${asm}_flye_simple.fa \
                ${reads}${indiv}/${indiv}_${asm}.fastq > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple.sam
        printf "FINISHED RUNNING MINIMAP STARTING RACON FOR ${indiv} AT "
        date
        racon -t 32 ${reads}${indiv}/${indiv}_${asm}.fastq \
                ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple.sam \
                ${dir}${indiv}/dorado_${asm}_simple/${indiv}_${asm}_flye_simple.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple_racon1.fa
        printf "FINISHED RACON, STARTING LONGSTITCH AT "
        date
        mkdir ${dir}${indiv}/${asm}_longstitch1
        cp ${reads}${indiv}/${indiv}_${asm}.fastq ${dir}${indiv}/${asm}_longstitch1/
        cp ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple_racon1.fa ${dir}${indiv}/${asm}_longstitch1/
        longstitch run -C ${dir}${indiv}/${asm}_longstitch1 \
                draft=${indiv}_${asm}_flye_simple_racon1 reads=${indiv}_${asm} t=32 G=1.2g --debug
        cat ${dir}${indiv}/${asm}_longstitch1/${indiv}_${asm}*.ntLink.scaffolds.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.fa
        printf "FINISHED ROUND 1 LONGSTITCH, STARTING ROUND 2 MINIMAP AT "
        date
        minimap2 -t 64 -ax map-ont ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.fa \
                ${reads}${indiv}/${indiv}_${asm}.fastq > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.sam
        printf "FINISHED ROUND 2 MINIMAP, STARTING ROUND 2 RACON AT "
        date
        racon -t 64 ${reads}${indiv}/${indiv}_${asm}.fastq \
                ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.sam \
                ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1_racon2.fa
        printf "FINISHED ROUND 2 RACON, STARTING ROUND 2 LONGSTITCH AT "
        date
        mkdir ${dir}${indiv}/${asm}_longstitch2
        mv ${dir}${indiv}/${asm}_longstitch1/${indiv}_${asm}.fastq ${dir}${indiv}/${asm}_longstitch2/
        cp ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1_racon2.fa ${dir}${indiv}/${asm}_longstitch2/
        longstitch run -C ${dir}${indiv}/${asm}_longstitch2 \
                draft=${indiv}_${asm}_longstitch1_racon2 reads=${indiv}_${asm} t=32 G=1.2g --debug
        cat ${dir}${indiv}/${asm}_longstitch2/${indiv}_${asm}*.ntLink.scaffolds.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch2.fa
        printf "\nFINISHED RUNNING 2 ROUNDS OF POLISHING AND SCAFFOLDING OF Q20 5kb ASSEMBLY FOR ${indiv} AT "
        date
done
```
Finally, assembly quality was assessed using [BUSCO](https://busco.ezlab.org/) v5.4.7 to assess how complete it may be, and [Quast](https://github.com/ablab/quast) v5.2.0 to assess contiguity (Table 2).  

**Table 2. Polished assembly summary statistics:**
| Read Inputs | % Complete BUSCO | % Missing BUSCO | # Scaffolds | N50 (Mbp) | L50 | Largest Scaffold Size (Mbp) | # N's per 100 kbp |
|:-----------:|:----------------:|:---------------:|:-----------:|:---------:|:---:|:---------------------------:|:-----------------:|
|   Q20, 1kb  |       97.7       |       1.9       |     433     |   48.3    |  8  |            100.3            |       25.9        |
|   Q20, 5kb  |       97.7       |       1.9       |     298     |   41.3    |  9  |            103.4            |       17.3        |
|  Q20, 10kb  |       97.5       |       2.1       |     283     |   35.4    | 10  |             86.9            |        6.2        |


[D-GENIES](https://dgenies.toulouse.inra.fr/) (figure below) to visualise structural differences with a HQ assembly for the common tern ([*Sterna hirundo*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009819605.1/)).  

In the end, the genome assembly generated using a minimum read length of 5kb was used for population read alignment and analyses as it represented a good balance of coverage (Table 1) and contiguity (Table 2).

<figure>
        <div style="text-align: center;">
        <img src="https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/Figures/Katie_q20_5kb_longstitch2_to_CommonTern.png"
             alt="Tara iti genome aligned against the common tern genome for comparing synteny and contiguity"
             width="600" height="600">
        </div>
        <figcaption>Tara iti draft assembly mapped to the VGP assembly for Common tern. In this dotplot, the tara iti assembly is represented along the y-axis and the common tern assembly is along the x-axis.</figcaption>
</figure>

Unsurprisingly, the *de novo* tara iti assembly was not chromosomally resolved. Because common terns and fairy terns are relatively related species and demonstrate high synteny, we used the common tern as a reference to scaffold the tara iti assembly with [RagTag](https://github.com/malonge/RagTag) v2.1.0 to maximise our ability to call structural variants.  
```
ragtag.py scaffold -o reference/SP01_ragtag/ reference/common_tern.fasta reference/SP01.fasta
```
**Table 3. Final assembly summary statistics:**
| Read Inputs | % Complete BUSCO | % Missing BUSCO | # Scaffolds | N50 (Mbp) | L50 | Largest Scaffold Size (Mbp) | # N's per 100 kbp |
|:-----------:|:----------------:|:---------------:|:-----------:|:---------:|:---:|:---------------------------:|:-----------------:|
|   Q20, 5kb  |       97.7       |       1.9       |     137     |    84.9   |  5  |            219.3            |       18.7        |

# Using SeqKit Properly
```
seqkit grep -p  "^[^]/|/(\D )(\w+) " -f autosome_names.txt ${REF} > ${OUT_REF}
```