**************************************************************************************
*   Genomic signatures of inbreeding and mutation load in a threatened rattlesnake   *
*                                        by                                          *
*                         Alexander Ochoa and H. Lisle Gibbs                         *
**************************************************************************************

#Please contact the authors (alexander.ochoa@ucf.edu; gibbs.128@osu.edu) for questions


-----------------------------------Trim Sequences-------------------------------------

#(Step 1) Trimming
trim_galore --paired --stringency 5 --nextera --nextseq 20 --length 100 --trim-n --fastqc R1.fastq R2.fastq
trim_galore --paired --stringency 5 --quality 20 --length 100 --trim-n --fastqc R1.fastq R2.fastq
###Used Trim Galore v.0.4.5 (https://github.com/FelixKrueger/TrimGalore) 
###'--paired' indicates paired-end reads
###'--stringency' indicates overlap INT with adapter sequence required to trim a read
###'--nextera' trims nextera (and Illumina) adapters from the reads
###'--nextseq' eliminates low quality bases (e.g., <20) on either side of the read and clips poly-G tails (only appicable to 2-dye systems)
###'--length' discards reads (and their partners) shorter than length INT
###'--trim-n' removes Ns on either side of the read
###'--fastqc' outputs FastQC report of your trimmed files
###The second command line is for samples sequenced with 4-dye systems


-----------------------------Mapping/Redundancy Removal-------------------------------

#(Step 2) Create index files for reference sequence
bwa index reference.fasta
samtools faidx reference.fasta
###Used BWA v.0.7.15 (https://github.com/lh3/bwa/releases)
###Used SAMtools v.1.3.1 (https://sourceforge.net/projects/samtools/files/samtools/)

#(Step 3) Mapping
bwa mem -M -t 15 -R "@RG\tID:1\tPL:illumina\tPU:GibbsLab\tLB:n/a\tSM:sampleID\tCN:OSU" -v 3 reference.fasta R1.trimmed.fastq R2.trimmed.fastq > output.sam
###'mem' uses an algorithm designed for short-read sequences
###'-t' inputs an INT number of threads
###'-R' links read group headers to the reads
###'-v' controls the verbose level of the output, 'normal messages' in this case
###Change 'sampleID' accordingly

#(Step 4) Convert to binary format
samtools view -Shu output.sam > output.bam
###'-Shu' includes the header to uncompressed bam output

#(Step 5) Sort BAM file
samtools sort output.bam -o output.sorted.bam
###This step arranges the mapped/unmapped reads with respect to scaffold name (as in the index files) and base position

#(Step 6) Remove redundancies
samtools rmdup output.sorted.bam - | samtools view -b -f 0x0002 -F 0x0004 -F 0x0008 - > output.sorted.rmdup.bam
###This step removes potential PCR duplicates, unmapped reads, and mapped reads with unmapped partners


----------------------------------Indel Realignment-----------------------------------

#(Step 7) Index reference sequence using Picard
java -jar picard.jar CreateSequenceDictionary REFERENCE=reference.fasta OUTPUT=reference.dict
###Download the picard.jar file from https://github.com/broadinstitute/picard/releases/tag/2.23.4

#(Step 8) Index BAM files using SAMtools
ls *.sorted.rmdup.bam | xargs -n1 -P5 samtools index

#(Step 9) Get coordinates for indel regions using GATK
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I output.sorted.rmdup.bam -o output.sorted.rmdup.intervals
###Used GATK v.4.0.1.2 (https://github.com/broadinstitute/gatk/releases)

#(Step 10) Realign indels using GATK
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R reference.fasta -I output.sorted.rmdup.bam -targetIntervals output.sorted.rmdup.intervals -o output.sorted.rmdup.realigned.bam


---------------------------------Coverage Estimates-----------------------------------

#(Step 11) Estimate coverage from BAM files using BEDTools
cut -f 1,2 reference.fasta.fai > bedtools_genome.txt
bedtools genomecov -ibam output.sorted.rmdup.realigned.bam -g bedtools_genome.txt | awk '$1=="genome" {tot += $2*$3; reads += $3} END {print tot/reads}' > output.coverage.txt
###Used BEDTools v.2.19.1 (https://bedtools.readthedocs.io/en/latest/)


-------------------------------------Subsampling--------------------------------------

#(Step 12) Subsample BAM files by depth using SAMtools
samtools view -s 1.00000000000 -b output.sorted.rmdup.realigned.bam > output.sorted.rmdup.realigned.sub.bam
###The '-s' parameter represents the fraction of mapped reads to be retained from the BAM file

#(Step 13) Subsample BAM files by scaffold using SAMtools
samtools view -bh output.sorted.rmdup.realigned.sub.bam jcf7180000038883 > output.sorted.rmdup.realigned.sub.scaff.bam
###As an example, only reads mapped to scaffold jcf7180000038883 will be retained
###Also, remove all other scaffolds from the header of the BAM file and from the FASTA reference assembly


-------------------------------Whole-Genome Alignment---------------------------------

#(Step 14) Run MashMap
mashmap -t 8 -r Cviridis.fasta -q reference.fasta -s 1000 --pi 95 -f one-to-one -o map1000.95.txt
generateDotPlot png large map1000.95.txt
###Used MashMap v.2.0 (https://github.com/marbl/MashMap)
###'-t' inputs an INT number of threads
###'-r' indicates the reference sequence (C. viridis, in this case)
###'-q' indicates the query sequence (S. catenatus, in this case)
###'-s' indicates the window size
###'--pi' indicates the minimum sequence identity
###This step is intended to identify and retain only S. catenatus autosomal scaffolds for downstream analyses (see Step 13)


-----------------------------Genetic Diversity Estimates------------------------------

#(Step 15) Run ANGSD
angsd -bam CEBO.filelist -minMapQ 20 -minQ 20 -doSaf 1 -anc reference.fasta -GL 2 -P 24 -out CEBO
realSFS CEBO.saf.idx -P 24 -fold 1 > CEBO.sfs
realSFS saf2theta CEBO.saf.idx -sfs CEBO.sfs -outname CEBO -fold 1
thetaStat do_stat CEBO.thetas.idx
###Used ANGSD v.0.930 (http://www.popgen.dk/angsd/index.php/ANGSD)
###'-bam' indicates a list of BAM files from a population (CEBO, in this case)
###'-minMapQ' indicates minimum mapping quality
###'-minQ' indicates minimum base quality
###'-doSaf 1' calculates the site allele frequency likelihood assuming HWE
###'-anc' indicates the reference sequence
###'-GL 2' uses GATK genotype likelihoods
###'-P' inputs an INT number of threads
###'-out' indicates the output file name
###'-fold 1' calculates the folded site frequency spectrum
###FASTA and BAM files were composed of autosomal scaffolds >2 Mb in length
###BAM files were subsampled to 5x sequence depth
###This step will estimate theta, pi, and D for each population

------------------------------------ROH Estimates-------------------------------------

#(Step 16) Run ROHan
rohan --size 50000 --rohmu 5e-4 -t 16 -o ROHout reference.fasta output.sorted.rmdup.realigned.sub.scaff.bam
###Used ROHan v.1.0 (https://github.com/grenaud/ROHan)
###'--size' indicates the window size
###'--rohmu' indicates the heterozygosity rate
###'-t' inputs an INT number of threads
###'-o' indicates the output file name
###FASTA and BAM files were composed of autosomal scaffolds >2 Mb in length
###BAM files were subsampled to 5x sequence depth


---------------------------------ROH-Size Clustering----------------------------------

#(Step 17) Run R
###Crate a one-column file (e.g., scat_rohs.txt) with ROH sizes across all S. catenatus samples
###These ROH sizes were obtained from Step 16 after excluding microchromosomal scaffolds
library(mclust)
scat_ROH<-read.csv("~/scat_rohs.txt", header = F)
test<-Mclust(scat_ROH)
scatbic<-mclustBIC(scat_ROH)
summary(scatbic, parameter=TRUE)
plot(scatbic)
summary(test, parameter=TRUE)

x<-seq(0,8000000, length=1000)
scag1<-dnorm(x, mean=test$parameters$mean[1], sd=sqrt(test$parameters$variance$sigmasq[1]))
scag2<-dnorm(x, mean=test$parameters$mean[2], sd=sqrt(test$parameters$variance$sigmasq[2]))
scag3<-dnorm(x, mean=test$parameters$mean[3], sd=sqrt(test$parameters$variance$sigmasq[3]))
hist(test$data, breaks=132, main="Distribution of ROH for inbred S. catenatus", xlab="Number of base pairs")
lines(x, (scag1*1800000000), col="red")
lines(x, (scag2*1800000000), col="blue")
lines(x, (scag3*1800000000), col="green")
legend("topright", title="Gaussian groups", legend=c("Short","Medium","Long"), lty=c(1,1,1), col=c("red", "blue", "green"))

hist(test$data, breaks=132, main="Distribution of ROH for S. catenatus", xlab="Number of base pairs")
abline(v=200000, lty=5, col="red")
abline(v=700000, lty=4, col="blue")
legend("topright", title="ROH boundries", legend=c("Short-Medium break","Medium-Long break"), lty=c(5,4), col=c("red", "blue"))

plot(x,scag1, col="red", type="l", main="ROH length distributions for S. catenatus",
     xlab="Number of base pairs", ylab="", yaxt="n")
lines(x, (scag2), col="blue")
lines(x, (scag3), col="green")
legend("topright", title="Gaussian groups", legend=c("Short","Medium","Long"), lty=c(1,1,1), col=c("red", "blue", "green"))


--------------------------------Relatedness Estimates---------------------------------

#(Step 18) Run ANGSD
angsd -b CEBO.filelist -minMapQ 20 -minQ 20 -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -out CEBO
###'-b' indicates a list of BAM files from a population (CEBO in this case)
###'-gl' and '-doGlf' uses and outputs genotype likelihoods
###'-domajorminor' and '-domaf' outputs minor allele frequencies
###'-minmaf' indicates the minor allele frequency threshold
###'-snp_pval' indicates the P-value threshold for SNPs
###BAM files were composed of autosomal scaffolds >2 Mb in length
###BAM files were subsampled to 5x sequence depth

#(Step 19) Obtain pairwise rxy estimates
zcat CEBO.mafs.gz | cut -f5 |sed 1d > CEBO.freq
ngsRelate -g CEBO.glf.gz -n 9 -f CEBO.freq -O CEBO.res
###Used NgsRelate v.2 (https://github.com/ANGSD/NgsRelate)
###'-g' indicates the genotype likelihoods file
###'-n' indicates number of samples from that population (CEBO, in this case)
###'-f' indicates the minor allele frequencies file


-------------------------------Demographic Simulations--------------------------------

#(Step 20) Annotate repeated elements in the reference sequence
RepeatMasker -pa 16 -gff -xsmall reference.fasta
###Used RepeatMasker v.4.0.7 (https://www.repeatmasker.org/)
###'-pa' inputs an INT number of threads
###'-gff' outputs annotations in GFF format
###'-xsmall' returns repeated elements in lower case
###FASTA file was composed of autosomal scaffolds of all lengths (not just those >2 Mb)

#(Step 21) Convert the resulting GFF file to BED format
bedtools sort -i reference.repeat.gff > reference.repeat.bed
bedtools merge -i reference.repeat.bed > reference.repeat.merged.bed
###Make sure the BED file is 0-based; if not, this can be done manually for each entry of the reference.repeat.bed file

#(Step 22) Call and filter SNPs in an individual sample
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fasta -I output.sorted.rmdup.realigned.bam -o output.sorted.rmdup.realigned.g.vcf -ERC GVCF --disableOptimizations --forceActive --dontTrimActiveRegions
java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R reference.fasta -V output.sorted.rmdup.realigned.g.vcf -o output.sorted.rmdup.realigned.vcf
java -jar GenomeAnalysisTK.jar --analysis_type SelectVariants -R reference.fasta -V output.sorted.rmdup.realigned.vcf --restrictAllelesTo BIALLELIC -o output.sorted.rmdup.realigned.bi.vcf
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fasta -V output.sorted.rmdup.realigned.bi.vcf -o output.sorted.rmdup.realigned.bi.filtered.vcf --filterExpression "DP < 8 || DP > 31 || FS > 60.0 || MQ < 20.0 || MQRankSum < -12.5 || QD < 2.0 || ReadPosRankSum < -8.0" --filterName "LOW_QUALITY"
grep -v "LOW_QUALITY" output.sorted.rmdup.realigned.bi.filtered.vcf > output.final.vcf
###BAM files were composed of autosomal scaffolds of all lengths (not just those >2 Mb)
###BAM files were not subsampled to 5x sequence depth

#(Step 23) Process data and run PSMC
bgzip -c output.final.vcf > output.final.vcf.gz
tabix -p vcf output.final.vcf.gz
zcat output.final.vcf.gz | vcf-subset -e -c output | grep -e "^#" -e "0/1" | bgzip > output.vcf.gz
tabix -p vcf output.vcf.gz
cat reference.fasta | seqtk seq -M reference.repeat.merged.bed -A - | vcf-consensus -i output.vcf.gz > output.fa
fq2psmcfa output.fa > output.psmcfa
splitfa output.psmcfa > output.split.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o output.psmc output.psmcfa
for i in {1..100}; do /data/ochoa.43/psmc-master/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o output.$i.psmc output.split.psmcfa; done
cat output.psmc output.{1..100}.psmc > output.boostrapped.psmc
psmc_plot.pl -P top -X3000000 -u7.2e-09 -p -g3 -T output -x10000 output.plot output.boostrapped.psmc
###Used VCFtools v.0.1.16 (http://vcftools.sourceforge.net/) and PSMC v.0.6.4 (https://github.com/lh3/psmc)


------------------------------------Mutation Load-------------------------------------

#(Step 24) Run BUSCO
run_BUSCO.py -i reference.fasta -o busco_results -l ~/tetrapoda_odb10/ -m genome -c 8
###Used BUSCO v.3.1.0 (https://busco.ezlab.org/)
###FASTA file was composed of autosomal scaffolds of all lengths (not just those >2 Mb)
###We then converted and merged the resulting GFF file into the BED format as outlined in Step 20
###For practical purposes, this BED file will be called BUSCOs.bed

#(Step 25) Index the BED file and obtain genotypes from each population using ANGSD
angsd sites index BUSCOs.bed
angsd -b CEBO.filelist -minMapQ 20 -minQ 20 -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -doGlf 3 -doGeno 4 -doPost 1 -postCutoff 0.95 -sites BUSCOs.bed -out CEBO
###BAM files were composed of autosomal scaffolds of all lengths (not just those >2 Mb)
###BAM files were subsampled to 5x sequence depth
###The resulting genotype file (CEBO.geno.gz) contains genotype calls across individuals from CEBO
###We then created a BED file (called BUSCOs.sites.bed) with a list of unique SNPs across all populations

#(Step 26) Call genotypes in each species (S. catenatus, S. tergeminus, and S. miliarius) separately
###As an example for S. catenatus:
GenomeAnalysisTK.jar -T GenotypeGVCFs -R reference.fasta -V SCAT.filelist -L BUSCOs.sites.bed -o SCAT.vcf --includeNonVariantSites
###BAM files were composed of autosomal scaffolds of all lengths (not just those >2 Mb)
###BAM files were not subsampled to 5x sequence depth
###We used WGS from S. miliarius as an outgroup to S. catenatus and S. tergeminus (see Section S2 of the Supplemental Information)
###This step is intended to identify monomorphic, fixed sites in each population and species since these sites were not called in ANGSD
###This step is also intended to identify 'ancestral' homozygous genotypes in S. miliarius
###We used similar parameters as in Step 21 for filtering the data (also see Section S2 of the Supplemental Information)
###We used information from Steps 24 and 25 to identify ancestral and derived alleles in S. catenatus and S. tergeminus
###From here, we used a similar 'vcf-consensus' approach as in Step 22 to create BUSCO consensus sequences based on fixed+ancestral and fixed+derived alleles for S. catenatus and S. tergeminus

#(Step 27) Annotate derived alleles into nonsynonymous and synonymous
###For each species-specific fixed+ancestral and fixed+derived FASTA files we ran the following command line:
bedtools getfasta -fo out -name -tab -s -fi input.fasta -bed BUSCOs.bed
###We then concatenated and translated the coding sequences (CDS) from each BUSCO and compared fixed+ancestral and fixed+derived files from each species to annotate derived alleles
###Please see MEC-21-0485/06_Mutation_load/MutationLoad.txt in Dryad for annotations

#(Step 28) Run PROVEAN on the nonsynonymous alleles from each species
provean.sh -q BUSCOlocus.fasta -v BUSCOlocus.var --num_threads 16 --save_supporting_set BUSCOlocus.sss
###Used PROVEAN v.1.1 (http://provean.jcvi.org/index.php)
###The FASTA file corresponds to a species-specific fixed+ancestral amino acid sequence for each BUSCO locus
###The VAR file contains a list (and positions) of species-specific amino acid substitutions for each BUSCO locus