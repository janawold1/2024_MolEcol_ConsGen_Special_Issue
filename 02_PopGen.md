# Population genomics with ANGSD
[ANGSDv0.940-dirty](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) (htslib:1.18) and ngsTools were used to estimate summary statistics for the *de novo* tara iti assembly and scaffolded as per the common tern reference (see [01_read_QC_and_alignment.md](github.com/janawold1/2023_EVOLAPP_Special_Issue/blob/main/01_read_QC_and_alignment.md) for details). The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. This was repeated alignments to all three reference geomes (Common tern, Tara iti scaffolded using common tern, and the unscaffolded tara iti reference). Sex chromosomes were excluded for alignments to the scaffolded assemblies (i.e., common tern, tara iti scaffolded using common tern). All analyses were performed under the SAMtools genotype likelihood model `-GL 1` unless otherwise denoted.  
```
angsd -P 8 -b GLOBAL.list -ref $ref -out qc/GLOBAL.qc \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -trim 0 -C 50 -baq 1 -minMapQ 15 -doQsdist 1 -doDepth 1 \
    -doCounts 1 -maxDepth 1000
```
Courtesy of scripts provided by [@mfumagalli](github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.

|      Population     |minimum mapQ|minimum Q|minimum depth|maximum depth|minimum individuals|
| ------------------- | ---------- | ------- | ----------- | ----------- | ----------------- |
|Australian Fairy Tern|     20     |    20   |     114     |     350     |        19         |
|      Tara iti       |     20     |    20   |     114     |     350     |        19         |
|       Global        |     20     |    20   |     300     |     630     |        38         |

A missingness threshold of 0% or 10% were used when appropriate for all subsequent analyses. Below is code run for analyses performed on alignments to the tara iti reference scaffolded using the common tern assembly. All analyses were performed in a similar manner, with the variables adjusted as per the table above.  
## Removing coding Regions
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.31.0 was used to find the complement of the UCSC annotation and exclude coding regions of the genome. An additional 1kb of sequence on either side of annotations were included to reduce linkage.  

First, the annotations for autosomal chromosomes were extracted and the window size for these regions increased by 1kb on either side for each annotation file.  
```
zcat reference/GCA_009819605.1/genes/GCA_009819605.1_bSteHir1.pri.xenoRefGene.gtf.gz | \
    sort | \
    grep -v WNM | \
    grep -v CM020500 | \
    grep -v CM020463 | \
    grep -v CM020462 | \
    awk '{print $1"_RagTag\t"$4-1000"\t"$5+1000"\t"$2}'> reference/GCA_009819605.1/genes/all_autosomal_annotations.bed

 zcat reference/GCA_009819605.1/genes/GCA_009819605.1_bSteHir1.pri.augustus.gtf.gz | \
    sort | \
    grep -v WNM | \
    grep -v CM020500 | \
    grep -v CM020463 | \
    grep -v CM020462 | \
    awk '{print $1"_RagTag\t"$4-1000"\t"$5+1000"\t"$2}' >> reference/GCA_009819605.1/genes/all_autosomal_annotations.bed

sort -k1,1 -k2,2n reference/GCA_009819605.1/genes/all_autosomal_annotations.bed > reference/GCA_009819605.1/genes/all_autosomal_annotations_sorted.bed
```
These files were then merged, leaving roughly 22,173 annotated genes.  
```
bedtools merge \
    -i reference/GCA_009819605.1/genes/all_autosomal_annotations_sorted.bed \
    -c 1,2,3 -o distinct > reference/GCA_009819605.1/genes/all_autosomal_annotations_merged.bed
```
The complement regions were extracted with BEDtools and used to filter regions for ANGSD.  
```
awk -v OFS='\t' {'print $1,$2'} reference/TI_scaffolded_as_CT.fasta.gz.fai > reference/TI_scaffolded_as_CT_bedtools_genome.txt

bedtools complement -i reference/GCA_009819605.1/genes/all_autosomal_annotations_sorted.bed \
    -g reference/TI_scaffolded_as_CT_bedtools_genome.txt | grep -v sf | grep -v WNM | \
    grep -v CM020500 | grep -v CM020463 | grep -v CM020462 > angsd_scaffolded/TI_scaffolded_neutral_regions.bed
```
No regions within 1kb of the start/end of the chromosome were retained to account for poor assembly at repetitive regions of the genome. This left 22,167 regions covering 615,256,912 bp for analyses. This corresponds to roughly 52% of the genome.  

This file was then indexed for use with the `angsd -sites` option.  
```
angsd sites index TI_scaffolded_neutral_regions.bed
```

We also extracted the autosomoal chromosomes from the common tern assembly, and renamed them with the same names in the tara iti reference assembly. This was so we could polarize the site frequency spectrum with common tern as the ancestral state. 

# Test for Sex Specific differences?

## Remember ANGSD method of reconstructing outgroup when TI ONT assembly complete
```
faSomeRecords common_tern.fa common_tern_autosomes.bed common_tern_autosomes.fasta

awk '{print $1"\t"$1"_RagTag"} common_tern_autosomes.bed > renames.txt

while read -r line
    do
    old=$(echo $line | awk '{print $1}')
    new=$(echo $line | awk '{print $2}')
    echo "RENAMING $old TO $new"
    sed -i "s/$old/$new/g" common_tern_autosomes.fasta
done < renames.txt
```
Once the autosomes were extracted and renamed, the common tern reference had to be subsetted to be usable for SFS analyses with ANGSD and GADMA. To ensure the appropriate regions of the genome were included for polarising the SFS, the common tern assembly was aligned to the scaffolded tara iti reference with [minimap2](https://github.com/lh3/minimap2) v2.26-r1175. This is by no means a perfect solution, but helped with maximising the likelihood that the appropriate region from the common tern assembly was used.  
```
minimap2 -x asm5 TI_scaffolded_as_CT.fasta.gz common_tern_autosomes.fasta.gz -a -o tern2tern.sam
samtools view -F tern2tern.sam | samtools fasta > common_tern_subset.fasta
```
```
bgzip common_tern_subset.fasta
samtools faidx common_tern_subset.fasta.gz
```
## Inbreeding Estimates
### F statistics
To estimate relative levels of inbreeding and account for the high liklihood that tara iti do not conform to Hardy-Weinberg Equilibrium (HWE), we first generated population specific genotype likelihoods for input into [ngsF](https://github.com/mfumagalli/ngsTools) v1.2.0-STD. The outputs of these analyses can also be used as a prior for populations that are not in hardy-weinburg equilibrium (HWE).  
```
angsd -P 8 -b AU.list -ref $ref -out inbreeding/AU -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 350 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 8 -b TI.list -ref $ref -out inbreeding/TI -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3

angsd -P 8 -b GLOBAL.list -ref $ref -out inbreeding/GLOBAL -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGlf 3
```
This left 980,578 sites in the `AU` dataset, 297,787 sites in the `TI` dataset, and 1,665,527 sites in the `GLOBAL` dataset upon the completion of ANGSD. Next, `ngsF` was used to estimate inbreeding. First, an initial search was performed. The below was perfomed for each of the `AU`, `TI`, and `GLOBAL` data sets upon the completion of ANGSD.  
```
zcat inbreeding/AU.glf.gz > inbreeding/AU.glf

NSITES=$(zcat inbreeding/AU.mafs.gz | tail -n+2 | wc -l)

ngsF --n_ind 19 --n_sites $NSITES --glf inbreeding/AU.glf --out inbreeding/AU_approx_indF \
    --approx_EM --init_values u --n_threads 8

ngsF --n_ind 19 --n_sites $NSITES --glf inbreeding/AU.glf --out inbreeding/AU.indF \
    --init_values inbreeding/AU_approx_indF.pars --n_threads 8

cat inbreeding/AU.indF
```
We also estimated per-individual inbreeding with [ngsF-HMM](https://github.com/mfumagalli/ngsTools) v1.1.0.
```

```
### Relatedness
We estimated relatedness among individuals with [ngsRelate](https://github.com/ANGSD/NgsRelate).  
```
for POP in AU TI
    do
    zcat inbreeding/${POP}.mafs.gz | cut -f 6 | sed 1d > inbreeding/${POP}_freqs
    NSITES=$(zcat inbreeding ${POP}.mafs.gz | tail -n +2 | wc -l)
    ngsRelate -g inbreeding/${POP}.glf -n 19 -f inbreeding/${POP}_freqs -O inbreeding/${POP}_relatedness
done
```

### Runs of Homozygosity
Two methods were used to estimate runs of homozygosity, [ROHAN](https://github.com/grenaud/rohan) and an approach using ANGSD.  

For ROHAN, all PCR-duplicates were removed as recommended. 
```
samtools markdup -@16 -r --write-index *_nodup_autosomes.cram *_rohan.bam
```
And the programme was run on for each individual.
```
for bam in *_rohan.bam
    do
    base=$(basename $bam _rohan.bam)
    echo "RUNNING ROHAN FOR ${base}"
    rohan --size 50000 --rohmu 4.6e-9 -t 16 --tstv 2.36 -o output/${base} $ref $bam
done
```
### ROH-Size Clustering
Run in R
Create a one-column file (e.g., scat_rohs.txt) with ROH sizes across all samples
These ROH sizes were obtained from above after excluding microchromosomal scaffolds
```
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
```

### Global Heterozygosity
Here, we implemented a global (genome-wide heterozygosity) method from ANGSD. Essentially, this estimate is a proportion of heterozygous genotypes / genome size (excluding regions of the genome with low confidence). Unlike other runs of ANGSD, individual CRAMs are used to estimate hetereozygosity, which is simply second value in the SFS/AFS.  
```
for cram in ${dir}*_nodup_autosomes.cram
    do
    base=$(basename $cram _nodup_autosomes.cram)
    angsd -i $cram -anc $ref -out heterozygosity/${base} -dosaf 1 -GL 1 -doCounts 1
    realSFS -fold 1 heterozygosity/${base}.saf.idx > heterozygosity/${base}_est.ml
done
```

## Population Structure
### PCA
To visualise population structure using a PCA, we first ran ANGSD using the `-doGeno 32` and `-doPost 1` options.  
```
angsd -P 8 -b GLOBAL.bamlist -ref $ref -out structure_PCA/GLOBAL_noMiss -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1
```
Then, [ngsCovar](https://github.com/mfumagalli/ngsPopGen) was run to estimate covariance.  
```
gunzip structure_PCA/GLOBAL_noMiss.geno.gz

NSITES=$(zcat structure_PCA/GLOBAL_noMiss.mafs.gz | tail -n+2 | wc -l)
ngsCovar -probfile structure_PCA/GLOBAL_noMiss.geno \
    -outfile structure_PCA/GLOBAL_noMiss.covar -nind 38 -nsites $NSITES -call 0 -norm 0
```
Before plotting with the Rscripts supplied by ngsTools.  
```
Rscript -e 'write.table(cbind(seq(1,38),rep(1,38),c(rep("AU",19),rep("TI",19))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="structure_PCA/GLOBAL_noMiss.clst", quote=F)'
Rscript plotPCA.R -i structure_PCA/GLOBAL_noMiss.covar -c 1-2 -a structure_PCA/GLOBAL_noMiss.clst -o structure_PCA/GLOBAL_noMiss.pca.pdf
```
### MDS
To construct a MDS for fairy tern populations, we first ran ANGSD as below. Notably, the only change is the `-doGeno 8` flag.  
```
angsd -P 8 -b GLOBAL.list -ref $ref -out structure_MDS/GLOBAL -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 38 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1
```
This command not only generates the required `*.mafs.gz`, but also the `*.geno.gz`, which contains posterior probabilities of all possible genotypes required for estimating genetic distance with [ngsDistv1.0.10](github.com/mfumagalli/ngsTools). In addition, a `pops.label` file denoting the population of origin (one entry on a new line for each samples) is necessary for estimating genetic distance.  
```
NSITES=$(zcat GLOBAL.mafs.gz | tail -n +2 | wc -l)
echo $NSITES

ngsDist -verbose 1 -geno distance/GLOBAL.geno.gz -probs \
    -n_ind 38 -n_sites $NSITES -labels pops.label -o distance/GLOBAL.dist
```
Extract and construct MDS.  
```
tail -n +3 distance/GLOBAL.dist | Rscript --vanilla --slave getMDS.R \
    --no_header --data_symm -n 4 -m "mds" -o distance/GLOBAL.mds

head distance/GLOBAL.mds
```
And finally plot the MDS.  
```
Rscript plotMDS.R -i distance/GLOBAL.mds -c 1-2 -a structure/GLOBAL_noMiss.clst -o distance/GLOBAL_mds.pdf
```
### Population Structure with Inbreeding
Initial inbreeding estimates for tara iti indicate that the population is likely out of hardy-weinburg equilibrium (HWE). To account for this, relative inbreeding levels were incorporated into assessments of population structure.  

### Geographic Population Structure (GPS)
Analyses of population structure using PCAs have been shown to exhibit bias ([Elhaik et al 2022](https://doi.org/10.1038/ncomms4513)). Admixture based appraoches, like those implemented in [GPS](https://github.com/arash-darzian/Geographic_Population_Structure_GPS/tree/main) have been proposed. This program leverages XXX may not want to use as it has a strong aDNA focus XXX

## Summary Statistics
### Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS. Here, we are using the common tern as the ancestral state.  
```
anc=reference/common_tern_autosomes.fasta.gz
region=TI_scaffolded_neutral_regions.bed
for POP in AU TI
    do
    angsd -P 8 -b ${POP}.list -ref $ref -anc $ref -out sfs/${POP}_fold -sites ${region} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
        -GL 1 -doSaf 1
done

angsd -P 8 -b GLOBAL.list -ref $ref -anc $anc -out sfs/GLOBAL -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 300 -setMaxDepth 630 -doCounts 1 \
    -GL 1 -doSaf 1
```

```
realSFS -fold 1 -bootstrap 100 -tole 1e-6 sfs/AU_fold.saf.idx > sfs/AU_fold_realSFS_folded_100boots.sfs

realSFS -fold 1 -bootstrap 100 -tole 1e-6 sfs/TI_fold.saf.idx > sfs/TI_fold_realSFS_folded_100boots.sfs

realSFS -fold 1 -bootstrap 100 -tole 1e-6 sfs/AU_fold.saf.idx sfs/TI_fold.saf.idx > sfs/total_fold_realSFS_folded_100boots.sfs
```
## Theta Statistics
Using the ANGSD `-doSaf 1` and realSFS outputs, theta neutrality statistics were estimated using realSFS and associated thetaStat programmes.  
```
for POP in AU TI total
    do
    tail -n1 sfs/${POP}_fold_realSFS_folded_100boots.sfs > diversity/${POP}_fold.sfs
    realSFS saf2theta sfs/${POP}_fold.saf.idx -outname diversity/${POP}_fold -sfs diversity/${POP}_fold.sfs -fold 1
    thetaStat do_stat diversity/${POP}_fold.thetas.idx
    thetaStat do_stat diversity/${POP}_fold.thetas.idx -win 50000 -step 10000 -outnames diversity/${POP}_fold_thetasWindow
done
```
Only the Watterson's Theta and Tajima's D can be estimated with a folded SFS. The mean value of Tajima's D was estimated by the below.  

Genetic differentiation.  
```
realSFS fst index sfs/AU.saf.idx sfs/TI.saf.idx -sfs diversity/total_fold.sfs -fstout distance/total_fst

realSFS fst stats distance/total_fst.idx \
    -tole 1e-6 -win 50000 -step 10000 -whichFst 1 > distance/total_fst_results.tsv
```
This yielded the result:  

```
-> Assuming idxname:distance/test.fst.idx
        -> Assuming .fst.gz file: distance/test.fst.gz
        -> args: tole:0.000001 nthreads:4 maxiter:100 nsites(block):0 start:(null) chr:(null) start:-1 stop:-1 fstout:(null) oldout:0 seed:-1 bootstrap:0 resample_chr:0 whichFst:1 fold:1 ref:(null) anc:reference/TI_scaffolded_as_CT.fasta.gz
        -> FST.Unweight[nObs:387393642]:0.002465 Fst.Weight:0.857685
```

We also ran:
```
realSFS fst stats2 distance/total_fst.idx \
    -tole 1e-6 -fold 1 -anc $ref -win 50000 -step 10000 -whichFST 1 > distance/total_fst2_results.tsv
```

```
angsd -P 8 -b AU.list -ref ${ref} -anc ${ref} -out diversity/AU_folded -sites ${region} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
    -minMapQ 20 -minQ 20 -minInd 19 -setMinDepth 114 -setMaxDepth 300 -doCounts 1 \
    -GL 1 -doSaf 1 -doThetas 1 -pest sfs/AU_fold.sfs
```
## Generating VCF
For other programs and to estimate a Ts/Tv rate for Rohan, a BCF was produced using ANGSD.  
```
angsd -P 8 -b GLOBAL.list -ref $ref -out angsd/genotypes/global_GATK_genotypes -uniqueOnly 1 \
    -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 -minMapQ 20 -minQ 20 -minInd 38 \
    -setMinDepth 300 -setMaxDepth 630 -doCounts 1 -skipTriallelic 1 -doBcf 1 -GL 2 \
    -doPost 1 -doMaf 1 -doGeno 10 --ignore-RG 0 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3

angsd -P 8 -b GLOBAL.list -ref $ref -out angsd/genotypes/global_samtools_genotypes -uniqueOnly 1 \
    -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 -minMapQ 20 -minQ 20 -minInd 38 \
    -setMinDepth 300 -setMaxDepth 630 -doCounts 1 -skipTriallelic 1 -doBcf 1 -GL 1 \
    -doPost 1 -doMaf 1 -doGeno 10 --ignore-RG 0 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3
```
