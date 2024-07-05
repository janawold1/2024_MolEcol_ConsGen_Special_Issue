## Historical N<sub>e</sub> Inference with PSMC
We created a consensus diploid sequence for three tara iti, Australian fairy tern and kakī samples. The conversion from FASTQ format to PSMCfa must use a `gzipped` file as input.  
```
for BAM in variants/fairy_bam/nodup/{AU01,AU27,AU28,SND06,SP01,TI37}_nodup_autosomes.bam
    do
    BASE=$(basename $BAM _nodup_autosomes.bam)
    printf "BEGAN GENERATING CONSENSUS SEQUENCE FOR $BASE AT "
    date
    bcftools mpileup --threads 26 -Ou -Q 30 -q 20 -f $REF $BAM | \
        bcftools call --threads 26 -c | \
        vcfutils.pl vcf2fq -d 8 -D 125 -Q 30 > ${PSMC}fastq/${BASE}.fq
    gzip ${PSMC}fastq/${BASE}_psmc.fq
    fq2psmcfa -q20 ${PSMC}fastq/${BASE}_psmc.fq.gz > ${PSMC}psmcfa/${BASE}.psmcfa
    printf "FINISHED GENERATING CONSENSUS SEQUENCE FOR $BASE AT "
    date
done
```
Once the file was prepared, we ran PSMC under the the following conditions... We assumed a generation time of 3 years for tara iti as the age of first breeding ranges between 2-4 years of age. For kakī, this was increased to 6 years. The avian mutation rate is estimated to range from 1.23 x 10<sup>-9</sup> - 2.21 x 10<sup>-9</sup>. Assuming the lower range of this estimate, a rate of 3.69 x 10<sup>-9</sup> and 7.38 x 10<sup>-9</sup> were used for fairy terns and kakī respectively.  
```
for FA in psmc/psmcfa/*.psmcfa
    do
    SAMP=$(basename $FA .psmcfa)
    splitfa $FA > psmc/psmcfa/${SAMP}_split.psmcfa
    printf "STARTED RUNNING PSMC FOR $SAMP AT "
    date
    psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o psmc/out/${SAMP}_diploid.psmc $FA
    seq 1 100 | parallel -j 10 psmc -N30 -t5 -r5 -b -p "4+30*2+4+6+10" -o psmc/out/${SAMP}_round-{}.psmc psmc/psmcfa/${SAMP}_split.psmcfa
    cat psmc/out/${SAMP}_*.psmc > psmc/out/${SAMP}_combined.psmc
    printf "FINISHED RUNNING PSMC FOR $SAMP AT "
    date
done

psmc_plot.pl -u 3.69e-09 -g 3 psmc/AU psmc/out/AU*_combined.psmc
psmc_plot.pl -u 3.69e-09 -g 3 psmc/TI psmc/out/TI*_combined.psmc
psmc_plot.pl -u 7.38e-09 -g 6 psmc/KI psmc/out/KI*_combined.psmc
```
## Contemporary N<sub>e</sub> Inference with StairwayPlot2
We extracted the variant sites from the SFS for Australian fairy tern, tara iti and kakī. These were then put into a blueprint file for StairwayPlot2. The nRand variable was run using the default settings of (nseq-2)(0.25), (nseq-2)(0.5), (nseq-2)(0.75), and (nseq-2). 

Once the blueprint file was prepared, the bash file was created with:
```
for POP in AU TI KI
    do
    java -cp stairway_plot_es Stairbuilder ${POP}_two-epoch.blueprint
    bash ${POP}_two-epoch.blueprint.sh
done
```
## Estimating fairy tern time of divergence with hPSMC


## Population Demography and Connectivity Inference with GADMA
[GADMA](https://gadma.readthedocs.io/en/latest/user_manual/input_data/snp_data_format.html) leverages the joint SFS to infer the demographic history of multiple populations. It can implement [dadi](https://bitbucket.org/gutenkunstlab/dadi/), [moments](https://github.com/MomentsLD/moments), [momi2](https://github.com/popgenmethods/momi2/), and [momentsLD](https://github.com/MomentsLD/moments).  

GADMA can take multiple input formats. Here we estimated the joint SFS in ANGSD and realSFS as above, with the exception that the `dadi` flag was on when running realSFS. This output was then converted with a perl script [realSFS2dadi.pl](https://github.com/z0on/2bRAD_denovo/blob/master/realsfs2dadi.pl).
```
realSFS dadi -ref $REF -anc $ANC AU.saf.idx TI.saf.ids -sfs AU.sfs -sfs TI.sfs > GLOBAL.dadi
realSFS2dadi.pl GLOBAL.dadi 19 15 > GLOBAL_GADMA_SNPformat.txt
```
GADMA is relatively straightforward and easy to run. A parameter file defining specific settings can be used as input. One important parameter is `sequence length`, which denotes the number of sites used to build the data (SFS in our case). Fortunately, the realSFS programme includes this (`nSites`) as part of its progress output. For the neutral dataset this was 512,061,918 sites while it was 976,823,680 sites for the whole-genome data set.  

Because the SNPs in either SFS were not filtered for linkage, we updated the paramfile to reflect this and provided a directory for bootstrapping. This has important implications for model selection methods, see [here](https://gadma.readthedocs.io/en/latest/user_manual/input_data/input_data.html#extra-information-about-data) for more information.  

As with PSMC and StairwayPlot2 above, we ran GADMA with a conservative mutation rate of 1.23e-9 and a generation time of 3 years. Finally, the Selection and Dominance options were set to true prior to running with `gadma -p GADMA_neutral.params -o GADMA_neutral_moments/`.  