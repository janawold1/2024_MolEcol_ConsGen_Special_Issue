# Population genomics with ANGSD
ANGSDvX.X and ngsTools were used to estimate summary statistics for the *de novo* tara iti assembly and scaffolded as per the common tern reference (see [01_read_QC_and_alignment.md](github.com/janawold1/2023_EVOLAPP_Special_Issue)).  
## Summary Statistics
Nucleotide diversity, inbreeding coefficients and Fst were estimated as per below.
```
for POP in AU TI
    do
    angsd -P 8 -b ${POP}.list -ref $ref -anc $ref \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \ 
        -trim 0 -C 50 -minMapQ 20 -minQ 20 -min Ind 11 \
        -skipTriallelic 1 -setMinDepth 66 -setMaxDepth 770 /
        -GL 1 -doMajorMinor 1 -doMaf -doCounts 1 -doSaf 1 XXXX
done
```
## Site Frequency Spectrum
The intermediate site frequency spectrum estimated in the example above were used to generate SFS files with realSFS.  
```
realSFS AU.saf.idx > AU.sfs
realSFS TI.saf.idx > TI.sfs
realSFS AU.saf.idx TI.saf.idx > AU_TI.sfs
```