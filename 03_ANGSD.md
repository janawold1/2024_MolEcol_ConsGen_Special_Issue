# Population genomics with ANGSD

## Site Frequency Spectrum
```
for POP in AU TI
    do
    echo "CALCULATING INITIAL SFS FILES FOR POPULATION $POP..."
    angsd -P 8 -b ${POP}.list -ref $ref -anc $ref \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 15 -minInd 11 -skipTriallelic 1 -setMinDepth 44 -setMaxDepth 770 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -doSaf 1 -doCounts 1 -out ${POP}
    echo "FINISHED INITIAL SFS FILES FOR POPULATION $POP..."
done

realSFS AU.saf.idx TI.saf.idx > AU_TI.sfs
```