# Tara iti repeat annotation, 18-25 September 2024

Running the repeat modelling and annotation for tara iti in response to feedback on the ms. Whole process takes around 46 hrs.

Input genome assembly file is `Katie_5kb_ragtag.fa`. 

Scripts to run a basic repeat modelling and annotation pipeline (based on https://darencard.net/blog/2022-07-09-genome-repeat-annotation/):

```
01-build-repmodDB.sh
02-repmod.sl
03-repclass.sh
03b-repclass.sh
04a-repmask.sl
04b-repmask.sl
04c-repmask.sl
05-consolidate.sh
06-summarise.sh
07-makeGFF.sh
08-maskfasta.sh
```

Software:

```
RepeatModeler v2.0.3
SeqKit v2.4.0 - only for minor data wrangling
RepeatMasker v4.1.0
RMBlast v2.10.0
repclassifier (https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier)
seqtk v1.4 - only for minor data wrangling
rmOutToGFF3custom (https://github.com/darencard/GenomeAnnotation/blob/master/rmOutToGFF3custom)
BEDTools v2.30.0
```

1. `01-build-repmodDB.sh` - build an assembly-specific database for RepeatModeler
2. `02-repmod.sl` - run RepeatModeler to identify de novo repeats
3. `03-repclass.sh` - split identified repeats into those that were successfully classified (known) and those that weren't (unknown), then use `repclassifier` that employs RepeatMasker to attempt to classify those unknown repeats
4. `03b-repclass.sh` - continue classifying unknown repeats via an iterative process
5. `04x-repmask.sl` - serially annotate the assembly with RepeatMasker, by annotating: a) simple repeats, b) well-curated repeats for birds, c) known species-specific repeats and then unknown species-specific repeats
6. `05-consolidate.sh` - consolidate across the rounds of 04x outputs to produce sets for simple repeats, complex repeats, and all repeats, then produce a summary of repeat composition with RepeatModeler
7. `06-summarise.sh` - produce an alternative (more informative) summary of repeat composition
8. `07-makeGFF.sh` - use `rmOutToGFF3custom` to convert outputs to GFF3 format
9. `08-maskfasta.sh` - use BEDTools to mask the genome assembly file in two ways: 1) simple repeats are soft-masked (replaced with lowercase), 2) complex repeats are hard masked (replaced with Ns)

## Notes:

During the iterative process of repeat classification, I retrieved 334 repeats in the  `tern-families.fa` file as output of `02-repmod.sl`. Of these, there are 196 known and 138 unclassified repeats. 
- Round 1 repclassifier: known = 213, unknown = 121
- Round 2 repclassifier: known = 220, unknown = 114
- Round 3 repclassifier: known = 221, unknown = 113
- Round 4 repclassifier: known = 224, unknown = 110
- Round 5 repclassifier: known = 225, unknown = 109
- Round 6 repclassifier: known = 225, unknown = 109.

This process plateaued at round 5, so I removed round 6 outputs and used round 5 outputs as input to `04c-repmask.sl`. 


Assembly fasta headers were too long for RepeatMasker (it can't handle headers over ~40 characters in length). So before running `04a-repmask.sl`, I have produced a version of the genome (`Katie_5kb_ragtag-renamed.fa`) with the headers renamed to format: `>taraIti1_1 >taraIti1_2` etc., using:

```bash
ml purge; ml SeqKit/2.4.0
seqkit replace -p '.+' -r 'taraIti1_{nr}' /nesi/nobackup/ga03186/taraiti-repeats/Katie_5kb_ragtag.fa > /nesi/nobackup/ga03186/taraiti-repeats/Katie_5kb_ragtag-renamed.fa
```

## Results 

Summarised output of `06-summarise.sh`:

```log
Repeat_group    Repeat_subgroup   TotalBP     %genome
ARTEFACT        NA                80          6.62162e-06
DNA             Crypton           85794       0.00710119
DNA             Crypton-A         35775       0.00296111
DNA             Kolobok           190247      0.0157468
DNA             Merlin            43684       0.00361574
DNA             NA                1116559     0.0924179
DNA             PIF-Harbinger     1578934     0.130689
DNA             TcMar             50875       0.00421094
DNA             TcMar-Pogo        13921       0.00115225
DNA             TcMar-Tigger      109377      0.00905317
DNA             hAT               175832      0.0145537
DNA             hAT-Ac            174327      0.0144291
DNA             hAT-Blackjack     106184      0.00878888
DNA             hAT-Charlie       104158      0.00862119
DNA             hAT-Tag1          54029       0.004472
DNA             hAT-Tip100        3943        0.000326363
DNA             hAT-hAT19         2426        0.000200801
LINE            CR1               72911494    6.0349
LINE            L1                22911       0.00189635
LINE            L2                457263      0.0378478
LINE            Penelope          32872       0.00272082
LTR             DIRS              7431        0.000615066
LTR             ERV               510703      0.042271
LTR             ERV1              1495054     0.123746
LTR             ERVK              2346723     0.194239
LTR             ERVL              12014231    0.994421
LTR             Gypsy             25297       0.00209384
LTR             NA                36822       0.00304777
Low_complexity  NA                2894190     0.239553
RC              Helitron          14087       0.00116598
SINE            5S-Deu-L2         346792      0.0287041
SINE            MIR               498840      0.0412891
SINE            tRNA              347298      0.028746
SINE            tRNA-CR1          12067       0.000998789
SINE            tRNA-Deu          25040       0.00207257
Satellite       NA                369839      0.0306117
Simple_repeat   NA                12099704    1.0015
Unknown         NA                8444372     0.698943
Unspecified     NA                7188128     0.594963
rRNA            NA                87213       0.00721864
scRNA           NA                355         2.93834e-05
snRNA           NA                17566       0.00145394
tRNA            NA                23974       0.00198433
```

And `05-consolidate.sh` RepeatModeler summary:

```log
==================================================
file name: Katie_5kb_ragtag-renamed.full_mask
sequences:           137
total length: 1208163107 bp  (1107585427 bp excl N/X-runs)
GC level:        Unknown %
bases masked:  122591736 bp ( 10.15 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements       209593     88989246 bp    7.37 %
   SINEs:             9612      1222839 bp    0.10 %
   Penelope            153        32817 bp    0.00 %
   LINEs:           171510     72325410 bp    5.99 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex     171263     72269805 bp    5.98 %
     R1/LOA/Jockey       0            0 bp    0.00 %
     R2/R4/NeSL          0            0 bp    0.00 %
     RTE/Bov-B           0            0 bp    0.00 %
     L1/CIN4            94        22788 bp    0.00 %
   LTR elements:     28471     15440997 bp    1.28 %
     BEL/Pao             0            0 bp    0.00 %
     Ty1/Copia           0            0 bp    0.00 %
     Gypsy/DIRS1       143        32653 bp    0.00 %
       Retroviral    28131     15371619 bp    1.27 %

DNA transposons      25729      3770430 bp    0.31 %
   hobo-Activator     3650       614764 bp    0.05 %
   Tc1-IS630-Pogo      960       172305 bp    0.01 %
   En-Spm                0            0 bp    0.00 %
   MuDR-IS905            0            0 bp    0.00 %
   PiggyBac              0            0 bp    0.00 %
   Tourist/Harbinger 11044      1518223 bp    0.13 %
   Other (Mirage,        0            0 bp    0.00 %
    P-element, Transib)

Rolling-circles         87        14001 bp    0.00 %

Unclassified:        52008     15011649 bp    1.24 %

Total interspersed repeats:   107771325 bp    8.92 %

Small RNA:            3322       470885 bp    0.04 %

Satellites:            763       333647 bp    0.03 %
Simple repeats:     284008     11559216 bp    0.96 %
Low complexity:      56922      2789307 bp    0.23 %
==================================================
```

These represent two slightly different ways of describing the repeat landscape identified. 

The primary outputs are masked assembly fasta files `Katie_5kb_ragtag-renamed.simple_mask.soft.fasta` and `Katie_5kb_ragtag-renamed.simple_mask.soft.complex_mask.hard.fasta`.

## Suggested text for the ms/supplementary

We followed the method of Card 2022 (https://darencard.net/blog/2022-07-09-genome-repeat-annotation/) to identify, classify, annotate, and mask repeats in the tara iti genome assembly, implementing RepeatModeler v2.0.3, RepeatMasker v4.1.0, RMBlast v2.10.0, BEDTools v2.30.0, and custom scripts repclassifier (https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier) and rmOutToGFF3custom (https://github.com/darencard/GenomeAnnotation/blob/master/rmOutToGFF3custom). Scripts for the full pipeline are available at ([GitHub repo]). In the final pipeline, we conducted five rounds of iterative classification of unknown repeats prior to full repeat annotation. 10.15% of the tara iti genome was masked as containing simple and complex repeats.

[Add `05-consolidate.sh` RepeatModeler summary as supplementary?]

