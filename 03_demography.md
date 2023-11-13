# Demographic Inference with GADMA
[GADMAvX.X](https://github.com/ctlab/GADMA/tree/master) was used for demographic inference. The *dadi* and *moments* engines were run as per the parameter file below.  

Data was first converted from the SAF output format of ANGSD.  
```
realSFS dadi -fold 1 sfs/AU_fold.saf.idx sfs/TI_fold.saf.idx \
    -sfs diversity/AU_fold.sfs -sfs diversity/TI_fold.sfs -ref ${ref} -anc ${ref} > gadma_folded/AU_TI_folded_dadi.tsv
./realsfs2dadi.pl AU_TI_folded_dadi.tsv 19 19 > AU_TI_folded_gadma.tsv
```
Below is an example of one parameter file used as input into GADMA. The main parameters varied included the `engine`, where both the `dadi` and `moments` engines were trialled. More information on each of these variables can be found [here]().  
```
Output directory: /media/jana/BigData/tara_iti_publication/gadma_folded
Input data: /media/jana/BigData/tara_iti_publication/AU_TI_folded_gadma.tsv
Population labels: [AU, TI]
Projections: [19, 19]
Outgroup: False
Sequence length: 1068161
Linked SNP's: True
Directory with bootstrap: /media/jana/BigData/tara_iti_publication/gadma_boot
Engine: moments
Pts: [19, 29, 39]
Theta0: Null
Mutation rate: 4.6e-9
Recombination rate: Null
Time for generation: 3
Custom filename: Null
Lower bound: Null
Upper bound: Null
Parameter identifiers: None
Initial structure: [1, 1]
Final structure: [1, 1]
Only sudden: False
Dynamics: [Sud, Lin, Exp]
No migrations: False
Symmetric migrations: False
Migration masks: Null
Selection: False
Dominance: False
Split fractions: True
Inbreeding: False
Ancestral size as parameter: False
Upper bound of first split: Null
Upper bound of second split: Null
Local optimizer: optimize_log_lbfgsb
Print models' code every N iteration: 0
Model plot engine: demes
Draw models every N iteration: 0
Units of time in drawing: years
Vmin: 1
Silence: False
Verbose: 1
Number of repeats: 3
Number of processes: 3
Resume from: Null
Only models: False
```

Finally, GADMA was run with `gadma --params gadma_folded_param_file`

## Historical Ne with StairwayPlot2
Inferring population demographic history from the SFS. Like with GADMA, [StairwayPlot2v2.1.1]() uses a parameter file outlining requirements for executing the programme.
```
# input setting
popid: AU
nseq: 38
L: Total number of observed nucleic sites
whether_folded: true
SFS:
#smallest_size_of_SFS_bin_used_for_estimation:
#largest_size_of_SFS_bin_used_for_estimation: 19
pct_training: 0.67
nrand: 7 15 22 28
project_dir:
stairway_plot_dir:
ninput: 
#random_seed: 
#output_setting
mu
year_per_generation: 3
#plot setting
plot_title: Two-epoch_fold
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12
```
Then run with `java -cp stairway_plot_es Stairbuilder blueprint_file`