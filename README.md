Package HMMancestry (dev)
=============

![alt tag](/FBgenotyper2.pdf)


The package fbgenotyper was designed to infer recombination breakpoints, hotspots, and coldspots 
in high-throughput, next-generation sequence data, even when sequencing coverage is relatively 
low. Currently, fbgenotyper handles three types of data:

1) Haploid yeast data where each of four spores of a tetrad are genotyped.
* Low throughput for the experimenter since it requires tetrad dissection.
* Can more precisely estimate different types of recombination events like crossover, non-crossover,
and gene conversion tracts.

2) Haploid yeast data where each spore comes from a different tetrad.
* This data type can be easily observed by the experimenter since spores can be liberated en masse 
(e.g., enzymatically breaking up tetrad and vortexing spores).

* Does not have the resolution to distinguish between different types of recombination since knowledge 
of other spores in the tetrad are unknown.

In either case, two haploid parent strains -- each sequenced at high coverage to accurately call 
bi-allelic snps between them -- are mated forming F1 zygotes. These diploids are sporulated to 
induce meiosis and the resulting haploid spores are isolated, grown in colonies, and DNA is 
extracted and sequenced using which ever method the experimenter desires. Read counts for each 
individually barcoded and sequenced spore are estimated from the parental snps using a mapping 
program like bowtie2 or bwa. These read counts are used in fbgenotyper to estimate recombination 
hotspots, average length of a gene conversion tracts, and more.

3) Diploid data with read counts for alternative snps. This is the newest part of fbgenotyper and is the part that is under the most development currently. Stay tuned for updates. 



To build:
```
make roxygen
make build
make install
```
