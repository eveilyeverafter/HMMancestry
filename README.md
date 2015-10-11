Package fbgenotyper (dev)
=============

![alt tag](https://github.com/tylerhether/fbgenotyper/blob/master/FBgenotyper2.pdf)


ATTN: fbgenotyper is in the middle of a complete overhaul as we optimize code for speed and improve functionality (e.g., added support for non-yeast species and diploids). For now please consider this package "experimental". 


The package fbgenotyper was designed to infer recombination breakpoints, hotspots, and coldspots 
in high-throughput, next-generation sequence data, even when sequencing coverage is relatively 
low. Currently, fbgenotyper handles two types of data, both from yeast matings:

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


To build:
```
make roxygen
make build
make install
```
