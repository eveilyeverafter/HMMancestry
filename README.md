Package HMMancestry
=============

![Image](/logo.pdf?raw=true)


The package HMMancestry was designed to infer recombination breakpoints, hotspots, and coldspots 
in high-throughput, next-generation sequence data, even when sequencing coverage is relatively 
low. Currently, HMMancestry handles three types of data:

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
program like bowtie2 or bwa. These read counts are used in HMMancestry to estimate recombination 
hotspots, average length of a gene conversion tracts, and more.

3) Diploid data with read counts for alternative snps. This is the newest part of HMMancestry and is the part that is under the most development currently.
* Can be used to infer local ancestry of a number of model organisms
* Advanced filial generations can be used to increase detection of hotspots
* NCO and CO-associated gene conversion tracts cannot be determined with diploids since the algorithm doesn't use phase information.

To install compiled version:
* download HMMancestry_*.tar.gz
* open the command line and navigate to the location of the package
* use R CMD INSTALL HMMancestry_*.tar.gz to install (replacing * with the actual version number) 

To build from source:
* Have the R package devtools and roxygen2 installed
* Have git installed on your machine
* Open Terminal, pull down the repo, and cd into the directory
* type the following:

```
make roxygen
make build
make install
make doc
```
