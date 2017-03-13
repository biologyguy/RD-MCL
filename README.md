[![Build Status](https://travis-ci.org/biologyguy/RD-MCL.svg?branch=master)](https://travis-ci.org/biologyguy/RD-MCL)
[![Coverage Status](https://img.shields.io/coveralls/biologyguy/RD-MCL/master.svg)](https://coveralls.io/github/biologyguy/RD-MCL?branch=master)
<p align="center"><a href="https://github.com/biologyguy/RD-MCL/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/RD-MCL/master/rdmcl/images/rdmcl-logo.png" width=70%/></a></p>

# Recursive Dynamic Markov Clustering
### A method for identifying hierarchical orthogroups among homologous sequences
'Orthology' is a term that was coined to describe 'homology via speciation'[¹](References), which is now a concept
 broadly used as a predictor of shared gene function among species[²](References). Orthology, like homology in general,
 is a binary state; either two genetic elements are orthologous or they are not — there is no in-between. To make life
 complicated, however, the relationship among orthologs is not necessarily transitive. A simple example of this is a
 case where gene A is orthologous to gene B, gene B is orthologous to gene C, but gene C is not orthologous to gene A.
 This is further amplified by technical constraints, because it can be difficult to reconstruct the evolutionary
 history of genetic elements, especially when deletions have occurred (these and other challenges are well reviewed by
 Koonin[³](References)). 
 


In essence, RD-MCL is an extension of conventional [Markov clustering](http://micans.org/mcl/)-based orthogroup
 prediction algorithms like [OrthoMCL](http://orthomcl.org/orthomcl/). 
 
 
 The key innovations MCL parameters (inflation 
 value and edge similarity threshold) are selected 
 by the software to maximize the size of clusters while minizing the
 number of paralogs in each cluster. Furthermore, orthogroups are 
 broken up recursively, allowing for a more fine-grained final result.
 
The software is undergoing extensive development at the moment, and not 
 all of the code it needs to run is available in the repo. I 
 won't write up a proper HowTo until RD-MCL moves out of alpha, but I'm
 happy to help you get up and running with if you're convinced
 you need RD-MCL in your life this very second. 

If you have any comments, suggestions, or concerns, feel free to create an issue in the issue tracker or to get in 
 touch with me directly at steve.bond@nih.gov

## Installation
RD-MCL is hosted on the Python Package Index, so the easiest way to get the software and most dependencies is via `pip`:

```bash
    $: pip install rdmcl
    $: rdmcl -setup
```

## Dependencies not installed by `pip`
The following will need to be installed separately and placed somewhere in your PATH. 

* [MAFFT](http://mafft.cbrc.jp/alignment/software/)
* [PSIPRED](http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/)

If you are using Conda, then the following should get you up and running:

```bash
    $: conda install -c biocore psipred mafft
```

## References
¹ Fitch, W. M. [Distinguishing homologous from analogous proteins](https://doi.org/10.2307/2412448).
 _Systemat. Zool._ **19**, 99–106 (1970).

² Gabaldón, T. and Koonin, E. V. [Functional and evolutionary implications of gene
 orthology](https://doi.org/10.1038/nrg3456). _Nature reviews. Genetics._ **14**, 360-366 (2013).

³ Koonin, E. V. [Orthologs, paralogs, and evolutionary genomics](https://doi.org/10.1146/annurev.genet.39.073003.114725).
 _Annual review of genetics._ **39**, 309-338 (2005).