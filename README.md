[![Build Status](https://travis-ci.org/biologyguy/RD-MCL.svg?branch=master)](https://travis-ci.org/biologyguy/RD-MCL)
[![Coverage Status](https://img.shields.io/coveralls/biologyguy/RD-MCL/master.svg)](https://coveralls.io/github/biologyguy/RD-MCL?branch=master)
# RD-MCL
### Recursive Dynamic Markov clustering

In essence, RD-MCL is an extension of conventional Markov 
 clustering-based orthogroup prediction algorithms like OrthoMCL. MCL
 parameters (inflation value and edge similarity threshold) are selected 
 by the software to maximize the size of clusters while minizing the
 number of paralogs in each cluster. Furthermore, orthogroups are 
 broken up recursively, allowing for a more fine-grained final result.
 
The software is undergoing extensive development at the moment, and not 
 all of the code it needs to run is available in the repo. I 
 won't write up a proper HowTo until RD-MCL moves out of alpha, but I'm
 happy to help you get up and running with if you're convinced
 you need RD-MCL in your life this very second. Just send me an email: 
 steve.bond@nih.gov

## Dependencies
* [MAFFT](http://mafft.cbrc.jp/alignment/software/)
* [PSIPRED](http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/)
* [BuddySuite](https://github.com/biologyguy/BuddySuite)
* [pandas](https://pandas-docs.github.io/pandas-docs-travis/install.html)

