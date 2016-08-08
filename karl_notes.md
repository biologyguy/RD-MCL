# Generating mock sequences for testing on OrthoMCL

## Tools
We first tried using the evolver tool in [PAML 4.8](http://abacus.gene.ucl.ac.uk/software/paml.html) to generate amino
acid sequences off of a base tree.

We called evolver with the following:

    $: evolver 7

using MCaa.dat as the datafile. Within MCaa.dat, we set the number of taxa to 5, the sequence length to 200, and 1
replicate. The model used was empirical (3), with the 'wag' model located in dat/wag.dat within the paml directory.
We ran it on a wide range of gamma (shape and # categories) parameters to find an optimal evolution rate for
reproducibility, producing sequences which produced trees with low robinson-foulds distance values when compared to the
base tree. However, we found that these "optimal values" were only optimal for the given tree, and differred vastly for
other trees. At the same time, other literature, like

    10.1093/gbe/evw005

and

    10.1371/journal.pone.0105015

used a tool called [indel-seq-gen](https://github.com/cstrope/indel-seq-gen) to simulate sequence evolution. However,
we had significant issues trying to compile the software and had to modify it in order to build it correctly. The
modified version is available [here](https://github.com/biologyguy/indel-seq-gen). We initially ran it using this as the
seed tree

```
(((C:1.2815,D:0.7178):1.1211,(F:0.7311,J:0.6326):2.839):1.0542,(H:1.5417,I:1.4833):2.7087);
```

with the following parameters

   $: indel-seq-gen -m JTT -l 200 -b 70 -o n -e test < ml_tree.newick

This took 10 hours to run, and produced sequences that looked pretty much random and that weren't particularly useful.
At this point we gave up and decided to write our own program...

Which we promptly halted after Steve tweeted out asking if anyone knew of any good software for simulating sequence
evolution, which prompted several responses. Stephanie J. Spielman, the author of one such piece of software suggested
we use a python package she co-wrote, entitled pyvolve

    10.1371/journal.pone.0139047

which enables the simulation of nucleotide and AA sequence evolution with a base tree and seed sequence. However, it
does not simulate indels, because Spielman claims:

```
At this time, I'm not convinced we understand how to model subs+indels together cohesively. Til then, sims are 2 unrealistic
```

Upon testing pyvolve with a Ctenophore tree, we got back sequence alignments that were very similar, and when run
through FastTree, produced a tree with very high similarity to the original tree.

## Generating trees for sequence evolution simulation
Using biopython's Phylo package, we wrote a script that generates a random species tree, and a random gene tree, then
replaces all of the gene tree's leaf nodes with copies of the species tree, allowing us to simulate "perfect"
orthogroups.