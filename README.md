[![Build Status](https://travis-ci.org/biologyguy/RD-MCL.svg?branch=master)](https://travis-ci.org/biologyguy/RD-MCL)
[![Coverage Status](https://img.shields.io/coveralls/biologyguy/RD-MCL/master.svg)](https://coveralls.io/github/biologyguy/RD-MCL?branch=master)
[![PyPi version](https://img.shields.io/pypi/v/rdmcl.svg)](https://pypi.python.org/pypi/rdmcl)
<p align="center"><a href="https://github.com/biologyguy/RD-MCL/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/RD-MCL/master/rdmcl/images/rdmcl-logo.png" width=70%/></a></p>

# Recursive Dynamic Markov Clustering
### A method for identifying hierarchical orthogroups among homologous sequences
'Orthology' is a term that was coined to describe 'homology via speciation'[¹](https://github.com/biologyguy/RD-MCL#references), which is now a concept broadly used as a predictor of shared gene function among species[²⁻³](https://github.com/biologyguy/RD-MCL#references). From a systematics perspective, orthology also represents a natural schema for classifying/naming genes coherently. As we move into the foreseeable future, when the genomes of all species on the planet have been sequences, it will be important to catalog the evolutionary history of all genes and label them in a rational way. Considerable effort has been made to programmatically identify orthologs, leading to excellent software solutions and several large public databases for genome-scale predictions. What is currently missing, however, is a convenient method for fine-grained analysis of specific gene families. 
  
In essence, RD-MCL is an extension of conventional [Markov clustering](http://micans.org/mcl/)-based orthogroup prediction algorithms like [OrthoMCL](http://orthomcl.org/orthomcl/), with three key differences:

1) The similarity metric used to describe the relatedness of sequences is based on multiple sequence alignments, not pair-wise sequence alignments or BLAST. This significantly improves the quality of the information available to the clustering algorithm.

2) The appropriate granularity of the Markov clustering algorithm, as is controlled by the 'inflation factor' and 'edge similarity threshold', is determined on the fly. This is in contrast to almost all other methods, where default parameters are selected at the outset and imposed indiscriminately on all datasets.

3) Differences in evolutionary rates among orthologous groups of sequences are accounted for by recursive rounds of clustering.

### Getting started

[Click here a full use-case tutorial](https://github.com/biologyguy/RD-MCL/wiki/Tutorial)

RD-MCL is hosted on the Python Package Index, so the easiest way to get the software and most dependencies is via `pip`:

```bash
$: pip install rdmcl
$: rdmcl -setup
```

The program will complain if you don't run '-setup' before the first time you use it, so make sure you do that.

The input for RD-MCL is a sequence file in [any of the many supported formats](https://github.com/biologyguy/BuddySuite/wiki/SB-Screw-formats#format--str-), where the name of each sequence is prefixed with an organism identifier. For example:

```fasta
>ath-At4g02970
MNVYIDTETGSSFSITIDFGETVLEIKEKIEKSQGIPVSKQILYLDGKALEDDLHKIDYM
ILFESRLLLRISPDADPNQSNEQTEQSKQIDDKKQEFCGIQDSSESKKITRVMARRVHNI
YSSLPAYSLDELLGPKYSATVAVGGRTNQVVQPTEQASTSGTAKEVLRDSDSPVEKKIKT
NPMKFTVHVKPYQEDTRMIHVEVNADDNVEELRKELVKMQERGELNLPHEAFHLLGLGSS
ETCPHQNRSEEPNQCPTILMSPHGLQAIVT
>cel-CE08215_2
QIFVKVLGVSYAFKIHREDTVFDIKNDIEHRHDIPQHSYWLSFSGKRLEDHCSIGDYNIQ
KSSTITMYFRSG
>cel-CE16986
MKATTVKENEVKDDRKLSLNEMLRKRCLQVKNTKMKNSSMPKFQYFVRLNGKTRTLNVNA
SDTVEQGKMQLCHNARSTRMSYGGKPLSDQITFGEYNISNNSTMDLHFRI
>hsa-Hs20473312
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQGKEGIPPDQQRLIFAGKQLEDGRTLSDYN
IQKESTLHLVLRLLVVLRKGRRSLTPLPRRISTRERRLSWLS
>sce-YDR139c
MIVKVKTLTGKEISVELKESDLVYHIKELLEEKEGIPPSQQRLIFQGKQIDDKLTVTDAH
LVEGMQLHLVLTLRGGN
```

The above is a few sequences from [KOG0001](https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=KOG0001), coming from Arabidopsis (ath), C. Elgans (cel), Human (hsa), and yeast (sce). Note the hyphen (-) separating each identifier from the gene name. This is important! Make sure there are no spurious hyphens in any of the gene names, and if you can't use a hyphen for some reason, set the delimiting character with the `-ts` flag.

Once you have your sequences named correctly, simply pass it into rdmcl:

```
$: rdmcl your_seq_file.fa
```

A new directory will be created which will contain all of the accoutrement associated with the run, including a 'final_clusters.txt' file, which is the result you'll probably be most interested in.

There are several parameters you can modify; use `:$ rdmcl -h` to get a listing of them. Things are still under development so I haven't written a wiki yet, but I'd be overjoyed to get feedback if you are confused by anything. Please do not hesitate to email me!

## Video of Evolution 2017 talk
I discuss the rationale and high level implementation details of RD-MCL

[![Alt text](https://img.youtube.com/vi/52STQpKv8j4/0.jpg)](https://www.youtube.com/watch?v=52STQpKv8j4)

## Running many RD-MCL jobs on a cluster (not relevant for single jobs)
RD-MCL will parallelize creation of all-by-all graphs while searching MCL parameter space. Once a graph has been
 created it is saved in a database, thus preventing repetition of the 'hard' work in case the same cluster is identified
 again at a later time. This means that the computational burden of a given run will tend to decrease with time and many
 cores may end up sitting idle as the job wears on. This could get you in trouble with your sys-admins and
 fellow users if running many RD-MCL jobs on shared resources, because it will look like you're tying up a bunch of
 CPUs without actually making use of them...
 
Instead, set up one or more dedicated worker nodes with the `launch_worker` script bundled with RD-MCL:
 
    $: launch_worker --workdb <path/to/desired/directory>
 
Make sure to sequester entire nodes for this script, as it will use all the cores it can find. Next, launch RD-MCL
 jobs with the `--max_cpus` flag set to 4, and the `--workdb` flag set to the same path you specified for 
 `launch_worker`:
 
    $: rdmcl --max_cpus 4 --workdb <path/to/same/directory/as/launch_worker>
 
RD-MCL will now send its expensive all-by-all work to a queue and wait around for one of the workers to do the 
 calculations. A good starting point is about one worker with ≥16 cores for every ten rdmcl jobs of <100 sequences.


## References
¹ Fitch, W. M. [Distinguishing homologous from analogous proteins](https://doi.org/10.2307/2412448).
 _Systemat. Zool._ **19**, 99–106 (1970).

² Gabaldón, T. and Koonin, E. V. [Functional and evolutionary implications of gene
 orthology](https://doi.org/10.1038/nrg3456). _Nature reviews. Genetics._ **14**, 360-366 (2013).

³ Koonin, E. V. [Orthologs, paralogs, and evolutionary genomics](https://doi.org/10.1146/annurev.genet.39.073003.114725).
 _Annual review of genetics._ **39**, 309-338 (2005).
 
## Contact
If you have any comments, suggestions, or concerns, feel free to create an issue in the issue tracker or to get in 
 touch with me directly at steve.bond@nih.gov