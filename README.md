[![Build Status](https://travis-ci.org/biologyguy/RD-MCL.svg?branch=master)](https://travis-ci.org/biologyguy/RD-MCL)
[![Coverage Status](https://img.shields.io/coveralls/biologyguy/RD-MCL/master.svg)](https://coveralls.io/github/biologyguy/RD-MCL?branch=master)
[![PyPi version](https://img.shields.io/pypi/v/rdmcl.svg)](https://pypi.python.org/pypi/rdmcl)
<p align="center"><a href="https://github.com/biologyguy/RD-MCL/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/RD-MCL/master/rdmcl/images/rdmcl-logo.png" width=70%/></a></p>

# Recursive Dynamic Markov Clustering
### A method for identifying hierarchical orthogroups among homologous sequences
'Orthology' is a term that was coined to describe 'homology via
 speciation'[¹](https://github.com/biologyguy/RD-MCL#references), which is now a concept broadly used as a predictor of 
 shared gene function among species[²⁻³](https://github.com/biologyguy/RD-MCL#references). From a systematics 
 perspective, orthology also represents a natural schema for classifying/naming genes coherently. As we move into the 
 foreseeable future, when the genomes of all species on the planet have been sequences, it will be important to catalog 
 the evolutionary history of all genes and label them in a rational way. Considerable effort has been made to 
 programmatically identify orthologs, leading to excellent software solutions and several large public databases for 
 genome-scale predictions. What is currently missing, however, is a convenient method for fine-grained analysis of 
 specific gene families. 
  
In essence, RD-MCL is an extension of conventional [Markov clustering](http://micans.org/mcl/)-based orthogroup
 prediction algorithms like [OrthoMCL](http://orthomcl.org/orthomcl/), with three key differences:

1) The similarity metric used to describe the relatedness of sequences is based on multiple sequence alignments, not
 pair-wise sequence alignments or BLAST. This significantly improves the quality of the information available to the 
 clustering algorithm.

2) The appropriate granularity of the Markov clustering algorithm, as is controlled by the 'inflation factor' and 'edge
 similarity threshold', is determined on the fly. This is in contrast to almost all other methods, where default 
 parameters are selected at the outset and imposed indiscriminately on all datasets.

3) Differences in evolutionary rates among orthologous groups of sequences are accounted for by recursive rounds of
 clustering.

### Getting started

[Click here for a full use-case tutorial](https://github.com/biologyguy/RD-MCL/wiki/Tutorial)

RD-MCL is hosted on the Python Package Index, so the easiest way to get the software and most dependencies is via `pip`:

```bash
$: pip install rdmcl
$: rdmcl -setup
```

The program will complain if you don't run '-setup' before the first time you use it, so make sure you do that.

The input for RD-MCL is a sequence file in 
[any of the many supported formats](https://github.com/biologyguy/BuddySuite/wiki/Supported-formats), where the name of
 each sequence is prefixed with an organism identifier. For example:

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

The above is a few sequences from [KOG0001](https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=KOG0001), coming
 from Arabidopsis (ath), C. Elgans (cel), Human (hsa), and yeast (sce). Note the hyphen (-) separating each identifier 
 from the gene name. This is important! Make sure there are no spurious hyphens in any of the gene names, and if you 
 can't use a hyphen for some reason, set the delimiting character with the `-ts` flag.

Once you have your sequences named correctly, simply pass it into rdmcl:

```
$: rdmcl your_seq_file.fa
```

A new directory will be created which will contain all of the accoutrement associated with the run, including a
 'final_clusters.txt' file, which is the result you'll probably be most interested in.

There are several parameters you can modify; use `:$ rdmcl -h` to get a listing of them. They are also individually
 described in the [wiki](wiki).

## Video of Evolution 2017 talk
I discuss the rationale and high level implementation details of RD-MCL

[![Alt text](https://img.youtube.com/vi/52STQpKv8j4/0.jpg)](https://www.youtube.com/watch?v=52STQpKv8j4)

## Distributing RD-MCL on a cluster
RD-MCL will parallelize creation of all-by-all graphs while searching MCL parameter space. Once a graph has been
 created it is saved in a database, thus preventing repetition of the 'hard' work if/when the same cluster is identified
 again at a later time. This means that the computational burden of a given run will tend to be high at the beginning of 
 that run and decrease with time.
 
To spread the work out across multiple nodes during the 'hard' part, launch workers with the 
 [launch_worker](https://github.com/biologyguy/RD-MCL/wiki/launch_worker) script
 bundled with RD-MCL:
 
    $: launch_worker --workdb <path/to/desired/directory>
 
By default, `launch_worker` will use all of the cores it can find, so either sequester the entire node or pass in the
 `--max_cpus` flag to restrict it. I have run as many as 100 workers at a time, but be aware that this sort of pressure 
 can lead to some instability (i.e., lost jobs from the queue and frozen master threads). Twenty workers is usually safe. 

Next, launch RD-MCL with the `--workdb` flag set to the same path you specified for `launch_worker`:
 
    $: rdmcl --workdb <path/to/same/directory/as/launch_worker>
 
RD-MCL will now send its expensive all-by-all work to a queue and wait around for one of the workers to do the 
 calculations. You can keep track of how busy the workers are by running the 
 [monitor script](https://github.com/biologyguy/RD-MCL/wiki/monitor_dbs) in the same directory as the workers:

``` bash
$: monitor_dbs

Press return to terminate.
#Master  AveMhb   #Worker  AveWhb   #queue   #subq   #proc   #subp   #comp   #HashWait #IdWait  ConnectTime
29       19.0     16       51.0     1        362     22      12      29      25        25       0.01
```

Also, you can send an arbitrary number of RD-MCL jobs to the same worker pool, no problem.

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