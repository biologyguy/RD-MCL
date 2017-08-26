|Build_Status| |Coverage_Status| |PyPi_version|

|RDMCL|

--------------

Recursive Dynamic Markov Clustering
===================================

A method for identifying hierarchical orthogroups among homologous sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'Orthology' is a term that was coined to describe 'homology via speciation'¹, which is now a concept broadly used as a predictor of shared gene function among species²⁻³. From a systematics perspective, orthology also represents a natural schema for classifying/naming genes coherently. As we move into the foreseeable future, when the genomes of all species on the planet have been sequences, it will be important to catalog the evolutionary history of all genes and label them in a rational way. Considerable effort has been made to programmatically identify orthologs, leading to excellent software solutions and several large public databases for genome-scale predictions. What is currently missing, however, is a convenient method for fine-grained analysis of specific gene families.

In essence, RD-MCL is an extension of conventional `Markov clustering <http://micans.org/mcl/>`_-based orthogroup prediction algorithms like `OrthoMCL <http://orthomcl.org/orthomcl/>`_, with three key differences:

1) The similarity metric used to describe the relatedness of sequences is based on multiple sequence alignments, not pair-wise sequence alignments or BLAST. This significantly improves the quality of the information available to the clustering algorithm.
2) The appropriate granularity of the Markov clustering algorithm, as is controlled by the 'inflation factor' and 'edge similarity threshold', is determined on the fly. This is in contrast to almost all other methods, where default parameters are selected at the outset and imposed indiscriminately on all datasets.
3) Differences in evolutionary rates among orthologous groups of sequences are accounted for by recursive rounds of clustering.


Getting started
~~~~~~~~~~~~~~~

`Click here a full use-case tutorial <https://github.com/biologyguy/RD-MCL/wiki/Tutorial>`_

RD-MCL is hosted on the Python Package Index, so the easiest way to get the software and most dependencies is via `pip`:

.. code:: text

  $: pip install rdmcl
  $: rdmcl -setup


The program will complain if you don't run '-setup' before the first time you use it, so make sure you do that.

The input for RD-MCL is a sequence file in `any of the many supported formats <https://github.com/biologyguy/BuddySuite/wiki/SB-Screw-formats#format--str->`_, where the name of each sequence is prefixed with an organism identifier. For example:

.. code:: text

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


The above is a few sequences from `KOG0001 <https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=KOG0001>`_, coming from Arabidopsis (ath), C. Elgans (cel), Human (hsa), and yeast (sce). Note the hyphen (-) separating each identifier from the gene name. This is important! Make sure there are no spurious hyphens in any of the gene names, and if you can't use a hyphen for some reason, set the delimiting character with the `-ts` flag.

Once you have your sequences named correctly, simply pass it into rdmcl:

.. code:: text

  $: rdmcl your_seq_file.fa


A new directory will be created which will contain all of the accoutrement associated with the run, including a 'final_clusters.txt' file, which is the result you'll probably be most interested in.

There are several parameters you can modify; use `:$ rdmcl -h` to get a listing of them. Things are still under development so I haven't written a wiki yet, but I'd be overjoyed to get feedback if you are confused by anything. Please do not hesitate to email me!

Contact
-------

Any comments you have would be really appreciated. Please feel free to
add issues in the GitHub issue tracker or contact Steve Bond (lead
developer) directly at steve.bond@nih.gov.

.. |Build_Status| image:: https://travis-ci.org/biologyguy/RD-MCL.svg?branch=master
   :target: https://travis-ci.org/biologyguy/RD-MCL
.. |Coverage_Status| image:: https://img.shields.io/coveralls/biologyguy/RD-MCL/master.svg
   :target: https://coveralls.io/github/biologyguy/RD-MCL?branch=master
.. |PyPi_version| image:: https://img.shields.io/pypi/v/rdmcl.svg
   :target: https://pypi.python.org/pypi/rdmcl
.. |RDMCL| image:: https://raw.githubusercontent.com/biologyguy/RD-MCL/master/rdmcl/images/rdmcl-logo.png
   :target: https://github.com/biologyguy/RD-MCL/wiki
   :height: 200 px