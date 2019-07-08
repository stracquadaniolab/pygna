# PyGNA: a Python framework for geneset network analysis

Current version: 0.4.2-dev

![build](https://circleci.com/gh/stracquadaniolab/baghera/tree/master.svg?style=svg)
![platform](https://anaconda.org/stracquadaniolab/baghera/badges/platforms.svg)
![anaconda](https://anaconda.org/stracquadaniolab/baghera/badges/version.svg)

PyGNA is a unified framework for network analysis of high-throughput experiment results. It can be used both as a standalone command line application or it can be included as a package in your own python code. 

For an overview of PyGNA functionalities check the infographic below, otherwise dive into the [Getting started](#getting-started) guide. For the complete API check the official [Documentation](#documentation)

![Infographic](docs/pygna_infographic-01.png)

## Installation

The easiest and fastest way to install `pygna` using `conda`:

    $ conda install -c stracquadaniolab -c bioconda -c conda-forge pygna

Alternatively you can install it through `pip`:

    $ pip install pygna

Please note, that `pip` will not install non Python requirements.

## Getting started

A typical `pygna` analysis consists of 3 steps:

1. Generate the RWR and SP matrices for the network you are using ( once they are generated, you won't need to repeat the same step again)
2. Make sure that the input genesets are in the right format. If a network uses entrez ID, and your file is in HUGO symbols, use the pygna utility for the name conversion.
3. Run the analysis you are interested into.
4. Once you have the output tables, you can choose to visualise one or more plots.

Otherwise you can check the snakemake pipeline for the full analysis of a geneset. 

We provide the data for a minimum working example in the zip folder named `min_working_example`.
The examples below show some basic analysis that can be carried out with pygna

### Example 1: Running pygna GNT analysis

Running `pygna` on this input as follows:

    $ cd ./your-path/min-working-example/
    $ pygna build-RWR-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5
    $ pygna analyse-RW barabasi.interactome.tsv  disgenet_cancer_groups_subset.gmt  interactome_RWR.hdf5  ../min-working-example/ interactome
    $ pygna pygna paint-datasets-stats interactome_table_RW.csv  ../min_working_example/ interactome

### Example 2: Running pygna GNA analysis
    
    $ cd ./your-path/min-working-example/
    $ #skip this step if the matrix is already computed
    $ pygna build-RWR-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5
    $ #The association analysis is run N x M times (N number of genesets, M number of pathways), we use only 100 permutations in this example to avoid long computations. Recommended is 1000
    $ pygna comparison-random-walk barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 ../min_working_example/ GO_cc_interactome -B GO_cc_subset.gmt -k --number-of-permutations 50 --show-results
    $ #If you don't include the --show-results flag at the comparison step, plot the matrix as follows
    $ pygna paint-comparison-RW GO_cc_interactome_table_association_RW.csv  ../min_working_example/ comparison_stats



## Documentation

The official documentation for `pygna` can be found on [readthedocs](https://pygna.readthedocs.io/).

## Authors

- Viola Fanfani (v.fanfani@sms.ed.ac.uk): lead developer.
- Giovanni Stracquadanio (giovanni.stracquadanio@ed.ac.uk)

## Citation

A unified framework for geneset network analysis.
Viola Fanfani and  Giovanni Stracquadanio
bioRxiv XX; doi: XX

## Issues

Please post an issue to report a bug or request new features.
