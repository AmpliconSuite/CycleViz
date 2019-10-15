## CycleViz

###Installation
Visualize outputs of AmpliconReconstructor in Circos-style images. 

Requires python2 and matplotlib version 2.0.0 or higher and intervaltree python module. 

To check your matplotlib version in python, type
```
import matplotlib
print matplotlib.__version__
```

To upgrade to latest matplotlib from command line, do 
```
pip install --upgrade matplotlib
```

To install intervaltree python package. 
```
pip install intervaltree
```

###Usage
There are two modes of visualization, circular and linear. Currently there are separate scripts to run each mode, but we intend to merge these in the future.

using the `--help` command will output a specific description of each argument below.

For circular visualizations do 

`CycleViz.py [-h] [--om_alignments] --cycles_file CYCLES_FILE --cycle
                   CYCLE [-c CONTIGS] [-s SEGS] [-g GRAPH] [-i PATH_ALIGNMENT]
                   [--sname SNAME] [--rot ROT] [--label_segs]
                   [--gene_subset_file GENE_SUBSET_FILE]
                   `

For linear visualizations do 
` LinearViz.py [-h] [--om_alignments] [-s SEGS] [-g GRAPH] [-c CONTIGS]
                    --cycles_file CYCLES_FILE --path PATH [-i PATH_ALIGNMENT]
                    [--sname SNAME] [--label_segs]
                    [--reduce_path REDUCE_PATH REDUCE_PATH]
                    [--gene_subset_file GENE_SUBSET_FILE | --gene_subset_list GENE_SUBSET_LIST [GENE_SUBSET_LIST
                    `


Note that the cycles file and cycle number (or "path number", in the linear case) are the only required arguments. It is highly recommended to use the Bushman oncogoene file for the gene_subset_file to make more readable plots