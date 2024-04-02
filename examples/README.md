# CycleViz examples

## Simple - GBM39 ecDNA
Data derived from GBM39 glioblastoma cell line published in Turner et al., https://pubmed.ncbi.nlm.nih.gov/28178237/

This simple example does not use a YAML file to specify plot parameters.

1. Enter GBM39 example directory
>`cd $CV_SRC/examples/GBM39`
 
2. Call CycleViz 
>`$CV_SRC/CycleViz.py --ref hg19 --cycles_file FF-4_amplicon1_cycles.txt --cycle 1 -g FF-4_amplicon1_graph.txt --rotate_to_min --figure_size_style small`

Recall from the documentation that `--figure_size_style small` rescales the image to enable better sizing of the image at reduced scales.

## Complex - T41 ecDNA
Data derived from oropharyngeal squamous cell carinoma with hybrid human-viral ecDNA from Pang et al., https://pubmed.ncbi.nlm.nih.gov/34548317/

1. Enter the T41 example directory
>`cd $CV_SRC/examples/T41`

2. Call CycleViz

**Note that the order in which the feature YAML files are provided specifies the order in which the interior tracks appear (first is innermost, last is outermost).**
>`$CV_SRC/CycleViz.py --input_yaml_file t41_CV.yaml --feature_yaml_list T41_link_track.yaml T41_wgs_track.yaml T41_rna_track.yaml`

