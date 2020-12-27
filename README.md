# CycleViz

Visualize outputs of [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect/) & [AmpliconReconstructor](https://github.com/jluebeck/AmpliconReconstructor) (AR) in Circos-style images. Supports hg19, hg38, and GRCh37. CycleViz is implemented in python and compatible with both python2 and python3. CycleViz has been tested on Ubuntu 16.04+ and MacOS 10+.

**Examples**: Left, a cycles file visualization without AR-reconstruction data. Right, a cycles file visualization with Bionano data. 

Seperate instructions [are available](#running-cycleviz-with-an-aa-generated-cycles-file) if using an AA cycles file without AR reconstruction.

<!---![AA example](images/exampleAA.png){:height="300px" width="300px"}
![AR example](images/exampleAR.png){:height="300px" width="300px"} --->

<img src="images/exampleAA.png" height="45%" width="45%"> <img src="images/exampleAR.png" height="45%" width="45%">


## Installation

Requires is matplotlib version 2.0.0 or higher and intervaltree python module, both of which are non-standard in Conda.

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

[optional] To get the Microsoft fonts on Ubuntu (CycleViz defaults to Arial font)
```
sudo apt-get install ttf-mscorefonts-installer
```

After cloning the CycleViz repo consider running the following to add a line to your `.bashrc` file:

`cd CycleViz/`

`echo export CV_SRC=$PWD >> ~/.bashrc`

`source ~/.bashrc`

## Usage

### Running CycleViz with an AA-generated cycles file
If using an AA-generated cycles file instead of AR, whose segments are merged across adjacent breakpoint graph segments, you will first need to unmerge the segments in the cycles file. We provide a script, `convert_cycles_file.py` which can be run on an AA cycles file and AA breakpoint graph. This script must be run to generate a BPG-converted cycles file, prior to running CycleViz. Many thanks to Siavash R. Dehkordi for providing this script.

### Command line arguments

There are two modes of visualization, circular and linear. Currently there are separate scripts to run each mode, but we intend to merge these in the future.

using the `--help` command will output a specific description of each argument below.

**The following arguments can be specied individually on the command line or put into a yaml file passed to CycleViz with `--input_yaml_file`.**

#### Specifying the genomic structure
The structure of the visualized regions can be specified with either an AA cycles file, or a bed file listing the regions (called a 'structure bed')


| Argument      | Description |
| :---        |    :----  | 
| `--cycles_file [filename]`      | Specify structure with an AA-formatted cycles file (converted to graph segment annotation using instructions above).       |
| `--structure_bed [filename]`   | Specify structure with a bed file listing the segments and their orientations |

**If `--cycles_file` is used for the structure, the following two additional arguments should be specified**

| Argument | Description |
| :---  |  :----  |
| `--cycle [int]` \[required\]  | Cycle ID number from AA cycles file to use for visualization |
| `-g `/`--graph [filename]` \[optional\] | AA graph file for the amplicon |

#### Annotation arguments

| Argument | Default | Description |
| :---  | :----: |  :----  |
| `--ref [hg19, GRCh37, hg38, GRCh38]` |  `hg19` | Reference coordinate system used |
| `--sname [filename prefix]` \[optional\] | | Output filename prefix |
| `--gene_subset_file [filename or "BUSHMAN"]` \[optional\] | | File giving a list of genes show in the plot (as GFF or textfile with RefGene name in last column), or use the [Bushman oncogene list](http://www.bushmanlab.org/links/genelists). If no file is specified, all overlapping RefGene genes will be shown. |
| `--gene_subset_list [string] [string] ...` \[optional\] | | Alternative to `--gene_subset_list`. Directly specify a list of RefGene gene names to show in plot.| 


#### Specifying optional properties related to plot appearance
| Argument      | Default | Description |
| :---        |    :----:   | :--- |
|`--figure_size_style ["normal", "small"]`| `"normal"` | Produce normally scaled figure or `small` figure rescaled for small image size. | 
| `--label_segs` | | Print the segment ID & direction under structure segment.|
| `--gene_fontsize [float]` | 7 | Gene name fontsize |
| `--gene_highlight_list [string] [string] ...` | | List of RefGene gene names to give alternate color (default red) |
| `--print_dup_names` | | Print gene name each time it is split across segments. Default, print gene name once if split across multiple segments.| 
| `--segment_end_ticks` | | Label exact coordinate endpoints of segments and do not show tick marks along segment. Default is off, and will print ticks of approx location (scaled by 10 kbp) along the segment. |
| `--tick_fontsize` | 7 | Fontsize for coordinate ticks or endpoint coordinates. |

#### Specifying properties related to interior data track features
| Argument      | Default | Description |
| :---        |    :----:   | :--- |
| `--feature_yaml_list [filename] [filename] ...` |  | List of interior feature track yaml files (**ordered from inner to outer**). |
| `--interior_segments_cycle [filename]` | | An AA cycles-formatted file indicating regions to draw an inlaid structure beneath the reference (e.g. a transcript structure from a rearranged genome). Subset segments should be ordered consistently with the outer structure. | 
| `--center_hole [float]` | 1.25 | Radius of center hole in CycleViz plot where no data appears (if interior features are used). |

#### Specifying properties related to Bionano data & AmpliconReconstructor output
| Argument      | Default | Description |
| :---        |    :----:   | :--- |
| `--om_alignments` | | Must be set to enable the Bionano/AR mode |
| `--om_segs [filename]` | | CMAP file of *in silico* reference segments |
| `-c`/`--contigs [filename]` | | CMAP file of the OM contigs | 
| `-i/--path_alignment [filename]` | | AR path alignment file |

#### Arguments specific to LinearViz
| Argument      | Default | Description |
| :---        |    :----:   | :--- |
| `--reduce_path [int] [int]` | `0 0` | Trim the following number of segments from path beginning and end, respectively. |


#### Examples
 
For circular visualizations

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
                    [--gene_subset_file GENE_SUBSET_FILE | --gene_subset_list GENE_SUBSET_LIST [GENE_SUBSET_LIST]
                    `

The `--gene_subset_file` or `--gene_subset_list` arguments can be used to reduce the genes shown in the visualization. By default all genes (except ncRNAs, lncRNAs, microRNAs, etc.) will be shown. An cancer-related genelist file is provided: `Bushman_group_allOnco_May2018.tsv`, provided on the [Bushman Lab website](http://www.bushmanlab.org/links/genelists). This can be specified simply by setting `--gene_subset_file Bushman`.

Note that the structure bed or the cycles file/cycle number (or "path number", in the linear case) are the only required arguments. It is highly recommended to use the Bushman oncogene file for the gene_subset_file to make more readable plots

