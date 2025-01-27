# CycleViz

Latest version: **0.2.0**.

Visualize outputs of [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect/), [CoRAL](https://github.com/AmpliconSuite/CoRAL), [ecSimulator](https://github.com/AmpliconSuite/ecSimulator), or [AmpliconReconstructor](https://github.com/jluebeck/AmpliconReconstructor)  in Circos-style images. 

CycleViz can also produce circular visualizations of general genomic regions using only a bed file. 
Supports hg19, GRCh37, hg38, and mm10. CycleViz is implemented in python and compatible with both python2 and python3. 
CycleViz has been tested on Ubuntu 16.04+, and macOS 10+.

**Examples**: Left, visualization of a simple AA-derived structure with multiple interior data tracks. Right, a cycles file visualization with Bionano data, generated from output of AmpliconReconstructor. 

<img src="images/T41_example.png" width=350px><img src="images/exampleAR.png" width=375px>


## Installation

Requires matplotlib version 2.0.0 or higher, intervaltree and pyyaml python modules.


To install relevant python packages. 
```bash
pip install intervaltree 'matplotlib>=2.0.0' numpy pyyaml
```
or
```bash
conda install intervaltree 'matplotlib-base>=2.0.0' numpy pyyaml
```

Then download the CycleViz repo
```
git clone https://github.com/AmpliconSuite/CycleViz
```

**Recommended: Setting a bash variable for CycleViz** 

After cloning the CycleViz repo consider running the following to add a line to your `.bashrc` file:

```bash
cd CycleViz/
echo export CV_SRC=$PWD >> ~/.bashrc
source ~/.bashrc
```

**Optional: Install Arial font in Linux**

CycleViz defaults to Arial font and will fall back to Deja Vu Sans if Arial is not present. Arial is already included on macOS. To get the Microsoft fonts on Ubuntu, do
```bash
sudo apt-get install ttf-mscorefonts-installer fontconfig
sudo fc-cache -f  # rebuilds the font cache
```
Then reset the matplotlib font cache using the instructions [here](https://stackoverflow.com/questions/37920935/matplotlib-cant-find-font-installed-in-my-linux-machine).

Alternatively, you can do this process using conda by running

`conda install mscorefonts` 

then launch python and do
```python
import matplotlib.font_manager
matplotlib.font_manager._load_fontmanager(try_read_cache=False)
```

## Usage

**The [examples directory](https://github.com/jluebeck/CycleViz/tree/master/examples) provides two worked examples and sample datasets as a jumping-off point for users.**

A very basic visualization would go like so:
>`$CV_SRC/CycleViz.py -g sample_amplicon1_graph.txt --cycles_file sample_amplicon1_cycles.txt --cycle 1 --ref GRCh38`

### Command line arguments

There are two modes of visualization, CycleViz and LinearViz.

using the `--help` command will output a specific description of each argument below.

**The following arguments can be specied individually on the command line or put into a yaml file passed to CycleViz with `--input_yaml_file`.**

#### Specifying the genomic structure
The structure of the visualized regions can be specified with either an AA cycles file, or a bed file listing the regions (called a 'structure bed')

| Argument | Description                                                                                                                                                                                |
| :--- |:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--cycles_file [filename]`   | Specify structure with an AA-formatted cycles file                                                                                                                                         |
| `--structure_bed [filename]` | Specify structure with a bed file listing the segments and their orientations. Some predefined reference genome structures are available by specifying `hg19`, `GRCh37`, `hg38`, `GRCh38`. |

**If `--cycles_file` is used for the structure, the following two additional arguments should be specified**

| Argument | Description |
| :---  |  :----  |
| `--cycle [int]` \[required\]  | Cycle ID number from AA cycles file to use for visualization |
| `-g `/`--graph [filename]` \[optional\] | AA graph file for the amplicon |

#### Annotation arguments

| Argument                                                  | Default | Description                                                                                                                                                                                                                                              |
|:----------------------------------------------------------| :----: |:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--ref [hg19, GRCh37, hg38, GRCh38, mm10, GRCh38_viral]`  | | Reference build                                                                                                                                                                                                                                          |
| `--sname [filename prefix]` \[optional\]                  | | Output filename prefix                                                                                                                                                                                                                                   |
| `--gene_subset_file [filename or "BUSHMAN"]` \[optional\] | | File giving a list of genes show in the plot (as GFF or textfile with RefGene name in last column), or use the [Bushman oncogene list](http://www.bushmanlab.org/links/genelists). If no file is specified, all overlapping RefGene genes will be shown. |
| `--gene_subset_list [string] [string] ...` \[optional\]   | | Alternative to `--gene_subset_list`. Directly specify a list of RefGene gene names to show in plot.                                                                                                                                                      | 


#### Specifying optional properties related to plot appearance
| Argument                                                       |  Default   | Description                                                                                                                                                                                               |
|:---------------------------------------------------------------|:----------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--annotate_structure ['genes' /path/to/chromosome/coloring/]` | `'genes'`  | Will either show locations of genes or user can give a filename of colors to assign to specific regions. There are cytoband coloring files in the `resources/` directory that can be set as the argument. |
| `--structure_color ['auto', matplotlib color]`                 |  `'auto'`  | Either use default coloring of chromosomes (`auto`), or specify a single matplotlib color to use on every chromosome. This will be plotted behind `annotate_structure`.                                   | 
| `--figure_size_style ["normal", "small"]`                      | `"normal"` | Produce normally scaled figure or `small` figure rescaled for small image size.                                                                                                                           | 
| `--label_segs ["names", "numbers"]`                            |            | Print the segment ID & direction (`numbers`) or chromosome name (`names`) or a custom list of labels (specified as multiple space separated arguments) under structure segment.                           |
| `--gene_fontsize [float]`                                      |     7      | Gene name fontsize                                                                                                                                                                                        |
| `--gene_highlight_list [string] [string] ...`                  |            | List of RefGene gene names to give alternate color (default red)                                                                                                                                          |
| `--print_dup_names`                                            |            | Print gene name each time it is split across segments. Default, print gene name once if split across multiple segments.                                                                                   | 
| `--tick_fontsize`                                              |     7      | Fontsize for coordinate ticks or endpoint coordinates.                                                                                                                                                    |
| `--tick_type`                                                  |   `'ends'`    | Represent ticks only at segment ends, or throughout ('ends' is default).                                                                                                                                  |
| `--hide_chrom_color_legend [True/False]`                       |  `False`   | Do not print a map of color to chromosome name on the left side. Perhaps set to `True` if showing more than ~10 chroms.                                                                                   |
| `--rotate_to_min`                                              |            | Rotate the plot such that the smallest genomic coordinate resides at the 0 degree position in the first quadrant (3 o'clock).                                                                             |
| `--no_PDF`                                                     |            | Do not save a PDF version of the plot.                                                                                                                                                                    |

#### Specifying properties related to interior data track features
| Argument                                                   | Default  | Description                                                                                                                                                                                                                       |
|:-----------------------------------------------------------|:--------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--feature_yaml_list [filename] [filename] ...`            |          | List of interior feature track yaml files (**ordered from inner to outer**).                                                                                                                                                      |
| `--interior_segments_cycle [filename]`                     |          | An AA cycles-formatted file indicating regions to draw an inlaid structure beneath the reference (e.g. a transcript structure from a rearranged genome). Subset segments should be ordered consistently with the outer structure. | 
| `--interior_segments_cycle_connect_width ["full", "thin", "none"]` | `'full'` | Defines how the segments should be connected visually in the interior cycle (full height (`full`) or `auto`, which is just a small connector line).                                                                               |
| `--interior_segments_cycle_matching ["exact", "relaxed"]`  |          | Specifying matching of interior cycle to outer cycle should be `exact` or `relaxed` (relaxed for local matches with shared start). Required if --interior_segments_cycle set                                                      |
| `--center_hole [float]`                                    |   1.25   | Radius of center hole in CycleViz plot where no data appears (if interior features are used).                                                                                                                                     |

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


### Formatting a feature yaml file
The arguments that can be specified inside a YAML file for each data track. 

| Argument                       |  Default   |                                                                                       Choices                                                                                       | Description                                                                                                                                                                                                                                                                                                                     |
|:-------------------------------|:----------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `name`                         |            |                                                                                                                                                                                     | Give a name to the feature you are plotting                                                                                                                                                                                                                                                                                     |
| `tracktype`                    | `standard` |                                                                          `['standard', 'links', 'rects']`                                                                           | The track will show numeric data (`standard`), links between points (`links`), or colored rectangles (`rects`).                                                                                                                                                                                                                 |
| `primary_feature_bedgraph`     |            |                                                           path to a file formatted as bed (or wig), with data in column 4                                                           | Required argument. This is the data that will be shown.                                                                                                                                                                                                                                                                         |
| `secondary_feature_bedgraph`   |            |                                                      path to a file formatted as bed (or wig), with position value in column 4                                                      | Put a second source of data on the same track.                                                                                                                                                                                                                                                                                  | 
| `primary_style`                |  `points`  |                                                                           `['points', 'lines', 'radial']`                                                                           | Plot primary data track as `points`, connected `lines`, or `radial` lines from the bottom of the track.                                                                                                                                                                                                                         |
| `secondary_style`              |  `points`  |                                                                           `['points', 'lines', 'radial']`                                                                           | Plot secondary data track as `points`, connected `lines`, or `radial` lines from the bottom of the track.                                                                                                                                                                                                                       |
| `rescale_by_secondary`         |  `False`   |                                                                        `[True, False, 'mean' (True), 'each']`                                                                        | Divide entries in the primary data track by entries in the secondary track's data, at the corresponding position. Uses mean if multiple overlapping.                                                                                                                                                                            |
| `rescale_secondary_to_primary` |  `False`   |                                                                                   `[True, False]`                                                                                   | Rescale secondary data track to use the same min/max height as primary data track. A separate axis in the legend will be created to show the corresponding true values at their scaled positions. Note this feature is experimental and may cause weird behaviors. We advise any data rescaling be done prior to using CycleViz |
| `rescale_by_count`             |  `False`   |                                                                                   `[True, False]`                                                                                   | Rescale the data track by the number of times that structure segment (identified by segment ID/name) appears.                                                                                                                                                                                                                   |
| `hide_secondary`               |  `False`   |                                                                                   `[True, False]`                                                                                   | Do not show the secondary data, however still applies all other rescaling operations specified.                                                                                                                                                                                                                                 |
| `indicate_zero`                |  `False`   |                                                                                   `[True, False]`                                                                                   | Plot a line in the data track corresponding to the value 0.                                                                                                                                                                                                                                                                     |
| `zero_facecolor`               |    `k`     |                                                                              a matplotlib named color                                                                               | Color of the line in the data track corresponding to the value 0.                                                                                                                                                                                                                                                               |
| `num_hlines`                   |    `5`     |                                                                                    integer >= 0                                                                                     | Number of horizontal grid lines in the track (only for `standard` tracktype).                                                                                                                                                                                                                                                   |
| `nice_hlines`                  |   `True`   |                                                                                   `[True, False]`                                                                                   | Automatically select reasonable min and max values for the gridlines. Turning off uses raw min and max from track as the hline min/max.                                                                                                                                                                                         |
| `grid_legend_fontsize`         |    `4`     |                                                                                     number > 0                                                                                      | Fontsize for gridline value ticks in track legend plot.                                                                                                                                                                                                                                                                         |
| `granularity`                  |    `0`     |                                                                                     number > 0                                                                                      | The frequency with which to draw points/breaks between lines. `0` indicates to use an automatic amount of granularity.                                                                                                                                                                                                          |
| `primary_smoothing`            |    `0`     |                                                                                     number > 0                                                                                      | Amount of smoothing (rolling average) to apply across points from primary data. This is applied after granularity is set. For mild smoothing try 50.                                                                                                                                                                            |
| `secondary_smoothing`          |    `0`     |                                                                                     number > 0                                                                                      | Amount of smoothing (rolling average) to apply across points from secondary data. This is applied after granularity is set. For mild smoothing try 50.                                                                                                                                                                          |
| `end_trim`                     |    `0`     |                                                                                     number > 0                                                                                      | Do not plot data within `[end_trim]` of the segment ends.                                                                                                                                                                                                                                                                       |
| `show_segment_copy_count`      |  `False`   |                                                                                   `[True, False]`                                                                                   | Plot a horizontal line indicating the number of copies this segment has in the structure.                                                                                                                                                                                                                                       |
| `segment_copy_count_scaling`   |    `1`     |                                                                                       number                                                                                        | (Used only if `show_segment_copy_count` is set) Since `show_segment_copy_count` may be indecipherable if the range is large, apply this value as a multiplicative constant to the copy number of the segment.                                                                                                                   |
| `primary_upper_cap`            |   `None`   |                                                                                       number                                                                                        | Set the following value as maximum value to allow in primary data. Assessed after rescaling. Values exceeding this parameter will be set to this parameter.                                                                                                                                                                     |
| `secondary_upper_cap`          |   `None`   |                                                                                       number                                                                                        | Set the following value as maximum value to allow in secondary data. Assessed after rescaling. Values exceeding this parameter will be set to this parameter.                                                                                                                                                                   |
| `primary_lower_cap`            |   `None`   |                                                                                       number                                                                                        | Set the following value as minimum value to allow in primary data. Assessed after rescaling. Values below this parameter will be set to this parameter.                                                                                                                                                                         |
| `secondary_lower_cap`          |   `None`   |                                                                                       number                                                                                        | Set the following value as minimum value to allow in secondary data. Assessed after rescaling. Values below this parameter will be set to this parameter.                                                                                                                                                                       |
| `linkpoint`                    |   `None`   |                                                                                `[None, 'midpoint']`                                                                                 | When drawing links, setting 'midpoint' will draw a line between the centers of the link bedpe entry start and end positions.                                                                                                                                                                                                    |
| `link_single_match`            |  `False`   |                                                                                   `[True, False]`                                                                                   | If a link object can go to multiple locations (because the segment is shown multiple times), draw it for all combinations of endpoints (`False`). If `True`, draw only closest pairing of start/end in the structure.                                                                                                           |
| `primary_kwargs`               |    `{}`    | [Line2D](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.lines.Line2D.html) or [Scatter](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.pyplot.scatter.html) **kwargs dict | **kwargs for the primary data. Will override CycleViz defaults. Requires user to know which **kwargs parent they are working with. `trackstyle: points` is Scatter, others are Line2D.                                                                                                                                          |
| `secondary_kwargs`             |    `{}`    | [Line2D](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.lines.Line2D.html) or [Scatter](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.pyplot.scatter.html) **kwargs dict | **kwargs for the secondary data. Will override CycleViz defaults. Requires user to know which **kwargs parent they are working with.                                                                                                                                                                                            |
| `link_kwargs`                  |    `{}`    |                               [Patch](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.patches.Patch.html#matplotlib.patches.Patch) **kwargs dict                                | **kwargs for link data. Will override CycleViz defaults. Requires user to know which **kwargs parent they are working with.                                                                                                                                                                                                     |
| `hline_kwargs`                 |    `{}`    |                               [Patch](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.patches.Patch.html#matplotlib.patches.Patch) **kwargs dict                                | **kwargs for horizontal gridlines. Will override CycleViz defaults.                                                                                                                                                                                                                                                             |
| `background_kwargs`            |    `{}`    |                                            [Patch](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.lines.Line2D.html) **kwargs dict                                             | **kwargs for background of the track. Will override CycleViz defaults.                                                                                                                                                                                                                                                          |


### Examples

**Please see the [examples folder](https://github.com/jluebeck/CycleViz/tree/master/examples) for worked examples with demo data.**
 
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

### Creating your own structure.bed file
If you would like to specify a collection of region of the genome to show please create a file formatted as follows

`chrom  pos1    pos2    strand  [connected to previous sequence (True/False)]`

e.g.

```
chr10  138300  850441  +   False
chr4    450220  602811  +   False    
```

The structure of the plot will be ordered by the structure you specified! 
