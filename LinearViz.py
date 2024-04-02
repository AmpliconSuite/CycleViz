#!/usr/bin/env python

import argparse
from collections import defaultdict
import copy
import os
import sys

from ast import literal_eval as make_tuple
import matplotlib
matplotlib.use('Agg')  # this import must happen immediately after importing matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np

from bionanoUtil import *
from convert_cycles_file import *
import VizUtil as vu
from _version import __version__


rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['pdf.fonttype'] = 42


CV_RESOURCES = os.path.dirname(os.path.abspath(__file__)) + "/resources/"
print("Using resources in " + CV_RESOURCES)

seg_spacing = 0.009
bar_width_scaling = 0.02
bar_drop_prop = 2
seg_bar_height = 0
contig_bar_height = 3
segment_bar_height = 2
gene_bar_height = 1
ref_bar_height = 0
bar_width = 1
gene_to_locations = defaultdict(list)
overlap_genes = []
all_relGenes = []
prev_start = 0
alternate = True
alternated = False


def plot_bpg_connection(ref_placements, prev_seg_index_is_adj, bpg_dict, seg_end_pos_d):
    connect_width = bar_width / 2.
    for ind, refObj in ref_placements.items():
        next_ind = (ind + 1) % len(ref_placements)
        next_refObj = ref_placements[next_ind]
        if not prev_seg_index_is_adj[next_ind]:  # or next_ind == 0 to try and close
            bpg_adjacency = vu.pair_is_edge(refObj.id, next_refObj.id, refObj.direction, next_refObj.direction,
                                            bpg_dict, seg_end_pos_d)

            if not bpg_adjacency or ind == len(ref_placements) - 1:
                continue

            bpg_connector_len = next_refObj.abs_start_pos - refObj.abs_end_pos
            # makes the reference genome wedges
            patches.append(mpatches.Rectangle((refObj.abs_end_pos, ref_bar_height + bar_width / 4.), bpg_connector_len,connect_width))
            f_color_v.append('grey')
            e_color_v.append('grey')
            lw_v.append(0.2)


# must return positions of transcribed regions of genes here
def plot_gene_track(currStart, currEnd, relGenes, pTup, total_length, seg_dir):
    global prev_start, alternate, alternated
    overlap_genes.append({})
    for gObj in relGenes:
        # e_posns is a list of tuples of exon (start,end)
        # these can be plotted similarly to how the coding region is marked
        gname, gstart, gend, e_posns = gObj.gname, gObj.gstart, gObj.gend, gObj.eposns
        seg_len = pTup[2] - pTup[1]
        if seg_dir == "+":
            normStart = currStart + max(0, gstart - pTup[1])
            normEnd = currStart + min(seg_len, gend - pTup[1])
        else:
            normEnd = currStart + min(seg_len, pTup[2] - gstart)
            normStart = currStart + max(0, pTup[2] - gend)

        gene_to_locations[gname].append((normStart, normEnd))
        box_len = normEnd - normStart
        # patches.append(mpatches.Wedge((0,0), seg_bar_height, start_angle, end_angle, width=bar_width/2.0))
        print(gname, "GBH",gene_bar_height,"BW",bar_width)
        patches.append(mpatches.Rectangle((normStart, gene_bar_height + bar_width), box_len, 0.6 * bar_width))
        f_color_v.append('k')
        e_color_v.append('k')
        lw_v.append(0)
        print(total_length, bar_width)
        # TODO:
        # draw some arrows over the black box
        # but first draw a white line in the box
        # then put some white arrow markers on it

        if gname not in overlap_genes[len(overlap_genes) - 2] or gstart > overlap_genes[len(overlap_genes) - 2].get(
                gname) or args.print_dup_genes:
            if abs(gstart - prev_start) < 8*bar_width and alternate and not alternated:
                ax.text(normStart + box_len / 2., gene_bar_height - 1.5*bar_width + 0.1 * bar_width, gname, style='italic', color='k',
                        ha="center", fontsize=11)
                alternate = False
            elif abs(gstart - prev_start) < 8*bar_width and not alternate and not alternated:
                ax.text(normStart + box_len / 2., gene_bar_height + 0.1 * bar_width, gname, style='italic',
                       color='k',
                       ha="center", fontsize=11)
                alternate = True
                alternated = True
            else:
                ax.text(normStart + box_len / 2., gene_bar_height - 0.7*bar_width + 0.1 * bar_width, gname, style='italic', color='k',
                        ha="center", fontsize=11)
                alternated = False
            prev_start = gstart

        if currEnd < gend:
            overlap_genes[len(overlap_genes) - 1][gname] = gend

        # TODO: add exon plotting
        # for exon in e_posns:
        #     if exon[1] > pTup[1] and exon[0] < pTup[2]:
        #         if strand == "+":
        #             normStart = currStart + max(1,exon[0]-pTup[1])
        #             normEnd = currStart + min(pTup[2]-pTup[1],exon[1]-pTup[1])

        #         else:
        #             normEnd = currStart + min(pTup[2]-pTup[1],pTup[2]-exon[0])
        #             normStart = currStart + max(1,pTup[2] - exon[1])

        #         # start_angle, end_angle = start_end_angle(normStart,normEnd,total_length)
        #         # patches.append(mpatches.Wedge((0,0), seg_bar_height-bar_width/2.0, start_angle, end_angle, width=bar_width/2.0))
        #         patches.append(mpatches.Rectangle((normStart,gene_bar_height + 0.4*bar_width),box_len,0.6*bar_width))
        #         f_color_v.append('r')
        #         e_color_v.append('r')
        #         lw_v.append(0)


# plot the reference genome
def plot_ref_genome(ref_placements, path, total_length, segSeqD, imputed_status, label_segs, color_map,  onco_set=None):
    if onco_set is None:
        onco_set = set()
    font0 = FontProperties()
    for ind, refObj in ref_placements.items():
        print(ind,refObj.to_string(),ref_bar_height)
        seg_coord_tup = segSeqD[path[ind][0]]
        # print(refObj.to_string())
        box_len = refObj.abs_end_pos - refObj.abs_start_pos
        # print start_angle,end_angle

        # makes the reference genome wedges
        patches.append(mpatches.Rectangle((refObj.abs_start_pos, ref_bar_height), box_len, bar_width))
        chrom = segSeqD[path[ind][0]][0]
        if color_map == "standard":
            try:
                c_col = chromosome_colors[chrom]
            except KeyError:
                print("Color not found for " + chrom + ". Using red.")
                chromosome_colors[chrom] = "red"
                c_col = chromosome_colors[chrom]

        else:
            c_col = color_map(float(refObj.id)/len(segSeqD))

        f_color_v.append(c_col)
        e_color_v.append(c_col)
        lw_v.append(0.2)

        # makes the ticks on the reference genome wedges
        if path[ind][1] == "+":
            ts = (seg_coord_tup[1], refObj.abs_start_pos, 0)
            te = (seg_coord_tup[2] + 1, refObj.abs_end_pos+1, 1)
            s = 1

        else:
            ts = (seg_coord_tup[2], refObj.abs_start_pos, 0)
            te = (seg_coord_tup[1] - 1, refObj.abs_end_pos + 1, -1)
            s = -1

        tick_freq = max(20000, 40000 * int(np.floor(total_length / 1000000)))
        print("Tick frequency set to " + str(tick_freq) + "bp")
        posns = []
        a = ts[0]
        b = te[0]
        # print(a, b, ts[0], te[0], step)
        for j in np.arange(a, b, s):
            if j % tick_freq == 0:
                rpos = vu.convert_gpos_to_ropos(j, refObj.abs_start_pos, refObj.abs_end_pos, seg_coord_tup[1],
                                                path[ind][1])
                posns.append((j, rpos, 0))

        for j in posns:
            if j[0] % tick_freq == 0:
                x_i, y_i = j[1], ref_bar_height
                x_f, y_f = j[1], ref_bar_height - bar_width * 0.3
                ax.plot([x_i, x_f], [y_i, y_f], color='grey', linewidth=1)
                txt = " " + str(int(round((j[0]) / 10000)))  # if ha == "left" else str(int(round((j[0])/10000))) + " "
                # txt = str(j[0])
                x_t, y_t = j[1], ref_bar_height - bar_width * 0.4
                ax.text(x_t, y_t, txt, color='grey', rotation=-90, rotation_mode="anchor",
                        ha="left", va="center", fontsize=7)

        # p_end = refObj.abs_end_pos
        # gene_tree = vu.parse_genes(seg_coord_tup[0], args.ref)
        relGenes = vu.rel_genes(gene_tree, seg_coord_tup, copy.copy(onco_set))

        # plot the gene track
        plot_gene_track(refObj.abs_start_pos, refObj.abs_end_pos, relGenes, seg_coord_tup, total_length, path[ind][1])

        # label the segments by number in path
        mid_sp = (refObj.abs_end_pos + refObj.abs_start_pos) / 2
        # text_angle = mid_sp/total_length*360.
        font = font0.copy()
        if imputed_status[ind]:
            font.set_style('italic')
        #     # font.set_weight('bold')

        # text_angle,ha = vu.correct_text_angle(text_angle)

        if label_segs:
            #     ax.text(x,y,cycle[ind][0]+cycle[ind][1],color='grey',rotation=text_angle,
            #         ha=ha,fontsize=5,fontproperties=font,rotation_mode='anchor')
            label_text = path[ind][1]
            if label_segs == "id":
                label_text = path[ind][0] + label_text

            ax.text(mid_sp, ref_bar_height + 0.25 * bar_width, label_text, color='grey', fontsize=8, fontproperties=font, ha='center')


# plot cmap track
def plot_cmap_track(seg_placements, total_length, unadj_bar_height, color, seg_id_labels=False):
    path_label_locs = defaultdict(list)
    for ind, segObj in seg_placements.items():
        bar_height = unadj_bar_height + segObj.track_height_shift
        print("cmap_plot", segObj.id)
        # print "cmap plotting abs end pos are"
        # print segObj.abs_start_pos, segObj.abs_end_pos
        box_len = segObj.abs_end_pos - segObj.abs_start_pos

        # Draw the box
        patches.append(mpatches.Rectangle((segObj.abs_start_pos, bar_height), box_len, bar_width))
        f_color_v.append(color)
        e_color_v.append('k')
        lw_v.append(0)

        linewidth = min(0.5 * 2000000 / total_length, 0.5)
        # Draw the labels in the box
        for i in segObj.label_posns:
            if i > segObj.abs_end_pos or i < segObj.abs_start_pos:
                continue

            y_i, y_f = bar_height, bar_height + bar_width
            ax.plot([i, i], [y_i, y_f], color='k', alpha=0.9, linewidth=linewidth)

        # TODO: fix for dense packing
        if seg_id_labels:
            mid_sp = (segObj.abs_end_pos + segObj.abs_start_pos) / 2
            text = segObj.id + segObj.direction
            ax.text(mid_sp, bar_height + 1.1 * bar_width, text, color='grey', fontsize=9, ha="center")

    return path_label_locs


# plot the connecting lines for the bionano track
def plot_alignment(contig_locs, segment_locs, total_length):
    linewidth = min(0.5 * 2000000 / total_length, 0.5)
    print("linewidth", linewidth, total_length)
    for a_d in aln_vect:
        c_id = a_d["contig_id"]
        c_num_dir = int(a_d["contig_dir"] + "1")
        contig_label_vect = contig_locs[c_id].label_posns
        seg_label_vect = segment_locs[a_d["seg_aln_number"]].label_posns
        clx = contig_label_vect[a_d["contig_label"] - 1]
        slx = seg_label_vect[a_d["seg_label"] - 1]
        # contig_top = seg_bar_height + contig_bar_height + contig_locs[c_id].track_height_shift + bar_width
        contig_bottom = seg_bar_height + contig_bar_height + contig_locs[c_id].track_height_shift
        ax.plot([slx, clx], [seg_bar_height + bar_width, contig_bottom], color="grey", linewidth=linewidth)


def construct_path_ref_placements(path, segSeqD, raw_path_length, prev_seg_index_is_adj, next_seg_index_is_adj,
                                  cycle_seg_counts, aln_vect=None):
    if aln_vect is None:
        aln_vect = []

    spacing_bp = seg_spacing * raw_path_length
    path_ref_placements = {}
    curr_start = 0.0
    for ind, i in enumerate(path):
        seg_id_count = cycle_seg_counts[i[0]]
        seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        seg_end = curr_start + seg_len
        padj, nadj = prev_seg_index_is_adj[ind], next_seg_index_is_adj[ind]
        curr_obj = vu.CycleVizElemObj(i[0], segSeqD[i[0]][0], segSeqD[i[0]][1], segSeqD[i[0]][2], i[1], curr_start,
                                      seg_end, seg_id_count, padj, nadj)

        path_ref_placements[ind] = curr_obj
        next_start = seg_end
        mod_ind = (ind + 1) % (len(prev_seg_index_is_adj))
        if not prev_seg_index_is_adj[mod_ind]:
            next_start += spacing_bp

        curr_start = next_start

    total_length = next_start
    return path_ref_placements, total_length


parser = argparse.ArgumentParser(description="Linear visualizations of AA & AR output")
parser.add_argument("--om_alignments",
                    help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
                    action='store_true')
parser.add_argument("-s", "--om_segs", help="segments cmap file")
parser.add_argument("-g", "--graph", help="breakpoint graph file")
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("--ref", help="reference genome", choices=["hg19", "hg38", "GRCh37", "GRCh38"], default="hg19")
parser.add_argument("--cycles_file", help="AA/AR cycles-formatted input file", required=True)
parser.add_argument("--path", help="path number to visualize", type=int, required=True)
parser.add_argument("--AR_path_alignment", help="AR path alignment file")
parser.add_argument("--outname", help="output prefix")
parser.add_argument("--label_segs", help="label segs with graph IDs and direction (id) or just direction (dir). Default: no labeling", choices=["id", "dir"], default=None)
parser.add_argument("--reduce_path", help="Number of path elements to remove from left and right ends. Must supply both values, \
                    default 0 0", nargs=2, type=int, default=[0, 0])
parser.add_argument("--print_dup_genes", help="If a gene appears multiple times print name every time.",
                    action='store_true', default=False)
parser.add_argument("--color_map", help="Set a matplotlib named color pallete (e.g. Blues, Purples), default is a custom chromosomal map", default="standard")
parser.add_argument("-v", "--version", action='version', version='LinearViz {version} \n Author: Jens Luebeck '
                    '(jluebeck [at] ucsd.edu)'.format(version=__version__))
group2 = parser.add_mutually_exclusive_group(required=False)
group2.add_argument("--gene_subset_file", help="File containing subset of genes to plot (e.g. oncogene genelist file)",
                    default="")
group2.add_argument("--gene_subset_list", help="List of genes to plot (e.g. MYC PVT1)", nargs="+", type=str)


# ----------------------
# handle arguments

args = parser.parse_args()
if args.ref == "GRCh38":
    args.ref = "hg38"

if not args.outname:
    args.outname = os.path.split(args.cycles_file)[1].split(".")[0]

outdir = os.path.dirname(args.outname)
if outdir and not os.path.exists(outdir):
    os.makedirs(outdir)

fname = args.outname + "_path_" + str(args.path) + "_trim_" + str(args.reduce_path[0]) + "_" + str(args.reduce_path[1])

print(args.reduce_path, "path reduction (L, R)")

print("Reading genes")
gene_tree = vu.parse_genes(args.ref, [])

print("Unaligned fraction cutoff set to " + str(vu.unaligned_cutoff_frac))

chromosome_colors = vu.get_chr_colors()
plt.clf()
fig, ax = plt.subplots(figsize=(10, 6))
patches = []
f_color_v = []
e_color_v = []
lw_v = []

# convert the cycles file and reset the arg cycles_file arg
bpg_cf = os.path.basename(args.cycles_file).rsplit("_cycles.txt")[0] + "_BPG_converted_cycles.txt"
print("Converting cycles file segment boundaries to graph file segment boundaries")
make_new_cycle(args.graph, args.cycles_file, bpg_cf)
print(bpg_cf)
args.cycles_file = bpg_cf
paths, segSeqD, circular_D = vu.parse_cycles_file(args.cycles_file)
path_num = args.path
path = paths[path_num]

if args.reduce_path != [0, 0]:
    isCycle = False
else:
    isCycle = circular_D[path_num]

prev_seg_index_is_adj = vu.adjacent_segs(path, segSeqD, isCycle)
print("PSIIA", prev_seg_index_is_adj)
raw_path_length = vu.get_raw_path_length(path, segSeqD)

bpg_dict, seg_end_pos_d = {}, {}
if args.graph:
    bpg_dict, seg_end_pos_d = vu.parse_BPG(args.graph)

gene_set = set()

if args.gene_subset_file.upper() == "BUSHMAN":
    sourceDir = os.path.dirname(os.path.abspath(__file__))
    args.gene_subset_file = sourceDir + "/resources/Bushman_group_allOnco_May2018.tsv"

if args.gene_subset_file:
    gff = True if args.gene_subset_file.endswith(".gff") else False
    gene_set = vu.parse_gene_subset_file(args.gene_subset_file, gff)

elif args.gene_subset_list:
    gene_set = set(args.gene_subset_list)

if not args.color_map == "standard":
    args.color_map = cm.get_cmap(args.color_map)

# ----------------------

if not args.om_alignments:
    if args.reduce_path != [0, 0]:
        # reduce alignments
        path, prev_seg_index_is_adj, _ = vu.reduce_path(path, prev_seg_index_is_adj, args.reduce_path)

    prev_seg_index_is_adj, next_seg_index_is_adj = vu.adjacent_segs(path, segSeqD, isCycle)
    cycle_seg_counts = vu.get_seg_amplicon_count(path)
    ref_placements, total_length = construct_path_ref_placements(path, segSeqD, raw_path_length, prev_seg_index_is_adj,
                                                                 next_seg_index_is_adj, cycle_seg_counts)

    imputed_status = [False] * len(path)
    # set heights
    # order of tracks goes:
    # contig_bar_height (level 3) #unused in this case
    # seg_bar_height (level 2) #unused in this case
    # gene_bar_height (level 1)
    # ref_bar_height (level 0)

    # scales the height based on the length. (bar_width = total_length * some proportion)
    # level is turned into absolute coordinates
    bar_width = total_length * bar_width_scaling

    #put the om contigs on top of the reference om segments
    contig_bar_height += bar_width * bar_drop_prop
    gene_bar_height = seg_bar_height - bar_width * bar_drop_prop + 0.7*bar_width
    ref_bar_height = seg_bar_height - (bar_width * 1.5 * bar_drop_prop) - 0.7*bar_width
    # the following is a holder point to make the plot height work when no OM data is present
    ax.plot(0,seg_bar_height + contig_bar_height, color='white', markersize=10)

else:
    seg_cmaps = parse_cmap(args.om_segs, True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    seg_cmap_lens = get_cmap_lens(args.om_segs)
    aln_vect, meta_dict = vu.parse_alnfile(args.AR_path_alignment)
    if args.reduce_path != [0, 0]:
        # reduce alignments
        path, prev_seg_index_is_adj, aln_vect = vu.reduce_path(path, prev_seg_index_is_adj, args.reduce_path, aln_vect)

    is_segdup, split_ind = vu.check_segdup(aln_vect, path, isCycle)
    if is_segdup:
        print("alignment shows simple segdup")
        path = [path[0]] * 2
        print(path)
        isCycle = False
        prev_seg_index_is_adj = [False, True]
        for a_ind in range(split_ind, len(aln_vect)):
            aln_vect[a_ind]["seg_aln_number"] = 1

    prev_seg_index_is_adj, next_seg_index_is_adj = vu.adjacent_segs(path, segSeqD, isCycle)
    cycle_seg_counts = vu.get_seg_amplicon_count(path)
    ref_placements, total_length = construct_path_ref_placements(path, segSeqD, raw_path_length, prev_seg_index_is_adj,
                                                                 next_seg_index_is_adj, cycle_seg_counts)

    path_seg_placements = vu.place_path_segs_and_labels(path, ref_placements, seg_cmap_vects)

    # set heights
    # this is same as non-om version, but total length is different, thus bar_width is different
    # order of tracks goes:
    # contig_bar_height (level 3)
    # seg_bar_height (level 2)
    # gene_bar_height (level 1)
    # ref_bar_height (level 0)

    # scales the height based on the length. (bar_width = total_length * some proportion)
    # level is turned into absolute coordinates
    bar_width = total_length * bar_width_scaling

    #put the om contigs on top of the reference om segments
    contig_bar_height += bar_width * bar_drop_prop
    gene_bar_height = seg_bar_height - bar_width * bar_drop_prop + 0.7*bar_width
    ref_bar_height = seg_bar_height - (bar_width * 1.5 * bar_drop_prop) - 0.7*bar_width

    contig_cmaps = parse_cmap(args.contigs, True)
    contig_cmap_vects = vectorize_cmaps(contig_cmaps)

    ###
    # TODO: TRIM REF SEGS
    ###

    contig_cmap_lens = get_cmap_lens(args.contigs)
    # path_seg_placements,aln_vect,total_length,contig_cmap_vects
    contig_placements, contig_list = vu.place_contigs_and_labels(path_seg_placements, aln_vect, total_length,
                                                                 contig_cmap_vects, isCycle, True, segSeqD)
    vu.decide_trim_contigs(contig_cmap_vects, contig_placements, total_length)

    # plot segs cmap
    print("SH", seg_bar_height + segment_bar_height)
    print("CH", seg_bar_height + contig_bar_height)
    plot_cmap_track(path_seg_placements, total_length, seg_bar_height + segment_bar_height, "darkorange")

    # check overlaps of contigs and adjust heights accordingly
    contig_height_shifts = vu.set_contig_height_shifts(contig_placements, contig_list, -bar_width)
    # plot contigs cmap
    plot_cmap_track(contig_placements, total_length, seg_bar_height + contig_bar_height, "cornflowerblue",
                    seg_id_labels=True)

    # plot alignments
    plot_alignment(contig_placements, path_seg_placements, total_length)

    imputed_status = vu.imputed_status_from_aln(aln_vect, len(path))

print("RH", ref_bar_height, "BW", bar_width)
plot_ref_genome(ref_placements, path, total_length, segSeqD, imputed_status, args.label_segs, args.color_map, gene_set)

if args.graph:
    plot_bpg_connection(ref_placements, prev_seg_index_is_adj, bpg_dict, seg_end_pos_d)

#ax.set_xlim(-(seg_bar_height+1.25), (seg_bar_height+1.25))
#ax.set_ylim(-(seg_bar_height+1.25), (seg_bar_height+1.25))
chrom_set = set()
for i in path:
    chrom_set.add(segSeqD[i[0]][0])

sorted_chrom = sorted(chrom_set, key=lambda x: x.rsplit("chr")[-1])
sorted_chrom_colors = [chromosome_colors[x] for x in sorted_chrom]
legend_patches = []
for chrom, color in zip(sorted_chrom, sorted_chrom_colors):
    legend_patches.append(mpatches.Patch(color=color, label=chrom))

# plt.legend(handles=legend_patches,fontsize=8,loc=3,bbox_to_anchor=(-.3,.15))
plt.legend(handles=legend_patches, fontsize=10,
           bbox_to_anchor=(0, 0))  # bbox_to_anchor=(0,-1.5))#,bbox_to_anchor=(.09,-1.5))

p = PatchCollection(patches)
p.set_facecolor(f_color_v)
p.set_edgecolor(e_color_v)
p.set_linewidth(lw_v)
ax.add_collection(p)
ax.set_aspect(1.0)
plt.axis('off')
plt.savefig(fname + '.png', dpi=600)
plt.savefig(fname + '.pdf', format='pdf')

plt.close()
print("finished")