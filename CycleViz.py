#!/usr/bin/env python

import argparse
import copy
import os

import matplotlib
matplotlib.use('Agg')  # this import must happen immediately after importing matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches
import numpy as np

from bionanoUtil import *
import VizUtil as vu

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

seg_spacing = 0.009
bar_width = 2.5 / 3
global_rot = 90.0
center_hole = 0.75
outer_bar = 10
bed_spacing = .5
contig_bar_height = -14 / 3
segment_bar_height = -8.0 / 3
gene_to_locations = defaultdict(list)
overlap_genes = []
all_relGenes = []


# get the start and end angle from the linear start and end
def start_end_angle(normStart, normEnd, total_length):
    start_angle = normStart / total_length * 360
    end_angle = normEnd / total_length * 360
    # does it cross the start boundary?
    if end_angle < 0 and start_angle > 0:
        end_angle += 360

    # handle circular contig
    if start_angle > 360 and start_angle % 360 > end_angle:
        start_angle, end_angle = 360, 0

    return start_angle, end_angle


def plot_bpg_connection(ref_placements, cycle, total_length, prev_seg_index_is_adj, bpg_dict, seg_end_pos_d):
    connect_width = bar_width / 2.
    for ind, refObj in ref_placements.items():
        next_ind = (ind + 1) % len(ref_placements)
        next_refObj = ref_placements[next_ind]
        if not prev_seg_index_is_adj[next_ind]:  # or next_ind == 0 to try and close
            bpg_adjacency = vu.pair_is_edge(refObj.id, next_refObj.id, refObj.direction, next_refObj.direction,
                                            bpg_dict, seg_end_pos_d)

            if not bpg_adjacency:
                continue

            start_angle, end_angle = start_end_angle(next_refObj.abs_start_pos, refObj.abs_end_pos, total_length)
            # makes the reference genome wedges
            patches.append(
                mpatches.Wedge((0, 0), outer_bar - bar_width / 4, end_angle, start_angle, width=connect_width))
            f_color_v.append('grey')
            e_color_v.append('grey')
            lw_v.append(0.2)


def plot_feature_track(currStart, currEnd, seg_dir, pTup, cfc, curr_chrom, total_length, seg_copies):
    gc, gs, ge = pTup
    granularity = cfc.track_props['granularity']
    if granularity == 0:
        granularity = max(1,int((ge - gs)/10000.0))

    # print(cfc.track_max, cfc.track_min, cfc.top, cfc.base)
    height_scale_factor = (cfc.top - cfc.base)/float(cfc.track_max - cfc.track_min)

    # plot the legens
    legend_points = np.linspace(currStart / total_length * 2 * np.pi, (currEnd + 1) / total_length * 2 * np.pi, 10000)
    lheights = list(np.linspace(cfc.base, cfc.top, 5))
    lheights.append(cfc.top)

    lvals = list(np.linspace(cfc.track_min, cfc.track_max, 5))
    lvals.append(cfc.track_max)

    print("TRACK LEGEND HEIGHTS", lvals)
    for lh in lheights:
        x_v, y_v = vu.polar_series_to_cartesians(legend_points, lh)
        plt.plot(x_v, y_v, color='lightgrey', linewidth=0.25, zorder=0)

    tertiary_data = []
    tertiary_style = 'lines'
    tertiary_color = 'mediumorchid'
    if cfc.track_props['show_segment_copy_count']:
        v = 2*seg_copies*cfc.track_props['segment_copy_count_scaling']
        tertiary_data = [[gs, ge, v]]


    hs = cfc.track_props['hide_secondary']
    if cfc.track_props['hide_secondary'] == "viral" and not (curr_chrom.startswith('chr') or len(curr_chrom) < 3):
        hs = True

    elif cfc.track_props['hide_secondary'] == "viral":
        hs = False

    for data_it, style, curr_color, elem_ind in zip([cfc.primary_data[curr_chrom], cfc.secondary_data[curr_chrom],tertiary_data],
                                          [cfc.track_props['primary_style'], cfc.track_props['secondary_style'], tertiary_style],
                                          [cfc.track_props['primary_color'], cfc.track_props['secondary_color'], tertiary_color],
                                          list(range(3))):

        datalist = [(max(x[0],gs), min(x[1], ge), height_scale_factor*x[2] + cfc.base) for x in data_it] #restrict to the coordinates of the region

        # convert the data into granular form
        sortrevdir = False if seg_dir == "+" else True
        datalist.sort(reverse=sortrevdir)
        point_data = []
        val_data = []
        for p in datalist:
            if p[1] - p[0] > granularity:
                plocs = np.linspace(p[0], p[1], (p[1] - p[0])/granularity)
                for newp in plocs:
                    point_data.append(newp)
                    val_data.append(p[2])

                point_data.append(p[1])
                val_data.append(p[2])

            else:
                point_data.append((p[0] + p[1])/2.0)
                val_data.append(p[2])

        # set the direction and convert to polars from proportional length
        if seg_dir == "+":
            normed_data = [(currStart + x - gs)/total_length * 2 * np.pi for x in point_data]
        else:
            normed_data = [(currStart + ge - x)/total_length * 2 * np.pi for x in point_data]

        # convert to cartesians
        x_v, y_v = vu.polar_series_to_cartesians(normed_data, val_data)

        zorder = 0 if elem_ind == 1 else 1

        # draw the points/lines
        if style == "points":
            if elem_ind != 1 or not hs:
                plt.scatter(x_v, y_v, s=cfc.track_props['pointsize'], color=curr_color, zorder=zorder)

        elif style == "lines":
            if elem_ind != 1 or not hs:
                plt.plot(x_v, y_v, linewidth=cfc.track_props['linewidth'], color=curr_color, zorder=zorder)

            # if seg_dir == "+":
            #     normeddata = [(currStart + x[0] - gs, currStart + x[1] - gs, x[2]) for x in datalist]
            # else:
            #     normeddata = [(currStart + ge - x[1], currStart + ge - x[0], x[2]) for x in datalist]

        else:
            print("feature_style must be either 'points' or 'lines'\n")


def plot_gene_direction_indicator(normStart, normEnd, drop):
    start_angle = normStart / total_length * 360
    end_angle = normEnd / total_length * 360

    patches.append(mpatches.Wedge((0, 0), outer_bar - drop, start_angle, end_angle, width=bar_width / 4.5))
    f_color_v.append('k')
    e_color_v.append('k')
    lw_v.append(0)


# must return positions of transcribed regions of genes here
def plot_gene_track(currStart, currEnd, relGenes, pTup, total_length, seg_dir, ind, plot_gene_direction=True):
    overlap_genes.append({})
    for gObj in relGenes:
        # e_posns is a list of tuples of exon (start,end)
        # these can be plotted similarly to how the coding region is marked
        gname, gstart, gend, e_posns = gObj.gname, gObj.gstart, gObj.gend, gObj.eposns
        # print(gname, gstart, gend, pTup, len(gObj.gdrops))
        seg_len = pTup[2] - pTup[1]
        hasStart = False
        hasEnd = False
        if seg_dir == "+":
            ts = max(0, gstart - pTup[1])
            te = min(seg_len, gend - pTup[1])
            if gObj.strand == "+":
                drop = 1.4 * bar_width
                if ts > 0: hasStart = True
                if te < seg_len: hasEnd = True
            else:
                drop = 2 * bar_width
                if ts > 0: hasEnd = True
                if te < seg_len: hasStart = True

            normStart = currStart + max(0, gstart - pTup[1])
            normEnd = currStart + min(seg_len, gend - pTup[1])

        else:
            te = min(seg_len, pTup[2] - gstart)
            ts = max(0, pTup[2] - gend)
            if gObj.strand == "+":
                drop = 2 * bar_width
                if te < seg_len: hasStart = True
                if ts > 0: hasEnd = True
            else:
                drop = 1.4 * bar_width
                if te < seg_len: hasEnd = True
                if ts > 0: hasStart = True

            normEnd = currStart + min(seg_len, pTup[2] - gstart)
            normStart = currStart + max(0, pTup[2] - gend)

        start_angle = normStart / total_length * 360
        end_angle = normEnd / total_length * 360
        text_angle = (start_angle + end_angle) / 2.0
        gene_to_locations[gname].append((start_angle / 360., end_angle / 360.))
        if end_angle < 0 and start_angle > 0:
            end_angle += 360

        patches.append(mpatches.Wedge((0, 0), outer_bar, start_angle, end_angle, width=bar_width / 2.0))
        f_color_v.append('k')
        e_color_v.append('k')
        lw_v.append(0)

        if gname not in overlap_genes[len(overlap_genes)-2] or not overlap_genes[len(overlap_genes)-2].get(gname)[0] or seg_dir != overlap_genes[len(overlap_genes)-2].get(gname)[1]:
            x_t, y_t = vu.pol2cart(outer_bar + bar_width + 1.7, (text_angle / 360 * 2 * np.pi))
            text_angle, ha = vu.correct_text_angle(text_angle)

            if gObj.highlight_name:
                ax.text(x_t, y_t, gname, style='italic', color='r', rotation=text_angle, ha=ha, va="center",
                        fontsize=gene_fontsize, rotation_mode='anchor')
            else:
                ax.text(x_t, y_t, gname, style='italic', color='k', rotation=text_angle, ha=ha, va="center",
                        fontsize=gene_fontsize, rotation_mode='anchor')

        # draw something to show direction and truncation status
        if plot_gene_direction:
            plot_gene_direction_indicator(normStart, normEnd, drop)
            gObj.gdrops.append((normStart, normEnd, total_length, seg_dir, currStart, currEnd, hasStart, hasEnd, ind, drop, pTup))
            # gObj.gdrops = [(normStart, normEnd, total_length, seg_dir, currStart, currEnd, pTup), ]

        print("PTUPCHECK",pTup,gstart,gend)
        if not (pTup[2] >= gend and pTup[1] <= gstart):
            overlap_genes[len(overlap_genes)-1][gname] = (True, seg_dir)

        for exon in e_posns:
            #fix exon orientation
            if exon[1] > pTup[1] and exon[0] < pTup[2]:
                if seg_dir == "+":
                    normStart = currStart + max(1, exon[0] - pTup[1])
                    normEnd = currStart + min(pTup[2] - pTup[1], exon[1] - pTup[1])

                else:
                    normEnd = currStart + min(pTup[2] - pTup[1], pTup[2] - exon[0])
                    normStart = currStart + max(1, pTup[2] - exon[1])

                start_angle, end_angle = start_end_angle(normStart, normEnd, total_length)
                patches.append(
                    mpatches.Wedge((0, 0), outer_bar - bar_width / 2.0, start_angle, end_angle, width=bar_width / 2.0))
                f_color_v.append('k')
                e_color_v.append('k')
                lw_v.append(0)


# plot the reference genome
def plot_ref_genome(ref_placements, cycle, total_length, segSeqD, imputed_status, label_segs, edge_ticks, onco_set=None):
    if onco_set is None:
        onco_set = set()

    font0 = FontProperties()
    p_end = 0
    # rot_sp = global_rot / 360. * total_length
    for ind, refObj in ref_placements.items():
        seg_coord_tup = segSeqD[cycle[ind][0]]
        start_angle, end_angle = start_end_angle(refObj.abs_end_pos, refObj.abs_start_pos, total_length)

        # makes the reference genome wedges
        patches.append(mpatches.Wedge((0, 0), outer_bar, end_angle, start_angle, width=bar_width))
        chrom = segSeqD[cycle[ind][0]][0]
        try:
            f_color_v.append(chromosome_colors[chrom])
        except KeyError:
            print("Color not found for " + chrom + ". Using red.")
            chromosome_colors[chrom] = "red"
            f_color_v.append("red")

        e_color_v.append(chromosome_colors[chrom])
        lw_v.append(0.2)

        # makes the ticks on the reference genome wedges
        if cycle[ind][1] == "+":
            posns = zip(range(seg_coord_tup[1], seg_coord_tup[2] + 1),
                        np.arange(refObj.abs_start_pos, refObj.abs_end_pos+1))
        else:
            posns = zip(np.arange(seg_coord_tup[2], seg_coord_tup[1] - 1, -1),
                        np.arange(refObj.abs_start_pos, refObj.abs_end_pos+1))

        tick_freq = max(10000, 30000 * int(np.floor(total_length / 800000)))
        # if refObj.abs_end_pos - refObj.abs_start_pos < 30000:
        # tick_freq = 25000

        # put the positions on the ends of the joined segs
        if edge_ticks:
            newposns = []
            tick_freq = 1
            text_trunc = 1
            # print(posns,seg_coord_tup,refObj.abs_start_pos, refObj.abs_end_pos)
            if not refObj.prev_is_adjacent:
                newposns.append(posns[0])

            if not refObj.next_is_adjacent:
                newposns.append(posns[-1])

            posns = newposns


        # put the positions not on the ends
        else:
            # if there are no labels present on the segment given the current frequency, AND this refobject is not adjacent
            # to the previous, get positions in this segment divisible by 10kbp, set the middle one as the labelling site
            # else just set it to 10000
            text_trunc = 10000

            if (not any(j[0] % tick_freq == 0 for j in posns)) and abs(refObj.abs_start_pos - p_end) > 1:
                tens = [j[0] for j in posns if j[0] % 10000 == 0]
                middleIndex = int((len(tens) - 1) / 2)
                if tens:
                    tick_freq = tens[middleIndex]
                else:
                    tick_freq = 10000


        p_end = refObj.abs_end_pos
        print("tick freq", tick_freq)
        for j in posns:
            if j[0] % tick_freq == 0:
                text_angle = j[1] / total_length * 360
                x, y = vu.pol2cart(outer_bar, (text_angle / 360 * 2 * np.pi))
                x_t, y_t = vu.pol2cart(outer_bar + 0.2, (text_angle / 360 * 2 * np.pi))
                ax.plot([x, x_t], [y, y_t], color='grey', linewidth=1)

                text_angle, ha = vu.correct_text_angle(text_angle)
                txt = " " + str(int(round((j[0]) / text_trunc))) if ha == "left" else str(int(round((j[0]) / text_trunc))) + " "

                ax.text(x_t, y_t, txt, color='grey', rotation=text_angle,
                        ha=ha, va="center", fontsize=tick_fontsize, rotation_mode='anchor')

        # gene_tree = vu.parse_genes(seg_coord_tup[0], args.ref)
        relGenes = vu.rel_genes(gene_tree, seg_coord_tup, copy.copy(onco_set))
        all_relGenes.extend(relGenes)
        # plot the gene track
        # print(ind, refObj.to_string(), len(relGenes))
        plot_gene_track(refObj.abs_start_pos, refObj.abs_end_pos, relGenes, seg_coord_tup, total_length, cycle[ind][1], ind)
        for index, cfc in enumerate(refObj.feature_tracks):
            plot_feature_track(refObj.abs_start_pos, refObj.abs_end_pos, refObj.direction, seg_coord_tup, cfc,
                               refObj.chrom, total_length, refObj.seg_count)

        # label the segments by number in cycle
        mid_sp = (refObj.abs_end_pos + refObj.abs_start_pos) / 2
        text_angle = mid_sp / total_length * 360.
        x, y = vu.pol2cart((outer_bar - 2 * bar_width), (text_angle / 360. * 2. * np.pi))
        font = font0.copy()
        if imputed_status[ind]:
            font.set_style('italic')
            # font.set_weight('bold')

        text_angle, ha = vu.correct_text_angle(text_angle)

        if label_segs:
            ax.text(x, y, cycle[ind][0] + cycle[ind][1], color='grey', rotation=text_angle,
                    ha=ha, fontsize=5, fontproperties=font, rotation_mode='anchor')


# set the heights of the bed track features
def get_feature_heights(ntracks, intertrack_spacing):
    if ntracks > 0:
        maxtop = outer_bar-(intertrack_spacing+2)
        bases = np.linspace(center_hole, maxtop, ntracks)
        tops = [x - intertrack_spacing for x in bases[1:]]
        tops.append(maxtop)
        return bases, tops

    return [], []


# plot cmap track for bionano
def plot_cmap_track(seg_placements, total_length, unadj_bar_height, color, seg_id_labels=False):
    cycle_label_locs = defaultdict(list)
    for ind, segObj in seg_placements.items():
        bar_height = unadj_bar_height + segObj.track_height_shift
        print("cmap_plot", segObj.id)
        start_angle, end_angle = start_end_angle(segObj.abs_end_pos, segObj.abs_start_pos, total_length)
        patches.append(mpatches.Wedge((0, 0), bar_height + bar_width, end_angle, start_angle, width=bar_width))
        f_color_v.append(color)
        e_color_v.append('k')
        lw_v.append(0)

        linewidth = min(0.25 * 2000000 / total_length, 0.25)
        for i in segObj.label_posns:
            if i > segObj.abs_end_pos or i < segObj.abs_start_pos:
                continue

            label_rads = i / total_length * 2 * np.pi
            x, y = vu.pol2cart(bar_height, label_rads)
            x_t, y_t = vu.pol2cart(bar_height + bar_width, label_rads)
            # linewidth = min(0.2*2000000/total_length,0.2)
            ax.plot([x, x_t], [y, y_t], color='k', alpha=0.9, linewidth=linewidth)

        if seg_id_labels:
            mid_sp = (segObj.abs_end_pos + segObj.abs_start_pos) / 2
            text_angle = mid_sp / total_length * 360.
            x, y = vu.pol2cart(bar_height - 1.2, (text_angle / 360. * 2. * np.pi))
            text_angle, ha = vu.correct_text_angle(text_angle)
            text = segObj.id + segObj.direction
            ax.text(x, y, text, color='grey', rotation=text_angle,
                    ha=ha, fontsize=5, rotation_mode='anchor')

    return cycle_label_locs


# plot the connecting lines for the bionano track
def plot_alignment(contig_locs, segment_locs, total_length):
    segs_base = outer_bar + segment_bar_height
    linewidth = min(0.25 * 2000000 / total_length, 0.25)
    for a_d in aln_vect:
        c_id = a_d["contig_id"]
        c_num_dir = int(a_d["contig_dir"] + "1")

        contig_label_vect = contig_locs[c_id].label_posns
        seg_label_vect = segment_locs[a_d["seg_aln_number"]].label_posns
        c_l_pos = contig_label_vect[a_d["contig_label"] - 1]
        c_l_loc = c_l_pos / total_length * 2. * np.pi
        # s_l_pos = seg_label_vect[s_num_dir*a_d["seg_label"]-(s_num_dir+1)/2]
        s_l_pos = seg_label_vect[a_d["seg_label"] - 1]
        s_l_loc = s_l_pos / total_length * 2. * np.pi
        contig_top = outer_bar + contig_bar_height + contig_locs[c_id].track_height_shift + bar_width
        x_c, y_c = vu.pol2cart(contig_top, c_l_loc)
        x_s, y_s = vu.pol2cart(segs_base, s_l_loc)
        ax.plot([x_c, x_s], [y_c, y_s], color="grey", linewidth=linewidth)


def construct_cycle_ref_placements(cycle, segSeqD, raw_cycle_length, prev_seg_index_is_adj, next_seg_index_is_adj,
                                   isCycle, cycle_seg_counts, aln_vect=None):
    if aln_vect is None:
        aln_vect = []

    spacing_bp = seg_spacing * raw_cycle_length
    cycle_ref_placements = {}
    curr_start = 0.0 if isCycle else spacing_bp
    for ind, i in enumerate(cycle):
        seg_id_count = cycle_seg_counts[i[0]]
        seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        seg_end = curr_start + seg_len
        padj, nadj = prev_seg_index_is_adj[ind], next_seg_index_is_adj[ind]
        curr_obj = vu.CycleVizElemObj(i[0], segSeqD[i[0]][0], segSeqD[i[0]][1], segSeqD[i[0]][2], i[1], curr_start,
                                      seg_end, seg_id_count, padj, nadj)
        cycle_ref_placements[ind] = curr_obj
        next_start = seg_end
        mod_ind = (ind + 1) % (len(prev_seg_index_is_adj))
        if not prev_seg_index_is_adj[mod_ind]:
            next_start += spacing_bp

        curr_start = next_start

    total_length = next_start
    return cycle_ref_placements, total_length


parser = argparse.ArgumentParser(description="Circular visualizations of AA & AR output")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--input_yaml_file", help="Specifiy all desired arguments in this file, OR use the options below\n")
group.add_argument("--cycles_file", help="AA/AR cycles-formatted input file [required]")
parser.add_argument("--cycle", help="cycle number to visualize [required]")
parser.add_argument("-g", "--graph", help="breakpoint graph file [required]")
parser.add_argument("--ref", help="reference genome", choices=["hg19", "hg38", "GRCh37", "GRCh38"], default="hg19")
parser.add_argument("--om_alignments",
                    help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
                    action='store_true')
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--rot", help="number of segments to rotate counterclockwise", type=int, default=0)
parser.add_argument("--label_segs", help="label segs with graph IDs", action='store_true')
parser.add_argument("--gene_subset_file", help="file containing subset of genes to plot (e.g. oncogene genelist file)",
                    default="")
parser.add_argument("--gene_subset_list", help="list of genes to plot (e.g. MYC PVT1)", nargs="+", type=str)
parser.add_argument("--print_dup_genes", help="if a gene appears multiple times print name every time.",
                    action='store_true', default=False)
parser.add_argument("--gene_highlight_list", help="list of gene names to highlight", nargs="+", type=str, default=[])
parser.add_argument("--gene_fontsize", help="font size for gene names", type=float, default=7)
parser.add_argument("--segment_end_ticks", help="Only place ticks at ends of non-contiguous segments",
                    action='store_true', default=False)
parser.add_argument("--tick_fontsize", help="font size for genomic position ticks", type=float, default=7)
parser.add_argument("--feature_yaml_list", nargs='+', help="list of the input yamls for bedgraph file feature "
                    "specifying additional data. Will be plotted from outside to inside given the order the filenames "
                    "appear in", default = [])

args = parser.parse_args()
if args.input_yaml_file:
    vu.parse_main_args_yaml(args)

if args.ref == "GRCh38":
    args.ref = "hg38"

print(args.ref)

if not args.sname:
    args.sname = os.path.split(args.cycles_file)[1].split(".")[0] + "_"

outdir = os.path.dirname(args.sname)
if outdir and not os.path.exists(outdir):
    os.makedirs(outdir)

fname = args.sname + "cycle_" + args.cycle

print("Reading genes")
gene_tree = vu.parse_genes(args.ref, args.gene_highlight_list)

print("Unaligned fraction cutoff set to " + str(vu.unaligned_cutoff_frac))
chromosome_colors = vu.get_chr_colors()
plt.clf()
fig, ax = plt.subplots()
patches = []
f_color_v = []
e_color_v = []
lw_v = []

cycles, segSeqD, circular_D = vu.parse_cycles_file(args.cycles_file)
cycle_num = args.cycle
isCycle = circular_D[cycle_num]
cycle = cycles[cycle_num]
cycle_seg_counts = vu.get_seg_amplicon_count(cycle)
prev_seg_index_is_adj, next_seg_index_is_adj = vu.adjacent_segs(cycle, segSeqD, isCycle)
raw_cycle_length = vu.get_raw_path_length(cycle, segSeqD)

gene_fontsize = args.gene_fontsize
tick_fontsize = args.tick_fontsize

bpg_dict, seg_end_pos_d = {}, {}
if args.graph:
    bpg_dict, seg_end_pos_d = vu.parse_BPG(args.graph)

gene_set = set()

if args.gene_subset_file.upper() == "BUSHMAN":
    sourceDir = os.path.dirname(os.path.abspath(__file__))
    args.gene_subset_file = sourceDir + "/Bushman_group_allOnco_May2018.tsv"

if args.gene_subset_file and not args.gene_subset_file == "NO":
    gff = True if args.gene_subset_file.endswith(".gff") else False
    gene_set = vu.parse_gene_subset_file(args.gene_subset_file, gff)

elif args.gene_subset_list:
    gene_set = set(args.gene_subset_list)


fbases, ftops = get_feature_heights(len(args.feature_yaml_list), 0.5)

if not args.om_alignments:
    ref_placements, total_length = construct_cycle_ref_placements(cycle, segSeqD, raw_cycle_length,
                                                                  prev_seg_index_is_adj, next_seg_index_is_adj,
                                                                  isCycle, cycle_seg_counts)
    imputed_status = [False] * len(cycle)

    # bedgraph
    if args.feature_yaml_list:
        for ind, yaml_file in enumerate(args.feature_yaml_list):
            cfc = vu.parse_feature_yaml(yaml_file, ind+1, len(args.feature_yaml_list))
            cfc.base, cfc.top = fbases[ind], ftops[ind]
            vu.store_bed_data(cfc, ref_placements, cfc.track_props['end_trim'])
            vu.reset_track_min_max(ref_placements, len(args.feature_yaml_list))

# only if bionano data present
else:
    print("Visualizing with alignments")
    print("Contig spacing set to " + str(vu.contig_spacing))
    seg_cmaps = parse_cmap(args.segs, True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    seg_cmap_lens = get_cmap_lens(args.segs)
    aln_vect, meta_dict = vu.parse_alnfile(args.path_alignment)
    is_segdup, split_ind = vu.check_segdup(aln_vect, cycle, isCycle)
    if is_segdup:
        print("alignment shows simple segdup")
        cycle = [cycle[0]] * 2
        print(cycle)
        isCycle = False
        prev_seg_index_is_adj = [False, True]
        next_seg_index_is_adj = [True, False]
        for a_ind in range(split_ind, len(aln_vect)):
            aln_vect[a_ind]["seg_aln_number"] = 1

    ref_placements, total_length = construct_cycle_ref_placements(cycle, segSeqD, raw_cycle_length,
                                                                  prev_seg_index_is_adj, next_seg_index_is_adj, isCycle,
                                                                  cycle_seg_counts, aln_vect)

    cycle_seg_placements = vu.place_path_segs_and_labels(cycle, ref_placements, seg_cmap_vects)

    # bedgraph
    if args.feature_yaml_list:
        for ind, yaml_file in enumerate(args.feature_yaml_list):
            cfc = vu.parse_feature_yaml(yaml_file, ind + 1, len(args.feature_yaml_list))
            cfc.base, cfc.top = fbases[ind], ftops[ind]
            vu.store_bed_data(cfc, ref_placements, cfc.track_props['end_trim'])
            vu.reset_track_min_max(ref_placements, len(args.feature_yaml_list))

    contig_cmaps = parse_cmap(args.contigs, True)
    contig_cmap_vects = vectorize_cmaps(contig_cmaps)

    ###
    # TODO: TRIM REF SEGS
    ###

    contig_cmap_lens = get_cmap_lens(args.contigs)
    contig_placements, contig_list = vu.place_contigs_and_labels(cycle_seg_placements, aln_vect, total_length,
                                                                 contig_cmap_vects, isCycle, True, segSeqD)
    vu.decide_trim_contigs(contig_cmap_vects, contig_placements, total_length)

    # plot cmap segs
    plot_cmap_track(cycle_seg_placements, total_length, outer_bar + segment_bar_height, "darkorange")

    # check overlaps of contigs and adjust heights accordingly
    contig_height_shifts = vu.set_contig_height_shifts(contig_placements, contig_list)
    # plot contigs
    plot_cmap_track(contig_placements, total_length, outer_bar + contig_bar_height, "cornflowerblue",
                    seg_id_labels=True)

    # plot alignments
    plot_alignment(contig_placements, cycle_seg_placements, total_length)
    imputed_status = vu.imputed_status_from_aln(aln_vect, len(cycle))

plot_ref_genome(ref_placements, cycle, total_length, segSeqD, imputed_status, args.label_segs, args.segment_end_ticks,
                gene_set)

for gObj in all_relGenes:
    gObj.draw_marker_ends(outer_bar)
    gObj.draw_seg_links(outer_bar, bar_width)
    gObj.draw_trunc_spots(outer_bar)

if args.graph:
    plot_bpg_connection(ref_placements, cycle, total_length, prev_seg_index_is_adj, bpg_dict, seg_end_pos_d)

ax.set_xlim(-(outer_bar + 1.25), (outer_bar + 1.25))
ax.set_ylim(-(outer_bar + 3.3), (outer_bar + 3.3))
chrom_set = set()
for i in cycle:
    chrom_set.add(segSeqD[i[0]][0])

sorted_chrom = sorted(chrom_set, key=lambda x: x.rsplit("chr")[-1])
sorted_chrom_colors = [chromosome_colors[x] for x in sorted_chrom]
legend_patches = []
for chrom, color in zip(sorted_chrom, sorted_chrom_colors):
    legend_patches.append(mpatches.Patch(color=color, label=chrom))

plt.legend(handles=legend_patches, fontsize=8, loc=3, bbox_to_anchor=(-.3, .15))

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
