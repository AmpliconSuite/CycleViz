#!/usr/bin/env python

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"
__version__ = "0.1.7"

import argparse
from collections import defaultdict
import copy
import os
import sys

from ast import literal_eval as make_tuple
from intervaltree import IntervalTree
import matplotlib
matplotlib.use('Agg')  # this import must happen immediately after importing matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np

from bionanoUtil import *
import VizUtil as vu

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['pdf.fonttype'] = 42

CV_RESOURCES = os.path.dirname(os.path.abspath(__file__)) + "/resources/"
print("Using resources in " + CV_RESOURCES)

seg_spacing = 0.009
bar_width = 2.1 / 3
global_rot = 90.0
center_hole = 1.25
outer_bar = 10
gene_spacing = 1.

contig_bar_height = -13.5 / 3
segment_bar_height = -8.0 / 3
gene_to_locations = defaultdict(list)
overlap_genes = []
all_relGenes = []


# get the start and end angle from the linear start and end
def start_end_angle(normStart, normEnd, total_length):
    start_angle = normStart / total_length * 360
    end_angle = normEnd / total_length * 360
    # does it cross the start boundary?
    if end_angle < 0 and start_angle > 0 and (start_angle - end_angle) > 360:
        # print(start_angle, end_angle, "X")
        end_angle += 360

    # handle circular contig
    if start_angle > 360 and start_angle % 360 > end_angle:
        start_angle, end_angle = 360, 0

    return start_angle, end_angle


def plot_bpg_connection(ref_placements, total_length, prev_seg_index_is_adj=None, bpg_dict=None, seg_end_pos_d=None,
                        manual_links=None, connect_width="auto"):

    if connect_width == "full":
        connect_width = bar_width
        ch = 0
        prev_seg_index_is_adj = defaultdict(bool)

    elif connect_width == "none":
        return

    else:
        if prev_seg_index_is_adj and bpg_dict and seg_end_pos_d:
            connect_width = bar_width / 2.
            ch = bar_width/4.

        else:
            connect_width = bar_width/15.0
            ch = bar_width/2.
            prev_seg_index_is_adj = defaultdict(bool)


    for ind, refObj in ref_placements.items():
        if refObj.custom_bh:
            curr_bh = refObj.custom_bh + refObj.track_height_shift
        else:
            curr_bh = outer_bar

        # this is for interior segment tracks
        if refObj.custom_face_color:
            connect_col = refObj.custom_face_color

        else:
            connect_col = 'grey'

        next_ind = (ind + 1) % len(ref_placements)
        next_refObj = ref_placements[next_ind]
        if not prev_seg_index_is_adj[next_ind]:  # or next_ind == 0 to try and close
            if not manual_links:
                bpg_adjacency = vu.pair_is_edge(refObj.id, next_refObj.id, refObj.direction, next_refObj.direction,
                                            bpg_dict, seg_end_pos_d)

            else:
                bpg_adjacency = manual_links[ind]

            if not bpg_adjacency:
                # print("No link", refObj.to_string(), next_refObj.to_string())
                continue

            start_angle, end_angle = start_end_angle(next_refObj.abs_start_pos, refObj.abs_end_pos, total_length)
            # makes the reference genome wedges
            ax.add_patch(
                mpatches.Wedge((0, 0), curr_bh - ch, end_angle, start_angle, edgecolor=connect_col,
                               facecolor=connect_col, linewidth=0, width=connect_width))
            # f_color_v.append(connect_col)
            # e_color_v.append(connect_col)
            # lw_v.append(0)


def plot_links(cfc):
    ig = cfc.base + intertrack_spacing  # radial location on the inside of the plot where the link passes over
    og = cfc.top + intertrack_spacing/3  # radial location on the edge of the plot where the link originates
    for currlinks in [cfc.primary_links, cfc.secondary_links]:
        for cLink in currlinks:
            for a_tup in cLink.posA_hits:
                acenter = sum(a_tup) / 2.0
                btuplist = []
                if cfc.track_props['link_single_match']:
                    btupmin = None
                    mindist = float('inf')
                    for b_tup in cLink.posB_hits:
                        bcenter = sum(b_tup)/2.0
                        tupdist1 = abs(bcenter - acenter)
                        tupdist2 = (total_length - bcenter) + acenter
                        tupdist = min(tupdist1, tupdist2)
                        if tupdist < mindist:
                            mindist = tupdist
                            btupmin = b_tup

                    if btupmin:
                        btuplist = [btupmin,]

                else:
                    btuplist = cLink.posB_hits

                for b_tup in btuplist:
                    bcenter = sum(b_tup) / 2.0
                    if a_tup[0] - a_tup[1] != 0 or b_tup[0] - b_tup[1] != 0 and not cfc.track_props['linkpoint'] == "midpoint":
                        aloc_0 = a_tup[0]
                        bloc_0 = b_tup[0]
                        aloc_1 = a_tup[1]
                        bloc_1 = b_tup[1]
                        # a0, am, a1, a1, b1, b1, bm, b0, b0, a0, a0
                        alocs = [aloc_0, acenter, aloc_1, aloc_1, bloc_0, bloc_0, bcenter, bloc_1, bloc_1, aloc_0, aloc_0, aloc_0]
                        aguides = [og, og, og, ig, ig, og, og, og, ig, ig, og, og]
                        aphis = np.multiply(alocs, ((1.0 / total_length) * 2 * np.pi))
                        point_zip = zip(aphis, aguides)
                        codes = [
                            Path.MOVETO,
                            Path.CURVE3,
                            Path.CURVE3,
                            Path.CURVE4,
                            Path.CURVE4,
                            Path.CURVE4,
                            Path.CURVE3,
                            Path.CURVE3,
                            Path.CURVE4,
                            Path.CURVE4,
                            Path.CURVE4,
                            Path.CLOSEPOLY
                        ]
                        fc = cLink.link_color
                        ec = 'lightgrey'

                    else:
                        aphi = acenter / total_length * 2 * np.pi
                        bphi = bcenter / total_length * 2 * np.pi
                        codes = [
                            Path.MOVETO,
                            Path.CURVE4,
                            Path.CURVE4,
                            Path.CURVE4,
                        ]
                        point_zip = zip([aphi, aphi, bphi, bphi], [og, ig, ig, og])
                        fc = 'none'
                        ec = cLink.link_color

                    verts = []
                    for phi, g in point_zip:
                        x, y = vu.pol2cart(g, phi)
                        verts.append((x, y))
                    # x_a, y_a = vu.pol2cart(outer_guide, aphi)
                    # x_b, y_b = vu.pol2cart(outer_guide, bphi)
                    # x_a_i, y_a_i = vu.pol2cart(inner_guide, aphi)
                    # x_b_i, y_b_i = vu.pol2cart(inner_guide, bphi)
                    # verts = [
                    #     (x_a, y_a),  # P0
                    #     (x_a_i, y_a_i),  # P1
                    #     (x_b_i, y_b_i),  # P2
                    #     (x_b, y_b),  # P3
                    # ]
                    lw_val = np.log2(cLink.score + 0.1) / 10
                    path = Path(verts, codes)
                    # patches.append(mpatches.PathPatch(path, facecolor='none', edgecolor=currcol, linewidth=lw_val,
                    #                                   alpha=0.5))
                    ax.add_patch(mpatches.PathPatch(path, facecolor=fc, edgecolor=ec, linewidth=lw_val, alpha=0.5))
                    # f_color_v.append('none')
                    # e_color_v.append(currcol)
                    # lw_v.append(np.log2(cLink.score + 0.1)/10)


def plot_rects(refObj, index):
    cfc = None
    for cfc_x in refObj.feature_tracks:
        if cfc_x.index == index:
            cfc = cfc_x

    currStart = refObj.abs_start_pos
    pTup = (refObj.chrom, refObj.ref_start, refObj.ref_end)
    for k, klist in cfc.primary_data.items():
        for x in klist:
            istart, iend = x[0], x[1]
            seg_len = pTup[2] - pTup[1]
            if refObj.direction == "+":
                normStart = currStart + max(0, istart - pTup[1])
                normEnd = currStart + min(seg_len, iend - pTup[1])

            else:
                normEnd = currStart + min(seg_len, pTup[2] - istart)
                normStart = currStart + max(0, pTup[2] - iend)

            start_angle, end_angle = start_end_angle(normStart, normEnd, total_length)
            text_angle = (start_angle + end_angle) / 2.0

            width = cfc.top - cfc.base  # height of the rectangle
            print(width, cfc.top, cfc.base)
            try:
                ctup = make_tuple("".join(x[2][2].split()))   #index | other_value | color_tuple
                if any([x > 1 for x in ctup]):
                    ctup = tuple([x / 255.0 for x in ctup])

            except ValueError:
                print("Warning - last column of bed file for rectangle is not a valid color tuple")
                ctup = None

            if ctup:
                cfc.track_props['primary_kwargs'].pop("facecolor", None)
                cfc.track_props['primary_kwargs'].pop("edgecolor", None)
                print(cfc.base, start_angle, end_angle)
                ax.add_patch(mpatches.Wedge((0, 0), cfc.base, start_angle, end_angle, facecolor=ctup, edgecolor=ctup,
                                        width=width, **cfc.track_props['primary_kwargs']))

            else:
                ax.add_patch(mpatches.Wedge((0, 0), cfc.base, start_angle, end_angle, width=width,
                                            **cfc.track_props['primary_kwargs']))


def plot_standard_IF_track(currStart, currEnd, seg_dir, pTup, cfc, curr_chrom, total_length, seg_copies, f_ind):
    print(currStart, currEnd)
    gc, gs, ge = pTup
    granularity = cfc.track_props['granularity']
    if granularity == 0:
        granularity = max(1, int((ge - gs)/10000.0))

    # print(cfc.track_max, cfc.track_min, cfc.top, cfc.base)
    height_scale_factor = (cfc.top - cfc.base)/float(cfc.track_max - cfc.track_min)

    # plot a background
    # print(cfc.track_props['background_kwargs'])
    if cfc.track_props['background_kwargs']['facecolor'] == 'auto':
        if f_ind % 2 == 1:
            cfc.track_props['background_kwargs']['facecolor'] = (0.9, 0.9, 0.9)
            if cfc.track_props['hline_kwargs']['facecolor'] == 'auto':
                cfc.track_props['hline_kwargs']['facecolor'] = 'white'

        else:
            cfc.track_props['background_kwargs']['facecolor'] = 'none'

    if cfc.track_props['hline_kwargs']['facecolor'] == 'auto':
        cfc.track_props['hline_kwargs']['facecolor'] = 'lightgrey'

    ax.add_patch(mpatches.Wedge((0, 0), cfc.top + intertrack_spacing / 2.0, 360, 0, zorder=-1,
                                width=cfc.top - cfc.base + intertrack_spacing,
                                **cfc.track_props['background_kwargs']))

    # plot the legends lines
    # legend_points = np.linspace(currStart / total_length * 2 * np.pi, (currEnd + 1) / total_length * 2 * np.pi, 10000)
    lheights = list(np.linspace(cfc.base, cfc.top, cfc.track_props['num_hlines']))
    legend_start_angle = currStart / total_length * 360
    legend_end_angle = currEnd / total_length * 360
    # if legend_end_angle < 0 and legend_start_angle > 0:
    #     legend_end_angle += 360
    # legend_start_angle, legend_end_angle = start_end_angle(currStart, currEnd, total_length)

    # print("TRACK LEGEND HEIGHTS", legend_ticks)
    for lh in lheights:
        ax.add_patch(mpatches.Wedge((0, 0), lh, legend_start_angle, legend_end_angle, zorder=1, width=0.07,
                                    **cfc.track_props['hline_kwargs']))

                                    #facecolor='k', edgecolor='k', linewidth=0,
                                    # width=bar_width / 6.0))

        # x_v, y_v = vu.polar_series_to_cartesians(legend_points, lh)
        # plt.plot(x_v, y_v, zorder=1, **cfc.track_props['hline_kwargs'])

    if cfc.track_props['indicate_zero']:
        zh = (-cfc.track_min)/(cfc.track_max - cfc.track_min) * (cfc.top - cfc.base) + cfc.base
        zeroline_kwargs = vu.create_kwargs(kwtype="Patch", facecolors=cfc.track_props['zero_facecolor'])
        ax.add_patch(mpatches.Wedge((0, 0), zh, legend_start_angle, legend_end_angle, zorder=1, width=0.12,
                                    **zeroline_kwargs))

        # plt.plot(x_v, y_v, color=cfc.track_props['indicate_zero'], linewidth=0.5, zorder=1)

    tertiary_data = []
    tertiary_style = 'lines'
    cfc.track_props['tertiary_kwargs'] = vu.create_kwargs(kwtype="Line2D", facecolors='mediumorchid',
                                                          edgecolors='mediumorchid')
    if cfc.track_props['show_segment_copy_count']:
        v = 2*seg_copies*cfc.track_props['segment_copy_count_scaling']
        tertiary_data = [[gs, ge, v]]


    hs = cfc.track_props['hide_secondary']
    if cfc.track_props['hide_secondary'] == "viral" and not (curr_chrom.startswith('chr') or len(curr_chrom) < 3):
        hs = True

    elif cfc.track_props['hide_secondary'] == "viral":
        hs = False

    for data_it, style, kwargs, smoothing, elem_ind in zip([cfc.primary_data[curr_chrom], cfc.secondary_data[curr_chrom], tertiary_data],
                                          [cfc.track_props['primary_style'], cfc.track_props['secondary_style'], tertiary_style],
                                          [cfc.track_props['primary_kwargs'], cfc.track_props['secondary_kwargs'], cfc.track_props['tertiary_kwargs']],
                                          [cfc.track_props['primary_smoothing'], cfc.track_props['secondary_smoothing'], 0],
                                          list(range(3))):

        #check if hiding secondary (hs)
        if elem_ind == 1 and hs:
            continue

        zorder = 3 if elem_ind == 1 else 2
        datalist = [(max(x[0], gs), min(x[1], ge), height_scale_factor*x[2] + cfc.base) for x in data_it]  # restrict to the coordinates of the region

        # convert the data into granular form
        sortrevdir = False if seg_dir == "+" else True
        datalist.sort(reverse=sortrevdir)
        point_data = []
        val_data = []
        for p in datalist:
            if p[1] - p[0] > granularity:
                plocs = np.linspace(p[0], p[1], int((p[1] - p[0])/granularity))
                for newp in plocs:
                    point_data.append(newp)
                    val_data.append(p[2])

                point_data.append(p[1])
                val_data.append(p[2])

            else:
                point_data.append((p[0] + p[1])/2.0)
                val_data.append(p[2])

        if smoothing > 0 and val_data and style == 'points':
            smvd = []
            for ind, h in enumerate(val_data):
                cr = val_data[max(0, ind-smoothing):min(len(val_data)-1, ind + smoothing + 1)]
                smvd.append(np.mean(cr))

            val_data = smvd

        # set the direction and convert to polars from proportional length
        if seg_dir == "+":
            normed_data = [(currStart + x - gs)/total_length * 2 * np.pi for x in point_data]
        else:
            normed_data = [(currStart + ge - x)/total_length * 2 * np.pi for x in point_data]

        # convert to cartesians
        x_v, y_v = vu.polar_series_to_cartesians(normed_data, val_data)

        # draw the points/lines
        if style == "points":
            # trying a kwargs-based method
            # print(kwargs)
            if len(x_v) > 0:
                print("Placing " + str(len(x_v)) + " points")

            plt.scatter(x_v, y_v, zorder=zorder, **kwargs)

            #plt.scatter(x_v, y_v, s=cfc.track_props['pointsize'], edgecolors='none', color=curr_color, marker='.',
            #            zorder=zorder)
            # plt.plot(x_v, y_v, linewidth=cfc.track_props['pointsize'], color=curr_color, zorder=zorder)

        elif style == "lines":
            if len(x_v) > 0:
                print("Placing " + str(len(x_v)) + " points")

            plt.plot(x_v, y_v, zorder=zorder, **kwargs)
            # plt.plot(x_v, y_v, linewidth=cfc.track_props['linewidth'], color=curr_color, zorder=zorder)
            # if seg_dir == "+":
            #     normeddata = [(currStart + x[0] - gs, currStart + x[1] - gs, x[2]) for x in datalist]
            # else:
            #     normeddata = [(currStart + ge - x[1], currStart + ge - x[0], x[2]) for x in datalist]

        elif style == "radial":
            x0_vect, y0_vect = vu.polar_series_to_cartesians(normed_data, [cfc.base]*len(normed_data))
            segs = []
            for x0, y0, x1, y1 in zip(x0_vect, y0_vect, x_v, y_v):
                segs.append([(x0, y0), (x1, y1)])

            # line_segments = LineCollection(segs, linewidths=cfc.track_props['linewidth'], colors=curr_color,
            #                                linestyle='solid')
            # line_segments = LineCollection(segs, linestyle='solid', zorder=zorder, **kwargs)
            line_segments = LineCollection(segs, linestyle="solid", zorder=zorder, colors='k', **kwargs)
            ax.add_collection(line_segments)

        else:
            print("feature_style must be either 'points', 'lines', or 'radial'\n")


def plot_interior_tracks(ref_placements):
    for ind, refObj in ref_placements.items():
        seg_coord_tup = (refObj.chrom, refObj.ref_start, refObj.ref_end)
        for cfc in refObj.feature_tracks:
            if cfc.index <= 0:
                continue

            if cfc.track_props['tracktype'] == 'standard':
                plot_standard_IF_track(refObj.abs_start_pos, refObj.abs_end_pos, refObj.direction, seg_coord_tup, cfc,
                                   refObj.chrom, total_length, refObj.seg_count, cfc.index)

            if cfc.track_props['tracktype'] == 'rects':
                plot_rects(refObj, cfc.index) # iterates over all and finds and plots


def plot_track_legend(refObj, ofpre, outer_bar, bar_width, noPDF):
    #create a figure that is ? tall
    legw = 1
    fig_l = plt.figure(figsize=(2, 8/3.0))
    ax_l = fig_l.add_subplot(111, aspect='equal')

    # plot the reference
    try:
        refcolor = chromosome_colors[refObj.chrom]
    except KeyError:
        # print("Color not found for " + chrom + ". Using red.")
        refcolor = "red"
    ax_l.add_patch(mpatches.Rectangle((0, outer_bar - bar_width), legw, bar_width/2, facecolor=refcolor,
                                      edgecolor=refcolor, lw=0.2))

    # TODO: plot the interior segments (if they exist)

    # plot the tracks (with ticks)
    for f_ind, cfc in enumerate(refObj.feature_tracks):
        # p_color = cfc.track_props['primary_color']
        # s_color = cfc.track_props['secondary_color']
        rh = cfc.top - cfc.base + intertrack_spacing
        rb = cfc.base - intertrack_spacing / 2
        ax_l.text(-2, rb + rh / 2.0, str(cfc.track_props["name"]), rotation=90, ha='center', va='center', color='k',
                  fontsize=cfc.track_props['grid_legend_fontsize'] + 2)

        if cfc.track_props['tracktype'] == 'standard':
            height_scale_factor = (cfc.top - cfc.base) / float(cfc.track_max - cfc.track_min)
            # plot a background
            # lcolor = 'lightgrey'
            ax_l.add_patch(mpatches.Rectangle((0, rb), legw, rh, zorder=-1, **cfc.track_props['background_kwargs']))

            # plot the legends lines
            lheights = list(np.linspace(cfc.base, cfc.top, cfc.track_props['num_hlines']))
            legend_ticks = list(np.linspace(cfc.track_min, cfc.track_max, cfc.track_props['num_hlines']))
            cfc.track_props['hline_kwargs']['linewidth'] *= 2
            for lh, lt in zip(lheights, legend_ticks):
                # sec_lt = (lt - cfc.minsec) * cfc.sec_rsf + cfc.sec_rss

                if cfc.track_props['rescale_secondary_to_primary']:
                    sec_lt = (lt - cfc.sec_rss) / cfc.sec_rsf + cfc.minsec  # inverse the previous transformation
                    sec_lt_str = '%s' % float('%.3g' % sec_lt)
                else:
                    sec_lt_str = str(lt)

                # ax_l.plot([0, legw], [lh, lh], zorder=1, color=lcolor, linewidth=0.5, zorder=1)
                # cfc.track_props['hline_kwargs']['color'] = cfc.track_props['hline_kwargs']['facecolor']
                # ax_l.plot([0, legw], [lh, lh], zorder=1, **cfc.track_props['hline_kwargs'])
                # print("legkw", cfc.track_props['hline_kwargs'])
                ax_l.add_patch(mpatches.Rectangle((0, lh), legw, 0.05, zorder=1, **cfc.track_props['hline_kwargs']))
                ax_l.text(-0.15, lh, str(lt), ha='right', va='center', fontsize=cfc.track_props['grid_legend_fontsize'],
                          color='k')

                if cfc.track_props['secondary_feature_bedgraph']:
                    ax_l.text(legw + 0.15, lh, sec_lt_str, ha='left', va='center', color='k',
                              fontsize=cfc.track_props['grid_legend_fontsize'])

            # background to text label will be the data color, so the font will need to be colored appropriately.
            p_color = 'k'
            for n in ['facecolors', 'markerfacecolor', 'facecolor', 'c']:
                if n in cfc.track_props['primary_kwargs']:
                    p_color = cfc.track_props['primary_kwargs'][n]

            s_color = 'k'
            for n in ['facecolors', 'markerfacecolor', 'facecolor', 'c']:
                if n in cfc.track_props['secondary_kwargs']:
                    s_color = cfc.track_props['secondary_kwargs'][n]

            r, g, b, a = matplotlib.colors.to_rgba(p_color, alpha=None)
            if r == g == b < 0.5:
                p_fc = 'white'
            else:
                p_fc = 'k'

            ax_l.text(-1.4, rb + rh/2.0, "Primary data", rotation=90, ha='center', va='center', color=p_fc,
                      fontsize=cfc.track_props['grid_legend_fontsize'] + 1, backgroundcolor=p_color)

            # print(cfc.secondary_data)
            # print("sec_check", [len(x) for x in cfc.primary_data.values()])
            # print("sec_check", [len(x) for x in cfc.secondary_data.values()])
            if cfc.track_props['secondary_feature_bedgraph']:
                r, g, b, a = matplotlib.colors.to_rgba(s_color, alpha=None)
                if r == g == b < 0.5:
                    s_fc = 'white'
                else:
                    s_fc = 'k'

                ax_l.text(legw + 1.4, rb + rh/2.0, "Secondary data", rotation=270, ha='center', va='center',
                          fontsize=cfc.track_props['grid_legend_fontsize'] + 1, color=s_fc, backgroundcolor=s_color)

            ax_l.plot([0, legw], [cfc.top + intertrack_spacing/2]*2, color='k', linewidth=0.75, zorder=1)

        # if it's links?
        # else:

    ax_l.axis('off')
    fig_l.savefig(ofpre + '.png', dpi=600)
    if noPDF is False:
        fig_l.savefig(ofpre + '.pdf', format='pdf')


def plot_gene_direction_indicator(s, e, total_length, drop, flanked, gInstance, shift):
    slant = 3.0
    if args.figure_size_style == "small":
        marker_freq = 0.015 * total_length
        clw = 0.8
    else:
        marker_freq = 0.007 * total_length
        clw = 0.4

    fullspace_a = np.arange(shift, shift + 2*total_length, marker_freq)
    fullspace_b = np.arange(shift + marker_freq/slant, shift + 2*total_length + marker_freq/slant, marker_freq)

    trim = drop / 4.0
    if drop < 0:
        # posns_a, posns_b = posns_b, posns_a
        fullspace_a, fullspace_b = fullspace_b, fullspace_a
        trim *= -1

    boolean_array = np.logical_and(fullspace_a >= s, fullspace_a <= e)
    in_range_indices = np.where(boolean_array)[0]
    # put one down if it's too skinny
    if len(in_range_indices) == 0 and not flanked:
        posns_a = [(e + s)/2.0]
        posns_b = [(e + s)/2.0 + marker_freq/slant]
        if drop < 0:
            posns_a, posns_b = posns_b, posns_a

    else:
        posns_a = [fullspace_a[x] for x in in_range_indices]
        posns_b = [fullspace_b[x] for x in in_range_indices]

    ttop = outer_bar - bar_width / 4.0 + drop - trim
    tbot = ttop - bar_width / 4.0 + trim

    btop = tbot
    bbot = tbot - bar_width/ 4.0 + trim

    for fpos, rpos in zip(posns_a, posns_b):
        pos_angle_a = fpos / total_length * 360
        pos_angle_b = rpos / total_length * 360


        x_b, y_b = vu.pol2cart(ttop, (pos_angle_a / 360 * 2 * np.pi))
        x_t, y_t = vu.pol2cart(tbot, (pos_angle_b / 360 * 2 * np.pi))
        plt.plot([x_b, x_t], [y_b, y_t], linewidth=clw, color='grey')

        x_b, y_b = vu.pol2cart(btop, (pos_angle_b / 360 * 2 * np.pi))
        x_t, y_t = vu.pol2cart(bbot, (pos_angle_a / 360 * 2 * np.pi))
        plt.plot([x_b, x_t], [y_b, y_t], linewidth=clw, color='grey')

    #draw marker starts and ends
    gInstance.draw_marker_ends(tbot)


def plot_gene_bars(currStart, currEnd, relGenes, pTup, total_length, seg_dir, ind, flanked, cycle, isCycle, shift,
                   plot_gene_direction=True):
    overlap_genes.append({})
    nhits = len(relGenes) + 1
    for gObj_ind, gObj in enumerate(sorted(relGenes, key=lambda x: x.gstart)):
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
                if ts > 0: hasStart = True
                if te < seg_len: hasEnd = True
            else:
                if ts > 0: hasEnd = True
                if te < seg_len: hasStart = True

            normStart = currStart + max(0, gstart - pTup[1])
            normEnd = currStart + min(seg_len, gend - pTup[1])

        else:
            te = min(seg_len, pTup[2] - gstart)
            ts = max(0, pTup[2] - gend)
            if gObj.strand == "+":
                if te < seg_len: hasStart = True
                if ts > 0: hasEnd = True
            else:
                if te < seg_len: hasEnd = True
                if ts > 0: hasStart = True

            normEnd = currStart + min(seg_len, pTup[2] - gstart)
            normStart = currStart + max(0, pTup[2] - gend)

        start_angle = normStart / total_length * 360
        end_angle = normEnd / total_length * 360
        # text_angle = (gObj_ind+1)*(start_angle + end_angle) / float(n)
        text_angle = start_angle + (gObj_ind+1)*(end_angle-start_angle)/float(nhits)
        gene_to_locations[gname].append((start_angle / 360., end_angle / 360.))
        if end_angle < 0 and start_angle > 0:
            end_angle += 360

        gsign = 1 if seg_dir == gObj.strand else -1
        drop = gsign * bar_width / 4.0
        gbh = outer_bar - 5.0*bar_width/12 + drop
        gObj.gdrops.append(gbh)
        ax.add_patch(mpatches.Wedge((0, 0), gbh, start_angle, end_angle, facecolor='k', edgecolor='k', linewidth=0,
                                      width=bar_width / 6.0))

        # TODO: REFACTOR TO OUTSIDE - put in the gParent
        if gname not in overlap_genes[len(overlap_genes)-2] or not overlap_genes[len(overlap_genes)-2].get(gname)[0] \
                or seg_dir != overlap_genes[len(overlap_genes)-2].get(gname)[1]:
            #handle wraparound
            if not (isCycle and ind == len(cycle)-1 and gname in overlap_genes[0]):
                x_t, y_t = vu.pol2cart(outer_bar + bar_width + gene_spacing, (text_angle / 360 * 2 * np.pi))
                text_angle, ha = vu.correct_text_angle(text_angle)
                if gObj.highlight_name:
                    ax.text(x_t, y_t, gname, style='italic', color='r', rotation=text_angle, ha=ha, va="center",
                            fontsize=gene_fontsize, rotation_mode='anchor')
                else:
                    ax.text(x_t, y_t, gname, style='italic', color='k', rotation=text_angle, ha=ha, va="center",
                            fontsize=gene_fontsize, rotation_mode='anchor')

        # draw something to show direction and truncation status
        if plot_gene_direction:
            # gene instance
            gInstance = vu.gene_viz_instance(gObj, normStart, normEnd, total_length, seg_dir, currStart,
                                             currEnd,
                                             hasStart, hasEnd, ind, pTup)

            gObj.gdrops.append(gInstance)
            plot_gene_direction_indicator(normStart, normEnd, total_length, drop, flanked, gInstance, shift)
            # gObj.gdrops = [(normStart, normEnd, total_length, seg_dir, currStart, currEnd, pTup), ]

        if not (pTup[2] >= gend and pTup[1] <= gstart):
            overlap_genes[len(overlap_genes)-1][gname] = (True, seg_dir)

        for exon in e_posns:
            if gObj.highlight_name:
                ecolor = 'r'
                lw = 0.3
            else:
                ecolor = 'k'
                lw = 0.3
            if exon[1] > pTup[1] and exon[0] < pTup[2]:
                if seg_dir == "+":
                    normStart = currStart + max(1, exon[0] - pTup[1])
                    normEnd = currStart + min(pTup[2] - pTup[1], exon[1] - pTup[1])

                else:
                    normEnd = currStart + min(pTup[2] - pTup[1], pTup[2] - exon[0])
                    normStart = currStart + max(1, pTup[2] - exon[1])

                start_angle, end_angle = start_end_angle(normStart, normEnd, total_length)

                ax.add_patch(
                    mpatches.Wedge((0, 0), outer_bar - bar_width / 4.0 + (drop), start_angle, end_angle,
                                   facecolor=ecolor, edgecolor=ecolor, linewidth=lw, width=bar_width / 2.0))


# Gene plotting
def plot_genes(ref_placements, cycle, isCycle, onco_set=None):
    if onco_set is None:
        onco_set = set()

    shift = min([x.abs_start_pos for x in ref_placements.values()])
    for ind, refObj in ref_placements.items():
        seg_coord_tup = (refObj.chrom, refObj.ref_start, refObj.ref_end)
        relGenes = vu.rel_genes(gene_tree, seg_coord_tup, copy.copy(onco_set))
        all_relGenes.extend(relGenes)
        # plot the gene track
        # print(ind, refObj.to_string(), len(relGenes))
        flanked = refObj.next_is_adjacent or refObj.prev_is_adjacent
        plot_gene_bars(refObj.abs_start_pos, refObj.abs_end_pos, relGenes, seg_coord_tup, total_length, cycle[ind][1],
                       ind, flanked, cycle, isCycle, shift)


# plot the reference genome
def plot_ref_genome(ref_placements, cycle, total_length, imputed_status, label_segs, edge_ticks, bw=bar_width):
    font0 = FontProperties()
    # rot_sp = global_rot / 360. * total_length
    prev_was_skinny = False
    length_since_last_bp = 0
    for ind, refObj in ref_placements.items():
        # print(ind, refObj.to_string())
        if refObj.custom_bh:
            curr_bh = refObj.custom_bh + refObj.track_height_shift
        else:
            curr_bh = outer_bar + refObj.track_height_shift

        # seg_coord_tup = segSeqD[cycle[ind][0]]
        # print(ind, refObj.to_string())
        chrom = refObj.chrom
        seg_coord_tup = (refObj.chrom, refObj.ref_start, refObj.ref_end)
        start_angle, end_angle = start_end_angle(refObj.abs_end_pos, refObj.abs_start_pos, total_length)

        # makes the reference genome wedges
        if not refObj.custom_face_color:
            if args.structure_color == "auto":
                if chrom not in chromosome_colors:
                    print("Color not found for " + chrom + ". Using red.")
                    chromosome_colors[chrom] = "red"

                f_color = chromosome_colors[chrom]
                e_color = chromosome_colors[chrom]
            else:
                f_color, e_color = args.structure_color, args.structure_color
                if e_color == f_color and (e_color == 'w' or e_color == 'white'):
                    e_color = 'k'

        # this happens with the interior track!
        else:
            f_color = refObj.custom_face_color
            e_color = refObj.custom_edge_color

        # lw_v.append(0.2)
        ax.add_patch(mpatches.Wedge((0, 0), curr_bh, end_angle, start_angle, facecolor=f_color, edgecolor=e_color,
                                    linewidth=0.2, width=bw))

        # makes the ticks on the reference genome wedges
        # TODO: Refactor outside
        if cycle[ind][1] == "+":
            ts = (seg_coord_tup[1], refObj.abs_start_pos, 0)
            te = (seg_coord_tup[2] + 1, refObj.abs_end_pos+1, 1)
            s = 1

        else:
            ts = (seg_coord_tup[2], refObj.abs_start_pos, 0)
            te = (seg_coord_tup[1] - 1, refObj.abs_end_pos + 1, -1)
            s = -1

        if not refObj.prev_is_adjacent:
            length_since_last_bp = te[1] - ts[1]
        else:
            length_since_last_bp+=(te[1] - ts[1])

        text_trunc = 1
        # put the positions on the ends of the joined segs
        posns = []
        if edge_ticks == "ends":
            newposns = []
            tick_freq = 1
            # if too narrow, just make one label in the middle
            # print(ts[0],te[0], (te[1] - ts[1])/total_length * 360)
            skinny = length_since_last_bp / total_length * 360 < 1.2
            # if it's very small, just put one marker in the center
            if skinny and not refObj.prev_is_adjacent and not refObj.next_is_adjacent:
                skinny = False
                tm = (str(ts[0]) + "-" + str(te[0]), np.mean((ts[1],te[1])), 0)
                newposns.append(tm)

            else:
                if not refObj.prev_is_adjacent:
                    newposns.append(ts)

                if not refObj.next_is_adjacent:
                    skinny = False
                    newposns.append(te)

                if refObj.prev_is_adjacent and refObj.next_is_adjacent:
                    skinny = False

            posns = newposns

        elif edge_ticks == "none":
            skinny = False
            tick_freq = float('inf')

        else:
            skinny = False
            text_trunc = 10000
            tick_freq = max(10000, 20000 * int(np.floor(total_length / 1000000)))
            # print("tick freq", tick_freq)
            step = int(tick_freq/10000)
            a = int(np.floor(ts[0] / 10000)) + 1
            b = int(np.floor(te[0] / 10000)) + 1
            # print(a, b, ts[0], te[0], step)
            for j in np.arange(a, b, s):
                if (j*10000) % tick_freq == 0:
                    sj = j*10000
                    rpos = vu.convert_gpos_to_ropos(sj, refObj.abs_start_pos, refObj.abs_end_pos, seg_coord_tup[1],
                                                    cycle[ind][1])
                    posns.append((sj, rpos, 0))

        for j in posns:
            text_angle = j[1] / total_length * 360
            x, y = vu.pol2cart(curr_bh, (text_angle / 360 * 2 * np.pi))
            x_t, y_t = vu.pol2cart(curr_bh + 0.2, (text_angle / 360 * 2 * np.pi))
            ax.plot([x, x_t], [y, y_t], color='grey', linewidth=1, zorder=-10)

            text_angle, ha = vu.correct_text_angle(text_angle)
            txtpad = " "
            if prev_was_skinny and len(posns) == 1 and skinny:
                txtpad += " "*(2*len(str(j[0])))
                txtpad = txtpad + "-" if ha == "left" else "-" + txtpad

            try:
                txt = txtpad + str(int(round((j[0]-j[2]) / text_trunc))) if ha == "left" else str(int(round(
                    (j[0]-j[2]) / text_trunc))) + txtpad
            except TypeError:
                txt = txtpad + str(j[0]) if ha == "left" else str(j[0]) + txtpad

            ax.text(x_t, y_t, txt, color='grey', rotation=text_angle,
                    ha=ha, va="center", fontsize=tick_fontsize, rotation_mode='anchor')

        prev_was_skinny = skinny
        # end ticking section

        # Segment labeling
        # TODO: Refactor to outside
        # label the segments by number in cycle
        if label_segs:
            # mid_sp = (refObj.abs_end_pos + refObj.abs_start_pos) / 2.0
            # centerpoint_angle = mid_sp / total_length * 360.
            centerpoint_angle = (start_angle + end_angle) / 2
            x, y = vu.pol2cart((curr_bh - ref_label_height_factor), (centerpoint_angle / 360. * 2. * np.pi))
            font = font0.copy()
            if imputed_status[ind]:
                font.set_style('italic')
                # font.set_weight('bold')

            if label_segs[0] == "numbers":
                t = str(cycle[ind][0]) + cycle[ind][1]
                text_angle, ha = vu.correct_text_angle(centerpoint_angle)
                va = 'baseline'

            elif label_segs[0] == "names":
                t = refObj.chrom.lstrip("chr")
                text_angle, temp = vu.correct_text_angle(centerpoint_angle + 90)
                ha = 'center'
                if temp == 'right':
                    va = 'bottom'
                else:
                    va = 'top'

            else:
                try:
                    t = label_segs[ind]
                except IndexError:
                    t = ""

                text_angle, temp = vu.correct_text_angle(centerpoint_angle + 90)
                ha = 'center'
                if temp == 'right':
                    va = 'bottom'
                else:
                    va = 'top'

            ax.text(x, y, t, color='grey', rotation=text_angle, ha=ha, va=va, fontsize=seg_label_fontsize,
                    fontproperties=font, rotation_mode='anchor')


# set the heights of the bed track features
def get_feature_heights(ntracks, intertrack_spacing, has_OM, has_IS):
    IS_height = segment_bar_height/2
    if ntracks > 0 or has_IS:
        maxtop = outer_bar-(intertrack_spacing + feature_ref_offset)
        if has_OM:
            maxtop += contig_bar_height
        if has_IS:
            maxtop += IS_height

        if maxtop <= center_hole:
            print("ERROR: om and segment height exceeds allowed height in track")
            sys.exit(1)

        divs = np.linspace(center_hole, maxtop, ntracks+1)
        bases = [divs[0] + intertrack_spacing, ]
        tops = []
        for p in divs[1:-1]:
            tops.append(p)
            bases.append(p + intertrack_spacing)

        tops.append(divs[-1])
        # establish the top of the segment bar
        smt = maxtop
        if has_IS:
            smt = maxtop - IS_height

        print("Intertrack spacing is ", bases, tops)
        return bases, tops, smt

    return [], [], 0


# plot cmap track for bionano
def plot_cmap_track(seg_placements, total_length, unadj_bar_height, color, seg_id_labels=False):
    cycle_label_locs = defaultdict(list)
    for ind, segObj in seg_placements.items():
        bar_height = unadj_bar_height + segObj.track_height_shift
        # print("cmap_plot", segObj.id)
        if segObj.abs_end_pos - segObj.abs_start_pos > total_length and segObj.ref_end - segObj.ref_start > total_length:
            start_angle, end_angle = 360, 0
            print("cmap went full circle")
            # print(segObj.to_string())

        else:
            start_angle, end_angle = start_end_angle(segObj.abs_end_pos, segObj.abs_start_pos, total_length)
        ax.add_patch(mpatches.Wedge((0, 0), bar_height + bar_width, end_angle, start_angle, facecolor=color,
                                      edgecolor='k', linewidth=0, width=bar_width))

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
            print(segObj.id, start_angle, end_angle, segObj.abs_end_pos, segObj.abs_start_pos)
            text_angle = (start_angle + end_angle) / 2
            print(segObj.to_string(), text_angle)
            x, y = vu.pol2cart(bar_height - 1.5, (text_angle / 360. * 2. * np.pi))
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
        # print(a_d["contig_id"],a_d["seg_id"],a_d["contig_label"],a_d["seg_label"])
        s_l_pos = seg_label_vect[a_d["seg_label"] - 1]
        s_l_loc = s_l_pos / total_length * 2. * np.pi
        contig_top = outer_bar + contig_bar_height + contig_locs[c_id].track_height_shift + bar_width
        x_c, y_c = vu.pol2cart(contig_top, c_l_loc)
        x_s, y_s = vu.pol2cart(segs_base, s_l_loc)
        ax.plot([x_c, x_s], [y_c, y_s], color="grey", linewidth=linewidth)


# identify the offset of the segment with minimum location
def compute_min_seg_offset(cycle, segSeqD, spacing_bp, prev_seg_index_is_adj, isCycle):
    minChrom = sys.maxsize
    chromLook = {"X": 23, "Y": 24, "M": 25}
    minLoc = float('inf')
    curr_start = 0.0 if isCycle else spacing_bp
    next_start = 0
    minSegCStart = curr_start
    for ind, i in enumerate(cycle):
        cc, ca, cb = segSeqD[i[0]]
        cc = cc.split("chr", 1)[-1]
        try:
            cc = int(cc)
            tc = cc

        except ValueError:
            if cc in chromLook:
                tc = chromLook[cc]
            else:
                tc = 0

        dir = i[1]
        if dir == "-":
            cs = cb
        else:
            cs = ca

        seg_len = cb - ca
        seg_end = curr_start + seg_len

        if tc < minChrom or (cc == minChrom and cs < minLoc):
            minChrom, minLoc, minSegCStart = tc, cs, curr_start

        next_start = seg_end
        mod_ind = (ind + 1) % (len(prev_seg_index_is_adj))
        if not prev_seg_index_is_adj[mod_ind]:
            next_start += spacing_bp

        curr_start = next_start

    total_length = next_start
    return minSegCStart, total_length


def construct_cycle_ref_placements(cycle, segSeqD, raw_cycle_length, prev_seg_index_is_adj, next_seg_index_is_adj,
                                   isCycle, cycle_seg_counts, rotate_to_min):
    spacing_bp = seg_spacing * raw_cycle_length
    cycle_ref_placements = {}
    curr_start = 0.0 if isCycle else spacing_bp
    cso, total_length = compute_min_seg_offset(cycle, segSeqD, spacing_bp, prev_seg_index_is_adj, isCycle)
    if rotate_to_min:
        curr_start-=cso

    next_start = 0
    for ind, i in enumerate(cycle):
        seg_id_count = cycle_seg_counts[i[0]]
        seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        # print(i, seg_len)
        seg_end = (curr_start + seg_len) % total_length
        padj, nadj = prev_seg_index_is_adj[ind], next_seg_index_is_adj[ind]
        curr_obj = vu.CycleVizElemObj(i[0], segSeqD[i[0]][0], segSeqD[i[0]][1], segSeqD[i[0]][2], i[1], curr_start,
                                      seg_end, seg_id_count, padj, nadj)
        cycle_ref_placements[ind] = curr_obj
        next_start = seg_end
        mod_ind = (ind + 1) % (len(prev_seg_index_is_adj))
        if not prev_seg_index_is_adj[mod_ind]:
            next_start += spacing_bp

        curr_start = next_start % total_length

    return cycle_ref_placements, total_length


parser = argparse.ArgumentParser(description="Circular visualizations of genome structures")
group = parser.add_mutually_exclusive_group()
group.add_argument("--input_yaml_file", help="Specifiy all desired arguments in this file, OR use the options below\n")
group.add_argument("--structure_bed", help="bed file specifying the structure of the regions to be plotted. To use a "
                                           "standard reference genome as the structure, specify 'hg19, GRCh37, hg38, "
                                           "GRCh38, mm10, or GRCm38")
parser.add_argument("--cycles_file", help="AA/AR cycles-formatted input file")
parser.add_argument("--cycle", help="cycle number to visualize [required with --cycles_file]", type=int)
parser.add_argument("-g", "--graph", help="breakpoint graph file [required with --cycles_file]")
parser.add_argument("--ref", help="reference genome", choices=["hg19", "hg38", "GRCh37", "GRCh38", "mm10", "GRCh38_viral"])
parser.add_argument("--om_alignments",
                    help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
                    action='store_true')
parser.add_argument("--interior_segments_cycle",
                    help="Enable visualization of an interior segments from an AA cycles file (e.g. long read", type=str, default="")
parser.add_argument("--interior_segments_cycle_connect_width", help="Width of drawn junction between segments of interior cycle",
                    choices=["full", "thin", "none"], default="full")
parser.add_argument("--interior_segments_cycle_matching", help="Specifying matching of interior cycle to outer cycle should be"
                    " exact or relaxed (local matches with shared start). Required if --interior_segments_cycle set", choices=["exact", "relaxed"])
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("--om_segs", help="segments cmap file")
parser.add_argument("--AR_path_alignment", help="AR path alignment file")
parser.add_argument("--outname", "-o", help="output prefix")
# parser.add_argument("--rot", help="number of segments to rotate counterclockwise", type=int, default=0)
parser.add_argument("--label_segs", help="label segs with segments number-direction ('numbers') or segment names "
                                         "('names') or a custom list", default="", type=str, nargs='+')
parser.add_argument("--seg_label_fontsize", help="reference segment label fontsize (default 4)", type=float, default=4)
parser.add_argument("--seg_label_height", help="Backbone segment label height (default 2.2)", type=float, default=2.2)
parser.add_argument("--gene_subset_file", help="file containing subset of genes to plot (e.g. oncogene genelist file)",
                    default="")
parser.add_argument("--gene_subset_list", help="list of genes to plot (e.g. MYC PVT1)", nargs="+", type=str)
parser.add_argument("--print_dup_genes", help="if a gene appears multiple times print name every time.",
                    action='store_true', default=False)
parser.add_argument("--gene_highlight_list", help="list of gene names to highlight", nargs="+", type=str, default=[])
parser.add_argument("--gene_fontsize", help="font size for gene names", type=float, default=7)
parser.add_argument("--gene_spacing", help="How far from reference to plot gene names. Default 1.7", type=float,
                    default=1.7)
parser.add_argument("--tick_type", help="Represent ticks only at segment ends, or throughout ('ends' is default)",
                    choices=["ends", "standard", "none"], default="ends")
parser.add_argument("--tick_fontsize", help="font size for genomic position ticks", type=float)
parser.add_argument("--feature_yaml_list", nargs='+', help="list of the input yamls for bedgraph file feature "
                    "specifying additional data. Will be plotted from inside to outside given the order the filenames "
                    "appear in", default=[])
parser.add_argument("--annotate_structure", help="What to plot on the outer structure indicator. Give a bed file or use"
                    " predefined 'genes' argument", type=str, default="genes")
parser.add_argument("--structure_color", help="Use 'auto' coloring, or specify a single color for everything",
                    type=str, default='auto')
parser.add_argument("--rotate_to_min", help="Rotate structure so the lowest genomic coordinate appears at 0-degrees in "
                    "along first quadrant (3 o'clock position)", action='store_true', default=False)
parser.add_argument("--center_hole", type=float, help="whitespace in center of plot (Default 1.25)", default=1.25)
parser.add_argument("--intertrack_spacing", help="spacing between feature tracks (default 0.7)", type=float, default=0.7)
parser.add_argument("--feature_ref_offset", type=float, help="Offset spacing between ref and feature tracks (Default 0.5)",
                    default=0.5)
parser.add_argument("--figure_size_style", choices=["normal", "small"], default="normal")
parser.add_argument("--noPDF", help="Do not generate PDF file of visualization (default unset)", action='store_true',
                    default=False)
parser.add_argument("--hide_chrom_color_legend", help="Do not show a legend of the chromosome colors",
                    action='store_true', default=False)
parser.add_argument("--trim_contigs", help="Trim unaligned OM contig regions from visualization (defauly unset)",
                    default='both', choices=['start', 'end', 'both', 'none'])
parser.add_argument("-v", "--version", action='version', version='CycleViz {version} \n Author: Jens Luebeck '
                    '(jluebeck [at] ucsd.edu)'.format(version=__version__))

args = parser.parse_args()
if args.input_yaml_file:
    vu.parse_main_args_yaml(args)
    if not args.cycles_file and not args.structure_bed:
        print("--input_yaml_file must provide a cycles_file or structure_bed key!")

if not args.cycles_file and not args.input_yaml_file and not args.structure_bed:
    print("One of --input_yaml_file, --cycles_file, --structure_bed required!")

if args.cycles_file and not args.cycle:
    print("Error: must specify --cycle [number] if providing --cycles_file")
    sys.exit(1)


ref_choices = ["hg19", "hg38", "GRCh37", "GRCh38", "mm10", "GRCh38_viral"]
if not args.ref:
    print("--ref is required")

elif args.ref not in ref_choices:
    print("--ref must be one of " + str(ref_choices))


if args.ref == "GRCh38" or args.ref == "GRCh38_viral":
    args.ref = "hg38"

print("Plots will use reference " + args.ref)

if not args.tick_fontsize:
    if not args.tick_type == 'ends':
        args.tick_fontsize = 3
    else:
        args.tick_fontsize = 4

if args.figure_size_style == "small":
    print("Using figure_size_style=small")
    bar_width *= 1.5
    args.gene_fontsize *= 2
    args.tick_fontsize *= 1.5
    print("Bar width scaled up 1.5x to " + str(bar_width))
    print("Gene fontsize scaled up 2x to " + str(args.gene_fontsize))
    print("Tick fontsize scaled up 1.5x to: " + str(args.tick_fontsize))

ref_label_height_factor = args.seg_label_height * bar_width
center_hole = args.center_hole
gene_spacing = args.gene_spacing
feature_ref_offset = args.feature_ref_offset
intertrack_spacing = args.intertrack_spacing
seg_label_fontsize = args.seg_label_fontsize

sourceDir = os.path.dirname(os.path.abspath(__file__)) + "/"

# print("Unaligned fraction cutoff set to " + str(vu.unaligned_cutoff_frac))
chromosome_colors = vu.get_chr_colors()
plt.clf()
fig, ax = plt.subplots()
patches = []
f_color_v = []
e_color_v = []
lw_v = []

gene_fontsize = args.gene_fontsize
tick_fontsize = args.tick_fontsize
bpg_dict, seg_end_pos_d = {}, {}

# use AA files to determine the structure
if args.cycles_file:
    if not args.outname:
        args.outname = os.path.splitext(os.path.basename(args.cycles_file))[0]
    fname = args.outname + "_cycle_" + str(args.cycle)
    cycles, segSeqD, circular_D = vu.parse_cycles_file(args.cycles_file)
    isCycle = circular_D[args.cycle]
    cycle = cycles[args.cycle]
    if args.graph:
        bpg_dict, seg_end_pos_d = vu.parse_BPG(args.graph)

# use the structure_bed format to determine the structure
else:
    if not args.outname:
        print("Must specify --outname with --structure-bed")
        sys.exit(1)
        # args.outname = os.path.splitext(os.path.basename(args.structure_bed))[0] + "_"
    fname = args.outname + "_cycle_1"
    if args.structure_bed in {"hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38"}:
        args.structure_bed = CV_RESOURCES + args.structure_bed + "_structure.bed"
    struct_data = vu.parse_bed(args.structure_bed, store_all_additional_fields=True)
    cycle, isCycle, segSeqD, seg_end_pos_d, bpg_dict = vu.handle_struct_bed_data(struct_data)

# bookkeeping
cycle_seg_counts = vu.get_seg_amplicon_count(cycle)
prev_seg_index_is_adj, next_seg_index_is_adj = vu.adjacent_segs(cycle, segSeqD, isCycle)
raw_cycle_length = vu.get_raw_path_length(cycle, segSeqD)

# determine which genes to show
gene_set = set()
if args.gene_subset_file.upper() == "BUSHMAN" and args.ref in {"hg19", "GRCh37", "hg38", "GRCh38"}:
    args.gene_subset_file = CV_RESOURCES + "Bushman_group_allOnco_May2018.tsv"

if args.gene_subset_file:
    gff = True if args.gene_subset_file.endswith(".gff") else False
    gene_set = vu.parse_gene_subset_file(args.gene_subset_file, gff)

elif args.gene_subset_list:
    gene_set = set(args.gene_subset_list)

fbases, ftops, IS_bh = get_feature_heights(len(args.feature_yaml_list), intertrack_spacing, args.om_alignments,
                                    args.interior_segments_cycle)

# standard case
if not args.om_alignments:
    ref_placements, total_length = construct_cycle_ref_placements(cycle, segSeqD, raw_cycle_length,
                                                                  prev_seg_index_is_adj, next_seg_index_is_adj,
                                                                  isCycle, cycle_seg_counts, args.rotate_to_min)
    imputed_status = [False] * len(cycle)

# only if bionano data present
else:
    print("Visualizing with alignments")
    print("Contig spacing set to " + str(vu.contig_spacing))
    seg_cmaps = parse_cmap(args.om_segs, True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    seg_cmap_lens = get_cmap_lens(args.om_segs)
    aln_vect, meta_dict = vu.parse_alnfile(args.AR_path_alignment)
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
                                                                  cycle_seg_counts, args.rotate_to_min)

    cycle_seg_placements = vu.place_path_segs_and_labels(cycle, ref_placements, seg_cmap_vects)
    contig_cmaps = parse_cmap(args.contigs, True)
    contig_cmap_vects = vectorize_cmaps(contig_cmaps)
    contig_cmap_lens = get_cmap_lens(args.contigs)
    contig_placements, contig_list = vu.place_contigs_and_labels(cycle_seg_placements, aln_vect, total_length,
                                                                 contig_cmap_vects, isCycle, True, segSeqD)

    if args.trim_contigs != 'none':
        trim_s, trim_e = False, False
        if args.trim_contigs == "both" or args.noTrim == "start":
            trim_s = True
        if args.trim_contigs == "both" or args.noTrim == "end":
            trim_e = True

        vu.decide_trim_contigs(contig_cmap_vects, contig_placements, total_length, trim_s, trim_e)

    # plot cmap segs
    print("Plotting graph segment CMAPs")
    plot_cmap_track(cycle_seg_placements, total_length, outer_bar + segment_bar_height, "darkorange")

    # check overlaps of contigs and adjust heights accordingly
    contig_height_shifts = vu.set_contig_height_shifts(contig_placements, contig_list)
    # plot contigs
    print("Plotting contig CMAPs")
    plot_cmap_track(contig_placements, total_length, outer_bar + contig_bar_height, "cornflowerblue",
                    seg_id_labels=True)

    # plot alignments
    plot_alignment(contig_placements, cycle_seg_placements, total_length)
    imputed_status = vu.imputed_status_from_aln(aln_vect, len(cycle))

print("plotting structure")
plot_ref_genome(ref_placements, cycle, total_length, imputed_status, args.label_segs, args.tick_type, bar_width)
if args.annotate_structure == 'genes' or gene_set:
    print("Reading genes")
    gene_tree = vu.parse_genes(args.ref, args.gene_highlight_list)
    print("plotting genes")
    plot_genes(ref_placements, cycle, isCycle, gene_set)


# Interior segments
if args.interior_segments_cycle:
    if not args.interior_segments_cycle_matching:
        print("Must set --interior_segments_cycle_matching [exact, relaxed] if providing --interior_segments_cycle_matching")
        sys.exit(1)

    covered_posns = IntervalTree()
    interior_cycles, interior_segSeqD, IS_circular_D = vu.parse_cycles_file(args.interior_segments_cycle)

    n_cycles = len(interior_cycles.keys())
    visualized_cycles = 0
    for IS_cnum in sorted(interior_cycles.keys()):
        interior_cycle, IS_isCircular = interior_cycles[IS_cnum], IS_circular_D[IS_cnum]

        # NEED TO ALSO CHECK REVERSE DIR
        hit = False
        for curr_interior_cycle in [interior_cycle, vu.reverse_cycle(interior_cycle)]:
            if args.interior_segments_cycle_matching == "relaxed":
                IS_rObj_placements, new_interior_cycle, new_IS_links = vu.handle_IS_data_relaxed(ref_placements,
                                                                                          curr_interior_cycle,
                                                                                          interior_segSeqD,
                                                                                          IS_isCircular, IS_bh)

                if IS_rObj_placements:
                    visualized_cycles+=1
                    plot_ref_genome(IS_rObj_placements, new_interior_cycle, total_length, [False] * len(new_interior_cycle),
                                    False, 'none',
                                    bar_width)

                    # print(IS_rObj_placements.keys())
                    plot_bpg_connection(IS_rObj_placements, total_length, manual_links=new_IS_links,
                                        connect_width=args.interior_segments_cycle_connect_width)

                else:
                    print(IS_cnum, "no hit")

            else: # strict mode
                print("\nTrying to match", curr_interior_cycle)
                for IS_rObj_placements, new_interior_cycle, new_IS_links in vu.handle_IS_data(ref_placements, cycle, isCycle,
                                                                        curr_interior_cycle, interior_segSeqD, IS_isCircular, IS_bh):
                    # print(new_interior_cycle)
                    # print(len(IS_rObj_placements), len(new_interior_cycle))
                    hit = True
                    max_olaps = 0
                    used_rows = set()
                    for currObj in IS_rObj_placements.values():
                        abs_start, abs_end = int(currObj.abs_start_pos)-1, int(currObj.abs_end_pos) + 1
                        olaps = covered_posns[abs_start:abs_end]
                        for o in olaps:
                            used_rows.add(o.data)
                        max_olaps = max(len(olaps), max_olaps)

                    openings = sorted(list(set(range(0, max_olaps+1)).difference(used_rows)))
                    # print(used_rows)
                    # print(max_olaps)
                    # print(openings)
                    if openings:
                        max_olaps = openings[0]

                    temp_interval_tree = IntervalTree()
                    for currObj in IS_rObj_placements.values():
                        currObj.track_height_shift = max_olaps * -1
                        abs_start, abs_end = int(currObj.abs_start_pos)-1, int(currObj.abs_end_pos) + 1
                        temp_interval_tree.addi(abs_start, abs_end, max_olaps)

                    temp_interval_tree.merge_overlaps()
                    new_interval_tree = IntervalTree()
                    for x in temp_interval_tree:
                        new_interval_tree.addi(x.begin, x.end, max_olaps)

                    # print(new_interval_tree)
                    covered_posns |= new_interval_tree

                    plot_ref_genome(IS_rObj_placements, new_interior_cycle, total_length, [False] * len(new_interior_cycle), False, 'none',
                                    bar_width)

                    # print(IS_rObj_placements.keys())
                    plot_bpg_connection(IS_rObj_placements, total_length, manual_links=new_IS_links,
                                        connect_width=args.interior_segments_cycle_connect_width)

        if hit:
            visualized_cycles+=1
        else:
            print(IS_cnum, "no hit")

        print(str(visualized_cycles) + " paths of the " + str(n_cycles) + " given interior paths were visualized")

# cytobanding
if args.annotate_structure != "genes":
    annot_struct_yaml = args.annotate_structure
    # if args.annotate_structure == 'cytoband':
    #     annot_struct_yaml = CV_RESOURCES + "cytoBand_" + args.ref + ".yaml"

    cfc = vu.parse_feature_yaml(annot_struct_yaml, 0, 1, CV_RESOURCES)
    cfc.base, cfc.top = outer_bar, outer_bar+bar_width
    vu.store_bed_data(cfc, ref_placements, cfc.track_props['end_trim'])
    for refObj in ref_placements.values():
        plot_rects(refObj, 0)

if args.feature_yaml_list:
    for ind, yaml_file in enumerate(args.feature_yaml_list):
        print("Processing " + yaml_file)
        cfc = vu.parse_feature_yaml(yaml_file, ind + 1, len(args.feature_yaml_list))
        cfc.base, cfc.top = fbases[ind], ftops[ind]
        vu.store_bed_data(cfc, ref_placements, cfc.track_props['end_trim'])
        if cfc.track_props['tracktype'] == 'standard' or cfc.track_props['tracktype'] == 'rects':
            vu.reset_track_min_max(ref_placements, ind, cfc)
            plot_interior_tracks(ref_placements)
        else:
            plot_links(cfc)

if bpg_dict:
    plot_bpg_connection(ref_placements, total_length, prev_seg_index_is_adj, bpg_dict, seg_end_pos_d)

ax.set_xlim(-(outer_bar + 1.25), (outer_bar + 1.25))
ax.set_ylim(-(outer_bar + 3.3), (outer_bar + 3.3))

if not args.hide_chrom_color_legend and args.structure_color == 'auto':
    chrom_set = set()
    for i in cycle:
        chrom_set.add(segSeqD[i[0]][0])

    sorted_chrom = sorted(chrom_set, key=lambda x: x.rsplit("chr")[-1])
    sorted_chrom_colors = [chromosome_colors[x] for x in sorted_chrom]
    legend_patches = []
    for chrom, color in zip(sorted_chrom, sorted_chrom_colors):
        legend_patches.append(mpatches.Patch(facecolor=color, label=chrom))

    plt.legend(handles=legend_patches, fontsize=8, loc=3, bbox_to_anchor=(-.3, .15), frameon=False)

ax.set_aspect(1.0)
plt.axis('off')

print("saving PNG")
plt.savefig(fname + '.png', dpi=600)
if args.noPDF is False:
    print("saving PDF")
    plt.savefig(fname + '.pdf', format='pdf')
    plt.close()

# make plots of the yaml tracks
if args.feature_yaml_list:
    print("saving legend")
    plot_track_legend(ref_placements[0], fname + "_legend", outer_bar, bar_width, args.noPDF)

print("finished\n")