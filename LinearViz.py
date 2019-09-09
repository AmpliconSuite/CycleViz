#!/usr/bin/env python

import sys
import os
import copy
import bisect
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
from collections import defaultdict
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
# import matplotlib._color_data as mcd
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
from intervaltree import Interval, IntervalTree
import VizUtil as vu

#non-interactive backend

try:
    sys.path.insert(0,os.environ['AR_SRC'])

except KeyError:
    sys.stderr.write("AmpliconReconstructor source directory bash variable ($AR_SRC) not found.\n")
    sys.stderr.write("Is AmpliconReconstructor configured?")
    sys.exit()

from bionanoUtil import *

seg_spacing = 0.009
bar_width_scaling = 0.02
bar_drop_prop = 2
seg_bar_base = 0 
bed_spacing = .5
contig_bar_height = 1
segment_bar_height = 0
ref_bar_height = -2
gene_bar_height = -1
bar_width = 1
gene_to_locations = defaultdict(list)
plotted_gene_names = set()

def plot_bpg_connection(ref_placements,path,total_length,prev_seg_index_is_adj,bpg_dict,seg_end_pos_d):
    connect_width = bar_width/2.
    for ind,refObj in ref_placements.iteritems():
        next_ind = (ind+1) % len(ref_placements)
        next_refObj = ref_placements[next_ind]
        if not prev_seg_index_is_adj[next_ind]: #or next_ind == 0 to try and close
            bpg_adjacency = vu.pair_is_edge(refObj.id, next_refObj.id,refObj.direction, next_refObj.direction, 
                bpg_dict, seg_end_pos_d)

            if not bpg_adjacency:
                continue

            # start_angle, end_angle = start_end_angle(next_refObj.abs_start_pos,refObj.abs_end_pos,total_length)
            bpg_connector_len = next_refObj.abs_start_pos - refObj.abs_end_pos
            #makes the reference genome wedges
            # patches.append(mpatches.Wedge((0,0), seg_bar_base - bar_width/4, end_angle, start_angle, width=connect_width))
            #draw the connecting triangle
            patches.append(mpatches.Rectangle((refObj.abs_end_pos,seg_bar_base + bar_width/4.),bpg_connector_len,connect_width))
            f_color_v.append('grey')
            e_color_v.append('grey')
            lw_v.append(0.2)

#must return positions of transcribed regions of genes here
def plot_gene_track(currStart, relGenes, pTup, total_length, strand):
    for ind,i in enumerate(relGenes):
        truncStart = False
        truncEnd = False
        #e_posns is a list of tuples of exon (start,end)
        #these can be plotted similarly to how the coding region is marked
        tstart,tend,e_posns = relGenes[i]
        seg_len = pTup[2] - pTup[1]
        if strand == "+":
            normStart = currStart + max(0,tstart-pTup[1])
            normEnd = currStart + min(seg_len,tend-pTup[1])
        else:
            normEnd = currStart + min(seg_len,pTup[2]-tstart)
            normStart = currStart + max(0,pTup[2] - tend)

        # print max(0,tstart-pTup[1]),min(seg_len,tend-pTup[1]),seg_len
        # print i,normStart,normEnd, currStart, currStart+seg_len,tstart,tend,strand

        gene_to_locations[i].append((normStart,normEnd))
        box_len = normEnd - normStart
        # patches.append(mpatches.Wedge((0,0), seg_bar_base, start_angle, end_angle, width=bar_width/2.0))
        patches.append(mpatches.Rectangle((normStart,gene_bar_height + bar_width),box_len,0.6*bar_width))
        f_color_v.append('k')
        e_color_v.append('k')
        lw_v.append(0)

        if i not in plotted_gene_names:
            # x,y = pol2cart(seg_bar_base+(bar_width/2.0),(text_angle/360*2*np.pi))
            # x_t,y_t = pol2cart(seg_bar_base + bar_width + 0.9,(text_angle/360*2*np.pi))
            #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.4)
            
            # text_angle,ha = vu.correct_text_angle(text_angle)
            # ax.text(x_t,y_t,i,color='grey',rotation=text_angle,ha=ha,fontsize=6.5,rotation_mode='anchor')
            ax.text(normStart + box_len/2.,gene_bar_height + 0.1*bar_width,i,color='grey',ha="center",fontsize=6)
            plotted_gene_names.add(i)

        # for exon in e_posns:
        #     if exon[1] > pTup[1] and exon[0] < pTup[2]:
        #         if strand == "+":
        #             normStart = currStart + max(1,exon[0]-pTup[1])
        #             normEnd = currStart + min(pTup[2]-pTup[1],exon[1]-pTup[1])

        #         else:
        #             normEnd = currStart + min(pTup[2]-pTup[1],pTup[2]-exon[0])
        #             normStart = currStart + max(1,pTup[2] - exon[1])

        #         # start_angle, end_angle = start_end_angle(normStart,normEnd,total_length)
        #         # patches.append(mpatches.Wedge((0,0), seg_bar_base-bar_width/2.0, start_angle, end_angle, width=bar_width/2.0))
        #         patches.append(mpatches.Rectangle((normStart,gene_bar_height + 0.4*bar_width),box_len,0.6*bar_width))
        #         f_color_v.append('r')
        #         e_color_v.append('r')
        #         lw_v.append(0)

#plot the reference genome
def plot_ref_genome(ref_placements,path,total_length,segSeqD,imputed_status,label_segs,onco_set=set()):
    font0 = FontProperties()
    p_end = 0
    for ind,refObj in ref_placements.iteritems():
        seg_coord_tup = segSeqD[path[ind][0]]
        # print(refObj.to_string())
        # start_angle, end_angle = start_end_angle(refObj.abs_end_pos,refObj.abs_start_pos,total_length)
        box_len = refObj.abs_end_pos - refObj.abs_start_pos
        # print start_angle,end_angle
        
        #makes the reference genome wedges    
        # patches.append(mpatches.Wedge((0,0), seg_bar_base, end_angle, start_angle, width=bar_width))
        patches.append(mpatches.Rectangle((refObj.abs_start_pos,ref_bar_height),box_len,bar_width))
        chrom = segSeqD[path[ind][0]][0]
        f_color_v.append(chromosome_colors[chrom])
        e_color_v.append('grey')
        lw_v.append(0.2)
        
        #makes the ticks on the reference genome wedges
        if path[ind][1] == "+":
            # posns = zip(range(seg_coord_tup[1],seg_coord_tup[2]+1),np.arange(refObj.abs_end_pos,refObj.abs_start_pos,-1))
            posns = zip(range(seg_coord_tup[1],seg_coord_tup[2]+1),np.arange(refObj.abs_start_pos,refObj.abs_end_pos))
        else:
            # posns = zip(np.arange(seg_coord_tup[2],seg_coord_tup[1]-1,-1),np.arange(refObj.abs_end_pos,refObj.abs_start_pos,-1))
            posns = zip(np.arange(seg_coord_tup[2],seg_coord_tup[1]-1,-1),np.arange(refObj.abs_start_pos,refObj.abs_end_pos))

        tick_freq = max(20000,40000*int(np.floor(total_length/1000000)))        
        #segment too small, nothing gets ticked 
        if (not any(j[0] % tick_freq == 0 for j in posns)) and abs(refObj.abs_start_pos - p_end) > 1:
            tens = [j[0] for j in posns if j[0] % 10000 == 0]
            middleIndex = (len(tens) - 1)/2
            tick_freq = tens[middleIndex]

        for j in posns:
            if j[0] % tick_freq == 0:
                x_i,y_i = j[1],ref_bar_height
                x_f,y_f = j[1],ref_bar_height-bar_width*0.3
                ax.plot([x_i,x_f],[y_i,y_f],color='grey',linewidth=0.5)
                txt = " " + str(int(round((j[0])/10000))) # if ha == "left" else str(int(round((j[0])/10000))) + " "
                # txt = str(j[0])
                x_t,y_t = j[1],ref_bar_height-bar_width*0.4
                ax.text(x_t,y_t,txt,color='grey',rotation=-90,rotation_mode="anchor",
                        ha="left",va="center",fontsize=6)

        p_end = refObj.abs_end_pos    
        gene_tree = vu.parse_genes(seg_coord_tup[0])
        relGenes = vu.rel_genes(gene_tree,seg_coord_tup,onco_set)
        # plot the gene track
        # TODO: IMPLEMENT
        plot_gene_track(refObj.abs_start_pos,relGenes,seg_coord_tup,total_length,path[ind][1])

        # #label the segments by number in path
        # mid_sp = (refObj.abs_end_pos + refObj.abs_start_pos)/2
        # text_angle = mid_sp/total_length*360.
        # x,y = pol2cart((seg_bar_base-2*bar_width),(text_angle/360.*2.*np.pi))
        # font = font0.copy()
        # if imputed_status[ind]:
        #     font.set_style('italic')
        #     # font.set_weight('bold')

        # text_angle,ha = vu.correct_text_angle(text_angle)

        # if label_segs:
        #     ax.text(x,y,cycle[ind][0]+cycle[ind][1],color='grey',rotation=text_angle,
        #         ha=ha,fontsize=5,fontproperties=font,rotation_mode='anchor')

#plot cmap track
def plot_cmap_track(seg_placements,total_length,unadj_bar_height,color,seg_id_labels = False):
    path_label_locs = defaultdict(list)
    for ind,segObj in seg_placements.iteritems():
        bar_height = unadj_bar_height + segObj.track_height_shift
        print "cmap_plot",segObj.id
        print "cmap plotting abs end pos are"
        print segObj.abs_start_pos, segObj.abs_end_pos
        box_len = segObj.abs_end_pos - segObj.abs_start_pos

        #Draw the box
        patches.append(mpatches.Rectangle((segObj.abs_start_pos,bar_height),box_len,bar_width))
        f_color_v.append(color)
        e_color_v.append('k')
        lw_v.append(0)

        linewidth = min(0.25*2000000/total_length,0.25)
        #Draw the labels in the box
        for i in segObj.label_posns:
            if i > segObj.abs_end_pos or i < segObj.abs_start_pos:
                continue

            y_i,y_f = bar_height,bar_height+bar_width
            ax.plot([i,i],[y_i,y_f],color='k',alpha=0.9,linewidth=linewidth)

        #TODO: fix for dense packing
        if seg_id_labels:
            mid_sp = (segObj.abs_end_pos + segObj.abs_start_pos)/2
            text = segObj.id + segObj.direction 
            ax.text(mid_sp,bar_height+1.1*bar_width,text,color='grey',fontsize=5,ha="center")

    return path_label_locs

#plot the connecting lines for the bionano track
def plot_alignment(contig_locs,segment_locs,total_length):
    linewidth = min(0.25*2000000/total_length,0.25)
    print "linewidth",linewidth,total_length
    for a_d in aln_vect:
        c_id = a_d["contig_id"]
        c_num_dir = int(a_d["contig_dir"]+"1")
        contig_label_vect = contig_locs[c_id].label_posns
        seg_label_vect = segment_locs[a_d["seg_aln_number"]].label_posns
        clx = contig_label_vect[a_d["contig_label"]-1]
        slx = seg_label_vect[a_d["seg_label"]-1]
        # contig_top = seg_bar_base + contig_bar_height + contig_locs[c_id].track_height_shift + bar_width
        contig_bottom = seg_bar_base + contig_bar_height + contig_locs[c_id].track_height_shift
        ax.plot([slx, clx], [seg_bar_base + bar_width, contig_bottom], color="grey",linewidth=linewidth)

def construct_path_ref_placements(path,segSeqD,raw_path_length,prev_seg_index_is_adj,aln_vect = []):
    spacing_bp = seg_spacing*raw_path_length
    path_ref_placements = {}
    curr_start = 0.0
    for ind,i in enumerate(path):
        seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        seg_end = curr_start+seg_len
        curr_obj = vu.CycleVizElemObj(i[0],i[1],curr_start,seg_end)
        path_ref_placements[ind] = curr_obj
        next_start = seg_end
        mod_ind = (ind+1) % (len(prev_seg_index_is_adj))
        if not prev_seg_index_is_adj[mod_ind]:
            next_start+=spacing_bp

        curr_start = next_start

    total_length = next_start
    return path_ref_placements,total_length

parser = argparse.ArgumentParser(description="CLinear visualizations of AA & AR output")
parser.add_argument("--om_alignments",help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
    action='store_true')
parser.add_argument("--cycles_file",help="AA/AR cycles-formatted input file",required=True)
parser.add_argument("--path",help="path number to visualize",required=True)
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("-g", "--graph", help="breakpoint graph file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--label_segs",help="label segs with graph IDs",action='store_true')
parser.add_argument("--gene_subset_file",help="File containing subset of genes to plot (e.g. oncogene genelist file)")
parser.add_argument("--reduce_path",help="Number of path elements to remove from left and right ends. Must supply both values, \
                    default 0 0",nargs=2,type=int,default=[0,0])

args = parser.parse_args()

if not args.sname:
    args.sname = os.path.split(args.cycles_file)[1].split(".")[0] + "_"

outdir = os.path.dirname(args.sname)
if outdir and not os.path.exists(outdir):
    os.makedirs(outdir)

fname = args.sname + "_path_" + args.path

print args.reduce_path

print("Unaligned fraction cutoff set to " + str(vu.unaligned_cutoff_frac))

chromosome_colors = vu.get_chr_colors()
plt.clf()
fig, ax = plt.subplots()
patches = []
f_color_v = []
e_color_v = []
lw_v = []

paths,segSeqD,circular_D = vu.parse_cycles_file(args.cycles_file)
path_num = args.path
path = paths[path_num]

if args.reduce_path != [0,0]:
    isCycle = False
else:
    isCycle = circular_D[path_num]

#debug
# vu.reduce_path(path,)
# print "PATH:",path
# print "SEGSEQD",segSeqD

prev_seg_index_is_adj = vu.adjacent_segs(path,segSeqD,isCycle)
print path
print "PSIIA",prev_seg_index_is_adj
raw_path_length = vu.get_raw_path_length(path,segSeqD,False,prev_seg_index_is_adj)

bpg_dict,seg_end_pos_d = {},{}
if args.graph:
    bpg_dict,seg_end_pos_d = vu.parse_BPG(args.graph)

gene_set = set()
if args.gene_subset_file:
    gene_set = vu.parse_gene_subset_file(args.gene_subset_file)

if not args.om_alignments:
    ref_placements,total_length = construct_path_ref_placements(cycle,segSeqD,raw_path_length,prev_seg_index_is_adj)
    if args.reduce_path != [0,0]:
        #reduce alignments
        path,prev_seg_index_is_adj,_ = vu.reduce_path(path,prev_seg_index_is_adj,args.reduce_path)

else:
    seg_cmaps = parse_cmap(args.segs,True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    seg_cmap_lens = get_cmap_lens(args.segs)
    aln_vect,meta_dict = vu.parse_alnfile(args.path_alignment)
    if args.reduce_path != [0,0]:
        #reduce alignments
        path,prev_seg_index_is_adj,aln_vect = vu.reduce_path(path,prev_seg_index_is_adj,args.reduce_path,aln_vect)

    is_segdup,split_ind = vu.check_segdup(aln_vect,path,isCycle)
    if is_segdup:
        print("alignment shows simple segdup")
        path = [path[0]]*2
        print path
        isCycle = False
        prev_seg_index_is_adj = [False,True]
        for a_ind in range(split_ind,len(aln_vect)):
            aln_vect[a_ind]["seg_aln_number"] = 1

    ref_placements,total_length = construct_path_ref_placements(path,segSeqD,raw_path_length,
                                                                prev_seg_index_is_adj,aln_vect)
    bar_width = total_length*bar_width_scaling
    contig_bar_height+=bar_width*bar_drop_prop
    path_seg_placements = vu.place_path_segs_and_labels(path,ref_placements,seg_cmap_vects)

    contig_cmaps = parse_cmap(args.contigs,True)
    contig_cmap_vects = vectorize_cmaps(contig_cmaps)
 
    ###
    #TODO: TRIM REF SEGS
    ###

    contig_cmap_lens = get_cmap_lens(args.contigs)
    #path_seg_placements,aln_vect,total_length,contig_cmap_vects
    contig_placements,contig_list = vu.place_contigs_and_labels(path_seg_placements,aln_vect,total_length,
                                                                contig_cmap_vects,isCycle,True)
    vu.decide_trim_contigs(contig_cmap_vects,contig_placements,total_length)

    #plot segs cmap
    print "SH",seg_bar_base+segment_bar_height
    print "CH",seg_bar_base+contig_bar_height
    plot_cmap_track(path_seg_placements,total_length,seg_bar_base+segment_bar_height,"darkorange")

    # check overlaps of contigs and adjust heights accordingly
    contig_height_shifts = vu.set_contig_height_shifts(contig_placements,contig_list,-bar_width)
    #plot contigs cmap
    plot_cmap_track(contig_placements,total_length,seg_bar_base+contig_bar_height,"cornflowerblue",seg_id_labels=True)

    #plot alignments
    plot_alignment(contig_placements,path_seg_placements,total_length)

imputed_status = vu.imputed_status_from_aln(aln_vect,len(path))

ref_bar_height=seg_bar_base
gene_bar_height=seg_bar_base
gene_bar_height-=bar_width*bar_drop_prop
ref_bar_height-=(bar_width*1.5*bar_drop_prop)
print "RH",seg_bar_base-ref_bar_height
plot_ref_genome(ref_placements,path,total_length,segSeqD,imputed_status,args.label_segs,gene_set)

if args.graph:
    plot_bpg_connection(ref_placements,path,total_length,prev_seg_index_is_adj,bpg_dict,seg_end_pos_d)

# ax.set_xlim(-(seg_bar_base+1.25), (seg_bar_base+1.25))
# ax.set_ylim(-(seg_bar_base+1.25), (seg_bar_base+1.25))
chrom_set = set()
for i in path:
    chrom_set.add(segSeqD[i[0]][0])

sorted_chrom = sorted(chrom_set,key=lambda x:x.rsplit("chr")[-1])
sorted_chrom_colors = [chromosome_colors[x] for x in sorted_chrom]
legend_patches = []
for chrom,color in zip(sorted_chrom,sorted_chrom_colors):
    legend_patches.append(mpatches.Patch(color=color,label=chrom))

# plt.legend(handles=legend_patches,fontsize=8,loc=3,bbox_to_anchor=(-.3,.15))
plt.legend(handles=legend_patches,fontsize=8,bbox_to_anchor=(.09,-1.5))

p = PatchCollection(patches)
p.set_facecolor(f_color_v)
p.set_edgecolor(e_color_v)
p.set_linewidth(lw_v)
ax.add_collection(p)
ax.set_aspect(1.0)
plt.axis('off')
plt.savefig(fname + '.png',dpi=600)
plt.savefig(fname + '.pdf',format='pdf')

plt.close()
print "finished"
