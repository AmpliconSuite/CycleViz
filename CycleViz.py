#!/usr/bin/env python

import sys
import os
import copy
import bisect
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
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

# try:
#     sys.path.insert(0,os.environ['AR_SRC'])

# except KeyError:
#     sys.stderr.write("AmpliconReconstructor source directory bash variable ($AR_SRC) not found.\n")
#     sys.stderr.write("Is AmpliconReconstructor configured?")
#     sys.exit()

from bionanoUtil import *
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

seg_spacing = 0.009
bar_width = 2.5/3
global_rot = 90.0
center_hole = 0.75
outer_bar = 10 
bed_spacing = .5
contig_bar_height = -14/3
segment_bar_height = -8.0/3
gene_to_locations = defaultdict(list)

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x) / (2. * np.pi) * 360 
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi) 
    y = rho * np.sin(phi)
    return(x, y)

def polar_series_to_cartesians(line_points,r):
    #homer brain
    # x_list = []
    # y_list = []
    # for i in line_points:
    #     x,y = pol2cart(r,i)
    #     x_list.append(x)
    #     y_list.append(y)

    # return x_list,y_list

    #better
    # xy_tups = [pol2cart(r,i) for i in line_points]
    # return zip(*xy_tups) 

    #"galaxy brain"
    return zip(*[pol2cart(r,i) for i in line_points])


def start_end_angle(normStart,normEnd,total_length):
    #print "init norm_e_s",normStart,normEnd
    start_angle = normStart/total_length*360
    end_angle = normEnd/total_length*360
    # print start_angle,end_angle
    # text_angle = (start_angle + end_angle)/2.0
    #print "init s_e_angle",start_angle,end_angle
    if end_angle < 0 and start_angle > 0:
        end_angle+=360
    ##    start_angle,end_angle = end_angle,start_angle

    #handle circularl contig
    if start_angle > 360 and start_angle % 360 > end_angle:
        start_angle,end_angle = 360,0
    
    #print "s_e_angle",start_angle,end_angle
    return start_angle,end_angle


def plot_bpg_connection(ref_placements,cycle,total_length,prev_seg_index_is_adj,bpg_dict,seg_end_pos_d):
    connect_width = bar_width/2.
    for ind,refObj in ref_placements.iteritems():
        next_ind = (ind+1) % len(ref_placements)
        next_refObj = ref_placements[next_ind]
        if not prev_seg_index_is_adj[next_ind]: #or next_ind == 0 to try and close
            bpg_adjacency = vu.pair_is_edge(refObj.id, next_refObj.id,refObj.direction, next_refObj.direction, 
                bpg_dict, seg_end_pos_d)

            if not bpg_adjacency:
                continue

            start_angle, end_angle = start_end_angle(next_refObj.abs_start_pos,refObj.abs_end_pos,total_length)
            #makes the reference genome wedges
            patches.append(mpatches.Wedge((0,0), outer_bar - bar_width/4, end_angle, start_angle, width=connect_width))
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

        start_angle = normStart/total_length*360
        end_angle = normEnd/total_length*360
        text_angle = (start_angle + end_angle)/2.0
        gene_to_locations[i].append((start_angle/360.,end_angle/360.))
        if end_angle < 0 and start_angle > 0:
            end_angle+=360
        
        patches.append(mpatches.Wedge((0,0), outer_bar, start_angle, end_angle, width=bar_width/2.0))
        f_color_v.append('k')
        e_color_v.append('k')
        lw_v.append(0)

        # x,y = pol2cart(outer_bar+(bar_width/2.0),(text_angle/360*2*np.pi))
        x_t,y_t = pol2cart(outer_bar + bar_width + 2,(text_angle/360*2*np.pi))
        #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.4)
        
        text_angle,ha = vu.correct_text_angle(text_angle)
        ax.text(x_t,y_t,i,color='k',rotation=text_angle,ha=ha,va="center",fontsize=9,rotation_mode='anchor')

        for exon in e_posns:
            if exon[1] > pTup[1] and exon[0] < pTup[2]:
                if strand == "+":
                    normStart = currStart + max(1,exon[0]-pTup[1])
                    normEnd = currStart + min(pTup[2]-pTup[1],exon[1]-pTup[1])

                else:
                    normEnd = currStart + min(pTup[2]-pTup[1],pTup[2]-exon[0])
                    normStart = currStart + max(1,pTup[2] - exon[1])

                start_angle, end_angle = start_end_angle(normStart,normEnd,total_length)
                patches.append(mpatches.Wedge((0,0), outer_bar-bar_width/2.0, start_angle, end_angle, width=bar_width/2.0))
                f_color_v.append('k')
                e_color_v.append('k')
                lw_v.append(0)

#plot the reference genome
def plot_ref_genome(ref_placements,cycle,total_length,segSeqD,imputed_status,label_segs,onco_set=set()):
    font0 = FontProperties()
    rot_sp = global_rot/360.*total_length
    for ind,refObj in ref_placements.iteritems():
        seg_coord_tup = segSeqD[cycle[ind][0]]
        # print(refObj.to_string())
        start_angle, end_angle = start_end_angle(refObj.abs_end_pos,refObj.abs_start_pos,total_length)
        # print start_angle,end_angle
        
        #makes the reference genome wedges    
        patches.append(mpatches.Wedge((0,0), outer_bar, end_angle, start_angle, width=bar_width))
        chrom = segSeqD[cycle[ind][0]][0]
        f_color_v.append(chromosome_colors[chrom])
        e_color_v.append('grey')
        lw_v.append(0.2)
        
        #makes the ticks on the reference genome wedges
        if cycle[ind][1] == "+":
            # posns = zip(range(seg_coord_tup[1],seg_coord_tup[2]+1),np.arange(refObj.abs_end_pos,refObj.abs_start_pos,-1))
            posns = zip(range(seg_coord_tup[1],seg_coord_tup[2]+1),np.arange(refObj.abs_start_pos,refObj.abs_end_pos))
        else:
            # posns = zip(np.arange(seg_coord_tup[2],seg_coord_tup[1]-1,-1),np.arange(refObj.abs_end_pos,refObj.abs_start_pos,-1))
            posns = zip(np.arange(seg_coord_tup[2],seg_coord_tup[1]-1,-1),np.arange(refObj.abs_start_pos,refObj.abs_end_pos))

        tick_freq = max(30000,30000*int(np.floor(total_length/1000000)))
        if refObj.abs_end_pos-refObj.abs_start_pos < 30000:
            tick_freq = 10000
        for j in posns:
            if j[0] % tick_freq == 0:
                text_angle = j[1]/total_length*360
                x,y = pol2cart(outer_bar,(text_angle/360*2*np.pi))
                x_t,y_t = pol2cart(outer_bar + 0.2,(text_angle/360*2*np.pi))
                ax.plot([x,x_t],[y,y_t],color='grey',linewidth=1)
                
                text_angle,ha = vu.correct_text_angle(text_angle)
                txt = " " + str(int(round((j[0])/10000))) if ha == "left" else str(int(round((j[0])/10000))) + " "

                ax.text(x_t,y_t,txt,color='grey',rotation=text_angle,
                ha=ha,va="center",fontsize=9,rotation_mode='anchor')
    
        gene_tree = vu.parse_genes(seg_coord_tup[0],args.ref)
        relGenes = vu.rel_genes(gene_tree,seg_coord_tup,copy.copy(onco_set))
        #plot the gene track
        plot_gene_track(refObj.abs_start_pos,relGenes,seg_coord_tup,total_length,cycle[ind][1])

        #label the segments by number in cycle
        mid_sp = (refObj.abs_end_pos + refObj.abs_start_pos)/2
        text_angle = mid_sp/total_length*360.
        x,y = pol2cart((outer_bar-2*bar_width),(text_angle/360.*2.*np.pi))
        font = font0.copy()
        if imputed_status[ind]:
            font.set_style('italic')
            # font.set_weight('bold')

        text_angle,ha = vu.correct_text_angle(text_angle)

        if label_segs:
            ax.text(x,y,cycle[ind][0]+cycle[ind][1],color='grey',rotation=text_angle,
                ha=ha,fontsize=5,fontproperties=font,rotation_mode='anchor')

#plot cmap track
def plot_cmap_track(seg_placements,total_length,unadj_bar_height,color,seg_id_labels = False):
    cycle_label_locs = defaultdict(list)
    for ind,segObj in seg_placements.iteritems():
        bar_height = unadj_bar_height + segObj.track_height_shift
        print "cmap_plot",segObj.id
        print "cmap plotting abs end pos are"
        print segObj.abs_end_pos,segObj.abs_start_pos
        start_angle, end_angle = start_end_angle(segObj.abs_end_pos,segObj.abs_start_pos,total_length)
        # print start_angle,end_angle
        patches.append(mpatches.Wedge((0,0), bar_height + bar_width, end_angle, start_angle, width=bar_width))
        f_color_v.append(color)
        e_color_v.append('k')
        lw_v.append(0)

        # print "linewidth alt",0.2*1000000/total_length
        linewidth = min(0.25*2000000/total_length,0.25)
        for i in segObj.label_posns:
            if i > segObj.abs_end_pos or i < segObj.abs_start_pos:
                continue

            label_rads = i/total_length*2*np.pi
            x,y = pol2cart(bar_height,label_rads)
            x_t,y_t = pol2cart(bar_height+bar_width,label_rads)
            # linewidth = min(0.2*2000000/total_length,0.2)
            ax.plot([x,x_t],[y,y_t],color='k',alpha=0.9,linewidth=linewidth)

        if seg_id_labels:
            mid_sp = (segObj.abs_end_pos + segObj.abs_start_pos)/2
            text_angle = mid_sp/total_length*360.
            x,y = pol2cart(bar_height-1.2,(text_angle/360.*2.*np.pi))
            text_angle,ha = vu.correct_text_angle(text_angle)
            text = segObj.id + segObj.direction 
            ax.text(x,y,text,color='grey',rotation=text_angle,
                    ha=ha,fontsize=5,rotation_mode='anchor')

    return cycle_label_locs

#plot the connecting lines for the bionano track
def plot_alignment(contig_locs,segment_locs,total_length):
    segs_base = outer_bar+segment_bar_height
    linewidth = min(0.25*2000000/total_length,0.25)
    print "linewidth",linewidth,total_length
    for a_d in aln_vect:
        c_id = a_d["contig_id"]
        c_num_dir = int(a_d["contig_dir"]+"1")
        # s_num_dir = int(a_d["seg_dir"]+"1")
        # s_num_dir_alt = int(segment_locs[a_d["seg_aln_number"]].direction+"1")

        contig_label_vect = contig_locs[c_id].label_posns
        seg_label_vect = segment_locs[a_d["seg_aln_number"]].label_posns
        # c_l_pos = contig_label_vect[a_d["contig_label"]-1]
        # c_l_pos = contig_label_vect[c_num_dir*a_d["contig_label"]-(c_num_dir+1)/2]
        # if s_num_dir == -1 and c_num_dir == -1:
        #     c_l_pos = contig_label_vect[-1*a_d["contig_label"]]
        # else:
        c_l_pos = contig_label_vect[a_d["contig_label"]-1]
        c_l_loc = c_l_pos/total_length*2.*np.pi
        # s_l_pos = seg_label_vect[s_num_dir*a_d["seg_label"]-(s_num_dir+1)/2]
        s_l_pos = seg_label_vect[a_d["seg_label"]-1]
        s_l_loc = s_l_pos/total_length*2.*np.pi
        contig_top = outer_bar + contig_bar_height + contig_locs[c_id].track_height_shift + bar_width
        x_c,y_c = pol2cart(contig_top,c_l_loc)
        x_s,y_s = pol2cart(segs_base,s_l_loc)
        ax.plot([x_c, x_s], [y_c, y_s], color="grey",linewidth=linewidth)

# def get_seg_locs_from_cycle(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj):
#     spacing_bp = seg_spacing*raw_cycle_length
#     start_points = [spacing_bp/2.]
#     print "start points"
#     for ind,i in enumerate(cycle[:-1]):
#         print ind,i
#         seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
#         # seg_len = seg_cmap_lens[i[0]]
#         if not prev_seg_index_is_adj[ind+1]:
#             print "adding spacing to ",cycle[ind+1]
#             seg_len+=spacing_bp

#         start_points.append(start_points[-1] + seg_len)
#         # print start_points[-1],seg_len

#     last_len = segSeqD[cycle[-1][0]][2] - segSeqD[cycle[-1][0]][1]
#     # last_len = lens[cycle[-1][0]]
#     endpoint = start_points[-1] + last_len
#     if not prev_seg_index_is_adj[0]:
#         endpoint+=spacing_bp

#     return start_points,endpoint

def construct_cycle_ref_placements(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj,isCycle,aln_vect = []):
    spacing_bp = seg_spacing*raw_cycle_length
    cycle_ref_placements = {}
    curr_start = 0.0 if isCycle else spacing_bp
    for ind,i in enumerate(cycle):
        seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        seg_end = curr_start+seg_len
        curr_obj = vu.CycleVizElemObj(i[0],i[1],curr_start,seg_end)
        cycle_ref_placements[ind] = curr_obj
        next_start = seg_end
        mod_ind = (ind+1) % (len(prev_seg_index_is_adj))
        if not prev_seg_index_is_adj[mod_ind]:
            next_start+=spacing_bp

        curr_start = next_start

    total_length = next_start
    return cycle_ref_placements,total_length


parser = argparse.ArgumentParser(description="Circular visualizations of AA & AR output")
parser.add_argument("--om_alignments",help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
    action='store_true')
parser.add_argument("--cycles_file",help="AA/AR cycles-formatted input file",required=True)
parser.add_argument("--cycle",help="cycle number to visualize",required=True)
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("-g", "--graph", help="breakpoint graph file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--ref",help="reference genome",choices=["hg19","hg38","GRCh37","GRCh38"],default="hg19")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--rot", help="number of segments to rotate counterclockwise",type=int,default=0)
parser.add_argument("--label_segs",help="label segs with graph IDs",action='store_true')
parser.add_argument("--gene_subset_file",help="File containing subset of genes to plot (e.g. oncogene genelist file)")

args = parser.parse_args()

if args.ref == "GRCh38":
    args.ref == "hg38"

print args.ref

if not args.sname:
    args.sname = os.path.split(args.cycles_file)[1].split(".")[0] + "_"

outdir = os.path.dirname(args.sname)
if outdir and not os.path.exists(outdir):
    os.makedirs(outdir)

fname = args.sname + "_cycle_" + args.cycle

print("Unaligned fraction cutoff set to " + str(vu.unaligned_cutoff_frac))
chromosome_colors = vu.get_chr_colors()
plt.clf()
fig, ax = plt.subplots()
patches = []
f_color_v = []
e_color_v = []
lw_v = []

cycles,segSeqD,circular_D = vu.parse_cycles_file(args.cycles_file)
cycle_num = args.cycle
isCycle = circular_D[cycle_num]
cycle = cycles[cycle_num]
prev_seg_index_is_adj = vu.adjacent_segs(cycle,segSeqD,isCycle)
raw_cycle_length = vu.get_raw_path_length(cycle,segSeqD,isCycle,prev_seg_index_is_adj)

bpg_dict,seg_end_pos_d = {},{}
if args.graph:
    bpg_dict,seg_end_pos_d = vu.parse_BPG(args.graph)

gene_set = set()
if args.gene_subset_file:
    gene_set = vu.parse_gene_subset_file(args.gene_subset_file)

# for i,k in seg_end_pos_d.iteritems():
#     print i,k

# print ""
# for i,k in bpg_dict.iteritems():
#     print i,k

# start_points,total_length = get_seg_locs_from_cycle(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj)
if not args.om_alignments:
    ref_placements,total_length = construct_cycle_ref_placements(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj,isCycle)
    imputed_status = [False]*len(cycle)

else:
    print("Visualizing with alignments")
    print("Contig spacing set to " + str(vu.contig_spacing))
    seg_cmaps = parse_cmap(args.segs,True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    seg_cmap_lens = get_cmap_lens(args.segs)
    aln_vect,meta_dict = vu.parse_alnfile(args.path_alignment)
    is_segdup,split_ind = vu.check_segdup(aln_vect,cycle,isCycle)
    if is_segdup:
        print("alignment shows simple segdup")
        cycle = [cycle[0]]*2
        print cycle
        isCycle = False
        prev_seg_index_is_adj = [False,True]
        for a_ind in range(split_ind,len(aln_vect)):
            aln_vect[a_ind]["seg_aln_number"] = 1

    ref_placements,total_length = construct_cycle_ref_placements(cycle,segSeqD,raw_cycle_length,
                                                                prev_seg_index_is_adj,isCycle,aln_vect)
    cycle_seg_placements = vu.place_path_segs_and_labels(cycle,ref_placements,seg_cmap_vects)

    contig_cmaps = parse_cmap(args.contigs,True)
    contig_cmap_vects = vectorize_cmaps(contig_cmaps)

    ###
    #TODO: TRIM REF SEGS
    ###

    contig_cmap_lens = get_cmap_lens(args.contigs)
    #cycle_seg_placements,aln_vect,total_length,contig_cmap_vects
    contig_placements,contig_list = vu.place_contigs_and_labels(cycle_seg_placements,aln_vect,total_length,
                                                                contig_cmap_vects,isCycle,True)
    vu.decide_trim_contigs(contig_cmap_vects,contig_placements,total_length)

    #plot segs

    #seg_height_shifts = {x:0 for x in seg_placements}
    #seg_placements,total_length,bar_height,color,seg_text = None
    plot_cmap_track(cycle_seg_placements,total_length,outer_bar+segment_bar_height,"darkorange")


    #plot contigs

    # check overlaps of contigs and adjust heights accordingly
    contig_height_shifts = vu.set_contig_height_shifts(contig_placements,contig_list)
    plot_cmap_track(contig_placements,total_length,outer_bar+contig_bar_height,"cornflowerblue",seg_id_labels=True)


    #plot alignments
    plot_alignment(contig_placements,cycle_seg_placements,total_length)

    imputed_status = vu.imputed_status_from_aln(aln_vect,len(cycle))

    # contig_lab_dict = {}
    # for c_tup,start_posn in contig_locs.iteritems():
    #     contig_cycle = [c_tup]
    #     contig_cmap = contig_cmaps[c_tup[0]]
    #     contig_vect = {c_tup[0]:[contig_cmap[i] for i in sorted(contig_cmap.keys())]}
    #     seg_text = c_tup[0]
    #     plot_cmap_track([start_posn],[0],total_length,contig_cycle,outer_bar+contig_bar_height,
    #         "cornflowerblue",contig_vect,seg_text)
    #     contig_lab_dict[c_tup[0]] = contig_label_posns[0]

    # #plot segs
    # aln_nums = range(0,len(cycle))
    # if len(cycle) > 1 and aln_vect[-1]["seg_aln_number"] == "0":
    #     aln_nums.append(0)


#plot_ref_genome(ref_placements,cycle,total_length,segSeqD,imputed_status,label_segs=True)
plot_ref_genome(ref_placements,cycle,total_length,segSeqD,imputed_status,args.label_segs,gene_set)

if args.graph:
    plot_bpg_connection(ref_placements,cycle,total_length,prev_seg_index_is_adj,bpg_dict,seg_end_pos_d)

ax.set_xlim(-(outer_bar+1.25), (outer_bar+1.25))
ax.set_ylim(-(outer_bar+1.25), (outer_bar+1.25))
chrom_set = set()
for i in cycle:
    chrom_set.add(segSeqD[i[0]][0])

sorted_chrom = sorted(chrom_set,key=lambda x:x.rsplit("chr")[-1])
sorted_chrom_colors = [chromosome_colors[x] for x in sorted_chrom]
legend_patches = []
for chrom,color in zip(sorted_chrom,sorted_chrom_colors):
    legend_patches.append(mpatches.Patch(color=color,label=chrom))

plt.legend(handles=legend_patches,fontsize=8,loc=3,bbox_to_anchor=(-.3,.15))

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
