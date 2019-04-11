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
import matplotlib._color_data as mcd
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
from intervaltree import Interval, IntervalTree

#non-interactive backend

try:
    sys.path.insert(0,os.environ['AR_SRC'])

except KeyError:
    sys.stderr.write("AmpliconReconstructor source directory bash variable ($AR_SRC) not found.\n")
    sys.stderr.write("Is AmpliconReconstructor configured?")
    sys.exit()

from bionanoUtil import *

contig_spacing = 0.015
seg_spacing = 0.009
bar_width = 2.5/3
global_rot = 90.0
center_hole = 0.75
outer_bar = 10 
bed_spacing = .5
contig_bar_height = -14/3
segment_bar_height = -8.0/3
unaligned_cutoff_frac = 1./12

gene_to_locations = defaultdict(list)

def round_to_1_sig(x):
    if x == 0:
        return 0.

    return round(x, -int(np.floor(np.log10(abs(x)))))

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x) / (2. * np.pi) * 360 
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi) 
    y = rho * np.sin(phi)
    return(x, y)

def polar_series_to_cartesians(line_points,r):
    x_list = []
    y_list = []
    for i in line_points:
        x,y = pol2cart(r,i)
        x_list.append(x)
        y_list.append(y)
        
    return x_list,y_list
    
def start_end_angle(normStart,normEnd,total_length_with_spacing):
    start_angle = normStart/total_length_with_spacing*360
    end_angle = normEnd/total_length_with_spacing*360
    # print start_angle,end_angle
    # text_angle = (start_angle + end_angle)/2.0
    if end_angle < 0 and start_angle > 0:
        end_angle+=360

    return start_angle,end_angle

class CycleVizElemObj(object):
    def __init__(self,m_id,direction,s,t,cmap_vect = []):
        self.id = m_id
        self.direction = direction
        self.abs_start_pos = s
        self.abs_end_pos = t
        self.scaling_factor = 1
        self.cmap_vect = cmap_vect
        self.aln_lab_ends = (None,None)
        self.aln_bound_posns = (None,None)
        self.label_posns = []
        self.track_height_shift = 0
        self.start_trim = False
        self.end_trim = False


    def compute_label_posns(self):
        if self.direction == "+":
            if self.abs_start_pos == None:
                self.abs_start_pos = self.aln_bound_posns[0] - self.scaling_factor*self.cmap_vect[self.aln_lab_ends[0]-1]
            
            if self.abs_end_pos == None:
                self.abs_end_pos = self.aln_bound_posns[1] + self.scaling_factor*(self.cmap_vect[-1] - self.cmap_vect[self.aln_lab_ends[1]-1])

            for i in self.cmap_vect[:-1]:
                self.label_posns.append(self.scaling_factor*i + self.abs_start_pos)

        else:
            if self.abs_start_pos == None:
                self.abs_start_pos = self.aln_bound_posns[0] - self.scaling_factor*(self.cmap_vect[-1] - self.cmap_vect[self.aln_lab_ends[1]-1])
            
            if self.abs_end_pos == None:
                self.abs_end_pos = self.aln_bound_posns[1] + self.scaling_factor*self.cmap_vect[self.aln_lab_ends[0]-1]

            rev_cmap_vect = [self.cmap_vect[-1] - x for x in self.cmap_vect[::-1][1:]] #does not include length value
            for i in rev_cmap_vect:
                self.label_posns.append(self.scaling_factor*i + self.abs_start_pos)

            self.label_posns = self.label_posns[::-1]

    def to_string(self):
        return "{}{} | Start: {} | End: {} | scaling {}".format(self.id,self.direction,str(self.abs_start_pos),str(self.abs_end_pos),str(self.scaling_factor))

def parse_genes(chrom):
    print("Building interval tree for chrom " + chrom)
    t = IntervalTree()
    with open(os.environ['AR_SRC'] + "/hg19_refGene.txt") as infile:
        for line in infile:
            fields = line.rsplit("\t")
            #gene = fields[-2]
            currChrom = fields[2]
            tstart = int(fields[4])
            tend = int(fields[5])
            if chrom == currChrom:
                t[tstart:tend] = fields

    return t

#rotate text to be legible on both sides of circle
def correct_text_angle(text_angle):
    if (abs(text_angle) > 90 and abs(text_angle) < 270):
        text_angle-=180
        ha = "right"
    else:
        ha = "left"

    return text_angle,ha

def rel_genes(chrIntTree,pTup):    
    relGenes = {}
    chrom = pTup[0]
    overlappingT = chrIntTree[pTup[1]:pTup[2]]
    for i in overlappingT:
        gene = i.data[-4].rsplit("-")[0]
        tstart = int(i.data[4])
        tend = int(i.data[5])
        e_s_posns = [int(x) for x in i.data[9].rsplit(",") if x]
        e_e_posns = [int(x) for x in i.data[10].rsplit(",") if x]
        if not gene.startswith("LOC"):
            if gene not in relGenes:
                relGenes[gene] = (tstart,tend,zip(e_s_posns,e_e_posns))

            else:
                oldTStart = relGenes[gene][0]
                oldTEnd = relGenes[gene][1]
                if tstart < oldTStart:
                    oldTStart = tstart
                if tend > oldTEnd:
                    oldTEnd = tend
                
                relGenes[gene] = (oldTStart,oldTEnd,relGenes[gene][2] + zip(e_s_posns,e_e_posns))
    
    return relGenes

#must return positions of transcribed regions of genes here
def plot_gene_track(currStart, relGenes, pTup, total_length_with_spacing, strand):
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

        start_angle = normStart/total_length_with_spacing*360
        end_angle = normEnd/total_length_with_spacing*360
        text_angle = (start_angle + end_angle)/2.0
        gene_to_locations[i].append((start_angle/360.,end_angle/360.))
        if end_angle < 0 and start_angle > 0:
            end_angle+=360
        
        patches.append(mpatches.Wedge((0,0), outer_bar, start_angle, end_angle, width=bar_width/2.0))
        f_color_v.append('k')
        e_color_v.append('k')
        lw_v.append(0)

        # x,y = pol2cart(outer_bar+(bar_width/2.0),(text_angle/360*2*np.pi))
        x_t,y_t = pol2cart(outer_bar + bar_width + 0.9,(text_angle/360*2*np.pi))
        #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.4)
        
        text_angle,ha = correct_text_angle(text_angle)
        ax.text(x_t,y_t,i,color='grey',rotation=text_angle,ha=ha,fontsize=6.5,rotation_mode='anchor')

        for exon in e_posns:
            if exon[1] > pTup[1] and exon[0] < pTup[2]:
                if strand == "+":
                    normStart = currStart + max(1,exon[0]-pTup[1])
                    normEnd = currStart + min(pTup[2]-pTup[1],exon[1]-pTup[1])

                else:
                    normEnd = currStart + min(pTup[2]-pTup[1],pTup[2]-exon[0])
                    normStart = currStart + max(1,pTup[2] - exon[1])

                start_angle, end_angle = start_end_angle(normStart,normEnd,total_length_with_spacing)
                patches.append(mpatches.Wedge((0,0), outer_bar-bar_width/2.0, start_angle, end_angle, width=bar_width/2.0))
                f_color_v.append('k')
                e_color_v.append('k')
                lw_v.append(0)

#plot the reference genome
def plot_ref_genome(ref_placements,cycle,total_length,segSeqD,imputed_status,label_segs=True):
    font0 = FontProperties()
    rot_sp = global_rot/360.*total_length
    for ind,refObj in ref_placements.iteritems():
        seg_coord_tup = segSeqD[cycle[ind][0]]
        print(refObj.to_string())
        start_angle, end_angle = start_end_angle(refObj.abs_end_pos,refObj.abs_start_pos,total_length)
        print start_angle,end_angle
        
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

        tick_freq = max(20000,20000*int(np.floor(total_length/1000000)))
        for j in posns:
            if j[0] % tick_freq == 0:
                text_angle = j[1]/total_length*360
                x,y = pol2cart(outer_bar,(text_angle/360*2*np.pi))
                x_t,y_t = pol2cart(outer_bar + 0.2,(text_angle/360*2*np.pi))
                ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.2)
                
                text_angle,ha = correct_text_angle(text_angle)
                txt = " " + str(int(round((j[0])/10000))) if ha == "left" else str(int(round((j[0])/10000))) + " "

                ax.text(x_t,y_t,txt,color='grey',rotation=text_angle,
                ha=ha,fontsize=6,rotation_mode='anchor')
    
        gene_tree = parse_genes(seg_coord_tup[0])
        relGenes = rel_genes(gene_tree,seg_coord_tup)
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

        text_angle,ha = correct_text_angle(text_angle)

        if label_segs:
            ax.text(x,y,cycle[ind][0]+cycle[ind][1],color='grey',rotation=text_angle,
                ha=ha,fontsize=5,fontproperties=font,rotation_mode='anchor')

#plot cmap track
def plot_cmap_track(seg_placements,total_length,unadj_bar_height,color,seg_id_labels = False):
    cycle_label_locs = defaultdict(list)
    for ind,segObj in seg_placements.iteritems():
        bar_height = unadj_bar_height + segObj.track_height_shift
        print "cmap_plot"
        print segObj.id
        print segObj.abs_end_pos,segObj.abs_start_pos
        start_angle, end_angle = start_end_angle(segObj.abs_end_pos,segObj.abs_start_pos,total_length)
        print start_angle,end_angle
        patches.append(mpatches.Wedge((0,0), bar_height + bar_width, end_angle, start_angle, width=bar_width))
        f_color_v.append(color)
        e_color_v.append('k')
        lw_v.append(0)

        for i in segObj.label_posns:
            label_rads = i/total_length*2*np.pi
            x,y = pol2cart(bar_height,label_rads)
            x_t,y_t = pol2cart(bar_height+bar_width,label_rads)
            ax.plot([x,x_t],[y,y_t],color='k',alpha=0.9,linewidth=0.2)

        if seg_id_labels:
            mid_sp = (segObj.abs_end_pos + segObj.abs_start_pos)/2
            text_angle = mid_sp/total_length*360.
            x,y = pol2cart(bar_height-1.2,(text_angle/360.*2.*np.pi))
            text_angle,ha = correct_text_angle(text_angle)
            text = segObj.id + segObj.direction 
            ax.text(x,y,text,color='grey',rotation=text_angle,
                    ha=ha,fontsize=5,rotation_mode='anchor')

    return cycle_label_locs

#plot the connecting lines for the bionano track
def plot_alignment(contig_locs,segment_locs,total_length):
    segs_base = outer_bar+segment_bar_height
    for a_d in aln_vect:
        contig_label_vect = contig_locs[a_d["contig_id"]].label_posns
        seg_label_vect = segment_locs[a_d["seg_aln_number"]].label_posns
        c_l_loc = contig_label_vect[a_d["contig_label"]-1]/total_length*2.*np.pi
        s_l_loc = seg_label_vect[a_d["seg_label"]-1]/total_length*2.*np.pi
        contig_top = outer_bar + contig_bar_height + contig_locs[a_d["contig_id"]].track_height_shift + bar_width
        x_c,y_c = pol2cart(contig_top,c_l_loc)
        x_s,y_s = pol2cart(segs_base,s_l_loc)
        ax.plot([x_c, x_s], [y_c, y_s], color="grey",linewidth=0.2)

def parse_cycles_file(cycles_file):
    cycles = {}
    segSeqD = {}
    circular_D = {}
    with open(cycles_file) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().split()
                lowerBound = int(fields[3])
                upperBound = int(fields[4])
                chrom = fields[2]
                segNum = fields[1]
                segSeqD[segNum] = (chrom,lowerBound,upperBound)

            elif "Cycle=" in line:
                isCircular = True
                curr_cycle = []
                fields = line.rstrip().rsplit(";")
                lineD = {x.rsplit("=")[0]:x.rsplit("=")[1] for x in fields}
                segs = lineD["Segments"].rsplit(",")
                for i in segs:
                    seg = i[:-1]
                    if seg != "0" and i:
                        strand = i[-1]
                        curr_cycle.append((seg,strand))

                    else:
                        isCircular = False

                cycles[lineD["Cycle"]] = curr_cycle
                circular_D[lineD["Cycle"]] = isCircular

    return cycles,segSeqD,circular_D

def parse_alnfile(path_aln_file):
    aln_vect = []
    with open(path_aln_file) as infile:
        meta_header = infile.next().rstrip()[1:].split()
        aln_metadata_fields = infile.next().rstrip()[1:].split()
        meta_dict = dict(zip(meta_header,aln_metadata_fields))
        aln_header = infile.next().rstrip()[1:].split()
        for line in infile:
            fields = line.rstrip().split()
            fields_dict = dict(zip(aln_header,fields))
            fields_dict["contig_label"] = int(fields_dict["contig_label"])
            fields_dict["seg_label"] = int(fields_dict["seg_label"])
            fields_dict["seg_aln_number"] = int(fields_dict["seg_aln_number"])
            aln_vect.append(fields_dict)

    return aln_vect,meta_dict

#determine segments linearly adjacent in ref genome
def adjacent_segs(cycle,segSeqD,isCircular):
    print "checking adjacency"
    print cycle
    prev_seg_index_is_adj = [False]*len(cycle)
    p_end = segSeqD[cycle[0][0]][2] if cycle[0][1] == "+" else segSeqD[cycle[0][0]][1]
    p_chrom = segSeqD[cycle[0][0]][0]
    for ind in range(1,len(cycle)):
        i = cycle[ind]
        curr_chrom = segSeqD[i[0]][0]
        curr_start = segSeqD[i[0]][2] if i[1] == "-" else segSeqD[i[0]][1]
        if curr_chrom == p_chrom and abs(curr_start - p_end) == 1:
            prev_seg_index_is_adj[ind] = True

        p_end = segSeqD[i[0]][2] if i[1] == "+" else segSeqD[i[0]][1]
        p_chrom = curr_chrom


    if isCircular and len(cycle) > 1:
        init_start = segSeqD[cycle[0][0]][2] if cycle[0][1] == "-" else segSeqD[cycle[0][0]][1]
        init_chr = segSeqD[cycle[0][0]][0]
        if p_chrom == curr_chrom and abs(init_start - p_end) == 1:
            prev_seg_index_is_adj[0] = True

    print prev_seg_index_is_adj
    return prev_seg_index_is_adj

def get_raw_cycle_length(cycle,segSeqD,isCircular,prev_seg_index_is_adj):
    raw_cycle_length= 0.0
    for i in cycle:
        s_tup = segSeqD[i[0]]
        s_len = s_tup[2] - s_tup[1]
        raw_cycle_length+=s_len

    return raw_cycle_length

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

#segment is imputed by AR or not
def imputed_status_from_aln(aln_vect,cycle_len):
    imputed_status = [int(aln_vect[0]["imputed"])]
    curr_seg_aln_number = 0
    for a_d in aln_vect:
        if a_d["seg_aln_number"] != curr_seg_aln_number:
            for i in range(curr_seg_aln_number+1,a_d["seg_aln_number"]):
                imputed_status.append(1)

            imputed_status.append(int(a_d["imputed"]))
            curr_seg_aln_number=a_d["seg_aln_number"]


    for i in range(curr_seg_aln_number+1,cycle_len):
        imputed_status.append(1)

    return imputed_status

#TODO: Implement contig end trimming
def trim_contigs(contig_cmap_vects,contig_placements,total_length):
    for cObj in contig_placements:
        cmap_vect = contig_cmap[vects[cObj.id]]
        first_lab,last_lab = cObj.aln_lab_ends
        if cmap_vect[first_lab-1] - cmap_vect[0] > unaligned_cutoff_frac*total_length:
            cObj.start_trim = True

        if cmap_vect[-1] - cmap_vect[last_lab-1] > unaligned_cutoff_frac*total_length:
            cObj.end_trim = True

#SET COLORS
def get_chr_colors():
    to_add = plt.cm.get_cmap(None, 4).colors[1:]
    color_vect = ["#f2bfca","indianred","salmon","burlywood",'#d5b60a',"xkcd:algae",to_add[0],"darkslateblue",
                 to_add[2],"#017374","#734a65","#bffe28","xkcd:darkgreen","#910951","xkcd:stone",
                 "xkcd:purpley","xkcd:topaz","lavender","darkseagreen","powderblue","#ff073a",to_add[1],"magenta"]

    chromosome_colors = dict(zip(["chr" + str(i) for i in range(1,24)],color_vect))
    return chromosome_colors

#TEMP SOLUTION (will break if too many consecutive overlaps)
def set_contig_height_shifts(contig_placements,contig_list):
    print "SETTING HEIGHTS"
    for ind,i in enumerate(contig_list[1:]):
        prevObj = contig_placements[contig_list[ind]]
        currObj = contig_placements[i]

        print "s",prevObj.abs_start_pos,prevObj.abs_end_pos
        print "t",currObj.abs_start_pos,currObj.abs_end_pos

        if currObj.abs_start_pos < prevObj.abs_end_pos:
            print "hit"
            #TODO: ROBUST FIX
            currObj.track_height_shift = prevObj.track_height_shift - 1.6

def construct_cycle_ref_placements(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj,isCircular):
    spacing_bp = seg_spacing*raw_cycle_length
    cycle_ref_placements = {}
    curr_start = 0.0 if isCircular else spacing_bp
    for ind,i in enumerate(cycle):
        seg_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        seg_end = curr_start+seg_len
        curr_obj = CycleVizElemObj(i[0],i[1],curr_start,seg_end)
        cycle_ref_placements[ind] = curr_obj
        next_start = seg_end
        mod_ind = ind % (len(prev_seg_index_is_adj)-1)
        if not prev_seg_index_is_adj[mod_ind]:
            next_start+=spacing_bp

        curr_start = next_start

    total_length = next_start
    return cycle_ref_placements,total_length

def place_cycle_segs_and_labels(cycle,ref_placements,seg_cmap_vects):
    cycle_seg_placements = {}
    for ind,i in enumerate(cycle):
        refObj = ref_placements[ind]
        segObj = copy.deepcopy(refObj)
        segObj.cmap_vect = seg_cmap_vects[i[0]]
        segObj.compute_label_posns()
        cycle_seg_placements[ind] = segObj

    return cycle_seg_placements

#create an object for each contig encoding variables such as position of start and end of contig (absolute ends)
#and positioning of contig labels
def place_contigs_and_labels(cycle_seg_placements,aln_vect,total_length,contig_cmap_vects):
    used_contigs = set()
    contig_aln_dict = defaultdict(list)
    wraparound = []
    contig_list = []
    for i in aln_vect:
        c_id = i["contig_id"]
        contig_aln_dict[c_id].append(i)
        if c_id not in contig_list: contig_list.append(c_id)

    contig_span_dict = {}
    for c_id,i_list in contig_aln_dict.iteritems():
        print "placing contigs computation step"
        print c_id
        cc_vect = contig_cmap_vects[c_id]
        san_f = i_list[0]["seg_aln_number"]
        sal_f = i_list[0]["seg_label"]
        cal_f = i_list[0]["contig_label"]
        san_l = i_list[-1]["seg_aln_number"]
        sal_l = i_list[-1]["seg_label"]
        cal_l = i_list[-1]["contig_label"]
        contig_dir = i_list[0]["contig_dir"]
        print san_f,sal_f,cal_f
        print san_l,sal_l,cal_l
        curr_contig_struct = CycleVizElemObj(c_id,contig_dir,None,None,cc_vect)

        #look up aln posns from cycle_seg_placements
        #look up position of first one
        segObj_start = cycle_seg_placements[san_f]
        seg_start_l_pos = segObj_start.label_posns[sal_f-1]

        #look up position of last one
        segObj_end = cycle_seg_placements[san_l]
        seg_end_l_pos = segObj_end.label_posns[sal_l-1]

        if seg_end_l_pos < seg_start_l_pos:
            print "hit add len"
            seg_end_l_pos = total_length + seg_end_l_pos

        #compute scaling
        scaled_seg_dist  = abs(seg_end_l_pos - seg_start_l_pos)*(1-contig_spacing/2)
        scaling_factor = scaled_seg_dist/(abs(cc_vect[cal_f-1] - cc_vect[cal_l-1]))
        #SET CONTIG SCALING FACTOR
        curr_contig_struct.scaling_factor = scaling_factor
        print scaling_factor,c_id

        abs_start_pos = seg_start_l_pos - (cc_vect[cal_f-1])*scaling_factor
        abs_end_pos = abs_start_pos + (cc_vect[-1])*scaling_factor
        
        # if contig_dir == "+":
        #     abs_start_pos = seg_start_l_pos - (cc_vect[cal_f-1])*scaling_factor
        #     # abs_end_pos = abs_start_pos + (cc_vect[-2] + (cc_vect[-1] - cc_vect[cal_l-1]))*scaling_factor
        #     abs_end_pos = abs_start_pos + (cc_vect[-1])*scaling_factor

        # else:
        #     # abs_start_pos = seg_end_l_pos - (cc_vect[-1] - cc_vect[cal_l-1])*scaling_factor
        #     abs_start_pos = seg_end_l_pos - (cc_vect[cal_f-1])*scaling_factor
        #     # abs_start_pos = seg_end_l_pos - (cc_vect[-1])*scaling_factor
        #     # abs_end_pos = seg_start_l_pos + (cc_vect[cal_f-1])*scaling_factor
        #     abs_end_pos = abs_start_pos + (cc_vect[-1])*scaling_factor
        #     # abs_end_pos = seg_start_l_pos + cc_vect[cal_f-1]*scaling_factor
        #     # abs_start_pos = (abs_end_pos - cc_vect[-2])*scaling_factor

        #s,t unset until label posns computed
        print "dists"
        print cc_vect[cal_f-1] - cc_vect[cal_l-1]
        print seg_end_l_pos - seg_start_l_pos
        print abs_end_pos - abs_start_pos

        curr_contig_struct.abs_start_pos = abs_start_pos
        curr_contig_struct.abs_end_pos = abs_end_pos

        #SET BOUNDARY ALN POSITIONS FROM TRACK
        curr_contig_struct.aln_bound_posns = (seg_start_l_pos,seg_end_l_pos)

        csl = min(i_list[-1]["contig_label"],i_list[0]["contig_label"])
        cel = max(i_list[-1]["contig_label"],i_list[0]["contig_label"])
        #SET FIRST AND LAST LABEL ALIGNED IN THE CONTIG
        curr_contig_struct.aln_lab_ends = (csl,cel)
        curr_contig_struct.compute_label_posns()
        contig_span_dict[c_id] = curr_contig_struct

    return contig_span_dict,contig_list


parser = argparse.ArgumentParser(description="Corrects and extends alignment paths to produce BioNano contig/AA segment scaffolds")
parser.add_argument("--om_alignments",help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
    action='store_true')
parser.add_argument("--cycles_file",help="cycles file",required=True)
parser.add_argument("--cycle",help="cycle number to visualize",required=True)
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("-g", "--graph", help="breakpoint graph file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--rot", help="number of segments to rotate counterclockwise",type=int,default=0)
parser.add_argument("--label_segs",help="label segs with graph IDs",action='store_true')
parser.add_argument("--feature_files",help="bed file list",nargs="+")
parser.add_argument("--feature_labels",help="bed feature names",nargs="+",default=[])
parser.add_argument("--log10_features",help="features which have been log10 scaled",nargs="+")
parser.add_argument("--transcript_features",help="features which are transcript info",nargs="+")
parser.add_argument("--oncogenes_only",help="only show oncogenes in the annotation track",action='store_true')

args = parser.parse_args()

if not args.sname:
    args.sname = os.path.split(args.cycles_file)[1].split(".")[0] + "_"

outdir = os.path.dirname(args.sname)
if outdir and not os.path.exists(outdir):
    os.makedirs(outdir)

fname = args.sname + "_cycle_" + args.cycle

chromosome_colors = get_chr_colors()
plt.clf()
fig, ax = plt.subplots()
patches = []
f_color_v = []
e_color_v = []
lw_v = []

cycles,segSeqD,circular_D = parse_cycles_file(args.cycles_file)
cycle_num = args.cycle
isCircular = circular_D[cycle_num]
cycle = cycles[cycle_num]
prev_seg_index_is_adj = adjacent_segs(cycle,segSeqD,isCircular)
raw_cycle_length = get_raw_cycle_length(cycle,segSeqD,isCircular,prev_seg_index_is_adj)
# start_points,total_length = get_seg_locs_from_cycle(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj)
ref_placements,total_length = construct_cycle_ref_placements(cycle,segSeqD,raw_cycle_length,prev_seg_index_is_adj,isCircular)

if args.om_alignments:
    seg_cmaps = parse_cmap(args.segs,True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    seg_cmap_lens = get_cmap_lens(args.segs)
    cycle_seg_placements = place_cycle_segs_and_labels(cycle,ref_placements,seg_cmap_vects)

    contig_cmaps = parse_cmap(args.contigs,True)
    contig_cmap_vects = vectorize_cmaps(contig_cmaps)

    aln_vect,meta_dict = parse_alnfile(args.path_alignment)
    contig_cmap_lens = get_cmap_lens(args.contigs)
    #cycle_seg_placements,aln_vect,total_length,contig_cmap_vects
    contig_placements,contig_list = place_contigs_and_labels(cycle_seg_placements,aln_vect,total_length,contig_cmap_vects)

    if not args.feature_labels:
        outside = False

    #plot segs

    #seg_height_shifts = {x:0 for x in seg_placements}
    #seg_placements,total_length,bar_height,color,seg_text = None
    plot_cmap_track(cycle_seg_placements,total_length,outer_bar+segment_bar_height,"darkorange")


    #plot contigs

    # check overlaps of contigs and adjust heights accordingly
    contig_height_shifts = set_contig_height_shifts(contig_placements,contig_list)
    plot_cmap_track(contig_placements,total_length,outer_bar+contig_bar_height,"cornflowerblue",seg_id_labels=True)


    #plot alignments
    plot_alignment(contig_placements,cycle_seg_placements,total_length)

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

imputed_status = imputed_status_from_aln(aln_vect,len(cycle))
#plot_ref_genome(ref_placements,cycle,total_length,segSeqD,imputed_status,label_segs=True)
plot_ref_genome(ref_placements,cycle,total_length,segSeqD,imputed_status,args.label_segs)

ax.set_xlim(-(outer_bar+1.25), (outer_bar+1.25))
ax.set_ylim(-(outer_bar+1.25), (outer_bar+1.25))
chrom_set = set()
for i in cycle:
    chrom_set.add(segSeqD[i[0]][0])

sorted_chrom = sorted(chrom_set,key=lambda x:int(x.rsplit("chr")[-1]))
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
print "done plotting bionano"
