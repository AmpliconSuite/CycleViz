#!/usr/bin/env python

import sys
sys.path.insert(0, "/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/BioNanoCycleBuilder")
import os
import bisect
import argparse
import numpy as np
from collections import defaultdict
# import matplotlib.path as mpath
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
import matplotlib._color_data as mcd
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
from intervaltree import Interval, IntervalTree
from bionanoUtil import *

contig_spacing = 0.01
seg_spacing = 0.0075
bar_width = 2.0/3
global_rot = 90.0
center_hole = 0.75
outer_bar = 10 
bed_spacing = .5

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
    text_angle = (start_angle + end_angle)/2.0
    if end_angle < 0 and start_angle > 0:
        end_angle+=360

    return start_angle,end_angle,text_angle

def parse_genes(chrom):
    print("Building interval tree for chrom " + chrom)
    t = IntervalTree()
    with open("refGene.txt") as infile:
        for line in infile:
            fields = line.rsplit("\t")
            #gene = fields[-2]
            currChrom = fields[2]
            tstart = int(fields[4])
            tend = int(fields[5])
            if chrom == currChrom:
                t[tstart:tend] = fields

    return t

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
            normStart = currStart - max(0,tstart-pTup[1])
            normEnd = currStart - min(pTup[2]-pTup[1],tend-pTup[1])

        else:
            normEnd = currStart - min(pTup[2]-pTup[1],pTup[2]-tstart)
            normStart = currStart - max(0,pTup[2] - tend)

        start_angle = normStart/total_length_with_spacing*360
        end_angle = normEnd/total_length_with_spacing*360
        text_angle = (start_angle + end_angle)/2.0
        gene_to_locations[i].append((start_angle/360.,end_angle/360.))
        if end_angle < 0 and start_angle > 0:
            end_angle+=360
        
    
        patches.append(mpatches.Wedge((0,0), outer_bar, end_angle, start_angle, width=bar_width/2.0))
        f_color_v.append('k')
        e_color_v.append('k')
        lw_v.append(0)

        # x,y = pol2cart(outer_bar+(bar_width/2.0),(text_angle/360*2*np.pi))
        x_t,y_t = pol2cart(outer_bar + bar_width + 0.75,(text_angle/360*2*np.pi))
        #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.4)
        
        if (abs(text_angle) > 90 and abs(text_angle) < 360):
            text_angle+=180 
            ax.text(x_t,y_t,i,color='grey',rotation=text_angle,
                ha='right',fontsize=6.5,rotation_mode='anchor')
        
        else:
            ax.text(x_t,y_t,i,color='grey',rotation=text_angle,
                ha='left',fontsize=6.5,rotation_mode='anchor')

        for exon in e_posns:
            if exon[1] > pTup[1] and exon[0] < pTup[2]:
                if strand == "+":
                    normStart = currStart - max(1,exon[0]-pTup[1])
                    normEnd = currStart - min(pTup[2]-pTup[1],exon[1]-pTup[1])

                else:
                    normEnd = currStart - min(pTup[2]-pTup[1],pTup[2]-exon[0])
                    normStart = currStart - max(1,pTup[2] - exon[1])

                start_angle, end_angle, _ = start_end_angle(normStart,normEnd,total_length_with_spacing)


                patches.append(mpatches.Wedge((0,0), outer_bar-bar_width/2.0, end_angle, start_angle, width=bar_width/2.0))
                f_color_v.append('k')
                e_color_v.append('k')
                lw_v.append(0)
        
#plot the reference genome
def plot_ref_genome(start_points,lens,cycle,total_length_with_spacing,segSeqD,imputed_status,label_segs=True):
    font0 = FontProperties()
    seg_posns = []
    global_start = total_length_with_spacing*(global_rot/360.0)
    for ind,sp in enumerate(start_points):
        start_point = int(global_start - sp)
        seg_coord_tup = segSeqD[cycle[ind][0]]
        start_angle, end_angle, _ = start_end_angle(start_point,start_point - lens[ind],total_length_with_spacing)
        
        #makes the reference genome wedges    
        patches.append(mpatches.Wedge((0,0), outer_bar, end_angle, start_angle, width=bar_width))
        chrom = segSeqD[cycle[ind][0]][0]
        f_color_v.append(chromosome_colors[chrom])
        e_color_v.append('grey')
        lw_v.append(0.2)
        
        #makes the ticks on the reference genome wedges
        if cycle[ind][1] == "+":
            posns = zip(range(seg_coord_tup[1],seg_coord_tup[2]+1),np.arange(start_point,start_point-lens[ind]-1,-1))
        else:
            posns = zip(np.arange(seg_coord_tup[2],seg_coord_tup[1]-1,-1),np.arange(start_point,start_point-lens[ind]-1,-1))

        tick_freq = max(10000,10000*int(np.floor(total_length_with_spacing/1000000)))
        for j in posns:
            if j[0] % tick_freq == 0:
                text_angle = j[1]/total_length_with_spacing*360
                x,y = pol2cart(outer_bar,(text_angle/360*2*np.pi))
                x_t,y_t = pol2cart(outer_bar + 0.2,(text_angle/360*2*np.pi))
                ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.2)
                
                if (abs(text_angle) > 90 and abs(text_angle) < 360):
                    text_angle-=180
                    ha = "right"
                    txt = str(int(round((j[0])/10000))) + " "

                else:
                    ha = "left"
                    txt = " " + str(int(round((j[0])/10000)))

                ax.text(x_t,y_t,txt,color='grey',rotation=text_angle,
                ha=ha,fontsize=4,rotation_mode='anchor')
    
        gene_tree = parse_genes(seg_coord_tup[0])
        relGenes = rel_genes(gene_tree,seg_coord_tup)
        #plot the gene track
        plot_gene_track(start_point,relGenes,seg_coord_tup,total_length_with_spacing,cycle[ind][1])

        #label the segments by number in cycle
        mid_sp = start_point - lens[ind]/2
        text_angle = mid_sp/total_length*360.
        x,y = pol2cart(outer_bar-0.9,(text_angle/360.*2.*np.pi))
        font = font0.copy()
        if imputed_status[ind]:
            font.set_style('italic')
            font.set_weight('bold')

        if (abs(text_angle) > 90 and abs(text_angle) < 360):
            text_angle-=180
            ha = "left"

        else:
            ha = "right"

        if label_segs:
            ax.text(x,y,cycle[ind][0]+cycle[ind][1],color='grey',rotation=text_angle,
                ha=ha,fontsize=5,fontproperties=font,rotation_mode='anchor')

#plot gene transcript info
#transcipt count info is log-2
def plot_gene_transcript(r,feat_min,feat_scaling,total_length_with_spacing):
    for gene,pos_list in gene_to_locations.iteritems():
        try:
            value = np.log2(transcript_dict[gene])
            for pos_tup in pos_list:
                curr_rad = r+(value-feat_min)/feat_scaling
                # x_s,y_s = pol2cart(curr_rad,pos_tup[0]*2*np.pi)
                # x_e,y_e = pol2cart(curr_rad,pos_tup[1]*2*np.pi)
                lstart,lend = tuple(sorted(pos_tup))
                line_points = np.arange(lstart*2*np.pi,lend*2*np.pi,0.001)
                x_vals,y_vals = polar_series_to_cartesians(line_points,curr_rad)
                ax.plot(x_vals,y_vals,color='k',linewidth=1)

        except KeyError:
            pass

#plot genomic features
def bed_feat_bounds(bed_feat_list,n_guides):
    #get min and max vals for each feature
    #set bar heights
    feat_to_min_max = {x:(0,0) for x in bed_feat_list}
    feat_to_bar_vals = {x:[] for x in bed_feat_list}
    for bar_level,curr_feat in enumerate(bed_feat_list):
        feat_min,feat_max = float('inf'),-float('inf')
        if curr_feat in transcript_set:
            for i in gene_to_locations:
                try:
                    val = np.log2(transcript_dict[i])
                    if val > feat_max:
                        feat_max = val
                    elif val < feat_min:
                        feat_min = val

                except KeyError:
                    print "couldn't find gene " + i
                    continue

        else:
            for chrom in bed_feat_dict[curr_feat]:
                for x in bed_feat_dict[curr_feat][chrom]:
                    if x[2] > feat_max:
                        feat_max = x[2]
                    elif x[2] < feat_min:
                        feat_min = x[2]

        #guide bars
        diff = (feat_max - feat_min)/n_guides
        top = int(round(feat_max))
        bottom = int(round(feat_min))
        if bottom < feat_min: feat_min = bottom
        if feat_max < top: feat_max = top
        for i in range(bottom,top+1):
            feat_to_bar_vals[curr_feat].append(i)

        feat_to_min_max[curr_feat] = (feat_min,feat_max)

    return feat_to_min_max,feat_to_bar_vals

#plot bed files of genomic features
def plot_bed_features(start_points,lens,cycle,total_length_with_spacing,segSeqD,bed_feat_list):
    n_guides = 5
    n_feats = float(len(bed_feat_list))
    interior_space = outer_bar - (bar_width+0.1) - center_hole
    space_per_bed = interior_space/(n_feats + bed_spacing*(n_feats-1))
    seg_posns = []
    global_start = total_length_with_spacing*(global_rot/360.0)
    top_bar_base = outer_bar - (bar_width+0.1)
    #get min and max vals for each feature and set bar heights
    feat_to_min_max,feat_to_bar_vals = bed_feat_bounds(bed_feat_list,n_guides)

    #plot transcribed stuff
    for i in transcript_set:
        feat_min,feat_max = feat_to_min_max[i]
        for ind,j in enumerate(bed_feat_list):
            if j == i: bar_level = ind

        r = top_bar_base - (bar_level+1)*(space_per_bed + bed_spacing)
        feat_scaling = float(feat_max-feat_min)/space_per_bed
        plot_gene_transcript(r,feat_min,feat_scaling,total_length_with_spacing)


    for ind,sp in enumerate(start_points):
        start_point = int(global_start - sp)
        seg_coord_tup = segSeqD[cycle[ind][0]]
        start_angle, end_angle, text_angle = start_end_angle(start_point,start_point - lens[ind],total_length_with_spacing)

        for bar_level,curr_feat in enumerate(bed_feat_list):
            feat_min,feat_max = feat_to_min_max[curr_feat]
            feat_scaling = float(feat_max-feat_min)/space_per_bed
            # r = curr_bar_base - diff - bed_spacing
            r = top_bar_base - (bar_level+1)*(space_per_bed + bed_spacing)

            if curr_feat not in transcript_set:
                full_feat_list = bed_feat_dict[curr_feat][seg_coord_tup[0]]
                s_posns = [x[0] for x in full_feat_list]
                left_start = bisect.bisect_left(s_posns,seg_coord_tup[1])
                right_end = bisect.bisect_right(s_posns,seg_coord_tup[2])
                #left overhang
                curr_feat_list = full_feat_list[left_start:right_end]
                if curr_feat_list[0][0] < seg_coord_tup[1]:
                    curr_feat_list[0][0] = seg_coord_tup[1]

                curr_feat_list[-1][1] = seg_coord_tup[2]

                x_values = []
                y_values = []
                colors = []
                #the value for a region in the bed is its midpoint
                for x_ind,x in enumerate(curr_feat_list):
                    lstart = start_point - (x[0] - seg_coord_tup[1])
                    lend = start_point - (x[1] - seg_coord_tup[1])
                    curr_rad = r+(x[2]-feat_min)/feat_scaling
                    x_s,y_s = pol2cart(curr_rad,lstart/total_length_with_spacing*2*np.pi)
                    x_e,y_e = pol2cart(curr_rad,lend/total_length_with_spacing*2*np.pi)
                    x_values.extend([x_s,x_e])
                    y_values.extend([y_s,y_e])

                ax.plot(x_values,y_values,color='k',linewidth=0.3)
            
            #plot guide bars
            for i in feat_to_bar_vals[curr_feat]:
                curr_rad = r+(i-feat_min)/feat_scaling
                patches.append(mpatches.Wedge((0,0), curr_rad, end_angle, start_angle, width=0.05))
                f_color_v.append('lightgrey')
                e_color_v.append('grey')
                lw_v.append(0)

        # break

    #write tick mark indicators
    for bar_level,curr_feat in enumerate(bed_feat_list):
        for i in feat_to_bar_vals[curr_feat]:
            feat_min,feat_max = feat_to_min_max[curr_feat]
            feat_scaling = float(feat_max-feat_min)/space_per_bed
            scaled_i = (i-feat_min)/feat_scaling
            r = top_bar_base - (bar_level+1)*(space_per_bed + bed_spacing)
            curr_rad = r+scaled_i
            x,y = pol2cart(curr_rad,global_rot/360*2*np.pi)
            bar_label = round_to_1_sig(10**i) if curr_feat in log10_set else 2**i
            if bar_label >= 1:
                bar_label = int(bar_label)

            ax.text(x,y,str(bar_label),color='grey',
                    ha="left",fontsize=4,rotation_mode='anchor')


#plot cmap track
def plot_cmap_track(start_points,aln_nums,total_length,cycle,bar_height,color,cmap_vects,seg_text = None):
    cycle_label_locs = defaultdict(list)
    global_start = total_length*(global_rot/360.0)
    for start_point,seg_tup,a_num in zip(start_points,cycle,aln_nums):
        #draw the tracks
        start_point = (global_start - start_point)
        curr_seg_cmap_vect = cmap_vects[seg_tup[0]]
        seg_cmap_lab_order = range(1,len(curr_seg_cmap_vect))
        start_angle, end_angle, _ = start_end_angle(start_point,start_point - curr_seg_cmap_vect[-1],total_length)

        patches.append(mpatches.Wedge((0,0), bar_height + bar_width, end_angle, start_angle, width=bar_width))
        f_color_v.append(color)
        e_color_v.append('k')
        lw_v.append(0)

        # draw the labels
        if seg_tup[1] == '-':
            seg_len = cmap_vects[seg_tup[0]][-1]
            curr_seg_cmap_vect = [seg_len - x for x in curr_seg_cmap_vect]

        lab_locs = {}
        for lab_id,i in zip(seg_cmap_lab_order,curr_seg_cmap_vect[:-1]):
            label_rads = (start_point - i)/total_length*2*np.pi
            x,y = pol2cart(bar_height,label_rads)
            x_t,y_t = pol2cart(bar_height+bar_width,label_rads)
            ax.plot([x,x_t],[y,y_t],color='k',alpha=0.9,linewidth=0.2)
            lab_locs[lab_id] = label_rads

        cycle_label_locs[a_num] = lab_locs

        if seg_text:
            mid_sp = start_point - curr_seg_cmap_vect[0]/2 if seg_tup[1] == "-" else start_point - curr_seg_cmap_vect[-1]/2
            text_angle = mid_sp/total_length*360.
            x,y = pol2cart(bar_height-0.7,(text_angle/360.*2.*np.pi))
            
            if (abs(text_angle) > 90 and abs(text_angle) < 360):
                text_angle-=180
                ha = "left"

            else:
                ha = "right"

            ax.text(x,y,seg_text,color='grey',rotation=text_angle,
                    ha=ha,fontsize=5,rotation_mode='anchor')


    return cycle_label_locs

#plot the connecting lines for the bionano track
def plot_alignment(segment_lab_locs,contig_lab_loc_dict,aln_vect,segs_base,contigs_top,total_length):
    for a_d in aln_vect:
        contig_label_dict = contig_lab_loc_dict[a_d["contig_id"]]
        seg_label_dict = segment_lab_locs[a_d["seg_aln_number"]]
        c_l_loc = contig_label_dict[a_d["contig_label"]]
        s_l_loc = seg_label_dict[a_d["seg_label"]]
        x_c,y_c = pol2cart(contigs_top,c_l_loc)
        x_s,y_s = pol2cart(segs_base,s_l_loc)
        ax.plot([x_c, x_s], [y_c, y_s], color="grey",linewidth=0.2)

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

#parse the alignment file
#refactor to util
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

    # while aln_vect[0]["seg_aln_number"] != 0:
    #     aln_vect = aln_vect[1:] + [aln_vect[0]]

    return aln_vect,meta_dict,

#parse cycles file
#refactor to util?
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

#format refGeneID ENSG value
def parse_transcript_data(transcript_file):
    transcript_data = {}
    with open(transcript_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            transcript_data[fields[0]] = float(fields[2])

    return transcript_data

#parse bed file
def parse_bed_file(bed_file):
    bed_list = []
    with open(bed_file) as infile:
        for line in infile:
            fields = line.rstrip().split()
            fields[1] = int(fields[1])
            fields[2] = int(fields[2])
            fields[3] = float(fields[3])
            bed_list.append(fields)

    return bed_list

#locations of each bionano contig on the figure
def get_contig_locs(aln_vect,start_points,cycle,total_length):
    # #go through the seg_seq and calculate total lengths of alignments
    contig_list = []
    first_aln = []
    prev_c_id = None
    for a_d in aln_vect:
        if a_d["contig_id"] != prev_c_id:
            prev_c_id = a_d["contig_id"]
            contig_list.append((prev_c_id,a_d["contig_dir"]))
            first_aln.append(a_d)

        prev_a_d = a_d

    total_contig_length = sum([contig_cmap_lens[i[0]] for i in contig_list])
    spacing = contig_spacing*total_contig_length*len(contig_list)
    total_contig_length+=spacing

    no_rot = False
    if total_contig_length < 0.9*total_length:
        print "no rotation"
        no_rot = True

    half_rot = seg_spacing*total_length*len(cycle)/2.0 if not no_rot else 0

    contig_start_points = []
    last_end = 0
    for ind,i in enumerate(contig_list):
        print i
        a_d = first_aln[ind]
        zero_shifted_start = start_points[a_d["seg_aln_number"]]
        if a_d["seg_dir"] == "+":
            zero_shifted_start+=seg_cmaps[a_d["seg_id"]][a_d["seg_label"]]
        else:
            zero_shifted_start+=(seg_cmap_vects[a_d["seg_id"]][-1] - seg_cmaps[a_d["seg_id"]][a_d["seg_label"]])

        if a_d["contig_dir"] == "+":
            zero_shifted_start-=contig_cmaps[i[0]][a_d["contig_label"]]
        else:
            zero_shifted_start+=(contig_cmap_vect[i[0]][-1]-contig_cmaps[i[0]][a_d["contig_label"]])

        curr_sp = zero_shifted_start
        print curr_sp,last_end
        if last_end:
            # last_end = contig_start_points[-1]+contig_cmap_vect[contig_list[ind-1][0]][-1]
            if curr_sp < last_end + total_length*contig_spacing:
                curr_sp = last_end + total_length*contig_spacing

        elif len(contig_list)==1:
            curr_sp+=half_rot

        elif ind == 0:
            if curr_sp > 0:
                curr_sp-=half_rot
            else:
                curr_sp+=half_rot

        contig_start_points.append(curr_sp)
        last_end = contig_start_points[-1] + contig_cmap_vect[i[0]][-1]

    print "circle " + str(isCircular)
    #don't want contig hanging over at the end and overlapping with first contig
    if (curr_sp + contig_cmap_vect[i[0]][-1] > total_contig_length) and not isCircular:
        total_contig_length = curr_sp + contig_cmap_vect[i[0]][-1]


    elif (curr_sp + contig_cmap_vect[i[0]][-1] > total_contig_length + contig_start_points[0]):
        total_contig_length = curr_sp + contig_cmap_vect[i[0]][-1]

    elif isCircular and (last_end < total_contig_length + contig_start_points[0]):
        total_contig_length = last_end - contig_start_points[0] - half_rot


    return dict(zip(contig_list,contig_start_points)),total_contig_length

#get start locations for a cycle
def get_seg_locs_from_cycle(cycle,segSeqD):
    lens = []
    consecutives = set()
    for i in cycle:
        curr_len = segSeqD[i[0]][2] - segSeqD[i[0]][1]
        lens.append(curr_len)

    for i in range(len(cycle)-1):
        k1 = cycle[i]
        k2 = cycle[i+1]
        if abs(segSeqD[k1[0]][2] - segSeqD[k2[0]][1]) == 1 and segSeqD[k1[0]][0] == segSeqD[k2[0]][0]:
            consecutives.add(i)

    if abs(segSeqD[k2[0]][2] - segSeqD[cycle[0][0]][1]) == 1 and segSeqD[k2[0]][0] == segSeqD[cycle[0][0]][0]:
            consecutives.add(i+1)


    total_seg_len = sum(lens)
    seg_padding = seg_spacing*(total_seg_len-len(consecutives))

    start_points = [0]
    total_broken = 0
    for ind,i in enumerate(lens[:-1]):
        #check if i and i+1 adjacent
        curr_padding = seg_padding if ind not in consecutives else 0
        start_points.append(start_points[-1] + i + curr_padding)

    curr_padding = seg_padding if ind+1 not in consecutives else 0 
    total_length = start_points[-1] + lens[-1] + curr_padding

    rot_len = start_points[args.rot]
    start_points = [x - rot_len for x in start_points]

    return start_points,lens,float(total_length)

def feat_bed_to_lookup(bed_list):
    bed_dict = defaultdict(list)
    for i in bed_list:
        if i[0].startswith("hs"):
            i[0] = i[0].replace("hs","chr")

        bed_dict[i[0]].append(i[1:])

    return bed_dict

# if __name__ == '__main__':
# Parses the command line arguments
parser = argparse.ArgumentParser(description="Corrects and extends alignment paths to produce BioNano contig/AA segment scaffolds")
parser.add_argument("--bionano_alignments",help="Enable Bionano visualizations (requires contigs,segs,key,path_alignment args)",
    action='store_true')
parser.add_argument("--cycles_file",help="cycles file",required=True)
parser.add_argument("--cycle",help="cycle number to visualize",required=True)
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--rot", help="number of segments to rotate counterclockwise",type=int,default=0)
parser.add_argument("--label_segs",help="label segs with graph IDs",action='store_true')
parser.add_argument("--feature_files",help="bed file list",nargs="+")
parser.add_argument("--feature_labels",help="bed feature names",nargs="+",default=[])
parser.add_argument("--log10_features",help="features which have been log10 scaled",nargs="+")
parser.add_argument("--transcript_features",help="features which are transcript info",nargs="+")

args = parser.parse_args()

if not args.sname:
    samp_name = os.path.split(args.cycles_file)[1].split(".")[0] + "_"
else:
    samp_name = args.sname.rsplit("/")[-1]

fname = samp_name + "_cycle_" + args.cycle

log10_set = set(args.log10_features) if args.log10_features else set()
transcript_set = set(args.transcript_features) if args.transcript_features else set()
transcript_dict = {}

bed_feat_dict = {}
transcript_data = {}
if args.feature_files:
    for i,j in zip(args.feature_files,args.feature_labels):
        #feature name -> chromosome -> ordered list of positions
        if i.endswith(".bed"):
            bed_list = parse_bed_file(i)
            bed_feat_dict[j] = feat_bed_to_lookup(bed_list)
        elif j in transcript_set:
            transcript_dict = parse_transcript_data(i)
            bed_feat_dict[j] = []

#SET COLORS
to_add = plt.cm.get_cmap(None, 4).colors[1:]
color_vect = ["#f2bfca","indianred","salmon","burlywood",'#d5b60a',"xkcd:algae",to_add[0],"darkslateblue",
             to_add[2],"#017374","#734a65","#bffe28","xkcd:darkgreen","#910951","xkcd:stone",
             "xkcd:purpley","xkcd:topaz","lavender","darkseagreen","powderblue","#ff073a",to_add[1],"magenta"]

chromosome_colors = dict(zip(["chr" + str(i) for i in range(1,24)],color_vect))

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

start_points,lens,total_length = get_seg_locs_from_cycle(cycle,segSeqD)
if args.bionano_alignments:
    #unfinished utility for bionano
    #parse cmaps
    # ref_cmaps = parse_cmap("/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/BioNanoCycleBuilder/hg19_BspQI.cmap")
    seg_cmaps = parse_cmap(args.segs,True)
    seg_cmap_vects = vectorize_cmaps(seg_cmaps)
    contig_cmaps = parse_cmap(args.contigs,True)

    contig_cmap_vect = vectorize_cmaps(contig_cmaps)
    #get cmap lens
    aln_vect,meta_dict = parse_alnfile(args.path_alignment)
    seg_cmap_lens = get_cmap_lens(args.segs)
    contig_cmap_lens = get_cmap_lens(args.contigs)
    contig_locs,contig_total_length = get_contig_locs(aln_vect,start_points,cycle,total_length)
    if contig_total_length > total_length:
        print "setting total length to contig total length"
        total_length = contig_total_length

    if not args.feature_labels:
        outside = False
    #plot contigs
    contig_bar_height = -13/3
    segment_bar_height = -8.0/3

    #outside
    #contig_bar_height = 4
    #segment_bar_height = 8.0/3

    contig_lab_dict = {}
    for c_tup,start_posn in contig_locs.iteritems():
        contig_cycle = [c_tup]
        contig_cmap = contig_cmaps[c_tup[0]]
        contig_vect = {c_tup[0]:[contig_cmap[i] for i in sorted(contig_cmap.keys())]}
        seg_text = c_tup[0]
        contig_label_posns = plot_cmap_track([start_posn],[0],total_length,contig_cycle,outer_bar+contig_bar_height,
            "cornflowerblue",contig_vect,seg_text)
        contig_lab_dict[c_tup[0]] = contig_label_posns[0]

    #plot segs
    aln_nums = range(0,len(cycle))
    if len(cycle) > 1 and aln_vect[-1]["seg_aln_number"] == "0":
        aln_nums.append(0)

    segment_label_list = plot_cmap_track(start_points,aln_nums,total_length,cycle,outer_bar+segment_bar_height,
        "darkorange",seg_cmap_vects)
    #plot alignment
    plot_alignment(segment_label_list,contig_lab_dict,aln_vect,outer_bar+segment_bar_height,
        outer_bar+contig_bar_height+bar_width,total_length)

    lens = [seg_cmap_vects[x[0]][-1] for x in cycle]


imputed_status = imputed_status_from_aln(aln_vect,len(cycle))

plot_ref_genome(start_points,lens,cycle,total_length,segSeqD,imputed_status,args.label_segs)

ax.set_xlim(-(outer_bar+1), (outer_bar+1))
ax.set_ylim(-(outer_bar+1), (outer_bar+1))
chrom_set = set()
for i in cycle:
    chrom_set.add(segSeqD[i[0]][0])

# for i in range(1,24):
#     chrom_set.add("chr" + str(i))

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

#handles basic (non-bionano case)
if args.feature_labels:
    #make a separate figure
    plt.clf()
    fig, ax = plt.subplots()
    patches = []
    f_color_v = []
    e_color_v = []
    lw_v = []
    ax.set_xlim(-(outer_bar+1), (outer_bar+1))
    ax.set_ylim(-(outer_bar+1), (outer_bar+1))

    plot_ref_genome(start_points,lens,cycle,total_length,segSeqD,imputed_status,args.label_segs)
    #plot the bed-files
    plot_bed_features(start_points,lens,cycle,total_length,segSeqD,args.feature_labels)


    plt.legend(handles=legend_patches,fontsize=8,loc=3,bbox_to_anchor=(-.3,.15))
    
    # plt.legend(handles=legend_patches,fontsize=5)
    p = PatchCollection(patches)
    p.set_facecolor(f_color_v)
    p.set_edgecolor(e_color_v)
    p.set_linewidth(lw_v)
    ax.add_collection(p)
    ax.set_aspect(1.0)
    plt.axis('off')
    plt.savefig(fname + '_genomic_features.png',dpi=600)
    plt.savefig(fname + '_genomic_features.pdf',format='pdf')
    plt.close()
    
