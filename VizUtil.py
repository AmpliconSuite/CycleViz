import sys
import os
import copy
import bisect
import numpy as np
import matplotlib
matplotlib.use('Agg')
from collections import defaultdict
from matplotlib import pyplot as plt
from intervaltree import Interval, IntervalTree

contig_spacing = 1./100
unaligned_cutoff_frac = 1./60

def round_to_1_sig(x):
    if x == 0:
        return 0.

    return round(x, -int(np.floor(np.log10(abs(x)))))

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
            # for i in self.cmap_vect[:-1]:
            #     self.label_posns.append(self.scaling_factor*i + self.abs_start_pos)

            self.label_posns = self.label_posns[::-1]
            # if self.id == "161":
            #     print self.label_posns[self.aln_lab_ends[0]:self.aln_lab_ends[-1] + 1]

    def to_string(self):
        return "{}{} | Start: {} | End: {} | scaling {}".format(self.id,self.direction,str(self.abs_start_pos),str(self.abs_end_pos),str(self.scaling_factor))

    def trim_obj_ends(self,total_length):
        #use the .end_trim,.start_trim to chop off.
        #overhang should go to 1/3 of the unaligned cutoff threshold

        if self.start_trim:
            p_abs_start = self.abs_start_pos
            print "DOING s TRIM on " + self.id
            self.abs_start_pos = self.aln_bound_posns[0] - unaligned_cutoff_frac*total_length/4.
            # if self.direction == "+":
            #     self.abs_start_pos = self.aln_bound_posns[0] - unaligned_cutoff_frac*total_length/4.

            # else:
            #     self.abs_start_pos = self.aln_bound_posns[0] - unaligned_cutoff_frac*total_length/4.

            if self.direction == "-":
                self.update_label_posns(p_abs_start - self.abs_start_pos)
            print "now",self.abs_start_pos



        if self.end_trim:
            p_abs_end = self.abs_end_pos
            print "DOING e TRIM on " + self.id
            print self.aln_bound_posns
            self.abs_end_pos = self.aln_bound_posns[-1] + unaligned_cutoff_frac*total_length/4.
            # if self.direction == "+":
            #     self.abs_end_pos = self.aln_bound_posns[-1] + unaligned_cutoff_frac*total_length/4.

            # else:
            #     self.abs_end_pos = self.aln_bound_posns[-1] + unaligned_cutoff_frac*total_length/4.
            if self.direction == "-":
                self.update_label_posns(p_abs_end - self.abs_end_pos)

            print "now",self.abs_end_pos

        
    #update label positions after trimming contigs
    def update_label_posns(self,s_diff):
        print "diff",s_diff
        for ind in range(len(self.label_posns)):
            self.label_posns[ind]-=s_diff

#SET COLORS
def get_chr_colors():
    to_add = plt.cm.get_cmap(None, 4).colors[1:]
    # color_vect = ["#ffe8ed","indianred","salmon","burlywood",'#d5b60a',"xkcd:algae",to_add[0],"darkslateblue",
    #              to_add[2],"#017374","#734a65","#bffe28","xkcd:darkgreen","#910951","xkcd:stone",
    #              "xkcd:purpley","xkcd:brown","lavender","darkseagreen","powderblue","#ff073a",to_add[1],
    #              "magenta","plum"]

    color_vect= ["aqua","gainsboro","salmon","bisque",'goldenrod',"xkcd:algae",to_add[0],"darkslateblue",
                 "yellow","sienna","purple","#bffe28","xkcd:darkgreen","#910951","xkcd:stone",
                 "xkcd:purpley","xkcd:brown","lavender","darkseagreen","powderblue","crimson",to_add[1],
                 "fuchsia","pink"]

    chrnames = [str(i) for i in (range(1,23) + ["X","Y"])]
    chromosome_colors = dict(zip(["chr" + i for i in chrnames],color_vect))
    for i in range(len(chrnames)):
        chromosome_colors[chrnames[i]] = color_vect[i]
    return chromosome_colors

#parse the breakpoint graph to indicate for two endpoints if there is an edge.
def parse_BPG(BPG_file):
    bidirectional_edge_dict = defaultdict(set)
    seg_end_pos_d = {}
    seqnum = 0
    with open(BPG_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if not fields:
                continue

            if fields[0] in ["concordant","discordant"]:
                e_rep = fields[1].rsplit("->")
                start = e_rep[0][:-1]
                end = e_rep[1][:-1]
                bidirectional_edge_dict[start].add(end)
                bidirectional_edge_dict[end].add(start)

            elif fields[0] == "sequence":
                seqnum+=1
                seg_end_pos_d[str(seqnum)] = (fields[1][:-1],fields[2][:-1])

    return bidirectional_edge_dict,seg_end_pos_d

#extract oncogenes from the bushman lab file, assumes refseq name in last column
def parse_gene_subset_file(gene_list_file):
    gene_set = set()
    with open(gene_list_file) as infile:
        for line in infile:
            fields = line.rstrip().split()
            if not fields:
                continue

            gene_set.add(fields[-1].strip("\""))

    return gene_set

def parse_genes(chrom,ref):
    # print("Building interval tree for chrom " + chrom)
    t = IntervalTree()
    # with open(os.environ['AR_SRC'] + "/hg19_refGene.txt") as infile:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    refGene_name = "refGene_" + ref + ".txt"
    with open(os.path.join(__location__, refGene_name)) as infile:
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

def rel_genes(chrIntTree,pTup, gene_set = set()):    
    relGenes = {}
    chrom = pTup[0]
    overlappingT = chrIntTree[pTup[1]:pTup[2]]
    gene_set_only = (len(gene_set) == 0)
    for i in overlappingT:
        gene = i.data[-4].rsplit("-")[0]
        tstart = int(i.data[4])
        tend = int(i.data[5])
        e_s_posns = [int(x) for x in i.data[9].rsplit(",") if x]
        e_e_posns = [int(x) for x in i.data[10].rsplit(",") if x] 
        is_other_feature = (gene.startswith("LOC") or gene.startswith("LINC") or gene.startswith("MIR"))
        if gene_set_only:
            gene_set.add(gene)

        if not is_other_feature and gene in gene_set:
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

def pair_is_edge(a_id,b_id,a_dir,b_dir,bpg_dict,seg_end_pos_d):
    rObj1_end = seg_end_pos_d[a_id][-1] if a_dir == "+" else seg_end_pos_d[a_id][0]
    rObj2_start = seg_end_pos_d[b_id][0] if b_dir == "+" else seg_end_pos_d[b_id][-1]
    return rObj1_end in bpg_dict[rObj2_start]

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
                isCycle = True
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
                        isCycle = False

                cycles[lineD["Cycle"]] = curr_cycle
                circular_D[lineD["Cycle"]] = isCycle

    return cycles,segSeqD,circular_D

def check_segdup(aln_vect,cycle,circular):
    print("Checking if segdup")
    #iterate over and delete the second half it's bad
    num_contigs = len(set([x["contig_id"] for x in aln_vect]))
    if num_contigs != 1:
        return False,-1

    if len(cycle) == 1 and circular:
        split_ind = -1
        direction = cycle[0][-1]
        first_set = set()
        second_set = set()
        first_label = aln_vect[0]["seg_label"]
        first_set.add(first_label)
        prev = first_label
        direction = "+" if aln_vect[1]["seg_label"] - first_label > 0 else "-"
        for ind,i in enumerate(aln_vect[1:]):
            curr_label = i["seg_label"]
            if not curr_label < prev and direction == "+":
                first_set.add(curr_label)
                prev = curr_label

            elif not curr_label > prev and direction == "-":
                first_set.add(curr_label)
                prev = curr_label
                
            else:
                second_set.add(curr_label)
                if split_ind == -1:
                    split_ind = ind+1

        s1,s2 = sorted([len(first_set),len(second_set)])
        print s1,s2,split_ind,s1/float(s2)
        return (s1/float(s2) > .25),split_ind

    return False,-1
                
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
def adjacent_segs(cycle,segSeqD,isCycle):
    print "checking adjacency"
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


    if isCycle and len(cycle) > 1:
        init_start = segSeqD[cycle[0][0]][2] if cycle[0][1] == "-" else segSeqD[cycle[0][0]][1]
        init_chr = segSeqD[cycle[0][0]][0]
        if p_chrom == curr_chrom and abs(init_start - p_end) == 1:
            prev_seg_index_is_adj[0] = True

    # print prev_seg_index_is_adj
    return prev_seg_index_is_adj

def get_raw_path_length(path,segSeqD,isCycle,prev_seg_index_is_adj):
    raw_path_length= 0.0
    for i in path:
        s_tup = segSeqD[i[0]]
        s_len = s_tup[2] - s_tup[1]
        raw_path_length+=s_len

    return raw_path_length

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

#check contig end trimming
def decide_trim_contigs(contig_cmap_vects,contig_placements,total_length):
    print "DECIDING TRIMMING"
    for cObj in contig_placements.values():
        print cObj.id
        cmap_vect = contig_cmap_vects[cObj.id]
        first_lab,last_lab = cObj.aln_lab_ends

        if (cmap_vect[first_lab-1] - cmap_vect[0])*cObj.scaling_factor > unaligned_cutoff_frac*total_length:
            cObj.start_trim = True
            print "start_trim true"

        if (cmap_vect[-1] - cmap_vect[last_lab-1])*cObj.scaling_factor > unaligned_cutoff_frac*total_length:
            cObj.end_trim = True
            print "end_trim true"

        if cObj.start_trim or cObj.end_trim:
            cObj.trim_obj_ends(total_length)

#TEMP SOLUTION (will break if too many consecutive overlaps)
def set_contig_height_shifts(contig_placements,contig_list,scale_mult=1):
    print "SETTING HEIGHTS"
    prev_offset = 0
    for ind,i in enumerate(contig_list[1:]):
        prevObj = contig_placements[contig_list[ind]]
        currObj = contig_placements[i]

        if currObj.abs_start_pos < prevObj.abs_end_pos:
            shift_mult = -1 if prev_offset == 0 else 0
            currObj.track_height_shift = shift_mult*1.5*scale_mult
            prev_offset = shift_mult

        else:
            prev_offset = 0




def place_path_segs_and_labels(path,ref_placements,seg_cmap_vects):
    path_seg_placements = {}
    for ind,i in enumerate(path):
        refObj = ref_placements[ind]
        segObj = copy.deepcopy(refObj)
        segObj.cmap_vect = seg_cmap_vects[i[0]]
        segObj.compute_label_posns()
        path_seg_placements[ind] = segObj

    return path_seg_placements

#create an object for each contig encoding variables such as position of start and end of contig (absolute ends)
#and positioning of contig labels
def place_contigs_and_labels(path_seg_placements,aln_vect,total_length,contig_cmap_vects,isCycle,circularViz):
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
        #print "placing contigs computation step"
        #print c_id
        cc_vect = contig_cmap_vects[c_id]
        san_f = i_list[0]["seg_aln_number"]
        sal_f = i_list[0]["seg_label"]
        cal_f = i_list[0]["contig_label"]
        san_l = i_list[-1]["seg_aln_number"]
        sal_l = i_list[-1]["seg_label"]
        cal_l = i_list[-1]["contig_label"]
        contig_dir = i_list[0]["contig_dir"]
        #print san_f,sal_f,cal_f
        #print san_l,sal_l,cal_l
        curr_contig_struct = CycleVizElemObj(c_id,contig_dir,None,None,cc_vect)

        #look up aln posns from path_seg_placements
        #look up position of first one
        segObj_start = path_seg_placements[san_f]
        seg_start_l_pos = segObj_start.label_posns[sal_f-1]

        #look up position of last one
        segObj_end = path_seg_placements[san_l]
        seg_end_l_pos = segObj_end.label_posns[sal_l-1]
       
        if seg_end_l_pos < seg_start_l_pos:
            seg_end_l_pos+= total_length

        #catch case where contig is overcircularized (e.g. circular assembly)
        if len(contig_aln_dict) == 1 and isCycle and len(i_list) > 2:
            san_s = i_list[1]["seg_aln_number"]
            segObj_second = path_seg_placements[san_s]
            second_seg_abs_end_pos = segObj_second.abs_end_pos
            if seg_end_l_pos < second_seg_abs_end_pos:
                seg_end_l_pos+=total_length

        #compute scaling
        scaling_factor = 1
        if circularViz:
            print c_id,"comp_scaling"
            scaled_seg_dist  = abs(seg_end_l_pos - seg_start_l_pos)*(1-contig_spacing)
            scaling_factor = scaled_seg_dist/(abs(cc_vect[cal_f-1] - cc_vect[cal_l-1]))
            print seg_start_l_pos,seg_end_l_pos,1-contig_spacing,scaled_seg_dist,total_length
            print scaled_seg_dist,scaling_factor
            #SET CONTIG SCALING FACTOR

        curr_contig_struct.scaling_factor = scaling_factor
        #print scaling_factor,c_id

        if contig_dir == "+":
            abs_start_pos = seg_start_l_pos - (cc_vect[cal_f-1])*scaling_factor
            abs_end_pos = abs_start_pos + (cc_vect[-1])*scaling_factor

        else:
            print "applying scaling to ends"
            abs_start_pos = seg_start_l_pos - (cc_vect[cal_l-1])*scaling_factor
            abs_end_pos = abs_start_pos + (cc_vect[-1])*scaling_factor
            print "now",abs_start_pos,abs_end_pos


        print "SEG PLACEMENT ",c_id
        print abs_start_pos,abs_end_pos
        print seg_start_l_pos,seg_end_l_pos,scaling_factor

        curr_contig_struct.abs_start_pos = abs_start_pos
        curr_contig_struct.abs_end_pos = abs_end_pos

        #SET BOUNDARY ALN POSITIONS FROM TRACK
        curr_contig_struct.aln_bound_posns = (seg_start_l_pos,seg_end_l_pos)

        csl = min(i_list[-1]["contig_label"],i_list[0]["contig_label"])
        cel = max(i_list[-1]["contig_label"],i_list[0]["contig_label"])
        print "CSL/CEL",csl,cel
        print ""
        #SET FIRST AND LAST LABEL ALIGNED IN THE CONTIG
        curr_contig_struct.aln_lab_ends = (csl,cel)
        curr_contig_struct.compute_label_posns()
        contig_span_dict[c_id] = curr_contig_struct

    return contig_span_dict,contig_list

def reduce_path(path,prev_seg_index_is_adj,inds,aln_vect=[]):
    #pass
    print "Reducing path by " + str(inds)
    print path
    left,right = inds
    path = path[left:]
    prev_seg_index_is_adj = prev_seg_index_is_adj[left:]
    prev_seg_index_is_adj[0] = False
    item_nums = [a_d["seg_aln_number"] for a_d in aln_vect]
    left_cut_position = bisect.bisect_left(item_nums,left)
    aln_vect = aln_vect[left_cut_position:]
    if right > 0:
        path = path[:-right]
        prev_seg_index_is_adj = prev_seg_index_is_adj[:-right]
        cut_val = len(path)+left
        item_nums = [a_d["seg_aln_number"] for a_d in aln_vect]
        right_cut_position = bisect.bisect_left(item_nums,cut_val)
        aln_vect = aln_vect[:right_cut_position]
    

    downshift = aln_vect[0]["seg_aln_number"]
    for a_ind, a_d in enumerate(aln_vect):
        aln_vect[a_ind]["seg_aln_number"] = aln_vect[a_ind]["seg_aln_number"] - downshift

    print path
    return path, prev_seg_index_is_adj, aln_vect
