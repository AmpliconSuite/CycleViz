import bisect
from collections import defaultdict
import copy
import os

from intervaltree import IntervalTree
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

matplotlib.use('Agg')

contig_spacing = 1. / 100
unaligned_cutoff_frac = 1. / 60


def cart2pol(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x) / (2. * np.pi) * 360
    return rho, phi


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def polar_series_to_cartesians(line_points, rad_vect):
    return zip(*[pol2cart(r, i) for r, i in zip(rad_vect, line_points)])


def round_to_1_sig(x):
    if x == 0:
        return 0.

    return round(x, -int(np.floor(np.log10(abs(x)))))


class CycleVizElemObj(object):
    def __init__(self, m_id, chrom, ref_start, ref_end, direction, s, t, cmap_vect=[]):
        self.id = m_id
        self.chrom = chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.direction = direction
        self.abs_start_pos = s
        self.abs_end_pos = t
        self.scaling_factor = 1
        self.cmap_vect = cmap_vect
        self.aln_lab_ends = (None, None)
        self.aln_bound_posns = (None, None)
        self.label_posns = []
        self.track_height_shift = 0
        self.start_trim = False
        self.end_trim = False
        self.bed_data = {}

    def compute_label_posns(self):
        if self.direction == "+":
            if self.abs_start_pos is None:
                self.abs_start_pos = self.aln_bound_posns[0] - self.scaling_factor * self.cmap_vect[
                    self.aln_lab_ends[0] - 1]

            if self.abs_end_pos is None:
                self.abs_end_pos = self.aln_bound_posns[1] + self.scaling_factor * (
                        self.cmap_vect[-1] - self.cmap_vect[self.aln_lab_ends[1] - 1])

            for i in self.cmap_vect[:-1]:
                self.label_posns.append(self.scaling_factor * i + self.abs_start_pos)

        else:
            if self.abs_start_pos is None:
                self.abs_start_pos = self.aln_bound_posns[0] - self.scaling_factor * (
                        self.cmap_vect[-1] - self.cmap_vect[self.aln_lab_ends[1] - 1])

            if self.abs_end_pos is None:
                self.abs_end_pos = self.aln_bound_posns[1] + self.scaling_factor * self.cmap_vect[
                    self.aln_lab_ends[0] - 1]

            rev_cmap_vect = [self.cmap_vect[-1] - x for x in self.cmap_vect[::-1][1:]]  # does not include length value
            for i in rev_cmap_vect:
                self.label_posns.append(self.scaling_factor * i + self.abs_start_pos)

            self.label_posns = self.label_posns[::-1]

    def to_string(self):
        return "{}{} | Start: {} | End: {} | scaling {}".format(self.id, self.direction, self.chrom,
                                                                str(self.abs_start_pos),
                                                                str(self.abs_end_pos), str(self.scaling_factor))

    # trim visualized contig if it's long and unaligned
    def trim_obj_ends(self, total_length):
        # use the .end_trim,.start_trim to chop off.
        # overhang should go to 1/3 of the unaligned cutoff threshold
        if self.start_trim:
            p_abs_start = self.abs_start_pos
            print("DOING s TRIM on " + self.id)
            self.abs_start_pos = self.aln_bound_posns[0] - unaligned_cutoff_frac * total_length / 4.
            if self.direction == "-":
                self.update_label_posns(p_abs_start - self.abs_start_pos)
            print("now", self.abs_start_pos)

        if self.end_trim:
            p_abs_end = self.abs_end_pos
            print("DOING e TRIM on " + self.id)
            print(self.aln_bound_posns)
            self.abs_end_pos = self.aln_bound_posns[-1] + unaligned_cutoff_frac * total_length / 4.
            if self.direction == "-":
                self.update_label_posns(p_abs_end - self.abs_end_pos)

            print("now", self.abs_end_pos)

    # update label positions after trimming contigs
    def update_label_posns(self, s_diff):
        print("diff", s_diff)
        for ind in range(len(self.label_posns)):
            self.label_posns[ind] -= s_diff


# makes a gene object from parsed refGene data
class gene(object):
    def __init__(self, gchrom, gstart, gend, gdata):
        self.gchrom = gchrom
        self.gstart = gstart
        self.gend = gend
        self.gname = gdata[-4]
        self.strand = gdata[3]
        estarts = [int(x) for x in gdata[9].rsplit(",") if x]
        eends = [int(x) for x in gdata[10].rsplit(",") if x]
        self.eposns = zip(estarts, eends)
        self.gdrops = []
        self.mdrop_shift = 1.07
        self.gdrops_go_to_link = set()

    def get_angles(self, seg_dir, normStart, normEnd, total_length):
        tm = "X"
        start_angle = normStart / total_length * 360
        end_angle = normEnd / total_length * 360

        if seg_dir == "+" and self.strand == "+":
            # mdrop_shift = 1.07
            s_ang = start_angle
            e_ang = end_angle
            sm = "<"
            em = "s"

        elif seg_dir == "+" and self.strand == "-":
            self.mdrop_shift = 1.05

            s_ang = end_angle
            e_ang = start_angle
            sm = ">"
            em = "s"

        elif seg_dir == "-" and self.strand == "+":
            self.mdrop_shift = 1.05
            s_ang = end_angle
            e_ang = start_angle
            sm = ">"
            em = "s"

        else:
            # mdrop_shift = 1.07
            s_ang = start_angle
            e_ang = end_angle
            sm = "<"
            em = "s"

        return s_ang, e_ang, sm, em, tm

    def draw_marker_ends(self, outer_bar):
        # iterate over gdrops and see how many times the gene appears.
        if self.gdrops:
            self.gdrops = sorted(self.gdrops, key=lambda x: x[-1])
            for gd in self.gdrops:
                normStart, normEnd, total_length, seg_dir, currStart, currEnd, hasStart, hasEnd, seg_ind, drop, pTup = gd
                if hasStart or hasEnd:  # has start
                    s_ang, e_ang, sm, em, tm = self.get_angles(seg_dir, normStart, normEnd, total_length)

                    if hasStart:
                        x_m, y_m = pol2cart(outer_bar - self.mdrop_shift * drop, (s_ang / 360 * 2 * np.pi))
                        t = matplotlib.markers.MarkerStyle(marker=sm)
                        t._transform = t.get_transform().rotate_deg(s_ang - 89)
                        plt.scatter(x_m, y_m, marker=t, s=20, color='k')

                    if hasEnd:
                        x_m, y_m = pol2cart(outer_bar - self.mdrop_shift * drop, (e_ang / 360 * 2 * np.pi))
                        t = matplotlib.markers.MarkerStyle(marker=em)
                        t._transform = t.get_transform().rotate_deg(e_ang - 91)
                        plt.scatter(x_m, y_m, marker=t, s=8, color='k')

    def draw_trunc_spots(self, outer_bar):
        if len(self.gdrops) > 1:
            print(self.gname)
            self.gdrops = sorted(self.gdrops, key=lambda x: x[-1])
            for ind, gd in enumerate(self.gdrops):
                normStart, normEnd, total_length, seg_dir, currStart, currEnd, hasStart, hasEnd, seg_ind, drop, pTup = gd
                if not hasEnd and ind not in self.gdrops_go_to_link:
                    s_ang, e_ang, sm, em, tm = self.get_angles(seg_dir, normStart, normEnd, total_length)
                    x_m, y_m = pol2cart(outer_bar - self.mdrop_shift * drop, (e_ang / 360 * 2 * np.pi))
                    t = matplotlib.markers.MarkerStyle(marker=tm)
                    t._transform = t.get_transform().rotate_deg(e_ang - 91)
                    plt.scatter(x_m, y_m, marker=t, s=12, color='r')

    def draw_seg_links(self, outer_bar, bar_width):
        if len(self.gdrops) > 1:
            print(self.gname)
            self.gdrops = sorted(self.gdrops, key=lambda x: x[-1])
            for ind, gd in enumerate(self.gdrops[1:]):
                normStart, normEnd, total_length, seg_dir, currStart, currEnd, hasStart, hasEnd, seg_ind, drop, pTup = gd
                pgd = self.gdrops[ind]
                pposTup = pgd[-1]
                pseg_ind = pgd[-3]
                if abs(seg_ind - pseg_ind) == 1:
                    print("RS")
                    if (pTup[0] != pposTup[0] or pTup[1] - pposTup[2] > 1):
                        if hasStart or hasEnd and (hasStart, hasEnd) == (pgd[-5], pgd[-4]):
                            continue

                        self.gdrops_go_to_link.add(ind)
                        if seg_dir == "+":
                            start_rad = pgd[1] / total_length * 2 * np.pi
                            end_rad = normStart / total_length * 2 * np.pi
                        else:
                            start_rad = pgd[0] / total_length * 2 * np.pi
                            end_rad = normEnd / total_length * 2 * np.pi

                        mid_rad = (start_rad + end_rad)/2.0
                        if self.mdrop_shift == 1.07:
                            bd_sign = 1

                        else:
                            drop *= self.mdrop_shift
                            bd_sign = -1


                        thetas1 = np.linspace(start_rad, mid_rad, 100)
                        rhos1 = np.linspace(outer_bar-drop, outer_bar-drop+bd_sign*bar_width/2.0, 100)
                        x1, y1 = polar_series_to_cartesians(thetas1, rhos1)

                        thetas2 = np.linspace(mid_rad, end_rad, 100)
                        rhos2 = np.linspace(outer_bar-drop+bd_sign*bar_width/2.0, outer_bar-drop, 100)
                        x2, y2 = polar_series_to_cartesians(thetas2, rhos2)

                        plt.plot(x1, y1, linewidth=1, color='grey')
                        plt.plot(x2, y2, linewidth=1, color='grey')

                elif pTup[0] == pposTup[0] and pTup[1] - pposTup[2] == 1:
                    self.gdrops_go_to_link.add(ind)

# SET COLORS
def get_chr_colors():
    to_add = plt.cm.get_cmap(None, 4).colors[1:]
    # color_vect = ["#ffe8ed","indianred","salmon","burlywood",'#d5b60a',"xkcd:algae",to_add[0],"darkslateblue",
    #              to_add[2],"#017374","#734a65","#bffe28","xkcd:darkgreen","#910951","xkcd:stone",
    #              "xkcd:purpley","xkcd:brown","lavender","darkseagreen","powderblue","#ff073a",to_add[1],
    #              "magenta","plum"]

    color_vect = ["aqua", "rosybrown", "salmon", "bisque", 'goldenrod', "xkcd:algae", to_add[0], "darkslateblue",
                  "yellow", "sienna", "purple", "#bffe28", "xkcd:darkgreen", "#910951", "xkcd:stone",
                  "xkcd:purpley", "xkcd:brown", "lavender", "darkseagreen", "powderblue", "crimson", to_add[1],
                  "fuchsia", "pink"]

    chrnames = [str(i) for i in (list(range(1, 23)))] + ["X", "Y"]
    chromosome_colors = dict(zip(["chr" + i for i in chrnames], color_vect))
    for i in range(len(chrnames)):
        chromosome_colors[chrnames[i]] = color_vect[i]

    return chromosome_colors


# parse the breakpoint graph to indicate for two endpoints if there is an edge.
def parse_BPG(BPG_file):
    bidirectional_edge_dict = defaultdict(set)
    seg_end_pos_d = {}
    seqnum = 0
    with open(BPG_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if not fields:
                continue

            if fields[0] in ["concordant", "discordant"]:
                e_rep = fields[1].rsplit("->")
                start = e_rep[0][:-1]
                end = e_rep[1][:-1]
                bidirectional_edge_dict[start].add(end)
                bidirectional_edge_dict[end].add(start)

            elif fields[0] == "sequence":
                seqnum += 1
                seg_end_pos_d[str(seqnum)] = (fields[1][:-1], fields[2][:-1])

    return bidirectional_edge_dict, seg_end_pos_d


# extract oncogenes from a file.
# Assumes refseq genome name in last column, or get the refseq name from a gff file
def parse_gene_subset_file(gene_list_file, gff=False):
    gene_set = set()
    with open(gene_list_file) as infile:
        for line in infile:
            fields = line.rstrip().split()
            if not fields:
                continue

            if not gff:
                gene_set.add(fields[-1].strip("\""))
            else:
                # parse the line and get the name
                propFields = {x.split("=")[0]: x.split("=")[1] for x in fields[-1].rstrip(";").split(";")}
                gene_set.add(propFields["Name"])

    return gene_set


def parse_genes(ref):
    t = defaultdict(IntervalTree)
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    if ref == "GRCh37" or ref == "hg19":
        refGene_name = "refGene_hg19.txt"
    else:
        refGene_name = "refGene_" + ref + ".txt"

    seenNames = set()
    with open(os.path.join(__location__, refGene_name)) as infile:
        for line in infile:
            fields = line.rsplit("\t")
            currChrom = fields[2]
            if ref == "GRCh37":
                currChrom = currChrom[3:]

            tstart = int(fields[4])
            tend = int(fields[5])
            gname = fields[-4]
            if gname not in seenNames:
                seenNames.add(gname)
                currGene = gene(currChrom, tstart, tend, fields)
                t[currChrom][tstart:tend] = currGene

    return t


# rotate text to be legible on both sides of circle
def correct_text_angle(text_angle):
    if abs(text_angle > 90 and abs(text_angle) < 270):
        text_angle -= 180
        ha = "right"
    else:
        ha = "left"

    return text_angle, ha


# return list of relevant genes sorted by starting position
def rel_genes(chrIntTree, pTup, gene_set=None):
    if gene_set is None:
        gene_set = set()

    currGenes = {}
    chrom = pTup[0]
    overlappingT = chrIntTree[chrom][pTup[1]:pTup[2]]
    gene_set_only = (len(gene_set) == 0)
    for i in overlappingT:
        gObj = i.data
        gname = gObj.gname
        is_other_feature = (gname.startswith("LOC") or gname.startswith("LINC") or gname.startswith("MIR"))
        if gene_set_only:
            gene_set.add(gname)

        if not is_other_feature and gname in gene_set:
            if gname not in currGenes:
                currGenes[gname] = gObj

            # gene appears in file twice, if one is larger, use it. else just use the widest endpoints
            else:
                oldTStart = currGenes[gname].gstart
                oldTEnd = currGenes[gname].gend
                if gObj.gend - gObj.gstart > oldTEnd - oldTStart:
                    currGenes[gname] = gObj

                else:
                    if gObj.gstart < oldTStart:
                        currGenes[gname].gstart = gObj.gstart
                    if gObj.gend > oldTEnd:
                        currGenes[gname].gend = gObj.gend

    relGenes = sorted(currGenes.values(), key=lambda x: (x.gstart, x.gend))
    return relGenes


def pair_is_edge(a_id, b_id, a_dir, b_dir, bpg_dict, seg_end_pos_d):
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
                segSeqD[segNum] = (chrom, lowerBound, upperBound)

            elif "Cycle=" in line:
                isCycle = True
                curr_cycle = []
                fields = line.rstrip().rsplit(";")
                lineD = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in fields}
                segs = lineD["Segments"].rsplit(",")
                for i in segs:
                    seg = i[:-1]
                    if seg != "0" and i:
                        strand = i[-1]
                        curr_cycle.append((seg, strand))

                    else:
                        isCycle = False

                cycles[lineD["Cycle"]] = curr_cycle
                circular_D[lineD["Cycle"]] = isCycle

    return cycles, segSeqD, circular_D


def check_segdup(aln_vect, cycle, circular):
    print("Checking if segdup")
    # iterate over and delete the second half it's bad
    num_contigs = len(set([x["contig_id"] for x in aln_vect]))
    if num_contigs != 1:
        return False, -1

    if len(cycle) == 1 and circular:
        split_ind = -1
        first_set = set()
        second_set = set()
        first_label = aln_vect[0]["seg_label"]
        first_set.add(first_label)
        prev = first_label
        direction = "+" if aln_vect[1]["seg_label"] - first_label > 0 else "-"
        for ind, i in enumerate(aln_vect[1:]):
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
                    split_ind = ind + 1

        s1, s2 = sorted([len(first_set), len(second_set)])
        print(s1, s2, split_ind, s1 / float(s2))
        return (s1 / float(s2) > .25), split_ind

    return False, -1


def parse_alnfile(path_aln_file):
    aln_vect = []
    with open(path_aln_file) as infile:
        # read a few special header lines directly by calling .next(). This will not work in python3!
        # Ideally there should be a way to read the next line from infile with a command that works in 2 & 3.
        meta_header = next(infile).rstrip()[1:].split()
        aln_metadata_fields = next(infile).rstrip()[1:].split()
        meta_dict = dict(zip(meta_header, aln_metadata_fields))
        aln_header = next(infile).rstrip()[1:].split()
        for line in infile:
            fields = line.rstrip().split()
            fields_dict = dict(zip(aln_header, fields))
            fields_dict["contig_label"] = int(fields_dict["contig_label"])
            fields_dict["seg_label"] = int(fields_dict["seg_label"])
            fields_dict["seg_aln_number"] = int(fields_dict["seg_aln_number"])
            aln_vect.append(fields_dict)

    return aln_vect, meta_dict


# determine segments linearly adjacent in ref genome
def adjacent_segs(cycle, segSeqD, isCycle):
    print("checking adjacency")
    prev_seg_index_is_adj = [False] * len(cycle)
    p_end = segSeqD[cycle[0][0]][2] if cycle[0][1] == "+" else segSeqD[cycle[0][0]][1]
    p_chrom = segSeqD[cycle[0][0]][0]
    for ind in range(1, len(cycle)):
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


def get_raw_path_length(path, segSeqD):
    raw_path_length = 0.0
    for i in path:
        s_tup = segSeqD[i[0]]
        s_len = s_tup[2] - s_tup[1]
        raw_path_length += s_len

    return raw_path_length


# segment is imputed by AR or not
def imputed_status_from_aln(aln_vect, cycle_len):
    imputed_status = [int(aln_vect[0]["imputed"])]
    curr_seg_aln_number = 0
    for a_d in aln_vect:
        if a_d["seg_aln_number"] != curr_seg_aln_number:
            for i in range(curr_seg_aln_number + 1, a_d["seg_aln_number"]):
                imputed_status.append(1)

            imputed_status.append(int(a_d["imputed"]))
            curr_seg_aln_number = a_d["seg_aln_number"]

    for i in range(curr_seg_aln_number + 1, cycle_len):
        imputed_status.append(1)

    return imputed_status


# check contig end trimming
def decide_trim_contigs(contig_cmap_vects, contig_placements, total_length):
    print("DECIDING TRIMMING")
    for cObj in contig_placements.values():
        print(cObj.id)
        cmap_vect = contig_cmap_vects[cObj.id]
        first_lab, last_lab = cObj.aln_lab_ends

        if (cmap_vect[first_lab - 1] - cmap_vect[0]) * cObj.scaling_factor > unaligned_cutoff_frac * total_length:
            cObj.start_trim = True
            print("start_trim true")

        if (cmap_vect[-1] - cmap_vect[last_lab - 1]) * cObj.scaling_factor > unaligned_cutoff_frac * total_length:
            cObj.end_trim = True
            print("end_trim true")

        if cObj.start_trim or cObj.end_trim:
            cObj.trim_obj_ends(total_length)


# TEMP SOLUTION (will break if too many consecutive overlaps)
def set_contig_height_shifts(contig_placements, contig_list, scale_mult=1):
    print("SETTING HEIGHTS")
    prev_offset = 0
    for ind, i in enumerate(contig_list[1:]):
        prevObj = contig_placements[contig_list[ind]]
        currObj = contig_placements[i]

        if currObj.abs_start_pos < prevObj.abs_end_pos:
            shift_mult = -1 if prev_offset == 0 else 0
            currObj.track_height_shift = shift_mult * 1.5 * scale_mult
            prev_offset = shift_mult

        else:
            prev_offset = 0


def place_path_segs_and_labels(path, ref_placements, seg_cmap_vects):
    path_seg_placements = {}
    for ind, i in enumerate(path):
        refObj = ref_placements[ind]
        segObj = copy.deepcopy(refObj)
        segObj.cmap_vect = seg_cmap_vects[i[0]]
        segObj.compute_label_posns()
        path_seg_placements[ind] = segObj

    return path_seg_placements


# create an object for each contig encoding variables such as position of start and end of contig (absolute ends)
# and positioning of contig labels
def place_contigs_and_labels(path_seg_placements, aln_vect, total_length, contig_cmap_vects, isCycle, circularViz,
                             segSeqD):
    contig_aln_dict = defaultdict(list)
    contig_list = []
    for i in aln_vect:
        c_id = i["contig_id"]
        contig_aln_dict[c_id].append(i)
        if c_id not in contig_list: contig_list.append(c_id)

    contig_span_dict = {}
    for c_id, i_list in contig_aln_dict.items():
        # print "placing contigs computation step"
        # print(c_id)
        cc_vect = contig_cmap_vects[c_id]
        san_f = i_list[0]["seg_aln_number"]
        sal_f = i_list[0]["seg_label"]
        cal_f = i_list[0]["contig_label"]
        san_l = i_list[-1]["seg_aln_number"]
        sal_l = i_list[-1]["seg_label"]
        cal_l = i_list[-1]["contig_label"]
        contig_dir = i_list[0]["contig_dir"]
        # print(san_f,sal_f,cal_f)
        # print(san_l,sal_l,cal_l)
        curr_contig_struct = CycleVizElemObj(c_id, segSeqD[c_id[0]][0], segSeqD[c_id[0]][1], segSeqD[c_id[0]][2],
                                             contig_dir, None, None, cc_vect)

        # look up aln posns from path_seg_placements
        # look up position of first one
        segObj_start = path_seg_placements[san_f]
        seg_start_l_pos = segObj_start.label_posns[sal_f - 1]

        # look up position of last one
        segObj_end = path_seg_placements[san_l]
        seg_end_l_pos = segObj_end.label_posns[sal_l - 1]

        if seg_end_l_pos < seg_start_l_pos:
            seg_end_l_pos += total_length

        # catch case where contig is overcircularized (e.g. circular assembly)
        if len(contig_aln_dict) == 1 and isCycle and len(i_list) > 2:
            san_s = i_list[1]["seg_aln_number"]
            segObj_second = path_seg_placements[san_s]
            second_seg_abs_end_pos = segObj_second.abs_end_pos
            if seg_end_l_pos < second_seg_abs_end_pos:
                seg_end_l_pos += total_length

        # compute scaling
        scaling_factor = 1
        if circularViz:
            print(c_id, "comp_scaling")
            scaled_seg_dist = abs(seg_end_l_pos - seg_start_l_pos) * (1 - contig_spacing)
            scaling_factor = scaled_seg_dist / (abs(cc_vect[cal_f - 1] - cc_vect[cal_l - 1]))
            print(seg_start_l_pos, seg_end_l_pos, 1 - contig_spacing, scaled_seg_dist, total_length)
            print(scaled_seg_dist, scaling_factor)
            # SET CONTIG SCALING FACTOR

        curr_contig_struct.scaling_factor = scaling_factor
        # print scaling_factor,c_id

        if contig_dir == "+":
            abs_start_pos = seg_start_l_pos - (cc_vect[cal_f - 1]) * scaling_factor
            abs_end_pos = abs_start_pos + (cc_vect[-1]) * scaling_factor

        else:
            print("applying scaling to ends")
            abs_start_pos = seg_start_l_pos - (cc_vect[cal_l - 1]) * scaling_factor
            abs_end_pos = abs_start_pos + (cc_vect[-1]) * scaling_factor
            print("now", abs_start_pos, abs_end_pos)

        print("SEG PLACEMENT ", c_id)
        print(abs_start_pos, abs_end_pos)
        print(seg_start_l_pos, seg_end_l_pos, scaling_factor)

        curr_contig_struct.abs_start_pos = abs_start_pos
        curr_contig_struct.abs_end_pos = abs_end_pos

        # SET BOUNDARY ALN POSITIONS FROM TRACK
        curr_contig_struct.aln_bound_posns = (seg_start_l_pos, seg_end_l_pos)

        csl = min(i_list[-1]["contig_label"], i_list[0]["contig_label"])
        cel = max(i_list[-1]["contig_label"], i_list[0]["contig_label"])
        print("CSL/CEL", csl, cel)
        print("")
        # SET FIRST AND LAST LABEL ALIGNED IN THE CONTIG
        curr_contig_struct.aln_lab_ends = (csl, cel)
        curr_contig_struct.compute_label_posns()
        contig_span_dict[c_id] = curr_contig_struct

    return contig_span_dict, contig_list


def reduce_path(path, prev_seg_index_is_adj, inds, aln_vect=None):
    if aln_vect is None:
        aln_vect = []

    print("Reducing path by " + str(inds))
    print(path)
    left, right = inds
    path = path[left:]
    prev_seg_index_is_adj = prev_seg_index_is_adj[left:]
    prev_seg_index_is_adj[0] = False
    item_nums = [a_d["seg_aln_number"] for a_d in aln_vect]
    left_cut_position = bisect.bisect_left(item_nums, left)
    aln_vect = aln_vect[left_cut_position:]
    if right > 0:
        path = path[:-right]
        prev_seg_index_is_adj = prev_seg_index_is_adj[:-right]
        cut_val = len(path) + left
        item_nums = [a_d["seg_aln_number"] for a_d in aln_vect]
        right_cut_position = bisect.bisect_left(item_nums, cut_val)
        aln_vect = aln_vect[:right_cut_position]

    if aln_vect:
        downshift = aln_vect[0]["seg_aln_number"]
        for a_ind, a_d in enumerate(aln_vect):
            aln_vect[a_ind]["seg_aln_number"] = aln_vect[a_ind]["seg_aln_number"] - downshift

    print(path)
    return path, prev_seg_index_is_adj, aln_vect


def parse_yaml(args):
    import yaml
    with open(args.yaml_file) as f:
        args.bedgraph = ''
        sample_data = yaml.safe_load(f)
        args.cycles_file = sample_data.get("cycles_file")
        print(args.cycles_file)
        args.cycle = str(sample_data.get("cycle"))
        if "om_alignments" in sample_data:
            args.om_alignments = sample_data.get("om_alignments")
        if "contigs" in sample_data:
            args.contigs = sample_data.get("contigs")
        if "segs" in sample_data:
            args.segs = sample_data.get("segs")
        if "graph" in sample_data:
            args.graph = sample_data.get("graph")
        if "i" in sample_data:
            args.path_alignment = sample_data.get("i")
        if "ref" in sample_data:
            args.ref = sample_data.get("ref")
        if "sname" in sample_data:
            args.sname = sample_data.get("sname")
        if "rot" in sample_data:
            args.rot = sample_data.get("rot")
        if "label_segs" in sample_data:
            args.label_segs = sample_data.get("label_segs")
        if "gene_subset_file" in sample_data:
            args.gene_subset_files = sample_data.get("gene_subset_file")
        if "gene_subset_list" in sample_data:
            args.gene_subset_list = sample_data.get("gene_subset_list")
            print(args.gene_subset_list)
        if "print_dup_genes" in sample_data:
            args.print_dup_genes = sample_data.get("print_dup_genes")
        if "gene_fontsize" in sample_data:
            args.gene_fontsize = sample_data.get("gene_fontsize")
        if "tick_fontsize" in sample_data:
            args.tick_fontsize = sample_data.get("tick_fontsize")
        if "bedgraph_file" in sample_data:
            args.bedgraph = sample_data.get("begraph_file")

    return args
