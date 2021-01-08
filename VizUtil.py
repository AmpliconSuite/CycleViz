import bisect
from collections import defaultdict
import copy
import os
import sys

from intervaltree import IntervalTree
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import yaml

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


# arguments: thetas, rhos, however rad_vect can be a single rho, and will get turned into a list of same length as theta
def polar_series_to_cartesians(line_points, rad_vect):
    if len(line_points) == 0:
        return [], []

    if not isinstance(rad_vect, list):
        rad_vect = [rad_vect,] * len(line_points)

    return zip(*[pol2cart(r, i) for r, i in zip(rad_vect, line_points)])


def round_to_1_sig(x):
    if x == 0:
        return 0.

    return round(x, -int(np.floor(np.log10(abs(x)))))


def convert_gpos_to_ropos(cgpos, ro_start, ro_end, g_start, dir):
    if dir == "+":
        return (cgpos - g_start) + ro_start
    else:
        return ro_end - (cgpos - g_start)


class CycleVizElemObj(object):
    def __init__(self, m_id, chrom, ref_start, ref_end, direction, s, t, seg_count, padj, nadj, cmap_vect=None):
        if cmap_vect is None: cmap_vect = []
        self.id = m_id
        self.chrom = chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.direction = direction
        self.abs_start_pos = s
        self.abs_end_pos = t
        self.scaling_factor = 1
        self.seg_count = seg_count
        self.cmap_vect = cmap_vect
        self.aln_lab_ends = (None, None)
        self.aln_bound_posns = (None, None)
        self.label_posns = []
        self.track_height_shift = 0
        self.start_trim = False
        self.end_trim = False
        self.feature_tracks = []
        self.prev_is_adjacent = padj
        self.next_is_adjacent = nadj
        self.custom_color = None
        self.custom_bh = None

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


# this stores the local properties of each gene's visualization
class gene_viz_instance(object):
    def __init__(self, gParent, normStart, normEnd, total_length, seg_dir, currStart, currEnd, hasStart, hasEnd, seg_ind, pTup):
        self.gParent = gParent
        self.normStart = normStart
        self.normEnd = normEnd
        self.total_length = total_length
        self.seg_dir = seg_dir
        self.currStart = currStart
        self.currEnd = currEnd
        self.hasStart = hasStart
        self.hasEnd = hasEnd
        self.seg_ind = seg_ind
        self.pTup = pTup

    def get_angles(self):
        tm = "X"
        start_angle = self.normStart / self.total_length * 360
        end_angle = self.normEnd / self.total_length * 360

        if self.seg_dir == "+" and self.gParent.strand == "+":
            s_ang = start_angle
            e_ang = end_angle
            sm = "<"
            em = "s"

        elif self.seg_dir == "+" and self.gParent.strand == "-":
            s_ang = end_angle
            e_ang = start_angle
            sm = ">"
            em = "s"

        elif self.seg_dir == "-" and self.gParent.strand == "+":
            s_ang = end_angle
            e_ang = start_angle
            sm = ">"
            em = "s"

        else:
            s_ang = start_angle
            e_ang = end_angle
            sm = "<"
            em = "s"

        return s_ang, e_ang, sm, em, tm

    def draw_marker_ends(self, gbh):
        # iterate over gdrops and see how many times the gene appears.
        # self.gdrops = sorted(self.gdrops, key=lambda x: x[-1])
        if self.hasStart or self.hasEnd:
            s_ang, e_ang, sm, em, tm = self.get_angles()

            if self.hasStart:
                x_m, y_m = pol2cart(gbh, (s_ang / 360 * 2 * np.pi))
                t = matplotlib.markers.MarkerStyle(marker=sm)
                t._transform = t.get_transform().rotate_deg(s_ang - 89)
                plt.scatter(x_m, y_m, marker=t, s=15, color='silver',zorder=3,alpha=0.8)

            if self.hasEnd:
                x_m, y_m = pol2cart(gbh, (e_ang / 360 * 2 * np.pi))
                t = matplotlib.markers.MarkerStyle(marker=em)
                t._transform = t.get_transform().rotate_deg(e_ang - 91)
                plt.scatter(x_m, y_m, marker=t, s=5, color='silver',zorder=3,alpha=0.8)


# makes a gene object from parsed refGene data
# this stores global properties for the gene
class gene(object):
    def __init__(self, gchrom, gstart, gend, gdata, highlight_name):
        self.gchrom = gchrom
        self.gstart = gstart
        self.gend = gend
        self.gname = gdata[-4]
        self.strand = gdata[3]
        self.highlight_name = highlight_name
        estarts = [int(x) for x in gdata[9].rsplit(",") if x]
        eends = [int(x) for x in gdata[10].rsplit(",") if x]
        self.eposns = zip(estarts, eends)
        self.gdrops = []
        self.gdrops_go_to_link = set()


class feature_track(object):
    def __init__(self, index, primary_data, secondary_data, dd, dv_min, dv_max):
        self.index = index
        self.primary_data = primary_data
        self.secondary_data = secondary_data
        self.track_props = dd
        self.track_min = dv_min
        self.track_max = dv_max
        self.sec_rsf = 1
        self.sec_rss = 0
        self.minsec = 0
        self.base = 0
        self.top = 0
        self.primary_links = []
        self.secondary_links = []

    class Link(object):
        def __init__(self, chromA, chromB, data_tup):
            self.chromA = chromA
            self.chromB = chromB
            self.startA, self.endA = data_tup[0], data_tup[1]
            self.startB, self.endB = data_tup[2], data_tup[3]
            self.score = data_tup[4][0]
            self.link_color = data_tup[4][1]
            self.posA_hits = []
            self.posB_hits = []


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


# rotate text to be legible on both sides of circle
def correct_text_angle(text_angle):
    if abs(text_angle > 90 and abs(text_angle) < 270):
        text_angle -= 180
        ha = "right"
    else:
        ha = "left"

    return text_angle, ha


def pair_is_edge(a_id, b_id, a_dir, b_dir, bpg_dict, seg_end_pos_d):
    rObj1_end = seg_end_pos_d[a_id][-1] if a_dir == "+" else seg_end_pos_d[a_id][0]
    rObj2_start = seg_end_pos_d[b_id][0] if b_dir == "+" else seg_end_pos_d[b_id][-1]
    return rObj1_end in bpg_dict[rObj2_start]


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


# Read an AA-formatted cycles file
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


def handle_struct_bed_data(struct_data):
    #sort by the index
    flatvals = [[chrom] + list(x) for chrom, cl in struct_data.items() for x in cl]
    flatvals.sort(key=lambda x: x[3][0])
    segSeqD = {}
    rev_segSeqD = {}
    seg_end_pos_d = {}
    bidirectional_edge_dict = defaultdict(set)
    cycle = []
    connections = []
    isCycle = True
    uind = 0
    for ind, i in enumerate(flatvals):
        print(i)
        # 0 is chrom, 1 is index, 2 is start, 3 is end, 4 is tuple of datastuff
        # datastuff is strand, connected
        pt = (i[0], i[1], i[2])
        if pt not in segSeqD.values():
            uind += 1
            segSeqD[uind] = pt
            rev_segSeqD[pt] = uind

        cind = rev_segSeqD[pt]
        cycle.append((cind, i[3][1]))
        seg_end_pos_d[cind] = (pt[0] + ':' + str(pt[1]), pt[0] + ':' + str(pt[2]))
        connections.append(i[3][2] == 'True')  # expecting 'True' or 'False' string in this column from the file

    # if non-cyclic path, put zeros on it
    if not all(connections):
        # put zeros on it -- no, those are ordinarily stripped!
        # cycle = [(0, '+')] + cycle + [(0, '+')]
        isCycle = False

    # make the bpg dict
    print(cycle)
    for a, b, conn in zip(cycle, cycle[1:] + [cycle[0]], connections):
        print(a,b,conn)
        if conn:
            ae = seg_end_pos_d[a[0]][1] if a[1] == '+' else seg_end_pos_d[a[0]][0]
            bs = seg_end_pos_d[b[0]][0] if b[1] == '+' else seg_end_pos_d[b[0]][1]
            bidirectional_edge_dict[ae].add(bs)
            bidirectional_edge_dict[bs].add(ae)

    print(bidirectional_edge_dict)
    return cycle, isCycle, segSeqD, seg_end_pos_d, bidirectional_edge_dict


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
                    currGenes[gname] = copy.copy(gObj)

                else:
                    if gObj.gstart < oldTStart:
                        currGenes[gname].gstart = gObj.gstart
                    if gObj.gend > oldTEnd:
                        currGenes[gname].gend = gObj.gend

    relGenes = sorted(currGenes.values(), key=lambda x: (x.gstart, x.gend))
    return relGenes


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


def parse_genes(ref, gene_highlight_list):
    t = defaultdict(IntervalTree)
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    if ref == "GRCh37" or ref == "hg19":
        refGene_name = "refGene_hg19.txt"
    else:
        refGene_name = "refGene_" + ref + ".txt"

    seenNames = set()
    with open(os.path.join(__location__, "resources", refGene_name)) as infile:
        for line in infile:
            fields = line.rsplit("\t")
            currChrom = fields[2]
            if ref == "GRCh37" and not currChrom.startswith("hpv"):
                currChrom = currChrom[3:]

            tstart = int(fields[4])
            tend = int(fields[5])
            gname = fields[-4]
            if gname not in seenNames:
                seenNames.add(gname)
                currGene = gene(currChrom, tstart, tend, fields, gname in gene_highlight_list)
                t[currChrom][tstart:tend] = currGene

    return t


def parse_bed(bedfile, store_all_additional_fields=False):
    data_dict = defaultdict(list)
    with open(bedfile) as infile:
        if bedfile.endswith(".bedpe"):
            for line in infile:
                line = line.rstrip()
                if not line.startswith("#") and line:
                    fields = line.rsplit("\t")
                    chromA = fields[0]
                    startA, endA = float(fields[1]), float(fields[2])
                    chromB = fields[3]
                    startB, endB = float(fields[4]), float(fields[5])
                    data = (float(fields[6]), fields[7])
                    data_dict[(chromA, chromB)].append((startA, endA, startB, endB, data))

        else:
            for entry_index, line in enumerate(infile):
                line = line.rstrip("\t")
                if not line.startswith("#") and line:
                    fields = line.rsplit()
                    chrom = fields[0]
                    begin, end = float(fields[1]), float(fields[2])
                    if not store_all_additional_fields:
                        if len(fields) > 3:
                            data = float(fields[-1])
                        else:
                            data = None
                    else:
                        data = tuple([entry_index] + fields[3:])

                    data_dict[chrom].append((begin, end, data))

    return data_dict


def rescale_by_secondary(primary_dset, secondary_dset, chrom, mode):
    # put secondary data into an intervaltree
    if mode == True:
        mode = "mean"

    if mode != "mean" and mode != "each":
        print("Incorrect norm by secondary mode selected, must be 'mean' or 'each'... using 'mean'")

    if len(secondary_dset) == 0:
        print("No secondary data! skipping normalization")
        return primary_dset, secondary_dset

    sit = IntervalTree()
    normed_primary = defaultdict(list)
    normed_secondary = defaultdict(list)
    if mode == "each":
        for point in secondary_dset[chrom]:
            sit.addi(point[0], point[1], point[2])
            normed_secondary[chrom].append([point[0], point[1], 2])

    elif mode == "mean":
        print("Normalizing 'mean' for secondary will update secondary")
        lscale = 100000.0
        runl = 0.
        runs = 0.
        c = 0
        #for everything in secondary, compute a mean
        for ival in secondary_dset.values():
            for point in ival:
                c+=1
                l = (point[1] - point[0])/lscale
                runl+=l
                runs+=(l*point[2])

        if c > 0:
            allmean = runs/runl
        else:
            print("no secondary data, setting scale to 1")
            allmean = 1.0

        #replace everything in secondary with that mean
        for point in secondary_dset[chrom]:
            sit.addi(point[0], point[1], allmean)
            normed_secondary[chrom].append([point[0], point[1], 2.0])


    for point in primary_dset[chrom]:
        hit_sec = list(sit[point[0]:point[1]])
        if not hit_sec:
            print("could not normalize " + str(point))

        elif len(hit_sec) > 1:
            print(str(point) + ": multiple secondary track hits for normalization, using first hit " + "(" + str(hit_sec[0]) + ")")

        else:
            normed_primary[chrom].append([point[0], point[1], point[2]/float(hit_sec[0].data)])

    return normed_primary, normed_secondary


# take the feature data (cfc) and extract only the regions overlapping the reference segment in question (obj)
# append the coordinate restricted feature (restricted_cfc) to a list of features kept by the reference object (obj)
def store_bed_data(cfc, ref_placements, primary_end_trim=0, secondary_end_trim=0):
    if cfc.track_props['tracktype'] == 'standard' or cfc.track_props['tracktype'] == 'rects':
        print("extracting features, ET", primary_end_trim)
        for obj in ref_placements.values():
            primeTrim = primary_end_trim
            if obj.ref_end - obj.ref_start <= primary_end_trim*2:
                primeTrim = max(0,(obj.ref_end - obj.ref_start)/2 - 2)
                print("reset ET ", primeTrim)

            secTrim = secondary_end_trim
            if obj.ref_end - obj.ref_start <= primary_end_trim * 2:
                secTrim = max(0, (obj.ref_end - obj.ref_start) / 2 - 2)
                print("reset ET ", secTrim)

            local_primary_data = defaultdict(list)
            local_secondary_data = defaultdict(list)

            #store primary data
            for dstore, currdata, currTrim, uc, lc, in zip([local_primary_data, local_secondary_data],
                                        [cfc.primary_data[obj.chrom], cfc.secondary_data[obj.chrom]],
                                        [primeTrim, secTrim], [cfc.track_props['primary_upper_cap'],
                                                               cfc.track_props['secondary_upper_cap']],
                                                           [cfc.track_props['primary_lower_cap'],
                                                            cfc.track_props['secondary_lower_cap']]):

                for init_point in currdata:
                    #open the interval by 1
                    point = [init_point[0], init_point[1]+1, init_point[2]]
                    if obj.ref_start+currTrim <= point[0] <= obj.ref_end-currTrim or \
                            obj.ref_start+currTrim <= point[1] <= obj.ref_end-currTrim:

                        if uc and point[2] > uc:
                            point[2] = uc

                        if lc and point[2] < lc:
                            point[2] = lc

                        dstore[obj.chrom].append(point)

                    elif point[0] < obj.ref_start+currTrim and point[1] > obj.ref_end-currTrim:
                        dstore[obj.chrom].append(point)


            restricted_cfc = copy.copy(cfc)
            if cfc.track_props['rescale_by_secondary']:
                print(obj.to_string(), "normalizing by secondary")
                normed_primary, normed_secondary = rescale_by_secondary(local_primary_data, local_secondary_data,
                                                                    obj.chrom, cfc.track_props['rescale_by_secondary'])
                restricted_cfc.primary_data = normed_primary
                restricted_cfc.secondary_data = normed_secondary

            elif cfc.track_props['rescale_by_count']:
                print(obj.to_string(), "recaling by count")
                normed_primary = defaultdict(list)
                for point in local_primary_data[obj.chrom]:
                    normed_primary[obj.chrom].append([point[0], point[1], point[2] / float(obj.seg_count)])

                restricted_cfc.primary_data = normed_primary
                restricted_cfc.secondary_data = local_secondary_data

            else:
                restricted_cfc.primary_data = local_primary_data
                restricted_cfc.secondary_data = local_secondary_data

            obj.feature_tracks.append(restricted_cfc)
            print(obj.to_string(),"now has",len(obj.feature_tracks),"tracks")

    elif cfc.track_props['tracktype'] == 'links':
        for cdat, c_link_store in zip([cfc.primary_data, cfc.secondary_data], [cfc.primary_links, cfc.secondary_links]):
            for cp, data_tup_list in cdat.items():
                for data_tup in data_tup_list:
                    cLink = cfc.Link(cp[0], cp[1], data_tup)
                    # hitA, hitB = False, False
                    for obj in ref_placements.values():
                        seg_len = obj.ref_end - obj.ref_start
                        # check if it's chromA
                        if obj.chrom == cLink.chromA:
                            # check if overlap coord
                            if obj.ref_start <= cLink.startA <= obj.ref_end or \
                                    obj.ref_start <= cLink.endA <= obj.ref_end:
                                # compute the normpos
                                # print("matched", cp, data_tup, "A")
                                if obj.direction == "+":
                                    normStart_A = obj.abs_start_pos + max(0, cLink.startA - obj.ref_start)
                                    normEnd_A = obj.abs_start_pos + min(seg_len, cLink.endA - obj.ref_start)

                                else:
                                    normEnd_A = obj.abs_start_pos + min(seg_len, obj.ref_end - cLink.startA)
                                    normStart_A = obj.abs_start_pos + max(0, obj.ref_end - cLink.endA)

                                cLink.posA_hits.append((normStart_A, normEnd_A))
                                # hitA = True

                        # check if it's chromB
                        if obj.chrom == cLink.chromB:
                            # check if overlap coord
                            if obj.ref_start <= cLink.startB <= obj.ref_end or \
                                    obj.ref_start <= cLink.endB <= obj.ref_end:
                                # print("matched", cp, data_tup, "B")
                                # compute the normpos
                                if obj.direction == "+":
                                    normStart_B = obj.abs_start_pos + max(0, cLink.startB - obj.ref_start)
                                    normEnd_B = obj.abs_start_pos + min(seg_len, cLink.endB - obj.ref_start)

                                else:
                                    normEnd_B = obj.abs_start_pos + min(seg_len, obj.ref_end - cLink.startB)
                                    normStart_B = obj.abs_start_pos + max(0, obj.ref_end - cLink.endB)

                                cLink.posB_hits.append((normStart_B, normEnd_B))
                                # hitB = True

                        # if cfc.track_props['link_single_match'] and hitA and hitB:
                        #     break

                    c_link_store.append(cLink)

        # holds the place of the feature
        for obj in ref_placements.values():
            obj.feature_tracks.append(cfc)
        #     print(obj.to_string(), "now has", len(obj.feature_tracks), "tracks")


# Assign interior track data to a location(s) in a refobj, and create a new refobj for that overlapping region,
# to represent the interior track data

# TODO: What if the segments in the inside track are scrambled w.r.t the reference track?
def handle_IS_data(ref_placements, IS_cycle, IS_segSeqD, IS_isCircular, IS_bh, cycleColor='lightskyblue'):
    cycle_seg_counts = get_seg_amplicon_count(IS_cycle)
    prev_seg_index_is_adj, next_seg_index_is_adj = adjacent_segs(IS_cycle, IS_segSeqD, IS_isCircular)
    # need to match to the ref_placements
    valid_link_pairs = set(zip(IS_cycle[:-1],IS_cycle[1:]))
    if IS_isCircular:
        valid_link_pairs.add((IS_cycle[-1],IS_cycle[0]))

    new_IS_cycle = []
    new_IS_links = []
    # lastIndHit = -1
    IS_rObj_placements = {}
    u_ind = 0
    for obj in ref_placements.values():
        s_to_add = []
        c_to_add = []
        seg_len = obj.ref_end - obj.ref_start
        for ind, (segID, segdir) in enumerate(IS_cycle):
            padj, nadj = prev_seg_index_is_adj[ind], next_seg_index_is_adj[ind]
            seg_id_count = cycle_seg_counts[segID]
            c, s, e = IS_segSeqD[segID]
            if c != obj.chrom:
                continue

            # overlaps the ref seg
            if obj.ref_start <= s <= obj.ref_end or obj.ref_start <= e <= obj.ref_end:
                # compute the normpos
                if obj.direction == "+":
                    normStart = obj.abs_start_pos + max(0, s - obj.ref_start)
                    normEnd = obj.abs_start_pos + min(seg_len, e - obj.ref_start)

                else:
                    normEnd = obj.abs_start_pos + min(seg_len, obj.ref_end - s)
                    normStart = obj.abs_start_pos + max(0, obj.ref_end - e)

                # make a refobj
                currObj = CycleVizElemObj(segID, c, s, e, obj.direction, normStart, normEnd, seg_id_count, padj, nadj)
                currObj.custom_color = cycleColor
                currObj.custom_bh = IS_bh
                s_to_add.append(currObj)
                c_to_add.append((segID, segdir))

        if not s_to_add:
            continue

        ssta, scta = zip(*sorted(zip(s_to_add, c_to_add), key=lambda x: x[0].abs_start_pos))
        new_IS_cycle.extend(scta)
        for rObj, ct in zip(ssta, scta):
            IS_rObj_placements[u_ind] = rObj
            u_ind+=1

    for a, b in zip(new_IS_cycle[:-1], new_IS_cycle[1:]):
        if (a,b) in valid_link_pairs or (b,a) in valid_link_pairs:
            new_IS_links.append(True)

        else:
            new_IS_links.append(False)

    if (new_IS_cycle[-1], new_IS_cycle[0]) in valid_link_pairs or (new_IS_cycle[0], new_IS_cycle[-1]) in valid_link_pairs:
        new_IS_links.append(True)
    else:
        new_IS_links.append(False)

    print(new_IS_links)
    print(new_IS_cycle)
    return IS_rObj_placements, new_IS_cycle, new_IS_links


# BIONANO PLOTTING FUNCTIONS
# -----------------------------------------


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
        return (s1 / float(s2) > .25), split_ind

    return False, -1


# for use with bionano data & AR output
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

# -----------------------------------------

# determine segments linearly adjacent in ref genome
def adjacent_segs(cycle, segSeqD, isCycle):
    print("checking adjacency")
    prev_seg_index_is_adj = [False] * len(cycle)
    next_seg_index_is_adj = [False] * len(cycle)
    p_end = segSeqD[cycle[0][0]][2] if cycle[0][1] == "+" else segSeqD[cycle[0][0]][1]
    p_chrom = segSeqD[cycle[0][0]][0]
    p_dir = cycle[0][1]
    for ind in range(1, len(cycle)):
        i = cycle[ind]
        curr_chrom = segSeqD[i[0]][0]
        curr_start = segSeqD[i[0]][2] if i[1] == "-" else segSeqD[i[0]][1]
        if curr_chrom == p_chrom and abs(curr_start - p_end) == 1 and p_dir == i[1]:
                prev_seg_index_is_adj[ind] = True
                next_seg_index_is_adj[ind-1] = True

        p_end = segSeqD[i[0]][2] if i[1] == "+" else segSeqD[i[0]][1]
        p_chrom = curr_chrom
        p_dir = i[1]

    if isCycle and len(cycle) > 1:
        init_start = segSeqD[cycle[0][0]][2] if cycle[0][1] == "-" else segSeqD[cycle[0][0]][1]
        init_chr = segSeqD[cycle[0][0]][0]
        if p_chrom == curr_chrom and abs(init_start - p_end) == 1 and p_dir == cycle[0][1]:
            prev_seg_index_is_adj[0] = True
            next_seg_index_is_adj[len(cycle)-1] = True

    # print prev_seg_index_is_adj
    return prev_seg_index_is_adj, next_seg_index_is_adj


# count the number of occurences of BPG segments in the cycle. Store in a dict.
def get_seg_amplicon_count(cycle):
    cycle_id_countd = defaultdict(int)
    for x, _ in cycle:
        cycle_id_countd[x]+=1

    return cycle_id_countd


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
                                             contig_dir, None, None, False, False, cc_vect)

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


def reset_track_min_max(ref_placements, index, primary_cfc):
    tmin, tmax = 0, 0
    for obj in ref_placements.values():
        cfc = obj.feature_tracks[index]
        hs = cfc.track_props['hide_secondary']
        if cfc.track_props['hide_secondary'] == "viral" and not (obj.chrom.startswith('chr') or len(obj.chrom) < 3):
            hs = True

        elif cfc.track_props['hide_secondary'] == "viral":
            hs = False

        curr_track_min, curr_track_max = track_min_max(cfc.primary_data, cfc.secondary_data, True,
                                                       hide_secondary=hs)
        tmin = min(tmin, curr_track_min)
        tmax = max(tmax, curr_track_max)
        if cfc.track_props['show_segment_copy_count']:
            tmin = min(tmin, obj.seg_count)
            tmax = max(tmax, obj.seg_count)

    for obj in ref_placements.values():
        obj.feature_tracks[index].track_min = tmin
        obj.feature_tracks[index].track_max = tmax

    primary_cfc.track_min = tmin
    primary_cfc.track_max = tmax


# go over the bedgraph data and find min and max values, if not specified. pad those values by 2.5% above
# and below for appearance.
def track_min_max(primary_data, secondary_data, nice_hlines, hide_secondary = False, pad_prop=0.025):
    dv = []
    iterlist = list(primary_data.values())
    if not hide_secondary:
        iterlist+=list(secondary_data.values())

    for ivallist in iterlist:
        for x in ivallist:
            if isinstance(x[-1], tuple):
                dv.append(x[-1][0])  # assume score is in first data column
            else:
                dv.append(x[-1])

    if dv:
        min_dv, max_dv = min(dv), max(dv)

    else:
        return 0, 0

    if not nice_hlines and max_dv > 10:
        spread = max_dv - min_dv
        pad = spread*pad_prop
        return max(0,min_dv - pad), max_dv + pad

    else:
        min_dv = 0
        om = np.floor(np.log10(max_dv))
        cap = 10.0**om
        newmax =  np.ceil(max_dv/cap)*cap
        if newmax - max_dv > cap/2:
            newmax-=cap/2

        return min_dv, newmax


def create_kwargs(kwtype="Collection", facecolors=None, edgecolors=None, marker='.', markersize=1., linewidth=1.,
                  alpha=1., fontsize=7.):
    if kwtype == "Scatter":
        curr_kwargs = {
            'edgecolors': edgecolors,
            'facecolors': facecolors,
            'marker': marker,
            's': markersize,
            'linewidth': linewidth
        }

    elif kwtype == "Line2D":
        if facecolors is None:
            facecolors = 'none'
        if edgecolors is None:
            edgecolors = 'none'
        curr_kwargs = {
            'markeredgecolor': edgecolors,
            'markerfacecolor': facecolors,
            'marker': marker,
            'markersize': markersize,
            'linewidth': linewidth
        }

    elif kwtype == "Patch":
        if facecolors is None:
            facecolors = 'none'
        if edgecolors is None:
            edgecolors = 'none'
        curr_kwargs = {
            "edgecolor": edgecolors,
            "facecolor": facecolors,
            "linewidth": linewidth
        }

    elif kwtype == "Text":
        curr_kwargs = {
            "fontsize": fontsize,
            "facecolor": facecolors
        }

    else:
        curr_kwargs = {}

    return curr_kwargs


# TODO: refactor to just use a dictionary (like parse_feature_yaml)
def parse_main_args_yaml(args):
    with open(args.input_yaml_file) as f:
        sample_data = yaml.safe_load(f)
        if "cycles_file" in sample_data:
            args.cycles_file = sample_data.get("cycles_file")
            print(args.cycles_file)
            args.cycle = str(sample_data.get("cycle"))
        else:
            args.structure_bed = sample_data.get("structure_bed")

        if "om_alignments" in sample_data:
            args.om_alignments = sample_data.get("om_alignments")
        if "contigs" in sample_data:
            args.contigs = sample_data.get("contigs")
        if "om_segs" in sample_data:
            args.om_segs = sample_data.get("om_segs")
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
        if "gene_highlight_list" in sample_data:
            args.gene_highlight_list = sample_data.get("gene_highlight_list")
        if "gene_fontsize" in sample_data:
            args.gene_fontsize = sample_data.get("gene_fontsize")
        if "tick_fontsize" in sample_data:
            args.tick_fontsize = sample_data.get("tick_fontsize")
        if "tick_type" in sample_data:
            args.tick_type = sample_data.get("tick_type")
        if "center_hole" in sample_data:
            args.center_hole = sample_data.get("center_hole")
        if "interior_segments_cycle" in sample_data:
            args.interior_segments_cycle = sample_data.get("interior_segments_cycle")
        if "hide_chrom_color_legend" in sample_data:
            args.hide_chrom_color_legend = sample_data["hide_chrom_color_legend"]
        if "structure_coloring" in sample_data:
            args.structure_coloring = sample_data["structure_coloring"]
        if "annotate_structure" in sample_data:
            args.annotate_structure = sample_data["annotate_structure"]


def parse_feature_yaml(yaml_file, index, totfiles):
    with open(yaml_file) as yf:
        # specifies the default track properties
        dd = {
            'tracktype': "standard",
            'primary_feature_bedgraph': "",
            'secondary_feature_bedgraph': "",
            'primary_style': 'points',
            'rescale_by_secondary': False,
            'rescale_secondary_to_primary': False,
            'rescale_by_count': False,
            'secondary_style': 'points',
            'hide_secondary': False,
            'indicate_zero': None,
            'num_hlines': 5,
            'nice_hlines': True,
            'grid_legend_fontsize': 4,
            'granularity': 0,
            'end_trim': 50,
            'show_segment_copy_count': False,
            'segment_copy_count_scaling': 1,
            #'background_color': 'auto',
            'primary_smoothing': 0,
            'secondary_smoothing': 0,
            'primary_upper_cap': None,
            'primary_lower_cap': None,
            'secondary_upper_cap': None,
            'secondary_lower_cap': None,
            'sec_resc_zero': 0, #this is a 'private' param. User setting it won't change anything.
            'linkpoint': "", #or 'midpoint'
            'link_single_match': False,
            'primary_kwargs': {},
            'secondary_kwargs': {},
            'link_kwargs': {},
            'hline_kwargs': {},
            'background_kwargs': {}
        }

        indd = yaml.safe_load(yf)
        print(indd)
        dd.update(indd)

        lkeys = ['tracktype', 'primary_style', 'secondary_style', 'linkpoint']
        for lkey in lkeys:
            dd[lkey] = dd[lkey].lower()

        # set kwargs
        # primary
        if dd['primary_style'] == "points":
            ktype = 'Scatter'
        else:
            ktype = 'Line2D'
        pkw = create_kwargs(kwtype=ktype, facecolors='k', markersize=1.0 / (2*totfiles), linewidth=2.0 / (totfiles))
        pkw.update(dd["primary_kwargs"])
        dd['primary_kwargs'] = pkw
        # secondary
        if dd['secondary_style'] == "points":
            ktype = 'Scatter'
        else:
            ktype = 'Line2D'
        skw = create_kwargs(kwtype=ktype, facecolors='lightgreen', markersize=1.0 / (2*totfiles), linewidth=2.0 / (totfiles))
        skw.update(dd["secondary_kwargs"])
        dd['secondary_kwargs'] = skw
        # feature 'links'
        linkw = create_kwargs(kwtype="Patch", facecolors='r', alpha=0.5)
        linkw.update(dd["link_kwargs"])
        dd['link_kwargs'] = linkw
        # legend lines
        llkw = create_kwargs(kwtype="Line2D", facecolors="auto", linewidth=0.25) #auto - should be lightgrey
        llkw.update(dd['hline_kwargs'])
        dd['hline_kwargs'] = llkw
        # background color for track
        bgkw = create_kwargs(kwtype="Patch", facecolors="auto", linewidth=0)
        bgkw.update(dd['background_kwargs'])
        dd['background_kwargs'] = bgkw

        # preprocess data.
        # TODO: refactor to function
        primary_data = defaultdict(list)
        secondary_data = defaultdict(list)

        if dd['tracktype'] == 'standard':
            if dd["primary_feature_bedgraph"]:
                primary_data = parse_bed(dd['primary_feature_bedgraph'])

            if dd["secondary_feature_bedgraph"]:
                secondary_data = parse_bed(dd['secondary_feature_bedgraph'])
            else:
                dd['rescale_secondary_to_primary'] = False

            minprimary = min([x[2] for y in primary_data for x in primary_data[y]])
            maxprimary = max([x[2] for y in primary_data for x in primary_data[y]])
            print("MAX PRIM, MINPRIM", maxprimary, minprimary)

            if len(secondary_data) > 0:
                minsecondary = min([x[2] for y in secondary_data for x in secondary_data[y]])
                maxsecondary = max([x[2] for y in secondary_data for x in secondary_data[y]])
            else:
                minsecondary, maxsecondary = 0, 0

            if dd['secondary_upper_cap']:
                # maxprimary = min(maxprimary, dd['upper_cap'])
                maxsecondary = min(maxsecondary, dd['secondary_upper_cap'])

            if dd['secondary_lower_cap']:
                # minprimary = max(minprimary, dd['lower_cap'])
                minsecondary= max(minsecondary, dd['secondary_lower_cap'])

            if dd['primary_upper_cap']:
                maxprimary = min(maxprimary, dd['primary_upper_cap'])
                # maxsecondary = min(maxsecondary, dd['upper_cap'])

            if dd['primary_lower_cap']:
                minprimary = max(minprimary, dd['primary_lower_cap'])
                # minsecondary= max(minsecondary, dd['lower_cap'])

            if dd['rescale_secondary_to_primary']:
                sec_rsf = (maxprimary - minprimary) / (maxsecondary - minsecondary)
                print("RESCALNG secondary TO primary")
                # print((maxsecondary - minsecondary),(maxprimary - minprimary))

                rs_sec = defaultdict(list)
                for chrom, ivall in secondary_data.items():
                    for ival in ivall:
                        cdat = ival[2]
                        if dd['secondary_upper_cap']:
                            cdat = min(cdat, dd['secondary_upper_cap'])
                        if dd['secondary_lower_cap']:
                            cdat = max(cdat, dd['secondary_lower_cap'])
                        rsdat = (cdat - minsecondary) * sec_rsf + minprimary
                        rs_sec[chrom].append([ival[0], ival[1], rsdat])

                secondary_data = rs_sec
                dd['secondary_upper_cap'], dd['secondary_lower_cap'] = None, None
                dd['sec_resc_zero'] = (-1.0*minsecondary) * sec_rsf + minprimary

                # print((-1.0*minsecondary) / (maxsecondary - minsecondary), dd['sec_resc_zero'],'sec_resc_zero')

            dv_min, dv_max = track_min_max(primary_data, secondary_data, dd['nice_hlines'],
                                           hide_secondary=dd['hide_secondary'], pad_prop=0.025)

        elif dd['tracktype'] == 'links' or dd['tracktype'] == 'link':
            dd['tracktype'] = 'links'

            if dd["primary_feature_bedgraph"]:
                primary_data = parse_bed(dd['primary_feature_bedgraph'])

            if dd["secondary_feature_bedgraph"]:
                secondary_data = parse_bed(dd['secondary_feature_bedgraph'])

            dv_min, dv_max = track_min_max(primary_data, secondary_data, nice_hlines=False,
                                           hide_secondary=dd['hide_secondary'], pad_prop=0)

        elif dd['tracktype'] == 'rects':
            if dd["primary_feature_bedgraph"]:
                primary_data = parse_bed(dd['primary_feature_bedgraph'], store_all_additional_fields=True)

            if dd["secondary_feature_bedgraph"]:
                secondary_data = parse_bed(dd['secondary_feature_bedgraph'], store_all_additional_fields=True)

            dv_min, dv_max = 0, 1

        else:
            print("ERROR: feature " + str(index) + ": Unrecognized track type - " + str(dd['tracktype']) + "\n")
            sys.exit(1)

    new_cfc = feature_track(index, primary_data, secondary_data, dd, dv_min, dv_max)
    if dd['rescale_secondary_to_primary']:
        new_cfc.minsec = minsecondary
        new_cfc.sec_rsf = sec_rsf
        new_cfc.sec_rss = minprimary

    return new_cfc