import sys, random
import os
import bisect
import argparse
import numpy as np
import hg19util as hg19
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.path as mpath
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib.collections import PatchCollection
from intervaltree import Interval, IntervalTree


#Create colormap for human chromosomes
chromosome_colors = dict([(hg19.chrName[i],c) for (i,c) in enumerate(plt.cm.get_cmap(None, 23).colors)])
(gene_list,gene_map) = load_gene_list()
contig_spacing = 0.01
seg_spacing = 0.01
bar_width = 2.0/3
global_rot = 90.0
tick_positions = 10000

outer_bar = 15
bar_width = 2.0/3
bar_spacing = 1.5
bed_track_height = 3
bed_track_spacing = 0.2

class BuildBedTrack:
  
  def __init__(self, track, bed_data, ymax = None, ymin = None, is_log = False, num_axis = 10, point_spacing = 100, color = 'k'):
    self.track = track
    self.track_rmin = outer_bar-bar_width-((track+1)*bed_track_height)+bed_track_spacing-1
    self.track_rmax = outer_bar-bar_width-((track)*bed_track_height)-bed_track_spacing-1
    self.is_log = is_log
    self.ymax = max([s.info['value'] for s in bed_data]) if ymax is None else ymax
    self.ymin = min([s.info['value'] for s in bed_data]) if ymin is None else ymin
    self.num_axis = num_axis
    self.point_spacing = point_spacing
    self.bed_patches = []
    self.color = color
    self.bed_data = bed_data
  
  def build_axis(self):
    i = self.num_axis    
    while (i >= 0):
      y_value = self.ymin+(self.ymax-self.ymin)/(self.num_axis)*i
      r_value = self.track_rmin+(self.track_rmax-self.track_rmin)/(self.num_axis)*i
      self.bed_patches.append(mpatches.Wedge((0, 0), r_value, 0, 360, width=0.01))
      x_t,y_t = pol2cart(r_value,np.pi*2*1/4)
      ax.text(x_t,y_t,int(y_value),color='grey', ha='right',fontsize=3,rotation_mode='anchor')
      i+=-1  
  
  def build_segments(self):
    points_x = []
    points_y = []
    previous_end = total_length_with_spacing*(global_rot/360.0)    
    for ind,sp in enumerate(start_points):
      start_point = int(previous_end - sp)
      start_angle = start_point/total_length_with_spacing*360
      end_angle = (start_point - lens[ind])/total_length_with_spacing*360
      
      #segseqD referenced as global variable here because I'm lazy
      segment = segSeqD[cycle[ind][0]] 
      strand = cycle[ind][1]
      hits = [h[0] for h in self.bed_data.intersection([segment])]
      for h in hits:
        for pos in xrange(h.start, h.end, self.point_spacing):
          if strand == "+":
            normStart = start_point - max(0,pos-segment.start)
            normEnd = start_point - min(segment.end-segment.start,pos-segment.start)
          else:
            normEnd = start_point - min(segment.end-segment.start,segment.end-pos)
            normStart = start_point - max(0,segment.end - pos)        
          r_scale_value = 1.*(h.info['value']-self.ymin)/(self.ymax-self.ymin)*(self.track_rmax-self.track_rmin)+self.track_rmin
          
          x_s,y_s = pol2cart(r_scale_value,normStart/total_length_with_spacing*2*np.pi)
          foo = points_x.append(x_s)
          foo = points_y.append(y_s)
    foo = ax.plot(points_x,points_y,'ro',color=self.color,markersize=0.1)    
  
  @staticmethod
  def load_bed(bed_file, log = False):
    bed_data = hg19.interval_list()
    for line in open(bed_file):
      res = line.split('\t')
      bed_data.append(hg19.interval(res[0], int(res[1]), int(res[2]), info={'value':float(res[3]) if not log else 10**float(res[3])}))
    bed_data.sort()
    return bed_data

def load_gene_list():
  #f = '/pedigree2/projects/namphuon/data/references/hg19/annotations/refseq.july_2017.bed'
  f = '/pedigree2/projects/namphuon/data/references/hg19/annotations/ensembl.txt'
  input = open(f,'r')
  lines = [line.strip().split('\t') for line in input]
  lines.pop(0)
  added = {}
  gene_list = hg19.interval_list()
  foo = [(gene_list.append(hg19.interval(res[2],int(res[4]),int(res[5]),-1 if res[3] == '-' else 1, info={'data':res})),added.setdefault("%s:%s-%s;%s" % (res[2],res[4],res[5],res[3]),res[12])) for res in lines if "%s:%s-%s;%s" % (res[2],res[4],res[5],res[3]) not in added]
  input.close()
  gene_list.sort()        
  
  input = open('/pedigree2/projects/namphuon/data/references/grch38/gene.annotation.gff','r')  
  gene_map = {}  
  for line in input:
    if line[0] == '#':
      continue
    res = line.split('\t')
    res_dict = dict([(a.strip().replace('"','').split(' ')[0],a.strip().replace('"','').split(' ')[1]) for a in res[-1].split(';') if len(a.strip().split(' ')) == 2])
    if 'gene_name' in res_dict:
      gene_map[res_dict['gene_id'].split('.')[0]]=res_dict['gene_name']
  return (gene_list,gene_map)

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

def parse_genes(chrom):
  print("Building interval tree for chrom " + chrom)
  t = IntervalTree()
  with open("/pedigree2/projects/namphuon/programs/docker/aa/data_repo//hg19/human_hg19_september_2011/Genes_July_2010_hg19.gff") as infile:
    for line in infile:
      fields = line.rsplit("\t")
      #gene = fields[-2]
      currChrom = fields[2]
      tstart = int(fields[4])
      tend = int(fields[5])
      if chrom == currChrom:
        t[tstart:tend] = fields
  return t

def rel_genes(genes):  
  relGenes = {}  
  for (gene, seg) in genes:
    if gene.info['data'][-4] not in gene_map:
      continue
    #print i.data
    tstart = gene.start
    tend = gene.end
    if not (gene.info['data'][-4].startswith("LOC") or gene.info['data'][-4].startswith("LINC")):
      if gene not in relGenes:
        relGenes[gene] = [tstart, tend, [(int(z[0]),int(z[1])) for z in zip(gene.info['data'][9].split(','), gene.info['data'][10].split(',')) if z[0] != '']]
      else:
        oldTStart = relGenes[gene][0]
        oldTEnd = relGenes[gene][1]
        if tstart < oldTStart:
          oldTStart = tstart
        if tend > oldTEnd:
          oldTEnd = tend        
        relGenes[gene] = (oldTStart,oldTEnd,relGenes[2]+[(int(z[0]),int(z[1])) for z in zip(gene.info['data'][9].split(','), gene.info['data'][10].split(',')) if z[0] != ''])
  return relGenes

#need to add functionality for plotting exon posns
def plot_gene_track(currStart, relGenes, segment, total_length_with_spacing, strand):
  for ind,i in enumerate(relGenes):    
    truncStart = False
    truncEnd = False
    #e_posns is a list of tuples of exon (start,end)
    #these can be plotted similarly to how the coding region is marked
    tstart,tend,e_posns = relGenes[i]
    seg_len = segment.end - segment.start
    if strand == "+":
      normStart = currStart - max(0,tstart-segment.start)
      normEnd = currStart - min(segment.end-segment.start,tend-segment.start)
    else:
      normEnd = currStart - min(segment.end-segment.start,segment.end-tstart)
      normStart = currStart - max(0,segment.end - tend)
    # if tstart < pTup[1]:
    #   truncStart = True
    # if tend > pTup[2]:
    #   truncEnd = True
    start_angle = normStart/total_length_with_spacing*360
    end_angle = normEnd/total_length_with_spacing*360
    
    text_angle = (start_angle + end_angle)/2.0
    if end_angle < 0 and start_angle > 0:
      end_angle+=360
    
    patches.append(mpatches.Wedge((0,0), outer_bar, end_angle, start_angle, width=bar_width/2.0))
    f_color_v.append('k')
    e_color_v.append('k')
    lw_v.append(0)
    
    # x,y = pol2cart(outer_bar+(bar_width/2.0),(text_angle/360*2*np.pi))
    x_t,y_t = pol2cart(outer_bar + bar_width + 0.5,(text_angle/360*2*np.pi))
    #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.4)
    
    if text_angle < -90 and text_angle > -360:
      text_angle+=180
      ax.text(x_t,y_t,gene_map[i.info['data'][-4]],color='grey',rotation=text_angle,
        ha='right',fontsize=4,rotation_mode='anchor')    
    else:
       ax.text(x_t,y_t,gene_map[i.info['data'][-4]],color='grey',rotation=text_angle,
        ha='left',fontsize=4,rotation_mode='anchor')

#g_locs,cycle[cycle_num],seg_padding
def plot_ref_genome(start_points,lens,cycle,total_length_with_spacing,seg_padding,bed_feat_list=[]):
  seg_posns = []
  previous_end = total_length_with_spacing*(global_rot/360.0)
  for ind,sp in enumerate(start_points):
    start_point = int(previous_end - sp)
    start_angle = start_point/total_length_with_spacing*360
    end_angle = (start_point - lens[ind])/total_length_with_spacing*360
    
    #segseqD referenced as global variable here because I'm lazy
    seg_coord_tup = segSeqD[cycle[ind][0]]
    text_angle = (start_angle + end_angle)/2.0
    # text_angle = start_angle
    x,y = pol2cart(outer_bar + bar_width*8,text_angle/360*2*np.pi)
    
    #Writes the really ugly label for the segment on the outside
    #segment_name = "Segment " + cycle[ind][0] + cycle[ind][1] + "\n" + seg_coord_tup.chrom + ":" + str(seg_coord_tup.start) + "-" + str(seg_coord_tup.end)
    segment_name = "Segment " + cycle[ind][0] + cycle[ind][1] + "\n" + seg_coord_tup.chrom
    if text_angle < -90 and text_angle > -360:
      text_angle-=180
      ha = "right"
    else:
      ha = "left"
    ax.text(x,y,segment_name,color='k',rotation = text_angle,ha=ha,fontsize=5,rotation_mode='anchor')
    
    #makes the reference genome wedges  
    if end_angle < 0 and start_angle > 0:
      end_angle+=360
      
    patches.append(mpatches.Wedge((0,0), outer_bar, end_angle, start_angle, width=bar_width))
    f_color_v.append(chromosome_colors[seg_coord_tup.chrom])
    e_color_v.append('grey')
    lw_v.append(0.2)
    
    #makes the ticks on the reference genome wedges
    if cycle[ind][1] == "+":
      posns = zip(range(seg_coord_tup.start,seg_coord_tup.end+1),np.arange(start_point,start_point-lens[ind]-1,-1))
    else:
      posns = zip(np.arange(seg_coord_tup.end,seg_coord_tup.start-1,-1),np.arange(start_point,start_point-lens[ind]-1,-1))
    
    for j in posns:
      if j[0] % tick_positions == 0:
        text_angle = j[1]/total_length_with_spacing*360
        x,y = pol2cart(outer_bar,(text_angle/360*2*np.pi))
        x_t,y_t = pol2cart(outer_bar + 0.2,(text_angle/360*2*np.pi))
        ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.2)
        
        if text_angle < -90 and text_angle > -360:
          text_angle-=180
          ha = "right"
          txt = str(int(round((j[0])/tick_positions))) + " "
        else:
          ha = "left"
          txt = " " + str(int(round((j[0])/tick_positions)))
        ax.text(x_t,y_t,txt,color='grey',rotation=text_angle,
        ha=ha,fontsize=2.5,rotation_mode='anchor')    
    #overhaul
    genes = gene_list.intersection([seg_coord_tup])
    relGenes = rel_genes(genes)
    plot_gene_track(start_point,relGenes,seg_coord_tup,total_length_with_spacing,cycle[ind][1])
    
#OVERHAUL
def plot_ref_cmaps(start_point,seg_cmap_vector):
  start_angle = start_point/total_length_with_spacing*360
  end_angle = (start_point - seg_cmap_vector[-1])/total_length_with_spacing*360
  if end_angle < 0 and start_angle >0:
    end_angle+=360  
  patches.append(mpatches.Wedge((0,0), mid_bar+bar_width, end_angle, start_angle, width=bar_width))
  f_color_v.append('darkorange')
  e_color_v.append('k')
  lw_v.append(0)
  
  lab_locs = []
  print seg_cmap_vector
  for i in seg_cmap_vector[:-1]:
    label_angle = (start_point - i)/total_length_with_spacing*2*np.pi
    x,y = pol2cart(mid_bar,label_angle)
    x_t,y_t = pol2cart(mid_bar+bar_width,label_angle)
    ax.plot([x,x_t],[y,y_t],color='k',alpha=0.9,linewidth=0.2)
    lab_locs.append(start_point - i)
  return lab_locs

#parse the alignment file
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
      aln_vect.append(fields_dict)
  return aln_vect,meta_dict

#parse cycles file
def parse_cycles_file(cycles_file):
  cycles = {}
  segSeqD = {}
  with open(cycles_file) as infile:
    for line in infile:
      if line.startswith("Segment"):
        fields = line.rstrip().split()
        lowerBound = int(fields[3])
        upperBound = int(fields[4])
        chrom = fields[2]
        segNum = fields[1]
        segSeqD[segNum] = hg19.interval(chrom,lowerBound,upperBound)
      elif "Cycle=" in line:
        curr_cycle = []
        fields = line.rstrip().rsplit(";")
        lineD = {x.rsplit("=")[0]:x.rsplit("=")[1] for x in fields}
        segs = lineD["Segments"].rsplit(",")
        #TODO: Need to rotate the segs in case there's a 0 within the path, or else path is incorrect
        for i in segs:
          seg = i[:-1]
          if seg != "0":
            strand = i[-1]
            curr_cycle.append((seg,strand))
        cycles[lineD["Cycle"]] = curr_cycle
  return cycles,segSeqD

#parse bed file
def parse_bed_file(bed_file):
  bed_list = []
  print "f",bed_file
  with open(bed_file) as infile:
    for line in infile:
      fields = line.rstrip().split()
      fields[1] = int(fields[1])
      fields[2] = int(fields[2])
      fields[3] = float(fields[3])
      bed_list.append(fields)
  return bed_list

def get_contig_locs(aln_vect,meta_dict):
  #go through the seg_seq and calculate total lengths of alignments
  seg_seq = meta_dict["seg_seq"].split(",")
  total_seg_length = [seg_cmap_lens[i] for i in seg_seq]
  contig_list = []
  prev = None
  for a_d in aln_vect:
    if a_d["contig_id"] != prev:
      prev = a_d["contig_id"]
      contig_list.append[prev]
  total_contig_length = sum([contig_cmap_lens[x] for i in contig_set])
  total_genomic_length = max(total_seg_length,total_contig_length)
  spacing = contig_spacing*total_contig_length
  total_genomic_length+=(spacing*len(contig_set))
  
  start_points = [0]
  for ind in len(contig_list)-1:
    start_points.append(start_points[-1] + spacing + contig_cmap_lens[contig_list[ind]]) #this gets the previous ind!
  return dict(zip(contig_list,start_points))

def get_seg_locs_from_contig_aln():
  pass

#get start locations for a cycle
def get_seg_locs_from_cycle(cycle,segSeqD):
  lens = []
  for i in cycle:
    curr_len = segSeqD[i[0]].end - segSeqD[i[0]].start
    lens.append(curr_len)
  total_seg_len = sum(lens)
  seg_padding = seg_spacing*total_seg_len
  start_points = [0]
  for i in lens[:-1]:
    start_points.append(start_points[-1] + i + seg_padding)
  total_length = start_points[-1] + lens[-1] + seg_padding
  return start_points,lens,float(total_length),seg_padding

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
parser.add_argument("--bionano_alignments",help="turns on Bionano visualization mode (requires contigs,segs,key,path_alignment args")
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("--cycles_file",help="cycles file")
parser.add_argument("--cycle",help="cycle number to visualize")
parser.add_argument("-k", "--keyfile", help="segments cmap key file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--bed_files",help="bed file list",nargs="+")
parser.add_argument("--feature_labels",help="bed feature names",nargs="+")
#args = parser.parse_args()
args.cycles_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/onco_amplicon1_cycles.txt'
args.cycle = '16'

if not args.sname:
  args.sname = args.segs.rsplit(".")[0] + "_"
else:
  samp_name = args.sname.rsplit("/")[-1]

fname = samp_name

bed_feat_dict = {}
if args.bed_files:
  for i,j in zip(args.bed_files,args.feature_labels):
    print j,i
    #feature name -> chromosome -> ordered list of positions
    bed_list = parse_bed_file(i)
    bed_feat_dict[j] = feat_bed_to_lookup(bed_list)

outer_bar = max(bed_track_height*(len(bed_feat_dict)+2),10)

if args.bionano_alignments:
  #unfinished utility for bionano
  #parse cmaps
  segment_locs = get_segment_locs_from_contig_aln(args.path_alignment)
  seg_cmaps = parse_cmap(args.segs,True)
  contig_cmaps = parse_cmap(args.contigs,True)
  #get cmap lens
  seg_cmap_lens = get_cmap_lens(args.segs)
  contig_cmap_lens = get_cmap_lens(args.contigs) 

#handles basic (non-bionano case)
else:
  cycles,segSeqD = parse_cycles_file(args.cycles_file)
  print cycles
  cycle_num = args.cycle
  print cycles[cycle_num]
  start_points,lens,total_length_with_spacing,seg_padding = get_seg_locs_from_cycle(cycles[cycle_num],segSeqD)
  plot_ref_genome(start_points, lens,cycles[cycle_num],total_length_with_spacing,seg_padding,args.feature_labels)


bed_data = hg19.interval_list([hg19.interval('chr8', 127638302, 127938302, info={'value':int(random.random()*100)}), hg19.interval('chr8', 128716346,128746346, info={'value':int(random.random()*100)})])
bed_data.sort()


plt.clf()
fig, ax = plt.subplots()

patches = []
f_color_v = []
e_color_v = []
lw_v = []
cycle = cycles[cycle_num]
plot_ref_genome(start_points, lens, cycle,total_length_with_spacing,seg_padding,args.feature_labels)

p = PatchCollection(patches)
p.set_facecolor(f_color_v)
p.set_edgecolor(e_color_v)
p.set_linewidth(lw_v)
ax.add_collection(p)
ax.set_xlim(-(outer_bar+5), (outer_bar+5))
ax.set_ylim(-(outer_bar+5), (outer_bar+5))
ax.set_aspect(1.0)

wgs_bed = BuildBedTrack(0, BuildBedTrack.load_bed('/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/colo320dm.wgs.1000.pileuplog.bed', log = True))
wgs_bed.build_axis()
axis = PatchCollection(wgs_bed.bed_patches)
ax.add_collection(axis)
wgs_bed.build_segments()

rna_bed = BuildBedTrack(1, BuildBedTrack.load_bed('/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/colo320dm_10.fpkm.bed', log = True))
rna_bed.build_axis()
axis = PatchCollection(rna_bed.bed_patches)
ax.add_collection(axis)
rna_bed.build_segments()

plt.axis('off')
plt.savefig(fname + '.png',dpi=600)
plt.close()
print "done"