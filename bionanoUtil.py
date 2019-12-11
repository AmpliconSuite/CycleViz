"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

#compute the median
def median(L):
    if L:
        L = sorted(L)
        n = len(L)
        m = n - 1
        return (L[n/2] + L[m/2]) / 2.0

    return None

#parse cmap into a dictionary, to maintain the 1-indexing used in this format. 
#specify keep_length to keep the length field of the cmap entry.
def parse_cmap(cmapf,keep_length = False):
    cmaps = {}
    #contigCovs = {}
    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                if fD["CMapId"] not in cmaps:
                    cmaps[fD["CMapId"]] = {}
                    #contigCovs[fD["CMapId"]] = {}

                #this is not a good way to parse label channel means color channel
                if fD["LabelChannel"] == "1":
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])
                    #contigCovs[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Coverage"])
                    
                elif fD["LabelChannel"] == "0" and keep_length:
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])

    return cmaps

def get_cmap_lens(cmapf):
    cmap_lens = {}
    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                if fD["CMapId"] not in cmap_lens:
                    cmap_lens[fD["CMapId"]] = float(fD["ContigLength"])

    return cmap_lens


#parse a bnx file into a vector. Can specify keep_length to keep the length value.
def parse_bnx(bnxF,keep_length = False):
    moleculeD = {}
    with open(bnxF) as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.rstrip().rsplit("\t")
                if line.startswith('0'):
                    currKey = fields[1]
                elif line.startswith('1'):
                    if keep_length:
                        moleculeD[currKey] = [float(x) for x in fields[1:]]
                    else:
                        moleculeD[currKey] = [float(x) for x in fields[1:-1]] 

    return moleculeD

#parse a bnx file and store map of mol_id -> mol_length
def get_mol_lens(bnxF):
    moleculeLenD = {}
    with open(bnxF) as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.rstrip().rsplit("\t")
                if line.startswith('0'):
                    currKey = fields[1]
                elif line.startswith('1'):
                    #gets the length of the molecule, note the position of the last label
                    moleculeLenD[currKey] = float(fields[-1])

    return moleculeLenD

#read in key file 
def parse_keyfile(keyF_name):
    keyCompD = {}
    with open(keyF_name) as infile:
        for line in infile:
            if not line.startswith("#"):
                if line.startswith("CompntId"):
                    head = line.rstrip().split()
                else:
                    fields = line.rstrip().split()
                    keyCompD[fields[0]] = (fields[1],float(fields[2]))
                    
    return keyCompD


#parse xmap
def parse_xmap(xmapf):
    detailFields = ["XmapEntryID","QryContigID","RefContigID","Orientation","Confidence","QryLen","RefLen",
    "QryStartPos","QryEndPos","RefStartPos","RefEndPos","Alignment"]
    
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                alnstring = ")" + fD["Alignment"] + "("
                #xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x:fD[x] for x in detailFields}

    return xmapPair

#Swap reference and query for a given xmap
def swap_xmap_RQ(xmapD):
    for xmap_id,fD in xmapD.iteritems():
        # fD["QryLen"],fD["RefLen"] = fD["RefLen"],fD["QryLen"] #do this in the parsing itself
        fD["QryStartPos"],fD["RefStartPos"] = fD["RefStartPos"],fD["QryStartPos"]
        fD["QryEndPos"],fD["RefEndPos"] = fD["RefEndPos"],fD["QryEndPos"]
        fD["QryContigID"],fD["RefContigID"] = fD["RefContigID"],fD["QryContigID"]
        aln_pairs = fD["Alignment"] 
        fD["Alignment"] = [(y,x) for x,y in aln_pairs]


#can handle poorly formatted .xmap files, such as those from OMBlast. Requires more inputs than parse_xmap.
def parse_generic_xmap(xmapf,qryLenD,refLenD,swap_Ref_Qry = False):
    detailFields = ["XmapEntryID","QryContigID","RefContigID","Orientation","Confidence","QryLen","RefLen",
    "QryStartPos","QryEndPos","RefStartPos","RefEndPos","HitEnum"]
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                
                #handle mis-capitalizations
                for i in fD.keys():
                    for dname in detailFields:
                        if i.lower() == dname.lower():
                            fD[dname] = fD[i]

                fD["Confidence"] = float(fD["Confidence"])

                #refactor to eliminate,reduce, or simplify need for extra checks
                try:
                    fD["QryLen"],fD["RefLen"] = float(fD["QryLen"]),float(fD["RefLen"])
                except KeyError:
                    if not swap_Ref_Qry:
                        fD["QryLen"],fD["RefLen"] = qryLenD[fD["QryContigID"]], refLenD[fD["RefContigID"]]
                    else:
                        fD["QryLen"],fD["RefLen"] = qryLenD[fD["RefContigID"]], refLenD[fD["QryContigID"]]


                fD["QryStartPos"],fD["QryEndPos"] = sorted([float(fD["QryStartPos"]),float(fD["QryEndPos"])])
                if fD["Orientation"] == "+":
                    fD["RefStartPos"],fD["RefEndPos"] = sorted([float(fD["RefStartPos"]),float(fD["RefEndPos"])])
                else:
                    fD["RefEndPos"],fD["RefStartPos"] = sorted([float(fD["RefStartPos"]),float(fD["RefEndPos"])])
                
                xmapPair[fD["XmapEntryID"]] = {x:fD[x] for x in detailFields}
                       
                try:
                    alnstring = ")" + fD["Alignment"] + "("
                    aln_pairs = [(int(x.rsplit(",")[0]),int(x.rsplit(",")[1])) for x in alnstring.rsplit(")(")[1:-1]]             
                    xmapPair[fD["XmapEntryID"]]["Alignment"] = aln_pairs
                except KeyError:
                    #xmap does not have Alignment field
                     xmapPair[fD["XmapEntryID"]]["Alignment"] = []


    #handle the case where the user wants to swap the reference and qry (e.g. segments aligned to contigs)
    if swap_Ref_Qry: swap_xmap_RQ(xmapPair)

    return xmapPair

#make cmap dictionary into 0-indexed vector
def vectorize_cmaps(cmap_d):
    vectorized_dict = {}
    for y in cmap_d:
        y_posns = [cmap_d[y][k] for k in sorted(cmap_d[y].keys())]
        vectorized_dict[y] = y_posns
            
    return vectorized_dict

def parse_bed(bedfile):
    bed_list = []
    with open(bedfile) as infile:
        for line in infile:
            if not line.startswith("#"):
                fields = line.rstrip().rsplit()
                bed_list.append([fields[0],int(fields[1]),int(fields[2])])
        
    return bed_list

def dict_from_bed_list(bed_list):
    bed_dict = {}
    for fields in bed_list:
        if fields[0] not in bed_dict:
            bed_dict[fields[0]] = []        
        bed_dict[fields[0]].append((int(fields[1]),int(fields[2])))

    for chrom in bed_dict:
        bed_dict[chrom] = sorted(bed_dict[chrom])

    return bed_dict


#take label number from a reversed segment and write it back to forward
#ASSUMES 0-BASED INDEX
def translate_reversed_label(cmaps,rev_seg_id,r_lab):
    seg_id = rev_seg_id.rsplit("_")[0]
    nlabs = len(cmaps[seg_id])
    return nlabs - r_lab

#cmaps input must include the segment artificial end label to work properly
def add_full_reverse_cmaps(cmaps,key_dict):
    #make a reverse keyfile
    iter_keys = cmaps.keys()
    for i in iter_keys:
        tot_labs = len(cmaps[i])-1
        cmap_len = cmaps[i][tot_labs+1]
        new_ID = i + "_r"
        #add new entry to key dict
        seg_rep = "|".join(key_dict[i][0].rsplit("|")[::-1])
        key_dict[new_ID] = (seg_rep,cmap_len)
        #add new entry to cmaps
        cmaps[new_ID] = {}
        for j in range(1,tot_labs+1):
            cmaps[new_ID][tot_labs - j + 1] = cmap_len - cmaps[i][j]
        
        cmaps[new_ID][tot_labs+1] = cmap_len

#binary search find the label corresponding to some position in a cmap dict.
#return the RIGHT index of bisect (cmap is 1 based)
def pos_to_label(x, item_cmap):
    import bisect
    arr = [item_cmap[k] for k in range(1,max(item_cmap.keys())+1)]
    return bisect.bisect(arr,x)

#convert XMAP format to SegAligner alignment format. OMPathFinder requires alignments in SA format.
def xmap_to_SA_aln(xmapD,outdir,fname_prefix,ref_cmaps,contig_cmaps):
    seg_contig_count = {}
    for xmap_id,fD in xmapD.iteritems():
        contig_id = fD["QryContigID"]
        seg_id = fD["RefContigID"]
        score = fD["Confidence"]
        orientation = fD["Orientation"]
        #update number of times segment has been aligned to this contig in a particular direction
        cso_key = (contig_id,seg_id,orientation)
        if cso_key not in seg_contig_count:
            seg_contig_count[cso_key] = 0

        seg_contig_count[cso_key]+=1

        outname = outdir + "/" + fname_prefix + "_" + contig_id + "_" + seg_id + "_"
        if orientation == "-":
            outname+="r_"

        outname+=(str(seg_contig_count[cso_key]) + "_aln.txt")

        with open(outname,'w') as outfile:
            outfile.write("#seg_seq\ttotal_score\tcircular\n")
            outfile.write("#" + seg_id + orientation + "\t" + str(score) + "\tFalse\n")
            outfile.write("#contig_id\tseg_id\tcontig_label\tseg_label\tcontig_dir\tseg_dir\tseg_aln_number\tscore\tscore_delta\n")
            
            #handle orientation
            if fD["Alignment"]:
                alist = fD["Alignment"]

            #if no alignment string given (incomplete XMAP), make a dummy alignment
            else: 
                ref_start_label = pos_to_label(fD["RefStartPos"],ref_cmaps[seg_id])
                contig_start_label = pos_to_label(fD["QryStartPos"],contig_cmaps[contig_id])
                ref_end_label = pos_to_label(fD["RefEndPos"],ref_cmaps[seg_id])
                contig_end_label = pos_to_label(fD["QryEndPos"],contig_cmaps[contig_id])
                alist= [(ref_start_label,contig_start_label),(ref_end_label,contig_end_label)]

            if orientation == "-":
                    alist = alist[::-1]

            #write converted alignment
            for i in alist[:-1]:
                outlist = [contig_id,seg_id,str(i[1]),str(i[0]),"+",orientation,"0","0","0"]
                outfile.write("\t".join(outlist) + "\n")

            #write the last one and include total score
            i = alist[-1]
            outlist = [contig_id,seg_id,str(i[1]),str(i[0]),"+",orientation,"0",str(score),"0"]
            outfile.write("\t".join(outlist)+"\n")

            
#takes vector of cmap vector of positions, including the length of the map
def write_cmap_from_vector(cmap_vector,fname):
    header_lines = "# hostname=BioNanoUtil\n"
    header_lines += "# $ BioNanoUtil.py\n# CMAP File Version:\t0.1\n# Label Channels:\t1\n# Nickase Recognition Site 1:\tunknown\n"
    header_lines += "# Number of Consensus Maps:\t"
    header_lines += str(len(cmap_vector))
    header_lines += "\n# Values corresponding to intervals (StdDev, HapDelta) refer to the interval between current site and next site\n#h\tCMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\tChimQuality\n#f\tint\tfloat\tint\tint\tint\tfloat\tfloat\tfloat\tfloat\tfloat\n"
    with open(fname,'w') as outfile:
        outfile.write(header_lines)
        for ind,cmap_posns in enumerate(cmap_vector):
            map_b_len = str(cmap_posns[-1])
            map_l_len = str(len(cmap_posns)-1)
            for p_i,pos in enumerate(cmap_posns):
                outfile.write("\t".join([str(ind+1),map_b_len,map_l_len,str(p_i+1),"1",str(pos),"1.0","1.0","1.0","0.0"]) + "\n")

#parses the output from SegAligner 
def parse_seg_alignment_file(alignfile):
    alignment = []
    tip_aln = True if "_tip_" in alignfile else False
    with open(alignfile) as infile:
        meta_head = infile.next().rstrip()[1:].rsplit()
        meta_vals = infile.next().rstrip()[1:].rsplit()
        meta_dict = dict(zip(meta_head,meta_vals))
        aln_head = infile.next().rstrip()[1:].rsplit()
        for line in infile:
            fields = line.rstrip().rsplit()
            alignment.append(dict(zip(aln_head,fields)))

        seg_id = meta_dict["seg_seq"][:-1]
        strand = meta_dict["seg_seq"][-1]
        tot_score = float(meta_dict["total_score"])
        seg_start = int(alignment[0]["seg_label"])
        seg_end = int(alignment[-1]["seg_label"])
        seg_ends = (seg_start,seg_end)
        contig_start = int(alignment[0]["contig_label"])
        contig_end = int(alignment[-1]["contig_label"])
        contig_ends = (contig_start,contig_end)

    return alignment[0]["contig_id"],[seg_id,seg_ends,contig_ends,strand,tot_score,alignment,tip_aln]