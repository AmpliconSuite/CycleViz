#!/usr/bin/env python

import argparse
from copy import copy
import os
from subprocess import call
import sys
import threading

try:
    CV_SRC = os.environ["CV_SRC"]

except KeyError:
    print("$CV_SRC environment variable must be set. See CycleViz README.")
    sys.exit(1)


# generic worker thread function
class workerThread(threading.Thread):
    def __init__(self, threadID, target, *args):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self._target = target
        self._args = args
        threading.Thread.__init__(self)

    def run(self):
        self._target(*self._args)


def run_cycles_conversion(CV_SRC, slist, updated_slist, ccf_loc, cyclic_only=False):
    while True:
        try:
            sname, cfile, gfile = slist.pop()

        except IndexError:
            break

        opath = "{}{}_BPG_converted".format(ccf_loc, cfile.rsplit("/")[-1].rsplit("_cycles.txt")[0])
        cmd = "{}/convert_cycles_file.py -c {} -g {} -o {}".format(CV_SRC, cfile, gfile, opath)
        call(cmd, shell=True)
        newcf = opath + "_cycles.txt"

        # open the newcf, get one entry per entry
        clist = []
        with open(newcf) as infile:
            for line in infile:
                if line.startswith("Cycle="):
                    cf = line.rstrip().rsplit(";")
                    cd = {x.rsplit("=")[0]:x.rsplit("=")[1] for x in cf}
                    if not cd["Segments"].startswith("0") or not cyclic_only:
                        clist.append((sname, newcf, gfile, cd["Cycle"]))

        updated_slist.extend(clist)


def run_CycleViz(CV_SRC, clist, CV_yaml, ref, cv_image_loc, feature_yamls=None):
    while True:
        try:
            sname, cfile, gfile, cnum = clist.pop()

        except IndexError:
            break

        opath = cv_image_loc + cfile.rsplit("/")[-1].rsplit("_cycles.txt")[0]
        cmd = "{}/CycleViz.py --input_yaml_file {} --ref {} -o {} --cycles_file {} -g {} --cycle {}".format(
            CV_SRC, CV_yaml, ref, opath, cfile, gfile, cnum)

        if feature_yamls:
            cmd += " --feature_yaml_list {}".format(" ".join(feature_yamls))

        print(cmd)
        call(cmd, shell=True)


parser = argparse.ArgumentParser(description="Batch script for circular visualizations of genome structures")
parser.add_argument("--CV_yaml_file", help="Specifiy all desired arguments for batch CycleViz plots in this file\n",
                    required=True)
parser.add_argument("--input", help="Three-column, mutli-row file formatted as: sample_name cycles_file graph_file",
                    required=True)
parser.add_argument("--ref", help="reference genome", choices=["hg19", "hg38", "GRCh37", "GRCh38", "mm10", "GRCm38"],
                    required=True)
parser.add_argument("--feature_yaml_list", nargs='+', help="list of the input yamls for bedgraph file feature "
                    "specifying additional data. Will be plotted from inside to outside given the order the filenames "
                    "appear in", default=[])
parser.add_argument("-t", "--nthreads", type=int, help="Number of threads to use. Default=1", default=1)
parser.add_argument("--cyclic_only", help="Only visualize entries marked as cyclic", action='store_true')
args = parser.parse_args()

datalist = []
updated_datalist = []
ccf_loc = "converted_cycles_files/"
if not os.path.exists(ccf_loc):
    os.makedirs(ccf_loc)

cv_image_loc = "cycleviz_images/"
if not os.path.exists(cv_image_loc):
    os.makedirs(cv_image_loc)

with open(args.input) as infile:
    for l in infile:
        fields = l.rstrip().rsplit()
        if len(fields) != 3:
            print("improperly formatted line: " + l)
            continue

        sname, original_cfile, gfile = fields
        datalist.append((sname, original_cfile, gfile))

# reformat the cycles files
print("Running CycleViz on " + str(len(updated_datalist)) + " with " + str(args.nthreads) + " threads...")
threadL = []
for i in range(int(args.nthreads)):
    threadL.append(workerThread(i, run_cycles_conversion, CV_SRC, datalist, updated_datalist, ccf_loc,
                                args.cyclic_only))
    threadL[i].start()

for t in threadL:
    t.join()

print(len(updated_datalist))

#now, make the visualizations!
threadL = []
for i in range(int(args.nthreads)):
    threadL.append(workerThread(i, run_CycleViz, CV_SRC, updated_datalist, args.CV_yaml_file, args.ref, cv_image_loc,
                                args.feature_yaml_list))
    threadL[i].start()

for t in threadL:
    t.join()

