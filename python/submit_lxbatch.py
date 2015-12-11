#!/bin/env python

# given a filelist, submit jobs for each file
#
# davide.gerbaudo@gmail.com
# Oct 2015

# for each file you need to run something like
# TruthPlot -v -f  /tmp/gerbaudo/dummy_list.txt -s samplename -o out/foo -n 10 2>&1
# the output will be in out/foo/hist-samplename.root

import os, sys, time
from optparse import OptionParser

def commonPrefix(list) : return os.path.commonprefix(list)
def commonSuffix(list) : return os.path.commonprefix([l[::-1] for l in list])[::-1]
def determine_prefix_suffix(input_files=[]):
    "given a list of files, determine the longest common prefix/suffix, but stop at separators"
    prefix = commonPrefix(input_files)
    suffix = commonSuffix(input_files)
    separators = [',','.','-','_','/']
    def ltruncate(w, separators=separators):
        start = min(w.find(s) for s in separators if w.find(s)!=-1)
        return w[start:]
    def rtruncate(w, separators=separators):
        stop = max(w.rfind(s) for s in separators if w.rfind(s)>-1)
        return w[:stop+1]
    prefix = rtruncate(prefix)
    suffix = ltruncate(suffix)
    return prefix, suffix

parser = OptionParser()
parser.add_option('-i', '--input-file', help='filelist.txt')
parser.add_option('-o', '--output-dir')
parser.add_option('-q', '--queue', default='8nm', help='see bqueues')
parser.add_option('-l', '--label', help='only used in jobname')
parser.add_option("-s", "--sample", help="sample name")
parser.add_option("-S", "--submit", action="store_true", default=False, help="Actually submit the jobs")

(options, args) = parser.parse_args()

if not options.input_file: parser.error('input file not given')
if not options.output_dir: parser.error('output dir not given')

thisScript = os.path.realpath(__file__)
packageName = 'print_postfit_yields' # will be postfit_yields
packageDirectory = thisScript[:thisScript.find(packageName+'/python/')]+packageName
batchScript=packageDirectory+'/script/run_lxbatch.sh'

output_dir = options.output_dir
if not os.path.isdir(output_dir):
    print 'making dir ',output_dir
    os.makedirs(output_dir)

input_files = [l.strip() for l in open(options.input_file).readlines() if l.strip() and '.root' in l]
prefix, suffix = determine_prefix_suffix(input_files)
jobIds = [f.replace(prefix, '').replace(suffix, '') for f in input_files]
output_dir = options.output_dir
queue = options.queue
sample = options.sample
submit = options.submit
for jobId, input_file in zip(jobIds, input_files):
    jobname = (sample+'_'+jobId) if sample else jobId
    jobname += (options.label if options.label else '')
    output_file = output_dir.rstrip('/ ')+'/'+jobId+'/'+'fit.root'
    bsubCmd = ("bsub "
               +" -q %s"%queue
               +" -J %s"%jobname
               +" %s"%batchScript
               +" %s"%input_file
               +" %s"%output_file
               +" %s"%jobId)
               
    print bsubCmd
    if submit :
        os.system(bsubCmd)
        time.sleep(1)

#

