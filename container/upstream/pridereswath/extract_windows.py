import argparse
import sys
parser = argparse.ArgumentParser(description='Extract and stratify windows.')
parser.add_argument('--in', dest='file', required=True, type=str, help='input mzML file')
args = parser.parse_args()
if not args.file:
	parser.print_help()
	import sys
	sys.exit()

import subprocess
import csv
from itertools import islice, tee
from io import StringIO
import os

def sliding_window_iter(iterable, size):
    iterables = tee(iter(iterable), size)
    window = zip(*(islice(t, n, None) for n,t in enumerate(iterables)))
    yield from window

#file = "/hps/nobackup2/proteomics/SwathMSdata/pride_datasets/PXD007130/test/work/07/2381fc52f2fad4b1c34a98815533cc/test1.mzML"
# window extraction shell command
mzmlfile = args.file
wesc = "cat {mzml} | grep -A 2 MS:1000827 | head -n 400 | cut -d '=' -f8 | cut -d '\"' -f2 | paste -d, - - - - | sed 's/.\{{3\}}$//' | awk -F ',' '{{$4=$1-$2;}} {{$5=$1+$3}} {{print $4,$5}}'".format(mzml=mzmlfile)  # double curly brackets, use ' only for shell comand escapes

process = subprocess.Popen([wesc], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)  # shell=True for piped commands
stdout, stderr = process.communicate()
windows = list()
window_circular = list(csv.reader(StringIO(stdout.decode()), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC ) )

custom_windows_file = '{}_windows.txt'.format(os.path.splitext(os.path.basename(mzmlfile))[0])    
dyn = False
overlaps = list()
windowcount = 0
last = window_circular[-1]
with open(custom_windows_file, 'w') as f:
	f.write('start\tend')
	for frst, nxt in sliding_window_iter(window_circular, 2):
		if frst[1] > nxt[1]:
			last = [frst[0], frst[1]]
			break
		overlaps.append(frst[1]-nxt[0])
		if frst[1]-frst[0] != nxt[1]-nxt[0]:
			dyn = True 
		f.write("\n{}\t{}".format(frst[0], nxt[0] if frst[1]>nxt[0] else frst[1]))
		windowcount += 1
	f.write("\n{}\t{}".format(last[0],last[1]))
	windowcount += 1

print("Window count: {}".format(windowcount))
if dyn:
	print("Attention: dynamic window size!")
else:
	print("Window size: {}".format(last[1]-last[0]))

if len(set(overlaps)) == 1:
	print("Overlap: {}".format(overlaps[0]))
else:
	print("Attention: varying overlap size!")
