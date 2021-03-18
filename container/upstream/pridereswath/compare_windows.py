import argparse
import sys
parser = argparse.ArgumentParser(description='Compare files.')
parser.add_argument('--in', nargs='+', dest='filepaths', required=True, type=str, help='input file paths')
args = parser.parse_args()
if not args.filepaths:
	parser.print_help()
	import sys
	sys.exit()

import filecmp
from itertools import combinations

cmps = [filecmp.cmp(pair[0],pair[1]) for pair in combinations(args.filepaths, 2)]
if not all(cmps):
	sys.exit(1)
