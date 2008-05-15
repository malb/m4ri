#!/usr/bin/python
import sys
import commands

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", dest="n", default=1000)
parser.add_option("-c", dest="cutoff", default=2048)
parser.add_option("-t", dest="typ", default="red")

options, args = parser.parse_args()
n = int(options.n)
cutoff = int(options.cutoff)
typ = options.typ

cycles = []

pattern = "cpu cycles:"
for i in range(17):
    if typ == "red":
	out = commands.getoutput("./bench_reduction %d"%n)
    elif typ == "mul":
	out = commands.getoutput("./bench_multiplication %d %d"%(n,cutoff))
    elif typ == "add":
	out = commands.getoutput("./bench_addition %d"%n)
    else:
	raise TypeError, "Value for parameter typ %s not understood."%typ
    ii = out.rfind(pattern)
    cycles.append(int(out[ii+len(pattern):]))

# median as adviced in the cpucycles docs
cycles = sorted(cycles)
l = len(cycles)
total_cycles = cycles[l/2]

print "n: %5d, cpu cycles:  %10d"%(n, total_cycles)
