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
for i in range(9):
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
cycles_med = cycles[l/2]
cycles_min = cycles[0]
cycles_max = cycles[-1]
cycles_avg = sum(cycles)/len(cycles)

print "n: %5d, min: %10d, med: %10d, avg: %10d, max: 10%d"%(n, cycles_min, cycles_med, cycles_avg, cycles_max)
