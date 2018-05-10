#!/usr/bin/python

from os.path import join
from math import sqrt, cos, sin
import sys

##
# Read in the commadn line options
##

from optparse import OptionParser
usage = "usage: %prog [options] <trigger file>"
parser = OptionParser(usage=usage)
parser.add_option("-p", "--momentum", default="10.0",      help="Momentum in GeV (%default)", metavar="P")
parser.add_option("-t", "--time",     default="0.0",       help="Start Time (%default)", metavar="T")
parser.add_option(      "--pdg",      default="13",        help="Particle PDG code (%default, mu-)")
parser.add_option("-n", "--number",   default="0",         help="Number of particles", metavar="N")
parser.add_option(      "--direction",default="y",         help="Direction of the scan")
parser.add_option(      "--min",      default="0",         help="Minimum value")
parser.add_option(      "--max",      default="220",       help="Maximum value")
(options, args) = parser.parse_args()


## Particle Masses
massref = { 11:0.000511, 13:0.105 }

##
# Function which creates events
##
def MakeEvent(outfile,
              status = 1,
              pdg = int(options.pdg),
              mother1 = 0, mother2 = 0, daughter1 = 0, daughter2 = 0,
              momentum = float(options.momentum),
              direction = (0.0, 0.0, 1.0),
              position = (50., 0.0, 0.0),
              time = float(options.time)  ):
    global eventno

    momx, momy, momz = tuple([ momentum*x for x in direction ])
    mass   = massref[abs(pdg)]
    energy = sqrt(momentum**2 + mass**2)
    posx, posy, posz = position
    
    print >>outfile, eventno, 1
    print >>outfile, status, pdg,
    print >>outfile, mother1, mother2, daughter1, daughter2,
    print >>outfile, momx, momy, momz, energy, mass,
    print >>outfile, posx, posy, posz, time
    eventno += 1



vecfile = open(options.direction+"scan.vec","w")

N = int(options.number)
minimum = float(options.min)
maximum = float(options.max)
step = (maximum - minimum)/(N - 1)


for eventno in range(N):
    val = minimum + eventno * step

    if options.direction == "x":
        pos = ( val, 0., 0.)
        dir = (0., 0., 1.)
    elif options.direction == "y":
        pos = (50., val, 0.)
        dir = (0., 0., 1.)
    elif options.direction == "z":
        pos = (50., 115., val)
        dir = (0., -1., 0.)
    else:
        print "Unknown direction", options.direction
        sys.exit(1)
    
    MakeEvent(vecfile, position=pos, direction=dir )

vecfile.close()
