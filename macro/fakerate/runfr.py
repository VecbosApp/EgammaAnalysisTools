#!/usr/bin/python
import sys
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--recompile",
                  action="store_true", dest="recompile", default=False,
                  help="recompile the fake rate estimator")

(options, args) = parser.parse_args()

if(options.recompile):
    print 'recompiling...'
    os.system('make clean')
    os.system('make fakerate')
    print 'done complining.'

os.system('./fakerate fr')
os.system('hadd -f fakerates.root fr-EleMisid*')

from ROOT import gROOT
gROOT.LoadMacro('drawFR.cc+')
from ROOT import drawIds
drawIds()


