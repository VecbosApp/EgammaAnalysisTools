#!/usr/bin/python
import sys
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--recompile",
                  action="store_true", dest="recompile", default=False,
                  help="recompile the fake rate estimator")
parser.add_option("-d", "--drawonly",
                  action="store_true", dest="drawonly", default=False,
                  help="do not run. Only draw the results based on existing ROOT files")

(options, args) = parser.parse_args()

if(options.drawonly == False):
    if(options.recompile):
        print 'recompiling...'
        os.system('make clean')
        os.system('make fakerate')
        print 'done complining.'

    os.system('./fakerate fr')
    os.system('hadd -f fakerates_trigger.root fr_trigger-EleMisid*')
    os.system('hadd -f fakerates_zee1fake.root fr_zee1fake-EleMisid*')

from ROOT import gROOT
gROOT.LoadMacro('drawFR.cc+')
from ROOT import drawIdsBiased
drawIdsBiased()
from ROOT import drawIdsUnbiased
drawIdsUnbiased()
from ROOT import drawIdsUnbiasedNewWPs
drawIdsUnbiasedNewWPs()
from ROOT import drawIdsUnbiasedIsoBins
drawIdsUnbiasedIsoBins()

# make the 2D map and save in a ROOT file
gROOT.LoadMacro('doFRMaps.cc+')
from ROOT import doFRMapsHzz4l
doFRMapsHzz4l(0)
doFRMapsHzz4l(1)
