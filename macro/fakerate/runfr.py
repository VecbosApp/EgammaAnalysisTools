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
parser.add_option("-e", "--electrons",
                  action="store_true", dest="electrons", default=True,
                  help="run the fake rate for electrons")
parser.add_option("-m", "--muons",
                  action="store_false", dest="electrons", default=False,
                  help="run the fake rate for muons")

(options, args) = parser.parse_args()

if(options.drawonly == False):
    if(options.recompile):
        print 'recompiling...'
        os.system('make clean')
        os.system('make fakerate')
        os.system('make fakeratehzz4l')
        print 'done complining.'

    if(options.electrons == True):
        os.system('./fakerate fr')
        os.system('hadd -f fakerates_trigger.root fr_trigger-EleMisid*')
        os.system('hadd -f fakerates_zee1fake.root fr_zee1fake-EleMisid*')
    if(options.electrons == False):
        os.system('./fakeratehzz4l fr')
        os.system('hadd -f fakerates_zll1m.root fr_zllm-MuonMisid*')

from ROOT import gROOT
gROOT.LoadMacro('drawFR.cc+')
if(options.electrons == True):
    from ROOT import drawIdsBiased
    drawIdsBiased()
    from ROOT import drawIdsUnbiased
    drawIdsUnbiased()
    from ROOT import drawIdsUnbiasedNewWPs
    drawIdsUnbiasedNewWPs()
    from ROOT import drawIdsUnbiasedIsoBins
    drawIdsUnbiasedIsoBins()
else:
    from ROOT import drawMuonIdsUnbiased
    drawMuonIdsUnbiased()

# make the 2D map and save in a ROOT file
gROOT.LoadMacro('doFRMaps.cc+')
if(options.electrons == True):
    from ROOT import doEleFRMapsHzz4l
    doEleFRMapsHzz4l(0)
    doEleFRMapsHzz4l(1)
else:
    from ROOT import doMuFRMapsHzz4l
    doMuFRMapsHzz4l(0)
    doMuFRMapsHzz4l(1)
