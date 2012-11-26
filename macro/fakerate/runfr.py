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

(options, args) = parser.parse_args()

if(options.drawonly == False):
    if(options.recompile):
        print 'recompiling...'
        os.system('make clean')
        os.system('make fakerate')
        print 'done complining.'

    if(options.electrons == True):
        os.system('./fakerate fakes  fr')
        os.system('hadd -f fakerates_trigger.root fr_trigger-EleMisid*')
        os.system('./fakerate fakes-ewksub-wlnu ewkw')
        os.system('hadd -f subewkw_trigger.root ewkw_trigger-EleMisid*')
        os.system('./fakerate fakes-ewksub-zee ewkz')
        os.system('hadd -f subewkz_trigger.root ewkz_trigger-EleMisid*')

from ROOT import gROOT
gROOT.LoadMacro('drawFR.cc+')
if(options.electrons == True):
    from ROOT import drawIdsBiased
    drawIdsBiased()

