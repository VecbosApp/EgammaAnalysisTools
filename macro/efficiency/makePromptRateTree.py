# usage: python makePromptRateTree.py
import os

print "starting..."
from ROOT import gROOT
gROOT.LoadMacro('makePromptRateTree.cc+')
from ROOT import makeList
makeList()
from ROOT import makeSmall
makeSmall()

os.system('rm elist.root')
print "DONE. Reduced file is ../results_data/electrons_hzzidbitsFriend_promptrate.root"
