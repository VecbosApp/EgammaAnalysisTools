# usage: python makePromptRateTree.py
import os

print "starting..."
from ROOT import gROOT
gROOT.LoadMacro('makePromptRateTree.cc+')
from ROOT import makeList
# cut is "denom==1" to do the real promptrate.
# in this case change the root file in input to take the data
makeList("mcmatch==1")
from ROOT import makeSmall
makeSmall()

os.system('rm elist.root')
print "DONE. Reduced file is ../results_data/electrons_zeemc_hzzidbitsFriend_promptrate.root"
