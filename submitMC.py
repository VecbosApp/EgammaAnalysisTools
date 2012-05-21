#! /usr/bin/env python
import os, sys, commands, time
#################################################
### usage: submitMC.py application prefix year
#################################################
if len(sys.argv) != 4:
    print "usage submitData.py application prefix year"
    sys.exit(1)
application = int(sys.argv[1])
prefix = sys.argv[2]
year = sys.argv[3]

# hardcoded production batch
production = "Summer12_V14_52X"

if application==8:
    os.system("python cmst3_submit_manyfilesperjob.py "+production+" DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball 5 EgammaAnalysis 4 8nh "+prefix+" 1 "+year)
else:
    print "check your application."

print "DONE."

