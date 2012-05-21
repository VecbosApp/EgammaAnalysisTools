#! /usr/bin/env python
import os, sys, commands, time
#################################################
### usage: submitData.py application prefix year
#################################################
if len(sys.argv) != 4:
    print "usage submitData.py application prefix year"
    sys.exit(1)
application = int(sys.argv[1])
prefix = sys.argv[2]
year = sys.argv[3]

# hardcoded production batch
production = "Collisions8TeV_V14_52X"

if application==8 or application==10:
    os.system("python cmst3_submit_manyfilesperjob.py "+production+" DoubleElectronRun2012A 5 EgammaAnalysis 4 8nh "+prefix+" 0 "+year)
elif application==11:
    os.system("python cmst3_submit_manyfilesperjob.py "+production+" SingleElectronRun2012A 5 EgammaAnalysis 4 8nh "+prefix+" 0 "+year)
    os.system("python cmst3_submit_manyfilesperjob.py "+production+" SingleMuonRun2012A 5 EgammaAnalysis 4 8nh "+prefix+" 0 "+year)
elif application==12:
    os.system("python cmst3_submit_manyfilesperjob.py "+production+" DoubleElectronRun2012A 5 EgammaAnalysis 4 8nh "+prefix+" 0 "+year)
    os.system("python cmst3_submit_manyfilesperjob.py "+production+" DoubleMuonRun2012A 5 EgammaAnalysis 4 8nh "+prefix+" 0 "+year)
else:
    print "check your application."

print "DONE."

