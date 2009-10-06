#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 5:
    print "usage cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue"
    sys.exit(1)
dataset = sys.argv[1]
inputlist = "../HiggsAnalysisTools/cmst3_32X/"+dataset+".list"
output = dataset
queue = sys.argv[4]
ijobmax = int(sys.argv[2])
application = sys.argv[3]
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
castordir = "/castor/cern.ch/user/XXX/"
diskoutputdir = "/cmsrm/pc21/crovelli/data/Like3.2.X/results1/"
outputmain = castordir+output
diskoutputmain = diskoutputdir+output
# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir "+output)
os.system("mkdir "+output+"/log/")
os.system("mkdir "+output+"/input/")
os.system("mkdir "+output+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
    os.system("rfmkdir "+outputmain)
    os.system("rfmkdir "+outputroot)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = pwd+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    if ijob != (ijobmax-1):
        for line in range(filesperjob):
            ntpfile = input.readline() 
            inputfile.write(ntpfile)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            inputfile.write(ntpfile)
            ntpfile = input.readline()
            continue
    inputfile.close()

    # prepare the script to run
    outputname = output+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    #outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    outputfile.write('cp -r '+pwd+'/config $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+" "+"\n")
    outputfile.write('rfcp *.root '+outputroot+'\n')
    outputfile.write('scp -o BatchMode=yes -o StrictHostKeyChecking=no *.root pccmsrm21:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    ijob = ijob+1
    continue
