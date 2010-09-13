#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 7:
    print "usage cmst3_submit_manyfilesperjob.py process dataset njobs applicationName sample queue"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
inputlist = "../cmst3_35X/MC/"+process+"/"+dataset+".list"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
queue = sys.argv[6]
ijobmax = int(sys.argv[3])
application = sys.argv[4]
sample = sys.argv[5]
# to write on local disks
################################################
castordir = "/castor/cern.ch/user/e/emanuele/VecBos33X/SCStudy/"
diskoutputdir = "/cmsrm/pc18/crovelli/data/Egamma.3.5.X/"
outputmain = castordir+"/"+process+"/"+output+"/"+sample
diskoutputmain = diskoutputdir+"/"+process+"/"+output+"/"+sample
# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir -p "+process+"/"+output+"/"+sample)
os.system("mkdir -p "+process+"/"+output+"/"+sample+"/log/")
os.system("mkdir -p "+process+"/"+output+"/"+sample+"/input/")
os.system("mkdir -p "+process+"/"+output+"/"+sample+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
    os.system("rfmkdir -p "+outputmain)
    os.system("rfmkdir -p "+outputroot)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputmain)

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
    inputfilename = pwd+"/"+process+"/"+output+"/"+sample+"/input/input_"+str(ijob)+".list"
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
    outputname = process+"/"+output+"/"+sample+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp -r '+pwd+'/config $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+sample+"\n")
#    if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root\n")
#    else:  
    outputfile.write('ls *.root | xargs -i rfcp {} '+outputroot+'\n')
    outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm21:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+process+"/"+output+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    ijob = ijob+1
    continue
