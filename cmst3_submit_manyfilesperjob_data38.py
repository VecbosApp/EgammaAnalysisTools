#! /usr/bin/env python
import os
import sys
import time
import commands
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 8:
    print "usage cmst3_submit_manyfilesperjob.py process dataset njobs applicationName sample queue dirname"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
inputlist = "../cmst3_38X/"+process+"/"+dataset+".list"
output = dataset
queue = sys.argv[6]
ijobmax = int(sys.argv[3])
application = sys.argv[4]
sample = sys.argv[5]
dirname = sys.argv[7]
################################################
# to write on local disks
################################################
castordir = "/castor/cern.ch/user/m/meridian/VecBos7TeV/"+dirname+"/"
diskoutputdir = "/cmsrm/pc18/crovelli/data/Egamma.3.8.X/"+dirname+"/"
outputmain = castordir+"/"+process+"/"+output+"/"+sample
diskoutputmain = diskoutputdir+"/"+process+"/"+output+"/"+sample
# prepare job to write on the cmst3 cluster disks
################################################
os.system("rm -rf "+process+"/"+output+"/"+dirname+"/config")
os.system("rm -rf "+process+"/"+output+"/"+dirname+"/"+sample)
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/"+sample)
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/"+sample+"/log/")
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/"+sample+"/input/")
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/"+sample+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
    os.system("rfrm -r "+castordir)
    os.system("rfmkdir -p "+castordir)
    os.system("rfmkdir -p "+outputmain)
    os.system("rfmkdir -p "+outputroot)
    os.system("rfchmod 777 "+castordir)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none":
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 rm -rf "+diskoutputdir)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputdir)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
#print inputlist
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0
os.system("cp -r config "+process+"/"+output+"/"+dirname)
while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+process+"/"+output+"/"+dirname+"/"+sample+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    
    command="echo "+inputfiles[0].strip()+" | awk -F 'vecbos_' '{print $2}' | awk -F '_' '{print $1}'"
    runNumberRef=commands.getoutput(command)
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop(0)
        if ntpfile != '':
            command="echo "+ntpfile.strip()+" | awk -F 'vecbos_' '{print $2}' | awk -F '_' '{print $1}'"
            runNumber=commands.getoutput(command)
            if (runNumber == runNumberRef):
                inputfile.write(ntpfile)
            else:
                inputfiles.insert(0,ntpfile)
                break
            
    inputfile.close()

    # prepare the script to run
    outputname = process+"/"+output+"/"+dirname+"/"+sample+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp -r '+pwd+"/"+process+"/"+output+"/"+dirname+'/config $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+" -signal="+sample+"\n")
#    if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root\n")
#    else:  
    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i rfcp {} '+outputroot+'\n')
    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm18:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+process+"/"+output+"/"+dirname+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    ijob = ijob+1
    time.sleep(3)
    continue
