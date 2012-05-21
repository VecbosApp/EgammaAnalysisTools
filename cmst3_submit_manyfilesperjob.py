#! /usr/bin/env python
import os, sys, commands, time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 10:
    print "usage cmst3_submit_manyfilesperjob.py process dataset njobs applicationName sample queue dirname isMC epoch"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
#settingfile = "config/RSZZsettings.txt"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = sys.argv[6]
#ijobmax = 40
ijobmax = int(sys.argv[3])
#application = "VecbosApp"
application = sys.argv[4]
sample = sys.argv[5]
dirname = sys.argv[7]
isMC = int(sys.argv[8])
epoch = int(sys.argv[9])
if isMC != 0:
    if epoch == 2011:
        inputlist = "cmst3_42X/MC/"+process+"/"+dataset+".list"
    else:
        inputlist = "cmst3_52X/MC/"+process+"/"+dataset+".list"
else:
    if epoch == 2011:
        inputlist = "cmst3_42X/"+process+"/"+dataset+".list"
    else:
        inputlist = "cmst3_52X/Data/"+process+"/"+dataset+".list"
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
castordir = "none"
if epoch==2011: 
    diskoutputdir = "/cmsrm/pc24_2/emanuele/data/Egamma4.2.X/"+dirname+"/"
else:
    diskoutputdir = "/cmsrm/pc24_2/emanuele/data/Egamma5.2.X/"+dirname+"/"
outputmain = castordir+"/"+process+"/"+output+"/"+sample
diskoutputmain = diskoutputdir+"/"+process+"/"+output+"/"+sample
# prepare job to write on the cmst3 cluster disks
################################################
os.system("rm -rf "+dirname+"/"+process+"/"+output+"/"+"/"+sample)
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample)
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/log/")
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/input/")
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
#    os.system("rfrm -r "+outputmain)
    os.system("rfmkdir -p "+castordir)
    os.system("rfmkdir -p "+outputmain)
    os.system("rfmkdir -p "+outputroot)
    os.system("rfchmod 777 "+castordir)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm24 rm -rf "+diskoutputmain)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm24 mkdir -p "+diskoutputdir)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm24 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
#print inputlist
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
os.system("cp -r config "+dirname)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)


    inputfile.close()

    # prepare the script to run
    outputname = dirname+"/"+process+"/"+output+"/"+"/"+sample+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
#    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp -r '+pwd+'/data/ $WORKDIR\n')
    outputfile.write('cp -r '+pwd+'/elebdtweights/ $WORKDIR\n')
    outputfile.write('cp '+pwd+'/pdfs_MC.root $WORKDIR\n')
    outputfile.write('cp -r '+pwd+"/"+dirname+'/config $WORKDIR\n')
    if epoch==2011:
        outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
        outputfile.write('cd /afs/cern.ch/user/e/emanuele/scratch0/higgs/CMSSW_4_2_8_patch7/\n')
    else:
        outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc462\n')
        outputfile.write('cd /afs/cern.ch/work/e/emanuele/releases/CMSSW_5_2_3/\n')        
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write('cd $WORKDIR\n')
    if isMC != 0:
        outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+" -signal="+sample+"\n")
    else:
        outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+" -signal="+sample+" --isData \n")
#    if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root\n")
#    else:  
#    outputfile.write('ls *.root | xargs -i rfcp {} '+outputroot+'\n')
    outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm24:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    time.sleep(2)
    ijob = ijob+1
    continue
