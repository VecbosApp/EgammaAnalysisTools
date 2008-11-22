#! /usr/bin/env python
import os
# set parameters to use cmst3 batch 
#######################################
#dataset = "WjetsMADGRAPH"
#dataset = "ZjetsMADGRAPH"
#dataset = "ttjetsMADGRAPH"
#dataset = "QCDEmEnrichedPt20to30"
#dataset = "QCDEmEnrichedPt30to80"
dataset = "QCDEmEnrichedPt80to170"
inputlist = "../VecBosApp/cmst3_21X/"+dataset+".list"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
queue = "cmst3"
ijobmax = 100
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
castordir = "none"
outputmain = output
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
#look for the current directory
#######################################
os.system("\\rm tmp.log")
os.system("echo $PWD > tmp.log")
tmp = open("tmp.log")
pwd = tmp.readline()
tmp.close()
os.system("\\rm tmp.log")
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = output+"/input/input_"+str(ijob)+".list"
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
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    outputfile.write('cd '+pwd)
    if castordir != "none": outputfile.write('./EgammaAnalysis '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root -start=0 -stop=10\n")
    else:  outputfile.write('./EgammaAnalysisTagProbe '+inputfilename+' jobQCD80to170jets'+str(ijob)+" \n")
#    else:  outputfile.write('./EgammaAnalysis '+inputfilename+' pdfs'+str(ijob)+".root\n") 
#    else:  outputfile.write('./EgammaAnalysisPdfTagAndProbe '+inputfilename+' pdfsZTaP'+str(ijob)+".root\n")
#    else:  outputfile.write('./EgammaAnalysisPdfWjets '+inputfilename+' pdfsJets'+str(ijob)+".root\n")
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd[:-1]+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd[:-1]+"/"+outputname)
    ijob = ijob+1
    continue
