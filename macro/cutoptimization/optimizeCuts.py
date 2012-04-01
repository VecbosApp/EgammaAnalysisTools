#!/usr/bin/python
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-b", "--biased",
                  dest="biased", default=0, type="int",
                  help="optimize for triggering electrons. Signal is Z->ee with MC match; background is fake rate triggers (odd events)")
(options, args) = parser.parse_args()

if(options.biased==1):
    print "===> Optimizing cuts for trigger biased electrons..."
    signalFile="../results_data/electrons_zeemc.root"
    backgroundFile='../results_data/fakes.root'
    signalFriendFile="../results_data/electrons_zeemc_hzzisoFriend.root"
    backgroundFriendFile='../results_data/fakes_hzzisoFriend.root'
    additionalCut='&& DenomFakeSmurf==1'
else:
    print "===> Optimizing cuts for unbiased electrons..."
    signalFile="../results_data/electrons_zeemc.root"
    backgroundFile="../results_data/fakes-unbiased-wlnu.root"
    signalFriendFile="../results_data/electrons_zeemc_hzzisoFriend.root"
    backgroundFriendFile="../results_data/fakes-unbiased-wlnu_hzzisoFriend.root"
    additionalCut=''


if(options.biased==1):
    cuts        = ['&& abs(eta)<0.8 && pt>10 && pt<20',
                   '&& abs(eta)>=0.8 && abs(eta)<1.479 && pt>10 && pt<20',
                   '&& abs(eta)>=1.479 && pt>10 && pt<20',
                   '&& abs(eta)<0.8 && pt>=20',
                   '&& abs(eta)>=0.8 && abs(eta)<1.479 && pt>=20',
                   '&& abs(eta)>=1.479 && pt>=20']
else:
    cuts        = ['&& abs(eta)<0.8 && pt>5 && pt<10',
                   '&& abs(eta)>=0.8 && abs(eta)<1.479 && pt>5 && pt<10',
                   '&& abs(eta)>=1.479 && pt>5 && pt<10',
                   '&& abs(eta)<0.8 && pt>=10',
                   '&& abs(eta)>=0.8 && abs(eta)<1.479 && pt>=10',
                   '&& abs(eta)>=1.479 && pt>=10']
    
if(options.biased==1):
    outputs     = ['TMVA_Trg_fullIso_inEB_Pt10To20.root',
                   'TMVA_Trg_fullIso_outEB_Pt10To20.root',
                   'TMVA_Trg_fullIso_EE_Pt10To20.root',
                   'TMVA_Trg_fullIso_inEB_Pt20.root',
                   'TMVA_Trg_fullIso_outEB_Pt20.root',
                   'TMVA_Trg_fullIso_EE_Pt20.root']
else:
    outputs     = ['TMVA_NoTrg_fullIso_inEB_Pt5To10.root',
                   'TMVA_NoTrg_fullIso_outEB_Pt5To10.root',
                   'TMVA_NoTrg_fullIso_EE_Pt5To10.root',
                   'TMVA_NoTrg_fullIso_inEB_Pt10.root',
                   'TMVA_NoTrg_fullIso_outEB_Pt10.root',
                   'TMVA_NoTrg_fullIso_EE_Pt10.root']

if(options.biased==1):
    weightfiles  = ['OptimizedCuts_Trg_fullIso_weights_inEB_Pt10To20.xml',
                    'OptimizedCuts_Trg_fullIso_weights_outEB_Pt10To20.xml',
                    'OptimizedCuts_Trg_fullIso_weights_EE_Pt10To20.xml',
                    'OptimizedCuts_Trg_fullIso_weights_inEB_Pt20.xml',
                    'OptimizedCuts_Trg_fullIso_weights_outEB_Pt20.xml',
                    'OptimizedCuts_Trg_fullIso_weights_EE_Pt20.xml']
else:
    weightfiles  = ['OptimizedCuts_NoTrg_fullIso_weights_inEB_Pt5To10.xml',
                    'OptimizedCuts_NoTrg_fullIso_weights_outEB_Pt5To10.xml',
                    'OptimizedCuts_NoTrg_fullIso_weights_EE_Pt5To10.xml',
                    'OptimizedCuts_NoTrg_fullIso_weights_inEB_Pt10.xml',
                    'OptimizedCuts_NoTrg_fullIso_weights_outEB_Pt10.xml',
                    'OptimizedCuts_NoTrg_fullIso_weights_EE_Pt10.xml']
    

for i in range(len(cuts)):
    if len(cuts) != len(outputs):
        raise RuntimeError('cut set does not correspond to output set!')
    fullcut=cuts[i]+additionalCut
    if(options.biased==1):
        fulloutput='biased_'+outputs[i]
        fullweightfile='weights/biased_'+weightfiles[i]
    else:
        fulloutput='unbiased_'+outputs[i]
        fullweightfile='weights/unbiased_'+weightfiles[i]
    print "Now running optimization with Cut = "+fullcut+". Output will be put in "+fulloutput+"  and xml file in "+fullweightfile
    print "now starting seriously........."
    os.system('python CutsOptimization.py -m Cuts -a "'+fullcut+'" -o '+fulloutput+' -i '+signalFile+' -j '+backgroundFile+' -f '+signalFriendFile+' -g '+backgroundFriendFile)
    os.system('mv weights/TMVAClassification_Cuts.weights.xml '+fullweightfile)
        
print "DONE!"
