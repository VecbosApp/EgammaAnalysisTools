#!/usr/bin/python
import os

cuts        = ['&& abs(eta)<1.0 && pt>10 && pt<20',
               '&& abs(eta)>=1.0 && abs(eta)<1.4442 && pt>10 && pt<20',
               '&& abs(eta)>=1.4442 && abs(eta)<1.566 && pt>10 && pt<20',
               '&& abs(eta)>=1.566 && abs(eta)<2.0 && pt>10 && pt<20',
               '&& abs(eta)>=2.0 && pt>10 && pt<20',
               '&& abs(eta)<1.0 && pt>=20',
               '&& abs(eta)>=1.0 && abs(eta)<1.4442 && pt>=20',
               '&& abs(eta)>=1.4442 && abs(eta)<1.566 && pt>=20',
               '&& abs(eta)>=1.566 && abs(eta)<2.0 && pt>=20',
               '&& abs(eta)>=2.0 && pt>=20']

outputs     = ['TMVA_eta10_Pt10To20.root',
               'TMVA_eta14442_Pt10To20.root',
               'TMVA_eta1566_Pt10To20.root',
               'TMVA_eta20_Pt10To20.root',
               'TMVA_eta25_Pt10To20.root',
               'TMVA_eta10_Pt20.root',
               'TMVA_eta14442_Pt20.root',
               'TMVA_eta1566_Pt20.root',
               'TMVA_eta20_Pt20.root',
               'TMVA_eta25_Pt20.root']

weightfile  = ['weights/OptimizedCuts_weights_eta10_Pt10To20.xml',
               'weights/OptimizedCuts_weights_eta14442_Pt10To20.xml',
               'weights/OptimizedCuts_weights_eta1566_Pt10To20.xml',
               'weights/OptimizedCuts_weights_eta20_Pt10To20.xml',
               'weights/OptimizedCuts_weights_eta25_Pt10To20.xml',
               'weights/OptimizedCuts_weights_eta10_Pt20.xml',
               'weights/OptimizedCuts_weights_eta14442_Pt20.xml',
               'weights/OptimizedCuts_weights_eta1566_Pt20.xml',
               'weights/OptimizedCuts_weights_eta20_Pt20.xml',
               'weights/OptimizedCuts_weights_eta25_Pt20.xml']


for i in range(len(cuts)):
    if len(cuts) != len(outputs):
        raise RuntimeError('cut set does not correspond to output set!')
    os.system('python CutsOptimization.py -m CutsGA -a "'+cuts[i]+'" -o '+outputs[i])
    os.system('mv weights/TMVAClassification_CutsGA.weights.xml '+weightfile[i])

print "DONE!"
