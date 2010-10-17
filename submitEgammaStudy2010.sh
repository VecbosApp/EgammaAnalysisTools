#!/bin/sh

echo "SUBMITTING PHOTON + JETS..."
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt0to15    1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt120to170 1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt15to20   1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt170to300 1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt20to30   1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt300to500 1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt30to50   1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt500toInf 1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt50to80   1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet PhotonJet_Pt80to120  1 EgammaAnalysis 4 8nh
echo "DONE PHOTON + JETS."

echo "submitting W+jets, Z+jets"
python cmst3_submit_manyfilesperjob.py WJetsMADGRAPH WJets-madgraph 20 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py ZJetsMADGRAPH ZJets-madgraph  2 EgammaAnalysis 4 8nh
echo "done with W+jets"

echo "submitting QCD di-jets..."
python cmst3_submit_manyfilesperjob.py QCD QCD_BCtoE_Pt20to30        5 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_BCtoE_Pt30to80        5 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_BCtoE_Pt80to170       5 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_EMEnriched_Pt20to30 106 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_EMEnriched_Pt30to80 121 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_EMEnriched_Pt80to170 10 EgammaAnalysis 4 8nh
echo "done QCD di-jets."

echo "submitting QCD not enriched..." 
python cmst3_submit_manyfilesperjob.py QCD QCD_Pt-20_TuneD6T 10 EgammaAnalysis 4 8nh
echo "done QCD not enriched" 

echo "SUBMITTING TTBAR and single top"
python cmst3_submit_manyfilesperjob.py TTbar TTbarJets-madgraph 3 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py SingleTop SingleTop_sChannel-madgraph  1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py SingleTop SingleTop_tChannel-madgraph  1 EgammaAnalysis 4 8nh
python cmst3_submit_manyfilesperjob.py SingleTop SingleTop_tWChannel-madgraph 1 EgammaAnalysis 4 8nh
echo "DONE TTBAR and SINGLE TOP."

