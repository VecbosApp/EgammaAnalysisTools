#!/bin/sh

echo "submitting signals..."
python cmst3_submit_manyfilesperjob.py WJetsMADGRAPH WJetsMADGRAPH 206 EgammaAnalysisTools 0 8nh
python cmst3_submit_manyfilesperjob.py WJetsMADGRAPH WJetsMADGRAPH 206 EgammaAnalysisTools 2 8nh
python cmst3_submit_manyfilesperjob.py ZJetsMADGRAPH ZJetsMADGRAPH 23 EgammaAnalysisTools 1 8nh
python cmst3_submit_manyfilesperjob.py ZJetsMADGRAPH ZJetsMADGRAPH 23 EgammaAnalysisTools 3 8nh
echo "done with signals."

echo "submitting QCD di-jets..."
python cmst3_submit_manyfilesperjob.py QCD QCD_BCtoE_Pt20to30 28 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_BCtoE_Pt30to80 63 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_BCtoE_Pt80to170 33 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_EMEnriched_Pt20to30 281 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_EMEnriched_Pt30to80 276 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py QCD QCD_EMEnriched_Pt80to170 134 EgammaAnalysisTools 4 8nh
echo "done QCD di-jets."

echo "SUBMITTING PHOTON + JETS..."
python cmst3_submit_manyfilesperjob.py PhotonJet Pt0to15 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt120to170 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt15to20 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt170to300 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt20to30 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt300to500 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt30to50 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt500toInf 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt50to80 5 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py PhotonJet Pt80to120 5 EgammaAnalysisTools 4 8nh
echo "DONE PHOTON + JETS."

echo "SUBMITTING TTBAR (MC@NLO)..."
python cmst3_submit_manyfilesperjob.py TTbar mcatnlo 99 EgammaAnalysisTools 4 8nh
echo "DONE TTBAR (MC@NLO)."

echo "SUBMITTING SINGLE TOP..."
python cmst3_submit_manyfilesperjob.py SingleTop SingleTop_sChannel_madgraph 10 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py SingleTop SingleTop_tChannel_madgraph 10 EgammaAnalysisTools 4 8nh
python cmst3_submit_manyfilesperjob.py SingleTop SingleTop_tWChannel_madgraph 10 EgammaAnalysisTools 4 8nh
echo "DONE SINGLE TOP."

