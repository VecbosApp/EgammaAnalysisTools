#!/bin/sh

prefix=$1

echo "SUBMITTING DATA..."
python cmst3_submit_manyfilesperjob_data38.py Data7TeV dataset_jetmettau 7 EgammaAnalysis 4 8nh dataJMet
echo "DONE WITH DATA."
