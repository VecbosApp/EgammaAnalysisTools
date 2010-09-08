#!/bin/sh

prefix=$1

echo "SUBMITTING DATA..."
python cmst3_submit_manyfilesperjob_data38.py Data7TeV dataset_eg 7 VecbosApp 4 8nh datiEgamma
echo "DONE WITH DATA."
