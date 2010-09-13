#! /bin/sh
# ./createFitDatasets.sh assumes that merged trees (one per species as defined in createFitDatasets.cc) exist in results/datasets_trees/
# creates the fit RooDataSets in directory results/datasets

mkdir -p results/datasets

root -l -b <<EOF

.include /afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.25.02-cms6/include
.L createFitDatasetsEG.cc++
createAll();

.q

EOF
