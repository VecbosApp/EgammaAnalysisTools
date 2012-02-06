#!/bin/bash

echo "merging data..."
hadd results_data/merged.root results_data/DoubleElectron*root
echo "DONE merging data rootples."

echo "merging MC (not weighted by xsec: binning in Pt is only important)..."
hadd results/merged.root results/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/4/DYToEE*root results/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia/4/DYToEE*root
echo "DONE merging MC rootples."
