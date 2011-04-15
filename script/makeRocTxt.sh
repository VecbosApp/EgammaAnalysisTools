#!/bin/tcsh

echo "++ Producing the roc.txt file ++"
rm -rf roc.txt
touch roc.txt
cat efficiencies-HighPt.txt >> roc.txt
cat efficiencies-LowPt.txt >> roc.txt
cat elefakes-HighPt.txt >> roc.txt
cat elefakes-LowPt.txt >> roc.txt
