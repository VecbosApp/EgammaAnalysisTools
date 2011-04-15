#!/bin/sh

dir=$1

mkdir $dir/merged

hadd $dir/merged/EleMisidEta.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidEta.root
hadd $dir/merged/EleMisidEtaHighPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidEtaHighPt.root
hadd $dir/merged/EleMisidEtaLowPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidEtaLowPt.root

hadd $dir/merged/EleMisidPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidPt.root
hadd $dir/merged/EleMisidPtBarrel.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidPtBarrel.root
hadd $dir/merged/EleMisidPtEndcap.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidPtEndcap.root

hadd $dir/merged/EleMisidPU.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidPU.root
hadd $dir/merged/EleMisidPUHighPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidPUHighPt.root
hadd $dir/merged/EleMisidPULowPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleMisidPULowPt.root
