#!/bin/sh

dir=$1

mkdir $dir/merged

hadd $dir/merged/EleEfficiencyEta.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyEta.root
hadd $dir/merged/EleEfficiencyEtaHighPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyEtaHighPt.root
hadd $dir/merged/EleEfficiencyEtaLowPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyEtaLowPt.root

hadd $dir/merged/EleEfficiencyPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyPt.root
hadd $dir/merged/EleEfficiencyPtBarrel.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyPtBarrel.root
hadd $dir/merged/EleEfficiencyPtEndcap.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyPtEndcap.root

hadd $dir/merged/EleEfficiencyPU.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyPU.root
hadd $dir/merged/EleEfficiencyPUHighPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyPUHighPt.root
hadd $dir/merged/EleEfficiencyPULowPt.root $dir/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_*_-EleEfficiencyPULowPt.root
