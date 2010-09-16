#! /bin/sh 
# ./mergeTrees.sh expects results of single sample files in results/merged/
# it creates a merged root file for each single specie of the fit in results/datasets

mkdir -p results/trees

# W alone
cp results/merged/WJetsMADGRAPH_tree.root results/trees

# Z alone
cp results/merged/ZJetsMADGRAPH_tree.root results/trees

# merging all the enriched QCD
hadd results/trees/QCD_tree.root  results/merged/QCD_*_tree.root 

# for a test: QCD not enriched pT>20 alone
cp results/merged/QCD_Pt-20_TuneD6T_tree.root results/trees

# merging all gamma+jets
hadd results/trees/GammaJets_tree.root results/merged/PhotonJet_*_tree.root 

# merging all the top samples
hadd results/trees/TOP_tree.root results/merged/TTbar_madgraph_tree.root  results/merged/SingleTop_*_tree.root  
