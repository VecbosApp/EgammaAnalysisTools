#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTrees.sh

echo "Adding weights..."
root -l -b <<EOF

.L addWeightsToEgammaTree.cc++

addWeights("results/merged/WJetsMADGRAPH_tree.root",  0.00320548);
addWeights("results/merged/ZJetsMADGRAPH_tree.root",  0.00294516);
addWeights("results/merged/TTbar_madgraph_tree.root", 0.000122721);

addWeights("results/merged/QCD_EMEnriched_Pt20to30_tree.root",  0.0567193);
addWeights("results/merged/QCD_EMEnriched_Pt30to80_tree.root",  0.104823);
addWeights("results/merged/QCD_EMEnriched_Pt80to170_tree.root", 0.0267058);

addWeights("results/merged/QCD_BCtoE_Pt20to30_tree.root",  0.047912);
addWeights("results/merged/QCD_BCtoE_Pt30to80_tree.root",  0.158477);
addWeights("results/merged/QCD_BCtoE_Pt80to170_tree.root", 1.08628);

addWeights("results/merged/PhotonJet_Pt0to15_tree.root",  732.746);
addWeights("results/merged/PhotonJet_Pt15to20_tree.root",   1.05656);
addWeights("results/merged/PhotonJet_Pt20to30_tree.root",   0.953);
addWeights("results/merged/PhotonJet_Pt30to50_tree.root",   0.150182);
addWeights("results/merged/PhotonJet_Pt50to80_tree.root",   0.0248155);
addWeights("results/merged/PhotonJet_Pt80to120_tree.root",  0.00402609);
addWeights("results/merged/PhotonJet_Pt120to170_tree.root", 0.000690459);
addWeights("results/merged/PhotonJet_Pt170to300_tree.root", 0.000180215);
addWeights("results/merged/PhotonJet_Pt300to500_tree.root", 1.43579e-05);
addWeights("results/merged/PhotonJet_Pt500toInf_tree.root", 1.62229e-06);

addWeights("results/merged/SingleTop_sChannel-madgraph_tree.root",  4.7748e-06);
addWeights("results/merged/SingleTop_tChannel-madgraph_tree.root",  4.59681e-05);
addWeights("results/merged/SingleTop_tWChannel-madgraph_tree.root", 0.000644887);

EOF

