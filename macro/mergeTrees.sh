#! /bin/sh
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results/merged

hadd results/merged/WJetsMADGRAPH_tree.root  results/WJetsMADGRAPH/WJets-madgraph/4/WJets-madgraph_*_tree.root 
hadd results/merged/ZJetsMADGRAPH_tree.root  results/ZJetsMADGRAPH/ZJets-madgraph/4/ZJets-madgraph_*_tree.root 
hadd results/merged/QCD_EMEnriched_Pt20to30_tree.root   results/QCD/QCD_EMEnriched_Pt20to30/4/QCD_EMEnriched_Pt20to30_*_tree.root
hadd results/merged/QCD_EMEnriched_Pt30to80_tree.root   results/QCD/QCD_EMEnriched_Pt30to80/4/QCD_EMEnriched_Pt30to80_*_tree.root
hadd results/merged/QCD_EMEnriched_Pt80to170_tree.root  results/QCD/QCD_EMEnriched_Pt80to170/4/QCD_EMEnriched_Pt80to170_*_tree.root
hadd results/merged/QCD_BCtoE_Pt20to30_tree.root        results/QCD/QCD_BCtoE_Pt20to30/4/QCD_BCtoE_Pt20to30_*_tree.root
hadd results/merged/QCD_BCtoE_Pt30to80_tree.root        results/QCD/QCD_BCtoE_Pt30to80/4/QCD_BCtoE_Pt30to80_*_tree.root
hadd results/merged/QCD_BCtoE_Pt80to170_tree.root       results/QCD/QCD_BCtoE_Pt80to170/4/QCD_BCtoE_Pt80to170_*_tree.root

hadd results/merged/QCD_Pt-20_TuneD6T_tree.root         results/QCD/QCD_Pt-20_TuneD6T/4/QCD_Pt-20_TuneD6T_*_tree.root

hadd results/merged/TTbar_madgraph_tree.root  results/TTbar/TTbarJets-madgraph/4/TTbarJets-madgraph_*_tree.root

hadd results/merged/PhotonJet_Pt0to15_tree.root    results/PhotonJet/PhotonJet_Pt0to15/4/PhotonJet_Pt0to15_*_tree.root
hadd results/merged/PhotonJet_Pt120to170_tree.root results/PhotonJet/PhotonJet_Pt120to170/4/PhotonJet_Pt120to170_*_tree.root
hadd results/merged/PhotonJet_Pt15to20_tree.root   results/PhotonJet/PhotonJet_Pt15to20/4/PhotonJet_Pt15to20_*_tree.root
hadd results/merged/PhotonJet_Pt170to300_tree.root results/PhotonJet/PhotonJet_Pt170to300/4/PhotonJet_Pt170to300_*_tree.root
hadd results/merged/PhotonJet_Pt20to30_tree.root   results/PhotonJet/PhotonJet_Pt20to30/4/PhotonJet_Pt20to30_*_tree.root
hadd results/merged/PhotonJet_Pt300to500_tree.root results/PhotonJet/PhotonJet_Pt300to500/4/PhotonJet_Pt300to500_*_tree.root
hadd results/merged/PhotonJet_Pt30to50_tree.root   results/PhotonJet/PhotonJet_Pt30to50/4/PhotonJet_Pt30to50_*_tree.root
hadd results/merged/PhotonJet_Pt500toInf_tree.root results/PhotonJet/PhotonJet_Pt500toInf/4/PhotonJet_Pt500toInf_*_tree.root
hadd results/merged/PhotonJet_Pt50to80_tree.root   results/PhotonJet/PhotonJet_Pt50to80/4/PhotonJet_Pt50to80_*_tree.root
hadd results/merged/PhotonJet_Pt80to120_tree.root  results/PhotonJet/PhotonJet_Pt80to120/4/PhotonJet_Pt80to120_*_tree.root

hadd results/merged/SingleTop_sChannel-madgraph_tree.root  results/SingleTop/SingleTop_sChannel-madgraph/4/SingleTop_sChannel-madgraph_*_tree.root
hadd results/merged/SingleTop_tChannel-madgraph_tree.root  results/SingleTop/SingleTop_tChannel-madgraph/4/SingleTop_tChannel-madgraph_*_tree.root
hadd results/merged/SingleTop_tWChannel-madgraph_tree.root results/SingleTop/SingleTop_tWChannel-madgraph/4/SingleTop_tWChannel-madgraph_*_tree.root

