#!/bin/bash                                                                                                                                                                        
echo "merging data kine..."
hadd results_data_fakes/merged_kine.root results_data_fakes/DoubleElectron_*FakeKineTree.root
echo "DONE merging data kine."

echo "merging data ID..."
hadd results_data_fakes/merged.root results_data_fakes/DoubleElectron_*FakeIDTree.root
echo "DONE merging data ID."
