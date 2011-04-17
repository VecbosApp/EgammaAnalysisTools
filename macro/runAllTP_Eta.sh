#!/bin/bash

echo "running WP95..."
cmsRun testTP_Eta.py >& eta_WP95.log
echo "done WP95."

echo "running WP90..."
sed -i 's/WP95/WP90/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_WP90.log
echo "done WP90."

echo "running WP85..."
sed -i 's/WP90/WP85/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_WP85.log
echo "done WP85."

echo "running WP80..."
sed -i 's/WP85/WP80/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_WP80.log
echo "done WP80."

echo "running WP70..."
sed -i 's/WP80/WP70/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_WP70.log
echo "done WP70."

echo "running LHVeryLoose..."
sed -i 's/WP70/LHVeryLoose/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_LHVeryLoose.log
echo "done LHVeryLoose."

echo "running LHLoose..."
sed -i 's/LHVeryLoose/LHLoose/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_LHLoose.log
echo "done LHLoose."

echo "running LHMedium..."
sed -i 's/LHLoose/LHMedium/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_LHMedium.log
echo "done LHMedium."

echo "running LHTight..."
sed -i 's/LHMedium/LHTight/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_LHTight.log
echo "done LHTight."

echo "running LHHyperTight..."
sed -i 's/LHTight/LHHyperTight/g' testTP_Eta.py
cmsRun testTP_Eta.py >& eta_LHHyperTight.log
echo "done LHHyperTight."

echo "now making the default python file..."
sed -i 's/LHHyperTight/WP95/g' testTP_Eta.py
echo "VERY DONE."
