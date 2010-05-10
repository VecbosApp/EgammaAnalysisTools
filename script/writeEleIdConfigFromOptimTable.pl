#! /usr/bin/perl

use Getopt::Std;

$table = "optimTable.txt";
open(CUTTABLE, "$table");

$ecal = "EB";

$optimtable = "optimTable.txt";
open(OPTIMFILE,"$optimtable") or die "cannot open $optimtable";

@rows=<OPTIMFILE>;
$nrows=$#rows;
$i=1;

for($i=1; $i<=$nrows; $i++) {
    if($rows[$i]=~/BARREL/) { $ecal = "EB"; }
    else { if($rows[$i]=~/ENDCAP/) { $ecal = "EE"; } }
    if($rows[$i]=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
        system("mkdir -p Vecbos$1");
        $file = "Vecbos"."$1"."/electronIDCuts"."$ecal".".txt";
        print "writing cuts for eff = $1 % in $ecal. File is: $file\n";
        open(CONFIGFILE,">$file"); 
        print CONFIGFILE "deta \t 0.0 \t $2 \n";
        print CONFIGFILE "dphiIn \t 0.0 \t $3 \n";
        print CONFIGFILE "sigmaEtaEta \t 0.0 \t $4 \n";
        print CONFIGFILE "hOverE \t 0.0 \t $5 \n";
        print CONFIGFILE "trackerIso \t -100 \t $6 \n";
        print CONFIGFILE "ecalIso \t -100 \t $7 \n";
        print CONFIGFILE "hcalIso \t -100 \t $8 \n";
        print CONFIGFILE "missingHits \t 0 \t 1 \n";
        print CONFIGFILE "eOverPin        0.0             999.\n";
        print CONFIGFILE "s9s25           0.0             1.0\n";
        print CONFIGFILE "eOverPout       0.0             1000\n";
        print CONFIGFILE "dphiOut         0.0             999.\n";
        print CONFIGFILE "invEminusInvP   0.0             1000\n";
        print CONFIGFILE "bremFraction    0.0             999.\n";
        print CONFIGFILE "sigmaPhiPhi     0.0             999.\n";
        print CONFIGFILE "likelihood      0.0             1.0\n";
    }
}
