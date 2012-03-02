#include "include/HZZEleIDSelector.hh"
#include "EgammaAnalysisTools/include/eIDMVACuts.h"
#include <math.h>

int HZZEleIDSelector::etabin(float eta) {
  if(fabs(eta)<1.0) return 0;
  else if(fabs(eta)<1.4442) return 1;
  else if(fabs(eta)<1.566) return 2;
  else if(fabs(eta)<2.0) return 3;
  else return 4;
}

int HZZEleIDSelector::ptbin(float pt) {
  if(pt<20) return 0;
  else return 1;
}

bool HZZEleIDSelector::output(float pt, float eta, float bdt, float iso, 
			      HZZEleIDSelector::wpfulliso WP) {
  int etab=etabin(eta);
  int ptb=ptbin(pt);
  float bdtcut=cutbdtfulliso[ptb][etab][WP];
  float isocut=cutfulliso[ptb][etab][WP];
  return (bdt>bdtcut && iso<isocut);
}

bool HZZEleIDSelector::output(float pt, float eta, float bdt, float iso, 
			      HZZEleIDSelector::wpchiso WP) {
  int etab=etabin(eta);
  int ptb=ptbin(pt);
  float bdtcut=cutbdtchiso[ptb][etab][WP];
  float isocut=cutchiso[ptb][etab][WP];
  return (bdt>bdtcut && iso<isocut);
}
