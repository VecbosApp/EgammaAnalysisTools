#ifndef LikelihoodMeasurements_h
#define LikelihoodMeasurements_h

struct LikelihoodMeasurements {
  int subdet; // 0=EB, 1=EE
  float pt;
  float deltaPhi,
    deltaEta,
    eSeedClusterOverPout,
    sSuperClusterOverP,
    hadronicOverEm,
    sigmaIEtaIEta,
    sigmaIPhiIPhi,
    fBrem;
  int nBremClusters;
};

#endif

