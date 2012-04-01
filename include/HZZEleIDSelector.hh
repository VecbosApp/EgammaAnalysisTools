#ifndef HZZEleIDSelector_H
#define HZZEleIDSelector_H

class HZZEleIDSelector {
public:

  HZZEleIDSelector() {}
  ~HZZEleIDSelector() {}
  enum wpfulliso {
    kWP95 = 0,
    kWP90, kWP85, kWP80, kWP70
  };

  enum wpchiso {
    kWP95ChIso = 0,
    kWP90ChIso, kWP85ChIso, kWP80ChIso, kWP70ChIso
  };

  enum mvatype {
    kMVABiased = 0,
    kMVAUnbiased
  };

  bool output(float pt, float eta, float bdt, float iso, 
	      HZZEleIDSelector::wpfulliso WP, HZZEleIDSelector::mvatype type);
  bool output(float pt, float eta, float bdt, float chiso, 
	      HZZEleIDSelector::wpchiso WP, HZZEleIDSelector::mvatype type);

private:
  int etabin(float eta);
  int ptbinTrg(float pt);
  int ptbinNoTrg(float pt);

};

#endif
