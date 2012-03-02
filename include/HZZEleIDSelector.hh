#ifndef HZZEleIDSelector_H
#define HZZEleIDSelector_H

class HZZEleIDSelector {
public:

  HZZEleIDSelector();
  ~HZZEleIDSelector();
  enum wpfulliso {
    kWP95 = 0,
    kWP90, kWP85, kWP80, kWP70
  };

  enum wpchiso {
    kWP95ChIso = 0,
    kWP90ChIso, kWP85ChIso
  };

  bool output(float pt, float eta, float bdt, float iso, 
	      HZZEleIDSelector::wpfulliso WP);
  bool output(float pt, float eta, float bdt, float chiso, 
	      HZZEleIDSelector::wpchiso WP);

private:
  int etabin(float eta);
  int ptbin(float pt);


};

#endif
