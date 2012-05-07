#ifndef ZbbAnalysis_AnalysisStep_AcceptanceCuts_h
#define ZbbAnalysis_AnalysisStep_AcceptanceCuts_h

struct AcceptanceCuts {
  
  double bjetEtaMax_,bjetPtMin_,muEtaMax_,muPtMin_,eleEtaMax_,elePtMin_;
  
  AcceptanceCuts(){;}
  
  void clear(){
    bjetEtaMax_ =  0.;
    bjetPtMin_  =  0.; 
    muEtaMax_   =  0.;
    muPtMin_    =  0.;
    eleEtaMax_  =  0.;
    elePtMin_   =  0.;
  }

  void set(double m_bjetEtaMax,double m_bjetPtMin,double m_muEtaMax,double muPtMin,double m_eleEtaMax,double m_elePtMin){     
    bjetEtaMax_ =  m_bjetEtaMax;
    bjetPtMin_  =  m_bjetPtMin; 
    muEtaMax_   =  m_muEtaMax;
    muPtMin_    =  muPtMin;
    eleEtaMax_  =  m_eleEtaMax;
    elePtMin_   =  m_elePtMin;
  }

  void setDefault(){
    bjetEtaMax_=2.1;
    bjetPtMin_=25;
    muEtaMax_=2.1;
    muPtMin_=20;
    eleEtaMax_=2.5;
    elePtMin_=25;
  }

  // setter methods
  void setbJetEtaMax(double m_bjetEtaMax){     
    bjetEtaMax_ =  m_bjetEtaMax;
  }
  
  void setbJetPtMin(double m_bjetPtMin){     
    bjetPtMin_  =  m_bjetPtMin; 
  }
  
  void setmuEtaMax(double m_muEtaMax){     
    muEtaMax_   =  m_muEtaMax;
  }
  
  void setmuPtMin(double muPtMin){     
    muPtMin_    =  muPtMin;
  }
  
  void seteleEtaMax(double m_eleEtaMax){     
    eleEtaMax_  =  m_eleEtaMax;
  }
  
  void setelePtMin(double m_elePtMin){     
    elePtMin_   =  m_elePtMin;
  }
  
  void print() {
    std::cout<<"* bjet eta max: "<<bjetEtaMax_ <<" bjet pt min: "<<bjetPtMin_  <<std::endl; 
    std::cout<<"* muon eta max: "<<muEtaMax_   <<" muon pt min: "<<muPtMin_    <<std::endl;
    std::cout<<"* ele  eta max: "<<eleEtaMax_  <<" ele  pt min: "<<elePtMin_   <<std::endl;
  }
};

#endif
