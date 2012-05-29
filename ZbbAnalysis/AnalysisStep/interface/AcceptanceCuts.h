#ifndef ZbbAnalysis_AnalysisStep_AcceptanceCuts_h
#define ZbbAnalysis_AnalysisStep_AcceptanceCuts_h

struct AcceptanceCuts {
  
  double bjetEtaMax_,bjetPtMin_,muEtaMax_,muPtMin_,eleEtaMax_,elePtMin_,betaCut_,betaStarCut_;
  
  AcceptanceCuts(){;}
  
  void clear(){
    bjetEtaMax_ =  0.;
    bjetPtMin_  =  0.; 
    muEtaMax_   =  0.;
    muPtMin_    =  0.;
    eleEtaMax_  =  0.;
    elePtMin_   =  0.;
    betaCut_    =  0.;
    betaStarCut_ = 0.;
  }

  void set(double m_bjetEtaMax,double m_bjetPtMin,double m_muEtaMax,double muPtMin,double m_eleEtaMax,double m_elePtMin,double m_betaCut,double m_betastarCut){     
    bjetEtaMax_ =  m_bjetEtaMax;
    bjetPtMin_  =  m_bjetPtMin; 
    muEtaMax_   =  m_muEtaMax;
    muPtMin_    =  muPtMin;
    eleEtaMax_  =  m_eleEtaMax;
    elePtMin_   =  m_elePtMin;
    betaCut_    =  m_betaCut;
    betaStarCut_=  m_betastarCut;
  }

  void setDefault(){
    bjetEtaMax_=2.1;
    bjetPtMin_=25;
    muEtaMax_=2.1;
    muPtMin_=20;
    eleEtaMax_=2.5;
    elePtMin_=25;
    betaCut_=0.15;	 
    betaStarCut_=0.85;
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

  void setBetaCut(double m_betaCut){     
    betaCut_   =  m_betaCut;
  }

  void setBetaStarCut(double m_betaStarCut){     
    betaStarCut_   =  m_betaStarCut;
  }
  
  void print() {
    std::cout<<"* bjet eta max: "<<bjetEtaMax_ <<" bjet pt min: "<<bjetPtMin_  <<std::endl; 
    std::cout<<"* muon eta max: "<<muEtaMax_   <<" muon pt min: "<<muPtMin_    <<std::endl;
    std::cout<<"* ele  eta max: "<<eleEtaMax_  <<" ele  pt min: "<<elePtMin_   <<std::endl;
    std::cout<<"* beta cut >  : "<<betaCut_    <<" beta* cut  : "<<betaStarCut_<<std::endl;
  }
};

#endif
