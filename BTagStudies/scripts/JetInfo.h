#ifndef JetInfo_h
#define JetInfo_h

#include "TObject.h"

// struct to store the relevant per-jet infos

// primary vertex related
class PVInfo {
 public:
  Float_t PVx,PVy,PVz,PVChiSq,PVndof,PVNormChiSq;  
  PVInfo() {}
  ClassDef(PVInfo,1)
};

// secondary vertex related
class SVInfo {
 public:
  Float_t SV3dDistance,SV3dDistanceError, SV2dDistance,SV2dDistanceError,
          SVChi2,SVDegreesOfFreedom,SVNormChi2,SVMass; 
  Int_t   SVtotCharge,SVnVertices,SVnVertexTracks,SVnVertexTracksAll;
  SVInfo() {}
  ClassDef(SVInfo,1)
}; 

// track related
class TRKInfo {
 public:
  Int_t quality;
  Float_t IP3d,IP3dError,IP3dDecayLength,pT,eta,phi;
  TRKInfo() {}
  ClassDef(TRKInfo,1)
};

// jet info related
class JetInfo : public TObject {

 public:
  Float_t pt,eta,phi;
  Int_t   nTracks;
  Int_t   jetId;
  Int_t   MCTrueFlavor;
  Float_t tche,tchp,ssvhe,ssvhp,csv,jp,jbp;
  PVInfo pv;
  SVInfo sv;
  TRKInfo trk[4];

  JetInfo();
  JetInfo(const JetInfo& j);
  ~JetInfo() {;}

  ClassDef(JetInfo,1)
}; 
#endif
