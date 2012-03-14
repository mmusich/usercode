#include "JetInfo.h"

ClassImp(JetInfo)

JetInfo::JetInfo() : TObject() {
}

JetInfo::JetInfo(const JetInfo& j) : TObject(j) {
  pt=j.pt;eta=j.eta;phi=j.phi;
  nTracks=j.nTracks;
  jetId=j.jetId;
  MCTrueFlavor=j.MCTrueFlavor;
  MCweight=j.MCweight;
  tche=j.tche;tchp=j.tchp;ssvhe=j.ssvhe;ssvhp=j.ssvhp;csv=j.csv;jp=j.jp;jbp=j.jbp;
  pv=j.pv;
  sv=j.sv;
  for (int i=0; i<4; i++)
    trk[i]=j.trk[i];
}


