#ifndef AuxCode_SLHCSimPhase2_RecHit_h
#define AuxCode_SLHCSimPhase2_RecHit_h

// STD
#include <memory>
#include <string>
#include <iostream>

static const int DIGIMAX = 1000000;
static const int CLUSTERMAX = 100000;
static const int DGPERCLMAX = 100;  

//--- Structures for ntupling
struct RecHit 
{
  int pdgid;
  int process;
  float q;
  float x;
  float y;
  float xx;
  float xy;
  float yy;
  float row;
  float col;
  float gx;
  float gy;
  float gz;
  int subid,module;
  int layer,ladder;           // BPix
  int disk,blade,panel,side;  // FPix
  int nsimhit;
  int spreadx,spready;
  float hx, hy;
  float tx, ty, tz;
  float theta, phi;
  
  // digis
  int fDgN;  
  int fDgRow[DIGIMAX], fDgCol[DIGIMAX];
  int fDgDetId[DIGIMAX];
  //   int fDgRoc[DIGIMAX], fDgRocR[DIGIMAX], fDgRocC[DIGIMAX];
  float fDgAdc[DIGIMAX], fDgCharge[DIGIMAX];
  int fDgClI[DIGIMAX];
  
  // clusters
  int fClN;  
  int fClDgN[CLUSTERMAX], fClDgI[CLUSTERMAX][DGPERCLMAX];
  
  void init() {
    float dummy_float = 9999.0;
    
    pdgid = 0;
    process = 0;
    q = dummy_float;
    x = dummy_float;
    y = dummy_float;
    xx = dummy_float;
    xy = dummy_float; 
    yy = dummy_float;
    row = dummy_float;
    col = dummy_float;
    gx = dummy_float;
    gy = dummy_float;
    gz = dummy_float;
    nsimhit = 0;
    subid=-99;
    module=-99;
    layer=-99;
    ladder=-99;
    disk=-99;
    blade=-99;
    panel=-99;
    side=-99;
    spreadx=0;
    spready=0;
    hx = dummy_float;
    hy = dummy_float;
    tx = dummy_float;
    ty = dummy_float;
    tz = dummy_float;
    theta = dummy_float;
    phi = dummy_float;
    
    fClN = CLUSTERMAX;  
    for (int i = 0; i < fClN; ++i) {
      for (int j = 0; j < DGPERCLMAX; ++j) {
	fClDgI[i][j] = -9999;
      }
      fClDgN[i] = 0;
    }  
    fClN = 0; 
    
    fDgN = DIGIMAX;
    for (int i = 0; i < fDgN; ++i) {
      fDgRow[i] = fDgCol[i] = -9999;
      fDgAdc[i] = fDgCharge[i] = -9999.;
      fDgClI[i] = -9999;
      //    fDgRoc[i] = fDgRocR[i] = fDgRocC[i] = -9999;
    }
    fDgN = 0;  
  }
};

#endif

