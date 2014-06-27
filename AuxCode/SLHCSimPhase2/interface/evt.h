#ifndef AuxCode_SLHCSimPhase2_evt_h
#define AuxCode_SLHCSimPhase2_evt_h

// STD
#include <memory>
#include <string>
#include <iostream>

 //--- Structures for ntupling:
struct evt
{
  int run;
  int evtnum;
  
  void init(){
    int dummy_int = 9999;
    run = dummy_int;
    evtnum = dummy_int;    
  }
};

#endif
