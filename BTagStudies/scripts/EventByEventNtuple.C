/*******************************************************************************
*
* Take two ntuples file1 and file2
* Store a TGraphErrors to a file matchList
* x,y in the TGraphErrors are the entries in file1 and file2, that match
* ex is the run number
* ey is the event number
*
* Based on an original script by N.Kypreos
*
*******************************************************************************/
#include <iostream>
#include <vector>
#include "TBranch.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLeafI.h"
#include "TTree.h"
#include "TObjString.h"
#include "TObjArray.h"

using namespace std; 

class MyEvent {

public:

  UInt_t entry;
  UInt_t lumi;
  UInt_t run;
  UInt_t event;

  MyEvent(UInt_t const _entry, UInt_t const _run, UInt_t const _lumi, UInt_t const _event)
  {
    entry = _entry;
    run   = _run;
    lumi  = _lumi;
    event = _event;
  };
};

// forward declarations
vector<MyEvent> GetListToSort(TTree *tree);
vector<MyEvent> RemoveDuplicates(vector<MyEvent> const &vec);

//------------------------------------------------------------------------------
// EventByEventNtuple
//------------------------------------------------------------------------------
void EventByEventNtuple(const TString& filename_and_label1 = "RunRangeIdeal-165364.root=IDEAL", const TString& filename_and_label2 = "RunRangeNEWKB-165364.root=KandB")
{

  TObjArray *name_and_label_pair1 = filename_and_label1.Tokenize("=");
  TFile *file1 = TFile::Open(name_and_label_pair1->At(0)->GetName());
  TTree *tree1 = (TTree*)file1->Get("t");

  TObjArray *name_and_label_pair2 = filename_and_label2.Tokenize("=");
  TFile *file2 = TFile::Open(name_and_label_pair2->At(0)->GetName());
  TTree *tree2 = (TTree*)file2->Get("t");

  TFile *outFile = new TFile("MatrixOfMatches.root", "RECREATE");

  if (file1->IsZombie() || file2->IsZombie() || outFile->IsZombie()) return;

  TGraphErrors *outGraph = new TGraphErrors();
	
  vector<MyEvent> list1 = GetListToSort(tree1);
  vector<MyEvent> list2 = GetListToSort(tree2);

  outGraph = new TGraphErrors();

  Int_t numFoundMatches = 0;

  for (vector<MyEvent>::const_iterator i1begin=list1.begin(), i1=i1begin, i1end=list1.end(); i1!=i1end; ++i1) {

    Int_t val1 	 = i1->entry;
    Int_t run1	 = i1->run;
    Int_t lumi1  = i1->lumi;
    Int_t event1 = i1->event;

    for (vector<MyEvent>::const_iterator i2begin=list2.begin(), i2=i2begin, i2end=list2.end(); i2!=i2end; ++i2) {

      Int_t val2   = i2->entry;
      Int_t run2   = i2->run;
      Int_t lumi2  = i2->lumi;
      Int_t event2 = i2->event;
			
      if (run1 == run2 && lumi1==lumi2 && event1 == event2) {

	// cout<<"match! run1:"<<run1<<" run2:"<<run2<<" event1:"<<event1<<" event2:"<<event2<<" lumi1:"<<lumi1<<" lumi2:"<<lumi2<<"entry1: "<<val1<<"entry2: "<<val2<<endl;

	Int_t ipop = i2 - list2.begin();
	
	outGraph->SetPoint(numFoundMatches, val1, val2);
	outGraph->SetPointError(numFoundMatches, run1, event1);
	
	list2.erase(list2.begin()+ipop);				

	numFoundMatches++;
	
	break;
      }
    }
  }

  cout << " Number of pairs: " << outGraph->GetN() << endl;
	
  outGraph->Write("matchList");

  // save the filenames of the two input files in the output root file to prevent mismatches in 
  // the next scripts

  TObjString *tos_file1  = new TObjString(name_and_label_pair1->At(0)->GetName());
  TObjString *tos_label1 = new TObjString(name_and_label_pair1->At(1)->GetName());
  TObjString *tos_file2  = new TObjString(name_and_label_pair2->At(0)->GetName());
  TObjString *tos_label2 = new TObjString(name_and_label_pair2->At(1)->GetName());
  tos_file1->Write("rootfile1");
  tos_label1->Write("label1");
  tos_file2->Write("rootfile2");
  tos_label2->Write("label2");

  outFile->Write();

  return;
}


//------------------------------------------------------------------------------
// Get list that contains the run and event number for an entry
//------------------------------------------------------------------------------
vector<MyEvent> GetListToSort(TTree* tree)
{
  Int_t numEntries = tree->GetEntries(); 

  vector<MyEvent> retVec;

  retVec.clear();

  for (Int_t jEntry=0; jEntry<numEntries; ++jEntry) {

    tree->GetEntry(jEntry);

    Int_t run   = (Int_t)tree->GetBranch("runNumber")->GetLeaf("runNumber")->GetValue();
    Int_t lumi  = (Int_t)tree->GetBranch("lumiBlockNumber")->GetLeaf("lumiBlockNumber")->GetValue();
    Int_t event = (Int_t)tree->GetBranch("eventNumber")->GetLeaf("eventNumber")->GetValue();
    
    retVec.push_back(MyEvent(jEntry, run, lumi, event));
  }

  RemoveDuplicates(retVec);

  return retVec;
}


//------------------------------------------------------------------------------
// Check for duplicates - this is NOT optimized
//------------------------------------------------------------------------------
vector<MyEvent> RemoveDuplicates(vector<MyEvent> const &vec)
{
  vector<MyEvent> retVec;

  retVec.clear();

  for (vector<MyEvent>::const_iterator it1=vec.begin(), vecEnd=vec.end(); it1!=vecEnd; ++it1) {

    Bool_t hasDouble = false;
    
    for (vector<MyEvent>::const_iterator it2=vec.begin(); it2!=vecEnd; ++it2) {
				
      Bool_t isDouble = true;

      isDouble &= (it1->entry != it2->entry);
      isDouble &= (it1->run   == it2->run);
      isDouble &= (it1->lumi  == it2->lumi);
      isDouble &= (it1->event == it2->event);

      hasDouble = isDouble;
      
      if (hasDouble) break;
    }
		
    if (hasDouble)
      continue;
    else
      retVec.push_back(*it1);
  }

  return retVec;
}


