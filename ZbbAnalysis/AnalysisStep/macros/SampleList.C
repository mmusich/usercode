#ifndef SAMPLELIST_H
#define SAMPLELIST_H

#include <iostream>

#include <TList.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TXMLAttr.h>
#include <TString.h>

#include "Sample.C"

using namespace std;

class SampleList {
 public:
  SampleList() { listOfSample=new TList();}
  Int_t ParseFile(TString); 
  void ParseSampleList (TXMLNode *node);
  Sample* ParseSample(TXMLNode *node, Bool_t isMC);
  TList* listOfSample;
  Int_t GetSize() const {return listOfSample->GetSize();}
};

Int_t SampleList::ParseFile(TString filename) {
  TDOMParser *domParser = new TDOMParser();
  domParser->SetValidate(true);  // validate with DTD
  Int_t parsecode = domParser->ParseFile(filename); // get XML document

  if (  parsecode < 0 ) { 
    cerr << domParser->GetParseCodeMessage(parsecode); 
    return -1; 
  }
  TXMLNode *node = domParser->GetXMLDocument()->GetRootNode();
  ParseSampleList(node);

  return 0; 
}


void SampleList::ParseSampleList (TXMLNode *node) {

  for (; node; node = node->GetNextNode()) {

    if (node->GetNodeType() == TXMLNode::kXMLElementNode) { // Element Node
      cout << node->GetNodeName() << ": ";

      if (strcmp(node->GetNodeName(), "Sample") == 0) {
	Bool_t isMC=false;
	if (node->HasAttributes()) {
	  TList *attrList = node->GetAttributes();
	  TXMLAttr *attr = 0;
	  TIter next(attrList);
	  while ((attr=(TXMLAttr*)next())) {  // loop on attributes
	    cout << attr->GetName() << ":" << attr->GetValue();
	    if (strcmp(attr->GetName(), "ID") == 0) 	      
	      cout<<"Reading sample: "<< atoi(attr->GetValue()) <<endl;
	    if (strcmp(attr->GetName(), "IsMC") == 0 )  
	      isMC = ( strcmp(attr->GetValue(), "true") == 0 ) ;	    	 
	  } // end loop on attributes
	}
	listOfSample->Add(ParseSample(node->GetChildren(), isMC));
      }
    } // end dealing with Element node

    if (node->GetNodeType() == TXMLNode::kXMLTextNode)  // Text node
      cout << node->GetContent();
    ParseSampleList(node->GetChildren());
  } // end loop on Element nodes
}

Sample* SampleList::ParseSample(TXMLNode *node, Bool_t isMC){  
  TString inputrootfile;
  TString histolabel("");
  Int_t histocolor(0);
  Float_t lumi(0.);
  Float_t xsection(0.);

  for ( ; node; node = node->GetNextNode()) {
    if (node->GetNodeType() == TXMLNode::kXMLElementNode) { // Element Node
      if (strcmp(node->GetNodeName(), "InputRootFile") == 0)
	inputrootfile=node->GetText();	
      if (strcmp(node->GetNodeName(), "Label") == 0)
	histolabel=node->GetText();
      if (strcmp(node->GetNodeName(), "Color") == 0)
	histocolor=atoi(node->GetText());
      if (strcmp(node->GetNodeName(), "Lumi") == 0)
	lumi=atof(node->GetText());
      if (strcmp(node->GetNodeName(), "XSection") == 0)
	xsection=atof(node->GetText());
    }
  }

  return new Sample(isMC, inputrootfile, histolabel, histocolor, lumi, xsection);
}

#endif
