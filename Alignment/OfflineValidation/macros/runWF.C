#include "TKey.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TText.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>

TString prepareString(const char *dirname=".", const char *ext=".root",TList *LabelList);
void list_files(const char *dirname=".", const char *ext=".root");

void runWF(TString namesandlabels,Int_t nOfFiles){

  // **********************************************
  // retrieve the input file names and the labels
  // **********************************************

  TList *FileList = new TList();
  TList *LabelList = new TList();
  
  TObjArray *nameandlabelpairs = namesandlabels.Tokenize(",");
  for (Int_t i = 0; i < nameandlabelpairs->GetEntries(); ++i) {
    TObjArray *aFileLegPair = TString(nameandlabelpairs->At(i)->GetName()).Tokenize("=");
    
    if(aFileLegPair->GetEntries() == 2) {
      FileList->Add(aFileLegPair->At(0));  
      LabelList->Add(aFileLegPair->At(1));
    }
    else {
      std::cout << "Please give file name and legend entry in the following form:\n" 
		<< " filename1=legendentry1,filename2=legendentry2\n";   
    }    
  }
  
  const Int_t NOfFiles = FileList->GetSize();
  
  if(NOfFiles!=nOfFiles){
    cout<<"&MSG-e: NOfFiles = "<<nOfFiles<<endl;
    return;
  }  
  
   
  TString files[NOfFiles];
  TString legs[NOfFiles]; 

  for(Int_t j=0; j < NOfFiles; j++) {   
    
    TObjString* file = (TObjString*)FileList->At(j);
    files[j] = file->String();

    TObjString* legend = (TObjString*)LabelList->At(j);
    legs[j] = legend->String();

    cout<<"LegLabels["<<j<<"]"<<legs[j]<<endl;
  }

 // **********************************************
  // load the analysis class
  // **********************************************

  gROOT->LoadMacro("./OfflineValidationTreeAnalysis.C++g");
  
  for(Int_t j=0; j < NOfFiles; j++) {   
    OfflineValidationTreeAnalysis* a = new OfflineValidationTreeAnalysis(files[j],legs[j]);
    a->Loop(); 
    delete a;
   
    gSystem->Sleep(100);
    
    // check if the output folder exists and if not creates
    const char* append = (legs[j].ReplaceAll(" ","_")).Data();
    TSystemDirectory directory(Form("plots_%s",append),Form("plots_%s",append));
    TList *l = directory.GetListOfFiles();
    if (l) {
      std::cout<< "&MSG-i: runWF(): >>>>>> directory plots_"<<append<<" exists already: deleting it before re-creating!" << std::endl;
      TString processline = Form(".! rm -fr plots_%s",append);
      gROOT->ProcessLine(processline.Data());
      processline.Clear();
      TString processline = Form(".! mkdir plots_%s",append);
      gROOT->ProcessLine(processline.Data());
      processline.Clear();
    } else {
      std::cout<< "&MSG-i: output directory does not exists. Creating it!" << std::endl;
      TString processline = Form(".! mkdir plots_%s",append);
      gROOT->ProcessLine(processline.Data());
      processline.Clear();
    }
  }
  
  // mv the output pdf files in the output folder
  for(Int_t j=0; j < NOfFiles; j++) {   
    const char* append = (legs[j].ReplaceAll(" ","_")).Data();
    TString processline2 = Form(".! mv cProf_*%s* cPull_*%s* plots_%s",append,append,append);
    gROOT->ProcessLine(processline2.Data());
    append = "";
  }    

  //list_files();

  gSystem->Sleep(100);
  
  // **********************************************
  // load the superimposition script
  // **********************************************

  gROOT->LoadMacro("./generalizedSuperImpose.C++g");

  // prepare the input string
  TString inputString=prepareString(".",".root",LabelList);

  // run the script
  generalizedSuperImpose(inputString,nOfFiles);

  // create the output folder name
  TString outCompDir = "Comparison";
  for(Int_t j=0; j < NOfFiles; j++) {
    const char* append = (legs[j].ReplaceAll(" ","_")).Data();
    //std::cout<<append<<std::endl;
    TString toAppend = Form("_%s",append);
    //std::cout<<toAppend<<std::endl;
    outCompDir.Append(toAppend.Data());
  }

  // check if the output folder exists and if not creates
  TSystemDirectory outdir(outCompDir,outCompDir);
  TList *m = outdir.GetListOfFiles();
  if (m) {
    std::cout<< "runWF(): >>>>>> directory "<<outCompDir<<" exists already: deleting it!" << std::endl;
    TString processline = Form(".! rm -fr %s",outCompDir.Data());
    gROOT->ProcessLine(processline.Data());
    processline.Clear();
    TString processline = Form(".! mkdir %s",outCompDir.Data());
    gROOT->ProcessLine(processline.Data());
    processline.Clear();
    } else {
    std::cout<<"runWF(): >>>>>> directory "<<outCompDir<<" does not exit yet: creating it!" << std::endl;
    TString processline = Form(".! mkdir %s",outCompDir.Data());
    gROOT->ProcessLine(processline.Data());
    processline.Clear();
  }

  // mv the final comparison output pdf files in the output folder  
  gROOT->ProcessLine((".! mv *.pdf ./"+outCompDir).Data());

}
/*--------------------------------------------------------------------*/
void list_files(const char *dirname, const char *ext)
/*--------------------------------------------------------------------*/
{
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
	cout << fname.Data() << endl;
      }
    }
  }
}

/*--------------------------------------------------------------------*/
TString prepareString(const char *dirname, const char *ext,TList *LabelList)
/*--------------------------------------------------------------------*/
{

  TString FilesList[10];
  TString LabelsList[10];

  TString outputString="";
  Int_t nfiles=0;

  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext) && fname.BeginsWith("TrackerOfflineValidationStandalone")) {
	bool isFileOk_=false;
	for(Int_t j=0; j < LabelList->GetEntries(); j++) {   
	  TObjString* legend = (TObjString*)LabelList->At(j);
	  TString leg = legend->String();
	  //std::cout<<"fname: "<<fname<<"| leg: "<<leg.ReplaceAll(" ","_")<<std::endl;
	  if(fname.Contains(leg.ReplaceAll(" ","_"))){
	    isFileOk_=true;
	  }
	}

	if(!isFileOk_) continue;

	cout <<"nfiles: "<< nfiles << " adding file: "<< fname.Data() << endl;
	FilesList[nfiles]=fname;
	TObjArray *all = fname.Tokenize(".");
	for (Int_t i = 0; i < all->GetEntries(); ++i) {
	  TString chunk = all->At(i)->GetName();
	  if (chunk.Contains("TrackerOfflineValidationStandalone")){
	    chunk.ReplaceAll("TrackerOfflineValidationStandalone_","");
	    chunk.ReplaceAll("_"," ");
	    LabelsList[nfiles]=chunk;
	  } // if .root in chunk
	} // loop on chunks in file name
	nfiles++;
      } // if correct file
    }
  }

  for(Int_t i=0;i<nfiles;i++){
    TString filestring  = FilesList[i];
    TString labelstring = LabelsList[i];
    outputString.Append(filestring);
    outputString.Append("=");
    outputString.Append(labelstring);
    if (i!=nfiles-1){
      outputString.Append(",");
    }
  }
  
  //cout << "the Output String is: "<<outputString<<endl;

  return outputString;

}






