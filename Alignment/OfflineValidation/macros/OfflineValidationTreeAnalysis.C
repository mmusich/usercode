#define OfflineValidationTreeAnalysis_cxx
#include "OfflineValidationTreeAnalysis.h"
#include <TH2.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TEntryList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TString.h>
#include "TArrow.h"
#include <TCanvas.h>

void OfflineValidationTreeAnalysis::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L OfflineValidationTreeAnalysis.C
  //      Root > OfflineValidationTreeAnalysis t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;
  
  Double_t cmToUm = 10000;

  Long64_t nentries = fChain->GetEntriesFast();
  const char * c_fChainName = fChainName.Data();

  // Pixel
  TCanvas *cProfBPixRMSx_DMR_Eta = new TCanvas(Form("%s_cProfXBPixEta",c_fChainName),"BPIX-x median vs Eta",700,700);
  TCanvas *cProfFPixRMSx_DMR_Eta = new TCanvas(Form("%s_cProfXFPixEta",c_fChainName),"FPIX-x median vs Eta",700,700);
  TCanvas *cProfBPixRMSy_DMR_Eta = new TCanvas(Form("%s_cProfYBPixEta",c_fChainName),"BPIX-y median vs Eta",700,700);
  TCanvas *cProfFPixRMSy_DMR_Eta = new TCanvas(Form("%s_cProfYFPixEta",c_fChainName),"FPIX-y median vs Eta",700,700);
										                         
  TCanvas *cProfBPixRMSx_DMR_Phi = new TCanvas(Form("%s_cProfXBPixPhi",c_fChainName),"BPIX-x median vs Phi",700,700);
  TCanvas *cProfFPixRMSx_DMR_Phi = new TCanvas(Form("%s_cProfXFPixPhi",c_fChainName),"FPIX-x median vs Phi",700,700);
  TCanvas *cProfBPixRMSy_DMR_Phi = new TCanvas(Form("%s_cProfYBPixPhi",c_fChainName),"BPIX-y median vs Phi",700,700);
  TCanvas *cProfFPixRMSy_DMR_Phi = new TCanvas(Form("%s_cProfYFPixPhi",c_fChainName),"FPIX-y median vs Phi",700,700);

  TCanvas *cProfBPixRMSx_DMR_Z = new TCanvas(Form("%s_cProfXBPixZ",c_fChainName),"BPIX-x median vs Z",700,700);
  TCanvas *cProfFPixRMSx_DMR_Z = new TCanvas(Form("%s_cProfXFPixZ",c_fChainName),"FPIX-x median vs Z",700,700);
  TCanvas *cProfBPixRMSy_DMR_Z = new TCanvas(Form("%s_cProfYBPixZ",c_fChainName),"BPIX-y median vs Z",700,700);
  TCanvas *cProfFPixRMSy_DMR_Z = new TCanvas(Form("%s_cProfYFPixZ",c_fChainName),"FPIX-y median vs Z",700,700);
									                           
  TCanvas *cPullBPixRMSx_Eta = new TCanvas(Form("%s_cPullXBPixEta",c_fChainName),"BPIX-x pull vs Eta",700,700);
  TCanvas *cPullFPixRMSx_Eta = new TCanvas(Form("%s_cPullXFPixEta",c_fChainName),"FPIX-x pull vs Eta",700,700);
  TCanvas *cPullBPixRMSy_Eta = new TCanvas(Form("%s_cPullYBPixEta",c_fChainName),"BPIX-y pull vs Eta",700,700);
  TCanvas *cPullFPixRMSy_Eta = new TCanvas(Form("%s_cPullYFPixEta",c_fChainName),"FPIX-y pull vs Eta",700,700);
									                           
  TCanvas *cPullBPixRMSx_Phi = new TCanvas(Form("%s_cPullXBPixPhi",c_fChainName),"BPIX-x pull vs Phi",700,700);
  TCanvas *cPullFPixRMSx_Phi = new TCanvas(Form("%s_cPullXFPixPhi",c_fChainName),"FPIX-x pull vs Phi",700,700);
  TCanvas *cPullBPixRMSy_Phi = new TCanvas(Form("%s_cPullYBPixPhi",c_fChainName),"BPIX-y pull vs Phi",700,700);
  TCanvas *cPullFPixRMSy_Phi = new TCanvas(Form("%s_cPullYFPixPhi",c_fChainName),"FPIX-y pull vs Phi",700,700);

  TCanvas *cPullBPixRMSx_Z = new TCanvas(Form("%s_cPullXBPixZ",c_fChainName),"BPIX-x pull vs Z",700,700);
  TCanvas *cPullFPixRMSx_Z = new TCanvas(Form("%s_cPullXFPixZ",c_fChainName),"FPIX-x pull vs Z",700,700);
  TCanvas *cPullBPixRMSy_Z = new TCanvas(Form("%s_cPullYBPixZ",c_fChainName),"BPIX-y pull vs Z",700,700);
  TCanvas *cPullFPixRMSy_Z = new TCanvas(Form("%s_cPullYFPixZ",c_fChainName),"FPIX-y pull vs Z",700,700);

  // Strips
  TCanvas *cProfTIBRMS_DMR_Eta   = new TCanvas(Form("%s_cProfTIBEta",c_fChainName),"TIB median vs Eta",700,700);
  TCanvas *cProfTIDRMS_DMR_Eta   = new TCanvas(Form("%s_cProfTIDEta",c_fChainName),"TID median vs Eta",700,700);
  TCanvas *cProfTOBRMS_DMR_Eta   = new TCanvas(Form("%s_cProfTOBEta",c_fChainName),"TOB median vs Eta",700,700);
  TCanvas *cProfTECRMS_DMR_Eta   = new TCanvas(Form("%s_cProfTECEta",c_fChainName),"TEC median vs Eta",700,700);
										                    
  TCanvas *cProfTIBRMS_DMR_Phi   = new TCanvas(Form("%s_cProfTIBPhi",c_fChainName),"TIB median vs Phi",700,700);
  TCanvas *cProfTIDRMS_DMR_Phi   = new TCanvas(Form("%s_cProfTIDPhi",c_fChainName),"TID median vs Phi",700,700);
  TCanvas *cProfTOBRMS_DMR_Phi   = new TCanvas(Form("%s_cProfTOBPhi",c_fChainName),"TOB median vs Phi",700,700);
  TCanvas *cProfTECRMS_DMR_Phi   = new TCanvas(Form("%s_cProfTECPhi",c_fChainName),"TEC median vs Phi",700,700);

  TCanvas *cProfTIBRMS_DMR_Z   = new TCanvas(Form("%s_cProfTIBZ",c_fChainName),"TIB median vs Z",700,700);
  TCanvas *cProfTIDRMS_DMR_Z   = new TCanvas(Form("%s_cProfTIDZ",c_fChainName),"TID median vs Z",700,700);
  TCanvas *cProfTOBRMS_DMR_Z   = new TCanvas(Form("%s_cProfTOBZ",c_fChainName),"TOB median vs Z",700,700);
  TCanvas *cProfTECRMS_DMR_Z   = new TCanvas(Form("%s_cProfTECZ",c_fChainName),"TEC median vs Z",700,700);
									                      
  TCanvas *cPullTIBRMS_Eta   = new TCanvas(Form("%s_cPullTIBEta",c_fChainName),"TIB pull vs Eta",700,700);
  TCanvas *cPullTIDRMS_Eta   = new TCanvas(Form("%s_cPullTIDEta",c_fChainName),"TID pull vs Eta",700,700);
  TCanvas *cPullTOBRMS_Eta   = new TCanvas(Form("%s_cPullTOBEta",c_fChainName),"TOB pull vs Eta",700,700);
  TCanvas *cPullTECRMS_Eta   = new TCanvas(Form("%s_cPullTECEta",c_fChainName),"TEC pull vs Eta",700,700);
									                      
  TCanvas *cPullTIBRMS_Phi   = new TCanvas(Form("%s_cPullTIBPhi",c_fChainName),"TIB pull vs Phi",700,700);
  TCanvas *cPullTIDRMS_Phi   = new TCanvas(Form("%s_cPullTIDPhi",c_fChainName),"TID pull vs Phi",700,700);
  TCanvas *cPullTOBRMS_Phi   = new TCanvas(Form("%s_cPullTOBPhi",c_fChainName),"TOB pull vs Phi",700,700);
  TCanvas *cPullTECRMS_Phi   = new TCanvas(Form("%s_cPullTECPhi",c_fChainName),"TEC pull vs Phi",700,700);

  TCanvas *cPullTIBRMS_Z   = new TCanvas(Form("%s_cPullTIBZ",c_fChainName),"TIB pull vs Z",700,700);
  TCanvas *cPullTIDRMS_Z   = new TCanvas(Form("%s_cPullTIDZ",c_fChainName),"TID pull vs Z",700,700);
  TCanvas *cPullTOBRMS_Z   = new TCanvas(Form("%s_cPullTOBZ",c_fChainName),"TOB pull vs Z",700,700);
  TCanvas *cPullTECRMS_Z   = new TCanvas(Form("%s_cPullTECZ",c_fChainName),"TEC pull vs Z",700,700);

  TCanvas* DMRCanvEta[8]  = {cProfBPixRMSx_DMR_Eta,
			     cProfBPixRMSy_DMR_Eta,
			     cProfFPixRMSx_DMR_Eta,
			     cProfFPixRMSy_DMR_Eta,
			     cProfTIBRMS_DMR_Eta,
			     cProfTIDRMS_DMR_Eta,
			     cProfTOBRMS_DMR_Eta,
			     cProfTECRMS_DMR_Eta};
  
  TCanvas* DMRCanvPhi[8]  = {cProfBPixRMSx_DMR_Phi,
			     cProfBPixRMSy_DMR_Phi,
			     cProfFPixRMSx_DMR_Phi,
			     cProfFPixRMSy_DMR_Phi,
			     cProfTIBRMS_DMR_Phi,
			     cProfTIDRMS_DMR_Phi,
			     cProfTOBRMS_DMR_Phi,
			     cProfTECRMS_DMR_Phi};

  TCanvas* DMRCanvZ[8]    = {cProfBPixRMSx_DMR_Z,
			     cProfBPixRMSy_DMR_Z,
			     cProfFPixRMSx_DMR_Z,
			     cProfFPixRMSy_DMR_Z,
			     cProfTIBRMS_DMR_Z,
			     cProfTIDRMS_DMR_Z,
			     cProfTOBRMS_DMR_Z,
			     cProfTECRMS_DMR_Z};

  TCanvas* PullCanvEta[8] = {cPullBPixRMSx_Eta,
			     cPullBPixRMSy_Eta,
			     cPullFPixRMSx_Eta,
			     cPullFPixRMSy_Eta,
			     cPullTIBRMS_Eta,
			     cPullTIDRMS_Eta,
			     cPullTOBRMS_Eta,
			     cPullTECRMS_Eta};

  TCanvas* PullCanvPhi[8] = {cPullBPixRMSx_Phi,
			     cPullBPixRMSy_Phi,
			     cPullFPixRMSx_Phi,
			     cPullFPixRMSy_Phi,
			     cPullTIBRMS_Phi,
			     cPullTIDRMS_Phi,
			     cPullTOBRMS_Phi,
			     cPullTECRMS_Phi};

  TCanvas* PullCanvZ[8] = {cPullBPixRMSx_Z,
			   cPullBPixRMSy_Z,
			   cPullFPixRMSx_Z,
			   cPullFPixRMSy_Z,
			   cPullTIBRMS_Z,
			   cPullTIDRMS_Z,
			   cPullTOBRMS_Z,
			   cPullTECRMS_Z};
  			  									  	      			 									      
  // get the maxima/minima
  //std::cout<<"overall min(medianX) == "<<fChain->GetMinimum("medianX")*cmToUm<<std::endl;
  //std::cout<<"overall max(medianX) == "<<fChain->GetMaximum("medianX")*cmToUm<<std::endl;
  //std::cout<<"overall min(medianY) == "<<fChain->GetMinimum("medianY")*cmToUm<<std::endl;
  //std::cout<<"overall max(medianY) == "<<fChain->GetMaximum("medianY")*cmToUm<<std::endl;

  TString dets[8]  = {"BPix-x","BPix-y","FPix-x","FPix-y","TIB","TID","TOB","TEC"};
  TString coord[8] = {"x'","y'","x'","y'","x'","x'","x'","x'"};
  Int_t theSubDetId[8] = {1,1,2,2,3,4,5,6};
  TString theCuts[8];
  Double_t theMaxima[8];
  Double_t theMinima[8];
  std::map<TString,std::pair<Double_t,Double_t > > mymap; 
  std::map<TString,Double_t> nbins; 

  for (UInt_t i=0;i<8;i++){
    theCuts[i]=Form("(subDetId==%i)&&(entries>30)",theSubDetId[i]);
    //std::cout<<"*********************************"<<std::endl;
    //std::cout<<"i="<<i<<", theCuts["<<i<<"] is:"<<theCuts[i]<<std::endl;
    fChain->Draw(">>eList",theCuts[i], "entrylist"); //Note: "entrylist" is a special option!
    TEntryList * const eList = (TEntryList*)gDirectory->Get("eList");
    if (eList) {
      fChain->SetEntryList(eList);
      //std::cout<<"After cut applied ("<<dets[i]<<") :\n";
      if((i==1)||(i==3)){
	//std::cout<<"min(y) == "<<fChain->GetMinimum("medianY")*cmToUm<<std::endl;
	//std::cout<<"max(y) == "<<fChain->GetMaximum("medianY")*cmToUm<<std::endl;
	theMinima[i]=fChain->GetMinimum("medianY")*cmToUm;
	theMaxima[i]=fChain->GetMaximum("medianY")*cmToUm;

	Double_t theExtreme(0.);
	if (TMath::Abs(theMinima[i])>theMaxima[i]){
	  theExtreme=TMath::Abs(theMinima[i]);
	} else {
	  theExtreme=TMath::Abs(theMaxima[i]);
	}

	std::pair <Double_t,Double_t> foo;
	foo = std::make_pair(-theExtreme*1.10,theExtreme*1.10);
	mymap[dets[i]] = foo;
	nbins[dets[i]] = (foo.second - foo.first)/0.2;

      } else {
	//std::cout<<"min(x) == "<<fChain->GetMinimum("medianX")*cmToUm<<std::endl;
	//std::cout<<"max(x) == "<<fChain->GetMaximum("medianX")*cmToUm<<std::endl;
	theMinima[i]=fChain->GetMinimum("medianX")*cmToUm;
	theMaxima[i]=fChain->GetMaximum("medianX")*cmToUm;
	
	Double_t theExtreme(0.);
	if (TMath::Abs(theMinima[i])>theMaxima[i]){
	  theExtreme=TMath::Abs(theMinima[i]);
	} else {
	  theExtreme=TMath::Abs(theMaxima[i]);
	}

	std::pair <Double_t,Double_t> foo;
	foo = std::make_pair(-theExtreme*1.10,theExtreme*1.10);
	mymap[dets[i]] = foo;
	nbins[dets[i]] = (foo.second - foo.first)/0.2;
      }
    }

    fChain->SetEntryList(0);
    delete eList;
  }


  for (UInt_t i=0;i<8;i++){
    // std::cout<<"mymap["<<dets[i]<<"]: "<<mymap[dets[i]].first<<","<<mymap[dets[i]].second<<std::endl;
  }

  /*
  fChain->Draw(">>eList", "(subDetId==1)&&(entries>30)", "entrylist"); //Note: "entrylist" is a special option!
  TEntryList * const eList = (TEntryList*)gDirectory->Get("eList");
  if (eList) {
    fChain->SetEntryList(eList);
    std::cout<<"After cut applied:\n";
    std::cout<<"min == "<<fChain->GetMinimum("medianX")<<std::endl;
    std::cout<<"max == "<<fChain->GetMaximum("medianX")<<std::endl;
  }
  */

  // rms of DMR plot maps

  std::map<TString,TProfile*> profMedianEta;
  std::map<TString,TH1F*> h1MedianEta;
  std::map<TString,TH2F*> h2MedianEta;

  std::map<TString,TProfile*> profMedianPhi;
  std::map<TString,TH1F*> h1MedianPhi;
  std::map<TString,TH2F*> h2MedianPhi;

  std::map<TString,TProfile*> profMedianZ;
  std::map<TString,TH1F*> h1MedianZ;
  std::map<TString,TH2F*> h2MedianZ;
  
  // mean of pull plot maps

  std::map<TString,TProfile*> profPullEta;
  std::map<TString,TH2F*> h2PullEta;

  std::map<TString,TProfile*> profPullPhi;
  std::map<TString,TH2F*> h2PullPhi;

  std::map<TString,TProfile*> profPullZ;
  std::map<TString,TH2F*> h2PullZ;

  for (UInt_t i=0;i<8;i++){
    
    // profMedian

    profMedianEta[dets[i]] = new TProfile(Form("prof_%s_VsEta",dets[i].Data()),
					  "profile x-residual vs |#eta|;detector pseudorapidity |#eta|;median of residual #mu_{1/2}(res_{x'})",
					  12,0.,2.8,mymap[dets[i]].first,mymap[dets[i]].second,"s");
    
    profMedianPhi[dets[i]] = new TProfile(Form("prof_%s_VsPhi",dets[i].Data()),
					  "profile x-residual vs #varphi;detector azimuth #varphi;median of residual #mu_{1/2}(res_{x'})",
					  12,-TMath::Pi(),TMath::Pi(),mymap[dets[i]].first,mymap[dets[i]].second,"s");

    profMedianZ[dets[i]]   = new TProfile(Form("prof_%s_VsZ",dets[i].Data()),
					  "profile x-residual vs z;detector z [cm];median of residual #mu_{1/2}(res_{x'})",
					  12,-280,280,mymap[dets[i]].first,mymap[dets[i]].second,"s");
    // profPull

    profPullEta[dets[i]]   = new TProfile(Form("profNormRes_%s_VsEta",dets[i].Data()),
					  Form("Normalized Residuals RMS vs |#eta| (%s);detector pseudorapidity |#eta|;r.m.s. of pull",dets[i].Data()),
					  12,0.,2.8,0.7,1.2,"i");
    
    profPullPhi[dets[i]]   = new TProfile(Form("profNormRes_%s_VsPhi",dets[i].Data()),
					  Form("Normalized Residuals RMS vs #varphi (%s);detector azimuth #phi;r.m.s. of pull",dets[i].Data()),
					  12,-TMath::Pi(),TMath::Pi(),0.7,1.2,"i");

    profPullZ[dets[i]]     = new TProfile(Form("profNormRes_%s_VsZ",dets[i].Data()),
					  Form("Normalized Residuals RMS vs z (%s);detector z [cm];r.m.s. of pull",dets[i].Data()),
					  12,-280.,280.,0.7,1.2,"i");
    
    // h2Median

    h2MedianEta[dets[i]]   = new TH2F(Form("h2_%s_VsEta",dets[i].Data()),
				      Form("profile #mu_{1/2}(res_{%s}) vs |#eta| (%s);detector pseudorapidity |#eta|;median of residual #mu_{1/2}(res_{x'})",coord[i].Data(),dets[i].Data()),
				      12,0.,2.8,nbins[dets[i]],mymap[dets[i]].first,mymap[dets[i]].second);

    h2MedianPhi[dets[i]]   = new TH2F(Form("h2_%s_VsPhi",dets[i].Data()),
				      Form("profile #mu_{1/2}(res_{%s}) vs #varphi (%s);detector azimuth #varphi;median of residual #mu_{1/2}(res_{x'})",coord[i].Data(),dets[i].Data()),
				      12,-TMath::Pi(),TMath::Pi(),nbins[dets[i]],mymap[dets[i]].first,mymap[dets[i]].second);

    h2MedianZ[dets[i]]     = new TH2F(Form("h2_%s_VsZ",dets[i].Data()),
				      Form("profile #mu_{1/2}(res_{%s}) vs z (%s);detector z [cm];median of residual #mu_{1/2}(res_{x'})",coord[i].Data(),dets[i].Data()),
				      12,-280.,280.,nbins[dets[i]],mymap[dets[i]].first,mymap[dets[i]].second);

    // h2Pull

    h2PullEta[dets[i]]     = new TH2F(Form("h2NormRes_%s_VsEta",dets[i].Data()),
				      Form("Normalized Residuals RMS vs |#eta| (%s);detector pseudorapidity |#eta|;r.m.s. of pull",dets[i].Data()),
				      12,0.,2.8,100,0.7,1.2);

    h2PullPhi[dets[i]]     = new TH2F(Form("h2NormRes_%s_VsPhi",dets[i].Data()),
				      Form("Normalized Residuals RMS vs #varphi (%s);detector azimuth #varphi;r.m.s. of pull",dets[i].Data()),
				      12,-TMath::Pi(),TMath::Pi(),100,0.7,1.2);  
    
    h2PullZ[dets[i]]       = new TH2F(Form("h2NormRes_%s_VsZ",dets[i].Data()),
				      Form("Normalized Residuals RMS vs z (%s);detector z [cm];r.m.s. of pull",dets[i].Data()),
				      12,-280.,280.,100,0.7,1.2);  

    // h1Median

    h1MedianEta[dets[i]]   = new TH1F(Form("h1_%sRMSEta",dets[i].Data()),
				      Form("RMS of DMR vs |#eta| (%s);detector pseudorapidity |#eta|;r.m.s. of DMR [#mum]",dets[i].Data()),
				      12,0.,2.8);

    h1MedianPhi[dets[i]]   = new TH1F(Form("h1_%sRMSPhi",dets[i].Data()),
				       Form("RMS of DMR vs #phi (%s);detector azimuth #varphi;r.m.s. of DMR [#mum]",dets[i].Data()),
				      12,-TMath::Pi(),TMath::Pi()); 

    h1MedianZ[dets[i]]     = new TH1F(Form("h1_%sRMSZ",dets[i].Data()),
				      Form("RMS of DMR vs z (%s);detector z [cm];r.m.s. of DMR [#mum]",dets[i].Data()),
				      12,-280.,280.); 
    
  }

  // std::cout<<"before event loop"<<std::endl;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // std::cout<<"ientry:"<<ientry<<std::endl;
      
      if(entries<30) continue;

      if(subDetId==1){

	// ***** vs x' *****

	h2MedianEta["BPix-x"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	profMedianEta["BPix-x"]->Fill(TMath::Abs(posEta),medianX*cmToUm);

	h2MedianPhi["BPix-x"]->Fill(posPhi,medianX*cmToUm);
	profMedianPhi["BPix-x"]->Fill(posPhi,medianX*cmToUm);
	
	h2MedianZ["BPix-x"]->Fill(posZ,medianX*cmToUm);
	profMedianZ["BPix-x"]->Fill(posZ,medianX*cmToUm);

	h2PullEta["BPix-x"]->Fill(TMath::Abs(posEta),rmsNormX);
	profPullEta["BPix-x"]->Fill(TMath::Abs(posEta),rmsNormX);
	
	h2PullPhi["BPix-x"]->Fill(posPhi,rmsNormX);
	profPullPhi["BPix-x"]->Fill(posPhi,rmsNormX);

	h2PullZ["BPix-x"]->Fill(posZ,rmsNormX);
	profPullZ["BPix-x"]->Fill(posZ,rmsNormX);

	// ***** vs y' *****

	h2MedianEta["BPix-y"]->Fill(TMath::Abs(posEta),medianY*cmToUm);
	profMedianEta["BPix-y"]->Fill(TMath::Abs(posEta),medianY*cmToUm);
	
	h2MedianPhi["BPix-y"]->Fill(posPhi,medianY*cmToUm);
	profMedianPhi["BPix-y"]->Fill(posPhi,medianY*cmToUm);

	h2MedianZ["BPix-y"]->Fill(posZ,medianY*cmToUm);
	profMedianZ["BPix-y"]->Fill(posZ,medianY*cmToUm);

	h2PullEta["BPix-y"]->Fill(TMath::Abs(posEta),rmsNormY);
	profPullEta["BPix-y"]->Fill(TMath::Abs(posEta),rmsNormY);

	h2PullPhi["BPix-y"]->Fill(posPhi,rmsNormY);
	profPullPhi["BPix-y"]->Fill(posPhi,rmsNormY);
	
	h2PullZ["BPix-y"]->Fill(posZ,rmsNormY);
	profPullZ["BPix-y"]->Fill(posZ,rmsNormY);
	
      } else if(subDetId==2){
	
	// ***** vs x' *****

	h2MedianEta["FPix-x"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	profMedianEta["FPix-x"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	
	h2MedianPhi["FPix-x"]->Fill(posPhi,medianX*cmToUm);
	profMedianPhi["FPix-x"]->Fill(posPhi,medianX*cmToUm);

	h2MedianZ["FPix-x"]->Fill(posZ,medianX*cmToUm);
	profMedianZ["FPix-x"]->Fill(posZ,medianX*cmToUm);

	h2PullEta["FPix-x"]->Fill(TMath::Abs(posEta),rmsNormX);
	profPullEta["FPix-x"]->Fill(TMath::Abs(posEta),rmsNormX);
	
	h2PullPhi["FPix-x"]->Fill(posPhi,rmsNormX);
	profPullPhi["FPix-x"]->Fill(posPhi,rmsNormX);
	
	h2PullZ["FPix-x"]->Fill(posZ,rmsNormX);
	profPullZ["FPix-x"]->Fill(posZ,rmsNormX);

	// ***** vs y' *****

	h2MedianEta["FPix-y"]->Fill(TMath::Abs(posEta),medianY*cmToUm);
	profMedianEta["FPix-y"]->Fill(TMath::Abs(posEta),medianY*cmToUm);

	h2MedianPhi["FPix-y"]->Fill(posPhi,medianY*cmToUm);
	profMedianPhi["FPix-y"]->Fill(posPhi,medianY*cmToUm);
	
	h2MedianZ["FPix-y"]->Fill(posZ,medianY*cmToUm);
	profMedianZ["FPix-y"]->Fill(posZ,medianY*cmToUm);

	h2PullEta["FPix-y"]->Fill(TMath::Abs(posEta),rmsNormY);
	profPullEta["FPix-y"]->Fill(TMath::Abs(posEta),rmsNormY);

	h2PullPhi["FPix-y"]->Fill(posPhi,rmsNormY);
	profPullPhi["FPix-y"]->Fill(posPhi,rmsNormY);

	h2PullZ["FPix-y"]->Fill(posZ,rmsNormY);
	profPullZ["FPix-y"]->Fill(posZ,rmsNormY);

      }  else if(subDetId==3){             

	h2MedianEta["TIB"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	profMedianEta["TIB"]->Fill(TMath::Abs(posEta),medianX*cmToUm);

	h2MedianPhi["TIB"]->Fill(posPhi,medianX*cmToUm);
	profMedianPhi["TIB"]->Fill(posPhi,medianX*cmToUm);

	h2MedianZ["TIB"]->Fill(posZ,medianX*cmToUm);
	profMedianZ["TIB"]->Fill(posZ,medianX*cmToUm);

	h2PullEta["TIB"]->Fill(TMath::Abs(posEta),rmsNormX);
	profPullEta["TIB"]->Fill(TMath::Abs(posEta),rmsNormX);

	h2PullPhi["TIB"]->Fill(posPhi,rmsNormX);
	profPullPhi["TIB"]->Fill(posPhi,rmsNormX);

	h2PullZ["TIB"]->Fill(posZ,rmsNormX);
	profPullZ["TIB"]->Fill(posZ,rmsNormX);
	
      } else if(subDetId==4){      
  
	h2MedianEta["TID"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	profMedianEta["TID"]->Fill(TMath::Abs(posEta),medianX*cmToUm);

	h2MedianPhi["TID"]->Fill(posPhi,medianX*cmToUm);
	profMedianPhi["TID"]->Fill(posPhi,medianX*cmToUm);

	h2MedianZ["TID"]->Fill(posZ,medianX*cmToUm);
	profMedianZ["TID"]->Fill(posZ,medianX*cmToUm);

	h2PullEta["TID"]->Fill(TMath::Abs(posEta),rmsNormX);
	profPullEta["TID"]->Fill(TMath::Abs(posEta),rmsNormX);

	h2PullPhi["TID"]->Fill(posPhi,rmsNormX);
	profPullPhi["TID"]->Fill(posPhi,rmsNormX);

	h2PullZ["TID"]->Fill(posZ,rmsNormX);
	profPullZ["TID"]->Fill(posZ,rmsNormX);
	
      }  else if(subDetId==5){

	h2MedianEta["TOB"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	profMedianEta["TOB"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	
	h2MedianPhi["TOB"]->Fill(posPhi,medianX*cmToUm);
	profMedianPhi["TOB"]->Fill(posPhi,medianX*cmToUm);

	h2MedianZ["TOB"]->Fill(posZ,medianX*cmToUm);
	profMedianZ["TOB"]->Fill(posZ,medianX*cmToUm);

	h2PullEta["TOB"]->Fill(TMath::Abs(posEta),rmsNormX);
	profPullEta["TOB"]->Fill(TMath::Abs(posEta),rmsNormX);
	
	h2PullPhi["TOB"]->Fill(posPhi,rmsNormX);
	profPullPhi["TOB"]->Fill(posPhi,rmsNormX);

	h2PullZ["TOB"]->Fill(posZ,rmsNormX);
	profPullZ["TOB"]->Fill(posZ,rmsNormX);

      } else if (subDetId==6){
   
	h2MedianEta["TEC"]->Fill(TMath::Abs(posEta),medianX*cmToUm);
	profMedianEta["TEC"]->Fill(TMath::Abs(posEta),medianX*cmToUm);

	h2MedianPhi["TEC"]->Fill(posPhi,medianX*cmToUm);
	profMedianPhi["TEC"]->Fill(posPhi,medianX*cmToUm);

	h2MedianZ["TEC"]->Fill(posZ,medianX*cmToUm);
	profMedianZ["TEC"]->Fill(posZ,medianX*cmToUm);

	h2PullEta["TEC"]->Fill(TMath::Abs(posEta),rmsNormX);
	profPullEta["TEC"]->Fill(TMath::Abs(posEta),rmsNormX);
     
	h2PullPhi["TEC"]->Fill(posPhi,rmsNormX);
	profPullPhi["TEC"]->Fill(posPhi,rmsNormX);

	h2PullZ["TEC"]->Fill(posZ,rmsNormX);
	profPullZ["TEC"]->Fill(posZ,rmsNormX);

      }
  }
  
  setStyle("fire");
 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for (UInt_t i=0;i<8;i++){
   
    // %%% - vs eta - %%%

    makeNiceCanv(DMRCanvEta[i]);
    DMRCanvEta[i]->cd();
    gPad->SetRightMargin(0.12);
    makeNice2DPlotStyle(h2MedianEta[dets[i]]);
    h2MedianEta[dets[i]]->Draw("box");
    profMedianEta[dets[i]]->SetLineColor(kRed);
    profMedianEta[dets[i]]->SetFillColor(kRed);
    profMedianEta[dets[i]]->SetFillStyle(3002);
    profMedianEta[dets[i]]->SetMarkerColor(kRed);
    profMedianEta[dets[i]]->SetLineWidth(2.);
    profMedianEta[dets[i]]->SetMarkerStyle(20);
    profMedianEta[dets[i]]->SetMarkerSize(1.4);
    profMedianEta[dets[i]]->Draw("E2same");
    DMRCanvEta[i]->Draw();

    TArrow *l0=new TArrow(DMRCanvEta[i]->GetUxmin(),0.,DMRCanvEta[i]->GetUxmax(),0.,0.001);
    l0->SetLineColor(kBlack);
    l0->SetLineWidth(3);
    l0->SetLineStyle(9);

    l0->Draw("same");
    DMRCanvEta[i]->SaveAs(Form("cProf_%s_RMS_DMR_Eta_",dets[i].Data())+fChainName+".pdf");

    fOutputFile->cd();
    DMRCanvZ[i]->Write();
   
    // %%% - vs phi - %%%

    makeNiceCanv(DMRCanvPhi[i]);
    DMRCanvPhi[i]->cd();
    gPad->SetRightMargin(0.12);
    makeNice2DPlotStyle(h2MedianPhi[dets[i]]);
    h2MedianPhi[dets[i]]->Draw("box");
    profMedianPhi[dets[i]]->SetLineColor(kRed);
    profMedianPhi[dets[i]]->SetFillColor(kRed);
    profMedianPhi[dets[i]]->SetFillStyle(3002);
    profMedianPhi[dets[i]]->SetMarkerColor(kRed);
    profMedianPhi[dets[i]]->SetLineWidth(2.);
    profMedianPhi[dets[i]]->SetMarkerStyle(20);
    profMedianPhi[dets[i]]->SetMarkerSize(1.4);
    profMedianPhi[dets[i]]->Draw("E2same");
    DMRCanvPhi[i]->Draw();

    TArrow *l1=new TArrow(DMRCanvPhi[i]->GetUxmin(),0.,DMRCanvPhi[i]->GetUxmax(),0.,0.001);
    l1->SetLineColor(kBlack);
    l1->SetLineWidth(3);
    l1->SetLineStyle(9);

    l1->Draw("same");
    DMRCanvPhi[i]->SaveAs(Form("cProf_%s_RMS_DMR_Phi_",dets[i].Data())+fChainName+".pdf");
  
    fOutputFile->cd();
    DMRCanvPhi[i]->Write();

    // %%% - vs z - %%%

    makeNiceCanv(DMRCanvZ[i]);
    DMRCanvZ[i]->cd();
    gPad->SetRightMargin(0.12);
    makeNice2DPlotStyle(h2MedianZ[dets[i]]);
    h2MedianZ[dets[i]]->Draw("box");
    profMedianZ[dets[i]]->SetLineColor(kRed);
    profMedianZ[dets[i]]->SetFillColor(kRed);
    profMedianZ[dets[i]]->SetFillStyle(3002);
    profMedianZ[dets[i]]->SetMarkerColor(kRed);
    profMedianZ[dets[i]]->SetLineWidth(2.);
    profMedianZ[dets[i]]->SetMarkerStyle(20);
    profMedianZ[dets[i]]->SetMarkerSize(1.4);
    profMedianZ[dets[i]]->Draw("E2same");
    DMRCanvZ[i]->Draw();

    TArrow *l2=new TArrow(DMRCanvZ[i]->GetUxmin(),0.,DMRCanvZ[i]->GetUxmax(),0.,0.001);
    l2->SetLineColor(kBlack);
    l2->SetLineWidth(3);
    l2->SetLineStyle(9);

    l2->Draw("same");
    DMRCanvZ[i]->SaveAs(Form("cProf_%s_RMS_DMR_Z_",dets[i].Data())+fChainName+".pdf");

    fOutputFile->cd();
    DMRCanvZ[i]->Write();

  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // DMR RMS VS ETA OF MODULE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for (Int_t i=0;i<=h1MedianEta["BPix-x"]->GetNbinsX();i++){
    h1MedianEta["BPix-x"]->SetBinContent(i,profMedianEta["BPix-x"]->GetBinError(i));
    h1MedianEta["BPix-y"]->SetBinContent(i,profMedianEta["BPix-y"]->GetBinError(i));
    h1MedianEta["FPix-x"]->SetBinContent(i,profMedianEta["FPix-x"]->GetBinError(i));
    h1MedianEta["FPix-y"]->SetBinContent(i,profMedianEta["FPix-y"]->GetBinError(i));
    h1MedianEta["TIB"]->SetBinContent(i,profMedianEta["TIB"]->GetBinError(i));
    h1MedianEta["TID"]->SetBinContent(i,profMedianEta["TID"]->GetBinError(i));
    h1MedianEta["TOB"]->SetBinContent(i,profMedianEta["TOB"]->GetBinError(i));
    h1MedianEta["TEC"]->SetBinContent(i,profMedianEta["TEC"]->GetBinError(i));

    h1MedianEta["BPix-x"]->SetBinError(i,0.0001);
    h1MedianEta["BPix-y"]->SetBinError(i,0.0001);
    h1MedianEta["FPix-x"]->SetBinError(i,0.0001);
    h1MedianEta["FPix-y"]->SetBinError(i,0.0001);
    h1MedianEta["TIB"]->SetBinError(i,0.0001);
    h1MedianEta["TID"]->SetBinError(i,0.0001);
    h1MedianEta["TOB"]->SetBinError(i,0.0001);
    h1MedianEta["TEC"]->SetBinError(i,0.0001);
  }
  

  makeNicePlotStyleAndColor(h1MedianEta["BPix-x"],0);
  makeNicePlotStyleAndColor(h1MedianEta["BPix-y"],1);
  makeNicePlotStyleAndColor(h1MedianEta["FPix-x"],2);
  makeNicePlotStyleAndColor(h1MedianEta["FPix-y"],3);
  makeNicePlotStyleAndColor(h1MedianEta["TIB"],4);
  makeNicePlotStyleAndColor(h1MedianEta["TID"],5);
  makeNicePlotStyleAndColor(h1MedianEta["TOB"],6);
  makeNicePlotStyleAndColor(h1MedianEta["TEC"],7);

  TCanvas *cAllRMSVsEta = new TCanvas("cAllRMSVsEta","RMS of DMR vs |#eta| (all Tracker)",700,700);
  makeNiceCanv(cAllRMSVsEta);
  cAllRMSVsEta->cd();
  //h1MedianEta["BPix-x"]->GetYaxis()->SetRangeUser(0.1,7.);
  h1MedianEta["BPix-x"]->Draw("LP");
  h1MedianEta["BPix-y"]->Draw("LPsame");
  h1MedianEta["FPix-x"]->Draw("LPsame");
  h1MedianEta["FPix-y"]->Draw("LPsame");
  h1MedianEta["TIB"]->Draw("LPsame");
  h1MedianEta["TID"]->Draw("LPsame");
  h1MedianEta["TOB"]->Draw("LPsame");
  h1MedianEta["TEC"]->Draw("LPsame");
 
  TLegend *leg = new TLegend(0.14,0.7,0.54,0.88);
  leg-> SetNColumns(2);
  leg->SetLineColor(10);
  leg->SetFillColor(10);
  leg->SetTextFont(42);
  leg->AddEntry(h1MedianEta["BPix-x"],"BPix (x)");
  leg->AddEntry(h1MedianEta["FPix-x"],"FPix (x)");
  leg->AddEntry(h1MedianEta["BPix-y"],"BPix (y)");
  leg->AddEntry(h1MedianEta["FPix-y"],"FPix (y)");
  leg->AddEntry(h1MedianEta["TIB"],"TIB");
  leg->AddEntry(h1MedianEta["TID"],"TID");
  leg->AddEntry(h1MedianEta["TOB"],"TOB");
  leg->AddEntry(h1MedianEta["TEC"],"TEC");
  	
  leg->Draw("same");

  //cAllRMSVsEta->SaveAs("cALLRMSVsEta_"+fChainName+".root");
  fOutputFile->cd("rmsDMR");
  cAllRMSVsEta->Write();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // DMR RMS VS PHI OF MODULE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for (Int_t i=0;i<=h1MedianPhi["BPix-x"]->GetNbinsX();i++){
    h1MedianPhi["BPix-x"]->SetBinContent(i,profMedianPhi["BPix-x"]->GetBinError(i));
    h1MedianPhi["BPix-y"]->SetBinContent(i,profMedianPhi["BPix-y"]->GetBinError(i));
    h1MedianPhi["FPix-x"]->SetBinContent(i,profMedianPhi["FPix-x"]->GetBinError(i));
    h1MedianPhi["FPix-y"]->SetBinContent(i,profMedianPhi["FPix-y"]->GetBinError(i));
    h1MedianPhi["TIB"]->SetBinContent(i,profMedianPhi["TIB"]->GetBinError(i));
    h1MedianPhi["TID"]->SetBinContent(i,profMedianPhi["TID"]->GetBinError(i));
    h1MedianPhi["TOB"]->SetBinContent(i,profMedianPhi["TOB"]->GetBinError(i));
    h1MedianPhi["TEC"]->SetBinContent(i,profMedianPhi["TEC"]->GetBinError(i));

    h1MedianPhi["BPix-x"]->SetBinError(i,0.0001);
    h1MedianPhi["BPix-y"]->SetBinError(i,0.0001);
    h1MedianPhi["FPix-x"]->SetBinError(i,0.0001);
    h1MedianPhi["FPix-y"]->SetBinError(i,0.0001);
    h1MedianPhi["TIB"]->SetBinError(i,0.0001);
    h1MedianPhi["TID"]->SetBinError(i,0.0001);
    h1MedianPhi["TOB"]->SetBinError(i,0.0001);
    h1MedianPhi["TEC"]->SetBinError(i,0.0001);
  }
  

  makeNicePlotStyleAndColor(h1MedianPhi["BPix-x"],0);
  makeNicePlotStyleAndColor(h1MedianPhi["BPix-y"],1);
  makeNicePlotStyleAndColor(h1MedianPhi["FPix-x"],2);
  makeNicePlotStyleAndColor(h1MedianPhi["FPix-y"],3);
  makeNicePlotStyleAndColor(h1MedianPhi["TIB"],4);
  makeNicePlotStyleAndColor(h1MedianPhi["TID"],5);
  makeNicePlotStyleAndColor(h1MedianPhi["TOB"],6);
  makeNicePlotStyleAndColor(h1MedianPhi["TEC"],7);

  TCanvas *cAllRMSVsPhi = new TCanvas("cAllRMSVsPhi","RMS of DMR vs |#eta| (all Tracker)",700,700);
  makeNiceCanv(cAllRMSVsPhi);
  cAllRMSVsPhi->cd();
  //h1MedianPhi["BPix-x"]->GetYaxis()->SetRangeUser(0.1,7.);
  h1MedianPhi["BPix-x"]->Draw("LP");
  h1MedianPhi["BPix-y"]->Draw("LPsame");
  h1MedianPhi["FPix-x"]->Draw("LPsame");
  h1MedianPhi["FPix-y"]->Draw("LPsame");
  h1MedianPhi["TIB"]->Draw("LPsame");
  h1MedianPhi["TID"]->Draw("LPsame");
  h1MedianPhi["TOB"]->Draw("LPsame");
  h1MedianPhi["TEC"]->Draw("LPsame");
 
  TLegend *legPhi = new TLegend(0.14,0.7,0.54,0.88);
  legPhi-> SetNColumns(2);
  legPhi->SetLineColor(10);
  legPhi->SetFillColor(10);
  legPhi->SetTextFont(42);
  legPhi->AddEntry(h1MedianPhi["BPix-x"],"BPix (x)");
  legPhi->AddEntry(h1MedianPhi["FPix-x"],"FPix (x)");
  legPhi->AddEntry(h1MedianPhi["BPix-y"],"BPix (y)");
  legPhi->AddEntry(h1MedianPhi["FPix-y"],"FPix (y)");
  legPhi->AddEntry(h1MedianPhi["TIB"],"TIB");
  legPhi->AddEntry(h1MedianPhi["TID"],"TID");
  legPhi->AddEntry(h1MedianPhi["TOB"],"TOB");
  legPhi->AddEntry(h1MedianPhi["TEC"],"TEC");
  	
  legPhi->Draw("same");

  //cAllRMSVsPhi->SaveAs("cALLRMSVsPhi_"+fChainName+".root");
  fOutputFile->cd("rmsDMR");
  cAllRMSVsPhi->Write();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // DMR RMS VS Z OF MODULE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for (Int_t i=0;i<=h1MedianZ["BPix-x"]->GetNbinsX();i++){
    h1MedianZ["BPix-x"]->SetBinContent(i,profMedianZ["BPix-x"]->GetBinError(i));
    h1MedianZ["BPix-y"]->SetBinContent(i,profMedianZ["BPix-y"]->GetBinError(i));
    h1MedianZ["FPix-x"]->SetBinContent(i,profMedianZ["FPix-x"]->GetBinError(i));
    h1MedianZ["FPix-y"]->SetBinContent(i,profMedianZ["FPix-y"]->GetBinError(i));
    h1MedianZ["TIB"]->SetBinContent(i,profMedianZ["TIB"]->GetBinError(i));
    h1MedianZ["TID"]->SetBinContent(i,profMedianZ["TID"]->GetBinError(i));
    h1MedianZ["TOB"]->SetBinContent(i,profMedianZ["TOB"]->GetBinError(i));
    h1MedianZ["TEC"]->SetBinContent(i,profMedianZ["TEC"]->GetBinError(i));

    h1MedianZ["BPix-x"]->SetBinError(i,0.0001);
    h1MedianZ["BPix-y"]->SetBinError(i,0.0001);
    h1MedianZ["FPix-x"]->SetBinError(i,0.0001);
    h1MedianZ["FPix-y"]->SetBinError(i,0.0001);
    h1MedianZ["TIB"]->SetBinError(i,0.0001);
    h1MedianZ["TID"]->SetBinError(i,0.0001);
    h1MedianZ["TOB"]->SetBinError(i,0.0001);
    h1MedianZ["TEC"]->SetBinError(i,0.0001);
  }
  

  makeNicePlotStyleAndColor(h1MedianZ["BPix-x"],0);
  makeNicePlotStyleAndColor(h1MedianZ["BPix-y"],1);
  makeNicePlotStyleAndColor(h1MedianZ["FPix-x"],2);
  makeNicePlotStyleAndColor(h1MedianZ["FPix-y"],3);
  makeNicePlotStyleAndColor(h1MedianZ["TIB"],4);
  makeNicePlotStyleAndColor(h1MedianZ["TID"],5);
  makeNicePlotStyleAndColor(h1MedianZ["TOB"],6);
  makeNicePlotStyleAndColor(h1MedianZ["TEC"],7);

  TCanvas *cAllRMSVsZ = new TCanvas("cAllRMSVsZ","RMS of DMR vs |#eta| (all Tracker)",700,700);
  makeNiceCanv(cAllRMSVsZ);
  cAllRMSVsZ->cd();
  //h1MedianZ["BPix-x"]->GetYaxis()->SetRangeUser(0.1,7.);
  h1MedianZ["BPix-x"]->Draw("LP");
  h1MedianZ["BPix-y"]->Draw("LPsame");
  h1MedianZ["FPix-x"]->Draw("LPsame");
  h1MedianZ["FPix-y"]->Draw("LPsame");
  h1MedianZ["TIB"]->Draw("LPsame");
  h1MedianZ["TID"]->Draw("LPsame");
  h1MedianZ["TOB"]->Draw("LPsame");
  h1MedianZ["TEC"]->Draw("LPsame");
 
  TLegend *legZ = new TLegend(0.14,0.7,0.54,0.88);
  legZ-> SetNColumns(2);
  legZ->SetLineColor(10);
  legZ->SetFillColor(10);
  legZ->SetTextFont(42);
  legZ->AddEntry(h1MedianZ["BPix-x"],"BPix (x)");
  legZ->AddEntry(h1MedianZ["FPix-x"],"FPix (x)");
  legZ->AddEntry(h1MedianZ["BPix-y"],"BPix (y)");
  legZ->AddEntry(h1MedianZ["FPix-y"],"FPix (y)");
  legZ->AddEntry(h1MedianZ["TIB"],"TIB");
  legZ->AddEntry(h1MedianZ["TID"],"TID");
  legZ->AddEntry(h1MedianZ["TOB"],"TOB");
  legZ->AddEntry(h1MedianZ["TEC"],"TEC");
  	
  legZ->Draw("same");

  //cAllRMSVsZ->SaveAs("cALLRMSVsZ_"+fChainName+".root");
  fOutputFile->cd("rmsDMR");
  cAllRMSVsZ->Write();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for (UInt_t i=0;i<8;i++){
    
    // %%% - vs eta - %%%

    makeNiceCanv(PullCanvEta[i]);
    PullCanvEta[i]->cd();
    gPad->SetRightMargin(0.12);
    makeNice2DPlotStyle(h2PullEta[dets[i]]);
    h2PullEta[dets[i]]->Draw("colz");
    profPullEta[dets[i]]->SetLineColor(kBlue);
    //profPullEta[dets[i]]->SetFillColor(kBlue);
    //profPullEta[dets[i]]->SetFillStyle(3002);
    profPullEta[dets[i]]->SetMarkerColor(kBlue);
    profPullEta[dets[i]]->SetLineWidth(2.);
    profPullEta[dets[i]]->SetMarkerStyle(20);
    profPullEta[dets[i]]->SetMarkerSize(1.8);
    profPullEta[dets[i]]->Draw("E1same");
    PullCanvEta[i]->Draw();
    TArrow *l0=new TArrow(PullCanvEta[i]->GetUxmin(),1.,PullCanvEta[i]->GetUxmax(),1.,0.001);
    l0->SetLineColor(kBlack);
    l0->SetLineWidth(3);
    l0->SetLineStyle(9);
    l0->Draw("same");
    PullCanvEta[i]->SaveAs(Form("cPull_%s_RMSNorm_Eta_",dets[i].Data())+fChainName+".pdf");
  
    fOutputFile->cd();
    PullCanvEta[i]->Write();

    // %%% - vs phi - %%%

    makeNiceCanv(PullCanvPhi[i]);
    PullCanvPhi[i]->cd();
    gPad->SetRightMargin(0.12);
    makeNice2DPlotStyle(h2PullPhi[dets[i]]);
    h2PullPhi[dets[i]]->Draw("colz");
    profPullPhi[dets[i]]->SetLineColor(kBlue);
    //profPullPhi[dets[i]]->SetFillColor(kBlue);
    //profPullPhi[dets[i]]->SetFillStyle(3002);
    profPullPhi[dets[i]]->SetMarkerColor(kBlue);
    profPullPhi[dets[i]]->SetLineWidth(2.);
    profPullPhi[dets[i]]->SetMarkerStyle(20);
    profPullPhi[dets[i]]->SetMarkerSize(1.8);
    profPullPhi[dets[i]]->Draw("E1same");
    PullCanvPhi[i]->Draw();
    TArrow *l1=new TArrow(PullCanvPhi[i]->GetUxmin(),1.,PullCanvPhi[i]->GetUxmax(),1.,0.001);
    l1->SetLineColor(kBlack);
    l1->SetLineWidth(3);
    l1->SetLineStyle(9);
    l1->Draw("same");
    PullCanvPhi[i]->SaveAs(Form("cPull_%s_RMSNorm_Phi_",dets[i].Data())+fChainName+".pdf");
    
    fOutputFile->cd();
    PullCanvPhi[i]->Write();

    // %%% - vs z - %%%

    makeNiceCanv(PullCanvZ[i]);
    PullCanvZ[i]->cd();
    gPad->SetRightMargin(0.12);
    makeNice2DPlotStyle(h2PullZ[dets[i]]);
    h2PullZ[dets[i]]->Draw("colz");
    profPullZ[dets[i]]->SetLineColor(kBlue);
    //profPullZ[dets[i]]->SetFillColor(kBlue);
    //profPullZ[dets[i]]->SetFillStyle(3002);
    profPullZ[dets[i]]->SetMarkerColor(kBlue);
    profPullZ[dets[i]]->SetLineWidth(2.);
    profPullZ[dets[i]]->SetMarkerStyle(20);
    profPullZ[dets[i]]->SetMarkerSize(1.8);
    profPullZ[dets[i]]->Draw("E1same");
    PullCanvZ[i]->Draw();
    TArrow *l2=new TArrow(PullCanvZ[i]->GetUxmin(),1.,PullCanvZ[i]->GetUxmax(),1.,0.001);
    l2->SetLineColor(kBlack);
    l2->SetLineWidth(3);
    l2->SetLineStyle(9);
    l2->Draw("same");
    PullCanvZ[i]->SaveAs(Form("cPull_%s_RMSNorm_Z_",dets[i].Data())+fChainName+".pdf");
   
    fOutputFile->cd();
    PullCanvZ[i]->Write();

  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //% PULL VS ETA OF MODULE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TProfile* profPullEtaxBPix = (TProfile*)profPullEta["BPix-x"]->Clone("profPullEtaBPix-x"); 
  TProfile* profPullEtayBPix = (TProfile*)profPullEta["BPix-y"]->Clone("profPullEtaBPix-y"); 
  TProfile* profPullEtaxFPix = (TProfile*)profPullEta["FPix-x"]->Clone("profPullEtaFPix-x"); 
  TProfile* profPullEtayFPix = (TProfile*)profPullEta["FPix-y"]->Clone("profPullEtaFPix-y"); 
  TProfile* profPullEtaTIB   = (TProfile*)profPullEta["TIB"]->Clone("profPullEtaTIB");	     
  TProfile* profPullEtaTID   = (TProfile*)profPullEta["TID"]->Clone("profPullEtaTID");	     
  TProfile* profPullEtaTOB   = (TProfile*)profPullEta["TOB"]->Clone("profPullEtaTOB");	     
  TProfile* profPullEtaTEC   = (TProfile*)profPullEta["TEC"]->Clone("profPullEtaTEC");      

  makeNiceProfileStyleAndColor(profPullEtaxBPix,0);
  makeNiceProfileStyleAndColor(profPullEtayBPix,1);
  makeNiceProfileStyleAndColor(profPullEtaxFPix,2);
  makeNiceProfileStyleAndColor(profPullEtayFPix,3);
  makeNiceProfileStyleAndColor(profPullEtaTIB,4);
  makeNiceProfileStyleAndColor(profPullEtaTID,5);
  makeNiceProfileStyleAndColor(profPullEtaTOB,6);
  makeNiceProfileStyleAndColor(profPullEtaTEC,7);

  TCanvas *cAllRMSNormVsEta = new TCanvas("cAllRMSNormVsEta","RMS of pull vs |#eta| (all Tracker)",700,700);
  makeNiceCanv(cAllRMSNormVsEta);
  cAllRMSNormVsEta->cd();
  
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.14);
  //profPullEtaxBPix->GetYaxis()->SetRangeUser(0.5,1.5);
  profPullEtaxBPix->Draw("LP");	
  profPullEtayBPix->Draw("LPsame");
  profPullEtaxFPix->Draw("LPsame");
  profPullEtayFPix->Draw("LPsame");
  profPullEtaTIB->Draw("LPsame");
  profPullEtaTID->Draw("LPsame");
  profPullEtaTOB->Draw("LPsame");
  profPullEtaTEC->Draw("LPsame");

  TLegend *leg2Eta = new TLegend(0.16,0.75,0.56,0.93);
  leg2Eta-> SetNColumns(2);
  leg2Eta->SetLineColor(10);
  leg2Eta->SetFillColor(10);
  leg2Eta->SetTextFont(42);
  leg2Eta->AddEntry(profPullEtaxBPix,"BPix (x)");
  leg2Eta->AddEntry(profPullEtayBPix,"FPix (x)");
  leg2Eta->AddEntry(profPullEtaxFPix,"BPix (y)");
  leg2Eta->AddEntry(profPullEtayFPix,"FPix (y)");
  leg2Eta->AddEntry(profPullEtaTIB,"TIB");
  leg2Eta->AddEntry(profPullEtaTID,"TID");
  leg2Eta->AddEntry(profPullEtaTOB,"TOB");
  leg2Eta->AddEntry(profPullEtaTEC,"TEC");
  
  leg2Eta->Draw("same");

  //cAllRMSNormVsEta->SaveAs("cALLRMSNormVsEta_"+fChainName+".root");
  fOutputFile->cd("Pulls");
  cAllRMSNormVsEta->Write();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //% PULL VS PHI OF MODULE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TProfile* profPullPhixBPix = (TProfile*)profPullPhi["BPix-x"]->Clone("profPullPhiBPix-x"); 
  TProfile* profPullPhiyBPix = (TProfile*)profPullPhi["BPix-y"]->Clone("profPullPhiBPix-y"); 
  TProfile* profPullPhixFPix = (TProfile*)profPullPhi["FPix-x"]->Clone("profPullPhiFPix-x"); 
  TProfile* profPullPhiyFPix = (TProfile*)profPullPhi["FPix-y"]->Clone("profPullPhiFPix-y"); 
  TProfile* profPullPhiTIB   = (TProfile*)profPullPhi["TIB"]->Clone("profPullPhiTIB");	     
  TProfile* profPullPhiTID   = (TProfile*)profPullPhi["TID"]->Clone("profPullPhiTID");	     
  TProfile* profPullPhiTOB   = (TProfile*)profPullPhi["TOB"]->Clone("profPullPhiTOB");	     
  TProfile* profPullPhiTEC   = (TProfile*)profPullPhi["TEC"]->Clone("profPullPhiTEC");      

  makeNiceProfileStyleAndColor(profPullPhixBPix,0);
  makeNiceProfileStyleAndColor(profPullPhiyBPix,1);
  makeNiceProfileStyleAndColor(profPullPhixFPix,2);
  makeNiceProfileStyleAndColor(profPullPhiyFPix,3);
  makeNiceProfileStyleAndColor(profPullPhiTIB,4);
  makeNiceProfileStyleAndColor(profPullPhiTID,5);
  makeNiceProfileStyleAndColor(profPullPhiTOB,6);
  makeNiceProfileStyleAndColor(profPullPhiTEC,7);

  TCanvas *cAllRMSNormVsPhi = new TCanvas("cAllRMSNormVsPhi","RMS of pull vs |#eta| (all Tracker)",700,700);
  makeNiceCanv(cAllRMSNormVsPhi);
  cAllRMSNormVsPhi->cd();
  
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.14);
  //profPullPhixBPix->GetYaxis()->SetRangeUser(0.8,1.2);
  profPullPhixBPix->Draw("LP");	
  profPullPhiyBPix->Draw("LPsame");
  profPullPhixFPix->Draw("LPsame");
  profPullPhiyFPix->Draw("LPsame");
  profPullPhiTIB->Draw("LPsame");
  profPullPhiTID->Draw("LPsame");
  profPullPhiTOB->Draw("LPsame");
  profPullPhiTEC->Draw("LPsame");

  TLegend *leg2Phi = new TLegend(0.16,0.75,0.56,0.93);
  leg2Phi-> SetNColumns(2);
  leg2Phi->SetLineColor(10);
  leg2Phi->SetFillColor(10);
  leg2Phi->SetTextFont(42);
  leg2Phi->AddEntry(profPullPhixBPix,"BPix (x)");
  leg2Phi->AddEntry(profPullPhiyBPix,"FPix (x)");
  leg2Phi->AddEntry(profPullPhixFPix,"BPix (y)");
  leg2Phi->AddEntry(profPullPhiyFPix,"FPix (y)");
  leg2Phi->AddEntry(profPullPhiTIB,"TIB");
  leg2Phi->AddEntry(profPullPhiTID,"TID");
  leg2Phi->AddEntry(profPullPhiTOB,"TOB");
  leg2Phi->AddEntry(profPullPhiTEC,"TEC");
  
  leg2Phi->Draw("same");

  //cAllRMSNormVsPhi->SaveAs("cALLRMSNormVsPhi_"+fChainName+".root");
  fOutputFile->cd("Pulls");
  cAllRMSNormVsPhi->Write();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //% PULL VS Z OF MODULE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TProfile* profPullZxBPix = (TProfile*)profPullZ["BPix-x"]->Clone("profPullZBPix-x"); 
  TProfile* profPullZyBPix = (TProfile*)profPullZ["BPix-y"]->Clone("profPullZBPix-y"); 
  TProfile* profPullZxFPix = (TProfile*)profPullZ["FPix-x"]->Clone("profPullZFPix-x"); 
  TProfile* profPullZyFPix = (TProfile*)profPullZ["FPix-y"]->Clone("profPullZFPix-y"); 
  TProfile* profPullZTIB   = (TProfile*)profPullZ["TIB"]->Clone("profPullZTIB");	     
  TProfile* profPullZTID   = (TProfile*)profPullZ["TID"]->Clone("profPullZTID");	     
  TProfile* profPullZTOB   = (TProfile*)profPullZ["TOB"]->Clone("profPullZTOB");	     
  TProfile* profPullZTEC   = (TProfile*)profPullZ["TEC"]->Clone("profPullZTEC");      

  makeNiceProfileStyleAndColor(profPullZxBPix,0);
  makeNiceProfileStyleAndColor(profPullZyBPix,1);
  makeNiceProfileStyleAndColor(profPullZxFPix,2);
  makeNiceProfileStyleAndColor(profPullZyFPix,3);
  makeNiceProfileStyleAndColor(profPullZTIB,4);
  makeNiceProfileStyleAndColor(profPullZTID,5);
  makeNiceProfileStyleAndColor(profPullZTOB,6);
  makeNiceProfileStyleAndColor(profPullZTEC,7);

  TCanvas *cAllRMSNormVsZ = new TCanvas("cAllRMSNormVsZ","RMS of pull vs |#eta| (all Tracker)",700,700);
  makeNiceCanv(cAllRMSNormVsZ);
  cAllRMSNormVsZ->cd();
  
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.14);
  //profPullZxBPix->GetYaxis()->SetRangeUser(0.8,1.2);
  profPullZxBPix->Draw("LP");	
  profPullZyBPix->Draw("LPsame");
  profPullZxFPix->Draw("LPsame");
  profPullZyFPix->Draw("LPsame");
  profPullZTIB->Draw("LPsame");
  profPullZTID->Draw("LPsame");
  profPullZTOB->Draw("LPsame");
  profPullZTEC->Draw("LPsame");

  TLegend *leg2Z = new TLegend(0.16,0.75,0.56,0.93);
  leg2Z-> SetNColumns(2);
  leg2Z->SetLineColor(10);
  leg2Z->SetFillColor(10);
  leg2Z->SetTextFont(42);
  leg2Z->AddEntry(profPullZxBPix,"BPix (x)");
  leg2Z->AddEntry(profPullZyBPix,"FPix (x)");
  leg2Z->AddEntry(profPullZxFPix,"BPix (y)");
  leg2Z->AddEntry(profPullZyFPix,"FPix (y)");
  leg2Z->AddEntry(profPullZTIB,"TIB");
  leg2Z->AddEntry(profPullZTID,"TID");
  leg2Z->AddEntry(profPullZTOB,"TOB");
  leg2Z->AddEntry(profPullZTEC,"TEC");
  
  leg2Z->Draw("same");

  //cAllRMSNormVsZ->SaveAs("cALLRMSNormVsZ_"+fChainName+".root");
  fOutputFile->cd("Pulls");
  cAllRMSNormVsZ->Write();

}

