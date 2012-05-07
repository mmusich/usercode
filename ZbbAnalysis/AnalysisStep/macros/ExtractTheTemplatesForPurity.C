#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "TSystem.h"
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TList.h>
#include <TTree.h>
#include <TF1.h>
#include <TFile.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TMath.h>
#include <Zbbstruct4JEC.h>
#include <Math/VectorUtil.h>
#include <fstream>

using namespace std;

void ExtractTheTemplatesForPurity(TString tagger_="ssvhpt"){

  // file in output
  TString fileoutname =" templates_"+tagger_+".root";
  TFile* fileout = new TFile(fileoutname,"RECREATE");

  TH1F::SetDefaultSumw2(kTRUE);
  
  const int n_MCSamples(7);

  // Open files 
  TFile *datafile = TFile::Open("analyzePAT_DATA2011.root","READ");

  TFile *mcfile[n_MCSamples];
  mcfile[0] = TFile::Open("analyzePAT_MC_Zb5fToLL_All.root","READ");
  mcfile[1] = TFile::Open("analyzePAT_MC_ZcToLL_All.root","READ");
  mcfile[2] = TFile::Open("analyzePAT_MC_ZlJets_All.root","READ");
  mcfile[3] = TFile::Open("analyzePAT_MC_ZtautauJets_All.root","READ");
  mcfile[4] = TFile::Open("analyzePAT_MC_TTJets_All.root","READ");
  mcfile[5] = TFile::Open("analyzePAT_MC_ZW_All.root","READ");
  mcfile[6] = TFile::Open("analyzePAT_MC_ZZ_All.root","READ");

  TString theTreeName_ ="finaldistros_"+tagger_+"/JetBTag/ZbbNtuple";
  cout<<"The tree name is: "<<theTreeName_<<endl;

  // Get the n-tuples
  TTree *ntupleDATA;
  ntupleDATA = (TTree*)datafile->Get(theTreeName_);

  TTree *ntupleMC[n_MCSamples];
  for (Int_t i=0; i<n_MCSamples; i++) ntupleMC[i] = (TTree*)mcfile[i]->Get(theTreeName_);

  // Declare the structures
  Event_info    *Event    = new Event_info[n_MCSamples+1];
  Z_info        *Z        = new Z_info[n_MCSamples+1];
  bjet_info     *bjet     = new bjet_info[n_MCSamples+1];
  jet2_info     *jet2     = new jet2_info[n_MCSamples+1];
  MET_info      *MET      = new MET_info[n_MCSamples+1];
  btagjets_info *btagjets = new btagjets_info[n_MCSamples+1];

  // Set the branches adress
  ntupleDATA->SetBranchAddress("Event",&Event[0]);
  ntupleDATA->SetBranchAddress("Z",&Z[0]);
  ntupleDATA->SetBranchAddress("bjet",&bjet[0]);
  ntupleDATA->SetBranchAddress("jet2",&jet2[0]);
  ntupleDATA->SetBranchAddress("MET",&MET[0]);
  ntupleDATA->SetBranchAddress("btagjets",&btagjets[0]);

  for (Int_t i=0; i<n_MCSamples; i++){
    ntupleMC[i]->SetBranchAddress("Event",&Event[i+1]);
    ntupleMC[i]->SetBranchAddress("Z",&Z[i+1]);
    ntupleMC[i]->SetBranchAddress("bjet",&bjet[i+1]);
    ntupleMC[i]->SetBranchAddress("jet2",&jet2[i+1]);
    ntupleMC[i]->SetBranchAddress("MET",&MET[i+1]);
    ntupleMC[i]->SetBranchAddress("btagjets",&btagjets[i+1]); 
  }

  // Calculation of the luminosity weights
  Double_t lumi = 2200;   // integrated luminosity in pb
  Double_t lumiweights[n_MCSamples],xsections[n_MCSamples],nevts[n_MCSamples];
  for (Int_t i = 0; i < n_MCSamples; i++) nevts[i] = ((TH1F*)mcfile[i]->Get("analyzePat/Selevents/SelectedEvts"))->GetBinContent(1);

  xsections[0] = 3048.;
  xsections[1] = 3048.;
  xsections[2] = 3048.;
  xsections[3] = 3048.;
  xsections[4] = 165.;
  xsections[5] = 19.79;
  xsections[6] = 6.4;
  
  for (Int_t i = 0; i < n_MCSamples; i++) lumiweights[i] = lumi*xsections[i]/nevts[i];
  for (Int_t i = 0; i < n_MCSamples; i++) cout<<"The weight ["<<i<<"]"<<" is: "<< lumiweights[i] <<endl;

  //**************************
  // SCAN NTUPLES
  //**************************

  // Declare the histograms all flavours
  TH1F* h_JPdiscDATA_ = new TH1F("JPdiscdata","JP discriminant DATA; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscDATA_calib_ = new TH1F("JPdiscCalibdata","JP discriminant DATA; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassDATA_ = new TH1F("SVmassdata","secondary vertex mass DATA; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_nbjetsDATA_ = new TH1F("nbjetsdata","number of b-jets DATA; number of b-tagged jets; events",5,-0.5,4.5);
  TH1F* h_TCHEdiscDATA_ = new TH1F("TCHEdiscdata","TCHE discriminant DATA ; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscDATA_ = new TH1F("TCHPdiscdata","TCHP discriminant DATA ; TCHP discriminant; jets",80,-20.,60.);
  
  TH1F* h_JPdiscMC_b_ = new TH1F("JPdisc_mc_b","JP discriminant MC - b flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_b_calib_ = new TH1F("JPdiscCalib_mc_b","JP discriminant MC - b flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_b_ = new TH1F("SVmass_mc_b","secondary vertex mass MC - b flavor ; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_b_ =new TH1F("TCHEdisc_mc_b","TCHE discriminant MC - b flavor; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_b_ =new TH1F("TCHPdisc_mc_b","TCHP discriminant MC - b flavor; TCHP discriminant; jets",80,-20.,60.);
 
  TH1F* h_JPdiscMC_c_ = new TH1F("JPdisc_mc_c","JP discriminant MC - c flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_c_calib_ = new TH1F("JPdiscCalib_mc_c","JP discriminant MC - c flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_c_ = new TH1F("SVmass_mc_c","secondary vertex mass MC - c flavor; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_c_ = new TH1F("TCHEdisc_mc_c","TCHE discriminant MC - c flavor; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_c_ = new TH1F("TCHPdisc_mc_c","TCHP discriminant MC - c flavor; TCHP discriminant; jets",80,-20.,60.);

  TH1F* h_JPdiscMC_l_ = new TH1F("JPdisc_mc_l","JP discriminant MC - light flavours; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_l_calib_ = new TH1F("JPdiscCalib_mc_l","JP discriminant MC - light flavours; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_l_ = new TH1F("SVmass_mc_l","secondary vertex mass MC - light flavours; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_l_ = new TH1F("TCHEdisc_mc_l","TCHE discriminant MC - light flavours; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_l_ = new TH1F("TCHPdisc_mc_l","TCHP discriminant MC - light flavours; TCHP discriminant; jets",80,-20.,60.);

  // muons
  TH1F* h_JPdiscDATA_mm_ = new TH1F("JPdiscdata_mm","JP discriminant DATA; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscDATA_calib_mm_ = new TH1F("JPdiscCalibdata_mm","JP discriminant DATA; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassDATA_mm_ = new TH1F("SVmassdata_mm","secondary vertex mass DATA; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_nbjetsDATA_mm_ = new TH1F("nbjetsdata_mm","number of b-jets DATA; number of b-tagged jets; events",5,-0.5,4.5);
  TH1F* h_TCHEdiscDATA_mm_ = new TH1F("TCHEdiscdata_mm","TCHE discriminant DATA muon; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscDATA_mm_ = new TH1F("TCHPdiscdata_mm","TCHP discriminant DATA muon; TCHP discriminant; jets",80,-20.,60.);

  TH1F* h_JPdiscMC_b_mm_ = new TH1F("JPdisc_mc_b_mm","JP discriminant MC - b flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_b_calib_mm_ = new TH1F("JPdiscCalib_mc_b_mm","JP discriminant MC - b flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_b_mm_ = new TH1F("SVmass_mc_b_mm","secondary vertex mass MC - b flavor ; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_b_mm_ =new TH1F("TCHEdisc_mc_b_mm","TCHE discriminant MC - b flavor; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_b_mm_ =new TH1F("TCHPdisc_mc_b_mm","TCHP discriminant MC - b flavor; TCHP discriminant; jets",80,-20.,60.);
   
  TH1F* h_JPdiscMC_c_mm_ = new TH1F("JPdisc_mc_c_mm","JP discriminant MC - c flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_c_calib_mm_ = new TH1F("JPdiscCalib_mc_c_mm","JP discriminant MC - c flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_c_mm_ = new TH1F("SVmass_mc_c_mm","secondary vertex mass MC - c flavor; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_c_mm_ = new TH1F("TCHEdisc_mc_c_mm","TCHE discriminant MC - c flavor; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_c_mm_ = new TH1F("TCHPdisc_mc_c_mm","TCHP discriminant MC - c flavor; TCHP discriminant; jets",80,-20.,60.);

  TH1F* h_JPdiscMC_l_mm_ = new TH1F("JPdisc_mc_l_mm","JP discriminant MC - light flavours; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_l_calib_mm_ = new TH1F("JPdiscCalib_mc_l_mm","JP discriminant MC - light flavours; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_l_mm_ = new TH1F("SVmass_mc_l_mm","secondary vertex mass MC - light flavours; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_l_mm_ = new TH1F("TCHEdisc_mc_l_mm","TCHE discriminant MC - light flavours; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_l_mm_ = new TH1F("TCHPdisc_mc_l_mm","TCHP discriminant MC - light flavours; TCHP discriminant; jets",80,-20.,60.);

  // electrons
  TH1F* h_JPdiscDATA_ee_ = new TH1F("JPdiscdata_ee","JP discriminant DATA; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscDATA_calib_ee_ = new TH1F("JPdiscCalibdata_ee","JP discriminant DATA; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassDATA_ee_ = new TH1F("SVmassdata_ee","secondary vertex mass DATA; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_nbjetsDATA_ee_ = new TH1F("nbjetsdata_ee","number of b-jets DATA; number of b-tagged jets; events",5,-0.5,4.5);
  TH1F* h_TCHEdiscDATA_ee_ = new TH1F("TCHEdiscdata_ee","TCHE discriminant DATA electron; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscDATA_ee_ = new TH1F("TCHPdiscdata_ee","TCHP discriminant DATA electron; TCHP discriminant; jets",80,-20.,60.);

  TH1F* h_JPdiscMC_b_ee_ = new TH1F("JPdisc_mc_b_ee","JP discriminant MC - b flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_b_calib_ee_ = new TH1F("JPdiscCalib_mc_b_ee","JP discriminant MC - b flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_b_ee_ = new TH1F("SVmass_mc_b_ee","secondary vertex mass MC - b flavor ; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_b_ee_ =new TH1F("TCHEdisc_mc_b_ee","TCHE discriminant MC - b flavor; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_b_ee_ =new TH1F("TCHPdisc_mc_b_ee","TCHP discriminant MC - b flavor; TCHP discriminant; jets",80,-20.,60.);
 
  TH1F* h_JPdiscMC_c_ee_ = new TH1F("JPdisc_mc_c_ee","JP discriminant MC - c flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_c_calib_ee_ = new TH1F("JPdiscCalib_mc_c_ee","JP discriminant MC - c flavor; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_c_ee_ = new TH1F("SVmass_mc_c_ee","secondary vertex mass MC - c flavor; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_c_ee_ = new TH1F("TCHEdisc_mc_c_ee","TCHE discriminant MC - c flavor; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_c_ee_ = new TH1F("TCHPdisc_mc_c_ee","TCHP discriminant MC - c flavor; TCHP discriminant; jets",80,-20.,60.);

  TH1F* h_JPdiscMC_l_ee_ = new TH1F("JPdisc_mc_l_ee","JP discriminant MC - light flavours; JP discriminant; jets",80,0.,3.);
  TH1F* h_JPdiscMC_l_calib_ee_ = new TH1F("JPdiscCalib_mc_l_ee","JP discriminant MC - light flavours; JP discriminant; jets",80,0.,3.);
  TH1F* h_SVmassMC_l_ee_ = new TH1F("SVmass_mc_l_ee","secondary vertex mass MC - light flavours; secondary vertex mass M_{SV} (GeV); jets",80,0.,5.);
  TH1F* h_TCHEdiscMC_l_ee_ = new TH1F("TCHEdisc_mc_l_ee","TCHE discriminant MC - light flavours; TCHE discriminant; jets",80,-20.,60.);
  TH1F* h_TCHPdiscMC_l_ee_ = new TH1F("TCHPdisc_mc_l_ee","TCHP discriminant MC - light flavours; TCHP discriminant; jets",80,-20.,60.);

  // Get the number of entries
  Int_t nEntries_DATA = ntupleDATA->GetEntries();
  Int_t nEntries_MC[n_MCSamples];
  for (Int_t i=0; i<n_MCSamples; i++) nEntries_MC[i] = ntupleMC[i]->GetEntries();

  //Scan DATA ntuple  
  for (Int_t iEntry=0; iEntry<nEntries_DATA; ++iEntry){ 
    ntupleDATA->GetEntry(iEntry);
    Int_t nBjets =  Event[0].nbjet_Event;
    h_nbjetsDATA_->Fill(nBjets);
    if(Event[0].isMuon){
      h_nbjetsDATA_mm_->Fill(nBjets);
    } else {
      h_nbjetsDATA_ee_->Fill(nBjets);
    }  
    for(Int_t bj=0; bj< nBjets; bj++){  
      h_JPdiscDATA_->Fill(btagjets[0].jpdiscr_bjets[bj]);
      h_JPdiscDATA_calib_->Fill(btagjets[0].jpdiscr_fromtags_bjets[bj]); 
      h_SVmassDATA_->Fill(btagjets[0].svmass_bjets[bj]);
      h_TCHEdiscDATA_->Fill(btagjets[0].tchediscr_bjets[bj]);
      h_TCHPdiscDATA_->Fill(btagjets[0].tchpdiscr_bjets[bj]);
      if(Event[0].isMuon){
	h_JPdiscDATA_mm_->Fill(btagjets[0].jpdiscr_bjets[bj]);
	h_JPdiscDATA_calib_mm_->Fill(btagjets[0].jpdiscr_fromtags_bjets[bj]); 
	h_SVmassDATA_mm_->Fill(btagjets[0].svmass_bjets[bj]);
	h_TCHEdiscDATA_mm_->Fill(btagjets[0].tchediscr_bjets[bj]);
        h_TCHPdiscDATA_mm_->Fill(btagjets[0].tchpdiscr_bjets[bj]);  
      } else {
	h_JPdiscDATA_ee_->Fill(btagjets[0].jpdiscr_bjets[bj]);
	h_JPdiscDATA_calib_ee_->Fill(btagjets[0].jpdiscr_fromtags_bjets[bj]); 
	h_SVmassDATA_ee_->Fill(btagjets[0].svmass_bjets[bj]);
	h_TCHEdiscDATA_ee_->Fill(btagjets[0].tchediscr_bjets[bj]);
        h_TCHPdiscDATA_ee_->Fill(btagjets[0].tchpdiscr_bjets[bj]);  
      }
    }
  }

  //Scan MC ntuple
  for (Int_t k=0; k<n_MCSamples; k++){ //loop over MC samples
    for (Int_t iEntry=0; iEntry<nEntries_MC[k]; ++iEntry){
      ntupleMC[k]->GetEntry(iEntry);
      Float_t w = Event[k+1].weight_Event*lumiweights[k];
      Int_t nBjets =  Event[k+1].nbjet_Event;
      for(Int_t bj=0; bj< nBjets; bj++){  
	if(fabs(btagjets[k+1].trueflavour_bjets[bj])==5) {                // b jets
	  h_JPdiscMC_b_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	  h_JPdiscMC_b_calib_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w); 
	  h_SVmassMC_b_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	  h_TCHEdiscMC_b_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	  h_TCHPdiscMC_b_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  if(Event[k+1].isMuon){
	    h_JPdiscMC_b_mm_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w);
	    h_JPdiscMC_b_calib_mm_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w);  
	    h_SVmassMC_b_mm_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	    h_TCHEdiscMC_b_mm_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	    h_TCHPdiscMC_b_mm_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  } else {
	    h_JPdiscMC_b_ee_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	    h_JPdiscMC_b_calib_ee_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w);  
	    h_SVmassMC_b_ee_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	    h_TCHEdiscMC_b_ee_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	    h_TCHPdiscMC_b_ee_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  }
	} else if (fabs(btagjets[k+1].trueflavour_bjets[bj])==4) {        // c jets
	  h_JPdiscMC_c_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	  h_JPdiscMC_c_calib_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w); 
	  h_SVmassMC_c_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	  h_TCHEdiscMC_c_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	  h_TCHPdiscMC_c_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  if(Event[k+1].isMuon){
	    h_JPdiscMC_c_mm_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	    h_JPdiscMC_c_calib_mm_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w);
	    h_SVmassMC_c_mm_->Fill(btagjets[k+1].svmass_bjets[bj],w); 
	    h_TCHEdiscMC_c_mm_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	    h_TCHPdiscMC_c_mm_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  } else {
	    h_JPdiscMC_c_ee_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	    h_JPdiscMC_c_calib_ee_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w);
	    h_SVmassMC_c_ee_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	    h_TCHEdiscMC_c_ee_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	    h_TCHPdiscMC_c_ee_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  }
	} else {                                                          // light jets
	  h_JPdiscMC_l_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	  h_JPdiscMC_l_calib_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w); 
	  h_SVmassMC_l_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	  h_TCHEdiscMC_l_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	  h_TCHPdiscMC_l_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  if(Event[k+1].isMuon){
	    h_JPdiscMC_l_mm_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w); 
	    h_JPdiscMC_l_calib_mm_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w); 
	    h_SVmassMC_l_mm_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	    h_TCHEdiscMC_l_mm_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	    h_TCHPdiscMC_l_mm_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);
	  } else {
	    h_JPdiscMC_l_ee_->Fill(btagjets[k+1].jpdiscr_bjets[bj],w);
	    h_JPdiscMC_l_calib_ee_->Fill(btagjets[k+1].jpdiscr_fromtags_bjets[bj],w); 
	    h_SVmassMC_l_ee_->Fill(btagjets[k+1].svmass_bjets[bj],w);
	    h_TCHEdiscMC_l_ee_->Fill(btagjets[k+1].tchediscr_bjets[bj],w);
	    h_TCHPdiscMC_l_ee_->Fill(btagjets[k+1].tchpdiscr_bjets[bj],w);  
	  }
	}
      }//ends loop on b tagged jets
    }// ends loop on tree
  }//ends loop on MC samples
  
  // Write the histograms into the output file
  
  fileout->cd();
  
  // sum of flavours

  h_JPdiscDATA_->Write(); 
  h_JPdiscDATA_calib_->Write(); 
  h_SVmassDATA_->Write();  
  h_nbjetsDATA_->Write();
  h_TCHEdiscDATA_->Write(); 
  h_TCHPdiscDATA_->Write(); 
  		      
  h_JPdiscMC_b_->Write();  
  h_JPdiscMC_b_calib_->Write(); 
  h_SVmassMC_b_->Write();
  h_TCHEdiscMC_b_->Write(); 
  h_TCHPdiscMC_b_->Write();   
 	      
  h_JPdiscMC_c_->Write();  
  h_JPdiscMC_c_calib_->Write(); 
  h_SVmassMC_c_->Write(); 
  h_TCHEdiscMC_c_->Write();
  h_TCHPdiscMC_c_->Write();
  		      
  h_JPdiscMC_l_->Write();  
  h_JPdiscMC_l_calib_->Write(); 
  h_SVmassMC_l_->Write(); 
  h_TCHEdiscMC_l_->Write();
  h_TCHPdiscMC_l_->Write();  
  
  // muon channel

  h_JPdiscDATA_mm_->Write(); 
  h_JPdiscDATA_calib_mm_->Write();
  h_SVmassDATA_mm_->Write();  
  h_nbjetsDATA_mm_->Write();
  h_TCHEdiscDATA_mm_->Write(); 
  h_TCHPdiscDATA_mm_->Write();

  h_JPdiscMC_b_mm_->Write(); 
  h_JPdiscMC_b_calib_mm_->Write(); 
  h_SVmassMC_b_mm_->Write();  
  h_TCHEdiscMC_b_mm_->Write(); 
  h_TCHPdiscMC_b_mm_->Write();   
  		      
  h_JPdiscMC_c_mm_->Write();  
  h_JPdiscMC_c_calib_mm_->Write();
  h_SVmassMC_c_mm_->Write();
  h_TCHEdiscMC_c_mm_->Write();
  h_TCHPdiscMC_c_mm_->Write();
  		      
  h_JPdiscMC_l_mm_->Write(); 
  h_JPdiscMC_l_calib_mm_->Write();
  h_SVmassMC_l_mm_->Write(); 
  h_TCHEdiscMC_l_mm_->Write();
  h_TCHPdiscMC_l_mm_->Write();  
  
  // electron channel

  h_JPdiscDATA_ee_->Write(); 
  h_JPdiscDATA_calib_ee_->Write();
  h_SVmassDATA_ee_->Write();  
  h_nbjetsDATA_ee_->Write();
  h_TCHEdiscDATA_ee_->Write(); 
  h_TCHPdiscDATA_ee_->Write();

  h_JPdiscMC_b_ee_->Write();  
  h_JPdiscMC_b_calib_ee_->Write(); 
  h_SVmassMC_b_ee_->Write(); 
  h_TCHEdiscMC_b_ee_->Write(); 
  h_TCHPdiscMC_b_ee_->Write();   
  		      
  h_JPdiscMC_c_ee_->Write();
  h_JPdiscMC_c_calib_ee_->Write(); 
  h_SVmassMC_c_ee_->Write();
  h_TCHEdiscMC_c_ee_->Write();
  h_TCHPdiscMC_c_ee_->Write();
  		      
  h_JPdiscMC_l_ee_->Write(); 
  h_JPdiscMC_l_calib_ee_->Write(); 
  h_SVmassMC_l_ee_->Write(); 
  h_TCHEdiscMC_l_ee_->Write();
  h_TCHPdiscMC_l_ee_->Write();  
  
  fileout->Close();
   
}
