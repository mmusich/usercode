#include <iostream>
#include <vector>
#include "TBranch.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"
#include "TLeafI.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include <fstream>
#include <iomanip>

using namespace std;

//---------------------------------------------------------------------------------------------
// struct for event efficiency calculation

struct myEff
{

  // in variables
  Double_t sumW_Num_,sumW_Den_,sumW2_Num_,sumW2_Den_;
  // out variables
  Double_t eff_, errCP_, errBin_, errImp_;

  myEff(Double_t sumW_Num,Double_t sumW_Den,Double_t sumW2_Num,Double_t sumW2_Den){

    // initialization
    sumW_Num_ =sumW_Num;
    sumW_Den_ =sumW_Den;
    sumW2_Num_=sumW2_Num; 
    sumW2_Den_=sumW2_Den;

    // calculation
    if(sumW_Den_==0 || sumW_Num_>=sumW_Den_){
      std::cout<<"myEff::myEff("<<sumW_Num<<","<< sumW_Den <<","<<sumW2_Num <<","<< sumW2_Den <<") : WARNING --- Something wrong is going on ---"<<std::endl;
    } else {
      this->calculate();
    }
  }
  
  void reset(){
    sumW_Num_=-1. ;
    sumW_Den_=-1. ;
    sumW2_Num_=-999.; 
    sumW2_Den_=-999.; 
    eff_=0.;errCP_=0.;errBin_=0.;errImp_=0.;    
  }

  void calculate(){
    eff_   = sumW_Num_/sumW_Den_;
    
    TEfficiency theEff;
    errCP_  = theEff.ClopperPearson(sumW_Den_+0.5,sumW_Num_+0.5,0.683,1)-eff_;
    errBin_ = TMath::Sqrt((eff_*(1-eff_))/sumW_Den_);
    
    if(sumW2_Num_!=0 && sumW2_Den_!=0){
      errImp_ = TMath::Sqrt(sumW2_Num_*(sumW_Den_-sumW_Num_)*(sumW_Den_-sumW_Num_)+(sumW2_Den_-sumW2_Num_)*sumW_Num_*sumW_Num_)/(sumW_Den_*sumW_Den_);
    }
  }
};


//---------------------------------------------------------------------------------------------
// method for nice graphics

void MakeNiceHistoStyle(TH1 *hist){
  
  // hist->SetTitleSize(0.09); 
  hist->SetMinimum(0.01);
  //hist->SetMaximum(theMaximum*1.10);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42);  
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetTitleOffset(1);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetLabelSize(0.05); 
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
}

//---------------------------------------------------------------------------------------------
// main method

void ScanTheTreeForCorrectionFactors(TString filename1 = "Acceptance_ZbbToEE_MADGRAPH.root",Bool_t useClopperPearsonErrors_=true, Bool_t partonLevel_=false,Bool_t showUnweighted_=false){

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //  Graphic settings
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //gStyle->SetOptStat(00000000);
  TH1::StatOverflows(kTRUE);
  //gROOT->SetStyle("Pub");
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetOptStat("emr");
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetNdivisions(303);  
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05); //---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1);        // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0); 
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);

  TH1F *h_GEN_mass_ZfromLeptons_ = new TH1F("GEN_mass_ZfromMuons","mass of ll (genParticles); GEN M_{l^{+}l^{-}} (GeV)",100.,0.,150.);

  TFile *file = new TFile(filename1,"open");

  TString ntp_addr_;

  if(partonLevel_){
    ntp_addr_=Form("calculate%sAcceptance/T","Parton");
  } else {
    ntp_addr_=Form("calculate%sAcceptance/T","Hadron");
  }

  TTree* myNTuple_ = NULL;
  if((TTree*)file->Get(ntp_addr_)){
    myNTuple_ = (TTree*)file->Get(ntp_addr_);
  } else {
    myNTuple_ = (TTree*)file->Get("calculateAcceptance/T");
  }

  Bool_t isRightFlavour_,isgenZyes_,isgenbyes_,isgenjetyes_,isgenZkinyes_,isrecZkinyes_,isrecbyes_,isrecjetyes_;
  Bool_t is_rec_lep_idiso_yes_,is_rec_b_HE_yes_,is_rec_b_HP_yes_; 
  Int_t nbgenJets_,nFilteredBgenJets_,nrecoJets_,nFilteredRecoJets_;
  Float_t llMass_,evtweight_,evtweightSq_,lepIdIsoWeight_,lepIdIsoWeightSq_,bTagWeightHE_,bTagWeightHESq_,bTagWeightHP_,bTagWeightHPSq_;;

  myNTuple_->SetBranchAddress("isRightFlavour",&isRightFlavour_);     
  myNTuple_->SetBranchAddress("isgenZyes",&isgenZyes_);          
  myNTuple_->SetBranchAddress("isgenbyes",&isgenbyes_);          
  myNTuple_->SetBranchAddress("isgenjetyes",&isgenjetyes_);          
  myNTuple_->SetBranchAddress("isgenZkinyes",&isgenZkinyes_);       
  myNTuple_->SetBranchAddress("isrecZkinyes",&isrecZkinyes_);       
  myNTuple_->SetBranchAddress("isrecbyes",&isrecbyes_);
  myNTuple_->SetBranchAddress("isrecjetyes",&isrecjetyes_);       
  myNTuple_->SetBranchAddress("nbgenJets",&nbgenJets_);          
  myNTuple_->SetBranchAddress("nFilteredBgenJets",&nFilteredBgenJets_);  
  myNTuple_->SetBranchAddress("nrecoJets",&nrecoJets_);          
  myNTuple_->SetBranchAddress("nFilteredRecoJets",&nFilteredRecoJets_);  
  myNTuple_->SetBranchAddress("llMass",&llMass_);             
  myNTuple_->SetBranchAddress("evtweight",&evtweight_);          
  myNTuple_->SetBranchAddress("evtweightSq",&evtweightSq_);      
  myNTuple_->SetBranchAddress("lepIdIsoWeight",&lepIdIsoWeight_);
  myNTuple_->SetBranchAddress("lepIdIsoWeightSq",&lepIdIsoWeightSq_);
  myNTuple_->SetBranchAddress("bTagWeightHE",&bTagWeightHE_);  
  myNTuple_->SetBranchAddress("bTagWeightHESq",&bTagWeightHESq_);
  myNTuple_->SetBranchAddress("bTagWeightHP",&bTagWeightHP_);  
  myNTuple_->SetBranchAddress("bTagWeightHPSq",&bTagWeightHPSq_);
  myNTuple_->SetBranchAddress("isreclepidiso_yes",&is_rec_lep_idiso_yes_);
  myNTuple_->SetBranchAddress("isrecbHE_yes",&is_rec_b_HE_yes_);
  myNTuple_->SetBranchAddress("isrecbHP_yes",&is_rec_b_HP_yes_);

  Int_t numEntries =  myNTuple_->GetEntries();

  //* variables for A_l and C_had
  Float_t N_genZkin_yes_genb_yes(0);
  Float_t sumW_genZkin_yes_genb_yes(0.);
  Float_t sumW2_genZkin_yes_genb_yes(0.);

  Float_t N_genZ_yes_genb_yes(0);
  Float_t sumW_genZ_yes_genb_yes(0.);
  Float_t sumW2_genZ_yes_genb_yes(0.);

  Float_t N_recoZkin_yes_recob_yes(0);
  Float_t sumW_recoZkin_yes_recob_yes(0.);
  Float_t sumW2_recoZkin_yes_recob_yes(0.);
  
  //* variables for C_reco (Z+b/Z+j)
  Float_t N_genZkin_yes_genjet_yes(0);
  Float_t sumW_genZkin_yes_genjet_yes(0.);
  Float_t sumW2_genZkin_yes_genjet_yes(0.);

  Float_t N_recoZkin_yes_recojet_yes(0);
  Float_t sumW_recoZkin_yes_recojet_yes(0.);
  Float_t sumW2_recoZkin_yes_recojet_yes(0.);

  Float_t N_recoZkinLepIdIso_yes_recojet_yes(0);
  Float_t sumW_recoZkinLepIdIso_yes_recojet_yes(0.);
  Float_t sumW2_recoZkinLepIdIso_yes_recojet_yes(0.);

  //* variables for Epsilon_l and Epsilon_b 
  Float_t N_recoZkinLepIdIso_yes_recob_yes(0);
  Float_t sumW_recoZkinLepIdIso_yes_recob_yes(0.);
  Float_t sumW2_recoZkinLepIdIso_yes_recob_yes(0.);

  Float_t N_recoZkinLepIdIso_yes_recobHE_yes(0);
  Float_t sumW_recoZkinLepIdIso_yes_recobHE_yes(0.);
  Float_t sumW2_recoZkinLepIdIso_yes_recobHE_yes(0.);

  Float_t N_recoZkinLepIdIso_yes_recobHP_yes(0);
  Float_t sumW_recoZkinLepIdIso_yes_recobHP_yes(0.);
  Float_t sumW2_recoZkinLepIdIso_yes_recobHP_yes(0.);

  for (Int_t jEntry=0; jEntry<numEntries; ++jEntry) {

    myNTuple_->GetEntry(jEntry);

    /* for A_l and C_had */
    if(isgenbyes_ && isgenZkinyes_) {
      N_genZkin_yes_genb_yes+=1.; 
      sumW_genZkin_yes_genb_yes+=evtweight_;
      sumW2_genZkin_yes_genb_yes+=pow(evtweight_,2);
    }
    
    if(isgenbyes_ && isgenZyes_){
      N_genZ_yes_genb_yes+=1.;
      sumW_genZ_yes_genb_yes+=evtweight_;
      sumW2_genZ_yes_genb_yes+=pow(evtweight_,2);
    }
    
    if(isrecbyes_ && isrecZkinyes_ ){
      N_recoZkin_yes_recob_yes+=1.;  
      sumW_recoZkin_yes_recob_yes+=evtweight_;
      sumW2_recoZkin_yes_recob_yes+=pow(evtweight_,2);
      
      if(is_rec_lep_idiso_yes_){
	N_recoZkinLepIdIso_yes_recob_yes+=1.;	   
	sumW_recoZkinLepIdIso_yes_recob_yes+=(evtweight_*lepIdIsoWeight_); 
	sumW2_recoZkinLepIdIso_yes_recob_yes+=pow(evtweight_*lepIdIsoWeight_,2);
	
	if(is_rec_b_HE_yes_){
	  N_recoZkinLepIdIso_yes_recobHE_yes+=1.;  
	  sumW_recoZkinLepIdIso_yes_recobHE_yes+=(evtweight_*lepIdIsoWeight_*bTagWeightHE_); 
	  sumW2_recoZkinLepIdIso_yes_recobHE_yes+=pow(evtweight_*lepIdIsoWeight_*bTagWeightHE_,2);
	}
	if(is_rec_b_HP_yes_){
	  N_recoZkinLepIdIso_yes_recobHP_yes+=1.;     
	  sumW_recoZkinLepIdIso_yes_recobHP_yes+=(evtweight_*lepIdIsoWeight_*bTagWeightHP_); 
	  sumW2_recoZkinLepIdIso_yes_recobHP_yes+=pow(evtweight_*lepIdIsoWeight_*bTagWeightHP_,2);
	}
      }
    }

    /* for C_reco (Z+b/Z+j)*/
    if(isgenjetyes_ && isgenZkinyes_){
      N_genZkin_yes_genjet_yes+=1.;
      sumW_genZkin_yes_genjet_yes+=evtweight_;
      sumW2_genZkin_yes_genjet_yes+=pow(evtweight_,2);
      h_GEN_mass_ZfromLeptons_->Fill(llMass_,evtweight_);
    }
    
    if(isrecjetyes_ && isrecZkinyes_){
      N_recoZkin_yes_recojet_yes+=1.;  
      sumW_recoZkin_yes_recojet_yes+=evtweight_;
      sumW2_recoZkin_yes_recojet_yes+=pow(evtweight_,2);

      if(is_rec_lep_idiso_yes_){
	N_recoZkinLepIdIso_yes_recojet_yes+=1.;
	sumW_recoZkinLepIdIso_yes_recojet_yes+=(evtweight_*lepIdIsoWeight_);
	sumW2_recoZkinLepIdIso_yes_recojet_yes+=pow(evtweight_*lepIdIsoWeight_,2);;
      }
    }
  }
  
  // Definition of efficiencies objects
  myEff Al_unw(N_genZkin_yes_genb_yes,N_genZ_yes_genb_yes,0.,0.);
  myEff Al_wgt(sumW_genZkin_yes_genb_yes,sumW_genZ_yes_genb_yes,sumW_genZkin_yes_genb_yes,sumW_genZ_yes_genb_yes);
  myEff Chad_unw(N_recoZkin_yes_recob_yes,N_genZkin_yes_genb_yes,0.,0.);
  myEff Chad_wgt(sumW_recoZkin_yes_recob_yes,sumW_genZkin_yes_genb_yes,sumW2_recoZkin_yes_recob_yes,sumW2_genZkin_yes_genb_yes);
  myEff ChadAny_unw(N_recoZkin_yes_recojet_yes,N_genZkin_yes_genjet_yes,0.,0.);
  myEff ChadAny_wgt(sumW_recoZkin_yes_recojet_yes,sumW_genZkin_yes_genjet_yes,sumW2_recoZkin_yes_recojet_yes,sumW2_genZkin_yes_genjet_yes);

  myEff Epsilonl_unw(N_recoZkinLepIdIso_yes_recob_yes,N_recoZkin_yes_recob_yes,0.,0.);    
  myEff EpsilonlAny_unw(N_recoZkinLepIdIso_yes_recojet_yes,N_recoZkin_yes_recojet_yes,0.,0.);     
  myEff EpsilonbHE_unw(N_recoZkinLepIdIso_yes_recobHE_yes,N_recoZkinLepIdIso_yes_recob_yes,0.,0.);
  myEff EpsilonbHP_unw(N_recoZkinLepIdIso_yes_recobHP_yes,N_recoZkinLepIdIso_yes_recob_yes,0.,0.);
  myEff Epsilonl_wgt(sumW_recoZkinLepIdIso_yes_recob_yes,sumW_recoZkin_yes_recob_yes,sumW2_recoZkinLepIdIso_yes_recob_yes,sumW2_recoZkin_yes_recob_yes);   
  myEff EpsilonlAny_wgt(sumW_recoZkinLepIdIso_yes_recojet_yes,sumW_recoZkin_yes_recojet_yes,sumW2_recoZkinLepIdIso_yes_recojet_yes,sumW2_recoZkin_yes_recojet_yes);    
  myEff EpsilonbHE_wgt(sumW_recoZkinLepIdIso_yes_recobHE_yes,sumW_recoZkinLepIdIso_yes_recob_yes,sumW2_recoZkinLepIdIso_yes_recobHE_yes,sumW2_recoZkinLepIdIso_yes_recob_yes);
  myEff EpsilonbHP_wgt(sumW_recoZkinLepIdIso_yes_recobHP_yes,sumW_recoZkinLepIdIso_yes_recob_yes,sumW2_recoZkinLepIdIso_yes_recobHP_yes,sumW2_recoZkinLepIdIso_yes_recob_yes);

  //   Definition of correction factors
  Double_t UnwAl          = Al_unw.eff_;
  Double_t UnwCHad        = Chad_unw.eff_;
  Double_t UnwCHad_Any    = ChadAny_unw.eff_;
  
  Double_t UnwEpsilon_l   = Epsilonl_unw.eff_;
  Double_t UnwEpsilon_l_Any = EpsilonlAny_unw.eff_;
  Double_t UnwEpsilon_bHE = EpsilonbHE_unw.eff_; 
  Double_t UnwEpsilon_bHP = EpsilonbHP_unw.eff_; 

  Double_t Al             = Al_wgt.eff_;
  Double_t CHad           = Chad_wgt.eff_;
  Double_t CHad_Any       = ChadAny_wgt.eff_;

  Double_t Epsilon_l      = Epsilonl_wgt.eff_;
  Double_t Epsilon_l_Any  = EpsilonlAny_wgt.eff_;
  Double_t Epsilon_bHE    = EpsilonbHE_wgt.eff_; 
  Double_t Epsilon_bHP    = EpsilonbHP_wgt.eff_; 

  Double_t ErrUnwAl(-1.),         ErrAl(-1.);  
  Double_t ErrUnwCHad(-1.),       ErrCHad(-1.);
  Double_t ErrUnwCHad_Any(-1.),   ErrCHad_Any(-1.);
  Double_t ErrUnwEpsilonl(-1.),   ErrEpsilonl(-1.);
  Double_t ErrUnwEpsilonlAny(-1.),ErrEpsilonlAny(-1.);
  Double_t ErrUnwEpsilonbHE(-1.), ErrEpsilonbHE(-1.);
  Double_t ErrUnwEpsilonbHP(-1.), ErrEpsilonbHP(-1.);

  if(useClopperPearsonErrors_){
    
    ErrUnwAl         = Al_unw.errCP_;
    ErrUnwCHad       = Chad_unw.errCP_;
    ErrUnwCHad_Any   = ChadAny_unw.errCP_;
    			 
    ErrAl            = Al_wgt.errCP_;
    ErrCHad          = Chad_wgt.errCP_;
    ErrCHad_Any      = ChadAny_wgt.errCP_;

    ErrUnwEpsilonl    = Epsilonl_unw.errCP_;	 
    ErrUnwEpsilonlAny = EpsilonlAny_unw.errCP_; 
    ErrUnwEpsilonbHE  = EpsilonbHE_unw.errCP_; 
    ErrUnwEpsilonbHP  = EpsilonbHP_unw.errCP_; 

    ErrEpsilonl      = Epsilonl_wgt.errCP_;
    ErrEpsilonlAny   = EpsilonlAny_wgt.errCP_;	
    ErrEpsilonbHE    = EpsilonbHE_wgt.errCP_;
    ErrEpsilonbHP    = EpsilonbHP_wgt.errCP_;

  } else{

    ErrUnwAl       = Al_unw.errBin_;	 
    ErrUnwCHad     = Chad_unw.errBin_;	 
    ErrUnwCHad_Any = ChadAny_unw.errBin_;
    		     	 	   
    ErrAl             = Al_wgt.errImp_;	 
    ErrUnwEpsilonlAny = EpsilonlAny_unw.errBin_; 
    ErrCHad           = Chad_wgt.errImp_;	 
    ErrCHad_Any       = ChadAny_wgt.errImp_;
    
    ErrUnwEpsilonl   = Epsilonl_unw.errBin_;	 
    ErrUnwEpsilonbHE = EpsilonbHE_unw.errBin_; 
    ErrUnwEpsilonbHP = EpsilonbHP_unw.errBin_; 

    ErrEpsilonl      = Epsilonl_wgt.errImp_;	
    ErrEpsilonlAny   = EpsilonlAny_wgt.errImp_;
    ErrEpsilonbHE    = EpsilonbHE_wgt.errImp_;
    ErrEpsilonbHP    = EpsilonbHP_wgt.errImp_;

  }
 
  MakeNiceHistoStyle(h_GEN_mass_ZfromLeptons_);

  // TCanvas *c1 = new TCanvas("c1","c1",800,600);
  //   c1->SetFillColor(0); 
  //   c1->cd();
  //   h_GEN_mass_ZfromLeptons_->Draw();
  
  Double_t xsection_unw(0.);
  Double_t xsection(0.);
  Double_t xsecCorrFactor(1.);

  if(filename1.Contains("MADGRAPH")){
    xsecCorrFactor= 3048./(14584158.); 
  } else if (filename1.Contains("SHERPA")){
    xsecCorrFactor= 3048./(369549.); 
  } else {
    xsecCorrFactor=1.;
  }

  // std::cout<<"xsecCorrFactor: "<<   xsecCorrFactor <<std::endl;

  xsection_unw = N_genZkin_yes_genb_yes*xsecCorrFactor;
  xsection     = sumW_genZkin_yes_genb_yes*xsecCorrFactor;

  /*  
      std::cout<<"============== print final report =================="<<std::endl;
      std::cout<<" N(ll+b):        "<< GEN_genZkin_yes_genb_yes<<" sumW(ll+b):   "<<sumW_GEN_genZkin_yes_genb_yes<<std::endl;
      std::cout<<" sigma(Z+b)(unw):"<< xsection_unw <<" pb  sigma(Z+b): "<<xsection<<" pb"<<std::endl;
      std::cout<<" A_l(unw):       "<< UnwAl<<"+/-"<<ErrUnwAl<<    "\t A_l(wgt):      "<<Al<<"+/-"<<ErrAl<<std::endl;
      std::cout<<" C_hadron(unw):  "<< UnwCHad<<"+/-"<<ErrUnwCHad<<"\t C_hadron(wgt): "<<CHad<<"+/-"<<ErrCHad<<std::endl;
      std::cout<<"============== print report for Z+b/Z+j =============="<<std::endl;  
      std::cout<<" epsilon_recoB(unw):       "<<UnwEpsReco_B     <<"+/-"<<  ErrUnwEpsReco_B   <<"\t epsilon_recoB(wgt):      "<<  EpsReco_B    <<"+/-"<< ErrEpsReco_B <<std::endl;
      std::cout<<" epsilon_recoAny(unw):  "<< UnwEpsReco_Any  <<"+/-"<<  ErrUnwEpsReco_Any  <<"\t epsilon_recoAny(wgt): "<< EpsReco_Any  <<"+/-"<< ErrEpsReco_Any <<std::endl;
      std::cout<<"=========================================="<<std::endl;
  */

  if ( !partonLevel_ ) {
    cout<<"================================================================="<<endl;
    cout<<"Final acceptance (A_l) & hadron correction factor (C_had) "<<endl;
    cout<<"================================================================="<<endl;

    if(showUnweighted_){
      cout<<"Unweighted A_l             :"<< setprecision (3) << UnwAl*100.           <<" +/- "<< setprecision (2) << ErrUnwAl*100.         <<"%"<<endl;
      cout<<"Unweighted C_had (any)     :"<< setprecision (3) << UnwCHad_Any*100.     <<" +/- "<< setprecision (2) << ErrUnwCHad_Any*100.   <<"%"<<endl;
      cout<<"Unweighted C_had           :"<< setprecision (3) << UnwCHad*100.         <<" +/- "<< setprecision (2) << ErrUnwCHad*100.       <<"%"<<endl;
      cout<<"Unweighted Epsilon_l (any) :"<< setprecision (3) << UnwEpsilon_l_Any*100.<<" +/- "<< setprecision (2) << ErrUnwEpsilonlAny*100.<<"%"<<endl;
      cout<<"Unweighted Epsilon_l       :"<< setprecision (3) << UnwEpsilon_l*100.    <<" +/- "<< setprecision (2) << ErrUnwEpsilonl*100.   <<"%"<<endl;
      cout<<"Unweighted Epsilon_b (HE)  :"<< setprecision (3) << UnwEpsilon_bHE*100.  <<" +/- "<< setprecision (2) << ErrUnwEpsilonbHE*100. <<"%"<<endl;
      cout<<"Unweighted Epsilon_b (HP)  :"<< setprecision (3) << UnwEpsilon_bHP*100.  <<" +/- "<< setprecision (2) << ErrUnwEpsilonbHP*100. <<"%"<<endl;
      cout<<"================================================================="<<endl;
    }
   
    cout<<"A_l                        :"<< setprecision (3) << Al*100.             <<" +/- "<< setprecision (2) << ErrAl*100.         <<"%"<<endl;
    cout<<"C_had(any)                 :"<< setprecision (3) << CHad_Any*100.       <<" +/- "<< setprecision (2) << ErrCHad_Any*100.   <<"%"<<endl;
    cout<<"C_had                      :"<< setprecision (3) << CHad*100.           <<" +/- "<< setprecision (2) << ErrCHad*100.       <<"%"<<endl;
    cout<<"Epsilon_l (any)            :"<< setprecision (3) << Epsilon_l_Any*100.  <<" +/- "<< setprecision (2) << ErrEpsilonlAny*100.<<"%"<<endl; 
    cout<<"Epsilon_l                  :"<< setprecision (3) << Epsilon_l*100.      <<" +/- "<< setprecision (2) << ErrEpsilonl*100.   <<"%"<<endl; 
    cout<<"Epsilon_b (HE)             :"<< setprecision (3) << Epsilon_bHE*100.    <<" +/- "<< setprecision (2) << ErrEpsilonbHE*100. <<"%"<<endl;
    cout<<"Epsilon_b (HP)             :"<< setprecision (3) << Epsilon_bHP*100.    <<" +/- "<< setprecision (2) << ErrEpsilonbHP*100. <<"%"<<endl; 
  } else {
    cout<<"================================================================="<<endl;
    cout<<"Final acceptance (A_l) & hadron correction factor (C_part) "<<endl;
    cout<<"================================================================="<<endl;
    
    if(showUnweighted_){
      cout<<"Unweighted A_l             :"<< setprecision (3) << UnwAl*100.           <<" +/- "<< setprecision (2) << ErrUnwAl*100.         <<"%"<<endl;
      cout<<"Unweighted C_part (any)    :"<< setprecision (3) << UnwCHad_Any*100.     <<" +/- "<< setprecision (2) << ErrUnwCHad_Any*100.   <<"%"<<endl;
      cout<<"Unweighted C_part          :"<< setprecision (3) << UnwCHad*100.         <<" +/- "<< setprecision (2) << ErrUnwCHad*100.       <<"%"<<endl;
      cout<<"Unweighted Epsilon_l (any) :"<< setprecision (3) << UnwEpsilon_l_Any*100.<<" +/- "<< setprecision (2) << ErrUnwEpsilonlAny*100.<<"%"<<endl;
      cout<<"Unweighted Epsilon_l       :"<< setprecision (3) << UnwEpsilon_l*100.    <<" +/- "<< setprecision (2) << ErrUnwEpsilonl*100.   <<"%"<<endl;
      cout<<"Unweighted Epsilon_b (HE)  :"<< setprecision (3) << UnwEpsilon_bHE*100.  <<" +/- "<< setprecision (2) << ErrUnwEpsilonbHE*100. <<"%"<<endl;
      cout<<"Unweighted Epsilon_b (HP)  :"<< setprecision (3) << UnwEpsilon_bHP*100.  <<" +/- "<< setprecision (2) << ErrUnwEpsilonbHP*100. <<"%"<<endl;   
      cout<<"================================================================="<<endl;
    }
    cout<<"A_l                        :"<< setprecision (3) << Al*100.             <<" +/- "<< setprecision (2) << ErrAl*100.           <<"%"<<endl;
    cout<<"C_part (any)               :"<< setprecision (3) << CHad_Any*100.       <<" +/- "<< setprecision (2) << ErrCHad_Any*100.     <<"%"<<endl;
    cout<<"C_part                     :"<< setprecision (3) << CHad*100.           <<" +/- "<< setprecision (2) << ErrCHad*100.         <<"%"<<endl;
    cout<<"Epsilon_l (any)            :"<< setprecision (3) << Epsilon_l_Any*100.  <<" +/- "<< setprecision (2) << ErrEpsilonlAny*100.  <<"%"<<endl;  
    cout<<"Epsilon_l                  :"<< setprecision (3) << Epsilon_l*100.      <<" +/- "<< setprecision (2) << ErrEpsilonl*100.     <<"%"<<endl; 
    cout<<"Epsilon_b (HE)             :"<< setprecision (3) << Epsilon_bHE*100.    <<" +/- "<< setprecision (2) << ErrEpsilonbHE*100.   <<"%"<<endl;
    cout<<"Epsilon_b (HP)             :"<< setprecision (3) << Epsilon_bHP*100.    <<" +/- "<< setprecision (2) << ErrEpsilonbHP*100.   <<"%"<<endl; 
  } 
  cout<<"================================================================="<<endl;
}

