#include "JetByJetComparisonHistos.h"
#include "JetInfo.h"

#include <TClass.h>
#include <TObjArray.h>
#include <iostream>

using namespace std;

const TString discriminators[5] = {"TCHE","TCHP","SSVHE","SSVHP","CSV"};     
const Float_t defaults[5] = {-100.,-100.,-1.,-1.,-1.};

JetByJetComparisonHistos::JetByJetComparisonHistos(const TString& s,TFile* fout) : dirname(s) {

  for(UInt_t disc=0; disc<5; disc++){
    defaultmap[discriminators[disc]] = defaults[disc];
  }

  //std::cout<<defaultmap.find("TCHE")->second <<" tagger "<< defaultmap.find("TCHE")->first<<std::endl;

  TObjArray *conditions_from_name = dirname.Tokenize("_"); 
  TString sa = conditions_from_name->At(0)->GetName();
  CondANotDef = sa.Contains("NotDefault");
  TString sb = conditions_from_name->At(1)->GetName();
  CondBNotDef = sb.Contains("NotDefault");
  h1vec.clear();
  h2vec.clear();
  hprofvec.clear();    
  fout->cd(); 
  fout->mkdir(dirname.Data()); 
  fout->cd(dirname.Data()); 
  addAllHistos();
}

void JetByJetComparisonHistos::addAllHistos() {
 
  //Histograms 1D 
  addHisto("hDeltaDiscrTCHE", "#Delta DiscrTCHE ;#Delta DiscrTCHE;jets",  1000,-200.,200.);
  addHisto("hDeltaDiscrTCHP", "#Delta DiscrTCHP ;#Delta DiscrTCHP;jets",  1000,-200.,200.);
  addHisto("hDeltaDiscrSSVHE","#Delta DiscrSSVHE ;#Delta DiscrSSVHE;jets",1000,-5.,5.    );
  addHisto("hDeltaDiscrSSVHP","#Delta DiscrSSVHP ;#Delta DiscrSSVHP;jets",1000,-5.,5.    );
  addHisto("hDeltaDiscrCSV",  "#Delta DiscrCSV ;#Delta DiscrCSV;jets",    1000,-5.,5.    );
 
  //Histograms 2D (scatter plot for differences)  
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2","#Delta discr TCHE vs IP3d2",200,-5.,5.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3"," #Delta discr TCHP vs IP3d3",200,-5.,5.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2Error","#Delta discr TCHE vs IP3d2Error",200,-0.5,0.5,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3Error"," #Delta discr TCHP vs IP3d3Error",200,-200.,200.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2/Error","#Delta discr TCHE vs IP3d2/Error",200,-200.,200.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3/Error"," #Delta discr TCHP vs IP3d3/Error",200,-200.,200.,200,-200.,200.);
  
  // Adding Profile plot   
  addProfile("ProfileDeltaDiscrTCHEvsIP3d2","#Delta discr TCHE vs IP3d2",200,-5.,5.,-200.,200.);
  addProfile("ProfileDeltaDiscrTCHPvsIP3d3"," #Delta discr TCHP vs IP3d3",200,-5.,5.,-200.,200.);
  addProfile("ProfileDeltaDiscrTCHEvsIP3d2Error","#Delta discr TCHE vs IP3d2Error",200,-0.5,0.5,-200.,200.);
  addProfile("ProfileDeltaDiscrTCHPvsIP3d3Error"," #Delta discr TCHP vs IP3d3Error",200,-0.5,0.5,-200.,200.);
  addProfile("ProfileDeltaDiscrTCHEvsIP3d2/Error","#Delta discr TCHE vs IP3d2/Error",200,-200.,200.,-200.,200.);
  addProfile("ProfileDeltaDiscrTCHPvsIP3d3/Error"," #Delta discr TCHP vs IP3d3/Error",200,-200.,200.,-200.,200.);
    
}

///////////////////////////////////  
void JetByJetComparisonHistos::addHisto(TString name, TString title,const int& nbins, const Float_t& min, const Float_t& max) {
  TH1F* h = new TH1F(name.Data(),title.Data(),nbins,min,max);       
    h->Sumw2();
    h1vec.push_back(h);
}

////////////////////////////////////
void JetByJetComparisonHistos::addHisto2D(TString name, TString title,const int& nbins,const Float_t& min,const Float_t& max,const int& nbinsy,const Float_t& miny, const Float_t& maxy) {
  TH2F* h2 = new TH2F(name.Data(),title.Data(),nbins,min,max,nbinsy,miny,maxy);       
    h2->Sumw2();
    h2vec.push_back(h2);
}

/////////////////////////////////////////////////////////////////////
void JetByJetComparisonHistos::addProfile(TString name, TString title,const int& nbins,const Float_t& min,const Float_t& max,const Float_t& miny, const Float_t& maxy) {
  TProfile* hprof = new TProfile(name.Data(),title.Data(),nbins,min,max,miny,maxy);       
    hprof->Sumw2();
    hprofvec.push_back(hprof);
}


/////////////////////////////////////////////////////////
TH1F* JetByJetComparisonHistos::findTH1(TString keyword){
    
    TH1F* htemp=h1vec[0];    
    UInt_t size = h1vec.size();
    for(UInt_t i=0; i< size; i++){
        if(((TString)h1vec[i]->GetName()).Contains(keyword)){
            htemp=h1vec[i];
            break;
        }
    }
    return htemp;
}

/////////////////////////////////////////////////////////////////////////////////////
TH2F* JetByJetComparisonHistos::findTH2(TString keyword){
    
    TH2F* htemp2=h2vec[0];    
    UInt_t size = h2vec.size();
    for(UInt_t i=0; i< size; i++){
        if(((TString)h2vec[i]->GetName()).Contains(keyword)){
            htemp2=h2vec[i];
            break;
        }
    }
    
    return htemp2;
    
}

//////////////////////////////////////////////////////////////////
TProfile* JetByJetComparisonHistos::findTProfile(TString keyword){
    
    TProfile* htempprofile=hprofvec[0];    
    UInt_t size = hprofvec.size();
    for(UInt_t i=0; i< size; i++){
        if(((TString)hprofvec[i]->GetName()).Contains(keyword)){
            htempprofile=hprofvec[i];
            break;
        }
    }
    
    return htempprofile;
    
}


////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparisonHistos::fillAllHistos(const JetInfo& ja, const JetInfo& jb, TFile* fout){
  fout->cd(dirname.Data());
  fillTH(findTH1("hDeltaDiscrTCHE"), ja.tche,  jb.tche);
  fillTH(findTH1("hDeltaDiscrTCHP"), ja.tchp,  jb.tchp);
  fillTH(findTH1("hDeltaDiscrSSVHE"),ja.ssvhe, jb.ssvhe);
  fillTH(findTH1("hDeltaDiscrSSVHP"),ja.ssvhp, jb.ssvhp);
  fillTH(findTH1("hDeltaDiscrCSV"),  ja.csv,   jb.csv);
  fillTH(findTH2("h2ScatDiscrTCHE"), ja.tche,jb.tche);
  fillTH(findTH2("h2ScatDiscrTCHP"), ja.tchp,jb.tchp);
  fillTH(findTH2("h2ScatDiscrCSV"), ja.csv,jb.csv);
  fillTH(findTH2("h2ScatDiscrJP"), ja.jp,jb.jp);
  fillTH(findTH2("h2ScatDiscrJBP"), ja.jbp,jb.jbp);
  fillTH(findTH2("h2ScatDiscrSSVHE"), ja.ssvhe,jb.ssvhp);
  fillTH(findTH2("h2ScatDiscrSSVHP"), ja.ssvhp,jb.ssvhp);
  
  /*  
      fillTH(findTH1("hDeltaDiscrTCHENotDefaultDefault"), ja.tche,  jb.tche);
      fillTH(findTH1("hDeltaDiscrTCHPNotDefaultDefault"), ja.tchp,  jb.tchp);
      fillTH(findTH1("hDeltaDiscrSSVHENotDefaultDefault"),ja.ssvhe, jb.ssvhe);
      fillTH(findTH1("hDeltaDiscrSSVHPNotDefaultDefault"),ja.ssvhp, jb.ssvhp);
      fillTH(findTH1("hDeltaDiscrCSVNotDefaultDefault"),  ja.csv,   jb.csv);
      fillTH(findTH1("hDeltaDiscrJPNotDefaultDefault"),  ja.jp,   jb.jp);
      fillTH(findTH1("hDeltaDiscrJBPNotDefaultDefault"),  ja.jbp,   jb.jbp);
      fillTH(findTH2("h2ScatDiscrTCHENotDefaultNotDefault"), ja.tche,jb.tche);
      fillTH(findTH2("h2ScatDiscrTCHPNotDefaultNotDefault"), ja.tchp,jb.tchp);
      fillTH(findTH2("h2ScatDiscrCSVNotDefaultNotDefault"), ja.csv,jb.csv);
      fillTH(findTH2("h2ScatDiscrJPNotDefaultNotDefault"), ja.jp,jb.jp);
      fillTH(findTH2("h2ScatDiscrJBPNotDefaultNotDefault"), ja.jbp,jb.jbp);
      fillTH(findTH2("h2ScatDiscrSSVHENotDefaultNotDefault"), ja.ssvhe,jb.ssvhp);
      fillTH(findTH2("h2ScatDiscrSSVHPNotDefaultNotDefault"), ja.ssvhp,jb.ssvhp)  
  */
  
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2"),ja.tche,jb.tche,ja.trk[1].IP3d);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3"),ja.tchp,jb.tchp,ja.trk[2].IP3d);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2Error"),ja.tche,jb.tche,ja.trk[1].IP3dError);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3Error"),ja.tchp,jb.tchp,ja.trk[2].IP3dError);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2/Error"),ja.tche,jb.tche,(ja.trk[1].IP3d)/(ja.trk[1].IP3dError));
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3/Error"),ja.tchp,jb.tchp,(ja.trk[2].IP3d)/(ja.trk[2].IP3dError));

  fillTH(findTProfile("ProfileDeltaDiscrTCHEvsIP3d2"),ja.tche,jb.tche,ja.trk[1].IP3d);
  fillTH(findTProfile("ProfileDeltaDiscrTCHPvsIP3d3"),ja.tchp,jb.tchp,ja.trk[2].IP3d);
  fillTH(findTProfile("ProfileDeltaDiscrTCHEvsIP3d2Error"),ja.tche,jb.tche,ja.trk[1].IP3dError);
  fillTH(findTProfile("ProfileDeltaDiscrTCHPvsIP3d3Error"),ja.tchp,jb.tchp,ja.trk[2].IP3dError);
  fillTH(findTProfile("ProfileDeltaDiscrTCHEvsIP3d2/Error"),ja.tche,jb.tche,(ja.trk[1].IP3d)/(ja.trk[1].IP3dError));
  fillTH(findTProfile("ProfileDeltaDiscrTCHPvsIP3d3/Error"),ja.tchp,jb.tchp,(ja.trk[2].IP3d)/(ja.trk[2].IP3dError));
   
}

////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparisonHistos::fillTH(TH1* p_h, float value1, float value2, float value3){ 
    
  Float_t default_(0.);
  TString tagger_;
  
  for(UInt_t disc=0; disc<5; disc++){                                    // opening loop on discriminant
    if(((TString)p_h->GetName()).Contains(discriminators[disc])) {       // match dicriminator/histogram name
      default_= defaultmap.find(discriminators[disc])->second;
      tagger_ = defaultmap.find(discriminators[disc])->first;
    }
  }                    
   
  //std::cout<<"tagger_="<<tagger_<<"default_="<<default_<<std::endl;
       
  if ( p_h->IsA()->InheritsFrom("TH1F") ) {                              // opening if type TH1F
    if (CondANotDef == true && CondBNotDef == false){                    // opening else if condition on case Default-Non Default
      if(value1 != default_ && value2 == default_){                   
	p_h->Fill(value1-value2);
      }                                                               
    } else if(CondANotDef == false && CondBNotDef == true){              // opening else if condition on case Not Default - Default
      if(value1 == default_ && value2 != default_){                    
	p_h->Fill(value1-value2);
      }                                                                
    }                                                                    // closing else if condition on case Not Default - Defaul
    else if(CondANotDef == true && CondBNotDef == true){                 // opening else if condition case NotDefault - Not Default
      if(value1 != default_ && value2 != default_){                      // opening if condition on values
	p_h->Fill(value1-value2);
      }                                                                  // closing if condition
    }                                                                    // closing else if condition case Not Default - Not Default 
  } else if ( p_h->IsA()->InheritsFrom("TH2F") ) {                       // opening if type TH2
    if (CondANotDef == true && CondBNotDef == false){                    // opening else if condition on case Default - Not Default 
      if(value1 != default_ && value2 == default_) {                  
	p_h->Fill(value1,value2);
      }                                                               
    } else if (CondANotDef == false && CondBNotDef == true) {            // opening else if condition on case Not Default - Default 
      if(value1 == default_ && value2 != default_) {                 
	p_h->Fill(value1,value2);
      }                                                                  // closing else if condition on case Not Default - Default
    }
    else if (CondANotDef == true && CondBNotDef == true){              
      if(value1 != default_ && value2 != default_){
	p_h->Fill(value1,value2);
      }
    }
  } else if ( p_h->IsA()->InheritsFrom("TProfile") ) {               //opening if type TProfile
    if (CondANotDef == false && CondBNotDef == true){               //opening else if condition on case Default-Non Default
      if(value1 == default_ && value2 != default_){                   
	p_h->Fill(value1-value2,value3);  
      }                          //opening if type TH1F 
    }
    else if (CondANotDef == true && CondBNotDef == false ) {        //opening else if condition on case Not Default - Default 
      if(value2 == default_ && value1 != default_) {                 
	p_h->Fill(value1-value2,value3);
      }                                                             // closing else if condition on case Not Default - Default
    }
    else if (CondANotDef == true && CondBNotDef == true){
      if(value1 != default_ && value2 != default_){
	p_h->Fill(value1-value2,value3);
      }
    }  
  } else {                                                               // opening if another TH* type
    cout << p_h->GetName() << " is unknown" << endl;  
  }                                                                      // closing if another type
}


