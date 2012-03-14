#include "JetByJetComparisonHistos.h"
#include "JetInfo.h"
#include <TColor.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TObjArray.h>
#include <iostream>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <TPaveStats.h>

using namespace std;

const TString discriminators[5] = {"TCHE","TCHP","SSVHE","SSVHP","CSV"};     
const Float_t defaults[5] = {-100.,-100.,-1.,-1.,-1.};

JetByJetComparisonHistos::JetByJetComparisonHistos(const TString& s,TFile* fout) : dirname(s) {

  for(UInt_t disc=0; disc<5; disc++){
    defaultmap[discriminators[disc]] = defaults[disc];
  }

  //  std::cout<<defaultmap.find("TCHE")->second <<" tagger "<< defaultmap.find("TCHE")->first<<std::endl;

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
  TString foutname = fout->GetName();

  TObjArray *obj_from_name= foutname.Tokenize("_");
  obj1name_ = obj_from_name->At(1)->GetName();
  obj2name_ = obj_from_name->At(3)->GetName();
  obj2name_.ReplaceAll(".root","");

  /*
    std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
    std::cout<<"obj name1: "<<obj1name_<<std::endl;
    std::cout<<"obj name2: "<<obj2name_<<std::endl;
    std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  */
  
  addAllHistos();

}

void JetByJetComparisonHistos::addAllHistos() {
 
  //Histograms 1D 
  addHisto("hDeltaDiscrTCHE","#Delta DiscrTCHE;#Delta D_{TCHE};jets",100,-100.,100.);
  addHisto("hDeltaDiscrTCHP","#Delta DiscrTCHP;#Delta D_{TCHP};jets",100,-100.,100.);
  addHisto("hDeltaDiscrSSVHE","#Delta DiscrSSVHE;#Delta D_{SSVHE};jets",100,-3.1,3.1);
  addHisto("hDeltaDiscrSSVHP","#Delta DiscrSSVHP;#Delta D_{SSVHP};jets",100,-3.1,3.1);
  addHisto("hDeltaDiscrCSV","#Delta DiscrCSV ;#Delta D_{CSV};jets", 100,-1.1,1.1);
  addHisto("hDeltaPT2TrackTCHE","#Delta PT 2 Track;#Deltap_{T} 2^{nd} track (GeV);jets",100,0.,150.);
  addHisto("hDeltaPT3TrackTCHP","#Delta PT 3 Track;#Deltap_{T} 3^{rd} track (GeV);jets",100,0.,150.);
  addHisto("hDeltaEta2TrackTCHE","#Delta Eta 2 Track ;#Delta#eta 2^{nd} track;jets",100,-1.5,1.5);
  addHisto("hDeltaEta3TrackTCHP","#Delta Eta 3 Track ;#Delta#eta 3^{rd} track;jets",100,-1.5,1.5);
  addHisto("hDeltaPhi2TrackTCHE","#Delta Phi 2 Track ;#Delta#phi 2^{nd} track;jets",100,-1.5,1.5);
  addHisto("hDeltaPhi3TrackTCHP","#Delta Phi 3 Track ;#Delta#phi 3^{rd} track;jets",100,-1.5,1.5);
  addHisto("hDeltaXYPVTCHE","#Delta XYPV ;#DeltaL_{xy}(PV) (cm);vertices",100,0,0.03);
  addHisto("hDeltaZPVTCHE","#Delta ZPV ;#DeltaZ_{PV} (cm);vertices",100,-0.3,0.3);

  //Histograms 2D (scatter plot for discriminants) ( for cross check)
  AddHisto2D("h2ScatDiscrTCHE","D_{TCHE}",obj1name_,obj2name_,100,-20.,20.,100,-20.,20.);
  AddHisto2D("h2ScatDiscrTCHP","D_{TCHP}",obj1name_,obj2name_,100,-20.,20.,100,-20.,20.);
  AddHisto2D("h2ScatDiscrSSVHE","D_{SSVHE}",obj1name_,obj2name_,100,-3.1,6.1,100,-3.1,6.1);
  AddHisto2D("h2ScatDiscrSSVHP","D_{SSVHP}",obj1name_,obj2name_,100,-3.1,6.1,100,-3.1,6.1);
  AddHisto2D("h2ScatDiscrCSV","D_{CSV}",obj1name_,obj2name_,100,-1.,4.,100,-1.,4.);
  AddHisto2D("h2ScatDiscrJP","D_{JP}",obj1name_,obj2name_,150,0.,4.,150,0.,4.);
  AddHisto2D("h2ScatDiscrJBP","D_{JBP}",obj1name_,obj2name_,150,0.,12.,150,0.,12.);

  // Histograms 2D (scatter plot for differences)  TCHE & TCHP DISCRIMINANTS
  // histogram of differences vs some third variable are always booked as:
  // valueA - valueB vs valueB ======> valueB is ALWAYS the reference!

  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2","#Delta discr TCHE vs IP3d2;IP_{3D}(2^{nd} track) (cm); #DeltaD_{TCHE}",80,-2.,6.,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3"," #Delta discr TCHP vs IP3d3;IP_{3D}(3^{rd} track) (cm); #DeltaD_{TCHP}",80,-2.,6,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2Error","#Delta discr TCHE vs IP3d2Error;#sigma_{IP3D}(2^{nd} track) (cm); #DeltaD_{TCHE}",80,0.,0.2,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3Error"," #Delta discr TCHP vs IP3d3Error;#sigma_{IP3D}(3^{rd} track) (cm); #DeltaD_{TCHP}",80,0.2,0.2,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2overError","#Delta discr TCHE vs IP3d2/Error;SIP_{3D}(2^{nd} track); #DeltaD_{TCHE}",80,-100.,100.,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3overError"," #Delta discr TCHP vs IP3d3/Error;SIP_{3D}(3^{rd} track); #DeltaD_{TCHP}",80,-100.,100.,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsEta","#Delta discr TCHE vs Eta;jet #eta; #DeltaD_{TCHE}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsEta"," #Delta discr TCHP vs Eta;jet #eta; #DeltaD_{TCHP}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsPhi","#Delta discr TCHE vs Phi;jet #phi; #DeltaD_{TCHE}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsPhi"," #Delta discr TCHP vs Phi;jet #phi; #DeltaD_{TCHP}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaPVz","#Delta discr TCHE vs DeltaPVz; #DeltaZ_{PV} (cm); #DeltaD_{TCHE}",60,-0.03,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaPVz","#Delta discr TCHP vs DeltaPVz; #DeltaZ_{PV} (cm); #DeltaD_{TCHP}",60,-0.03,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDelta3PV","#Delta discr TCHE vs Delta3PV; #DeltaL_{3D}(PV) (cm); #DeltaD_{TCHE}",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDelta3PV","#Delta discr TCHP vs Delta3PV; #DeltaL_{3D}(PV) (cm); #DeltaD_{TCHP}",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaXYPV","#Delta discr TCHE vs DeltaXYPV; #DeltaL_{xy}(PV) (cm); #DeltaD_{TCHE}",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaXYPV","#Delta discr TCHP vs DeltaXYPV; #DeltaL_{xy}(PV) (cm); #DeltaD_{TCHP}",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track","#Delta discr TCHE vs DeltaIP3dEta 2nd Track; #Delta#eta(2^{nd} track); #DeltaD_{TCHE}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track","#Delta discr TCHP vs DeltaIP3dEta 3rd Track; #Delta#eta(3^{rd} track); #DeltaD_{TCHP}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track","#Delta discr TCHE vs DeltaIP3dPhi 2nd Track; #Delta#phi(2^{nd} track); #DeltaD_{TCHE}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track","#Delta discr TCHP vs DeltaIP3dPhi 3rd Track; #Delta#phi(3^{rd} track); #DeltaD_{TCHP}",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track","#Delta discr TCHE vs DeltaIP3dPt 2nd Track; #Deltap_{T}(2^{nd} track) (GeV); #DeltaD_{TCHE}",50,0,150,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track","#Delta discr TCHP vs DeltaIP3dPt 3rd Track; #Deltap_{T}(3^{rd} track) (GeV); #DeltaD_{TCHP}",50,0,150,80,-20.,20.);
  
  // Histograms 2D (scatter plot for differences) SSVHE & SSVHP discriminants  
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSV3dDistance","#Delta discr SSVHE vs SV3dDistance;L_{3D}(SV) (cm); #DeltaD_{SSVHE}",150,-2.,14.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSV3dDistanceError","#Delta discr SSVHE vs SV3dDistanceError;#sigma_{L_{3D}}(SV) (cm);#DeltaD_{SSVHE}",50,0.,2.5,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSV3dDistanceoverSV3dDistanceError","#Delta discr SSVHE vs SV3dDistanceoverSV3dDistanceError;SIP_{L_{3D}}(SV);#DeltaD_{SSVHE}",50,0.,100.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSVMass","#Delta discr SSVHE vs SVMass;SV mass (GeV);#DeltaD_{SSVHE}",100,0.,15.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSVtotCharge","#Delta discr SSVHE vs SVtotCharge; SV tot charge;#DeltaD_{SSVHE}",21,-10.5,10.5,100,-3.1,3.1);
  
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSV3dDistance","#Delta discr SSVHP vs SV3dDistance;L_{3D}(SV) (cm); #DeltaD_{SSVHP}",150,-2.,14.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSV3dDistanceError","#Delta discr SSVHP vs SV3dDistanceError;#sigma_{L_{3D}}(SV) (cm);#DeltaD_{SSVHP}",50,0.,2.5,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSV3dDistanceoverSV3dDistanceError","#Delta discr SSVHP vs SV3dDistanceoverSV3dDistanceError;SIP_{L_{3D}}(SV);#DeltaD_{SSVHP}",100,0.,100.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSVMass","#Delta discr SSVHP vs SVMass;SV mass (GeV);#DeltaD_{SSVHP}",30,0.,15.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSVtotCharge","#Delta discr SSVHE vs SVtotCharge; SV tot charge;#DeltaD_{SSVHE}",21,-10.5,10.5,100,-3.1,3.1);
    
  //Histograms 2D (scatter plot for differences) CSV discriminant
  addHisto2D("h2ScatDeltaDiscrCSVvsSV3dDistance","#Delta discr CSV vs SV3dDistance;L_{3D}(SV) (cm);#DeltaD_{CSV}",150,-2.,14.,100,-1.1,1.1);
  addHisto2D("h2ScatDeltaDiscrCSVvsSV3dDistanceError","#Delta discr CSV vs SV3dDistanceError;#sigma_{L_{3D}}(SV) (cm);#DeltaD_{CVS}",50,0.,2.5,100,-1.1,1.1);
  addHisto2D("h2ScatDeltaDiscrCSVvsSV3dDistanceoverSV3dDistanceError","#Delta discr CSV vs SV3dDistanceoverSV3dDistanceError;SIP_{L_{3D}}(SV);#DeltaD_{SSVHP}",50,0.,100.,100,-1.1,1.1);
  addHisto2D("h2ScatDeltaDiscrCSVvsDeltaPVz","#Delta discr CSV vs DeltaPVz;#DeltaZ_{PV} (cm);#DeltaD_{CVS}",60,-0.03,0.03,100,-1.1,1.1);
  addHisto2D("h2ScatDeltaDiscrCSVvsDelta3PV","#Delta discr CSV vs Delta3PV;#DeltaL_{3D}(PV) (cm);#DeltaD_{CSV}",60,0,0.03,100,-1.1,1.1);
  addHisto2D("h2ScatDeltaDiscrCSVvsDeltaXYPV","#Delta discr CSV vs DeltaXYPV;#DeltaL_{XY}(PV) (cm);#DeltaD_{CSV}",60,0.,0.03,100,-1.1,1.1);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparisonHistos::AddHisto2D(std::string name, std::string title,TString firstCond,TString secondCond ,const int& nbins, const Float_t& min, const Float_t& max, const int& nbinsy, const Float_t& miny, const Float_t& maxy)  {
        
  char titlehisto[80];
  sprintf(titlehisto,"%s %s vs %s;%s - %s; %s -  %s",title.c_str(),firstCond.Data(),secondCond.Data(),title.c_str(),firstCond.Data(),title.c_str(),secondCond.Data());

  std::cout<<"titlehisto:"<<titlehisto<<std::endl;

  TH2F* h2 = new TH2F(name.c_str(),titlehisto,nbins,min,max,nbinsy,miny,maxy);

  h2->Sumw2();
  h2vec.push_back(h2);
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

  fillTH(findTH1("hDeltaDiscrTCHE"),ja.tche,jb.tche,ja.tche-jb.tche);
  fillTH(findTH1("hDeltaDiscrTCHP"),ja.tchp,jb.tchp,ja.tchp-jb.tchp);
  fillTH(findTH1("hDeltaDiscrSSVHE"),ja.ssvhe, jb.ssvhe,ja.ssvhe-jb.ssvhe);
  fillTH(findTH1("hDeltaDiscrSSVHP"),ja.ssvhp, jb.ssvhp,ja.ssvhp-jb.ssvhp);
  fillTH(findTH1("hDeltaDiscrCSV"),ja.csv,jb.csv,ja.csv-jb.csv);
  fillTH(findTH1("hDeltaPT2TrackTCHE"),ja.tche,jb.tche,(ja.trk[1].pT)-(jb.trk[1].pT));
  fillTH(findTH1("hDeltaPT3TrackTCHP"),ja.tchp,jb.tchp,(ja.trk[2].pT)-(jb.trk[2].pT));
  fillTH(findTH1("hDeltaEta2TrackTCHE"),ja.tche,jb.tche,(ja.trk[1].eta)-(jb.trk[1].eta));
  fillTH(findTH1("hDeltaEta3TrackTCHP"),ja.tchp,jb.tchp,(ja.trk[2].eta)-(jb.trk[2].eta));
  fillTH(findTH1("hDeltaPhi2TrackTCHE"),ja.tche,jb.tche,(ja.trk[1].phi)-(jb.trk[1].phi));
  fillTH(findTH1("hDeltaPhi3TrackTCHP"),jb.tchp,jb.tchp,(ja.trk[2].phi)-(jb.trk[2].phi));
     
  float deltax = (ja.pv.PVx-jb.pv.PVx);
  float deltay = (ja.pv.PVy-jb.pv.PVy);
  float deltaz = (ja.pv.PVz-jb.pv.PVz);

  fillTH(findTH1("hDeltaZPVTCHE"),ja.tche,jb.tche,deltaz);
  fillTH(findTH1("hDeltaXYPVTCHE"),jb.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)));

  fillTH(findTH2("h2ScatDiscrTCHE"),ja.tche,jb.tche,ja.tche,jb.tche);
  fillTH(findTH2("h2ScatDiscrTCHP"),ja.tchp,jb.tchp,ja.tchp,jb.tchp);
  fillTH(findTH2("h2ScatDiscrCSV"),ja.csv,jb.csv,ja.csv,jb.csv);
  fillTH(findTH2("h2ScatDiscrJP"),ja.jp,jb.jp,ja.jp,jb.jp);
  fillTH(findTH2("h2ScatDiscrJBP"),ja.jbp,jb.jbp,ja.jbp,jb.jbp);
  fillTH(findTH2("h2ScatDiscrSSVHE"),ja.ssvhe,jb.ssvhe,ja.ssvhe,jb.ssvhe);
  fillTH(findTH2("h2ScatDiscrSSVHP"),ja.ssvhp,jb.ssvhp,ja.ssvhp,jb.ssvhp);
  
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2"),ja.tche,jb.tche,jb.trk[1].IP3d,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3"),ja.tchp,jb.tchp,jb.trk[2].IP3d,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2Error"),ja.tche,jb.tche,jb.trk[1].IP3dError,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3Error"),ja.tchp,jb.tchp,jb.trk[2].IP3dError,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2overError"),ja.tche,jb.tche,(jb.trk[1].IP3d)/(jb.trk[1].IP3dError),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3overError"),ja.tchp,jb.tchp,(jb.trk[2].IP3d)/(jb.trk[2].IP3dError),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsEta"),ja.tche,jb.tche,jb.eta,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsEta"),ja.tchp,jb.tchp,jb.eta,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsPhi"),ja.tche,jb.tche,jb.phi,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsPhi"),ja.tchp,jb.tchp,jb.phi,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaPVz"),ja.tche,jb.tche,(ja.pv.PVz-jb.pv.PVz),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaPVz"),ja.tchp,jb.tchp,(ja.pv.PVz-jb.pv.PVz),ja.tchp-jb.tchp);
  
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDelta3PV"),ja.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDelta3PV"),ja.tchp,jb.tchp,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaXYPV"),ja.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaXYPV"),ja.tchp,jb.tchp,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track"),ja.tche,jb.tche,(ja.trk[1].eta-jb.trk[1].eta),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track"),ja.tchp,jb.tchp,(ja.trk[2].eta-jb.trk[2].eta),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track"),ja.tche,jb.tche,(ja.trk[1].phi-jb.trk[1].phi),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track"),ja.tchp,jb.tchp,(ja.trk[2].phi-jb.trk[2].phi),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track"),ja.tche,jb.tche,(ja.trk[1].pT-jb.trk[1].pT),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track"),ja.tchp,jb.tchp,(ja.trk[2].pT-jb.trk[2].pT),ja.tchp-jb.tchp);           

  fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSV3dDistance"),ja.ssvhe,jb.ssvhe,jb.sv.SV3dDistance,ja.ssvhe-jb.ssvhe);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSV3dDistanceError"),ja.ssvhe,jb.ssvhe,jb.sv.SV3dDistanceError,ja.ssvhe-jb.ssvhe);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSV3dDistanceoverSV3dDistanceError"),ja.ssvhe,jb.ssvhe,(jb.sv.SV3dDistance)/(jb.sv.SV3dDistanceError),ja.ssvhe-jb.ssvhe);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSVMass"),ja.ssvhe,jb.ssvhe,jb.sv.SVMass,ja.ssvhe-jb.ssvhe);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSVtotCharge"),ja.ssvhe,jb.ssvhe,jb.sv.SVtotCharge,ja.ssvhe-jb.ssvhe);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSV3dDistance"),ja.ssvhp,jb.ssvhp,jb.sv.SV3dDistance,ja.ssvhp-jb.ssvhp);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSV3dDistanceError"),ja.ssvhp,jb.ssvhp,jb.sv.SV3dDistanceError,ja.ssvhp-jb.ssvhp);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSV3dDistanceoverSV3dDistanceError"),ja.ssvhp,jb.ssvhp,(jb.sv.SV3dDistance)/(jb.sv.SV3dDistanceError),ja.ssvhp-jb.ssvhp);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSVMass"),ja.ssvhp,jb.ssvhp,jb.sv.SVMass,ja.ssvhp-jb.ssvhp);
  fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSVtotCharge"),ja.ssvhp,jb.ssvhp,jb.sv.SVtotCharge,ja.ssvhp-jb.ssvhp);
  
  // filling CSV plots
  fillTH(findTH2("h2ScatDeltaDiscrCSVvsSV3dDistance"),ja.csv,jb.csv,jb.sv.SV3dDistance,ja.csv-jb.csv);
  fillTH(findTH2("h2ScatDeltaDiscrCSVvsSV3dDistanceError"),ja.csv,jb.csv,jb.sv.SV3dDistanceError,ja.csv-jb.csv);
  fillTH(findTH2("h2ScatDeltaDiscrCSVvsSV3dDistanceoverSV3dDistanceError"),ja.csv,jb.csv,(jb.sv.SV3dDistance)/(jb.sv.SV3dDistanceError),ja.csv-jb.csv);
  fillTH(findTH2("h2ScatDeltaDiscrCSVvsDelta3PV"),ja.csv,jb.csv,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.csv-jb.csv);
  fillTH(findTH2("h2ScatDeltaDiscrCSVvsDeltaXYPV"),ja.csv,jb.csv,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.csv-jb.csv);
  fillTH(findTH2("h2ScatDeltaDiscrCSVvsDeltaPVz"),ja.csv,jb.csv,(ja.pv.PVz-jb.pv.PVz),ja.csv-jb.csv);
     
}

////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparisonHistos::fillTH(TH1* p_h, float value1, float value2, float valueX, float valueY){ 
  // value1, value2 needed for...
  // valueX
  // valueY ...    
  Float_t default_(0.);
  TString tagger_;
  
  for(UInt_t disc=0; disc<5; disc++){                                    // opening loop on discriminant
    if(((TString)p_h->GetName()).Contains(discriminators[disc])) {       // match dicriminator/histogram name
      default_= defaultmap.find(discriminators[disc])->second;
      tagger_ = defaultmap.find(discriminators[disc])->first;
    }
  }                    
  
  // std::cout<<"tagger_="<<tagger_<<"default_="<<default_<<std::endl;
  
  if ( p_h->IsA()->InheritsFrom("TH1F") ) {                              // opening if type TH1F
    if (CondANotDef == true && CondBNotDef == false){                    // opening else if condition on case Default-Non Default
      if(value1 != default_ && value2 == default_){                   
	p_h->Fill(valueX);
      }                                                               
    } else if(CondANotDef == false && CondBNotDef == true){              // opening else if condition on case Not Default - Default
      if(value1 == default_ && value2 != default_){                    
	p_h->Fill(valueX);
      }                                                                
    }                                                                    // closing else if condition on case Not Default - Defaul
    else if(CondANotDef == true && CondBNotDef == true){                 // opening else if condition case NotDefault - Not Default
      if(value1 != default_ && value2 != default_){                      // opening if condition on values
	p_h->Fill(valueX);
      }                                                                  // closing if condition
    }                                                                    // closing else if condition case Not Default - Not Default 
  } else if ( p_h->IsA()->InheritsFrom("TH2F") ) {                       // opening if type TH2
    if (CondANotDef == true && CondBNotDef == false){                    // opening else if condition on case Default - Not Default 
      if(value1 != default_ && value2 == default_) {                  
	p_h->Fill(valueX,valueY);
      }                                                               
    } else if (CondANotDef == false && CondBNotDef == true) {            // opening else if condition on case Not Default - Default 
      if(value1 == default_ && value2 != default_) {                 
	p_h->Fill(valueY,valueX);
      }                                                                  // closing else if condition on case Not Default - Default
    }
    else if (CondANotDef == true && CondBNotDef == true){              
      if(value1 != default_ && value2 != default_){
	p_h->Fill(valueX,valueY);
      }
    }
  } else if ( p_h->IsA()->InheritsFrom("TProfile") ) {               //opening if type TProfile
    if (CondANotDef == false && CondBNotDef == true){               //opening else if condition on case Default-Non Default
      if(value1 == default_ && value2 != default_){                   
	p_h->Fill(valueX,valueY);  
      }                          //opening if type TH1F 
    }
    else if (CondANotDef == true && CondBNotDef == false ) {        //opening else if condition on case Not Default - Default 
      if(value2 == default_ && value1 != default_) {                 
	p_h->Fill(valueX,valueY);
      }                                                             // closing else if condition on case Not Default - Default
    }
    else if (CondANotDef == true && CondBNotDef == true){
      if(value1 != default_ && value2 != default_){
	p_h->Fill(valueX,valueY);
      }
    }  
  } else {                                                               // opening if another TH* type
    cout << p_h->GetName() << " is unknown" << endl;  
  }                                                                      // closing if another type
}

void JetByJetComparisonHistos::drawNice2dHistos(TFile* fout)
{
  
  // LOOP on the 1D histograms
  UInt_t nOfHistos = h1vec.size();
  TObject    *statObj[nOfHistos];
  TPaveStats *stats[nOfHistos];
  
  for(UInt_t h=0; h<h1vec.size(); h++){
    TCanvas *c = new TCanvas(h1vec[h]->GetName()+dirname,h1vec[h]->GetName(),600,600);
    c->cd()->SetLogy();
    h1vec[h]->Draw();
    c->Draw();
    
    statObj[h] = h1vec[h]->GetListOfFunctions()->FindObject("stats");
    stats[h]= static_cast<TPaveStats*>(statObj[h]);
    stats[h]->SetFillColor(10);
    stats[h]->SetLineWidth(1);
    stats[h]->SetShadowColor(0);
    stats[h]->SetTextFont(42);
    stats[h]->SetTextSize(0.025);
    //stats[h]->SetLineColor(LineColors[h]);
    //stats[h]->SetTextColor(LineColors[h]);
    stats[h]->SetX1NDC(0.75);
    stats[h]->SetY1NDC(0.72);
    stats[h]->SetX2NDC(0.97);
    stats[h]->SetY2NDC(0.92);
    stats[h]->Draw("same"); 
   
    cmsPrel(60.);
    TString canvName = h1vec[h]->GetName()+dirname;
    c->SaveAs(canvName+".png");
  }
  
  // LOOP on the 2D histograms
  for(UInt_t h=0; h< h2vec.size(); h++){
    TCanvas *c = new TCanvas(h2vec[h]->GetName()+dirname,h2vec[h]->GetName(),800,600);
    c->cd();
    gPad->SetTopMargin(0.08);
    gPad->SetRightMargin(0.15);
    h2vec[h]->SetStats(kFALSE);
    h2vec[h]->Draw("colz");
    // c->Draw();

    TProfile *hpfx_tmp = (TProfile*) h2vec[h]->ProfileX("_pfx",1,-1,"o");
    hpfx_tmp->SetStats(kFALSE);
    hpfx_tmp->SetMarkerColor(kRed); 
    hpfx_tmp->SetMarkerSize(1); 
    hpfx_tmp->SetMarkerStyle(20); 
    hpfx_tmp->Draw("psame");

    // c->Draw();
    fout->cd(dirname.Data()); 
    h2vec[h]->Write();
    hpfx_tmp->Write();
    cmsPrel(60.);
    
    TString canvName = h2vec[h]->GetName()+dirname;
    c->cd()->SetLogz();
    c->SaveAs(canvName+".png");
    if ( hpfx_tmp!=0) delete hpfx_tmp;
    if ( h2vec[h]!=0) delete h2vec[h];
  }
}

void JetByJetComparisonHistos::cmsPrel(const double& intLumi) {
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22

  latex->SetTextAlign(13);
  //  latex->DrawLatex(0.12, 0.99, Form("CMS Preliminary 2011,     #sqrt{s} = 7 TeV,  L = %.2g pb^{ -1}",intLumi));
  latex->DrawLatex(0.2,0.99,"CMS Simulation 2011, top-top #sqrt{s} = 7 TeV");


  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.045);
  latex2->SetTextFont(72); //22
  latex2->DrawLatex(0.20,0.85,obj1name_+" vs "+obj2name_);

}


