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


  addAllHistos();

}

void JetByJetComparisonHistos::addAllHistos() {
 
  //Histograms 1D 
  addHisto("hDeltaDiscrTCHE", "#Delta DiscrTCHE ;#Delta DiscrTCHE;jets",  100,-100.,100.);
  addHisto("hDeltaDiscrTCHP", "#Delta DiscrTCHP ;#Delta DiscrTCHP;jets",  100,-100.,100.);
  addHisto("hDeltaDiscrSSVHE","#Delta DiscrSSVHE ;#Delta DiscrSSVHE;jets",100,-3.1,3.1    );
  addHisto("hDeltaDiscrSSVHP","#Delta DiscrSSVHP ;#Delta DiscrSSVHP;jets",100,-3.1,3.1    );
  addHisto("hDeltaDiscrCSV",  "#Delta DiscrCSV ;#Delta DiscrCSV;jets",    100,-1.1,1.1    );
  addHisto("hDeltaPT2TrackTCHE",  "#Delta PT 2 Track ;#Delta PT 2 Track;jets",    100,0.,150.    );
  addHisto("hDeltaPT3TrackTCHP",  "#Delta PT 3 Track;#Delta PT 3 Track;jets",    100,0.,150.    );
  addHisto("hDeltaEta2TrackTCHE",  "#Delta Eta 2 Track ;#Delta Eta 2 Track;jets",    100,-1.5,1.5    );
  addHisto("hDeltaEta3TrackTCHP",  "#Delta Eta 3 Track ;#Delta Eta 3 Track;jets",    100,-1.5,1.5    );
  addHisto("hDeltaPhi2TrackTCHE",  "#Delta Phi 2 Track ;#Delta Phi 2 Track;jets",    100,-1.5,1.5    );
  addHisto("hDeltaPhi3TrackTCHP",  "#Delta Phi 3 Track ;#Delta Phi 3 Track;jets",    100,-1.5,1.5    );
  addHisto("hDeltaXYPVTCHE", "#Delta XYPV ;#Delta XYPV;jets",    100, 0,0.03    );
  addHisto("hDeltaZPVTCHE",  "#Delta ZPV ;#Delta ZPV;jets",      100,-0.3,0.3    );
 

    
  //Histograms 2D (scatter plot for discriminants) ( for cross check)
  addHisto2D("h2ScatDiscrTCHE","Discr TCHE",obj1name_,obj2name_,100,-20.,20.,100,-20.,20.);
  addHisto2D("h2ScatDiscrTCHP","Discr TCHP",obj1name_,obj2name_,100,-20.,20.,100,-20.,20.);
  addHisto2D("h2ScatDiscrSSVHE","Discr SSVHE",obj1name_,obj2name_,100,-3.1,6.1,100,-3.1,6.1);
  addHisto2D("h2ScatDiscrSSVHP","Discr SSVHP",obj1name_,obj2name_,100,-3.1,6.1,100,-3.1,6.1);
  addHisto2D("h2ScatDiscrCSV","Discr CSV",obj1name_,obj2name_,100,-1.,1.4,100,-1.,1.4);
  addHisto2D("h2ScatDiscrJP","Discr JP",obj1name_,obj2name_,150,0.,4.,150,0.,4.);
  addHisto2D("h2ScatDiscrJBP","Discr JBP",obj1name_,obj2name_,150,0.,12.,150,0.,12.);

  //Histograms 2D (scatter plot for differences)  TCHE & TCHP DISCRIMINANTS
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2","#Delta discr TCHE vs IP3d2",80,-2.,6.,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3"," #Delta discr TCHP vs IP3d3",80,-2.,6,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2Error","#Delta discr TCHE vs IP3d2Error",80,0,0.2,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3Error"," #Delta discr TCHP vs IP3d3Error",80,0,0.2,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2overError","#Delta discr TCHE vs IP3d2overError",80,-100.,100.,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3overError"," #Delta discr TCHP vs IP3d3overError",80,-100.,100.,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsEta","#Delta discr TCHE vs Eta",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsEta"," #Delta discr TCHP vs Eta",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsPhi","#Delta discr TCHE vs Phi",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsPhi"," #Delta discr TCHP vs Phi",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaPVz","#Delta discr TCHE vs DeltaPVz",60,-0.03,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaPVz","#Delta discr TCHP vs DeltaPVz",60,-0.03,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDelta3PV","#Delta discr TCHE vs Delta3PV",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDelta3PV","#Delta discr TCHP vs Delta3PV",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaXYPV","#Delta discr TCHE vs DeltaXYPV",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaXYPV","#Delta discr TCHP vs DeltaXYPV",60,0.,0.03,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track","#Delta discr TCHE vs DeltaIP3dEta 2nd Track",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track","#Delta discr TCHP vs DeltaIP3dEta 3rd Track",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track","#Delta discr TCHE vs DeltaIP3dPhi 2nd Track",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track","#Delta discr TCHP vs DeltaIP3dPhi 3rd Track",30,-1.5,1.5,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track","#Delta discr TCHE vs DeltaIP3dPt 2nd Track",50,0,150,80,-20.,20.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track","#Delta discr TCHP vs DeltaIP3dPt 3rd Track",50,0,150,80,-20.,20.);
  
        // Histograms 2D (scatter plot for differences) SSVHE & SSVHP discriminants
    
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSV3dDistance","#Delta discr SSVHE vs SV3dDistance",150,-2.,14.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSV3dDistanceError"," #Delta discr SSVHE vs SV3dDistanceError",50,0.,2.5,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSV3dDistanceoverSV3dDistanceError"," #Delta discr SSVHE vs SV3dDistanceoverSV3dDistanceError",50,0.,100.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSVMass"," #Delta discr SSVHE vs SVMass",30,0.,15.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHEvsSVtotCharge"," #Delta discr SSVHE vs SVtotCharge",20,-10.,10.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSV3dDistance","#Delta discr SSVHP vs SV3dDistance",150,-2.,14.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSV3dDistanceError"," #Delta discr SSVHP vs SV3dDistanceError",50,0.,2.5,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSV3dDistanceoverSV3dDistanceError"," #Delta discr SSVHP vs SV3dDistanceoverSV3dDistanceError",100,0.,100.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSVMass"," #Delta discr SSVHP vs SVMass",30,0.,15.,100,-3.1,3.1);
  addHisto2D("h2ScatDeltaDiscrSSVHPvsSVtotCharge"," #Delta discr SSVHE vs SVtotCharge",20,-10.,10.,100,-3.1,3.1);
    
        //Histograms 2D (scatter plot for differences) CSV discriminant
    
    
 addHisto2D("h2ScatDeltaDiscrCSVvsSV3dDistance","#Delta discr CSV  vs SV3dDistance",150,-2.,14.,100,-1.1,1.1);
 addHisto2D("h2ScatDeltaDiscrCSVvsSV3dDistanceError"," #Delta discr CSV vs SV3dDistanceError",50,0.,2.5,100,-1.1,1.1);
 addHisto2D("h2ScatDeltaDiscrCSVvsSV3dDistanceoverSV3dDistanceError"," #Delta discr CSV vs SV3dDistanceoverSV3dDistanceError",50,0.,100.,100,-1.1,1.1);
 addHisto2D("h2ScatDeltaDiscrCSVvsDeltaPVz","#Delta discr CSV vs DeltaPVz",60,-0.03,0.03,100,-1.1,1.1);
 addHisto2D("h2ScatDeltaDiscrCSVvsDelta3PV","#Delta discr CSV vs Delta3PV",60,0,0.03,100,-1.1,1.1);
 addHisto2D("h2ScatDeltaDiscrCSVvsDeltaXYPV","#Delta discr CSV  vs DeltaXYPV",60,0.,0.03,100,-1.1,1.1);








    

    
    
    
    
  
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparisonHistos::addHisto2D(TString name,TString title,TString firstCond,TString secondCond ,const int& nbins, const Float_t& min, const Float_t& max, const int& nbinsy, const Float_t& miny, const Float_t& maxy)  {

  TString titlehisto;
  Form(titlehisto,"%s %s vs %s;%s - %s; %s - %s",title.Data(),firstCond.Data(),secondCond.Data(),title.Data(),firstCond.Data(),title.Data(),secondCond.Data());
  
  TH2F* h2 = new TH2F(name.Data(),titlehisto.Data(),nbins,min,max,nbinsy,miny,maxy);

 
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
  fillTH(findTH1("hDeltaDiscrTCHE"),  ja.tche,  jb.tche, ja.tche-jb.tche);
  fillTH(findTH1("hDeltaDiscrTCHP"),  ja.tchp,  jb.tchp, ja.tchp-jb.tchp);
  fillTH(findTH1("hDeltaDiscrSSVHE"), ja.ssvhe, jb.ssvhe,ja.ssvhe-jb.ssvhe);
  fillTH(findTH1("hDeltaDiscrSSVHP"), ja.ssvhp, jb.ssvhp,ja.ssvhp-jb.ssvhp);
  fillTH(findTH1("hDeltaDiscrCSV"),   ja.csv,   jb.csv,  ja.csv-jb.csv);
  fillTH(findTH1("hDeltaPT2TrackTCHE"),    ja.tche,   jb.tche,  (ja.trk[1].pT)-(jb.trk[1].pT));
  fillTH(findTH1("hDeltaPT3TrackTCHP"),    ja.tchp,   jb.tchp,  (ja.trk[2].pT)-(jb.trk[2].pT));
  fillTH(findTH1("hDeltaEta2TrackTCHE"),   ja.tche,   jb.tche,  (ja.trk[1].eta)-(jb.trk[1].eta));
  fillTH(findTH1("hDeltaEta3TrackTCHP"),   ja.tchp,   jb.tchp,  (ja.trk[2].eta)-(jb.trk[2].eta));
  fillTH(findTH1("hDeltaPhi2TrackTCHE"),   ja.tche,   jb.tche,  (ja.trk[1].phi)-(jb.trk[1].phi));
  fillTH(findTH1("hDeltaPhi3TrackTCHP"),   jb.tchp,   jb.tchp,  (ja.trk[2].phi)-(jb.trk[2].phi));


float deltax = (ja.pv.PVx-jb.pv.PVx) ;
  float deltay = (ja.pv.PVy-jb.pv.PVy) ;
  float deltaz = (ja.pv.PVz-jb.pv.PVz) ;


  fillTH(findTH1("hDeltaZPVTCHE"),   ja.tche,   jb.tche,  deltaz);


  fillTH(findTH1("hDeltaXYPVTCHE"),   jb.tche,   jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)));

    
    
  fillTH(findTH2("h2ScatDiscrTCHE"),  ja.tche,  jb.tche, ja.tche, jb.tche);
  fillTH(findTH2("h2ScatDiscrTCHP"),  ja.tchp,  jb.tchp, ja.tchp, jb.tchp);
  fillTH(findTH2("h2ScatDiscrCSV"),   ja.csv,   jb.csv,  ja.csv,  jb.csv);
  fillTH(findTH2("h2ScatDiscrJP"),    ja.jp,    jb.jp,   ja.jp,   jb.jp);
  fillTH(findTH2("h2ScatDiscrJBP"),   ja.jbp,   jb.jbp,  ja.jbp,  jb.jbp);
  fillTH(findTH2("h2ScatDiscrSSVHE"), ja.ssvhe, jb.ssvhe,ja.ssvhe,jb.ssvhe);
  fillTH(findTH2("h2ScatDiscrSSVHP"), ja.ssvhp, jb.ssvhp,ja.ssvhp,jb.ssvhp);
  
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
  
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2"),      ja.tche,jb.tche,ja.trk[1].IP3d,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3"),      ja.tchp,jb.tchp,ja.trk[2].IP3d,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2Error"), ja.tche,jb.tche,ja.trk[1].IP3dError,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3Error"), ja.tchp,jb.tchp,ja.trk[2].IP3dError,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2overError"),ja.tche,jb.tche,(ja.trk[1].IP3d)/(ja.trk[1].IP3dError),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3overError"),ja.tchp,jb.tchp,(ja.trk[2].IP3d)/(ja.trk[2].IP3dError),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsEta"),        ja.tche,jb.tche,ja.eta,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsEta"),        ja.tchp,jb.tchp,ja.eta,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsPhi"),        ja.tche,jb.tche,ja.phi,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsPhi"),        ja.tchp,jb.tchp,ja.phi,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaPVz"),        ja.tche,jb.tche,(ja.pv.PVz-jb.pv.PVz),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaPVz"),        ja.tchp,jb.tchp,(ja.pv.PVz-jb.pv.PVz),ja.tchp-jb.tchp);
  
  // float deltax = (ja.pv.PVx-jb.pv.PVx) ;
  //  float deltay = (ja.pv.PVy-jb.pv.PVy) ;
  //  float deltaz = (ja.pv.PVz-jb.pv.PVz) ;
 
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDelta3PV"),          ja.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDelta3PV"),          ja.tchp,jb.tchp,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaXYPV"),         ja.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaXYPV"),         ja.tchp,jb.tchp,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track"),ja.tche,jb.tche,(ja.trk[1].eta-jb.trk[1].eta),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track"),ja.tchp,jb.tchp,(ja.trk[2].eta-jb.trk[2].eta),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track"),ja.tche,jb.tche,(ja.trk[1].phi-jb.trk[1].phi),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track"),ja.tchp,jb.tchp,(ja.trk[2].phi-jb.trk[2].phi),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track"), ja.tche,jb.tche,(ja.trk[1].pT-jb.trk[1].pT),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track"), ja.tchp,jb.tchp,(ja.trk[2].pT-jb.trk[2].pT),ja.tchp-jb.tchp);           


  
    fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSV3dDistance"),ja.ssvhe,jb.ssvhe,ja.sv.SV3dDistance,ja.ssvhe-jb.ssvhe);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSV3dDistanceError"),ja.ssvhe,jb.ssvhe,ja.sv.SV3dDistanceError,ja.ssvhe-jb.ssvhe);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSV3dDistanceoverSV3dDistanceError"),ja.ssvhe,jb.ssvhe,(ja.sv.SV3dDistance)/(ja.sv.SV3dDistanceError),ja.ssvhe-jb.ssvhe);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSVMass"),ja.ssvhe,jb.ssvhe,ja.sv.SVMass,ja.ssvhe-jb.ssvhe);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHEvsSVtotCharge"),ja.ssvhe,jb.ssvhe,ja.sv.SVtotCharge,ja.ssvhe-jb.ssvhe);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSV3dDistance"),ja.ssvhp,jb.ssvhp,ja.sv.SV3dDistance,ja.ssvhp-jb.ssvhp);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSV3dDistanceError"),ja.ssvhp,jb.ssvhp,ja.sv.SV3dDistanceError,ja.ssvhp-jb.ssvhp);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSV3dDistanceoverSV3dDistanceError"),ja.ssvhp,jb.ssvhp,(ja.sv.SV3dDistance)/(ja.sv.SV3dDistanceError),ja.ssvhp-jb.ssvhp);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSVMass"),ja.ssvhp,jb.ssvhp,ja.sv.SVMass,ja.ssvhp-jb.ssvhp);
    fillTH(findTH2("h2ScatDeltaDiscrSSVHPvsSVtotCharge"),ja.ssvhp,jb.ssvhp,ja.sv.SVtotCharge,ja.ssvhp-jb.ssvhp);

    // filling CSV Plot

    fillTH(findTH2("h2ScatDeltaDiscrCSVvsSV3dDistance"),ja.csv,jb.csv,ja.sv.SV3dDistance,ja.csv-jb.csv);
    fillTH(findTH2("h2ScatDeltaDiscrCSVvsSV3dDistanceError"),ja.csv,jb.csv,ja.sv.SV3dDistanceError,ja.csv-jb.csv);
    fillTH(findTH2("h2ScatDeltaDiscrCSVvsSV3dDistanceoverSV3dDistanceError"),ja.csv,jb.csv,(ja.sv.SV3dDistance)/(ja.sv.SV3dDistanceError),ja.csv-jb.csv);
    fillTH(findTH2("h2ScatDeltaDiscrCSVvsDelta3PV"),          ja.csv,jb.csv,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.csv-jb.csv);
    fillTH(findTH2("h2ScatDeltaDiscrCSVvsDeltaXYPV"),         ja.csv,jb.csv,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.csv-jb.csv);
    fillTH(findTH2("h2ScatDeltaDiscrCSVvsDeltaPVz"),        ja.csv,jb.csv,(ja.pv.PVz-jb.pv.PVz),ja.csv-jb.csv);




    
    
    
       
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
   
    //  cmsPrel(60.);
    TString canvName = h1vec[h]->GetName()+dirname;
    c->SaveAs(canvName+".png");
  }
  
  // LOOP on the 2D histograms
  for(UInt_t h=0; h< h2vec.size(); h++){
    TCanvas *c = new TCanvas(h2vec[h]->GetName()+dirname,h2vec[h]->GetName(),800,600);
    c->cd();
    gPad->SetTopMargin(0.07);
    gPad->SetRightMargin(0.15);
    h2vec[h]->SetStats(kFALSE);
    h2vec[h]->Draw("colz");

    TProfile *hpfx_tmp = (TProfile*) h2vec[h]->ProfileX("_pfx",1,-1,"o");
    hpfx_tmp->SetStats(kFALSE);
    hpfx_tmp->SetMarkerColor(kBlack); 
    hpfx_tmp->SetMarkerSize(0.75); 
    hpfx_tmp->SetMarkerStyle(20); 
    hpfx_tmp->Draw("psame");

    //    c->Draw();
    
    fout->cd(dirname.Data()); 
    h2vec[h]->Write();
    hpfx_tmp->Write();
    
    //   cmsPrel(60.);
    
    TString canvName = h2vec[h]->GetName()+dirname;
    c->SaveAs(canvName+".png");
    if ( hpfx_tmp!=0) delete hpfx_tmp;
    if ( h2vec[h]!=0) delete h2vec[h];
  }
}



