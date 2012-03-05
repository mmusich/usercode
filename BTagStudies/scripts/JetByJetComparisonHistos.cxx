#include "JetByJetComparisonHistos.h"
#include "JetInfo.h"

#include <TColor.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TObjArray.h>
#include <iostream>
#include <TStyle.h>
#include <TMath.h>

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
  setTDRStyle();
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
    
  //Histograms 2D (scatter plot for discriminants) ( for cross check)
  addHisto2D("h2ScatDiscrTCHE","Discr TCHE",200,-100.,100.,200,-100.,100.);
  addHisto2D("h2ScatDiscrTCHP","Discr TCHP",200,-100.,100.,200,-100.,100.);
  addHisto2D("h2ScatDiscrSSVHE","Discr SSVHE",100,-5.,5.,100,-5.,5.);
  addHisto2D("h2ScatDiscrSSVHP","Discr SSVHP",100,-5.,5.,100,-5.,5.);
  addHisto2D("h2ScatDiscrCSV","Discr CSV",100,-5.,5.,100,-5.,5.);
  addHisto2D("h2ScatDiscrJP","Discr JP",100,0.,4.,100,0.,4.);
  addHisto2D("h2ScatDiscrJBP","Discr JBP",100,0.,12.,100,0.,12.);

  //Histograms 2D (scatter plot for differences)  
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2","#Delta discr TCHE vs IP3d2",200,-5.,5.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3"," #Delta discr TCHP vs IP3d3",200,-5.,5.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2Error","#Delta discr TCHE vs IP3d2Error",200,-0.5,0.5,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3Error"," #Delta discr TCHP vs IP3d3Error",200,-200.,200.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsIP3d2/Error","#Delta discr TCHE vs IP3d2/Error",200,-200.,200.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsIP3d3/Error"," #Delta discr TCHP vs IP3d3/Error",200,-200.,200.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsEta","#Delta discr TCHE vs Eta",200,-4.,4.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsEta"," #Delta discr TCHP vs Eta",200,-4.,4.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsPhi","#Delta discr TCHE vs Phi",200,-4,4,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsPhi"," #Delta discr TCHP vs Phi",200,-4.,4.,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaPVz","#Delta discr TCHE vs DeltaPVz",200,-1,1,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaPVz","#Delta discr TCHP vs DeltaPVz",200,-1,1,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDelta3PV","#Delta discr TCHE vs Delta3PV",200,-1,1,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDelta3PV","#Delta discr TCHP vs Delta3PV",200,-1,1,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaXYPV","#Delta discr TCHE vs DeltaXYPV",200,-1,1,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaXYPV","#Delta discr TCHP vs DeltaXYPV",200,-1,1,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track","#Delta discr TCHE vs DeltaIP3dEta 2nd Track",200,-5,5,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track","#Delta discr TCHP vs DeltaIP3dEta 3rd Track",200,-5,5,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track","#Delta discr TCHE vs DeltaIP3dPhi 2nd Track",200,-5,5,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track","#Delta discr TCHP vs DeltaIP3dPhi 3rd Track",200,-5,5,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track","#Delta discr TCHE vs DeltaIP3dPt 2nd Track",200,-5000,5000,200,-200.,200.);
  addHisto2D("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track","#Delta discr TCHP vs DeltaIP3dPt 3rd Track",200,-5000,5000,200,-200.,200.);
  
    
  // Adding Profile plot   
//   addProfile("h2ScatDeltaDiscrTCHEvsIP3d2_pfx","#Delta discr TCHE vs IP3d2",200,-5.,5.,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHPvsIP3d3_pfx"," #Delta discr TCHP vs IP3d3",200,-5.,5.,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHEvsIP3d2Error_pfx","#Delta discr TCHE vs IP3d2Error",200,-0.5,0.5,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHPvsIP3d3Error_pfx"," #Delta discr TCHP vs IP3d3Error",200,-0.5,0.5,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHEvsIP3d2/Error_pfx","#Delta discr TCHE vs IP3d2/Error",200,-200.,200.,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHPvsIP3d3/Error_pfx"," #Delta discr TCHP vs IP3d3/Error",200,-200.,200.,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHEvsEta_pfx"," #Delta discr TCHE vs Eta",200,-4,4,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHPvsEta_pfx"," #Delta discr TCHP vs Eta",200,-4,4,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHEvsPhi_pfx"," #Delta discr TCHE vs Phi",200,-4,4,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHPvsPhi_pfx"," #Delta discr TCHP vs Phi",200,-4,4,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHEvsDeltaPVz_pfx", "#Delta discr TCHE vs DeltaPVz",200,-1,1,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHPvsDeltaPVz_pfx", "#Delta discr TCHP vs DeltaPVz",200,-1,1,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHEvsDelta3PV_pfx", "#Delta discr TCHE vs Delta3PV",200,-1,1,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHPvsDelta3PV_pfx", "#Delta discr TCHP vs Delta3PV",200,-1,1,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHEvsDeltaXYPV_pfx", "#Delta discr TCHE vs DeltaXYPV",200,-1,1,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHPvsDeltaXYPV_pfx", "#Delta discr TCHE vs DeltaXYPV",200,-1,1,-200,200);
//   addProfile("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track_pfx","#Delta discr TCHE vs DeltaIP3dEta 2nd Track",200,-5,5,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track_pfx","#Delta discr TCHP vs DeltaIP3dEta 3rd Track",200,-5,5,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track_pfx","#Delta discr TCHE vs DeltaIP3dPhi 2nd Track",200,-5,5,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track_pfx","#Delta discr TCHP vs DeltaIP3dPhi 3rd Track",200,-5,5,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track_pfx","#Delta discr TCHE vs DeltaIP3dPt 2nd Track",200,-5000,5000,-200.,200.);
//   addProfile("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track_pfx","#Delta discr TCHP vs DeltaIP3dPt 3rd Track",200,-5000,5000,-200.,200.);
  
    
}

    ////////////////////////////////////////////

void JetByJetComparisonHistos::setTDRStyle(){
    
    TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

    //    tdrStyle->SetGreyscale();
    
     static const UInt_t Number = 3;
     Double_t Red[Number]   = { 1.00, 0.50, 0.00};
     Double_t Green[Number] = { 1.00, 0.50, 0.00};
     Double_t Blue[Number]  = { 1.00, 0.50, 0.00};
     Double_t Stops[Number] = { 0.00, 0.90, 1.00};
     Int_t nb=10;
     TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);

        // For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(600); //Height of canvas
    tdrStyle->SetCanvasDefW(600); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //POsition on screen
    tdrStyle->SetCanvasDefY(0);
    
        // For the Pad:
    tdrStyle->SetPadBorderMode(0);
        // tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadColor(kWhite);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);
    
        // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);
    
        // For the histo:
        // tdrStyle->SetHistFillColor(1);
        // tdrStyle->SetHistFillStyle(0);
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
        // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
        // tdrStyle->SetNumberContours(Int_t number = 20);
    
    tdrStyle->SetEndErrorSize(2);
        //  tdrStyle->SetErrorMarker(20);
    tdrStyle->SetErrorX(0.);
    
    tdrStyle->SetMarkerStyle(20);
    
        //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);
    
        //For the date:
    tdrStyle->SetOptDate(0);
        // tdrStyle->SetDateX(Float_t x = 0.01);
        // tdrStyle->SetDateY(Float_t y = 0.01);
    
        // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat("emruo"); // To display the mean and RMS:   SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.025);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.4g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.1);
    tdrStyle->SetStatW(0.15);
        // tdrStyle->SetStatStyle(Style_t style = 1001);
        // tdrStyle->SetStatX(Float_t x = 0);
        // tdrStyle->SetStatY(Float_t y = 0);
    
        // Margins:
    tdrStyle->SetPadTopMargin(0.07);
    tdrStyle->SetPadBottomMargin(0.15);
    tdrStyle->SetPadLeftMargin(0.16);
    tdrStyle->SetPadRightMargin(0.02);
    
        // For the Global title:
    tdrStyle->SetOptTitle(0);
    tdrStyle->SetTitleFont(42);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.05);
        // tdrStyle->SetTitleH(0); // Set the height of the title box
        // tdrStyle->SetTitleW(0); // Set the width of the title box
        // tdrStyle->SetTitleX(0); // Set the position of the title box
        // tdrStyle->SetTitleY(0.985); // Set the position of the title box
        // tdrStyle->SetTitleStyle(Style_t style = 1001);
        // tdrStyle->SetTitleBorderSize(2);
    
        // For the axis titles:
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");
        // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
        // tdrStyle->SetTitleYSize(Float_t size = 0.02);
    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);
        // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
        // For the axis labels:
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.05, "XYZ");
    
        // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);
    
        // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);
    
        // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
        // tdrStyle->SetLineScalePS(Float_t scale = 3);
        // tdrStyle->SetLineStyleString(Int_t i, const char* text);
        // tdrStyle->SetHeaderPS(const char* header);
        // tdrStyle->SetTitlePS(const char* pstitle);
    
        // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
        // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
        // tdrStyle->SetPaintTextFormat(const char* format = "g");
        // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
        // tdrStyle->SetTimeOffset(Double_t toffset);
        // tdrStyle->SetHistMinimumZero(kTRUE);
    
    tdrStyle->cd();
    
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
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsIP3d2/Error"),ja.tche,jb.tche,(ja.trk[1].IP3d)/(ja.trk[1].IP3dError),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsIP3d3/Error"),ja.tchp,jb.tchp,(ja.trk[2].IP3d)/(ja.trk[2].IP3dError),ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsEta"),        ja.tche,jb.tche,ja.eta,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsEta"),        ja.tchp,jb.tchp,ja.eta,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsPhi"),        ja.tche,jb.tche,ja.phi,ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsPhi"),        ja.tchp,jb.tchp,ja.phi,ja.tchp-jb.tchp);
  fillTH(findTH2("h2ScatDeltaDiscrTCHEvsDeltaPVz"),        ja.tche,jb.tche,(ja.pv.PVz-jb.pv.PVz),ja.tche-jb.tche);
  fillTH(findTH2("h2ScatDeltaDiscrTCHPvsDeltaPVz"),        ja.tchp,jb.tchp,(ja.pv.PVz-jb.pv.PVz),ja.tchp-jb.tchp);
  
  float deltax = (ja.pv.PVx-jb.pv.PVx) ;
  float deltay = (ja.pv.PVy-jb.pv.PVy) ;
  float deltaz = (ja.pv.PVz-jb.pv.PVz) ;
 
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

//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsIP3d2"),      ja.tche,jb.tche,ja.trk[1].IP3d,ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsIP3d3"),      ja.tchp,jb.tchp,ja.trk[2].IP3d,ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsIP3d2Error"), ja.tche,jb.tche,ja.trk[1].IP3dError,ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsIP3d3Error"), ja.tchp,jb.tchp,ja.trk[2].IP3dError,ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsIP3d2/Error"),ja.tche,jb.tche,(ja.trk[1].IP3d)/(ja.trk[1].IP3dError),ja.tche-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsIP3d3/Error"),ja.tchp,jb.tchp,(ja.trk[2].IP3d)/(ja.trk[2].IP3dError),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsDelta3PV"),        ja.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsDelta3PV"),        ja.tchp,jb.tchp,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)+TMath::Power(deltaz,2)),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsDeltaXYPV"),        ja.tche,jb.tche,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsDeltaXYPV"),        ja.tchp,jb.tchp,TMath::Sqrt(TMath::Power(deltax,2)+TMath::Power(deltay,2)),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsDeltaPVz"),        ja.tche,jb.tche,(ja.pv.PVz-jb.pv.PVz),ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsDeltaPVz"),        ja.tchp,jb.tchp,(ja.pv.PVz-jb.pv.PVz),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsDeltaIP3dEta2Track"),        ja.tche,jb.tche,(ja.trk[1].eta-jb.trk[1].eta),ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsDeltaIP3dEta3Track"),        ja.tchp,jb.tchp,(ja.trk[2].eta-jb.trk[2].eta),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsDeltaIP3dPhi2Track"),        ja.tche,jb.tche,(ja.trk[1].phi-jb.trk[1].phi),ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsDeltaIP3dPhi3Track"),        ja.tchp,jb.tchp,(ja.trk[2].phi-jb.trk[2].phi),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsDeltaIP3dPt2Track"),        ja.tche,jb.tche,(ja.trk[1].pT-jb.trk[1].pT),ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsDeltaIP3dPt3Track"),        ja.tchp,jb.tchp,(ja.trk[2].pT-jb.trk[2].pT),ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsEta"),        ja.tche,jb.tche,ja.eta,ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsEta"),        ja.tchp,jb.tchp,ja.eta,ja.tchp-jb.tchp);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHEvsPhi"),        ja.tche,jb.tche,ja.phi,ja.tche-jb.tche);
//   fillTH(findTProfile("h2ScatDeltaDiscrTCHPvsPhi"),        ja.tchp,jb.tchp,ja.phi,ja.tchp-jb.tchp);
       
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

 // LOOP on the 2D histograms
  for(UInt_t h=0; h< h2vec.size(); h++){
    TCanvas *c = new TCanvas(h2vec[h]->GetName(),h2vec[h]->GetName(),800,600);
    c->cd();
    gPad->SetTopMargin(0.07);
    gPad->SetRightMargin(0.15);
    h2vec[h]->SetStats(kFALSE);

    h2vec[h]->Draw("colz");

    TProfile *hpfx_tmp = (TProfile*) h2vec[h]->ProfileX("_pfx",1,-1);
    hpfx_tmp->SetStats(kFALSE);
    hpfx_tmp->SetMarkerColor(2); 
    hpfx_tmp->SetMarkerSize(0.3); 
    hpfx_tmp->SetMarkerStyle(21); 
    hpfx_tmp->Draw("psame");

    //    c->Draw();
    
    fout->cd(dirname.Data()); 
    h2vec[h]->Write();
    hpfx_tmp->Write();
    TString canvName = h2vec[h]->GetName();
    c->SaveAs(dirname+canvName+".png");
      
  }
  
}



