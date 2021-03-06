#ifndef JetByJetComparisonHistos_h
#define JetByJetComparisonHistos_h

#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <vector>
#include <map>

class JetInfo;

class JetByJetComparisonHistos{

 public:
  JetByJetComparisonHistos(const TString& s, TFile* fout);
    void addAllHistos();
    void addMigrPlotsHistos();
    void addHisto(TString name, TString title,const int& nbins, const Float_t& min, const Float_t& max);
    void AddHisto2D(std::string name, std::string title,TString firstCond,TString secondCond ,const int& nbins, const Float_t& min, const Float_t& max, const int& nbinsy, const Float_t& miny, const Float_t& maxy);
    void addHisto2D(TString name,TString title, const int& nbins,const Float_t& min,const Float_t& max,const int& nbinsy,const Float_t& miny,const Float_t& maxy);
    void addProfile(TString name, TString title,const int& nbins, const Float_t& min, const Float_t& max, const Float_t& miny, const Float_t& maxy);
    TH1F* findTH1(TString keyword);
    TH2F* findTH2(TString keyword);
    TProfile* findTProfile(TString keyword);
    void fillAllHistos(const JetInfo& ja, const JetInfo& jb, TFile* fout);
    void fillAllMigrationMatrices(const JetInfo& ja, const JetInfo& jb, TFile* fout);
    void fillTH(TH1* h, float value1, float value2, float valueX, float valueY=0);
    void fillMigrationMatrix(TH1* p_h, float value1, float value2);
    void drawNice2dHistos(TFile*);
    void cmsPrel(const double& intLumi);
    //   void setTDRStyle(TString palettename);
    //   void cmsPrel(const double& intLumi);

 private:
  std::map<TString,Float_t> defaultmap;
  TString dirname; 
  bool CondANotDef, CondBNotDef;
TString obj1name_, obj2name_;
  
  // store vector for dealing with canvas
  std::vector<TH1F*> h1vec;
  std::vector<TH2F*> h2vec;
  std::vector<TProfile*> hprofvec;

};

#endif
