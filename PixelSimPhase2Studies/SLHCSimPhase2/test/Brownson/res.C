#include <TFile.h>
#include <TTree.h>
//#include <TH1F.h>
//#include <TH1I.h>
#include <TCut.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <utility>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCut.h>
#include <TGraphAsymmErrors.h>
#include <TBranch.h>
#include <math.h>
#include <map>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include <TLine.h>
#include "TProfile.h"
#include "TF1.h"

using namespace std;
int main (){
void TH1_SQRT(TH1D*);
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./stdgrechitfullph1g_ntuple.root");
      if (!f) {
         f = new TFile("./stdgrechitfullph1g_ntuple.root");
      }
      TTree* tree = (TTree*)gDirectory->Get("PixelNtuple");
    TCanvas MyCanvas("MyCanvas", "New Geometry");
    gROOT->SetStyle("Plain");
    MyCanvas.UseCurrentStyle();
    gStyle->SetErrorX(0);

Long64_t nentries = tree->GetEntriesFast();
cout <<"\nThere are NEntries = "<<nentries<<"\n";

// *****************************************************************************
// Save some things here to copy and paste later...
// alphaPrime =  abs(abs(alpha)-(3.14159/2));
const double LorentzAngle = 0.3829;
const double HALF_PI = 1.57079632679489656;
// *****************************************************************************
// Microns of cm ... do everything in one or the other.
int units = 1; 	//     1 == cm
			// 10000 == Microns
if (units==1){cout <<"\nUsing units of cm\n";}
if (units==10000){cout <<"\nUsing units of microns\n";}
// *****************************************************************************
// Fill basic histograms
int MaxBXSpread = 3; // BPix spreadX==4+ will be lumped togeather
int MaxBYSpread = 9; 
int MaxFXSpread = 3;
int MaxFYSpread = 3;

std::vector<TH1F*> BPix_1x  ,BPix_1y;	// BPix Layer ==1
std::vector<TH1F*> BPix_234x,BPix_234y;	// BPix Layer ==2+
std::vector<TH1F*> FPix_x  ,FPix_y;     // FPix All disks

int color =0;
// *****************************************************************************
// Fill them
for (int i=1;i<=MaxBXSpread+1;i++){
	color++;    if (color==5) { color++; }	if (color==10) { color++; }
	std::ostringstream i_;i_<<i;
	string title = "BPix Layer1 spreadX=="+i_.str();
	string def   = "BPix_1x_"+i_.str();
	BPix_1x.push_back(new TH1F(def.c_str(),title.c_str(),80,units*-0.02,units*0.02));
	string what = "(pixel_recHit.hx-pixel_recHit.x)>>"+def;
	if (units == 10000) {what="10000*"+what;}
	string cuts = "pixel_recHit.subid==1&&pixel_recHit.layer==1&&pixel_recHit.spreadx=="+i_.str();
	if (i>MaxBXSpread){cuts = "pixel_recHit.subid==1&&pixel_recHit.layer==1&&pixel_recHit.spreadx>="+i_.str();}
	tree->Draw(what.c_str(),cuts.c_str(),0,nentries);
	BPix_1x[i-1]->SetStats(0);    BPix_1x[i-1]->SetLineColor(color);
	BPix_1x[i-1]->SetXTitle("RecHit X Residual");

        title = "BPix Layer2-4 spreadX=="+i_.str();
        def   = "BPix_234x_"+i_.str();
        BPix_234x.push_back(new TH1F(def.c_str(),title.c_str(),80,units*-0.02,units*0.02));
        what = "(pixel_recHit.hx-pixel_recHit.x)>>"+def;
        if (units == 10000) {what="10000*"+what;}
        cuts = "pixel_recHit.subid==1&&pixel_recHit.layer>1&&pixel_recHit.spreadx=="+i_.str();
        if (i>MaxBXSpread){cuts = "pixel_recHit.subid==1&&pixel_recHit.layer>1&&pixel_recHit.spreadx>="+i_.str();}
        tree->Draw(what.c_str(),cuts.c_str(),0,nentries);
        BPix_234x[i-1]->SetStats(0);    BPix_234x[i-1]->SetLineColor(color);
        BPix_234x[i-1]->SetXTitle("RecHit X Residual");


    } // End loop over MaxBXSpread
// Draw BPix_1x
std::vector<float> BPix_1x_Mean     , BPix_1x_Meanerr   , BPix_1x_RMS  , BPix_1x_RMSerr	 ;
std::vector<float> BPix_1x_FitMean  , BPix_1x_FitSigma  				 ;
std::vector<float> BPix_234x_Mean   , BPix_234x_Meanerr , BPix_234x_RMS, BPix_234x_RMSerr;
std::vector<float> BPix_234x_FitMean, BPix_234x_FitSigma				 ;

std::cout <<"\n\n";
for (int i=1;i<=MaxBXSpread+1;i++){
	BPix_1x_Mean.push_back(BPix_1x[i-1]->GetMean());
	BPix_1x_Meanerr.push_back(BPix_1x[i-1]->GetMeanError());
	BPix_1x_RMS.push_back(BPix_1x[i-1]->GetRMS());
	BPix_1x_RMSerr.push_back(BPix_1x[i-1]->GetRMSError());
	BPix_1x[i-1]->Fit("gaus","Q");
	BPix_1x[i-1]->GetFunction("gaus")->SetLineWidth(1);
	BPix_1x_FitMean.push_back(BPix_1x[i-1]->GetFunction("gaus")->GetParameter(1));
	BPix_1x_FitSigma.push_back(BPix_1x[i-1]->GetFunction("gaus")->GetParameter(2));
	std::cout <<"\nBPix X L1 Mean,RMS = \t\t else if ( sizex == "<<i<<" ) "	<<BPix_1x_Mean[i-1]	<<" , "
							<<BPix_1x_RMS[i-1]	<<"  "
							<<"\n\t\tFit Mean,Sigma \telse if ( sizex == "<<i<<" ) "
							<<BPix_1x_FitMean[i-1]  <<" , "
                                                	<<BPix_1x_FitSigma[i-1] ;
    } // End loop over MaxBXSpread
std::cout <<"\n\n";
for (int i=1;i<=MaxBXSpread+1;i++){
        BPix_234x_Mean.push_back(BPix_234x[i-1]->GetMean());
        BPix_234x_Meanerr.push_back(BPix_234x[i-1]->GetMeanError());
        BPix_234x_RMS.push_back(BPix_234x[i-1]->GetRMS());
        BPix_234x_RMSerr.push_back(BPix_234x[i-1]->GetRMSError());
        BPix_234x[i-1]->Fit("gaus","Q");
        BPix_234x[i-1]->GetFunction("gaus")->SetLineWidth(1);
        BPix_234x_FitMean.push_back(BPix_234x[i-1]->GetFunction("gaus")->GetParameter(1));
        BPix_234x_FitSigma.push_back(BPix_234x[i-1]->GetFunction("gaus")->GetParameter(2));
        std::cout <<"\nBPix X L234 Mean,RMS = \t\telse if ( sizex == "<<i<<" ) "       <<BPix_234x_Mean[i-1]     <<" , "
                                                        <<BPix_234x_RMS[i-1]      <<"  "
                                                        <<"\n\t\tFit Mean,Sigma \telse if ( sizex == "<<i<<" ) "
                                                        <<BPix_234x_FitMean[i-1]  <<" , "
                                                        <<BPix_234x_FitSigma[i-1] ;
    } // End loop over MaxBXSpread
//

cout << "\nTaking the ratio of Layer 1 to Layers 2-4\n";

for (int i=1;i<=MaxBXSpread+1;i++){
	cout <<setprecision(3)<< BPix_1x_FitSigma[i-1]/BPix_234x_FitSigma[i-1]<<"\t";
    } // End loop over MaxBXSpread
//

cout <<"\n\n";
MyCanvas.Clear();
MyCanvas.Divide(3,2);
TLegend  *legend_1xB  = new TLegend(.05,.05,.95,.95);
	{std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_1x_FitMean[0]<<" , #sigma="<<BPix_1x_FitSigma[0]; legend_1xB->AddEntry(BPix_1x[0],lab.str().c_str(),"pl");}
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_1x_FitMean[1]<<" , #sigma="<<BPix_1x_FitSigma[1]; legend_1xB->AddEntry(BPix_1x[1],lab.str().c_str(),"pl");}
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_1x_FitMean[2]<<" , #sigma="<<BPix_1x_FitSigma[2]; legend_1xB->AddEntry(BPix_1x[2],lab.str().c_str(),"pl");}
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_1x_FitMean[3]<<" , #sigma="<<BPix_1x_FitSigma[3]; legend_1xB->AddEntry(BPix_1x[3],lab.str().c_str(),"pl");}
for (int i=1;i<=MaxBXSpread+1;i++) {MyCanvas.cd(i);BPix_1x[i-1]->Draw();}
MyCanvas.cd(MaxBXSpread+2);
legend_1xB->Draw("Same");
MyCanvas.Print("BPix_Layer1_LocalX.eps");
MyCanvas.Clear();

MyCanvas.Divide(3,2);
TLegend  *legend_234xB  = new TLegend(.05,.05,.95,.95);
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_234x_FitMean[0]<<" , #sigma="<<BPix_234x_FitSigma[0]; legend_234xB->AddEntry(BPix_234x[0],lab.str().c_str(),"pl");}
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_234x_FitMean[1]<<" , #sigma="<<BPix_234x_FitSigma[1]; legend_234xB->AddEntry(BPix_234x[1],lab.str().c_str(),"pl");}
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_234x_FitMean[2]<<" , #sigma="<<BPix_234x_FitSigma[2]; legend_234xB->AddEntry(BPix_234x[2],lab.str().c_str(),"pl");}
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_234x_FitMean[3]<<" , #sigma="<<BPix_234x_FitSigma[3]; legend_234xB->AddEntry(BPix_234x[3],lab.str().c_str(),"pl");}
for (int i=1;i<=MaxBXSpread+1;i++) {MyCanvas.cd(i);BPix_234x[i-1]->Draw();}
MyCanvas.cd(MaxBXSpread+2);
legend_234xB->Draw("Same");
MyCanvas.Print("BPix_Layer234_LocalX.eps");
MyCanvas.Clear();
// *****************************************************************************
// Fill them
color = 0;
for (int i=1;i<=MaxBYSpread+1;i++){
        color++;    if (color==5) { color++; }  if (color==10) { color++; }
        std::ostringstream i_;i_<<i;
        string title = "BPix Layer1 spreadY=="+i_.str();
        string def   = "BPix_1y_"+i_.str();
        BPix_1y.push_back(new TH1F(def.c_str(),title.c_str(),80,units*-0.02,units*0.02));
        string what = "(pixel_recHit.hy-pixel_recHit.y)>>"+def;
        if (units == 10000) {what="10000*"+what;}
        string cuts = "pixel_recHit.subid==1&&pixel_recHit.layer==1&&pixel_recHit.spready=="+i_.str();
        if (i>MaxBYSpread){cuts = "pixel_recHit.subid==1&&pixel_recHit.layer==1&&pixel_recHit.spready>="+i_.str();}
        tree->Draw(what.c_str(),cuts.c_str(),0,nentries);
        BPix_1y[i-1]->SetStats(0);    BPix_1y[i-1]->SetLineColor(color);
        BPix_1y[i-1]->SetXTitle("RecHit y Residual");

        title = "BPix Layer2-4 spreadY=="+i_.str();
        def   = "BPix_234y_"+i_.str();
        BPix_234y.push_back(new TH1F(def.c_str(),title.c_str(),80,units*-0.02,units*0.02));
        what = "(pixel_recHit.hy-pixel_recHit.y)>>"+def;
        if (units == 10000) {what="10000*"+what;}
        cuts = "pixel_recHit.subid==1&&pixel_recHit.layer>1&&pixel_recHit.spready=="+i_.str();
        if (i>MaxBYSpread){cuts = "pixel_recHit.subid==1&&pixel_recHit.layer>1&&pixel_recHit.spready>="+i_.str();}
        tree->Draw(what.c_str(),cuts.c_str(),0,nentries);
        BPix_234y[i-1]->SetStats(0);    BPix_234y[i-1]->SetLineColor(color);
        BPix_234y[i-1]->SetXTitle("RecHit y Residual");
    } // End loop over MaxBYSpread
// Draw BPix_1y
std::vector<float> BPix_1y_Mean, BPix_1y_Meanerr, BPix_1y_RMS, BPix_1y_RMSerr;
std::vector<float> BPix_1y_FitMean,BPix_1y_FitSigma;
std::vector<float> BPix_234y_Mean, BPix_234y_Meanerr, BPix_234y_RMS, BPix_234y_RMSerr;
std::vector<float> BPix_234y_FitMean,BPix_234y_FitSigma;

std::cout <<"\n\n";
for (int i=1;i<=MaxBYSpread+1;i++){
        BPix_1y_Mean.push_back(BPix_1y[i-1]->GetMean());
        BPix_1y_Meanerr.push_back(BPix_1y[i-1]->GetMeanError());
        BPix_1y_RMS.push_back(BPix_1y[i-1]->GetRMS());
        BPix_1y_RMSerr.push_back(BPix_1y[i-1]->GetRMSError());
        BPix_1y[i-1]->Fit("gaus","Q");
        BPix_1y[i-1]->GetFunction("gaus")->SetLineWidth(1);
        BPix_1y_FitMean.push_back(BPix_1y[i-1]->GetFunction("gaus")->GetParameter(1));
        BPix_1y_FitSigma.push_back(BPix_1y[i-1]->GetFunction("gaus")->GetParameter(2));
        std::cout <<"\nBPix Y L1 Mean,RMS from plot = \t\telse if ( sizey == "<<i<<" ) "       <<BPix_1y_Mean[i-1]     <<" , "
                                                        <<BPix_1y_RMS[i-1]      <<"  "
                                                        <<"\n\t\tFit Mean,Sigma \telse if ( sizey == "<<i<<" ) "
                                                        <<BPix_1y_FitMean[i-1]  <<" , "
                                                        <<BPix_1y_FitSigma[i-1] ;
    } // End loop over MaxBYSpread
std::cout <<"\n\n";
for (int i=1;i<=MaxBYSpread+1;i++){
        BPix_234y_Mean.push_back(BPix_234y[i-1]->GetMean());
        BPix_234y_Meanerr.push_back(BPix_234y[i-1]->GetMeanError());
        BPix_234y_RMS.push_back(BPix_234y[i-1]->GetRMS());
        BPix_234y_RMSerr.push_back(BPix_234y[i-1]->GetRMSError());
        BPix_234y[i-1]->Fit("gaus","Q");
        BPix_234y[i-1]->GetFunction("gaus")->SetLineWidth(1);
        BPix_234y_FitMean.push_back(BPix_234y[i-1]->GetFunction("gaus")->GetParameter(1));
        BPix_234y_FitSigma.push_back(BPix_234y[i-1]->GetFunction("gaus")->GetParameter(2));
        std::cout <<"\nBPix Y L234 Mean,RMS from plot = \t\telse if ( sizey == "<<i<<" ) "       <<BPix_234y_Mean[i-1]     <<" , "
                                                        <<BPix_234y_RMS[i-1]      <<"  "
                                                        <<"\n\t\tFit Mean,Sigma \telse if ( sizey == "<<i<<" ) "
                                                        <<BPix_234y_FitMean[i-1]  <<" , "
                                                        <<BPix_234y_FitSigma[i-1] ;
    } // End loop over MaxBYSpread
//
cout << "\nTaking the ratio of Layer 1 to Layers 2-4\n";
for (int i=1;i<=MaxBYSpread+1;i++){
        cout <<setprecision(3)<< BPix_1y_FitSigma[i-1]/BPix_234y_FitSigma[i-1]<<"\t";
    } // End loop over MaxBXSpread
//
cout <<"\n\n";
MyCanvas.Clear();
MyCanvas.Divide(4,3);
TLegend  *legend_1yB  = new TLegend(.05,.05,.95,.95);
for (int i=0;i<=MaxBYSpread ;i++){
        std::ostringstream i_;i_<<i+1;
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_1y_FitMean[i]<<" , #sigma="<<BPix_1y_FitSigma[i]; legend_1yB->AddEntry(BPix_1y[i],lab.str().c_str(),"pl");}
}
for (int i=1;i<=MaxBYSpread+1;i++) {MyCanvas.cd(i);BPix_1y[i-1]->Draw();}
MyCanvas.cd(MaxBYSpread+2);
legend_1yB->Draw("Same");
MyCanvas.Print("BPix_Layer1_LocalY.eps");
MyCanvas.Clear();
//
cout <<"\n\n";
MyCanvas.Clear();
MyCanvas.Divide(4,3);
TLegend  *legend_234yB  = new TLegend(.05,.05,.95,.95);
for (int i=0;i<=MaxBYSpread ;i++){
        std::ostringstream i_;i_<<i+1;
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<BPix_234y_FitMean[i]<<" , #sigma="<<BPix_234y_FitSigma[i]; legend_234yB->AddEntry(BPix_234y[i],lab.str().c_str(),"pl");}
}
for (int i=1;i<=MaxBYSpread+1;i++) {MyCanvas.cd(i);BPix_234y[i-1]->Draw();}
MyCanvas.cd(MaxBYSpread+2);
legend_234yB->Draw("Same");
MyCanvas.Print("BPix_Layer234_LocalY.eps");
MyCanvas.Clear();

// *****************************************************************************
// Fill them
color = 0;
for (int i=1;i<=MaxFXSpread+1;i++){
        color++;    if (color==5) { color++; }  if (color==10) { color++; }
        std::ostringstream i_;i_<<i;
        string title = "FPix spreadX=="+i_.str();
        string def   = "FPix_x_"+i_.str();
        FPix_x.push_back(new TH1F(def.c_str(),title.c_str(),80,units*-0.02,units*0.02));
        string what = "(pixel_recHit.hx-pixel_recHit.x)>>"+def;
        if (units == 10000) {what="10000*"+what;}
        string cuts = "pixel_recHit.subid==2&&pixel_recHit.spreadx=="+i_.str();
        if (i>MaxFXSpread){cuts = "pixel_recHit.subid==2&&pixel_recHit.spreadx>="+i_.str();}
        tree->Draw(what.c_str(),cuts.c_str(),0,nentries);
        FPix_x[i-1]->SetStats(0);    FPix_x[i-1]->SetLineColor(color);
        FPix_x[i-1]->SetXTitle("RecHit X Residual");
    } // End loop over MaxFXSpread
// Output FPix_x
std::vector<float> FPix_x_Mean, FPix_x_Meanerr, FPix_x_RMS, FPix_x_RMSerr;
std::vector<float> FPix_x_FitMean,FPix_x_FitSigma;
std::cout <<"\n\n";
for (int i=1;i<=MaxFXSpread+1;i++){
        FPix_x_Mean.push_back(FPix_x[i-1]->GetMean());
        FPix_x_Meanerr.push_back(FPix_x[i-1]->GetMeanError());
        FPix_x_RMS.push_back(FPix_x[i-1]->GetRMS());
        FPix_x_RMSerr.push_back(FPix_x[i-1]->GetRMSError());
        FPix_x[i-1]->Fit("gaus","Q");
        FPix_x[i-1]->GetFunction("gaus")->SetLineWidth(1);
        FPix_x_FitMean.push_back(FPix_x[i-1]->GetFunction("gaus")->GetParameter(1));
        FPix_x_FitSigma.push_back(FPix_x[i-1]->GetFunction("gaus")->GetParameter(2));
        std::cout <<"\nFPix X Mean,RMS from plot= \t\telse if ( sizex == "<<i<<" ) "          <<FPix_x_Mean[i-1]     <<" , "
                                                        <<FPix_x_RMS[i-1]      <<"  "
                                                        <<"\n\t\tFit Mean,Sigma \telse if ( sizex == "<<i<<" ) "
                                                        <<FPix_x_FitMean[i-1]  <<" , "
                                                        <<FPix_x_FitSigma[i-1] ;
//        if (fabs(FPix_x_Mean[i-1])>fabs(2*FPix_x_Meanerr[i-1])) {std::cout << "\nThe X Mean for FPix spreadX=="<<i<<" is over 2 Standard Deviations of the Mean away from zero\n";}
    } // End loop over MaxBYSpread
//
cout <<"\n\n";
MyCanvas.Clear();
MyCanvas.Divide(3,2);
TLegend  *legend_fxB  = new TLegend(.05,.05,.95,.95);
for (int i=0;i<=MaxFXSpread ;i++){
        std::ostringstream i_;i_<<i+1;
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<FPix_x_FitMean[i]<<" , #sigma="<<FPix_x_FitSigma[i]; legend_fxB->AddEntry(FPix_x[i],lab.str().c_str(),"pl");}
}
for (int i=1;i<=MaxFXSpread+1;i++) {MyCanvas.cd(i);FPix_x[i-1]->Draw();}
MyCanvas.cd(MaxFXSpread+2);
legend_fxB->Draw("Same");
MyCanvas.Print("FPix_LocalX.eps");
MyCanvas.Clear();

// *****************************************************************************
// Fill them
color = 0;
for (int i=1;i<=MaxFYSpread+1;i++){
        color++;    if (color==5) { color++; }  if (color==10) { color++; }
        std::ostringstream i_;i_<<i;
        string title = "FPix spreadY=="+i_.str();
        string def   = "FPix_y_"+i_.str();
        FPix_y.push_back(new TH1F(def.c_str(),title.c_str(),80,units*-0.02,units*0.02));
        string what = "(pixel_recHit.hy-pixel_recHit.y)>>"+def;
        if (units == 10000) {what="10000*"+what;}
        string cuts = "pixel_recHit.subid==2&&pixel_recHit.spready=="+i_.str();
        if (i>MaxFYSpread){cuts = "pixel_recHit.subid==2&&pixel_recHit.spready>="+i_.str();}
        tree->Draw(what.c_str(),cuts.c_str(),0,nentries);
        FPix_y[i-1]->SetStats(0);    FPix_y[i-1]->SetLineColor(color);
        FPix_y[i-1]->SetXTitle("RecHit Y Residual");
    } // End loop over MaxFYSpread
// Output FPix_y
std::vector<float> FPix_y_Mean, FPix_y_Meanerr, FPix_y_RMS, FPix_y_RMSerr;
std::vector<float> FPix_y_FitMean, FPix_y_FitSigma;
std::cout <<"\n\n";
for (int i=1;i<=MaxFYSpread+1;i++){
        FPix_y_Mean.push_back(FPix_y[i-1]->GetMean());
        FPix_y_Meanerr.push_back(FPix_y[i-1]->GetMeanError());
        FPix_y_RMS.push_back(FPix_y[i-1]->GetRMS());
        FPix_y_RMSerr.push_back(FPix_y[i-1]->GetRMSError());
        FPix_y[i-1]->Fit("gaus","Q");
        FPix_y[i-1]->GetFunction("gaus")->SetLineWidth(1);
        FPix_y_FitMean.push_back(FPix_y[i-1]->GetFunction("gaus")->GetParameter(1));
        FPix_y_FitSigma.push_back(FPix_y[i-1]->GetFunction("gaus")->GetParameter(2));
        std::cout <<"\nFPix Y Mean,RMS from plot = \t\telse if ( sizey == "<<i<<" ) "          <<FPix_y_Mean[i-1]     <<" , "
                                                        <<FPix_y_RMS[i-1]      <<"  "
                                                        <<"\n\t\tFit Mean,Sigma \telse if ( sizey == "<<i<<" ) "
                                                        <<FPix_y_FitMean[i-1]  <<" , "
                                                        <<FPix_y_FitSigma[i-1] ;
//        if (fabs(FPix_y_Mean[i-1])>fabs(2*FPix_y_Meanerr[i-1])) {std::cout << "\nThe Y Mean for FPix spreadX=="<<i<<" is over 2 Standard Deviations of the Mean away from zero\n";}
    } // End loop over MaxBYSpread
//
cout <<"\n\n";
MyCanvas.Clear();
MyCanvas.Divide(3,2);
TLegend  *legend_fyB  = new TLegend(.05,.05,.95,.95);
for (int i=0;i<=MaxFYSpread ;i++){
        std::ostringstream i_;i_<<i+1;
        {std::ostringstream lab; lab << "#mu="<<setprecision(4)<<FPix_y_FitMean[i]<<" , #sigma="<<FPix_y_FitSigma[i]; legend_fyB->AddEntry(FPix_y[i],lab.str().c_str(),"pl");}
}
for (int i=1;i<=MaxFYSpread+1;i++) {MyCanvas.cd(i);FPix_y[i-1]->Draw();}
MyCanvas.cd(MaxFYSpread+2);
legend_fyB->Draw("Same");
MyCanvas.Print("FPix_LocalY.eps");
MyCanvas.Clear();
// *****************************************************************************
/*
    xerr_barrel_l1_= cms.untracked.vdouble(0.00115, 0.00120, 0.00088),
    xerr_barrel_l1_def_=cms.untracked.double(0.01030),
    yerr_barrel_l1_= cms.untracked.vdouble(0.00375,0.00230,0.00250,0.00250,0.00230,0.00230,0.00210,0.00210,0.00240),
    yerr_barrel_l1_def_=cms.untracked.double(0.00210),
    xerr_barrel_ln_= cms.untracked.vdouble(0.00115, 0.00120, 0.00088),
    xerr_barrel_ln_def_=cms.untracked.double(0.01030),
    yerr_barrel_ln_= cms.untracked.vdouble(0.00375,0.00230,0.00250,0.00250,0.00230,0.00230,0.00210,0.00210,0.00240),
    yerr_barrel_ln_def_=cms.untracked.double(0.00210),
    xerr_endcap_= cms.untracked.vdouble(0.0020, 0.0020),
    xerr_endcap_def_=cms.untracked.double(0.0020),
    yerr_endcap_= cms.untracked.double(0.00210),
    yerr_endcap_def_=cms.untracked.double(0.00075)

*/
std::cout 	<< "\n\nOutputing the sigmas from the gaus fits to be inserted into the code."
		<< "\n    xerr_barrel_ln_ = cms.untracked.vdouble("	
                                                <<BPix_234x_FitSigma[0]/units<<","
						<<BPix_234x_FitSigma[1]/units<<","
                                                <<BPix_234x_FitSigma[2]/units
						<<"),"
		<< "\n    xerr_barrel_ln_def_ = cms.untracked.double("<<BPix_234x_FitSigma[3]<<"),"
                << "\n    yerr_barrel_ln_ = cms.untracked.vdouble("	
                                                <<BPix_234y_FitSigma[0]/units<<","
                                                <<BPix_234y_FitSigma[1]/units<<","
                                                <<BPix_234y_FitSigma[2]/units<<","
                                                <<BPix_234y_FitSigma[3]/units<<","
                                                <<BPix_234y_FitSigma[4]/units<<","
                                                <<BPix_234y_FitSigma[5]/units<<","
                                                <<BPix_234y_FitSigma[6]/units<<","
                                                <<BPix_234y_FitSigma[7]/units<<","
                                                <<BPix_234y_FitSigma[8]/units
						<<"),"
                << "\n    yerr_barrel_ln_def_ = cms.untracked.double("<<BPix_234y_FitSigma[9]<<"),"
                << "\n    xerr_endcap_ = cms.untracked.vdouble("	
                                                <<FPix_x_FitSigma[0]/units<<","
						<<FPix_x_FitSigma[1]/units<<","
                                                <<FPix_x_FitSigma[2]/units
						<<"),"
                << "\n    xerr_endcap_def_ = cms.untracked.double("<<FPix_x_FitSigma[3]<<"),"
                << "\n    yerr_endcap_ = cms.untracked.vdouble("	
                                                <<FPix_y_FitSigma[0]/units<<","
                                                <<FPix_y_FitSigma[1]/units<<","
                                                <<FPix_y_FitSigma[2]/units
						<<"),"
                << "\n    yerr_endcap_def_ = cms.untracked.double("<<FPix_y_FitSigma[3]<<"),"
		<< "\n\n"
                << "\n    xerr_barrel_l1_ = cms.untracked.vdouble("   
                                                <<BPix_1x_FitSigma[0]/units<<","
                                                <<BPix_1x_FitSigma[1]/units<<","
                                                <<BPix_1x_FitSigma[2]/units
                                                <<"),"
                << "\n    xerr_barrel_l1_def_ = cms.untracked.double("<<BPix_1x_FitSigma[3]<<"),"
                << "\n    yerr_barrel_l1_ = cms.untracked.vdouble("   
      <<BPix_1y_FitSigma[0]/units<<","
                                                <<BPix_1y_FitSigma[1]/units<<","
                                                <<BPix_1y_FitSigma[2]/units<<","
                                                <<BPix_1y_FitSigma[3]/units<<","
                                                <<BPix_1y_FitSigma[4]/units<<","
                                                <<BPix_1y_FitSigma[5]/units<<","
                                                <<BPix_1y_FitSigma[6]/units<<","
                                                <<BPix_1y_FitSigma[7]/units<<","
                                                <<BPix_1y_FitSigma[8]/units
                                                <<"),"
                << "\n    yerr_barrel_l1_def_ = cms.untracked.double("<<BPix_1y_FitSigma[9]<<")"
                << "\n\n";












std::cout << "\nThe program ran\n\a";
return 0;
}


void TH1_SQRT(TH1D *Hist)
{ // Takes a TH1F and returns the sqrt of everything
  Int_t nbinsx = Hist->GetNbinsX();
  for (int i = 0 ; i < nbinsx ; i++){
      double e1 = Hist->GetBinError(i);
      double b1 = Hist->GetBinContent(i);
	//cout << "\nLooking at "<<Hist->GetTitle()<<" bin "<<i<<", with "<<b1<<" +/- "<<e1;
      Hist->SetBinContent(i,sqrt(Hist->GetBinContent(i)));
      Hist->SetBinError(i,(sqrt(b1+e1)-sqrt(b1-e1))/2.0 );
  }
  return;
}
