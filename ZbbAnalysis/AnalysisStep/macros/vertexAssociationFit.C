#include <iostream>
#include <TH1.h>
#include <TPaveStats.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <cmath> 
#include <TLatex.h>
#include <TStyle.h>

// constants                
double xminfit[16] =  {0.1, 0.2 ,0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3,  0.3,  0.3,  0.3, 0.3, 0.3, 0.3,  0.2};   
double xmaxfit[16] =  {0.9, 0.9 ,0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.75, 0.75, 0.75, 0.8, 0.8, 0.8,  0.9}; 
unsigned int nbins = 15;
unsigned int nSteps = 20;
float cutmin = 0.;
float cutmax = 0.5;
float defaultCut = 0.15;
//float defaultCut=0.5;
float defaultPur = 0.99;
float defaultEff = 0.95;
bool constrainSlope = false;
double contrainedSlopeParameters[2] = { 4.443, -0.1367 };

Double_t getConstrainedMaximum(TH1* histo, double xmin, double xmax){

  int xfirst = histo->FindBin(xmin);
  int xlast  = histo->FindBin(xmax);
  
  double MAXVAL_ = -99999999.;
  int theBin     = -100;

  for(int bin=xfirst; bin<=xlast; bin++){
    double value = histo->GetBinContent(bin);
    if(value > MAXVAL_ ) {
      MAXVAL_=value;
      theBin = bin;
    }
  }
  return histo->GetBinCenter(theBin);
}

Double_t getConstrainedMinimum(TH1* histo, double xmin, double xmax){

  int xfirst = histo->FindBin(xmin);
  int xlast  = histo->FindBin(xmax);
  
  double MINVAL_ = 99999999.;
  int theBin     = -100;
  
  for(int bin=xfirst; bin<=xlast; bin++){
    double value = histo->GetBinContent(bin);
    if(value < MINVAL_ ) {
      MINVAL_=value;
      theBin = bin;
    }
  }
  return histo->GetBinCenter(theBin);
}

void MakeNiceStatBox(TH1* histo){

  TObject    *statObj;
  TPaveStats *stats;

  statObj = histo->GetListOfFunctions()->FindObject("stats");
  stats= static_cast<TPaveStats*>(statObj);
  stats->SetFillColor(10);
  stats->SetLineWidth(1);
  stats->SetShadowColor(0);
  stats->SetTextFont(42);
  stats->SetTextSize(0.05);
  //stats->SetLineColor(LineColors);
  //stats->SetTextColor(LineColors);
  stats->SetX1NDC(0.615);
  stats->SetY1NDC(0.747);
  stats->SetX2NDC(0.979);
  stats->SetY2NDC(0.986);
}

void MakeNiceHistoStyle(TH1 *hist){
  //hist->SetTitleSize(0.09); 
  //hist->SetMinimum(0.01);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42);  
  hist->SetLineWidth(2);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->CenterTitle(); 
  hist->GetXaxis()->SetLabelSize(0.07); 
  hist->GetYaxis()->SetLabelSize(0.07);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
}

void MakeNiceGraphStyle(TGraphErrors *graph){

  double x_,y_,maximum_(-1.), minimum_(1000.);

  for(UInt_t n=0;n<(UInt_t)graph->GetN();++n) {
    graph->GetPoint(n,x_,y_);
    if(y_>maximum_) maximum_ = y_;
    if(y_<minimum_) minimum_ = y_;
  }

  double delta = (maximum_ -minimum_);
  
  graph->SetMarkerColor(kBlack);
  graph->SetMarkerStyle(20);
  graph->SetLineWidth(2);
  graph->GetYaxis()->SetRangeUser(minimum_-delta,maximum_+delta);
  //graph->SetLineColor(kRed);
  graph->GetYaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetTitleOffset(1.0);
  graph->GetYaxis()->SetTitleSize(0.07);
  graph->GetXaxis()->SetTitleSize(0.07);
  graph->GetXaxis()->SetTitleFont(42); 
  graph->GetYaxis()->SetTitleFont(42);  
  graph->GetXaxis()->CenterTitle(); 
  graph->GetYaxis()->CenterTitle(); 
  graph->GetXaxis()->SetLabelSize(0.07); 
  graph->GetYaxis()->SetLabelSize(0.07);
  graph->GetXaxis()->SetLabelFont(42);
  graph->GetYaxis()->SetLabelFont(42);
}

void findWorkingPoint(TGraphErrors* graph, unsigned int& n, double& x, double& y, double& xe, double& ye,unsigned int known){
  switch(known) {
    // known 1 : fixed n
    case 1: {
      graph->GetPoint(n,x,y);
      xe = graph->GetErrorX(n);
      ye = graph->GetErrorY(n);
      break;
    }
    // known 2 : fixed x
    case 2: {
      /*
      double dist = 1000.;
      int best = n;
      double target = x;
      for(int i=0;i<graph->GetN();++i) {
        graph->GetPoint(i,x,y);
        if (fabs(target-x)<dist) {
          dist = fabs(target-x);
          best = i;
        }
      }
      n = best;
      graph->GetPoint(n,x,y);
      */
      double target = x;
      for(n=0;n<(UInt_t)graph->GetN();++n) {
        graph->GetPoint(n,x,y);
	xe = graph->GetErrorX(n);
	ye = graph->GetErrorY(n);
	if(x>target) break;
      }
      break;
    }
    // known 3 : fixed y
    case 3: {
      /*
      double dist = 1000.;
      int best = n;
      double target = y;
      for(int i=0;i<graph->GetN();++i) {
        graph->GetPoint(i,x,y);
        if (fabs(target-y)<dist) {
          dist = fabs(target-y);
          best = i;
        }
      }
      n = best;
      graph->GetPoint(n,x,y);
      */
      double target = y;
      for(n=0;n<(UInt_t)graph->GetN();++n) {
        graph->GetPoint(n,x,y);
	xe = graph->GetErrorX(n);
	ye = graph->GetErrorY(n);
	if(y<target) break;
      }
      break;
    }
  }
}

void estimateEfficiencyAndPurity(TH1* fraction,double cut,double& efficiency, double& purity, double& efferror, double& purerror, TFitResultPtr fitRes) {
  // efficiency loss is defined as the estimated signal below the cut over the estimated total signal
  double rangeLow, rangeHigh;
  TF1* expo = fraction->GetFunction("expo");
  expo->GetRange(rangeLow, rangeHigh);

  double signalLoss         = expo->Integral(0,cut)/fraction->GetBinWidth(1);
  double signalLossError = expo->IntegralError(0,cut,fitRes->GetParams(),fitRes->GetCovarianceMatrix().GetMatrixArray())/fraction->GetBinWidth(1);
  
  double signal_expopart         = expo->Integral(cut,fraction->GetBinLowEdge(fraction->FindBin(rangeHigh)+1))/fraction->GetBinWidth(1);
  double signal_expopartError = expo->IntegralError(cut, fraction->GetBinLowEdge(fraction->FindBin(rangeHigh)+1),fitRes->GetParams(),fitRes->GetCovarianceMatrix().GetMatrixArray())/fraction->GetBinWidth(1);

  double signal_toppoart = fraction->Integral(fraction->FindBin(rangeHigh)+1,fraction->FindBin(1.)+1);
  
  efficiency = (signal_expopart+signal_toppoart)/(signal_expopart+signal_toppoart+signalLoss);
  efferror = TMath::Sqrt( pow(signalLoss,2)*pow(signal_expopartError,2) + pow(signal_toppoart+signal_expopart,2)*pow(signalLossError,2) )/pow(signal_expopart+signal_toppoart+signalLoss,2);
  
  // purity is defined as the signal above the cut (signal_expopart+signal_toppoart) 
  // over the total data above that cut
  
  double data_kept = fraction->Integral(fraction->FindBin(cut),fraction->FindBin(1.)+1);
  double allBeforeRangeHigh = fraction->Integral(fraction->FindBin(cut),fraction->FindBin(rangeHigh));
  purity = (signal_expopart+signal_toppoart)/data_kept;
  purerror = signal_expopartError/data_kept; 
 
  if(purity > 1.){
    std::cout << "#################################################"<<std::endl;
    std::cout << "cut: "<<cut<<" purity:" << purity <<std::endl;
    std::cout << "signalLoss: " << signalLoss << " signalExtrapolated: " << signal_expopart << " signalTop: " << signal_toppoart << std::endl;
    std::cout << "all b.r.h.: " << allBeforeRangeHigh << std::endl;
    std::cout << "total signal with expo: "<< (signal_expopart+signal_toppoart)<< std::endl;
    std::cout << "total signal integrating: "<< (allBeforeRangeHigh+signal_toppoart)<< std::endl;
    std::cout << "data kept:  " << data_kept  << std::endl;
    std::cout << "ratio:      " << allBeforeRangeHigh/signal_expopart << std::endl;
    std::cout << "#################################################"<<std::endl;
  }
 
  //   std::cout << "estimateEfficiencyAndPurity for cut=" << cut << std::endl;
  //   std::cout << "signal: " << signalLoss << " " << signal_expopart << " " << signal_toppoart << std::endl;
  //   std::cout << "data kept: " << data_kept << std::endl;
  //   std::cout << "eff: " << efficiency << "+/- " << efferror << " pur: " << purity << " +/- " << purerror << std::endl;
}

void vertexAssociationAnalysis(TH1* fraction, unsigned int index, double& slope,double& slopeError,double& constNorm,double& constError,TGraphErrors*& graph,TFitResultPtr fitRes) {
  // do the fit
  if(constrainSlope) {
    fraction->Fit("expo","LS","",xminfit[index],xmaxfit[index]);
    TF1* expo = (TF1*)gROOT->FindObjectAny("expo");
    expo->FixParameter(1,contrainedSlopeParameters[0]+contrainedSlopeParameters[1]*(index+1));
    fitRes = fraction->Fit("expo","LSB","",xminfit[index],xmaxfit[index]);
    expo->SetParLimits(1,-10000,10000);
  } else {
    // fitRes = fraction->Fit("expo","LS","",xminfit[index],xmaxfit[index]);  
    fitRes = fraction->Fit("expo","LS","",getConstrainedMinimum(fraction,0.1,0.3),getConstrainedMaximum(fraction,0.7,0.9));
  }
  // estimate efficiency and purity as a function of the cut
  double *efficiency, *purity, *efferror, *purerror;
  efficiency = new double[nSteps];
  purity     = new double[nSteps];
  efferror   = new double[nSteps];
  purerror   = new double[nSteps];
  unsigned int step = 0;
  for(float cut=cutmin; cut<cutmax; cut += (cutmax-cutmin)/float(nSteps), ++step) {
    estimateEfficiencyAndPurity(fraction,cut,efficiency[step],purity[step],efferror[step],purerror[step],fitRes);
  }
  // output: the graph with efficiency and purity vs cut and the expo slope
  slope      = fraction->GetFunction("expo")->GetParameter(1);
  slopeError = fraction->GetFunction("expo")->GetParError(1);
  
  constNorm  = fraction->GetFunction("expo")->GetParameter(0);
  constError = fraction->GetFunction("expo")->GetParError(0);

  //slope = result->Parameter(1);
  graph = new TGraphErrors(nSteps,purity,efficiency,purerror,efferror);
  //graph = new TGraph(nSteps,purity,efficiency);   
  delete[] efficiency;
  delete[] purity;
  delete[] efferror;
  delete[] purerror;
}

void vertexAssociationFit(TFile* input, TFile* output) {
  
  // style
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetOptTitle(1);
  // margins
  gStyle->SetStatFormat("6.4g"); 
  gStyle->SetFitFormat("5.2g"); 
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.02);

  // prepare the output
  // * canvas to put the fit per bin
  TCanvas* canvas_0 = new TCanvas("FitPerBin","Fit for each vertex bin",1200,900);
  canvas_0->Divide(4,4);
  // * canvas to put the plots per bin
  TCanvas* canvas_1 = new TCanvas("EffPerBin","Efficiency vs purity for each vertex bin",1200,900);
  canvas_1->Divide(4,4);
  // * canvas with other plots: slope vs nvertices, efficiency and purity for fixed cut, etc.
  TCanvas* canvas_2 = new TCanvas("jvAssociation","Efficiency vs purity summary plots",1000,600);
  canvas_2->Divide(4,2);
  // summary graphs
  TGraphErrors* slopeGraph  = new TGraphErrors(nbins);
  slopeGraph->SetNameTitle("slopeGraph","expo slope in each n_PU bin");
  TGraphErrors* constGraph  = new TGraphErrors(nbins);
  constGraph->SetNameTitle("constGraph","expo normalization in each n_PU bin");
  TGraphErrors* cuteffGraph = new TGraphErrors(nbins);
  cuteffGraph->SetNameTitle("cuteffGraph","Efficiency in each n_PU bin");
  TGraphErrors* cutpurGraph = new TGraphErrors(nbins);
  cutpurGraph->SetNameTitle("cutpurGraph","Purity in each n_PU bin");
  TGraphErrors* effcutGraph = new TGraphErrors(nbins);
  effcutGraph->SetNameTitle("effcutGraph","Cut for fixed efficiency in each n_PU bin");
  TGraphErrors* effpurGraph = new TGraphErrors(nbins);
  effpurGraph->SetNameTitle("effpurGraph","Purity for fixed efficiency in each n_PU bin");
  TGraphErrors* purcutGraph = new TGraphErrors(nbins);
  purcutGraph->SetNameTitle("purcutGraph","Cut for fixed purity in each n_PU bin");
  TGraphErrors* pureffGraph = new TGraphErrors(nbins);
  pureffGraph->SetNameTitle("pureffGraph","Efficiency for fixed purity in each n_PU bin");
  // loop over the directories
  double slope = 0.;
  double slopeError = -999.;

  double constnorm = 0.;
  double constError = -999.;

  TGraphErrors* result = NULL;
  TFitResultPtr fitRes(-1);
  double x,y, xe,ye;
  unsigned int index;
  TString histonames[nbins+1];
  TString labelnames[nbins+1];


  // get the histogram    
  TH1* evts_per_PUbin = (TH1*)input->Get("finaldistros_ssvhpt/GoodJet/GoodJet_nvertices");
  double sumTot(0), effTot(0);

  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.05);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();
   
  TLatex *latexLabel2 = new TLatex();
  latexLabel2->SetTextSize(0.09);
  latexLabel2->SetTextFont(42);
  latexLabel2->SetLineWidth(2);
  latexLabel2->SetNDC();
   
  for(unsigned int i=0; i<=nbins; ++i) {
    // get the histogram  
    labelnames[i]=Form("n_PU = %d bin",i);
    std::cout<<"labelnames["<<i<<"]"<<labelnames[i]<<std::endl;
    histonames[i]=Form("finaldistros_ssvhpt/GoodJet/GoodJet_jetvtx_ratio2b_%d_PUvtx",i);
    std::cout<<"seizing: histonames["<<i<<"]:"<<histonames[i]<<std::endl;
    TH1* fraction = (TH1*)input->Get(histonames[i]);
    //fraction->Rebin(2);
    // do it
    canvas_0->cd(i+1);
    fraction->SetTitle("");
    MakeNiceHistoStyle(fraction);
    if(fraction) vertexAssociationAnalysis(fraction,i,slope,slopeError,constnorm,constError,result,fitRes);
    latexLabel2->DrawLatex(0.15, 0.85,labelnames[i]);
    latexLabel2->DrawLatex(0.15, 0.75,"CMS #sqrt{s}= 7 TeV");
    gPad->Draw();
    gPad->SetLogy();
    MakeNiceStatBox(fraction);
    output->cd();
    fraction->Write();
    // plot
    canvas_1->cd(i+1);
    result->SetTitle("");
    MakeNiceGraphStyle(result);
    result->Draw("alP");
    latexLabel->DrawLatex(0.05, 0.95,labelnames[i]);
    result->GetXaxis()->SetTitle("Purity");
    result->GetYaxis()->SetTitle("Efficiency");
    result->GetXaxis()->SetRangeUser(0.8,1.2);
    result->GetYaxis()->SetRangeUser(0.85,1.0);
    output->cd();
    result->Write(); 
    // extract other quantities
    slopeGraph->SetPoint(i,i,slope);
    slopeGraph->SetPointError(i,0,slopeError);
    slopeGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    constGraph->SetPoint(i,i,constnorm);
    constGraph->SetPointError(i,0,constError);
    constGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    index = unsigned(((defaultCut-cutmin)/(cutmax-cutmin)*nSteps)+0.5);
    findWorkingPoint(result,index,x,y,xe,ye,1);

    // to compute the weighted efficiency
    int ivtx_bin=evts_per_PUbin->FindBin(i+1); 
    sumTot+=evts_per_PUbin->GetBinContent(ivtx_bin);
    effTot+=y*evts_per_PUbin->GetBinContent(ivtx_bin);

    cuteffGraph->SetPoint(i,i,y);
    cuteffGraph->SetPointError(i,0,ye);
    cuteffGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    cutpurGraph->SetPoint(i,i,x);
    cutpurGraph->SetPointError(i,0,xe);
    cutpurGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    x = defaultPur;
    findWorkingPoint(result,index,x,y,xe,ye,2);
    purcutGraph->SetPoint(i,i,cutmin+(cutmax-cutmin)/nSteps*index);
    purcutGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    pureffGraph->SetPoint(i,i,y);
    pureffGraph->SetPointError(i,0,ye);
    pureffGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    y = defaultEff;
    findWorkingPoint(result,index,x,y,xe,ye,3);
    effcutGraph->SetPoint(i,i,cutmin+(cutmax-cutmin)/nSteps*index);
    effcutGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    effpurGraph->SetPoint(i,i,x);
    effpurGraph->SetPointError(i,0,xe);
    effpurGraph->GetXaxis()->SetTitle("n_{PU} vertices");
    std::cout<<"end of other quantities: bin["<<i<<"]"<<std::endl;
  }

  std::cout << "Weighted Efficiency" << effTot/sumTot << std::endl;

  // plot summary graphs
  TF1* f=NULL;
  canvas_2->cd(1);
  MakeNiceGraphStyle(slopeGraph);
  slopeGraph->Draw("alP");
  slopeGraph->Fit("pol1");
  f = slopeGraph->GetFunction("pol1");
  if(f->GetParError(1)>fabs(f->GetParameter(1))) slopeGraph->Fit("pol0");
  
  canvas_2->cd(2);
  MakeNiceGraphStyle(constGraph);
  constGraph->Draw("alP");
  constGraph->Fit("pol2");
  f = constGraph->GetFunction("pol2");
  
  canvas_2->cd(3);
  MakeNiceGraphStyle(cuteffGraph);
  cuteffGraph->Draw("alP");
  cuteffGraph->Fit("pol1");
  f = cuteffGraph->GetFunction("pol1");
  if(f->GetParError(1)>fabs(f->GetParameter(1))) cuteffGraph->Fit("pol0");
  latexLabel2->DrawLatex(0.15, 0.15,Form("Jet Fraction cut: %.2f",defaultCut));
  
  canvas_2->cd(4);
  MakeNiceGraphStyle(cutpurGraph);
  cutpurGraph->Draw("alP");
  cutpurGraph->Fit("pol1");
  f = cutpurGraph->GetFunction("pol1");
  if(f->GetParError(1)>fabs(f->GetParameter(1))) cutpurGraph->Fit("pol0");
  latexLabel2->DrawLatex(0.15, 0.15,Form("Jet Fraction cut: %.2f",defaultCut));
  
  canvas_2->cd(5);
  MakeNiceGraphStyle(effcutGraph);
  effcutGraph->Draw("alP");
  effcutGraph->Fit("expo");
  latexLabel2->DrawLatex(0.15, 0.15,Form("Fixed efficiency: %.2f",defaultEff));
  
  canvas_2->cd(6);
  MakeNiceGraphStyle(effpurGraph);
  effpurGraph->Draw("alP");
  effpurGraph->Fit("pol1");
  f = effpurGraph->GetFunction("pol1");
  if(f->GetParError(1)>fabs(f->GetParameter(1))) effpurGraph->Fit("pol0");
  latexLabel2->DrawLatex(0.15, 0.15,Form("Fixed efficiency: %.2f",defaultEff));
  
  canvas_2->cd(7);
  gPad->SetLeftMargin(0.12);
  MakeNiceGraphStyle(purcutGraph);
  purcutGraph->Draw("alP");
  purcutGraph->Fit("expo");
  latexLabel2->DrawLatex(0.15, 0.15,Form("Fixed purity: %.3f",defaultPur));
  
  canvas_2->cd(8);
  MakeNiceGraphStyle(pureffGraph);
  pureffGraph->Draw("alP");
  pureffGraph->Fit("pol1");
  f = pureffGraph->GetFunction("pol1");
  if(f->GetParError(1)>fabs(f->GetParameter(1))) pureffGraph->Fit("pol0");
  latexLabel2->DrawLatex(0.15, 0.15,Form("Fixed purity: %.3f",defaultPur));
  
  // save
  output->cd();
  canvas_0->Write();
  canvas_1->Write();
  canvas_2->Write();

  canvas_0->SaveAs((TString)(canvas_0->GetName())+".png");
  canvas_1->SaveAs((TString)(canvas_1->GetName())+".png");
  canvas_2->SaveAs((TString)(canvas_2->GetName())+".png");

}

void vertexAssociationFit(const char* filename, const char* output="output.root") {
  TFile* file = TFile::Open(filename);
  TFile* fileOut = TFile::Open(output,"recreate");
  vertexAssociationFit(file,fileOut);
  //file->Close();
  //fileOut->Close();
}


