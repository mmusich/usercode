#define JetByJetComparison_cxx
#include "JetByJetComparison.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TPaveStats.h>
#include <TLatex.h>

void JetByJetComparison::Loop()
{

  //gROOT->SetStyle("Plain");
  setTDRStyle();
  
  TH1F* hDeltapt		  = new TH1F("hDeltapt","#Delta p_{T};jet #Deltap_{T} (GeV);jets",100,-5.,5.);		
  TH1F* hDeltaeta		  = new TH1F("hDeltaeta","#Delta #eta;jet #Delta#eta;jets",100,-5.,5.);		
  TH1F* hDeltaphi		  = new TH1F("hDeltaphi","#Delta #phi;jet #Delta#phi (rad);jets",100,-5.,5.);		
  TH1F* hDeltajetId		  = new TH1F("hDeltajetId","#Delta jetId;jet #Delta(jetId);jets",8,-4,4);		
  TH1F* hDeltaMCTrueFlavor	  = new TH1F("hDeltaMCTrueFlavor","#Delta MC true flavour;jet #Delta(MC true flavour);jets",40,-20,20);	
  TH1F* hDeltaSV3dDistance	  = new TH1F("hDeltaSV3dDistance","#Delta SV 3d distance;#Delta SV 3d distance (cm);secondary vertices",100,-5.,5.);	
  TH1F* hDeltaSV3dDistanceError   = new TH1F("hDeltaSV3dDistanceError","#Delta SV 3d distance error;#Delta SV 3d distance error (cm);secondary vertices",100,-5.,5.); 
  TH1F* hDeltaSV2dDistance	  = new TH1F("hDeltaSV2dDistance","#Delta SV 2d distance;#Delta SV 2d distance (cm);secondary vertices",100,-5.,5.);	
  TH1F* hDeltaSV2dDistanceError   = new TH1F("hDeltaSV2dDistanceError","#Delta SV 2d distance error;#Delta SV 2d distance (cm);secondary vertices",100,-5.,5.);
  TH1F* hDeltaSVChi2		  = new TH1F("hDeltaSVChi2","#Delta SV #chi^{2};#Delta SV #chi^{2};secondary vertices",100,-5.,5.);	
  TH1F* hDeltaSVDegreesOfFreedom  = new TH1F("hDeltaSVDegreesOfFreedom","#Delta SV d.o.f.;#Delta SV d.o.f.;secondary vertices",10,-5.,5.);
  TH1F* hDeltaSVNormChi2	  = new TH1F("hDeltaSVNormChi2","#Delta SV #chi^{2}/ndof;#Delta SV #chi^{2}/ndof;secondary vertices",100,-5.,5.);		
  TH1F* hDeltaSVMass		  = new TH1F("hDeltaSVMass","#Delta SV mass;#Delta SV mass (GeV);secondary vertices",100,-5.,5.);			
  TH1F* hDeltaSVtotCharge	  = new TH1F("hDeltaSVtotCharge","#Delta SV total charge;#Delta SV total charge;secondary vertices",10,-5.,5.);		        
  TH1F* hDeltaSVnVertices	  = new TH1F("hDeltaSVnVertices","#Delta SV number of vertices;#Delta SV N_{vertices};secondary vertices",10,-5.,5.);		
  TH1F* hDeltaSVnVertexTracks	  = new TH1F("hDeltaSVnVertexTracks","#Delta SV number of tracks;#Delta SV N_{tracks};secondary vertices",10,-5.,5.);		
  TH1F* hDeltaSVnVertexTracksAll  = new TH1F("hDeltaSVnVertexTracksAll","#Delta SV number of all tracks;#Delta SV N^{all}_{tracks};secondary vertices",10,-5.,5.);	
  TH1F* hDeltaDiscrTCHE           = new TH1F("hDeltaDiscrTCHE","#Delta DiscrTCHE ;#Delta DiscrTCHE;jets",1000,-5.,5.);
  TH1F* hDeltaDiscrTCHP           = new TH1F("hDeltaDiscrTCHP","#Delta DiscrTCHP ;#Delta DiscrTCHP;jets",1000,-5.,5.);
  TH1F* hDeltaDiscrSSVHE          = new TH1F("hDeltaDiscrSSVHE","#Delta DiscrSSVHE ;#Delta DiscrSSVHE;jets",1000,-5.,5.);
  TH1F* hDeltaDiscrSSVHP          = new TH1F("hDeltaDiscrSSVHP","#Delta DiscrSSVHP ;#Delta DiscrSSVHP;jets",1000,-5.,5.);
  TH1F* hDeltaDiscrCSV            = new TH1F("hDeltaDiscrCSV","#Delta DiscrCSV ;#Delta DiscrCSV;jets",1000,-5.,5.);
  TH1F* hDeltaDiscrJP             = new TH1F("hDeltaDiscrJP","#Delta DiscrJP ;#Delta DiscrJP;jets",1000,-5.,5.);
  TH1F* hDeltaDiscrJBP            = new TH1F("hDeltaDiscrJBP","#Delta DiscrJBP ;#Delta DiscrJBP;jets",1000,-5.,5.);

  std::vector<TH1F*> histosvec;

  //   In a ROOT session, you can do:
  //      Root > .L JetByJetComparison.C
  //      Root > JetByJetComparison t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //  This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //  To read only selected branches, Insert statements like:
  //  METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  //  METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //    by  b_branchname->GetEntry(ientry); //read only this branch
  
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // Fill histograms
      hDeltaDiscrTCHE->Fill(tche[0]-tche[1]);
      hDeltaDiscrTCHP->Fill(tchp[0]-tchp[1]);
      hDeltaDiscrSSVHE->Fill(ssvhe[0]-ssvhe[1]);
      hDeltaDiscrSSVHP->Fill(ssvhp[0]-ssvhp[1]); 
      hDeltaDiscrCSV->Fill(csv[0]-csv[1]);
      hDeltaDiscrJP->Fill(jp[0]-jp[1]);
      hDeltaDiscrJBP->Fill(jbp[0]-jbp[1]);
      hDeltapt->Fill(pt[0]-pt[1]);		 
      hDeltaeta->Fill(eta[0]-eta[1]);		 
      hDeltaphi->Fill(phi[0]-phi[1]);		 
      hDeltajetId->Fill(jetId[0]-jetId[1]);		 
      hDeltaMCTrueFlavor->Fill(MCTrueFlavor[0]-MCTrueFlavor[1]);
      hDeltaSV3dDistance->Fill(SV3dDistance[0]-SV3dDistance[1]);
      hDeltaSV3dDistanceError->Fill(SV3dDistanceError[0]-SV3dDistanceError[1]);
      hDeltaSV2dDistance->Fill(SV2dDistance[0]-SV2dDistance[1]);
      hDeltaSV2dDistanceError->Fill(SV2dDistanceError[0]-SV2dDistanceError[1]);
      hDeltaSVChi2->Fill(SVChi2[0]-SVChi2[1]);
      hDeltaSVDegreesOfFreedom->Fill(SVDegreesOfFreedom[0]-SVDegreesOfFreedom[1]);
      hDeltaSVNormChi2->Fill(SVNormChi2[0]-SVNormChi2[1]);
      hDeltaSVMass->Fill(SVMass[0]-SVMass[1]);
      hDeltaSVtotCharge->Fill(SVtotCharge[0]-SVtotCharge[1]);
      hDeltaSVnVertices->Fill(SVnVertices[0]-SVnVertices[1]);
      hDeltaSVnVertexTracks->Fill(SVnVertexTracks[0]-SVnVertexTracks[1]);
      hDeltaSVnVertexTracksAll->Fill(SVnVertexTracksAll[0]-SVnVertexTracksAll[1]);

   }

   histosvec.push_back(hDeltaDiscrTCHE);
   histosvec.push_back(hDeltaDiscrTCHP); 
   histosvec.push_back(hDeltaDiscrSSVHE);
   histosvec.push_back(hDeltaDiscrSSVHP);
   histosvec.push_back(hDeltaDiscrCSV);  
   histosvec.push_back(hDeltaDiscrJP);   
   histosvec.push_back(hDeltaDiscrJBP); 
   histosvec.push_back(hDeltapt);		 
   histosvec.push_back(hDeltaeta);		 
   histosvec.push_back(hDeltaphi);		 
   histosvec.push_back(hDeltajetId);		 
   histosvec.push_back(hDeltaMCTrueFlavor);	 
   histosvec.push_back(hDeltaSV3dDistance);	 
   histosvec.push_back(hDeltaSV3dDistanceError);  
   histosvec.push_back(hDeltaSV2dDistance);	 
   histosvec.push_back(hDeltaSV2dDistanceError);  
   histosvec.push_back(hDeltaSVChi2);		 
   histosvec.push_back(hDeltaSVDegreesOfFreedom); 
   histosvec.push_back(hDeltaSVNormChi2);	 
   histosvec.push_back(hDeltaSVMass);		 
   histosvec.push_back(hDeltaSVtotCharge);	 
   histosvec.push_back(hDeltaSVnVertices);	 
   histosvec.push_back(hDeltaSVnVertexTracks);	 
   histosvec.push_back(hDeltaSVnVertexTracksAll); 

   TFile *file_out=new TFile("JetByJetComparisonPlots.root","recreate");  
   file_out->cd();
   
   UInt_t nOfHistos = histosvec.size();

   TObject *statObj[nOfHistos];
   TPaveStats *stats[nOfHistos];

   for(UInt_t h=0; h<histosvec.size(); h++){
     TCanvas *c = new TCanvas(histosvec[h]->GetName(),histosvec[h]->GetName(),600,600);
     c->cd()->SetLogy();
     histosvec[h]->Draw();
     c->Draw();
   
     statObj[h] = histosvec[h]->GetListOfFunctions()->FindObject("stats");
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

     file_out->cd();
     cmsPrel(0.);
     histosvec[h]->Write();

   }
   
   file_out->cd();
   file_out->Write();
   file_out->Close();
   
}

void JetByJetComparison::setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

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

void JetByJetComparison::cmsPrel(const double& intLumi) {

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22

  latex->SetTextAlign(13);
  latex->DrawLatex(0.12, 0.99, Form("CMS Preliminary 2011,     #sqrt{s} = 7 TeV,  L = %.2g pb^{ -1}",intLumi));
  latex->DrawLatex(0.20, 0.90, CompNames[0]+" vs "+CompNames[1]);

}
